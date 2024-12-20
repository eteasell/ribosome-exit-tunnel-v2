import subprocess
from Bio import AlignIO
from Bio.PDB.MMCIFParser import MMCIFParser
from pathlib import Path
from protocol.data_access import *
from protocol.domain.landmark import *
from protocol.domain.types import PROTOTYPES
from protocol.taxonomy import *
from protocol.domain.sequence import *
from protocol.locate_residues import *

def process_list(list_rcsb_id: list[str], universal: bool = False):
    processed_ids = []
    changed_files = []
    for rcsb_id in list_rcsb_id:
        if universal is True:
            kingdom = None
            prototype = UNIVERSAL_PROTOTYPE
        else:
            kingdom = find_kingdom(rcsb_id)
            prototype = PROTOTYPES[kingdom]
        
        path = f"../ribosome-exit-tunnel-v2/data/mmcif/{rcsb_id}.cif"
        if not Path(path).is_file():
            success = get_mmcif(rcsb_id)
            if success is False:
                print(f"Cannot access {rcsb_id} file.")
                continue
            
        profile = None
        for polymer in prototype.keys():
            add_to_processed = True
            find_seq = check_fasta_for_rcsb_id(rcsb_id, polymer, kingdom)
            if find_seq is None:
                if profile is None:
                    profile = get_profile(rcsb_id)
                
                if profile is None:
                    add_to_processed = False
                    continue
                    
                changed_file = add_polymer_to_fasta_list(rcsb_id, polymer, profile, kingdom)
                
                if changed_file is None:
                    add_to_processed = False
                elif changed_file not in changed_files:
                    changed_files.append(changed_file)
        
        if add_to_processed:
            processed_ids.append(rcsb_id)
                               
    return processed_ids, changed_files
    
def align(files: list[str]):
    
    for file in files:
        input = f"data/fasta/{file}.fasta"
        output = f"data/output/fasta/aligned_{file}.fasta"
        mafft_command = f"mafft --auto --leavegappyregion {input} > {output}"
        subprocess.call(mafft_command, shell=True)


def select_landmarks(conservation_threshold: float, distance_threshold: float, kingdom: str = None):
    
    if kingdom is None:
        prototype = UNIVERSAL_PROTOTYPE
    else:
        prototype = PROTOTYPES[kingdom]
    proto_id = prototype['uL4']['parent_id']
    
    parser = MMCIFParser(QUIET=True)
    path = f"../ribosome-exit-tunnel-v2/data/mmcif/{proto_id}.cif"
    structure = parser.get_structure(proto_id, path)
    
    total_conserved = []
    
    for polymer in prototype.keys():
        print(f"Starting polymer {polymer}...")
        
        asym_id = prototype[polymer]['auth_asym_id']
        
        if kingdom is None:
            path = f"data/output/fasta/aligned_sequences_{polymer}.fasta"
        else:
            path = f"data/output/fasta/aligned_sequences_{kingdom}_{polymer}.fasta"

        alignment = AlignIO.read(path, "fasta")

        conserved_positions = []

        for i in range(alignment.get_alignment_length()):
            most_common = find_conserved(alignment[:,i], conservation_threshold)
    
            if (most_common is not None):
                conserved_positions.append(Landmark(i,most_common))
           
        proto_seq = None
         
        # Find prototype in alignment
        for i, seq in enumerate(alignment):
            if proto_id.upper() in seq.name:
                proto_seq = seq
                
        if proto_seq is None:
            print(f"Cannot find {polymer} prototype {proto_id}")
            continue
        
        chain = structure.child_dict[0][asym_id]
        chain_container = SequenceMappingContainer(chain)
        flat_seq = chain_container.flat_sequence
        
        conserved_positions = cherry_pick(polymer, proto_seq, conserved_positions, distance_threshold, chain, flat_seq, kingdom)

        total_conserved = total_conserved + conserved_positions
        
    return total_conserved

# given a full list of landmarks (for all polymers), and a list of those polymers,
# returns a list of landmark location for this rcsb_id accross all landmarks
def locate_landmarks(rcsb_id: str, conserved_positions: list[Landmark], polymers: list[str], kingdom : str = None):
    rows = []
    
    path = f"../ribosome-exit-tunnel-v2/data/mmcif/{rcsb_id}.cif"
    if not Path(path).is_file():
        get_mmcif(rcsb_id)
    
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(rcsb_id, path)
    except:
        print(f"BioPython MMCIFParser cannot handle {rcsb_id}.cif file.")
        return
    
    for polymer in polymers:
        
        if kingdom is None:
            path = f"data/output/fasta/aligned_sequences_{polymer}.fasta"
        else:
            path = f"data/output/fasta/aligned_sequences_{kingdom}_{polymer}.fasta"
            
        alignment = AlignIO.read(path, "fasta")
    
        seq = get_rcsb_in_alignment(alignment, rcsb_id)
        if seq is None:
            print("rcsb_id not found in alignment.")
            continue
    
        name_arr = seq.name.split('_')
        parent = name_arr[1]
        asym_id = name_arr[2]
        
        try:
            chain = structure.child_dict[0][asym_id]
        except:
            print(f"Chain {asym_id} not found in {rcsb_id} PDB file")
            continue
        
        chain_container = SequenceMappingContainer(chain)
        flat_seq = chain_container.flat_sequence
    
        i = 0
        for pos in conserved_positions:
            if polymer not in pos.name: continue
        
            landmark = Landmark(pos.position, pos.residue, pos.name)

            coords = locate_residues(landmark, polymer, asym_id, parent, chain, flat_seq, kingdom)
        
            if coords is not None:
                rows.append(coords)
                
            i += 1
            
    return rows
    