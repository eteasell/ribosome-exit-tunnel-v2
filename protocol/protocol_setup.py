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


def align(kingdom: str):
    
    prototype = PROTOTYPES[kingdom]
    
    for polymer in prototype.keys():
        print(f"Starting polymer {polymer}...")
        # Perform multiple sequence alignment   
        input = f"data/fasta/sequences_{kingdom}_{polymer}.fasta"
        output = f"data/output/fasta/aligned_sequences_{kingdom}_{polymer}.fasta"
        mafft_command = f"mafft --auto --leavegappyregion {input} > {output}"
        subprocess.call(mafft_command, shell=True)
    

def select_landmarks(kingdom: str, conservation_threshold: float, distance_threshold: float):
    
    prototype = PROTOTYPES[kingdom]
    proto_id = prototype['uL4']['parent_id']
    
    parser = MMCIFParser(QUIET=True)
    path = f"../ribosome-exit-tunnel/data/mmcif/{proto_id}.cif"
    structure = parser.get_structure(proto_id, path)
    
    total_conserved = []
    
    for polymer in prototype.keys():
        print(f"Starting polymer {polymer}...")
        
        asym_id = prototype[polymer]['auth_asym_id']

        alignment = AlignIO.read(f"data/output/fasta/aligned_sequences_{kingdom}_{polymer}.fasta", "fasta")

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
        
        # TODO: may want to be careful with model 0 choice here
        chain = structure.child_dict[0][asym_id]
        chain_container = SequenceMappingContainer(chain)
        flat_seq = chain_container.flat_sequence
        
        conserved_positions = cherry_pick(polymer, proto_seq, conserved_positions, distance_threshold, kingdom, chain, flat_seq)

        total_conserved = total_conserved + conserved_positions
        
    return total_conserved

# given a full list of landmarks (for all polymers), and a list of those polymers,
# returns a list of landmark location for this rcsb_id accross all landmarks
def locate_landmarks(rcsb_id: str, kingdom: str, conserved_positions: list[Landmark], polymers: list[str]):
    rows = []
    
    path = f"../ribosome-exit-tunnel/data/mmcif/{rcsb_id}.cif"
    if not Path(path).is_file():
        get_mmcif(rcsb_id)
    
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(rcsb_id, path)
    
    for polymer in polymers:
        
        alignment = AlignIO.read(f"data/output/fasta/aligned_sequences_{kingdom}_{polymer}.fasta", "fasta")
    
        seq = get_rcsb_in_alignment(alignment, rcsb_id)
        if seq is None:
            print("rcsb_id not found in alignment.")
            return None 
    
        name_arr = seq.name.split('_')
        parent = name_arr[1]
        asym_id = name_arr[2]
        
        chain = structure.child_dict[0][asym_id]
        chain_container = SequenceMappingContainer(chain)
        flat_seq = chain_container.flat_sequence
    
        i = 0
        for pos in conserved_positions:
            if polymer not in pos.name: continue
        
            landmark = Landmark(pos.position, pos.residue, f'{polymer}-{i}')

            coords = locate_residues(landmark, polymer, asym_id, parent, kingdom, chain, flat_seq)
        
            if coords is not None:
                rows.append(coords)
                
            i += 1
            
    return rows
    