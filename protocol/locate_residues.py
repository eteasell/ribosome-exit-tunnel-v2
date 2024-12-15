from protocol.data_access import *
from protocol.domain.sequence import *
from protocol.domain.landmark import *
from protocol.tunnel_coordinates import find_closest_point
import subprocess
from Bio import AlignIO, Seq
from Bio.PDB import Structure
from collections import Counter
from io import StringIO
import numpy as np

pairwise_alignment_cache = {}

def locate_residues(landmark: Landmark, 
                    polymer: str, 
                    polymer_id: str, 
                    rcsb_id: str,  
                    chain: Structure, 
                    flat_seq,
                    kingdom: str = None) -> dict:
    
    '''
    This method takes a landmark centered on the alignment, and finds this residue on the given rcsb_id's polymer.
    Returns the residues position and coordinates.
    
    landmark: the landmark to be located
    polymer: the poymer on which this landmark lies
    polymer_id: the polymer id specific to this rcsb_id
    rcsb_id: the id of the ribosome instance
    chain: the biopython Chain object holding the sequence
    flat_seq: from SequenceMappingContainer, tuple holding (seq, flat_index_to_residue_map, auth_seq_id_to_flat_index_map)
    kingdom: kingdom to which this rcsb_id belongs, or none if being called from main_universal.py
    '''
    
    # access aligned sequence from alignment files
    if kingdom is None:
        path = f"data/output/fasta/aligned_sequences_{polymer}.fasta"
    else:
        path = f"data/output/fasta/aligned_sequences_{kingdom}_{polymer}.fasta"
    alignment = AlignIO.read(path, "fasta")
    aligned_seq = get_rcsb_in_alignment(alignment, rcsb_id)
    
    # find the position of the landmark on the original riboXYZ seq
    alignment_position = map_to_original(aligned_seq, landmark.position) 
    
    # access riboXYZ sequence (pre alignment)
    orig_seq = check_fasta_for_rcsb_id(rcsb_id, polymer, kingdom)

    if orig_seq is None:
        print("Cannot access uniprot sequence")
        return
    
    # run pairwise alignment on the riboXYZ sequence and the flattened PDB sequence
    alignment = run_pairwise_alignment(rcsb_id, polymer_id, orig_seq, flat_seq[0])
    
    if alignment is None:
        return None
        
    # map the alignment_position from the original riboXYZ sequence to the pairwise-aligned flattened PDB sequence
    flattened_seq_aligned = alignment[1]
    flat_aligned_position = None
    if alignment_position is not None:  
        flat_aligned_position = map_to_original(flattened_seq_aligned, alignment_position)
    
    if flat_aligned_position is None:
        print(f"Cannot find {landmark} on {rcsb_id} {polymer}")
        return None 
        
    # use the MappingSequenceContainer flat_index_to_residue_map to access to residue in the PDB sequence
    resi_id = flat_seq[1][flat_aligned_position].get_id()
    residue = chain[resi_id]
    
    # check that the located residue is the same as the landmark
    landmark_1_letter = landmark.residue.upper()
    landmark_3_letter = ResidueSummary.one_letter_code_to_three(landmark_1_letter)
    if (residue.get_resname() != landmark_1_letter and residue.get_resname() != landmark_3_letter):
        return None
        
    # find atomic coordinates for the selected residue
    atom_coords = [atom.coord for atom in residue]
    if (len(atom_coords) == 0): 
       return None
    
    # take the mean coordinate for the atoms in residue
    vec = np.zeros(3)
    for coord in atom_coords:
        tmp_arr = np.array([coord[0], coord[1], coord[2]])
        vec += tmp_arr
    vec = vec / len(atom_coords)
    vec = vec.astype(np.int32)
        
    return {    
                "parent_id": rcsb_id, 
                "landmark": landmark.name, 
                "residue": landmark.residue, 
                "position": resi_id[1],
                "x": vec[0], "y": vec[1], "z": vec[2]
            }
    
def run_pairwise_alignment(rcsb_id, polymer_id, seq1, seq2):
    
    cache_key = f"{rcsb_id}-{polymer_id}"
    if cache_key in pairwise_alignment_cache:
        return pairwise_alignment_cache[cache_key]
        
    fasta = f">seq1\n{seq1}\n>seq2\n{seq2}"

    # Run MAFFT command
    process = subprocess.run(
    ["mafft", "--auto", "-"],
    input=fasta,
    text=True,
    capture_output=True
    )

    # Read the MAFFT output directly from stdout into a MultipleSeqAlignment object
    try:
        alignment =  AlignIO.read(StringIO(process.stdout), "fasta")
    except:
        print(f"Issue with pairwise sequence alignment on {rcsb_id}")
        pairwise_alignment_cache[cache_key] = None
        return None
    
    pairwise_alignment_cache[cache_key] = alignment
    
    return alignment

    
def cherry_pick(polymer: str, 
                seq: Seq, 
                conserved_positions: list[Landmark], 
                threshold: float,
                chain: Structure,
                flat_seq,
                kingdom : str = None
                ) -> list[Landmark]:
    
    '''
    This method choses landmarks based on distance to the MOLE model tunnel centerline on the prototype specified by seq
    
    polymer: the chain of interest (eg. uL4)
    seq: a Seq object which contains the protoype aligned sequence
    conserved_positions: list of Landmark candidates
    threshold: conservation threshold between [0,1]
    chain: the chain sequence from the structural PDB file (may contain gaps)
    flat_seq: outputted from SequenceMappingContainer(chain), chain with all gaps removed, but original indices stored
    kingdom: one of eukaryota, bacteria, archaea, or None
    '''
    
    name_arr            = seq.name.split('_')
    parent              = name_arr[1]
    chain_id            = name_arr[2]
    mapped_conserved    = []
    landmark_number     = 0
    
    for i, pos in enumerate(conserved_positions):
        
        coords = locate_residues(pos, polymer, chain_id, parent, chain, flat_seq, kingdom)
        
        if coords is None: continue
        xyz = [coords['x'], coords['y'], coords['z']]
        
        dist = find_closest_point(xyz, parent) 
        print(f"{pos} = {dist}")
        
        if dist <= threshold:
            
            if kingdom is not None:
                tag = kingdom_tag(kingdom)
                pos.name = f'{polymer}-{landmark_number}-{tag}'
            else:
                pos.name = f'{polymer}-{landmark_number}'
                
            landmark_number += 1
            mapped_conserved.append(pos)
        
    return mapped_conserved    
    
# Map conserved residue locations to orignal sequence positions
def map_to_original(sequence: Seq, position: int) -> int:
    ungapped_position = 0
    
    # iterate through the aligned sequence
    for i, residue in enumerate(sequence):
        if residue != "-":
            if i == position:
                return ungapped_position
            ungapped_position += 1
           
    return None

# Given an MSA column, return the most common element if it is at least as frequent as threshold
def find_conserved(column, threshold):
    counter = Counter(column)
    mode = counter.most_common(1)[0]
    
    if (mode[0] != '-' and mode[1] / len(column) >= threshold):
        return mode[0]
    
    return None
    
def get_rcsb_in_alignment(alignment, rcsb_id):
    for i, seq in enumerate(alignment):
        parent = seq.name.split('_')[1]
        if parent == rcsb_id: return alignment[i]
    return None

def kingdom_tag(kingdom: str):
    if kingdom == "eukaryota":
        return 'e'
    elif kingdom == 'bacteria':
        return 'b'
    elif kingdom == 'archaea':
        return 'a'
    else:
        return None