from ribosome_exit_tunnel.tunnel_coordinates import *
from ribosome_exit_tunnel.landmark import *
from Bio.Seq import Seq
from collections import Counter

# Given a list of potential landmarks in the form: [(31, 'F'), (43, 'V'), ..]
def cherry_pick(polymer: str, 
                seq: Seq, 
                conserved_positions: list[Landmark], 
                threshold: float,
                rna: bool = False,
                ) -> list[Landmark]:
    
    name_arr = seq.name.split('_')
    parent = name_arr[1]
    chain = name_arr[2]
    mapped_conserved = []
    landmark_number = 0
    
    for i, pos in enumerate(conserved_positions):
        val = map_to_original(seq, pos.position)
        if val is None: continue
        
        landmark = Landmark(val, pos.residue)
        coords = landmark.get_landmark_coordinates(chain, parent, rna)
        if coords is None: continue
        xyz = [coords['x'], coords['y'], coords['z']]
        
        dist = find_closest_point(xyz, parent) 
        print(f"{pos} = {dist}")
        
        if dist <= threshold:
            pos.name = f'{polymer}-{landmark_number}'
            landmark_number += 1
            mapped_conserved.append(pos)
        
    return mapped_conserved

# Given an MSA column, return the most common element if it is at least as frequent as threshold
def find_conserved(column, threshold):
    counter = Counter(column)
    mode = counter.most_common(1)[0]
    
    if (mode[0] != '-' and mode[1] / len(column) >= threshold):
        return mode[0]
    
    return None

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