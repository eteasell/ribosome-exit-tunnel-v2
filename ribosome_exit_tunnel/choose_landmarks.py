from ribosome_exit_tunnel.get_coorinates import *
from ribosome_exit_tunnel.tunnel_coordinates import *

# Given a list of potential landmarks in the form: [(31, 'F'), (43, 'V'), ..]
def cherry_pick(seq, conserved_positions, threshold):
    name_arr = seq.name.split('_')
    parent = name_arr[1]
    chain = name_arr[2]
    mapped_conserved = []
    
    for i, pos in enumerate(conserved_positions):
        val = map_to_orignal(seq, pos[0])
        
        coords = get_landmark_coordinates((f'uL4-{i}', pos[1], val+1), chain, parent)
        xyz = [coords['x'], coords['y'], coords['z']]
        
        dist = find_closest_point(xyz, '4ug0') 
        print(f"{pos} = {dist}")
        
        if dist <= threshold:
            mapped_conserved.append(pos)
        
    return mapped_conserved