from pymol import cmd
import numpy as np
from Bio.SeqUtils import seq3
from data_access import *
from pathlib import Path

# Example parameters:
# landmark = ('uL4-1', 'G', '56')
# chain = 'C'
# parent = '6vwl'
def get_landmark_coordinates(landmark, chain, parent):
    
    try:
        parent_chain = f'{parent}_{chain}'
        file = f'data/mmcif/{parent_chain}.cif'
        
        if Path(file).is_file() is False:
            get_mmcif(parent, chain)
               
        if parent_chain not in cmd.get_names():
            cmd.load(f'data/mmcif/{parent_chain}.cif', object=f'{parent_chain}')
        
    except:
        return None
    
    select = f"resi {landmark[2]} and chain {chain}"
    
    atom_coords = []
    cmd.iterate_state(1, select, 'atom_coords.append((chain, resn, x, y, z))', space={'atom_coords': atom_coords})
    
    if (atom_coords[0][1] != seq3(landmark[1]).upper()):
        return None
    
    vec = np.zeros(3)
    for coord in atom_coords:
        tmp_arr = np.array([coord[2], coord[3], coord[4]])
        vec += tmp_arr

    vec = vec / len(atom_coords)
    vec = vec.astype(np.int32)

    return {"parent_id": parent, "landmark": landmark[0], "residue": f"{landmark[1]}{landmark[2]}", "x": vec[0], "y": vec[1], "z": vec[2]}

# Map conserved residue locations to orignal sequence positions
def map_to_orignal(sequence, position):
    ungapped_position = 0
    
    for i, residue in enumerate(sequence):
        if residue != "-":
            if i == position:
                return ungapped_position
            ungapped_position += 1
            
    return None