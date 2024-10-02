from pymol import cmd
import numpy as np
import csv
from Bio.SeqUtils import seq3

# Example parameters:
# landmark = ('uL4-1', 'G', '56')
# chain = 'C'
# parent = '6vwl'
def get_landmark_coordinates(landmark, chain, parent):
    
    try:
        cmd.load(f'data/mmcif/{parent}.cif')
    except:
        return None
    
    select = f"resi {landmark[2]} and chain {chain}"
    
    atom_coords = []
    cmd.iterate_state(1, select, 'atom_coords.append((resn, x, y, z))', space={'atom_coords': atom_coords})
    cmd.delete(f'{parent}')
    
    if (atom_coords[0][0] != seq3(landmark[1]).upper()):
        return None
    
    vec = np.zeros(3)
    for coord in atom_coords:
        tmp_arr = np.array([coord[1], coord[2], coord[3]])
        vec += tmp_arr

    vec = vec / len(atom_coords)
    vec = vec.astype(np.int32)

    return {"parent_id": parent, "landmark": landmark[0], "x": vec[0], "y": vec[1], "z": vec[2]}