from pymol import cmd
import numpy as np
from Bio.SeqUtils import seq3
from protocol.data_access import *
from pathlib import Path

class Landmark:
    def __init__(self, position: int, residue: str, name: str = None):
        self.name = name
        self.position = position
        self.residue = residue
        
    def __str__(self):
        return f"{self.residue}{self.position}"

    def get_landmark_coordinates(self, chain: str, parent: str, rna: bool = False) -> dict:
        '''
        Note that this method is not used in the main code, as it was replaced by "locate_residues", which is more robust to
        inconsistent PDB residue numbering.
        '''
    
        try:
            file = f'data/mmcif/{parent}.cif'
        
            if Path(file).is_file() is False:
                get_mmcif(parent)
               
            print(cmd.get_names())
            if f'{parent}_{chain}' not in cmd.get_names():
                cmd.load(f'data/mmcif/{parent}.cif', object=f'{parent}_{chain}')
                cmd.remove(f'not chain {chain}')
            cmd.select(f'{parent}_{chain}')
        
        except:
            return None
    
        select = f"resi {self.position + 1}"
    
        atom_coords = []
        cmd.iterate_state(1, select, 'atom_coords.append((chain, resn, x, y, z))', space={'atom_coords': atom_coords})
    
        if (len(atom_coords) == 0): 
            return None
        elif (rna and atom_coords[0][1] != self.residue.upper()): 
            return None
        elif (not rna and atom_coords[0][1] != seq3(self.residue).upper()):
            return None
    
        vec = np.zeros(3)
        for coord in atom_coords:
            tmp_arr = np.array([coord[2], coord[3], coord[4]])
            vec += tmp_arr

        vec = vec / len(atom_coords)
        vec = vec.astype(np.int32)
        
        return {"parent_id": parent, 
                "landmark": self, 
                "x": vec[0], "y": vec[1], "z": vec[2]}