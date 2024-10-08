from Bio import AlignIO
from get_coorinates import *

# Simply a tester to run functions
# This is on the first protein in the alignment 6vwl (e coli)

conserved_positions = [(31, 'F'), (43, 'V'), (52, 'R'), (53, 'Q'), (65, 'T'), 
(68, 'E'), (71, 'G'), (73, 'G'), (84, 'G'), (86, 'G'), 
(89, 'R'), (91, 'G'), (101, 'G'), (102, 'G'), (105, 'F'), 
(126, 'A'), (129, 'S'), (157, 'V'), (169, 'K'), (170, 'T'), 
(171, 'K'), (178, 'K'), (231, 'A'), (233, 'R'), (234, 'N'), 
(238, 'V'), (267, 'A')]

alignment = AlignIO.read("data/aligned_sequences_uL4.fasta", "fasta")
seq = alignment[0]

name_arr = seq.name.split('_')
parent = name_arr[1]
chain = name_arr[2]
mapped_conserved = []
    
for i, pos in enumerate(conserved_positions):
    val = map_to_orignal(seq, pos[0])
    mapped_conserved.append((val, pos[1]))
    
print(mapped_conserved)

'''
[(18, 'F'), (30, 'V'), (39, 'R'), (40, 'Q'), (47, 'T'), (50, 'E'), 
(53, 'G'), (55, 'G'), (63, 'G'), (65, 'G'), (68, 'R'), (70, 'G'), 
(80, 'G'), (81, 'G'), (84, 'F'), (103, 'A'), (106, 'S'), (120, 'V'), 
(129, 'K'), (130, 'T'), (131, 'K'), (138, 'K'), (159, 'A'), (161, 'R'), 
(162, 'N'), (166, 'V'), (191, 'A')]

The idea for the manual cherry-picking is:

1. Use a conserved residue on uL4 to start the algo
2. run the MOLE 2.0 algo to get the tunnel cooridinates
    a) 
3. for each potential landmark, 
    a) find the closest point on the tunnel centerline to the given residue.
    b) find the distance between the residue and the centerline radius at that given point
    c) decide some exlusion threshold for the maximal allowed distance to the tunnel
    
* is it even possible to run the full MOLE algo from the script?
'''

# lets choose G65 as our starter residue

cmd.fetch('6vwl', path= "data")

# Select residues 65 from chain C
cmd.select('my_selection', 'resi 65 and chain C')

# Select atoms within 80 Å of the selection
cmd.select('within_80A', 'br. 6vwl w. 80 of my_selection')

# Create a new object from the selection
cmd.create('6vwl_edited', 'within_80A')

# Save the new object as a .cif file
cmd.save('output_structure.cif', '6vwl_edited')

## Need to locate the PTC

# from riboxyz
DORIS_ET_AL = {
    "SITE_6": "AAGACCC",
    "SITE_8": "GGAUAAC",
    "SITE_9": "GAGCUGGGUUUA",
}

# or instead use the known locations U4452 in euks and U2585 in proks

select = 'resi 2585 and chain 2'

atom_coords = []
cmd.iterate_state(1, select, 'atom_coords.append((chain, resn, x, y, z))', space={'atom_coords': atom_coords})
    
vec = np.zeros(3)
for coord in atom_coords:
    tmp_arr = np.array([coord[2], coord[3], coord[4]])
    vec += tmp_arr

vec = vec / len(atom_coords)
vec = vec.astype(np.int32)

# now vec has the coordinates for the ptc
print(vec)

