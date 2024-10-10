'''
The following is a script to be run in PyMol to extract and save the tunnel coordinates after running MOLE 2.0
'''
from pymol import cmd

# Get the coordinates of the tunnel (assuming 'Tunnel5' is the object or selection name)
xyz = cmd.get_coords('Tunnel5', 1)

# Initialize an empty list for the radii
r = []

# Iterate over the atoms in 'Tunnel' and append the van der Waals radii to the list
cmd.iterate('Tunnel5', 'r.append(vdw)', space=locals())

# Save the coordinates and radius data to text files
import numpy as np
np.savetxt('../data/tunnel/tunnel_coordinates_4ug0.txt', xyz, fmt='%.2f')
np.savetxt('../data/tunnel/tunnel_radius_4ug0.txt', r, fmt='%.2f')