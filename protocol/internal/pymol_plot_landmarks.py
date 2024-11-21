# This is a script to plot the landmarks into Pymol

import csv 
from pymol import cmd

kingdom = "bacteria"
rcsb_id = '6YST'
coordinates = []

path = f'/Users/ellateasell/ISCI448B/code/repos/-v2/data/output/landmarks/{kingdom}/landmarks_{rcsb_id}.csv'

# get list of coordinates
with open(path, mode='r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row['parent_id'] == rcsb_id:
                coordinates.append((int(row['x']),int(row['y']),int(row['z']), row['landmark']))

# Create a new object to store the points
obj_name = 'points'
obj_name2 = 'labels'

# Loop over the coordinates and create pseudo-atoms at each point
for i, (x, y, z, landmark) in enumerate(coordinates):
    cmd.pseudoatom(obj_name, pos=[x, y, z], name=landmark)
    cmd.pseudoatom(obj_name2, pos=[x, y, z], name=landmark)
    cmd.label(obj_name2, "name")


# Show the points as spheres
cmd.show('spheres', obj_name)
cmd.set('sphere_scale', 1.25)  # Adjust the size of the points

# Optionally, you can change the color
cmd.color('red', obj_name)

# Refresh the display
cmd.zoom(obj_name)

