import numpy as np
from protocol.data_access import find_kingdom
import csv

cache = {}

# Instance is a RCSB ID (for example 4ug0)
def get_tunnel_coordinates_eukaryota(instance: str) -> dict[int,list[float]]:
    
    if instance not in cache:
        xyz = open(f"data/tunnel/tunnel_coordinates_{instance}.txt", mode='r')
        xyz_lines = xyz.readlines()
        xyz.close()
    
        r = open(f"data/tunnel/tunnel_radius_{instance}.txt", mode='r')
        r_lines = r.readlines()
        r.close()
    
        coords = {}
    
        for i, line in enumerate(xyz_lines):
            if (i >= len(r_lines)): break
        
            content = line.split(" ")
            content.append(r_lines[i])
        
            cleaned = []
            for str in content:
                str.strip()
                try:
                    val = float(str)
                    cleaned.append(val)
                except:
                    None
        
            coords[i] = cleaned
        cache[instance] = coords
        
    # Each value in coords is of the form [x, y, z, r]
    return cache[instance]


def get_tunnel_coordinates(instance: str) -> dict[int,list[float]]:
    
    if instance not in cache:
        
        coords = {}
        with open(f"data/tunnel/tunnel_{instance}.csv", mode='r', newline='') as file:
            reader = csv.DictReader(file)
            i = 0
            for row in reader:
                content = [float(row['X']), float(row['Y']), float(row['Z']), float(row['FreeRadius'])]
                coords[i] = content
                i+=1                
        cache[instance] = coords
        
    # Each value in coords is of the form [x, y, z, r]
    return cache[instance]


# p is a list [x,y,z]
# instance is RCSB_ID code
def find_closest_point(p, instance):
    
    if find_kingdom(instance) == "eukaryota":
        coords = get_tunnel_coordinates_eukaryota(instance)
    else:
        coords = get_tunnel_coordinates(instance)
        
    dist = np.inf
    r = 0
    p = np.array(p)
    
    for coord in coords.values():
        xyz = np.array(coord[0:3])
        euc_dist = np.sqrt(np.sum(np.square(xyz - p))) - coord[3]
        if euc_dist < dist:
            dist = euc_dist
    
    return dist