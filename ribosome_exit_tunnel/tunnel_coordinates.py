import numpy as np

# Instance is a RCSB ID (for example 4ug0)
def get_tunnel_coordinates(instance: str) -> dict[int,list[float]]:
    
    if instance not in get_tunnel_coordinates.cache:
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
        get_tunnel_coordinates.cache[instance] = coords
        
    # Each value in coords is of the form [x, y, z, r]
    return get_tunnel_coordinates.cache[instance]

get_tunnel_coordinates.cache = {}

# p is a list [x,y,z]
# instance is RCSB_ID code
def find_closest_point(p, instance):
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