from ribosome_exit_tunnel.protocol_setup import *
import csv

# This is the script to run to get landmarks for a given rcsb_id
rcsb_id = '4UG0'
rerun_alignment = False
reselect_landmarks = False

CONSERVATION= 0.9
DISTANCE = 7.5

kingdom = find_kingdom(rcsb_id)
if kingdom is not None and PROTOTYPES[kingdom] is not None:
    
    if rerun_alignment:
        align(kingdom)
        
    if reselect_landmarks:
        conserved = select_landmarks(kingdom, CONSERVATION, DISTANCE)
        
        with open(f"data/output/conserved/alignment_conserved_{kingdom}.csv", mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["chain", "residue", "position"])
    
            if file.tell() == 0:
                writer.writeheader()
    
            writer.writerows([{'chain': obj.name, 'residue':obj.residue, 'position': obj.position} for i, obj in enumerate(conserved)])

        
    polymers = PROTOTYPES[kingdom].keys()

    conserved_positions = []
    with open("data/output/conserved/alignment_conserved_eukaryota.csv", mode='r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            landmark = Landmark(int(row['position']), row['residue'], row['chain'])
            conserved_positions.append(landmark)
        
    rows = locate_landmarks(rcsb_id, 'eukaryota', conserved_positions, polymers)
                 
    with open(f"data/output/landmarks/landmarks_{rcsb_id}.csv", mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "residue", "position", "x", "y", "z"])
    
        if file.tell() == 0:
            writer.writeheader()
    
        writer.writerows(rows)
else:
    if kingdom is None: 
        print(f"Cannot access data for specimen {rcsb_id}")
    else:
        print(f"System not yet calibrated to {kingdom}.")