from protocol.protocol_setup import *
import csv
import time
import argparse

'''
This is the script to run to get landmarks for a given rcsb_id.

To run this script, navigate to 'ribosome-exit-tunnel, and execute the following prompt:

    python3 -m protocol.main [rcsb_id]
    
where [rcsb_id] is the ID of the ribosome of interest.

Note that for debugging, the rcsb_id parameter is set in the launch.json file.

To rerun align or landmark selection, change the boolean flags to true. 
NOTE: Rerunning the alignment will take a long time (~ 10 minutes). The only reason for rerun alignment is if you have
changed the input fasta files to include new specimens that were not previously present. Otherwise, there should be no need.

To see changes in CONSERVATION and DISTANCE parameters, you must change 'reselect_landmarks' to True.
'''
    
reselect_landmarks = False

CONSERVATION= 0.9
DISTANCE = 7.5

def main(list_rcsb_id):  

    t1 = time.time()
    
    processed_ids, changed_files = process_list(list_rcsb_id)
    
    align(changed_files)
    
    for rcsb_id in processed_ids:
        
        t2 = time.time()
    
        kingdom = find_kingdom(rcsb_id)
        
        if kingdom is not None and PROTOTYPES[kingdom] is not None:
            
            if reselect_landmarks:
                conserved = select_landmarks(kingdom, CONSERVATION, DISTANCE)
            
                with open(f"data/output/conserved/alignment_close_conserved_{kingdom}.csv", mode='w', newline='') as file:
                    writer = csv.DictWriter(file, fieldnames=["chain", "residue", "position"])
        
                    if file.tell() == 0:
                        writer.writeheader()
        
                    writer.writerows([{'chain': obj.name, 'residue':obj.residue, 'position': obj.position} for i, obj in enumerate(conserved)])

            
            polymers = PROTOTYPES[kingdom].keys()

            conserved_positions = []
            with open(f"data/output/conserved/alignment_close_conserved_{kingdom}.csv", mode='r', newline='') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    landmark = Landmark(int(row['position']), row['residue'], row['chain'])
                    conserved_positions.append(landmark)
            
            rows = locate_landmarks(rcsb_id, kingdom, conserved_positions, polymers)
            
            if rows is None: continue
                    
            with open(f"data/output/landmarks/{kingdom}/landmarks_{rcsb_id}.csv", mode='w', newline='') as file:
                writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "residue", "position", "x", "y", "z"])
                if file.tell() == 0:
                    writer.writeheader()
                writer.writerows(rows)
        else:
            if kingdom is None: 
                print(f"Cannot access data for specimen {rcsb_id}")
            else:
                print(f"System not yet calibrated to {kingdom}.")
                
        print(f"Completed {rcsb_id}. Duration: {time.time() - t2}")
        
    duration = time.time() - t1
    print(f"Duration: {duration}")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process multiple ribosome instances.")
    parser.add_argument("rcsb_ids", type=str, nargs="+", help="The RCSB IDs to process")
    args = parser.parse_args()
    
    main(args.rcsb_ids)