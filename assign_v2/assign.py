from library.protocol_setup import *
import csv
import time
import argparse
from assign_v2 import LANDMARKS_DIR

'''
This is a new version of the 'main.py' and 'main_universal.py' files. This version of the landmark assignment protocol
addresses the issue of biased phylogeny being used in conservation analysis. There is only one version of the
taxonomy-adusted protocol (i.e. no seperation by kingdom or universal options). For simplicity, this script takes only one
structure at a time, which means it will be much slower than the other two scripts for assigning landmarks to a large number
of structures.

To run this script, navigate to 'ribosome-exit-tunnel-v2', and execute the following prompt:

    python -m assign_v2.assign rcsb_id
    
where rcsb_id is the ID of the ribosome of interest.

Note that for debugging, the rcsb_id parameter is set in the launch.json file.

NOTE: This script assumes that the landmarks have already been selected, which is done by assign_v2/select_landmarks.py. 
This script is only for assigning previously chosen landmarks.

'''

def main(rcsb_id):  

    t1 = time.time()
    
    # NOTE: changing the below params only changes which file the landmarks are chosen from. Changing these numbers to ones
    # that have not been run using select_landmarks.py will cause the code to break
    CONS_PARAM = 0.8
    DIST_PARAM = 7.5
    
    polymers = UNIVERSAL_PROTOTYPE.keys()

    conserved_positions = []
    path = ASSIGN_DIR / f"close_conserved_{CONS_PARAM}_{DIST_PARAM}.csv"
    with open(path, mode='r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            landmark = Landmark(int(row['position']), row['residue'], row['chain'])
            conserved_positions.append(landmark)
    
    rows = locate_landmarks(rcsb_id, conserved_positions, polymers, taxv2=True)
                
    with open(LANDMARKS_DIR / f"assign_v2/landmarks_{rcsb_id}.csv", mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "residue", "position", "x", "y", "z"])
        if file.tell() == 0:
            writer.writeheader()
        writer.writerows(rows)
    
    duration = time.time() - t1
    print(f"Duration: {duration}")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assign landmarks to a ribosome structure")
    parser.add_argument("rcsb_id", type=str, help="The RCSB ID to process")

    args = parser.parse_args()
    main(args.rcsb_id)