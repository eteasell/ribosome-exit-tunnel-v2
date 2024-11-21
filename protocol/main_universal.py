from protocol.protocol_setup import *
import csv
import time
import argparse

'''
This is the script to run to get UNIVERSAL landmarks for a given list of rcsb_id's. This means that the landmarks 
are not selected separately by kingdom.

To run this script, navigate to 'ribosome-exit-tunnel', and execute the following prompt:

    python3 -m protocol.main_universal [rcsb_id]
    
where [rcsb_id] is the list of IDs of the ribosomes of interest.

Note that for debugging, the rcsb_id parameter is set in the launch.json file.

To see changes in CONSERVATION and DISTANCE parameters, you must change 'reselect_landmarks' to True.

'''
    
reselect_landmarks = False

CONSERVATION= 0.9
DISTANCE = 7.5

def main(list_rcsb_id):  

    t1 = time.time()
    
    '''
    Processed id's are the filtered list of id's for which the structure and sequence files were successfully accessed
    Changed files is the list of fasta files which have had additions to include the given ids. 
    '''
    processed_ids, changed_files = process_list(list_rcsb_id, universal=True)
    
    print(f"Starting alignment on {len(changed_files)} files.")
    
    align(changed_files)
    
    print(f"Alignments completed in {time.time() - t1} seconds.")
    
    if len(changed_files) > 0 or reselect_landmarks:
        conserved = select_landmarks(CONSERVATION, DISTANCE)
    
        with open(f"data/output/conserved/alignment_close_conserved.csv", mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["chain", "residue", "position"])

            if file.tell() == 0:
                writer.writeheader()

            writer.writerows([{'chain': obj.name, 'residue':obj.residue, 'position': obj.position} for i, obj in enumerate(conserved)])
        
    for rcsb_id in processed_ids:
        
        t2 = time.time()
            
        polymers = UNIVERSAL_PROTOTYPE.keys()

        conserved_positions = []
        with open(f"data/output/conserved/alignment_close_conserved.csv", mode='r', newline='') as file:
            reader = csv.DictReader(file)
            for row in reader:
                landmark = Landmark(int(row['position']), row['residue'], row['chain'])
                conserved_positions.append(landmark)
        
        rows = locate_landmarks(rcsb_id, conserved_positions, polymers)
        
        if rows is None: continue
                
        with open(f"data/output/landmarks/universal/landmarks_{rcsb_id}.csv", mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "residue", "position", "x", "y", "z"])
            if file.tell() == 0:
                writer.writeheader()
            writer.writerows(rows)
                
        print(f"Completed {rcsb_id}. Duration: {time.time() - t2}")
        
    duration = time.time() - t1
    print(f"Duration: {duration}")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process multiple ribosome instances.")
    parser.add_argument("rcsb_ids", type=str, nargs="+", help="The RCSB IDs to process")
    args = parser.parse_args()
    
    main(args.rcsb_ids)