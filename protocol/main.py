from protocol.protocol_setup import *
import csv
import time
import argparse

'''
This is the script to run to get landmarks for a given list of rcsb_id's.

To run this script, navigate to 'ribosome-exit-tunnel-v2', and execute the following prompt:

    python -m protocol.main [rcsb_id]
    
where [rcsb_id] is the list of IDs of the ribosomes of interest.

Note that for debugging, the rcsb_id parameter is set in the launch.json file.

To see changes in CONSERVATION and DISTANCE parameters, you must change 'reselect_landmarks' to True.

The FASTA files stored in /data/fasta contain sequences for a large number of ribosome specimens. However, if you are 
running this code on a large number of new ribosomes, you should set 'reselect_landmarks' to True. This is because there is a
chance that conserved residues will hvae changed do to the new sequences in th alignments, so reselection and reassignment
should be done for best accuracy.

NOTE: This code will automatically run the alignements when the input fasta files have been changed; however, the alignments
can also be run with MAFFT online, which will be fasta than through this code.
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
    processed_ids, changed_files = process_list(list_rcsb_id)
    
    print(f"Starting alignment on {len(changed_files)} files.")
    
    align(changed_files)
    
    print(f"Alignments completed in {time.time() - t1} seconds.")
    
    # Keep track of which kingdoms we have already reselected landmarks for in this runtime 
    reselected = []
    
    for rcsb_id in processed_ids:
        
        t2 = time.time()
    
        kingdom = find_kingdom(rcsb_id)
        
        if kingdom is not None and PROTOTYPES[kingdom] is not None:
            
            if reselect_landmarks and kingdom not in reselected:
                conserved = select_landmarks(CONSERVATION, DISTANCE, kingdom)
            
                with open(f"data/output/conserved/alignment_close_conserved_{kingdom}.csv", mode='w', newline='') as file:
                    writer = csv.DictWriter(file, fieldnames=["chain", "residue", "position"])
        
                    if file.tell() == 0:
                        writer.writeheader()
        
                    writer.writerows([{'chain': obj.name, 'residue':obj.residue, 'position': obj.position} for i, obj in enumerate(conserved)])
                
                reselected.append(kingdom)
            
            polymers = PROTOTYPES[kingdom].keys()

            conserved_positions = []
            with open(f"data/output/conserved/alignment_close_conserved_{kingdom}.csv", mode='r', newline='') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    landmark = Landmark(int(row['position']), row['residue'], row['chain'])
                    conserved_positions.append(landmark)
            
            rows = locate_landmarks(rcsb_id, conserved_positions, polymers, kingdom)
            
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
