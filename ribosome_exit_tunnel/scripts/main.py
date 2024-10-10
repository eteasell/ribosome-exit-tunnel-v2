from Bio.Seq import Seq
import subprocess
from Bio import AlignIO
from ribosome_exit_tunnel.data_access import *
import csv
from ribosome_exit_tunnel.get_coorinates import *
from ribosome_exit_tunnel.choose_landmarks import *
import time

# This script is for uL4 only - for now
t1 = time.time()

conserved_threshold = 0.9
distance_threshold = 10

# Itialize with 4ug0 (prototype)
proto_uL4 = {'parent_id': '4ug0', 
             'auth_asym_id': 'LC', 
             'seq': 'MACARPLISVYSEKGESSGKNVTLPAVFKAPIRPDIVNFVHTNLRKNNRQPYAVSELAGHQTSAESWGTGRAVARIPRVRGGGTHRSGQGAFGNMCRGGRMFAPTKTWRRWHRRVNTTQKRYAICSALAASALPALVMSKGHRIEEVPELPLVVEDKVEGYKKTKEAVLLLKKLKAWNDIKKVYASQRMRAGKGKMRNRRRIQRRGPCIIYNEDNGIIKAFRNIPGITLLNVSKLNILKLAPGGHVGRFCIWTESAFRKLDELYGTWRKAASLKSNYNLPMHKMINTDLSRILKSPEIQRALRAPRKKIHRRVLKKNPLKNLRIMLKLNPYAKTMRRNTILRQARNHKLRVDKAAAAAAALQAKSDEKAAVAGKKPVVGKKGKKAAVGVKKQKKPLVGKKAAATKKPAPEKKPAEKKPTTEEKKPAA'}
seqs =[proto_uL4]
# Call RibosomeXYZ API to get all uL4 proteins
seqs = seqs + read_polymers_file("uL4")
# Saves to "input_sequences_uL4.fasta"
save_fasta(seqs, "uL4")
      
# Perform multiple sequence alignment    
mafft_command = f"mafft --auto --leavegappyregion {"data/output/input_sequences_uL4.fasta"} > {"data/output/aligned_sequences_uL4.fasta"}"
subprocess.call(mafft_command, shell=True)
alignment = AlignIO.read("data/output/aligned_sequences_uL4.fasta", "fasta")

conserved_positions = []

for i in range(alignment.get_alignment_length()):
    column = alignment[:,i]
    most_common = find_conserved(alignment[:,i], conserved_threshold)
    
    if (most_common is not None):
        conserved_positions.append((i,most_common))
        
conserved_positions = cherry_pick(alignment[0], conserved_positions, distance_threshold)

with open("data/output/alignment_conserved.csv", mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["chain", "position"])
    
    if file.tell() == 0:
        writer.writeheader()
    
    writer.writerows([{'chain': 'uL4', 'position': obj} for obj in conserved_positions])

rows = []

for i, seq in enumerate(alignment[0:50]):
    print(i)
    name_arr = seq.name.split('_')
    parent = name_arr[1]
    chain = name_arr[2]
    
    for i, pos in enumerate(conserved_positions):
        val = map_to_orignal(seq, pos[0])
        
        if val is not None:
            coords = get_landmark_coordinates((f'uL4-{i}', pos[1], val+1), chain, parent)
        
        if coords is not None:
            rows.append(coords)
        
with open("data/output/landmarks.csv", mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "residue", "x", "y", "z"])
    
    if file.tell() == 0:
        writer.writeheader()
    
    writer.writerows(rows)
    
duration = time.time() - t1
print("Script Time: " + str(duration))

# Script Time on 100: 561.4538388252258
# drastically slowed down as it went on