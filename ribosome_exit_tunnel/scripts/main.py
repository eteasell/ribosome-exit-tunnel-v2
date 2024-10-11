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
        
conserved_positions = cherry_pick('uL4' ,alignment[0], conserved_positions, distance_threshold)

with open("data/output/alignment_conserved.csv", mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["chain", "position"])
    
    if file.tell() == 0:
        writer.writeheader()
    
    writer.writerows([{'chain': 'uL4', 'position': obj} for obj in conserved_positions])

rows = []
parents = []

for i, seq in enumerate(alignment[0:50]):
    print(i)
    name_arr = seq.name.split('_')
    parent = name_arr[1]
    chain = name_arr[2]
    
    parents.append(parent)
    
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
print("Script Tim Part 1: " + str(duration))

# Script Time under 1 min on 50 ribosomes, with pymol call refactored

# Now adding uL22...

proto_uL22 = {'parent_id': '4ug0', 
             'auth_asym_id': 'LP', 
             'seq': 'MVRYSLDPENPTKSCKSRGSNLRVHFKNTRETAQAIKGMHIRKATKYLKDVTLQKQCVPFRRYNGGVGRCAQAKQWGWTQGRWPKKSAEFLLHMLKNAESNAELKGLDVDSLVIEHIQVNKAPKMRRRTYRAHGRINPYMSSPCHIEMILTEKEQIVPKPEEEVAQKKKISQKKLKKQKLMARE'}

seqs =[proto_uL22]

polymers = read_polymers_file("uL22")
polymers_to_include = list(filter(lambda x: x['parent_id'] in parents, polymers))

seqs = seqs + polymers_to_include
# Saves to "input_sequences_uL4.fasta"
save_fasta(seqs, "uL22")

# Perform multiple sequence alignment    
mafft_command = f"mafft --auto --leavegappyregion {"data/output/input_sequences_uL22.fasta"} > {"data/output/aligned_sequences_uL22.fasta"}"
subprocess.call(mafft_command, shell=True)
alignment = AlignIO.read("data/output/aligned_sequences_uL22.fasta", "fasta")

conserved_positions = []

for i in range(alignment.get_alignment_length()):
    column = alignment[:,i]
    most_common = find_conserved(alignment[:,i], conserved_threshold)
    
    if (most_common is not None):
        conserved_positions.append((i,most_common))
    
        
conserved_positions = cherry_pick('uL22', alignment[0], conserved_positions, distance_threshold)

with open("data/output/alignment_conserved.csv", mode='a', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["chain", "position"])
    
    if file.tell() == 0:
        writer.writeheader()
    
    writer.writerows([{'chain': 'uL22', 'position': obj} for obj in conserved_positions])
    
rows = []

for i, seq in enumerate(alignment):
    print(i)
    name_arr = seq.name.split('_')
    parent = name_arr[1]
    chain = name_arr[2]
    
    parents.append(parent)
    
    for i, pos in enumerate(conserved_positions):
        val = map_to_orignal(seq, pos[0])
        
        if val is not None:
            coords = get_landmark_coordinates((f'uL22-{i}', pos[1], val+1), chain, parent)
        
        if coords is not None:
            rows.append(coords)
            
with open("data/output/landmarks.csv", mode='a', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "residue", "x", "y", "z"])
    
    if file.tell() == 0:
        writer.writeheader()
    
    writer.writerows(rows)
    
duration = time.time() - t1
print("Script Time: " + str(duration))