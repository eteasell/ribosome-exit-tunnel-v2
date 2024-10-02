from Bio.Seq import Seq
import subprocess
from Bio import AlignIO
from data_access import *
import csv
from get_coorinates import *

# This script is for uL4 only
chain = 'C'
# Call RibosomeXYZ API to get all uL4 proteins
seqs = get_proteins("uL4")
# Saves to "input_sequences_uL4.fasta"
save_fasta(seqs, "uL4")
      
# Perform multiple sequence alignment    
mafft_command = f"mafft --auto --leavegappyregion {"data/input_sequences_uL4.fasta"} > {"data/aligned_sequences_uL4.fasta"}"
subprocess.call(mafft_command, shell=True)
alignment = AlignIO.read("data/aligned_sequences_uL4.fasta", "fasta")

num_sequences = len(alignment)

conserved_positions = []

# Collapse each column in the MSA, if length is 1, sequence is fully conserved
# To Do: should maybe lower the conservation criteria to 90% or something
for i in range(alignment.get_alignment_length()):
    column = alignment[:,i]
    if (len(set(column)) == 1):
        conserved_positions.append((i,column[0]))

rows = []

for seq in alignment:
    parent = seq.name.split('_')[1]
    mapped_conserved = []
    
    for i, pos in enumerate(conserved_positions[0:3]):
        val = map_to_orignal(seq, pos[0])
        mapped_conserved.append((val, pos[1]))
        
        coords = get_landmark_coordinates((f'uL4-{i}', pos[1], val+1), chain, parent)
        
        if coords is not None:
            rows.append(coords)
        
with open("data/output/landmarks.csv", mode='w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=["parent_id", "landmark", "x", "y", "z"])
    
    if file.tell() == 0:
        writer.writeheader()
    
    writer.writerows(rows)