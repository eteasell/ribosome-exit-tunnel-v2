from Bio.Seq import Seq
import subprocess
from Bio import AlignIO
import requests
from data_access import *
import pymol

# Call RibosomeXYZ API to get all uL4 proteins
seqs = get_proteins("uL4")
# Saves to "input_sequences_uL4.fasta"
save_fasta(seqs, "uL4")
          
mafft_command = f"mafft --auto --leavegappyregion {"data/input_sequences_uL4.fasta"} > {"data/aligned_sequences_uL4.fasta"}"
subprocess.call(mafft_command, shell=True)

alignment = AlignIO.read("data/aligned_sequences_uL4.fasta", "fasta")

num_sequences = len(alignment)

conserved_positions = []

for i in range(alignment.get_alignment_length()):
    column = alignment[:,i]
    
    if (len(set(column)) == 1):
        conserved_positions.append((i,column[0]))
  
      
print(conserved_positions)

# Find conserved for first example G73

seq1 = alignment[0]
mapped_conserved = []

for pos in conserved_positions:
    print(pos[0])
    val = map_to_orignal(seq1, pos[0])
    mapped_conserved.append((val, pos[1]))
    
print(mapped_conserved)

'''
For first example only
using pymol

G55 at (296,274,335)
R68 at (300,286,319)
G81 at (289,273,325)
'''