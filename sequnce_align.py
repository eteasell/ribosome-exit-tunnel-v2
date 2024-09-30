from Bio.Seq import Seq
import subprocess
from Bio import AlignIO
import requests
from data_access import *

# Call RibosomeXYZ API to get all uL4 proteins
seqs = get_proteins("uL4")
# Saves to "input_sequences_uL4.fasta"
save_fasta(seqs, "uL4")
          
mafft_command = f"mafft --auto {"data/input_sequences_uL4.fasta"} > {"data/aligned_sequences_uL4.fasta"}"
subprocess.call(mafft_command, shell=True)

alignment = AlignIO.read("data/aligned_sequences_uL4.fasta", "fasta")
print(alignment)