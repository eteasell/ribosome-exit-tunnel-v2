from protocol.data_access import *
from protocol.protocol_setup import *

chain = "23S"
kingdom = 'archaea'

get_mmcif('4V6U')

seqs = read_polymers_file_filter_by_kingdom(chain, kingdom)

save_fasta_by_kingdom(seqs, chain, kingdom)

input = f"data/fasta/sequences_{kingdom}_{chain}.fasta"
output = f"data/output/fasta/aligned_sequences_{kingdom}_{chain}.fasta"
mafft_command = f"mafft --auto --leavegappyregion {input} > {output}"
subprocess.call(mafft_command, shell=True)