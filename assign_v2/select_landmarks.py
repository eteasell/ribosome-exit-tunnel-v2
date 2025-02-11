from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB.MMCIFParser import MMCIFParser
from pathlib import Path
from library.data_access import *
from library.protocol_setup import select_landmarks
from library.types import Landmark
from library.types import PROTOTYPES
from library.taxonomy import *
from library.sequence import *
from library.locate_residues import *
from assign_v2 import ASSIGN_DIR, FASTA_DIR
import csv

def main(conservation_threshold: float, distance_threshold: float):
    
    conserved = select_landmarks(conservation_threshold, distance_threshold, kingdom=None, taxv2=True)
    
    path = ASSIGN_DIR / f"close_conserved_{conservation_threshold}_{distance_threshold}.csv"
    with open(path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=["chain", "residue", "position"])

        if file.tell() == 0:
            writer.writeheader()

        writer.writerows([{'chain': obj.name, 'residue':obj.residue, 'position': obj.position} for i, obj in enumerate(conserved)])
    
    
if __name__ == "__main__":
    CONSERVATION = 0.8
    DISTANCE = 7.5
    
    main(CONSERVATION, DISTANCE)