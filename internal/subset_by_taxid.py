from library.data_access import get_taxid_from_profile
from internal import ASSIGN_DIR, FASTA_DIR
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from library.taxonomy import *

class RcsbTaxIdConverter:
    
    def __init__(self):
        self.rcsb_to_tax_map: dict[str, int] = {}
        self.tax_to_rcsb_map: dict[int, list[str]] = {}
    
    def get_taxid(self, rcsb: str, rank: PhylogenyRank | None = None) -> int:
        if self.rcsb_to_tax_map.get(rcsb) is None:
            tax_id = get_taxid_from_profile(rcsb)

            if rank is not None:
                tax_id = TaxId.coerce_to_rank(tax_id, rank)
            
            self.rcsb_to_tax_map[rcsb] = tax_id
            if tax_id not in self.tax_to_rcsb_map:
                self.tax_to_rcsb_map[tax_id] = []
            self.tax_to_rcsb_map[tax_id].append(rcsb)
        return self.rcsb_to_tax_map[rcsb]

    def get_rcsbs(self, taxid: int) -> list[str] | None:
        return self.tax_to_rcsb_map.get(taxid)


def read_fasta_to_tax_dicts(file: str):
    converter = RcsbTaxIdConverter()
    
    path = FASTA_DIR / file
    with open(path, "r") as fasta_in:
        records: list[SeqRecord] = [*SeqIO.parse(fasta_in, "fasta")]
        
    for record in records:
        rcsb = record.id.split("_")[1]
        converter.get_taxid(rcsb, "phylum")
        
    return converter


def choose_rcsbs_from_tax_list(taxids: list[int], converter: RcsbTaxIdConverter, rank: PhylogenyRank | None = None):
    selected_rcsbs = set()
    for id in taxids:
        id = TaxId.coerce_to_rank(id, rank)
        rcsbs = converter.get_rcsbs(id)
        selected_rcsbs.add(rcsbs.pop())
    return selected_rcsbs

def read_tax_file(file: str):
    path = ASSIGN_DIR / file
    with open(path, "r") as tax_file:
        return [int(line.strip()) for line in tax_file]


if __name__ == "__main__":
    tax_file = "selected_taxids.txt"
    fasta_file = "sequences_uL4.fasta"
    
    taxids = read_tax_file(tax_file)
    converter = read_fasta_to_tax_dicts(fasta_file)
    selected = choose_rcsbs_from_tax_list(taxids, converter, "phylum")
    print(len(selected))
    
    outfile = ASSIGN_DIR / "selected_rcsb_ids.txt"
    with open(outfile, "w") as file:
        for id in selected:
            file.write(str(id) + "\n")
    