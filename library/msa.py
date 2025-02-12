from pprint import pprint
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import typing
from functools import reduce
from typing import Callable
from library.taxonomy import  TaxId, ncbi
from typing import Dict, List, Set
from collections import defaultdict
from ete3 import NCBITaxa
import random
from rich.console import Console
from rich.tree import Tree as RichTree
from rich import print as rprint
from rich.panel import Panel
from rich.table import Table
from library import ASSIGN_DIR

'''
Code in this class is taken from 'riboxyz':
A. Kushner, A.S. Petrov, K. Dao Duc, (2022) “RiboXYZ: A comprehensive database for ribosome structures”,  Nucleic Acids Research, gkac939
(https://github.com/rtviii/riboxyz/blob/b50cb8426677c884580568e2889033191343b60a/ribctl/lib/libmsa.py)
'''

class Fasta:
    records: list[SeqRecord]
    def __init__( self, file: str | None = None, records: list[SeqRecord] | None = None) -> None:

        if file is None:
            path = None
        else:
            path = ASSIGN_DIR / file
            
        if path is not None:
            try:
                with open(path, "r") as fasta_in:
                    self.records: list[SeqRecord] = [*SeqIO.parse(fasta_in, "fasta")]
            except FileNotFoundError:
                print(f"File not found: {path}")

            except Exception as e:
                print(f"An error occurred: {str(e)}")
                exit(-1)
                
        elif records is not None:
            self.records = records
                            
        
    def _yield_subset(self, predicate: Callable[[SeqRecord], bool]) -> list[SeqRecord]:
        return [*filter(predicate, self.records)]

    def pick_descendants_of_taxid(self, taxid: int) -> list[SeqRecord]:
        """Given a taxid, return the subset of records that are descendants of that taxid"""
        return self._yield_subset( lambda record: TaxId.is_descendant_of(taxid, int(record.id)))

    @staticmethod
    def write_fasta(seqrecords: list[SeqRecord], outfile: str):
        with open(outfile, "w") as fasta_out:
            SeqIO.write(seqrecords, fasta_out, "fasta")
            
    @staticmethod
    def write_record_ids(seqids: set, outfile: str):
        with open(outfile, "w") as file:
            for id in seqids:
                file.write(str(id) + "\n")

    @staticmethod
    def fasta_display_species(taxids: list[int]):
        taxids = TaxId.coerce_all_to_rank(taxids, "species")
        tree   = ncbi.get_topology(taxids)

        for node in tree.traverse():
            taxid           = int(node.name)
            scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
            node.name       = scientific_name

    def pick_taxids(self, _taxids_int: list[int]) -> list[SeqRecord]:
         
        taxids = list(map(lambda x: str(x), _taxids_int))
        for taxid in set(taxids):
            if taxid not in self.all_taxids(return_as='str'):
                raise Exception( f"Taxid {taxid} not found in records. Violated assumption. Did you alter the fast archives recently? There might be a cached HMM Scanner with taxid that is no longer present. " )

        def filter_duplicates(acc, record):
            if record.id not in acc["seen"]:
                acc["seen"].add(record.id)
                acc["result"].append(record)
            return acc

        initial_accumulator = {"seen": set(), "result": []}
        filtered_records = reduce(
            filter_duplicates,
            filter(lambda record: record.id in taxids, self.records),
            initial_accumulator,
        )

        return list(filtered_records["result"])

    def all_taxids(self, return_as:typing.Literal["int", "str"]="str") -> list[int] | list[str]:
        """Given a fasta file, return all the taxids present in it
        With the assumption that the tax id is the id of each seq record."""
        taxids = []
        for record in self.records:
            taxids = [*taxids, record.id]
        if return_as == "str":
            return taxids
        elif return_as == "int":
            return list(map(lambda _: int(_), taxids))
        else:
            raise Exception("Invalid type passed to all_taxids")


class FastaBalancer:
    def __init__(self, fasta: Fasta):
        """Initialize with a Fasta instance containing sequences with taxids as headers."""
        self.fasta = fasta
        self.ncbi = NCBITaxa()
        self.sequences = self._load_sequences()
        self.console = Console()

    def _load_sequences(self) -> Dict[int, str]:
        """Load sequences from Fasta instance."""
        return {int(record.id): str(record.seq) for record in self.fasta.records}

    def get_taxonomic_hierarchy(self) -> Dict[str, Dict[str, List[int]]]:
        """Build complete taxonomic hierarchy of sequences."""
        hierarchy = defaultdict(lambda: defaultdict(list))

        for taxid in self.sequences.keys():
            try:
                kingdom = TaxId.superkingdom(taxid)
                lineage = TaxId.get_lineage(taxid)
                phylum_id = next(
                    (tid for tid in lineage if TaxId.rank(tid) == "phylum"), None
                )

                if phylum_id:
                    phylum_name = TaxId.get_name(phylum_id)
                    hierarchy[kingdom][phylum_name].append(taxid)
            except Exception as e:
                print(f"Warning: Could not process taxid {taxid}: {str(e)}")

        return hierarchy

    def balance_dataset_sparse(self, target_size: int) -> Fasta:
        """
        Balance dataset by maximizing phylogenetic diversity for a given target size.
        Returns a new Fasta instance with only selected sequences.

        Args:
            target_size: Desired total number of sequences
        """
        original_taxids = set(self.sequences.keys())
        hierarchy = self.get_taxonomic_hierarchy()

        self.console.print(
            Panel(
                f"[bold]Starting sparse dataset balancing[/bold]\nTarget total sequences: {target_size}"
            )
        )

        selected_taxids = set()

        total_kingdoms = len(hierarchy)
        total_phyla = sum(len(phyla) for phyla in hierarchy.values())

        if target_size <= total_kingdoms:
            for kingdom in list(hierarchy.keys())[:target_size]:
                phylum_sizes = [
                    (p, len(taxa)) for p, taxa in hierarchy[kingdom].items()
                ]
                largest_phylum = max(phylum_sizes, key=lambda x: x[1])[0]
                selected = random.sample(hierarchy[kingdom][largest_phylum], 1)
                selected_taxids.update(selected)
        else:
            seqs_per_kingdom = target_size // total_kingdoms
            remaining = target_size % total_kingdoms

            for kingdom, phyla in hierarchy.items():
                kingdom_target = seqs_per_kingdom + (1 if remaining > 0 else 0)
                remaining = max(0, remaining - 1)

                if kingdom_target == 0:
                    continue

                n_phyla = len(phyla)
                seqs_per_phylum = max(1, kingdom_target // n_phyla)
                phylum_remainder = kingdom_target % n_phyla

                self.console.print(f"\n[bold]{kingdom}[/bold]:")
                self.console.print(
                    f"Found {n_phyla} phyla, allocating {kingdom_target} sequences"
                )

                for phylum_name, phylum_taxids in phyla.items():
                    phylum_target = min(
                        seqs_per_phylum + (1 if phylum_remainder > 0 else 0),
                        len(phylum_taxids),
                    )
                    phylum_remainder = max(0, phylum_remainder - 1)

                    selected = random.sample(phylum_taxids, phylum_target)
                    selected_taxids.update(selected)
                    self.console.print(
                        f"  • {phylum_name}: selected {len(selected)}/{len(phylum_taxids)} sequences"
                    )

        self.print_filtering_summary(original_taxids, selected_taxids)

        # Create new Fasta instance with only selected sequences

        selected_records = [record for record in self.fasta.records 
                           if int(record.id) in selected_taxids]
        Fasta.write_record_ids(selected_taxids, ASSIGN_DIR / "selected_taxids")
        return Fasta(records=selected_records)

    # [Previous visualization methods remain unchanged]
    def print_taxonomic_tree(self, taxids: List[int], title: str):
        """Print a hierarchical view of the taxonomic distribution."""
        rich_tree = RichTree(f"[bold blue]{title}[/bold blue]")

        kingdom_groups = defaultdict(list)
        for taxid in taxids:
            try:
                kingdom = TaxId.superkingdom(taxid)
                kingdom_groups[kingdom].append(taxid)
            except Exception:
                continue

        for kingdom, kingdom_taxids in kingdom_groups.items():
            kingdom_branch = rich_tree.add(
                f"[bold green]{kingdom}[/bold green] ({len(kingdom_taxids)} sequences)"
            )

            phylum_groups = defaultdict(list)
            for taxid in kingdom_taxids:
                try:
                    lineage = TaxId.get_lineage(taxid)
                    phylum_id = next(
                        (tid for tid in lineage if TaxId.rank(tid) == "phylum"), None
                    )
                    if phylum_id:
                        phylum_groups[TaxId.get_name(phylum_id)].append(taxid)
                except Exception:
                    continue

            for phylum, phylum_taxids in phylum_groups.items():
                phylum_branch = kingdom_branch.add(
                    f"[yellow]{phylum}[/yellow] ({len(phylum_taxids)} sequences)"
                )

        self.console.print(rich_tree)

    def print_filtering_summary(
        self, original_taxids: Set[int], selected_taxids: Set[int]
    ):
        """Print detailed summary of filtering process."""
        filtered_taxids = original_taxids - selected_taxids

        table = Table(title="Filtering Summary")
        table.add_column("Category", style="cyan")
        table.add_column("Count", justify="right", style="green")
        table.add_column("Percentage", justify="right", style="yellow")

        total = len(original_taxids)
        selected = len(selected_taxids)
        filtered = len(filtered_taxids)
        print(">>>>>>>>>> SEEING {} SEQS".format(total))
        table.add_row("Original sequences", str(total), "100%")
        table.add_row(
            "Selected sequences", str(selected), f"{(selected/total)*100:.1f}%"
        )
        table.add_row(
            "Filtered sequences", str(filtered), f"{(filtered/total)*100:.1f}%"
        )

        self.console.print("\n")
        self.console.print(table)

        self.console.print("\n[bold]Taxonomic Distribution Before Filtering:[/bold]")
        self.print_taxonomic_tree(list(original_taxids), "Original Dataset")

        self.console.print("\n[bold]Taxonomic Distribution After Filtering:[/bold]")
        self.print_taxonomic_tree(list(selected_taxids), "Balanced Dataset")

if __name__ == "__main__":
    file = ASSIGN_DIR / "renamed.fasta"
    
    fasta          = Fasta(file)
    balancer       = FastaBalancer(fasta)
    filtered_fasta = balancer.balance_dataset_sparse(target_size=20)