from ete3 import NCBITaxa
from typing_extensions import Literal
import typing

'''
Code in this class is taken from 'riboxyz':
A. Kushner, A.S. Petrov, K. Dao Duc, (2022) “RiboXYZ: A comprehensive database for ribosome structures”,  Nucleic Acids Research, gkac939
(https://github.com/rtviii/riboxyz/blob/331dbef5adf87f8393bdde0ee770414ea3c28e23/ribctl/lib/libtax.py)
'''
TAXID_BACTERIA  = 2
TAXID_EUKARYOTA = 2759
TAXID_ARCHAEA   = 2157

PhylogenyRank   = Literal[
    "superkingdom",
    "phylum",
    "class",
    "order",
    "clade",
    "family",
    "genus",
    "species",
    "strain",
    "isolate",
    "subspecies",
    "no rank",
    "suborder",
    "kingdom",
    "subfamily",
    "subgenus",
    "subphylum",
    "infraorder",
    "superorder",
    "superclass",
    "superfamily",
    "parvorder",
    "cohort",
    "infraclass",
    "subclass",
    "subkingdom",
    "species group",
    "tribe",
    "species subgroup",
    "subcohort",
    "subtribe",
]

ncbi = NCBITaxa()

class TaxId:

    @staticmethod
    def get_lineage(
        taxid, include_only: None | list[PhylogenyRank] = None
    ) -> list[int]:
        """Return ncbi lineage, except filter out the ranks that are not among the @PhylogenyRank."""
        lin = ncbi.get_lineage(taxid)
        if include_only is not None:
            return list(filter(lambda x: TaxId.rank(x) in include_only, lin))
        return lin if lin is not None else []

    @staticmethod
    def is_descendant_of(parent_taxid: int, target_taxid: int) -> bool:
        lineage = ncbi.get_lineage(target_taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        return False if parent_taxid not in lineage else True
    
    @staticmethod
    def superkingdom(
        taxid: int,
    ) -> typing.Literal["bacteria", "eukaryota", "archaea", "virus"]:
        match (
            TaxId.is_descendant_of(TAXID_EUKARYOTA, taxid),
            TaxId.is_descendant_of(TAXID_BACTERIA, taxid),
            TaxId.is_descendant_of(TAXID_ARCHAEA, taxid),
        ):
            case (False, False, True):
                return "archaea"
            case (False, True, False):
                return "bacteria"
            case (True, False, False):
                return "eukaryota"
            case (False, False, False):
                print("Probably a virus")
                return "virus"
            case _:
                raise ValueError(
                    "Taxid {} is not a descendant of any of the three domains".format(
                        taxid
                    )
                )
                
    @staticmethod
    def rank(taxid: int) -> PhylogenyRank:
        """Given a @taxid, return the rank of the taxid"""
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage)[taxid]
    
    @staticmethod
    def get_name(taxid):
        return list(ncbi.get_taxid_translator([taxid]).values())[0]
    
    @staticmethod
    def coerce_to_rank(taxid: int, target_rank: PhylogenyRank) -> int | None:
        """Given a @taxid and a @rank, return the taxid of the first ancestor of @taxid that is at @rank"""
        lineage = ncbi.get_lineage(taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        for item in lineage:
            rank = ncbi.get_rank([item])[item]
            if rank == target_rank:
                return item

        raise IndexError("Taxid {} does not have a {} level".format(taxid, target_rank))
    
    @staticmethod
    def coerce_all_to_rank(taxids: list[int], level: PhylogenyRank) -> list[int]:
        """Given a list of taxids, return a list of the same taxids but coerced to the given lineage level(rank)."""
        new = []
        for taxid in taxids:
            try:
                new.append(TaxId.coerce_to_rank(taxid, level))
            except Exception as e:
                print(e)
                raise Exception(e)
        return new