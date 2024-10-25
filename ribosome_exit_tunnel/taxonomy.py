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
        # lin = list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )
        # lin = list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )
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