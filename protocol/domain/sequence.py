from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from protocol.domain.types import *

class SequenceMappingContainer(Chain):
    """ 
    This class is taken from RibosomeXYZ:
    A. Kushner, A.S. Petrov, K. Dao Duc, (2022) “RiboXYZ: A comprehensive database for ribosome structures”,  Nucleic Acids Research, gkac939
    (https://github.com/rtviii/riboxyz/blob/master/ribctl/lib/libseq.py#L20)
    
    A container for keeping track of the correspondence between the structural and sequence indices within a give polymer chain.
    \nStructural data(`mmcif`) frequently has unresolved and modified residues adheres to the author's numbering (which is arbitrary for all intents and purposes, ex. starts at, say, 7).
    More here:https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
    
    This container keeps pointers between the naive (convenient) indices and the structural residues. 
    I refer to the following in the code:

    - **auth_seq_id** is the author-assigned residue number and is frequently used to refer to structural components. (ex. in Molstar)
    - **flat_index** is the the straightforward arithmetic (0-start) numbering. I produce it by removing all the modified residues from the "primary sequence"
    - **primary_sequence** is the sequence of structural residues (as packed into the BioPython `Chain` object) represented as a string.
    - **flat_sequence** is the **primary_sequence** with the modified residues removed, represented as a string.
    
    There is lots to optimize in this code (it builds index->Residue<object> maps by enumeration),
    but ideally this is taken care of at the parser level or at the deposition level.

    Again, see more: 
    - https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
    - https://bioinformatics.stackexchange.com/questions/14210/pdb-residue-numbering
    - https://bioinformatics.stackexchange.com/questions/20458/how-is-the-canonical-version-entity-poly-pdbx-seq-one-letter-code-obtaine
    - https://www.biostars.org/p/9588718/
    - https://stackoverflow.com/questions/45466408/biopython-resseq-doesnt-match-pdb-file
    
    """

    chain                        : Chain

    flat_index_to_residue_map    : dict[int, Residue]
    auth_seq_id_to_flat_index_map: dict[int, int]

    @property
    def primary_sequence(self, represent_noncanonical_as:str=".") -> tuple[str, dict[int,int]]:
        seq = ""
        auth_seq_id_to_primary_ix = {}
        for ix, residue in enumerate(self.chain.get_residues()):
            if residue.resname in [*AMINO_ACIDS.keys()]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
            elif residue.resname in [*NUCLEOTIDES]:
                seq = seq + residue.resname
            else:
                seq = seq + represent_noncanonical_as
            auth_seq_id_to_primary_ix[residue.get_id()[1]] = ix

        return seq, auth_seq_id_to_primary_ix

    @property
    def flat_sequence(self) -> tuple[str, dict[int,Residue], dict[int,int]]:
        res: list[Residue] = [*self.chain.get_residues()]

        flat_index_to_residue_map     = {}
        auth_seq_id_to_flat_index_map = {}
        seq                           = ""
        flat_index                    = 0

        for residue in res:
            if residue.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
                flat_index_to_residue_map[flat_index] = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = flat_index
                flat_index += 1
            else:
                continue
        return seq, flat_index_to_residue_map, auth_seq_id_to_flat_index_map

    def __init__(self, chain: Chain):
        self.chain = chain