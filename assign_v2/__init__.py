import pathlib

DATA_DIR = pathlib.Path(__file__).parent.parent / "data"

FASTA_DIR       = DATA_DIR / "fasta"
MMCIF_DIR       = DATA_DIR / "mmcif"
OUTPUT_DIR      = DATA_DIR / "output"
POLYMERS_DIR    = DATA_DIR / "polymers"
TUNNEL_DIR      = DATA_DIR / "tunnel"
ASSIGN_DIR         = DATA_DIR / "assign_v2"

MMCIF_DIR.mkdir(parents=True, exist_ok=True)