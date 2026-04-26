"""
build_synthetic_reference.py

Build a synthetic reference for FIP1L1::PDGFRA fusion detection.

Two contigs are produced:

  chr4_normal — A wild-type-like contig representing the unfused
    architecture. The FIP1L1 region and the PDGFRA region are joined
    via a short representative intervening sequence (the actual ~864 kb
    deleted interval is too large to include in full).

  chr4_fusion — A fusion contig representing what fusion-positive
    genomic DNA looks like at the locus. FIP1L1 sequence runs from the
    start of our cropped region through the 5' breakpoint
    (chr4:53,411,065). PDGFRA sequence picks up at the 3' breakpoint
    (chr4:54,274,884) and runs through the end of our PDGFRA crop.

Both contigs use GRCh38 reference sequence content. Sample 16 from Walz
et al. (2009) carries one or more variants relative to GRCh38 within
the breakpoint zone (one SNV in PDGFRA exon 12, position 46; possible
SNVs in the FIP1L1 intron 10 anchor region). These are not reproduced
in the synthetic reference; the breakpoint positions are anchored on
sample 16's published configuration but the surrounding sequence is the
GRCh38 reference.

Output: synthetic_reference/synthetic_reference.fasta (multi-record)
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
FIP1L1_FILE = WORKING_DIR / "FIP1L1.gb"
PDGFRA_FILE = WORKING_DIR / "PDGFRA.gb"

# Where to write the synthetic reference
OUTPUT_DIR = WORKING_DIR / "synthetic_reference"
OUTPUT_FILE = OUTPUT_DIR / "synthetic_reference.fasta"

# Chromosomal start of each cropped region on NC_000004.12
FIP1L1_OFFSET = 53_370_000
PDGFRA_OFFSET = 54_225_000

# Breakpoint coordinates derived from Walz et al. (2009) Sup. Table 4 sample 16
FIP1L1_BREAKPOINT_CHROM = 53_411_065  # last retained base on FIP1L1 side
PDGFRA_BREAKPOINT_CHROM = 54_274_885  # first retained base on PDGFRA side

# For the chr4_normal contig: how much "intervening" sequence to include
# between the FIP1L1 and PDGFRA portions. The real deletion spans ~864 kb,
# which is too long to include in a synthetic test reference. We use a
# short stand-in (1 kb of N's) to keep the contig structure parallel to
# chr4_fusion but with the deletion-equivalent region intact.
INTERVENING_LENGTH = 1000


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def chrom_to_local(chrom_pos: int, file_offset: int) -> int:
    """Convert 1-based chromosomal coordinate to 0-based record-local."""
    return chrom_pos - file_offset   # was: chrom_pos - file_offset - 1


def load_sequence(filepath: Path) -> str:
    """Load the raw sequence from a GenBank file as an uppercase string."""
    record = SeqIO.read(filepath, "genbank")
    return str(record.seq).upper()


# -----------------------------------------------------------------------------
# Contig construction
# -----------------------------------------------------------------------------

def build_normal_contig(fip1l1_seq: str, pdgfra_seq: str) -> str:
    """
    Build the wild-type-equivalent contig.

    Structure:
        [full FIP1L1 region]  +  [N-padded intervening]  +  [full PDGFRA region]

    The N-padded intervening section is a placeholder for the ~864 kb
    deleted interval. We don't try to represent the real intervening
    genomic sequence (CHIC2, LNX1, etc.) because including 864 kb of real
    sequence in a synthetic test reference would defeat the purpose of
    keeping the file lightweight and tractable for read simulation.
    """
    intervening = "N" * INTERVENING_LENGTH
    return fip1l1_seq + intervening + pdgfra_seq


def build_fusion_contig(fip1l1_seq: str, pdgfra_seq: str) -> tuple[str, int]:
    """
    Build the fusion contig.

    Structure:
        [FIP1L1 from start of crop through breakpoint]  +
        [PDGFRA from breakpoint through end of crop]

    The junction position (0-based, in the resulting fusion contig) is
    returned so we can verify the splice point and document the
    junction coordinate.
    """
    # Convert chromosomal breakpoint coordinates to record-local
    fip1l1_breakpoint_local = chrom_to_local(
        FIP1L1_BREAKPOINT_CHROM, FIP1L1_OFFSET
    )
    pdgfra_breakpoint_local = chrom_to_local(
        PDGFRA_BREAKPOINT_CHROM, PDGFRA_OFFSET
    )

    # FIP1L1 contribution: positions [0, breakpoint] inclusive.
    # We include the breakpoint base itself as the last retained FIP1L1 base,
    # so we slice [:breakpoint + 1] to get bases 0 through breakpoint inclusive.
    fip1l1_contribution = fip1l1_seq[:fip1l1_breakpoint_local + 1]

    # PDGFRA contribution: positions [breakpoint, end). The breakpoint base
    # is the first retained PDGFRA base, so we slice from there to the end.
    pdgfra_contribution = pdgfra_seq[pdgfra_breakpoint_local:]

    fusion_contig = fip1l1_contribution + pdgfra_contribution
    junction_position = len(fip1l1_contribution)  # 0-based index of first PDGFRA base

    return fusion_contig, junction_position


# -----------------------------------------------------------------------------
# Verification
# -----------------------------------------------------------------------------

def verify_fusion_junction(fusion_seq: str, junction_pos: int) -> None:
    """
    Print the sequence on either side of the fusion junction so we can
    eyeball it against expectations.

    Expected: FIP1L1 intron 10 sequence on the left, PDGFRA exon 12
    sequence (starting near 'CCCAGA' — note: GRCh38 reference, with the
    sample-16 SNV resolved as the reference allele) on the right.
    """
    print(f"\n--- Fusion junction verification ---")
    print(f"  Junction at position {junction_pos:,} of fusion contig "
          f"({len(fusion_seq):,} bp total)")
    print(f"\n  Sequence context (50 bp each side):")
    print(f"    FIP1L1 side: ...{fusion_seq[junction_pos-50:junction_pos]}|")
    print(f"    PDGFRA side: |{fusion_seq[junction_pos:junction_pos+50]}...")
    print(f"\n  (The | marks the splice point. FIP1L1 intronic sequence "
          f"should be on the left, PDGFRA exon 12 sequence on the right.)")


# -----------------------------------------------------------------------------
# FASTA writing
# -----------------------------------------------------------------------------

def write_fasta(records: list, output_path: Path) -> None:
    """Write multiple SeqRecords to a single multi-record FASTA file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as handle:
        SeqIO.write(records, handle, "fasta")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("Building synthetic reference: chr4_normal + chr4_fusion")
    print("=" * 70)

    # Load source sequences
    print(f"\nLoading source sequences...")
    fip1l1_seq = load_sequence(FIP1L1_FILE)
    pdgfra_seq = load_sequence(PDGFRA_FILE)
    print(f"  FIP1L1 region: {len(fip1l1_seq):,} bp "
          f"(chr4:{FIP1L1_OFFSET+1:,}-{FIP1L1_OFFSET+len(fip1l1_seq):,})")
    print(f"  PDGFRA region: {len(pdgfra_seq):,} bp "
          f"(chr4:{PDGFRA_OFFSET+1:,}-{PDGFRA_OFFSET+len(pdgfra_seq):,})")

    # Build chr4_normal
    print(f"\nBuilding chr4_normal contig...")
    normal_seq = build_normal_contig(fip1l1_seq, pdgfra_seq)
    print(f"  Total length: {len(normal_seq):,} bp")
    print(f"    FIP1L1 portion:    {len(fip1l1_seq):,} bp")
    print(f"    Intervening (N's): {INTERVENING_LENGTH:,} bp")
    print(f"    PDGFRA portion:    {len(pdgfra_seq):,} bp")

    # Build chr4_fusion
    print(f"\nBuilding chr4_fusion contig...")
    fusion_seq, junction_position = build_fusion_contig(fip1l1_seq, pdgfra_seq)
    fip1l1_breakpoint_local = chrom_to_local(
        FIP1L1_BREAKPOINT_CHROM, FIP1L1_OFFSET
    )
    pdgfra_breakpoint_local = chrom_to_local(
        PDGFRA_BREAKPOINT_CHROM, PDGFRA_OFFSET
    )
    print(f"  Total length: {len(fusion_seq):,} bp")
    print(f"    FIP1L1 contribution: {fip1l1_breakpoint_local + 1:,} bp "
          f"(crop start through chr4:{FIP1L1_BREAKPOINT_CHROM:,})")
    print(f"    PDGFRA contribution: {len(pdgfra_seq) - pdgfra_breakpoint_local:,} bp "
          f"(chr4:{PDGFRA_BREAKPOINT_CHROM:,} through crop end)")
    print(f"    Junction at fusion-contig position: {junction_position:,}")

    # Verify the junction looks right
    verify_fusion_junction(fusion_seq, junction_position)

    # Build SeqRecord objects with informative descriptions
    normal_record = SeqRecord(
        Seq(normal_seq),
        id="chr4_normal",
        description=(
            f"Wild-type-equivalent: FIP1L1 region (chr4:{FIP1L1_OFFSET+1:,}-"
            f"{FIP1L1_OFFSET+len(fip1l1_seq):,}) + {INTERVENING_LENGTH} bp N-padded "
            f"intervening + PDGFRA region (chr4:{PDGFRA_OFFSET+1:,}-"
            f"{PDGFRA_OFFSET+len(pdgfra_seq):,})"
        ),
    )

    fusion_record = SeqRecord(
        Seq(fusion_seq),
        id="chr4_fusion",
        description=(
            f"FIP1L1::PDGFRA fusion at Walz et al. (2009) sample 16 "
            f"breakpoint (FIP1L1 5' bp chr4:{FIP1L1_BREAKPOINT_CHROM:,}; "
            f"PDGFRA 3' bp chr4:{PDGFRA_BREAKPOINT_CHROM:,}; "
            f"junction at fusion-contig position {junction_position:,})"
        ),
    )

    # Write to FASTA
    write_fasta([normal_record, fusion_record], OUTPUT_FILE)
    print(f"\n--- Output ---")
    print(f"  Wrote 2 records to: {OUTPUT_FILE}")
    print(f"    chr4_normal: {len(normal_seq):,} bp")
    print(f"    chr4_fusion: {len(fusion_seq):,} bp")
    print(f"\n  Total file size on disk:")

    # Print file size
    size_bytes = OUTPUT_FILE.stat().st_size
    size_mb = size_bytes / 1024 / 1024
    print(f"    {size_bytes:,} bytes ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()