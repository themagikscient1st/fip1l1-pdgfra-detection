"""
build_deletion_reference.py

Build a "deletion-bearing" version of the FIP1L1 region: take the
wild-type FIP1L1 sequence and remove a defined 5 kb interval. Simulated
reads from this contig, aligned back to the wild-type FIP1L1 sequence,
will exhibit split alignments that Sniffles can cluster into a clean
deletion call.

The deletion is placed at chr4:53,400,000-53,405,000 (record-local
positions 30,000-35,000 in our cropped GenBank file). This puts it
inside FIP1L1 intron 7 (NCBI numbering) -- well away from annotated
exons and any other features that could complicate interpretation.

Outputs:
    sniffles_test/wildtype.fasta              -- the unmodified FIP1L1 region
    sniffles_test/deletion_bearing.fasta      -- same region with 5 kb removed
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
FIP1L1_FILE = WORKING_DIR / "FIP1L1.gb"

OUTPUT_DIR = WORKING_DIR / "sniffles_test"
WILDTYPE_FASTA = OUTPUT_DIR / "wildtype.fasta"
DELETION_FASTA = OUTPUT_DIR / "deletion_bearing.fasta"

# Chromosomal coordinates of the deletion
DELETION_CHROM_START = 53_400_000  # first deleted base (1-based, inclusive)
DELETION_CHROM_END = 53_405_000    # last deleted base (1-based, inclusive)
FIP1L1_OFFSET = 53_370_000


def main() -> None:
    print("=" * 70)
    print("Building wild-type and deletion-bearing FIP1L1 references")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load the FIP1L1 sequence
    record = SeqIO.read(FIP1L1_FILE, "genbank")
    sequence = str(record.seq).upper()
    seq_length = len(sequence)
    print(f"\nLoaded FIP1L1 region: {seq_length:,} bp "
          f"(chr4:{FIP1L1_OFFSET+1:,}-{FIP1L1_OFFSET+seq_length:,})")

    # Convert deletion coordinates to 0-based record-local
    # Chromosomal pos N corresponds to record-local pos (N - offset - 1)
    del_start_local = DELETION_CHROM_START - FIP1L1_OFFSET - 1  # 0-based, inclusive
    del_end_local = DELETION_CHROM_END - FIP1L1_OFFSET           # 0-based, exclusive
    deletion_size = del_end_local - del_start_local
    print(f"Deletion coordinates: chr4:{DELETION_CHROM_START:,}-"
          f"{DELETION_CHROM_END:,} ({deletion_size:,} bp)")
    print(f"Record-local: {del_start_local:,}-{del_end_local:,}")

    # Build wild-type record (just the unmodified sequence)
    wildtype_record = SeqRecord(
        Seq(sequence),
        id="chr4_fip1l1_wt",
        description=(
            f"FIP1L1 wild-type region (chr4:{FIP1L1_OFFSET+1:,}-"
            f"{FIP1L1_OFFSET+seq_length:,}, GRCh38)"
        ),
    )

    # Build deletion-bearing record: concatenate everything before the deletion
    # with everything after the deletion
    deletion_seq = sequence[:del_start_local] + sequence[del_end_local:]
    deletion_record = SeqRecord(
        Seq(deletion_seq),
        id="chr4_fip1l1_del",
        description=(
            f"FIP1L1 region with synthetic 5 kb deletion at "
            f"chr4:{DELETION_CHROM_START:,}-{DELETION_CHROM_END:,}"
        ),
    )

    # Write outputs
    with open(WILDTYPE_FASTA, "w") as fh:
        SeqIO.write([wildtype_record], fh, "fasta")
    with open(DELETION_FASTA, "w") as fh:
        SeqIO.write([deletion_record], fh, "fasta")

    print(f"\nWrote:")
    print(f"  Wild-type:        {WILDTYPE_FASTA} ({len(sequence):,} bp)")
    print(f"  Deletion-bearing: {DELETION_FASTA} ({len(deletion_seq):,} bp)")
    print(f"  Length difference: {len(sequence) - len(deletion_seq):,} bp")

    # Sanity check: the junction sequence in the deletion-bearing reference
    # should be the bases just before the deletion concatenated to the bases
    # just after.
    junction_pos = del_start_local
    print(f"\nJunction context (50 bp each side of the splice point):")
    print(f"  Before:  ...{deletion_seq[junction_pos-50:junction_pos]}|")
    print(f"  After:   |{deletion_seq[junction_pos:junction_pos+50]}...")

    print(f"\n{'=' * 70}")
    print(f"DONE")
    print(f"{'=' * 70}")
    print(f"Next: simulate reads from {DELETION_FASTA.name},")
    print(f"      align to {WILDTYPE_FASTA.name},")
    print(f"      run Sniffles to detect the {deletion_size:,} bp deletion.")


if __name__ == "__main__":
    main()