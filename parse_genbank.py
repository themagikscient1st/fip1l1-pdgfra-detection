"""
parse_genbank.py

Parse the FIP1L1 and PDGFRA GenBank files and report the coordinates
of each exon, with particular attention to the breakpoint zones:
    FIP1L1 exon 10 / intron 10 (5' breakpoint zone)
    PDGFRA exon 12              (3' breakpoint zone, always)

This is a read-only inspection script. It does not extract sequence
or make any decisions about breakpoint placement; it just tells us
what's in the GenBank files so we can make informed choices later.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# The directory containing the .gb files. Path() handles spaces in directory
# names automatically, so we don't need to escape "FIP1L1_PDGFRA FASTAs".
WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"

# Each file maps to a gene name and a target exon we want to flag.
# The chromosomal_offset is the start coordinate of the cropped region on
# NC_000004.12; we add it to record-local coordinates to get chromosomal ones.
FILES = {
    "FIP1L1.gb": {
        "gene": "FIP1L1",
        "target_exon": 10,
        "chromosomal_offset": 53_370_000,
    },
    "PDGFRA.gb": {
        "gene": "PDGFRA",
        "target_exon": 12,
        "chromosomal_offset": 54_225_000,
    },
}


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def get_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    """
    Pull a single qualifier value from a feature.

    GenBank qualifiers are stored as lists (because a feature can have
    multiple values for the same qualifier), so feature.qualifiers["gene"]
    returns ["FIP1L1"] rather than "FIP1L1". This helper unwraps the list
    and returns the first value, or a default if the qualifier is absent.
    """
    values = feature.qualifiers.get(key, [default])
    return values[0] if values else default


def feature_belongs_to_gene(feature: SeqFeature, gene_name: str) -> bool:
    """
    Check whether a feature's /gene= qualifier matches our target gene.

    The cropped region may contain features from neighboring genes (within
    our flanking sequence), so we filter to just the gene of interest.
    """
    return get_qualifier(feature, "gene") == gene_name


def format_coordinates(
    start: int,
    end: int,
    offset: int,
) -> str:
    """
    Format a coordinate range in both record-local and chromosomal systems.

    Returns something like:
        "7,641-90,862 (chr4: 53,377,641-53,460,862)"

    Note: Biopython uses 0-based, half-open coordinates internally
    (Python convention), so we add 1 to start positions when reporting
    in 1-based GenBank style.
    """
    local_start = start + 1  # convert 0-based to 1-based
    local_end = end           # half-open end is already correct for 1-based inclusive
    chrom_start = local_start + offset
    chrom_end = local_end + offset
    return (
        f"{local_start:,}-{local_end:,} "
        f"(chr4: {chrom_start:,}-{chrom_end:,})"
    )


# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

def analyze_gene_record(
    filepath: Path,
    gene_name: str,
    target_exon: int,
    chromosomal_offset: int,
) -> None:
    """
    Read one GenBank file and print a structured report of its contents.
    """
    print(f"\n{'=' * 70}")
    print(f"FILE: {filepath.name}")
    print(f"GENE: {gene_name}")
    print(f"{'=' * 70}")

    # SeqIO.read expects exactly one record in the file, which is what we have.
    # If there were multiple records we'd use SeqIO.parse instead.
    record: SeqRecord = SeqIO.read(filepath, "genbank")

    print(f"Record ID:        {record.id}")
    print(f"Description:      {record.description}")
    print(f"Length:           {len(record.seq):,} bp")
    print(f"Total features:   {len(record.features)}")
    print(f"Chromosomal range: chr4:{chromosomal_offset + 1:,}-"
          f"{chromosomal_offset + len(record.seq):,}")

    # -------------------------------------------------------------------------
    # Step 1: Find the gene feature for our target gene
    # -------------------------------------------------------------------------
    gene_features = [
        f for f in record.features
        if f.type == "gene" and feature_belongs_to_gene(f, gene_name)
    ]

    if not gene_features:
        print(f"\n  WARNING: No gene feature found for {gene_name}.")
        print(f"  Check the .gb file — the gene may be annotated under a synonym.")
        return

    gene_feature = gene_features[0]
    print(f"\n--- Gene feature: {gene_name} ---")
    print(f"  Boundaries: {format_coordinates(int(gene_feature.location.start), int(gene_feature.location.end), chromosomal_offset)}")
    print(f"  Strand:     {'+' if gene_feature.location.strand == 1 else '-'}")

    # -------------------------------------------------------------------------
    # Step 2: Find all exon features for our target gene
    # -------------------------------------------------------------------------
    # Each exon feature has a /number= qualifier indicating which exon it is.
    # Note: GenBank may include exons from multiple transcripts. For now we
    # collect all exons for our gene and report them; if there are duplicates
    # we'll see it in the output and can decide how to handle it.
    exon_features = [
        f for f in record.features
        if f.type == "exon" and feature_belongs_to_gene(f, gene_name)
    ]

    print(f"\n--- Exon features ({len(exon_features)} total) ---")

    if not exon_features:
        print(f"  WARNING: No exon features found for {gene_name}.")
        print(f"  This may mean the GenBank record only annotates mRNA features,")
        print(f"  in which case we'd need to extract exon coordinates from the")
        print(f"  mRNA feature's join() statement instead.")
        return

    # Report each exon's coordinates and flag the target
    target_exon_feature = None
    for exon in exon_features:
        exon_number_str = get_qualifier(exon, "number", "?")
        try:
            exon_number = int(exon_number_str)
        except ValueError:
            exon_number = None

        is_target = (exon_number == target_exon)
        marker = "  >>> TARGET" if is_target else ""

        coords = format_coordinates(
            int(exon.location.start),
            int(exon.location.end),
            chromosomal_offset,
        )
        length = int(exon.location.end) - int(exon.location.start)

        print(f"  Exon {exon_number_str:>3}: {coords}  [{length:,} bp]{marker}")

        if is_target:
            target_exon_feature = exon

    # -------------------------------------------------------------------------
    # Step 3: Report the breakpoint zone in detail
    # -------------------------------------------------------------------------
    if target_exon_feature is None:
        print(f"\n  WARNING: Target exon {target_exon} not found in the annotations.")
        return

    print(f"\n--- Breakpoint zone: {gene_name} exon {target_exon} ---")
    target_start = int(target_exon_feature.location.start)
    target_end = int(target_exon_feature.location.end)
    print(f"  Exon {target_exon} coordinates: "
          f"{format_coordinates(target_start, target_end, chromosomal_offset)}")
    print(f"  Exon {target_exon} length:      {target_end - target_start:,} bp")

    # If this is FIP1L1, also locate intron 10 (between exons 10 and 11).
    # The breakpoint sits inside this intron in Type A fusions.
    if gene_name == "FIP1L1":
        next_exon_features = [
            f for f in exon_features
            if get_qualifier(f, "number") == str(target_exon + 1)
        ]
        if next_exon_features:
            next_exon = next_exon_features[0]
            intron_start = target_end
            intron_end = int(next_exon.location.start)
            print(f"\n  Intron {target_exon} (between exons {target_exon} "
                  f"and {target_exon + 1}):")
            print(f"    Coordinates: "
                  f"{format_coordinates(intron_start, intron_end, chromosomal_offset)}")
            print(f"    Length:      {intron_end - intron_start:,} bp")
            print(f"    --> This is where Type A breakpoints sit "
                  f"(Walz et al., 2009)")


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

def main() -> None:
    if not WORKING_DIR.exists():
        print(f"ERROR: Working directory not found: {WORKING_DIR}")
        return

    for filename, info in FILES.items():
        filepath = WORKING_DIR / filename
        if not filepath.exists():
            print(f"ERROR: File not found: {filepath}")
            continue
        analyze_gene_record(
            filepath=filepath,
            gene_name=info["gene"],
            target_exon=info["target_exon"],
            chromosomal_offset=info["chromosomal_offset"],
        )


if __name__ == "__main__":
    main()