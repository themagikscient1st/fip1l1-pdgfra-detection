"""
detect_breakpoint.py

Custom breakpoint detection from aligned BAM, designed as an analytical
demonstration for the FIP1L1::PDGFRA assay capstone. Operates on either
the deletion-validation BAM (5 kb synthetic FIP1L1 intronic deletion)
or the fusion BAM (FIP1L1::PDGFRA chimeric construct).

Two evidence channels:
  1. Within-read CIGAR deletions: reads whose alignment includes a
     large 'D' operation indicate the deletion sits within their
     spanned interval. Position, size, and read identity are recorded.
  2. Split-read alignments: reads whose primary + supplementary
     alignment endpoints flank a defined gap indicate the read
     bridges a structural rearrangement.

Both channels are clustered to identify breakpoint coordinates with
multi-read support, and a summary is reported with confidence metrics.

This is the analytical core of what Sniffles2 does internally; we
reproduce it explicitly here to demonstrate that the synthetic
deletion is detectable by the same evidence pattern Sniffles uses.

Usage:
    python detect_breakpoint.py
"""

from collections import defaultdict
from pathlib import Path
import pysam


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
BAM_PATH = WORKING_DIR / "sniffles_test" / "deletion_aligned.sorted.bam"

# Detection parameters
MIN_DELETION_SIZE = 1000  # ignore CIGAR D operations smaller than this (within-read indels)
CLUSTER_WINDOW = 50       # bp tolerance for clustering breakpoint coordinates
MIN_SUPPORT = 3           # minimum reads supporting a breakpoint to call it


# -----------------------------------------------------------------------------
# Channel 1: CIGAR-based deletions
# -----------------------------------------------------------------------------

def find_cigar_deletions(bam_path: Path) -> list:
    """
    Scan every primary alignment for large D operations in its CIGAR.

    Returns a list of (chrom, start, end, size, read_name) tuples for
    each large deletion event found.

    pysam.AlignmentFile lets us iterate over BAM records as objects,
    where each .cigartuples is a list of (op_code, length) pairs:
      0=M (match), 1=I (insert), 2=D (delete), 4=S (soft-clip), etc.
    """
    deletions = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            # Skip secondary and supplementary alignments; we want to scan
            # primary alignments for within-read deletions
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                continue
            if read.cigartuples is None:
                continue

            # Walk through the CIGAR, tracking position on the reference
            # so we know where each D operation sits in chromosomal coords
            ref_pos = read.reference_start  # 0-based start of alignment

            for op_code, op_length in read.cigartuples:
                # M = 0, =/X = 7/8: consume both query and reference
                # D = 2, N = 3: consume reference only
                # I = 1, S = 4, H = 5: consume query only

                if op_code == 2 and op_length >= MIN_DELETION_SIZE:
                    # Found a D operation big enough to record
                    deletions.append((
                        read.reference_name,
                        ref_pos,                 # 0-based start of deletion
                        ref_pos + op_length,     # 0-based end (exclusive)
                        op_length,
                        read.query_name,
                    ))

                # Advance ref_pos for ops that consume reference
                if op_code in (0, 2, 3, 7, 8):
                    ref_pos += op_length

    return deletions


# -----------------------------------------------------------------------------
# Channel 2: Split-read alignments
# -----------------------------------------------------------------------------

def find_split_alignments(bam_path: Path) -> list:
    """
    Identify pairs of (primary, supplementary) alignments that flank
    a putative breakpoint.

    Returns a list of (chrom, breakpoint_pos_5p, breakpoint_pos_3p,
    gap_size, read_name) tuples where the gap_size is the distance
    between the end of the primary alignment and the start of the
    supplementary alignment.
    """
    # First pass: collect all primary alignments by read name
    primaries = {}
    supplementaries = defaultdict(list)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            if read.is_secondary:
                continue

            entry = (
                read.reference_name,
                read.reference_start,
                read.reference_end,  # 0-based exclusive
            )

            if read.is_supplementary:
                supplementaries[read.query_name].append(entry)
            else:
                primaries[read.query_name] = entry

    # Second pass: match up primaries with their supplementaries
    splits = []
    for read_name, prim in primaries.items():
        for supp in supplementaries.get(read_name, []):
            # Only consider same-chromosome pairs (we're looking for deletions)
            if prim[0] != supp[0]:
                continue

            # Determine which alignment is the 5' end and which is the 3' end
            # based on reference start position
            if prim[1] < supp[1]:
                five_prime_end = prim[2]      # end of primary
                three_prime_start = supp[1]   # start of supplementary
            else:
                five_prime_end = supp[2]
                three_prime_start = prim[1]

            gap = three_prime_start - five_prime_end
            if gap > MIN_DELETION_SIZE:
                splits.append((
                    prim[0],
                    five_prime_end,
                    three_prime_start,
                    gap,
                    read_name,
                ))

    return splits


# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------

def cluster_breakpoints(events: list, label: str) -> list:
    """
    Cluster events by breakpoint position using a simple window-based
    scheme: events whose breakpoints fall within CLUSTER_WINDOW of each
    other are grouped.

    Events tuple format: (chrom, start, end, size, read_name)

    Returns a list of clusters, each a dict with consensus breakpoint
    coordinates and supporting reads.
    """
    if not events:
        return []

    # Sort events by start position
    sorted_events = sorted(events, key=lambda e: (e[0], e[1]))

    clusters = []
    current_cluster = [sorted_events[0]]

    for event in sorted_events[1:]:
        last_event = current_cluster[-1]
        # Same chromosome and start positions within the cluster window?
        if (event[0] == last_event[0]
                and abs(event[1] - last_event[1]) <= CLUSTER_WINDOW
                and abs(event[2] - last_event[2]) <= CLUSTER_WINDOW):
            current_cluster.append(event)
        else:
            clusters.append(_summarize_cluster(current_cluster, label))
            current_cluster = [event]

    clusters.append(_summarize_cluster(current_cluster, label))
    return clusters


def _summarize_cluster(events: list, label: str) -> dict:
    """Compute consensus stats for a cluster of breakpoint events."""
    starts = [e[1] for e in events]
    ends = [e[2] for e in events]
    sizes = [e[3] for e in events]

    return {
        "evidence": label,
        "chrom": events[0][0],
        "start_consensus": sum(starts) // len(starts),
        "start_min": min(starts),
        "start_max": max(starts),
        "end_consensus": sum(ends) // len(ends),
        "size_consensus": sum(sizes) // len(sizes),
        "n_supporting_reads": len(events),
        "supporting_read_names": [e[4] for e in events],
    }


# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------

def print_cluster(cluster: dict) -> None:
    """Pretty-print one called breakpoint."""
    if cluster["n_supporting_reads"] < MIN_SUPPORT:
        return

    print(f"\n  --- Breakpoint detected ---")
    print(f"    Evidence type:    {cluster['evidence']}")
    print(f"    Reference contig: {cluster['chrom']}")
    print(f"    Start (consensus): {cluster['start_consensus']:,} "
          f"(range: {cluster['start_min']:,}-{cluster['start_max']:,})")
    print(f"    End (consensus):   {cluster['end_consensus']:,}")
    print(f"    Size:              {cluster['size_consensus']:,} bp")
    print(f"    Supporting reads:  {cluster['n_supporting_reads']}")


def main() -> None:
    print("=" * 70)
    print("Custom breakpoint detection on synthetic deletion BAM")
    print("=" * 70)
    print(f"\nInput BAM: {BAM_PATH}")
    print(f"Detection parameters:")
    print(f"  Minimum deletion size: {MIN_DELETION_SIZE:,} bp")
    print(f"  Cluster window:        {CLUSTER_WINDOW} bp")
    print(f"  Minimum support:       {MIN_SUPPORT} reads")

    # Channel 1: CIGAR-based deletions
    print(f"\n--- Channel 1: scanning CIGAR strings for D operations ---")
    cigar_deletions = find_cigar_deletions(BAM_PATH)
    print(f"  Found {len(cigar_deletions)} large deletion events in primary alignments")

    cigar_clusters = cluster_breakpoints(cigar_deletions, "CIGAR deletion")
    confident_cigar = [c for c in cigar_clusters if c["n_supporting_reads"] >= MIN_SUPPORT]
    print(f"  Clustered into {len(cigar_clusters)} candidate breakpoints, "
          f"{len(confident_cigar)} passing minimum support threshold")

    # Channel 2: Split-read alignments
    print(f"\n--- Channel 2: scanning for split-read alignments ---")
    split_alignments = find_split_alignments(BAM_PATH)
    print(f"  Found {len(split_alignments)} primary+supplementary pairs with same-chromosome gap")

    split_clusters = cluster_breakpoints(split_alignments, "Split alignment")
    confident_splits = [c for c in split_clusters if c["n_supporting_reads"] >= MIN_SUPPORT]
    print(f"  Clustered into {len(split_clusters)} candidate breakpoints, "
          f"{len(confident_splits)} passing minimum support threshold")

    # Combined report
    print(f"\n" + "=" * 70)
    print(f"FINAL CALLS (passing minimum support of {MIN_SUPPORT} reads)")
    print(f"=" * 70)

    all_confident = confident_cigar + confident_splits

    if not all_confident:
        print(f"\n  No breakpoints detected with sufficient support.")
        return

    for cluster in all_confident:
        print_cluster(cluster)

    print(f"\n" + "=" * 70)
    print(f"Total breakpoint calls: {len(all_confident)}")
    print(f"=" * 70)


if __name__ == "__main__":
    main()