"""
plot_validation.py

Generate a three-panel figure showing the synthetic deletion detected
by all three evidence channels:

  Panel A: Per-base coverage from mosdepth
  Panel B: Per-read alignment intervals (primary + supplementary)
  Panel C: Breakpoint calls from detect_breakpoint.py overlaid on the
           CIGAR-deletion positions

All three panels share the x-axis (chr4_fip1l1_wt record-local
position), making the orthogonal-evidence story visually direct.

Output:
    sniffles_test/figure_deletion_validation.png
"""

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pysam


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
BAM_PATH = WORKING_DIR / "sniffles_test" / "deletion_aligned.sorted.bam"
COVERAGE_BED = WORKING_DIR / "sniffles_test" / "coverage" / "coverage.per-base.bed.gz"
OUTPUT_FIG = WORKING_DIR / "sniffles_test" / "figure_deletion_validation.png"

# Ground truth from build_deletion_reference.py
DELETION_START_LOCAL = 29_999  # 0-based, inclusive
DELETION_END_LOCAL = 35_000    # 0-based, exclusive


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

def load_coverage(bed_path: Path) -> tuple:
    """
    Read mosdepth's per-base BED.gz and expand its run-length encoding
    into per-position arrays.

    Returns (positions, coverages) as parallel lists. Each entry in
    positions is the chromosomal coordinate (1-based, inclusive); each
    entry in coverages is the depth at that position.
    """
    positions = []
    coverages = []

    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            chrom, start, end, depth = line.strip().split("\t")
            start = int(start)
            end = int(end)
            depth = int(depth)

            # mosdepth BED: 0-based start, exclusive end. Convert to a
            # range of 1-based inclusive positions for plotting.
            for pos in range(start + 1, end + 1):
                positions.append(pos)
                coverages.append(depth)

    return positions, coverages


def load_alignment_intervals(bam_path: Path) -> dict:
    """
    Pull alignment intervals from the BAM, organized by alignment type.

    Returns a dict with three keys:
      'primary':       list of (start, end, has_cigar_deletion)
      'supplementary': list of (start, end)
      'cigar_dels':    list of (start, end) where a large D operation
                       was found within a primary alignment
    """
    primaries = []
    supplementaries = []
    cigar_dels = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary:
                continue

            # Determine if this primary alignment carries a large D
            has_big_del = False
            del_start = None
            del_end = None

            if read.cigartuples and not read.is_supplementary:
                ref_pos = read.reference_start
                for op_code, op_length in read.cigartuples:
                    if op_code == 2 and op_length >= 1000:  # D >= 1 kb
                        has_big_del = True
                        del_start = ref_pos
                        del_end = ref_pos + op_length
                        cigar_dels.append((del_start, del_end))
                    if op_code in (0, 2, 3, 7, 8):
                        ref_pos += op_length

            entry = (
                read.reference_start,
                read.reference_end,
                has_big_del,
            )
            if read.is_supplementary:
                supplementaries.append((read.reference_start, read.reference_end))
            else:
                primaries.append(entry)

    return {
        "primary": primaries,
        "supplementary": supplementaries,
        "cigar_dels": cigar_dels,
    }


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def plot_panel_a_coverage(ax, positions, coverages):
    """Per-base coverage as a filled area chart."""
    ax.fill_between(positions, coverages, step="post", alpha=0.6, color="steelblue")
    ax.plot(positions, coverages, color="steelblue", linewidth=0.5)
    ax.set_ylabel("Read depth")
    ax.set_title("Panel A: Per-base coverage (mosdepth)")
    ax.grid(True, alpha=0.3)
    ax.axvspan(
        DELETION_START_LOCAL,
        DELETION_END_LOCAL,
        color="red",
        alpha=0.15,
        label="Ground-truth deletion",
    )
    ax.legend(loc="upper right", fontsize=8)


def plot_panel_b_alignments(ax, alignment_data):
    """Per-read alignment intervals as horizontal bars."""
    primaries = alignment_data["primary"]
    supplementaries = alignment_data["supplementary"]

    # Sort primary alignments by start position so the y-axis ordering
    # is deterministic and visually informative
    primaries_sorted = sorted(primaries, key=lambda x: x[0])

    # Create y-positions: each primary read gets a row, supplementaries
    # share row positions with their parent primary (we don't try to
    # match them up; we just stack supplementaries below the primaries)
    y_offset = 0
    for start, end, has_del in primaries_sorted:
        color = "tomato" if has_del else "lightsteelblue"
        ax.barh(y_offset, end - start, left=start, height=0.7,
                color=color, edgecolor="none")
        y_offset += 1

    # Plot supplementaries as a separate band below
    supp_y_start = y_offset + 5  # gap to separate
    for start, end in sorted(supplementaries, key=lambda x: x[0]):
        ax.barh(supp_y_start, end - start, left=start, height=0.7,
                color="darkorange", edgecolor="none")
        supp_y_start += 1

    # Annotate the deletion region
    ax.axvspan(
        DELETION_START_LOCAL,
        DELETION_END_LOCAL,
        color="red",
        alpha=0.15,
    )

    # Manual legend
    legend_handles = [
        mpatches.Patch(color="lightsteelblue", label="Primary alignment (no large indel)"),
        mpatches.Patch(color="tomato", label="Primary with CIGAR deletion ≥1 kb"),
        mpatches.Patch(color="darkorange", label="Supplementary alignment"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
    ax.set_ylabel("Read index")
    ax.set_title("Panel B: Per-read alignment intervals")
    ax.set_yticks([])  # individual read indices aren't meaningful


def plot_panel_c_breakpoints(ax, alignment_data):
    """
    CIGAR-deletion intervals plotted as horizontal bars at unique
    y-positions, with the consensus breakpoint annotated.
    """
    cigar_dels = alignment_data["cigar_dels"]
    cigar_dels_sorted = sorted(cigar_dels, key=lambda x: x[0])

    for i, (start, end) in enumerate(cigar_dels_sorted):
        ax.barh(i, end - start, left=start, height=0.7,
                color="tomato", edgecolor="darkred", linewidth=0.5)

    # Compute consensus breakpoint
    if cigar_dels:
        starts = [d[0] for d in cigar_dels]
        ends = [d[1] for d in cigar_dels]
        consensus_start = sum(starts) // len(starts)
        consensus_end = sum(ends) // len(ends)

        # Annotate consensus position as vertical lines
        ax.axvline(consensus_start, color="darkred", linestyle="--",
                   linewidth=1.5, label=f"Consensus 5' breakpoint: {consensus_start:,}")
        ax.axvline(consensus_end, color="darkred", linestyle="--",
                   linewidth=1.5, label=f"Consensus 3' breakpoint: {consensus_end:,}")

    # Annotate ground truth as solid vertical lines
    ax.axvline(DELETION_START_LOCAL, color="green", linestyle="-",
               linewidth=1.5, alpha=0.8,
               label=f"Ground-truth 5': {DELETION_START_LOCAL:,}")
    ax.axvline(DELETION_END_LOCAL, color="green", linestyle="-",
               linewidth=1.5, alpha=0.8,
               label=f"Ground-truth 3': {DELETION_END_LOCAL:,}")

    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel("chr4_fip1l1_wt position (record-local, 0-based)")
    ax.set_ylabel("Supporting read index")
    ax.set_title(f"Panel C: CIGAR-deletion evidence ({len(cigar_dels)} reads)")
    ax.set_yticks([])


def main() -> None:
    print("=" * 70)
    print("Generating validation figure")
    print("=" * 70)

    # Load all data
    print(f"\nLoading coverage track from {COVERAGE_BED.name}...")
    positions, coverages = load_coverage(COVERAGE_BED)
    print(f"  Loaded {len(positions):,} per-base entries")
    print(f"  Coverage range: {min(coverages)}x to {max(coverages)}x")

    print(f"\nLoading alignment data from {BAM_PATH.name}...")
    alignment_data = load_alignment_intervals(BAM_PATH)
    print(f"  Primary alignments: {len(alignment_data['primary'])}")
    print(f"    With large CIGAR deletion: {len(alignment_data['cigar_dels'])}")
    print(f"  Supplementary alignments: {len(alignment_data['supplementary'])}")

    # Build figure with three vertically stacked panels
    print(f"\nRendering figure...")
    fig, axes = plt.subplots(
        3, 1,
        figsize=(14, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 2, 1.5]},
    )

    plot_panel_a_coverage(axes[0], positions, coverages)
    plot_panel_b_alignments(axes[1], alignment_data)
    plot_panel_c_breakpoints(axes[2], alignment_data)

    fig.suptitle(
        "Synthetic deletion validation: 5 kb deletion in FIP1L1 intron 7\n"
        "Three orthogonal evidence channels confirming breakpoint coordinates",
        fontsize=12,
    )

    fig.tight_layout()
    fig.savefig(OUTPUT_FIG, dpi=150, bbox_inches="tight")

    print(f"\n  Figure saved: {OUTPUT_FIG}")
    print(f"  Resolution: 150 DPI, 14×10 in")


if __name__ == "__main__":
    main()