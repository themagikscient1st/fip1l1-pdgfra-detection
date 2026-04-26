"""
build_fusion_figure.py

End-to-end FIP1L1::PDGFRA pipeline run for the clinical-target figure:
  1. Simulate reads from chr4_fusion at 30x with badread (chimeras off)
  2. Align to chr4_normal-only reference (the clinical-equivalent setup)
  3. Compute coverage with mosdepth
  4. Run breakpoint detection
  5. Render a three-panel figure analogous to figure_deletion_validation.png

The signal in this figure is the FIP1L1::PDGFRA junction. Reads from
the fusion contig spanning chr4_fusion position ~41,066 will produce
split alignments at the FIP1L1↔PDGFRA boundary in chr4_normal
(positions ~100,001 and ~101,002, separated by the 1 kb N-padded gap).

Output:
    fusion_figure/figure_fip1l1_pdgfra_target.png
"""

import gzip
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pysam
from Bio import SeqIO


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
SYNTHETIC_REFERENCE = WORKING_DIR / "synthetic_reference" / "synthetic_reference.fasta"

OUTPUT_DIR = WORKING_DIR / "fusion_figure"
CHR4_NORMAL_FASTA = OUTPUT_DIR / "chr4_normal.fasta"
CHR4_FUSION_FASTA = OUTPUT_DIR / "chr4_fusion.fasta"
READS_FASTQ = OUTPUT_DIR / "fusion_reads.fastq.gz"
READS_LOG = OUTPUT_DIR / "fusion_reads.log"
ALIGNED_BAM = OUTPUT_DIR / "fusion_aligned.sorted.bam"
ALIGN_LOG = OUTPUT_DIR / "alignment.log"
COVERAGE_DIR = OUTPUT_DIR / "coverage"
COVERAGE_PREFIX = COVERAGE_DIR / "coverage"
COVERAGE_BED = COVERAGE_DIR / "coverage.per-base.bed.gz"
OUTPUT_FIG = OUTPUT_DIR / "figure_fip1l1_pdgfra_target.png"

# Ground-truth coordinates in chr4_normal (181,002 bp contig)
FIP1L1_REGION_END = 100_001     # last base of FIP1L1 region in chr4_normal
PDGFRA_REGION_START = 101_002   # first base of PDGFRA region in chr4_normal

# Effective deletion coordinates from the chr4_fusion-vs-chr4_normal comparison.
# The chr4_fusion contig contains FIP1L1[1:41,066] joined directly to
# PDGFRA[~150,888:181,002]. When fusion-derived reads align to chr4_normal,
# the apparent deletion spans these positions.
EFFECTIVE_DEL_START = 41_066
EFFECTIVE_DEL_END = 150_888


# -----------------------------------------------------------------------------
# Pipeline steps
# -----------------------------------------------------------------------------

def split_reference() -> None:
    """Extract chr4_normal and chr4_fusion as single-contig FASTAs."""
    print(f"\n--- Step 1: extracting per-contig FASTAs ---")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for record in SeqIO.parse(SYNTHETIC_REFERENCE, "fasta"):
        if record.id == "chr4_normal":
            with open(CHR4_NORMAL_FASTA, "w") as fh:
                SeqIO.write([record], fh, "fasta")
            print(f"  chr4_normal: {len(record.seq):,} bp")
        elif record.id == "chr4_fusion":
            with open(CHR4_FUSION_FASTA, "w") as fh:
                SeqIO.write([record], fh, "fasta")
            print(f"  chr4_fusion: {len(record.seq):,} bp")

    subprocess.run(["samtools", "faidx", str(CHR4_NORMAL_FASTA)], check=True)


def simulate_fusion_reads() -> None:
    """Generate 30x reads from chr4_fusion with chimeras and junk disabled."""
    print(f"\n--- Step 2: simulating reads from chr4_fusion ---")

    cmd = [
        "badread", "simulate",
        "--reference", str(CHR4_FUSION_FASTA),
        "--quantity", "30x",
        "--length", "15000,13000",
        "--error_model", "nanopore2023",
        "--qscore_model", "nanopore2023",
        "--chimeras", "0",
        "--junk_reads", "0",
        "--random_reads", "0",
    ]

    with open(READS_FASTQ, "wb") as out_handle, open(READS_LOG, "w") as log_handle:
        gzip_proc = subprocess.Popen(["gzip", "-c"], stdin=subprocess.PIPE, stdout=out_handle)
        badread_proc = subprocess.Popen(cmd, stdout=gzip_proc.stdin, stderr=log_handle)
        badread_proc.wait()
        if gzip_proc.stdin:
            gzip_proc.stdin.close()
        gzip_proc.wait()

    if badread_proc.returncode != 0:
        with open(READS_LOG) as f:
            raise RuntimeError(f"badread failed:\n{f.read()}")

    size_mb = READS_FASTQ.stat().st_size / 1024 / 1024
    print(f"  Reads: {size_mb:.2f} MB")


def align_to_chr4_normal() -> None:
    """Align fusion reads to chr4_normal only."""
    print(f"\n--- Step 3: aligning to chr4_normal ---")

    minimap_cmd = [
        "minimap2", "-ax", "map-ont", "-t", "4", "-Y", "--MD",
        str(CHR4_NORMAL_FASTA), str(READS_FASTQ),
    ]
    samtools_cmd = ["samtools", "sort", "-O", "bam", "-@", "4", "-o", "-"]

    with open(ALIGN_LOG, "w") as log_handle, open(ALIGNED_BAM, "wb") as out_handle:
        minimap_proc = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=log_handle)
        samtools_proc = subprocess.Popen(samtools_cmd, stdin=minimap_proc.stdout,
                                         stdout=out_handle, stderr=log_handle)
        if minimap_proc.stdout:
            minimap_proc.stdout.close()
        samtools_proc.wait()
        minimap_proc.wait()

    if minimap_proc.returncode != 0 or samtools_proc.returncode != 0:
        with open(ALIGN_LOG) as f:
            raise RuntimeError(f"Alignment failed:\n{f.read()}")

    subprocess.run(["samtools", "index", str(ALIGNED_BAM)], check=True)

    # Print summary
    flagstat = subprocess.run(["samtools", "flagstat", str(ALIGNED_BAM)],
                              capture_output=True, text=True, check=True)
    print(f"  Alignment summary (first 5 lines):")
    for line in flagstat.stdout.splitlines()[:5]:
        print(f"    {line}")


def compute_coverage() -> None:
    """Run mosdepth for per-base coverage track."""
    print(f"\n--- Step 4: computing per-base coverage ---")

    COVERAGE_DIR.mkdir(parents=True, exist_ok=True)
    cmd = ["mosdepth", "--threads", "4", "--fast-mode",
           str(COVERAGE_PREFIX), str(ALIGNED_BAM)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mosdepth failed: {result.stderr}")

    print(f"  Coverage written: {COVERAGE_BED.name}")


# -----------------------------------------------------------------------------
# Plotting (reuses logic from plot_validation.py)
# -----------------------------------------------------------------------------

def load_coverage(bed_path: Path) -> tuple:
    positions, coverages = [], []
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            chrom, start, end, depth = line.strip().split("\t")
            start, end, depth = int(start), int(end), int(depth)
            for pos in range(start + 1, end + 1):
                positions.append(pos)
                coverages.append(depth)
    return positions, coverages


def load_alignment_intervals(bam_path: Path) -> dict:
    primaries = []
    supplementaries = []
    cigar_dels = []  # we won't have many of these for the fusion data,
                     # since the deletion is too large to encode as a CIGAR D op

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary:
                continue

            has_big_del = False
            if read.cigartuples and not read.is_supplementary:
                ref_pos = read.reference_start
                for op_code, op_length in read.cigartuples:
                    if op_code == 2 and op_length >= 1000:
                        has_big_del = True
                        cigar_dels.append((ref_pos, ref_pos + op_length))
                    if op_code in (0, 2, 3, 7, 8):
                        ref_pos += op_length

            if read.is_supplementary:
                supplementaries.append((read.reference_start, read.reference_end))
            else:
                primaries.append((read.reference_start, read.reference_end, has_big_del))

    return {"primary": primaries, "supplementary": supplementaries, "cigar_dels": cigar_dels}


def build_split_read_pairs(bam_path: Path) -> list:
    """For each read with a primary+supplementary alignment, return their endpoint pair."""
    primaries = {}
    supplementaries = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary:
                continue
            entry = (read.reference_start, read.reference_end)
            if read.is_supplementary:
                supplementaries.setdefault(read.query_name, []).append(entry)
            else:
                primaries[read.query_name] = entry

    pairs = []
    for read_name, prim in primaries.items():
        for supp in supplementaries.get(read_name, []):
            pairs.append((prim[0], prim[1], supp[0], supp[1], read_name))
    return pairs


def plot_panel_a_coverage(ax, positions, coverages):
    ax.fill_between(positions, coverages, step="post", alpha=0.6, color="steelblue")
    ax.plot(positions, coverages, color="steelblue", linewidth=0.5)
    ax.set_ylabel("Read depth")
    ax.set_title("Panel A: Per-base coverage on chr4_normal (mosdepth)")
    ax.grid(True, alpha=0.3)
    # Highlight the N-padded gap (the synthetic 'deletion' equivalent)
    ax.axvspan(
        EFFECTIVE_DEL_START, EFFECTIVE_DEL_END,
        color="red", alpha=0.15,
        label=f"Effective deletion (chr4_normal:{EFFECTIVE_DEL_START:,}-"
              f"{EFFECTIVE_DEL_END:,}; "
              f"~{(EFFECTIVE_DEL_END - EFFECTIVE_DEL_START) // 1000} kb)",
    )
    ax.legend(loc="upper right", fontsize=8)


def plot_panel_b_alignments(ax, alignment_data):
    primaries = sorted(alignment_data["primary"], key=lambda x: x[0])
    supplementaries = sorted(alignment_data["supplementary"], key=lambda x: x[0])

    y = 0
    for start, end, has_del in primaries:
        color = "tomato" if has_del else "lightsteelblue"
        ax.barh(y, end - start, left=start, height=0.7, color=color, edgecolor="none")
        y += 1

    supp_y = y + 5
    for start, end in supplementaries:
        ax.barh(supp_y, end - start, left=start, height=0.7, color="darkorange", edgecolor="none")
        supp_y += 1

    ax.axvspan(EFFECTIVE_DEL_START, EFFECTIVE_DEL_END, color="red", alpha=0.15)

    legend_handles = [
        mpatches.Patch(color="lightsteelblue", label="Primary alignment"),
        mpatches.Patch(color="darkorange", label="Supplementary alignment"),
    ]
    if any(p[2] for p in primaries):
        legend_handles.insert(1, mpatches.Patch(color="tomato",
                                                 label="Primary with CIGAR del ≥1 kb"))

    ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
    ax.set_ylabel("Read index")
    ax.set_title("Panel B: Per-read alignment intervals (chr4_normal)")
    ax.set_yticks([])


def plot_panel_c_split_reads(ax, split_pairs):
    """
    Plot each junction-spanning read as a pair of bars (primary on top,
    supplementary connected by a thin dashed line) to make the
    split-read structure visible.
    """
    # Sort by primary start position so the y-axis is deterministic
    pairs_sorted = sorted(split_pairs, key=lambda x: x[0])

    for i, (prim_start, prim_end, supp_start, supp_end, _name) in enumerate(pairs_sorted):
        # Primary alignment (lighter)
        ax.barh(i, prim_end - prim_start, left=prim_start, height=0.7,
                color="lightsteelblue", edgecolor="none")
        # Supplementary alignment (orange, same row)
        ax.barh(i, supp_end - supp_start, left=supp_start, height=0.7,
                color="darkorange", edgecolor="none")
        # Connecting line between them (dashed, indicates the split)
        ax.plot([min(prim_end, supp_end), max(prim_start, supp_start)],
                [i, i], color="gray", linewidth=0.5, linestyle=":", alpha=0.6)

    ax.axvspan(EFFECTIVE_DEL_START, EFFECTIVE_DEL_END, color="red", alpha=0.15)

    legend_handles = [
        mpatches.Patch(color="lightsteelblue", label="Primary alignment piece"),
        mpatches.Patch(color="darkorange", label="Supplementary alignment piece"),
        mpatches.Patch(color="red", alpha=0.3,
                       label=f"Effective deletion (~{(EFFECTIVE_DEL_END - EFFECTIVE_DEL_START) // 1000} kb)"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
    ax.set_xlabel("chr4_normal position (record-local, 0-based)")
    ax.set_ylabel("Junction-spanning read index")
    ax.set_title(f"Panel C: Split-read evidence ({len(split_pairs)} junction-spanning reads)")
    ax.set_yticks([])


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def render_figure() -> None:
    print(f"\n--- Step 5: rendering figure ---")

    positions, coverages = load_coverage(COVERAGE_BED)
    print(f"  Loaded coverage: {len(positions):,} positions, "
          f"range {min(coverages)}x to {max(coverages)}x")

    alignment_data = load_alignment_intervals(ALIGNED_BAM)
    print(f"  Primary alignments: {len(alignment_data['primary'])}")
    print(f"  Supplementary alignments: {len(alignment_data['supplementary'])}")

    split_pairs = build_split_read_pairs(ALIGNED_BAM)
    print(f"  Split-read pairs: {len(split_pairs)}")

    fig, axes = plt.subplots(
        3, 1, figsize=(14, 10), sharex=True,
        gridspec_kw={"height_ratios": [1, 2, 1.5]},
    )

    plot_panel_a_coverage(axes[0], positions, coverages)
    plot_panel_b_alignments(axes[1], alignment_data)
    plot_panel_c_split_reads(axes[2], split_pairs)

    fig.suptitle(
        "FIP1L1::PDGFRA fusion detection: synthetic clinical-target validation\n"
        "Reads simulated from chr4_fusion, aligned to chr4_normal (wild-type-equivalent)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG, dpi=150, bbox_inches="tight")
    print(f"\n  Figure saved: {OUTPUT_FIG}")


def main() -> None:
    print("=" * 70)
    print("Building FIP1L1::PDGFRA clinical-target figure")
    print("=" * 70)

    #split_reference()
    #simulate_fusion_reads()
    #align_to_chr4_normal()
    #compute_coverage()
    render_figure()

    print(f"\n{'=' * 70}")
    print(f"DONE")
    print(f"{'=' * 70}")
    print(f"\nOutput: {OUTPUT_FIG}")


if __name__ == "__main__":
    main()