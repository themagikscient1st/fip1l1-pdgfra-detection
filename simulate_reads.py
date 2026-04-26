"""
simulate_reads.py

Simulate nanopore reads from the synthetic FIP1L1::PDGFRA reference
using badread. Outputs a single FASTQ file containing reads from both
the wild-type-equivalent and fusion contigs at a defined ratio.

Configuration:
  - 15x total coverage, 1:1 fusion:normal ratio (~7.5x each contig)
  - Badread "nanopore2023" preset (R10.4.1-equivalent error profile)
  - Reads from both contigs are merged into a single FASTQ output

Output: simulated_reads/simulated_reads.fastq.gz
"""

import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
SYNTHETIC_REFERENCE = WORKING_DIR / "synthetic_reference" / "synthetic_reference.fasta"

OUTPUT_DIR = WORKING_DIR / "simulated_reads"
TEMP_DIR = OUTPUT_DIR / "tmp"
FINAL_OUTPUT = OUTPUT_DIR / "simulated_reads.fastq.gz"

# Per-contig coverage. Half the total 15x goes to each contig.
COVERAGE_PER_CONTIG = "30x"

# Badread chemistry preset
BADREAD_PRESET = "nanopore2023"

# Read length distribution (badread's defaults are sensible for nanopore;
# we override to make reads more uniform around 5-20 kb, which is realistic
# for adaptive-sampling output)
MEAN_LENGTH = "15000"
LENGTH_STDEV = "13000"


# -----------------------------------------------------------------------------
# Helper: split the multi-record reference into per-contig files
# -----------------------------------------------------------------------------

def split_reference_by_contig(reference: Path, output_dir: Path) -> dict:
    """
    Badread takes a single-contig reference at a time. Split our two-record
    synthetic reference into per-contig FASTAs.

    Returns a dict mapping contig name -> path to its single-contig FASTA.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    contig_files = {}

    for record in SeqIO.parse(reference, "fasta"):
        contig_path = output_dir / f"{record.id}.fasta"
        # Write a single-record FASTA
        with open(contig_path, "w") as fh:
            SeqIO.write([record], fh, "fasta")
        contig_files[record.id] = contig_path
        print(f"  Split {record.id}: {len(record.seq):,} bp -> {contig_path.name}")

    return contig_files


# -----------------------------------------------------------------------------
# Helper: run badread for one contig
# -----------------------------------------------------------------------------

def run_badread(
    contig_fasta: Path,
    output_fastq: Path,
    coverage: str,
    preset: str,
    mean_length: str,
    length_stdev: str,
) -> None:
    """
    Run badread on a single contig and write reads to a gzipped FASTQ.

    Stderr is written to a sibling .log file rather than captured in a
    pipe, which avoids OS-level pipe-buffer deadlocks when badread
    produces a lot of progress output.
    """
    print(f"\n  Simulating reads for {contig_fasta.stem}...")
    print(f"    Coverage:    {coverage}")
    print(f"    Preset:      {preset}")
    print(f"    Length:      mean {mean_length} bp (sd {length_stdev})")

    cmd = [
        "badread", "simulate",
        "--reference", str(contig_fasta),
        "--quantity", coverage,
        "--length", f"{mean_length},{length_stdev}",
        "--error_model", preset,
        "--qscore_model", preset,
        "--chimeras", "0",       # Disable chimeric reads (they create false-positive SVs)
        "--junk_reads", "0",      # Disable junk reads (they were inflating the unmapped count)
        "--random_reads", "0",    # Same idea
    ]

    log_file = output_fastq.with_suffix(".log")

    with open(output_fastq, "wb") as out_handle, open(log_file, "w") as log_handle:
        gzip_proc = subprocess.Popen(
            ["gzip", "-c"],
            stdin=subprocess.PIPE,
            stdout=out_handle,
        )
        badread_proc = subprocess.Popen(
            cmd,
            stdout=gzip_proc.stdin,
            stderr=log_handle,
        )
        badread_proc.wait()
        if gzip_proc.stdin:
            gzip_proc.stdin.close()
        gzip_proc.wait()

        if badread_proc.returncode != 0:
            log_handle.flush()
            with open(log_file) as f:
                error_text = f.read()
            raise RuntimeError(
                f"badread failed (exit code {badread_proc.returncode}) for "
                f"{contig_fasta.name}.\n--- badread stderr ---\n{error_text}"
            )

    with open(log_file) as f:
        for line in f:
            if any(key in line for key in ["coverage", "yield", "reads", "length"]):
                print(f"    {line.strip()}")

    size_bytes = output_fastq.stat().st_size
    if size_bytes < 100:
        raise RuntimeError(
            f"badread produced a near-empty output ({size_bytes} bytes) "
            f"for {contig_fasta.name}. Check {log_file} for details."
        )
    print(f"    Output: {size_bytes:,} bytes written to {output_fastq.name}")


# -----------------------------------------------------------------------------
# Helper: concatenate multiple gzipped FASTQs into one
# -----------------------------------------------------------------------------

def concatenate_fastqs(input_files: list, output_file: Path) -> None:
    """
    Merge per-contig FASTQ files into a single combined output.

    Gzipped FASTQs can be concatenated directly at the byte level (gzip
    is designed to support this), so we just cat them together.
    """
    print(f"\n  Concatenating {len(input_files)} FASTQs into {output_file.name}...")
    with open(output_file, "wb") as out_handle:
        for input_file in input_files:
            with open(input_file, "rb") as in_handle:
                out_handle.write(in_handle.read())
    size_mb = output_file.stat().st_size / 1024 / 1024
    print(f"    Combined output: {size_mb:.2f} MB")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("Simulating nanopore reads from synthetic FIP1L1::PDGFRA reference")
    print("=" * 70)

    # Verify the reference exists
    if not SYNTHETIC_REFERENCE.exists():
        print(f"ERROR: Synthetic reference not found at {SYNTHETIC_REFERENCE}")
        print(f"Run build_synthetic_reference.py first.")
        return

    # Set up output directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    # Split the reference into per-contig files
    print(f"\nSplitting reference into per-contig files:")
    contig_files = split_reference_by_contig(SYNTHETIC_REFERENCE, TEMP_DIR)

    # Run badread on each contig
    per_contig_fastqs = []
    for contig_name, contig_fasta in contig_files.items():
        per_contig_fastq = TEMP_DIR / f"{contig_name}.fastq.gz"
        run_badread(
            contig_fasta=contig_fasta,
            output_fastq=per_contig_fastq,
            coverage=COVERAGE_PER_CONTIG,
            preset=BADREAD_PRESET,
            mean_length=MEAN_LENGTH,
            length_stdev=LENGTH_STDEV,
        )
        per_contig_fastqs.append(per_contig_fastq)

    # Merge into final output
    concatenate_fastqs(per_contig_fastqs, FINAL_OUTPUT)

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Output: {FINAL_OUTPUT}")
    print(f"Per-contig outputs preserved in: {TEMP_DIR}")
    print(f"\nNext step: align these reads to GRCh38 and run SV calling.")


if __name__ == "__main__":
    main()