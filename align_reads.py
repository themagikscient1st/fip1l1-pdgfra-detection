"""
align_reads.py

Align simulated nanopore reads back to the synthetic reference using
minimap2, then sort and index the output BAM with samtools.

We align to the synthetic reference (rather than to GRCh38) because:
    1. It's faster — no need to load a 3 GB human genome index
    2. It's analytically equivalent for testing — reads simulated from
       chr4_normal and chr4_fusion will map to those contigs, and
       junction-spanning reads on chr4_fusion will produce the
       split-alignment signature that SV callers detect
    3. Ground-truth verification is direct — we know exactly which
       reads should map to which contigs and where the junction is

Output: aligned_reads/aligned_reads.sorted.bam (+ .bai index)
"""

import subprocess
from pathlib import Path


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
SYNTHETIC_REFERENCE = WORKING_DIR / "synthetic_reference" / "synthetic_reference.fasta"
SIMULATED_READS = WORKING_DIR / "simulated_reads" / "simulated_reads.fastq.gz"

OUTPUT_DIR = WORKING_DIR / "aligned_reads"
SORTED_BAM = OUTPUT_DIR / "aligned_reads.sorted.bam"
LOG_FILE = OUTPUT_DIR / "alignment.log"


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def run_command(cmd: list, log_handle, description: str) -> None:
    """
    Run a shell command, logging stderr to a file and raising on failure.

    Stdout is allowed to flow through (which is needed for piped commands
    where one tool's stdout is another tool's stdin), but stderr goes to
    a log file so we can inspect it after the fact without flooding the
    terminal.
    """
    print(f"\n  {description}")
    print(f"    Command: {' '.join(cmd)}")
    log_handle.write(f"\n=== {description} ===\n")
    log_handle.write(f"Command: {' '.join(cmd)}\n")
    log_handle.flush()

    result = subprocess.run(cmd, stderr=log_handle, check=False)
    if result.returncode != 0:
        log_handle.flush()
        with open(LOG_FILE) as f:
            log_text = f.read()
        raise RuntimeError(
            f"Command failed (exit code {result.returncode}): {' '.join(cmd)}\n"
            f"--- log file contents ---\n{log_text}"
        )


def run_pipeline(
    minimap_cmd: list,
    samtools_cmd: list,
    output_path: Path,
    log_handle,
    description: str,
) -> None:
    """
    Run a piped command pipeline: minimap2 | samtools sort > output.bam

    minimap2 outputs SAM to stdout; samtools sort reads SAM from stdin
    and writes sorted BAM to stdout (which we redirect to the output file).
    """
    print(f"\n  {description}")
    print(f"    Pipeline: {' '.join(minimap_cmd)} | {' '.join(samtools_cmd)} > {output_path.name}")
    log_handle.write(f"\n=== {description} ===\n")
    log_handle.write(f"Pipeline: {' '.join(minimap_cmd)} | {' '.join(samtools_cmd)} > {output_path}\n")
    log_handle.flush()

    with open(output_path, "wb") as out_handle:
        minimap_proc = subprocess.Popen(
            minimap_cmd,
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        samtools_proc = subprocess.Popen(
            samtools_cmd,
            stdin=minimap_proc.stdout,
            stdout=out_handle,
            stderr=log_handle,
        )
        # Allow minimap to receive SIGPIPE if samtools exits early
        if minimap_proc.stdout:
            minimap_proc.stdout.close()
        samtools_proc.wait()
        minimap_proc.wait()

    if minimap_proc.returncode != 0 or samtools_proc.returncode != 0:
        log_handle.flush()
        with open(LOG_FILE) as f:
            log_text = f.read()
        raise RuntimeError(
            f"Pipeline failed (minimap2={minimap_proc.returncode}, "
            f"samtools={samtools_proc.returncode}).\n"
            f"--- log file contents ---\n{log_text}"
        )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("Aligning simulated reads to synthetic reference")
    print("=" * 70)

    # Verify inputs
    if not SYNTHETIC_REFERENCE.exists():
        raise FileNotFoundError(f"Reference not found: {SYNTHETIC_REFERENCE}")
    if not SIMULATED_READS.exists():
        raise FileNotFoundError(f"Reads not found: {SIMULATED_READS}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(LOG_FILE, "w") as log_handle:

        # ---------------------------------------------------------------------
        # Step 1: Align reads with minimap2 and sort with samtools, in a pipeline
        # ---------------------------------------------------------------------
        # minimap2 flags:
        #   -a               output SAM format (default is PAF)
        #   -x map-ont       preset for Oxford Nanopore long reads
        #   -t 4             use 4 threads (tune for your system)
        #   -Y               use soft-clipping for supplementary alignments
        #                    (preserves clipped sequence in BAM, important for SV calling)
        # samtools sort flags:
        #   -O bam           output BAM format
        #   -o -             write to stdout (which we redirect to the output file)
        minimap_cmd = [
            "minimap2",
            "-ax", "map-ont",
            "-t", "4",
            "-Y",
            str(SYNTHETIC_REFERENCE),
            str(SIMULATED_READS),
        ]
        samtools_sort_cmd = [
            "samtools", "sort",
            "-O", "bam",
            "-@", "4",
            "-o", "-",
        ]
        run_pipeline(
            minimap_cmd,
            samtools_sort_cmd,
            SORTED_BAM,
            log_handle,
            description="Aligning reads with minimap2 and sorting with samtools",
        )

        # ---------------------------------------------------------------------
        # Step 2: Index the sorted BAM
        # ---------------------------------------------------------------------
        run_command(
            ["samtools", "index", str(SORTED_BAM)],
            log_handle,
            description="Indexing the sorted BAM",
        )

        # ---------------------------------------------------------------------
        # Step 3: Print alignment summary
        # ---------------------------------------------------------------------
        print(f"\n  Alignment summary:")
        log_handle.write(f"\n=== Alignment summary (samtools flagstat) ===\n")
        log_handle.flush()

        flagstat = subprocess.run(
            ["samtools", "flagstat", str(SORTED_BAM)],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in flagstat.stdout.splitlines():
            print(f"    {line}")
            log_handle.write(line + "\n")

        # ---------------------------------------------------------------------
        # Step 4: Per-contig read counts
        # ---------------------------------------------------------------------
        print(f"\n  Reads per contig:")
        log_handle.write(f"\n=== Reads per contig (samtools idxstats) ===\n")
        log_handle.flush()

        idxstats = subprocess.run(
            ["samtools", "idxstats", str(SORTED_BAM)],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in idxstats.stdout.strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                contig, length, mapped, unmapped = parts[0], parts[1], parts[2], parts[3]
                print(f"    {contig:>15s}  length={int(length):>8,}  "
                      f"mapped={int(mapped):>4}  unmapped={int(unmapped):>4}")
                log_handle.write(line + "\n")

    print(f"\n" + "=" * 70)
    print(f"DONE")
    print(f"=" * 70)
    print(f"Sorted BAM: {SORTED_BAM}")
    print(f"BAM index:  {SORTED_BAM}.bai")
    print(f"Log file:   {LOG_FILE}")
    print(f"\nNext step: SV calling with Sniffles2.")


if __name__ == "__main__":
    main()