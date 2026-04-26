"""
align_to_normal.py

Re-align simulated reads to ONLY the chr4_normal contig of the
synthetic reference. This mimics the clinical workflow: a fusion-
positive sample is aligned against the wild-type genome, and the
fusion is detected as a structural variant (a large deletion between
the FIP1L1 and PDGFRA breakpoints).

Aligning against chr4_normal alone removes the alignment ambiguity we
saw when both contigs were available simultaneously, and forces all
junction-spanning reads into split alignments across the deletion gap
- which is exactly the signal Sniffles is built to detect.

Output: aligned_to_normal/aligned_to_normal.sorted.bam
"""

import subprocess
from pathlib import Path
from Bio import SeqIO


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
SYNTHETIC_REFERENCE = WORKING_DIR / "synthetic_reference" / "synthetic_reference.fasta"
SIMULATED_READS = WORKING_DIR / "simulated_reads" / "simulated_reads.fastq.gz"

OUTPUT_DIR = WORKING_DIR / "aligned_to_normal"
NORMAL_ONLY_REFERENCE = OUTPUT_DIR / "chr4_normal.fasta"
SORTED_BAM = OUTPUT_DIR / "aligned_to_normal.sorted.bam"
LOG_FILE = OUTPUT_DIR / "alignment.log"


def extract_chr4_normal(input_fasta: Path, output_fasta: Path) -> None:
    """Pull just the chr4_normal record from the multi-record reference."""
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id == "chr4_normal":
            records.append(record)
    if not records:
        raise ValueError("chr4_normal not found in reference FASTA")
    with open(output_fasta, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    print(f"  Extracted chr4_normal: {len(records[0].seq):,} bp -> {output_fasta.name}")


def main() -> None:
    print("=" * 70)
    print("Re-aligning to chr4_normal only (clinical-workflow simulation)")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Extract chr4_normal as a single-record FASTA
    print(f"\nExtracting chr4_normal as single-contig reference:")
    extract_chr4_normal(SYNTHETIC_REFERENCE, NORMAL_ONLY_REFERENCE)

    # Step 2: Index the new reference
    print(f"\nIndexing the chr4_normal reference:")
    subprocess.run(["samtools", "faidx", str(NORMAL_ONLY_REFERENCE)], check=True)
    print(f"  Index created: {NORMAL_ONLY_REFERENCE}.fai")

    # Step 3: Align with minimap2 + sort with samtools
    print(f"\nAligning reads to chr4_normal:")
    minimap_cmd = [
        "minimap2",
        "-ax", "map-ont",
        "-t", "4",
        "-Y",
        "--MD",  # add MD tags (helps Sniffles refine breakpoints)
        str(NORMAL_ONLY_REFERENCE),
        str(SIMULATED_READS),
    ]
    samtools_sort_cmd = [
        "samtools", "sort",
        "-O", "bam",
        "-@", "4",
        "-o", "-",
    ]

    with open(LOG_FILE, "w") as log_handle, open(SORTED_BAM, "wb") as out_handle:
        minimap_proc = subprocess.Popen(
            minimap_cmd,
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        samtools_proc = subprocess.Popen(
            samtools_sort_cmd,
            stdin=minimap_proc.stdout,
            stdout=out_handle,
            stderr=log_handle,
        )
        if minimap_proc.stdout:
            minimap_proc.stdout.close()
        samtools_proc.wait()
        minimap_proc.wait()

    if minimap_proc.returncode != 0 or samtools_proc.returncode != 0:
        with open(LOG_FILE) as f:
            log_text = f.read()
        raise RuntimeError(
            f"Pipeline failed (minimap2={minimap_proc.returncode}, "
            f"samtools={samtools_proc.returncode}).\n--- log ---\n{log_text}"
        )

    # Step 4: Index the BAM
    print(f"\nIndexing BAM:")
    subprocess.run(["samtools", "index", str(SORTED_BAM)], check=True)

    # Step 5: Print summary
    print(f"\n  Alignment summary:")
    flagstat = subprocess.run(
        ["samtools", "flagstat", str(SORTED_BAM)],
        capture_output=True, text=True, check=True,
    )
    for line in flagstat.stdout.splitlines():
        print(f"    {line}")

    print(f"\n  Reads per contig:")
    idxstats = subprocess.run(
        ["samtools", "idxstats", str(SORTED_BAM)],
        capture_output=True, text=True, check=True,
    )
    for line in idxstats.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 4:
            contig, length, mapped, unmapped = parts[0], parts[1], parts[2], parts[3]
            print(f"    {contig:>15s}  length={int(length):>8,}  "
                  f"mapped={int(mapped):>4}  unmapped={int(unmapped):>4}")

    print(f"\n" + "=" * 70)
    print(f"DONE")
    print(f"=" * 70)
    print(f"BAM:     {SORTED_BAM}")
    print(f"BAM.bai: {SORTED_BAM}.bai")
    print(f"\nNext: re-run Sniffles against this BAM.")


if __name__ == "__main__":
    main()