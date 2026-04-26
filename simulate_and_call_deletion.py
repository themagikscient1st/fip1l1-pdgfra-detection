"""
simulate_and_call_deletion.py

Simulate reads from the deletion-bearing reference, align to wild-type,
run Sniffles2. The expected output is a single DEL variant call at
chr4_fip1l1_wt:30,000-35,000 (record-local) corresponding to the 5 kb
synthetic deletion.

Outputs:
    sniffles_test/deletion_reads.fastq.gz
    sniffles_test/deletion_aligned.sorted.bam
    sniffles_test/sniffles_deletion.vcf
"""

import subprocess
from pathlib import Path


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
TEST_DIR = WORKING_DIR / "sniffles_test"

WILDTYPE_FASTA = TEST_DIR / "wildtype.fasta"
DELETION_FASTA = TEST_DIR / "deletion_bearing.fasta"

READS_FASTQ = TEST_DIR / "deletion_reads.fastq.gz"
READS_LOG = TEST_DIR / "deletion_reads.log"

ALIGNED_BAM = TEST_DIR / "deletion_aligned.sorted.bam"
ALIGN_LOG = TEST_DIR / "alignment.log"

VCF_OUTPUT = TEST_DIR / "sniffles_deletion.vcf"
SNIFFLES_LOG = TEST_DIR / "sniffles.log"


def run_simulation() -> None:
    """Simulate 30x coverage from the deletion-bearing reference."""
    print(f"\n--- Step 1: Simulating reads from deletion-bearing reference ---")

    cmd = [
        "badread", "simulate",
        "--reference", str(DELETION_FASTA),
        "--quantity", "30x",
        "--length", "15000,13000",
        "--error_model", "nanopore2023",
        "--qscore_model", "nanopore2023",
        "--chimeras", "0",
        "--junk_reads", "0",
        "--random_reads", "0",
    ]

    with open(READS_FASTQ, "wb") as out_handle, open(READS_LOG, "w") as log_handle:
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
        with open(READS_LOG) as f:
            raise RuntimeError(f"badread failed:\n{f.read()}")

    size_mb = READS_FASTQ.stat().st_size / 1024 / 1024
    print(f"  Reads written: {READS_FASTQ.name} ({size_mb:.2f} MB)")


def run_alignment() -> None:
    """Align deletion reads to wild-type reference."""
    print(f"\n--- Step 2: Aligning to wild-type reference ---")

    # Index the wild-type reference if not already indexed
    if not (WILDTYPE_FASTA.parent / (WILDTYPE_FASTA.name + ".fai")).exists():
        subprocess.run(["samtools", "faidx", str(WILDTYPE_FASTA)], check=True)
        print(f"  Indexed wild-type reference")

    minimap_cmd = [
        "minimap2",
        "-ax", "map-ont",
        "-t", "4",
        "-Y",
        "--MD",
        str(WILDTYPE_FASTA),
        str(READS_FASTQ),
    ]
    samtools_sort_cmd = [
        "samtools", "sort",
        "-O", "bam",
        "-@", "4",
        "-o", "-",
    ]

    with open(ALIGN_LOG, "w") as log_handle, open(ALIGNED_BAM, "wb") as out_handle:
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
        with open(ALIGN_LOG) as f:
            raise RuntimeError(f"Alignment failed:\n{f.read()}")

    subprocess.run(["samtools", "index", str(ALIGNED_BAM)], check=True)

    # Print summary
    flagstat = subprocess.run(
        ["samtools", "flagstat", str(ALIGNED_BAM)],
        capture_output=True, text=True, check=True,
    )
    print(f"  Alignment summary:")
    for line in flagstat.stdout.splitlines()[:5]:
        print(f"    {line}")


def run_sniffles() -> None:
    """Call SVs with Sniffles2."""
    print(f"\n--- Step 3: Calling SVs with Sniffles2 ---")

    cmd = [
        "sniffles",
        "--input", str(ALIGNED_BAM),
        "--reference", str(WILDTYPE_FASTA),
        "--vcf", str(VCF_OUTPUT),
        "--minsupport", "3",
        "--minsvlen", "30",
        "--allow-overwrite",
    ]

    with open(SNIFFLES_LOG, "w") as log_handle:
        result = subprocess.run(cmd, stderr=log_handle, stdout=log_handle)

    if result.returncode != 0:
        with open(SNIFFLES_LOG) as f:
            raise RuntimeError(f"Sniffles failed:\n{f.read()}")

    # Parse the VCF and report calls
    print(f"  --- SV calls ---")
    n_records = 0
    with open(VCF_OUTPUT) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            n_records += 1
            fields = line.strip().split("\t")
            chrom, pos, sv_id, ref, alt, qual, filt, info = fields[:8]

            info_dict = dict(
                f.split("=", 1) if "=" in f else (f, "True")
                for f in info.split(";")
            )
            svtype = info_dict.get("SVTYPE", "?")
            svlen = info_dict.get("SVLEN", "?")
            end = info_dict.get("END", "?")
            support = info_dict.get("SUPPORT", "?")

            print(f"\n  Variant {n_records}:")
            print(f"    Position: {chrom}:{pos}-{end}")
            print(f"    Type:     {svtype}")
            print(f"    Length:   {svlen} bp")
            print(f"    Quality:  {qual}")
            print(f"    Filter:   {filt}")
            print(f"    Support:  {support} reads")

    if n_records == 0:
        print(f"  No SVs called. Check {SNIFFLES_LOG}")
    else:
        print(f"\n  Total: {n_records} variant(s) called")


def main() -> None:
    print("=" * 70)
    print("Sniffles2 validation: 5 kb synthetic deletion in FIP1L1")
    print("=" * 70)

    if not WILDTYPE_FASTA.exists() or not DELETION_FASTA.exists():
        raise FileNotFoundError(
            "Wild-type and/or deletion-bearing FASTAs not found. "
            "Run build_deletion_reference.py first."
        )

    run_simulation()
    run_alignment()
    run_sniffles()

    print(f"\n{'=' * 70}")
    print(f"DONE")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()