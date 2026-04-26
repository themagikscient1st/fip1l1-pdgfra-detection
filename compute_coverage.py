"""
compute_coverage.py

Run mosdepth on the deletion-validation BAM to produce per-base
coverage tracks. The output BED file is the input to the visualization
script.

Output:
    sniffles_test/coverage/coverage.per-base.bed.gz
"""

import subprocess
from pathlib import Path


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
BAM_PATH = WORKING_DIR / "sniffles_test" / "deletion_aligned.sorted.bam"

OUTPUT_DIR = WORKING_DIR / "sniffles_test" / "coverage"
OUTPUT_PREFIX = OUTPUT_DIR / "coverage"


def main() -> None:
    print("=" * 70)
    print("Computing per-base coverage with mosdepth")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # mosdepth flags:
    #   --by 1                  per-base resolution (rather than windowed)
    #   --no-per-base           DON'T disable per-base output
    #                           (we explicitly want it; default is windowed)
    #   --threads 4             parallelize
    #   --fast-mode             skip computing some per-base statistics we don't need
    #                           (significantly faster, output is still per-base coverage)
    #   <prefix>                output files share this prefix
    #   <bam>                   input BAM
    cmd = [
        "mosdepth",
        "--threads", "4",
        "--fast-mode",
        str(OUTPUT_PREFIX),
        str(BAM_PATH),
    ]

    print(f"\n  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"\n  STDOUT: {result.stdout}")
        print(f"  STDERR: {result.stderr}")
        raise RuntimeError(f"mosdepth failed with exit code {result.returncode}")

    # mosdepth writes several output files; the per-base track is the
    # gzipped BED file with .per-base.bed.gz suffix
    per_base_bed = OUTPUT_DIR / "coverage.per-base.bed.gz"
    if not per_base_bed.exists():
        raise FileNotFoundError(
            f"Expected output not found: {per_base_bed}\n"
            f"mosdepth wrote: {list(OUTPUT_DIR.iterdir())}"
        )

    size_kb = per_base_bed.stat().st_size / 1024
    print(f"\n  Per-base coverage written: {per_base_bed.name} ({size_kb:.1f} KB)")
    print(f"  Other mosdepth outputs in {OUTPUT_DIR}:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        print(f"    {f.name}")


if __name__ == "__main__":
    main()