"""
call_svs_v2.py

Run Sniffles2 against the chr4_normal-only alignment, where
junction-spanning reads have been forced into split alignments across
the deletion gap.

Output: variants/sniffles_chr4_normal.vcf
"""

import subprocess
from pathlib import Path


WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
NORMAL_ONLY_REFERENCE = WORKING_DIR / "aligned_to_normal" / "chr4_normal.fasta"
SORTED_BAM = WORKING_DIR / "aligned_to_normal" / "aligned_to_normal.sorted.bam"

OUTPUT_DIR = WORKING_DIR / "variants"
OUTPUT_VCF = OUTPUT_DIR / "sniffles_chr4_normal.vcf"
LOG_FILE = OUTPUT_DIR / "sniffles_chr4_normal.log"


def main() -> None:
    print("=" * 70)
    print("Calling SVs against chr4_normal-only alignment")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cmd = [
        "sniffles",
        "--input", str(SORTED_BAM),
        "--reference", str(NORMAL_ONLY_REFERENCE),
        "--vcf", str(OUTPUT_VCF),
        "--minsupport", "3",
        "--minsvlen", "30",
        "--allow-overwrite",
    ]

    print(f"\n  Command: {' '.join(cmd)}")

    with open(LOG_FILE, "w") as log_handle:
        result = subprocess.run(cmd, stderr=log_handle, stdout=log_handle)

    if result.returncode != 0:
        with open(LOG_FILE) as f:
            log_text = f.read()
        raise RuntimeError(f"Sniffles failed.\n--- log ---\n{log_text}")

    # Show called SVs
    print(f"\n  --- SV calls ---")
    n_records = 0
    with open(OUTPUT_VCF) as vcf:
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
            print(f"    Length:   {svlen}")
            print(f"    Quality:  {qual}")
            print(f"    Filter:   {filt}")
            print(f"    Support:  {support} reads")

    if n_records == 0:
        print(f"  No SVs called. Check {LOG_FILE} for details.")

    print(f"\nVCF:     {OUTPUT_VCF}")
    print(f"Log:     {LOG_FILE}")


if __name__ == "__main__":
    main()