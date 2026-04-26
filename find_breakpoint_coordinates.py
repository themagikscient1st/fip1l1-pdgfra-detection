"""
find_breakpoint_coordinates.py

Locate the Type A breakpoint coordinates (FIP1L1 exon 10 / PDGFRA exon 12
position 44) in our cropped GenBank files, anchored on Walz et al. (2009)
Supplementary Table 4 sample 16. Sample 16 is selected because:
    - It carries the most common Type A configuration
      (FIP1L1 exon 10 / PDGFRA exon 12 position 44)
    - The DNA junction is precisely nucleotide-resolved
      (no microhomology between the two genes at the breakpoint)
    - The full DNA junction sequence is available

Strategy:
    1. For FIP1L1: search for sample 16's intronic anchor sequence in
       FIP1L1.gb. The position immediately after the anchor is the
       FIP1L1-side breakpoint.
    2. For PDGFRA: locate exon 12 of the canonical PDGFRA mRNA
       (NM_006206.6) via the mRNA feature, then position 44 of exon 12
       is exon_12_start + 43.
    3. Sanity-check by confirming the sequence at position 44 matches
       the Walz-published sequence (CCCGGA).

Output: chromosomal coordinates for both breakpoints, plus diagnostic
context.
"""

from pathlib import Path
from Bio import SeqIO


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

WORKING_DIR = Path.home() / "miniconda3" / "FIP1L1_PDGFRA FASTAs"
FIP1L1_FILE = WORKING_DIR / "FIP1L1.gb"
PDGFRA_FILE = WORKING_DIR / "PDGFRA.gb"

# Chromosomal start of each cropped region on NC_000004.12.
# Used to convert record-local coordinates to chromosomal ones.
FIP1L1_OFFSET = 53_370_000
PDGFRA_OFFSET = 54_225_000

# Walz et al. (2009) Supplementary Table 4, sample 16:
#     FIP1L1 exon 10 / PDGFRA exon 12 position 44, Type A
#     mRNA junction: AGCAGCCCGGA
#     DNA junction:  TTTTTGGGTGTTTCTCATCAGCCCGGATGGACATGAATATA
# Split:
#     FIP1L1 intronic side: TTTTTGGGTGTTTCTCATCAG (21 bp)
#     PDGFRA exon 12 side:  CCCGGATGGACATGAATATA (starts at exon 12 pos 44)
WALZ_FIP1L1_INTRONIC_ANCHOR = "TTTTTGGGTGTTTCTCATCAG"
WALZ_PDGFRA_EXON12_POS44_SEQ = "CCCGGA"

# Canonical PDGFRA mRNA RefSeq accession
PDGFRA_CANONICAL_MRNA = "NM_006206.6"


# -----------------------------------------------------------------------------
# FIP1L1 lookup
# -----------------------------------------------------------------------------

def find_fip1l1_breakpoint() -> dict:
    """
    Find the FIP1L1-side breakpoint by searching for the Walz sample 16
    intronic anchor in FIP1L1.gb.

    Returns a dictionary with the breakpoint coordinate and diagnostic info,
    or an empty dictionary if the lookup fails.
    """
    print(f"\n--- FIP1L1 breakpoint lookup ---")
    print(f"  Source:  Walz et al. (2009) Supplementary Table 4, sample 16")
    print(f"  Anchor:  {WALZ_FIP1L1_INTRONIC_ANCHOR}")

    record = SeqIO.read(FIP1L1_FILE, "genbank")
    sequence = str(record.seq).upper()

    # Search the forward strand
    hit_position = sequence.find(WALZ_FIP1L1_INTRONIC_ANCHOR)

    if hit_position == -1:
        # Not found on forward strand; try reverse complement.
        # If the anchor is on the reverse strand, our breakpoint analysis
        # gets more complicated, so we just report the situation and stop.
        from Bio.Seq import Seq
        rc_anchor = str(Seq(WALZ_FIP1L1_INTRONIC_ANCHOR).reverse_complement())
        rc_hit = sequence.find(rc_anchor)
        if rc_hit != -1:
            print(f"  ERROR: anchor not found on forward strand, but reverse")
            print(f"  complement was found at position {rc_hit:,}.")
            print(f"  This shouldn't happen — FIP1L1 is on the forward strand")
            print(f"  per the Gene record. Investigate before proceeding.")
        else:
            print(f"  ERROR: anchor not found on either strand of FIP1L1.gb.")
            print(f"  The cropped region (chr4:{FIP1L1_OFFSET+1:,}-"
                  f"{FIP1L1_OFFSET+len(sequence):,}) may not contain")
            print(f"  intron 10, or the anchor sequence may not be in the")
            print(f"  reference genome (sample 16's intron 10 may differ).")
        return {}

    # Verify uniqueness — if the anchor appears more than once, we have an
    # ambiguous lookup and need to investigate.
    second_hit = sequence.find(WALZ_FIP1L1_INTRONIC_ANCHOR, hit_position + 1)
    if second_hit != -1:
        print(f"  WARNING: anchor sequence appears more than once in the file.")
        print(f"    First occurrence:  position {hit_position:,}")
        print(f"    Second occurrence: position {second_hit:,}")
        print(f"  This is unexpected for a 21 bp sequence. Investigate.")
        return {}

    # The breakpoint sits at the first base AFTER the FIP1L1 portion.
    # In our DNA junction, that's the 'C' of CCCGGA — but in the genome,
    # it's the first base of intron sequence that gets deleted. Position
    # in our genomic file = end of the anchor.
    anchor_length = len(WALZ_FIP1L1_INTRONIC_ANCHOR)
    breakpoint_local = hit_position + anchor_length  # 0-based; first deleted base
    breakpoint_chrom = breakpoint_local + FIP1L1_OFFSET + 1  # 1-based

    print(f"  Found at:        position {hit_position:,} (0-based, record-local)")
    print(f"                   chr4:{hit_position + FIP1L1_OFFSET + 1:,} (1-based, chromosomal)")
    print(f"  Breakpoint:      chr4:{breakpoint_chrom:,}")
    print(f"                   (first base of the deleted interval on the FIP1L1 side)")
    print(f"\n  Sequence context (50 bp each side of breakpoint):")
    print(f"    Upstream:   ...{sequence[breakpoint_local-50:breakpoint_local]}|")
    print(f"    Downstream: |{sequence[breakpoint_local:breakpoint_local+50]}...")
    print(f"  (The | marks the breakpoint position.)")

    return {
        "breakpoint_chrom": breakpoint_chrom,
        "breakpoint_local": breakpoint_local,
        "anchor_position": hit_position,
    }


# -----------------------------------------------------------------------------
# PDGFRA lookup
# -----------------------------------------------------------------------------

def find_pdgfra_breakpoint() -> dict:
    """
    Find the PDGFRA-side breakpoint via two-step lookup:
      1. Locate the canonical PDGFRA mRNA in PDGFRA.gb feature list
      2. Pull exon 12 from its CompoundLocation parts
      3. Position 44 of exon 12 = exon_12_start + 43 (1-based to 0-based)
      4. Sanity-check that the sequence at position 44 matches the
         Walz-published 'CCCGGA' anchor.
    """
    print(f"\n--- PDGFRA breakpoint lookup ---")
    print(f"  Strategy: locate exon 12 of {PDGFRA_CANONICAL_MRNA},")
    print(f"            then position 44 within it (Walz exon 12 numbering)")

    record = SeqIO.read(PDGFRA_FILE, "genbank")
    sequence = str(record.seq).upper()

    # Find the canonical mRNA feature
    canonical_mrna = None
    for feature in record.features:
        if feature.type != "mRNA":
            continue
        transcript_id = feature.qualifiers.get("transcript_id", [""])[0]
        if transcript_id == PDGFRA_CANONICAL_MRNA:
            canonical_mrna = feature
            break

    if canonical_mrna is None:
        print(f"  ERROR: {PDGFRA_CANONICAL_MRNA} not found in PDGFRA.gb.")
        print(f"  Available PDGFRA mRNA transcripts:")
        for feature in record.features:
            if (feature.type == "mRNA"
                    and feature.qualifiers.get("gene", [""])[0] == "PDGFRA"):
                tid = feature.qualifiers.get("transcript_id", ["?"])[0]
                n_exons = len(feature.location.parts)
                print(f"    {tid}: {n_exons} exons")
        return {}

    # mRNA features have a CompoundLocation whose parts list the exons in
    # transcription order (5' to 3').
    parts = canonical_mrna.location.parts
    print(f"  Found {PDGFRA_CANONICAL_MRNA} with {len(parts)} exons.")

    if len(parts) < 12:
        print(f"  ERROR: only {len(parts)} exons present; need exon 12.")
        return {}

    exon_12 = parts[11]  # 0-based index, so parts[11] is the 12th exon
    exon_12_start_local = int(exon_12.start)  # 0-based, inclusive
    exon_12_end_local = int(exon_12.end)      # 0-based, exclusive
    exon_12_length = exon_12_end_local - exon_12_start_local

    print(f"  Exon 12 record-local: positions {exon_12_start_local:,}-"
          f"{exon_12_end_local:,} (0-based)")
    print(f"  Exon 12 length:       {exon_12_length:,} bp")
    print(f"  Exon 12 chromosomal:  chr4:{exon_12_start_local + PDGFRA_OFFSET + 1:,}-"
          f"{exon_12_end_local + PDGFRA_OFFSET:,}")

    # Position 44 within exon 12 (1-based) is exon_12_start + 43 (0-based)
    breakpoint_local = exon_12_start_local + 43
    breakpoint_chrom = breakpoint_local + PDGFRA_OFFSET + 1

    # Sanity check: confirm the sequence at position 44 matches Walz's CCCGGA
    observed = sequence[breakpoint_local:breakpoint_local + len(WALZ_PDGFRA_EXON12_POS44_SEQ)]

    print(f"\n  Position 44 of exon 12: chr4:{breakpoint_chrom:,}")
    print(f"  Sequence at position 44:  {observed}")
    print(f"  Walz-published expected:  {WALZ_PDGFRA_EXON12_POS44_SEQ}")
    if observed == WALZ_PDGFRA_EXON12_POS44_SEQ:
        print(f"  --> MATCH: position 44 confirmed by sequence.")
    else:
        print(f"  --> MISMATCH: investigate before using this coordinate.")
        return {}

    print(f"\n  Sequence context (50 bp each side of position 44):")
    print(f"    Upstream:   ...{sequence[breakpoint_local-50:breakpoint_local]}|")
    print(f"    Downstream: |{sequence[breakpoint_local:breakpoint_local+50]}...")

    return {
        "breakpoint_chrom": breakpoint_chrom,
        "breakpoint_local": breakpoint_local,
        "exon_12_start_chrom": exon_12_start_local + PDGFRA_OFFSET + 1,
        "exon_12_end_chrom": exon_12_end_local + PDGFRA_OFFSET,
    }


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

def main() -> None:
    print("=" * 70)
    print("Locating Type A FIP1L1::PDGFRA breakpoint coordinates")
    print("Anchored on Walz et al. (2009) Supp. Table 4, sample 16")
    print("Configuration: FIP1L1 exon 10 / PDGFRA exon 12 position 44")
    print("=" * 70)

    fip1l1_result = find_fip1l1_breakpoint()
    pdgfra_result = find_pdgfra_breakpoint()

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if not (fip1l1_result and pdgfra_result):
        print("One or both lookups failed; see diagnostic output above.")
        return

    fip1l1_bp = fip1l1_result['breakpoint_chrom']
    pdgfra_bp = pdgfra_result['breakpoint_chrom']
    deletion_size = pdgfra_bp - fip1l1_bp

    print(f"FIP1L1 5' breakpoint:    chr4:{fip1l1_bp:,}")
    print(f"PDGFRA 3' breakpoint:    chr4:{pdgfra_bp:,}")
    print(f"Implied deletion size:   {deletion_size:,} bp ({deletion_size/1000:.1f} kb)")
    print(f"\n(Walz et al. report ~800 kb canonical deletion size; ")
    print(f" our derived figure should be in the 700-900 kb ballpark.)")


if __name__ == "__main__":
    main()