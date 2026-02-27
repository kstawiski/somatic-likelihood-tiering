#!/usr/bin/env python3
"""Somatic Likelihood Tiering (SLT) — 4-tier classification for tumor-only WES.

Classifies somatic variants from tumor-only whole-exome sequencing using
four orthogonal evidence layers (POPAF, gnomAD, GERMQ, COSMIC) and PureCN
posterior somatic probabilities, with integrated CHIP detection.

All thresholds frozen from development cohort. No post-hoc optimization.

Usage:
    python slt_classify.py --input variants.tsv --output slt_tiers.tsv
    python slt_classify.py --input mc3_annotated.maf --output mc3_slt.tsv --degraded
"""

import argparse
import csv
import math
import sys
from typing import Optional


# =============================================================================
# PureCN Posterior Thresholds (Riester et al., 2016 recommended levels)
# =============================================================================
POSTERIOR_A = 0.8   # SLT-A eligible gate
POSTERIOR_B = 0.5   # SLT-B eligible gate
POSTERIOR_C = 0.2   # SLT-C eligible gate

# =============================================================================
# Evidence Layer Thresholds (frozen from development cohort)
# =============================================================================
POPAF_THRESHOLD = 5.0          # Layer 1: Phred-scaled pop AF from Mutect2
GNOMAD_AF_THRESHOLD = 0.001    # Layer 2: gnomAD exome AF
GERMQ_THRESHOLD = 30           # Layer 3: Mutect2 germline quality score
COSMIC_CONFIRMED_MIN = 5       # Layer 4a: COSMIC confirmed somatic count
COSMIC_HOTSPOT_MIN = 10        # Layer 4b: COSMIC hotspot count
COSMIC_CGC_SAMPLE_MIN = 2      # Layer 4c: CGC gene + COSMIC sample count

# =============================================================================
# CHIP Gene Lists (Jaiswal et al., 2014; Genovese et al., 2014; Bolton et al., 2020)
# =============================================================================
CHIP_TIER1_GENES = frozenset([
    "DNMT3A", "TET2", "ASXL1", "PPM1D", "JAK2", "SF3B1", "SRSF2",
    "U2AF1", "IDH1", "IDH2", "ZBTB33", "GNB1", "CBL",
])

CHIP_TIER2_GENES = frozenset([
    "TP53", "GNAS", "STAT3", "BCOR", "BCORL1", "CALR", "CEBPA", "CSF3R",
    "CUX1", "DDX41", "ETV6", "EZH2", "ETNK1", "FLT3", "KIT", "KRAS",
    "MYD88", "NF1", "NOTCH1", "NPM1", "NRAS", "PHF6", "RAD21", "RUNX1",
    "SETBP1", "STAG2", "SUZ12", "ZRSR2",
])

# CHIP Hotspot Patterns
CHIP_TIER1_HOTSPOTS = {
    "DNMT3A": ["R882"],
    "JAK2": ["V617"],
    "SF3B1": ["K700", "K666"],
    "SRSF2": ["P95"],
    "U2AF1": ["S34", "Q157"],
    "IDH1": ["R132"],
    "IDH2": ["R140", "R172"],
    "CBL": ["Y371", "C384"],
    "GNB1": ["K57"],
}

CHIP_TIER2_HOTSPOTS = {
    "TP53": ["R175", "R248", "R273"],
    "KRAS": ["G12", "G13"],
}

# Truncating variant consequences (LOF)
TRUNCATING_CONSEQUENCES = frozenset([
    "nonsense", "frameshift_del", "frameshift_ins", "splice_site",
    "translation_start_site", "nonstop_mutation",
    "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
])

# CHIP VAF Thresholds
CHIP_TIER1_VAF = 0.15
CHIP_TP53_VAF = 0.15
CHIP_KRAS_VAF = 0.10


def safe_float(value, default=None):
    """Parse a float from string, returning default for missing/invalid values."""
    if value is None or value == "" or value == "." or value == "NA" or value == "nan":
        return default
    try:
        f = float(value)
        return default if math.isnan(f) else f
    except (ValueError, TypeError):
        return default


def safe_int(value, default=None):
    """Parse an int from string, returning default for missing/invalid."""
    if value is None or value == "" or value == "." or value == "NA":
        return default
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return default


# =============================================================================
# Evidence Layer Evaluation
# =============================================================================

def evaluate_layer1_popaf(popaf: Optional[float]) -> bool:
    """Layer 1: POPAF ≥ 5.0 (somatic-supporting)."""
    if popaf is None:
        return False
    return popaf >= POPAF_THRESHOLD


def evaluate_layer2_gnomad(gnomad_af: Optional[float]) -> bool:
    """Layer 2: gnomAD_AF < 0.001 OR missing (somatic-supporting)."""
    if gnomad_af is None:
        return True  # Missing gnomAD_AF is somatic-supporting by design
    return gnomad_af < GNOMAD_AF_THRESHOLD


def evaluate_layer3_germq(germq: Optional[float]) -> bool:
    """Layer 3: GERMQ ≥ 30 (somatic-supporting)."""
    if germq is None:
        return False
    return germq >= GERMQ_THRESHOLD


def evaluate_layer4_cosmic(
    cosmic_confirmed: Optional[int] = None,
    cosmic_hotspot: Optional[int] = None,
    is_cgc_gene: bool = False,
    cosmic_sample: Optional[int] = None,
) -> bool:
    """Layer 4: COSMIC composite criterion (somatic-supporting)."""
    if cosmic_confirmed is not None and cosmic_confirmed >= COSMIC_CONFIRMED_MIN:
        return True
    if cosmic_hotspot is not None and cosmic_hotspot >= COSMIC_HOTSPOT_MIN:
        return True
    if is_cgc_gene and cosmic_sample is not None and cosmic_sample >= COSMIC_CGC_SAMPLE_MIN:
        return True
    return False


def count_somatic_layers(layer1: bool, layer2: bool, layer3: bool, layer4: bool) -> int:
    """Count active somatic-supporting evidence layers."""
    return sum([layer1, layer2, layer3, layer4])


# =============================================================================
# Evidence Level Assignment
# =============================================================================

def assign_evidence_level(n_somatic_layers: int, layer4_cosmic: bool) -> str:
    """Assign evidence level: high / medium / low."""
    if n_somatic_layers >= 3 and layer4_cosmic:
        return "high"
    elif n_somatic_layers >= 2:
        return "medium"
    else:
        return "low"


# =============================================================================
# CHIP Classification
# =============================================================================

def is_hotspot_match(gene: str, protein_change: str, hotspot_dict: dict) -> bool:
    """Check if a variant matches a CHIP hotspot pattern."""
    if gene not in hotspot_dict or not protein_change:
        return False
    for pattern in hotspot_dict[gene]:
        if protein_change.startswith(pattern) or protein_change.startswith("p." + pattern):
            return True
    return False


def is_truncating(consequence: str) -> bool:
    """Check if variant consequence is truncating (LOF)."""
    if not consequence:
        return False
    return any(c in TRUNCATING_CONSEQUENCES for c in consequence.split(","))


def classify_chip(
    gene: str,
    vaf: Optional[float],
    protein_change: str,
    consequence: str,
) -> str:
    """CHIP 3-state classification per SLT classification rules.

    Returns: 'chip_likely', 'chip_possible', or 'no_chip'
    """
    if not gene:
        return "no_chip"

    is_tier1 = gene in CHIP_TIER1_GENES
    is_tier2 = gene in CHIP_TIER2_GENES

    if not is_tier1 and not is_tier2:
        return "no_chip"

    # LOF genes: TET2, ASXL1, PPM1D — any truncating = possible_CH
    lof_genes = {"TET2", "ASXL1", "PPM1D"}

    if is_tier1:
        is_hotspot = is_hotspot_match(gene, protein_change, CHIP_TIER1_HOTSPOTS)
        is_trunc = is_truncating(consequence)
        is_lof_gene = gene in lof_genes

        # possible_CH: Tier1 AND (VAF < 0.15 OR hotspot/truncating)
        if vaf is not None and vaf < CHIP_TIER1_VAF:
            return "chip_likely"
        if is_hotspot:
            return "chip_likely"
        if is_lof_gene and is_trunc:
            return "chip_likely"
        # Tier1 gene but conditions not met
        return "chip_possible"

    if is_tier2:
        # TP53 hotspot AND VAF < 0.15
        if gene == "TP53":
            is_hotspot = is_hotspot_match(gene, protein_change, CHIP_TIER2_HOTSPOTS)
            if is_hotspot and vaf is not None and vaf < CHIP_TP53_VAF:
                return "chip_likely"

        # KRAS G12/G13 AND VAF < 0.10
        if gene == "KRAS":
            is_hotspot = is_hotspot_match(gene, protein_change, CHIP_TIER2_HOTSPOTS)
            if is_hotspot and vaf is not None and vaf < CHIP_KRAS_VAF:
                return "chip_likely"

        # Tier2 gene but conditions not met
        return "chip_possible"

    return "no_chip"


# =============================================================================
# SLT Tier Assignment Cascade
# =============================================================================

def assign_slt_tier(
    posterior: Optional[float],
    evidence_level: str,
    n_somatic_layers: int,
    chip_status: str,
) -> str:
    """Assign SLT tier via deterministic cascade.

    Evaluated in order; first matching condition wins.
    """
    # Handle missing posterior: NaN/None → 0.0
    p = posterior if posterior is not None else 0.0

    # Handle missing chip_status
    cs = chip_status if chip_status else "no_chip"

    # SLT-A: posterior ≥ 0.8 AND evidence ∈ {high, medium} AND chip ≠ chip_likely
    if p >= POSTERIOR_A and evidence_level in ("high", "medium") and cs != "chip_likely":
        return "SLT-A"

    # SLT-B: posterior ≥ 0.5 OR (evidence = high AND n_layers ≥ 3)
    if p >= POSTERIOR_B:
        return "SLT-B"
    if evidence_level == "high" and n_somatic_layers >= 3:
        return "SLT-B"

    # SLT-C: posterior ≥ 0.2 OR evidence = medium
    if p >= POSTERIOR_C:
        return "SLT-C"
    if evidence_level == "medium":
        return "SLT-C"

    # SLT-D: everything else
    return "SLT-D"


# =============================================================================
# Full SLT Pipeline
# =============================================================================

def classify_variant(row: dict, degraded: bool = False) -> dict:
    """Apply full SLT classification to a single variant.

    Args:
        row: dict with variant fields
        degraded: if True, disable layers 1/3 and posterior (MAF-level analysis)

    Returns:
        dict with SLT classification fields added
    """
    # Parse input fields
    posterior = None if degraded else safe_float(row.get("POSTERIOR.SOMATIC") or row.get("POSTERIOR_SOMATIC"))
    popaf = None if degraded else safe_float(row.get("POPAF"))
    gnomad_af = safe_float(row.get("gnomAD_AF") or row.get("gnomAD_exome_AF"))
    germq = None if degraded else safe_float(row.get("GERMQ"))
    cosmic_confirmed = safe_int(row.get("COSMIC_CONFIRMED_SOMATIC"))
    cosmic_hotspot = safe_int(row.get("COSMIC_HOTSPOT"))
    cosmic_sample = safe_int(row.get("COSMIC_SAMPLE") or row.get("COSMIC_TOTAL_OCC"))
    is_cgc = str(row.get("CGC_GENE", "")).upper() in ("TRUE", "YES", "1")
    gene = row.get("Hugo_Symbol") or row.get("GENE") or row.get("gene") or ""
    vaf = safe_float(row.get("t_alt_freq") or row.get("AF") or row.get("VAF"))
    protein_change = row.get("HGVSp_Short") or row.get("Protein_Change") or row.get("AAChange") or ""
    consequence = row.get("Variant_Classification") or row.get("Consequence") or ""

    # Evaluate evidence layers
    l1 = evaluate_layer1_popaf(popaf)
    l2 = evaluate_layer2_gnomad(gnomad_af)
    l3 = evaluate_layer3_germq(germq)
    l4 = evaluate_layer4_cosmic(cosmic_confirmed, cosmic_hotspot, is_cgc, cosmic_sample)

    n_layers = count_somatic_layers(l1, l2, l3, l4)
    ev_level = assign_evidence_level(n_layers, l4)

    # CHIP classification
    chip = classify_chip(gene, vaf, protein_change, consequence)

    # Tier assignment
    tier = assign_slt_tier(posterior, ev_level, n_layers, chip)

    # Build output
    result = dict(row)
    result["slt_layer1_popaf"] = l1
    result["slt_layer2_gnomad"] = l2
    result["slt_layer3_germq"] = l3
    result["slt_layer4_cosmic"] = l4
    result["slt_n_somatic_layers"] = n_layers
    result["slt_evidence_level"] = ev_level
    result["slt_chip_status"] = chip
    result["slt_posterior_used"] = posterior if posterior is not None else "NA"
    result["slt_tier"] = tier
    result["slt_mode"] = "degraded" if degraded else "full"

    return result


def process_file(input_path: str, output_path: str, degraded: bool = False):
    """Process a TSV/MAF file and add SLT classifications."""
    with open(input_path, "r") as fin:
        # Skip comment lines
        header_line = None
        for line in fin:
            if not line.startswith("#"):
                header_line = line
                break

        if not header_line:
            print("ERROR: No header found", file=sys.stderr)
            sys.exit(1)

        reader = csv.DictReader(fin, fieldnames=header_line.strip().split("\t"), delimiter="\t")

        # Determine output fields
        first_row = next(reader, None)
        if first_row is None:
            print("ERROR: No data rows", file=sys.stderr)
            sys.exit(1)

        classified = classify_variant(first_row, degraded)
        out_fields = list(classified.keys())

        with open(output_path, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerow(classified)

            n = 1
            tier_counts = {classified["slt_tier"]: 1}

            for row in reader:
                c = classify_variant(row, degraded)
                writer.writerow(c)
                n += 1
                tier_counts[c["slt_tier"]] = tier_counts.get(c["slt_tier"], 0) + 1

    # Summary
    mode_label = "DEGRADED (no PureCN, 2/4 layers)" if degraded else "FULL (4 layers)"
    print(f"SLT Classification Complete [{mode_label}]")
    print(f"  Total variants: {n}")
    for tier in ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]:
        count = tier_counts.get(tier, 0)
        pct = 100 * count / n if n > 0 else 0
        print(f"  {tier}: {count} ({pct:.1f}%)")
    print(f"  Output: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="SLT: Somatic Likelihood Tiering for tumor-only WES",
        epilog="See README.md for full documentation",
    )
    parser.add_argument("--input", "-i", required=True, help="Input TSV/MAF file")
    parser.add_argument("--output", "-o", required=True, help="Output TSV with SLT tiers")
    parser.add_argument("--degraded", action="store_true",
                        help="Degraded mode: no PureCN posterior, no POPAF/GERMQ (MAF-level)")
    args = parser.parse_args()

    process_file(args.input, args.output, args.degraded)


if __name__ == "__main__":
    main()
