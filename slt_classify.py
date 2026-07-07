#!/usr/bin/env python3
"""R1 SLT classifier with callability-aware gnomAD Layer 2.

This revision module preserves the submitted SLT cascade for full Mutect2/PureCN
inputs while fixing the reviewer-identified annotation-only contradiction and
preventing missing gnomAD annotation from being treated as rarity evidence.

`classify_variant` expects pre-annotated gnomAD AF/AN/state fields. The allele
matching helpers in this module are small utilities and regression-test fixtures
for the upstream annotation contract; they are not a replacement for cohort-scale
reference FASTA-aware normalization before final metric generation.
"""

import argparse
import csv
import math
import sys
from typing import Optional


POSTERIOR_A = 0.8
POSTERIOR_B = 0.5
POSTERIOR_C = 0.2

POPAF_THRESHOLD = 5.0
GNOMAD_AF_THRESHOLD = 0.001
GNOMAD_MIN_AN = 10000
GERMQ_THRESHOLD = 30
COSMIC_CONFIRMED_MIN = 5
COSMIC_HOTSPOT_MIN = 10
COSMIC_CGC_SAMPLE_MIN = 2

GNOMAD_RARE_CALLABLE = "rare_callable"
GNOMAD_COMMON = "common"
GNOMAD_UNEVALUABLE = "unevaluable"
GNOMAD_STATES = {GNOMAD_RARE_CALLABLE, GNOMAD_COMMON, GNOMAD_UNEVALUABLE}
POSTERIOR_FIELDS = ("POSTERIOR.SOMATIC", "POSTERIOR_SOMATIC")

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

TRUNCATING_CONSEQUENCES = frozenset([
    "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
    "nonsense", "frameshift_del", "frameshift_ins", "splice_site",
    "translation_start_site", "nonstop_mutation",
    "stop_gained", "frameshift_variant",
    "splice_donor_variant", "splice_acceptor_variant",
    "start_lost", "stop_lost",
])

CHIP_TIER1_VAF = 0.15
CHIP_TP53_VAF = 0.15
CHIP_KRAS_VAF = 0.10


def safe_float(value, default=None):
    if value is None or value == "" or value == "." or value == "NA" or value == "nan":
        return default
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return default if math.isnan(parsed) else parsed


def safe_int(value, default=None):
    if value is None or value == "" or value == "." or value == "NA":
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def safe_bool(value) -> Optional[bool]:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "callable", "evaluable"}:
        return True
    if text in {"0", "false", "f", "no", "n", "uncallable", "unevaluable"}:
        return False
    return None


def first_present(row, names):
    for name in names:
        value = row.get(name)
        if value not in (None, "", ".", "NA"):
            return value
    return None


def normalize_chrom(chrom) -> str:
    text = str(chrom).strip()
    return text[3:] if text.lower().startswith("chr") else text


def normalize_variant_key(chrom, pos, ref, alt) -> tuple[str, str, str, str]:
    """Return a deterministic minimal key for allele-specific matching.

    This trims shared sequence where possible without requiring a reference
    genome. It does not replace full left-alignment in repetitive indel contexts.
    """
    chrom_norm = normalize_chrom(chrom)
    pos_int = safe_int(pos)
    if pos_int is None:
        raise ValueError(f"Invalid variant position: {pos!r}")
    ref_norm = str(ref).strip().upper()
    alt_norm = str(alt).strip().upper()

    while len(ref_norm) > 1 and len(alt_norm) > 1 and ref_norm[-1] == alt_norm[-1]:
        ref_norm = ref_norm[:-1]
        alt_norm = alt_norm[:-1]

    while len(ref_norm) > 1 and len(alt_norm) > 1 and ref_norm[0] == alt_norm[0]:
        pos_int += 1
        ref_norm = ref_norm[1:]
        alt_norm = alt_norm[1:]

    return chrom_norm, str(pos_int), ref_norm, alt_norm


def iter_normalized_variant_keys(chrom, pos, ref, alt_field):
    for alt in str(alt_field).split(","):
        alt = alt.strip()
        if alt:
            yield normalize_variant_key(chrom, pos, ref, alt)


def variant_keys_intersect(candidate, gnomad_record) -> bool:
    candidate_keys = set(iter_normalized_variant_keys(*candidate))
    gnomad_keys = set(iter_normalized_variant_keys(*gnomad_record))
    return bool(candidate_keys & gnomad_keys)


def evaluate_layer1_popaf(popaf: Optional[float]) -> bool:
    return popaf is not None and popaf >= POPAF_THRESHOLD


def derive_gnomad_state(row: dict, min_an: int = GNOMAD_MIN_AN) -> tuple[str, str]:
    """Derive three-state gnomAD evidence without conflating absence and missingness."""
    af = safe_float(first_present(row, ["gnomAD_AF", "gnomAD_exome_AF", "gnomad_af"]))
    an = safe_int(first_present(row, ["gnomAD_AN", "gnomAD_exome_AN", "gnomad_an"]))

    explicit = first_present(row, ["gnomAD_state", "gnomad_state"])
    explicit_state = None
    if explicit:
        state = str(explicit).strip().lower().replace("-", "_")
        if state in GNOMAD_STATES:
            explicit_state = state

    if af is not None and af >= GNOMAD_AF_THRESHOLD:
        if explicit_state in {GNOMAD_RARE_CALLABLE, GNOMAD_UNEVALUABLE}:
            return GNOMAD_UNEVALUABLE, "conflicting_gnomad_state_af_an"
        return GNOMAD_COMMON, "af_at_or_above_threshold"

    if af is not None and af < GNOMAD_AF_THRESHOLD and an is not None and an >= min_an:
        if explicit_state in {GNOMAD_COMMON, GNOMAD_UNEVALUABLE}:
            return GNOMAD_UNEVALUABLE, "conflicting_gnomad_state_af_an"
        return GNOMAD_RARE_CALLABLE, "rare_af_with_adequate_an"

    if explicit_state == GNOMAD_RARE_CALLABLE:
        if an is not None and an < min_an:
            return GNOMAD_UNEVALUABLE, "conflicting_gnomad_state_af_an"
        return GNOMAD_RARE_CALLABLE, "explicit_state"

    if explicit_state == GNOMAD_COMMON:
        return GNOMAD_UNEVALUABLE, "missing_observed_common_af"

    if explicit_state == GNOMAD_UNEVALUABLE:
        return GNOMAD_UNEVALUABLE, "explicit_state"

    if af is None and safe_bool(first_present(row, [
        "gnomAD_callable", "gnomAD_evaluable", "gnomad_callable", "gnomad_evaluable",
    ])) is True:
        return GNOMAD_UNEVALUABLE, "missing_observed_af_for_sites_vcf"

    return GNOMAD_UNEVALUABLE, "missing_or_insufficient_callability"


def evaluate_layer2_gnomad_state(gnomad_state: str) -> bool:
    return gnomad_state == GNOMAD_RARE_CALLABLE


def evaluate_layer3_germq(germq: Optional[float]) -> bool:
    return germq is not None and germq >= GERMQ_THRESHOLD


def evaluate_layer4_cosmic(
    cosmic_confirmed: Optional[int] = None,
    cosmic_hotspot: Optional[int] = None,
    is_cgc_gene: bool = False,
    cosmic_sample: Optional[int] = None,
) -> bool:
    if cosmic_confirmed is not None and cosmic_confirmed >= COSMIC_CONFIRMED_MIN:
        return True
    if cosmic_hotspot is not None and cosmic_hotspot >= COSMIC_HOTSPOT_MIN:
        return True
    if is_cgc_gene and cosmic_sample is not None and cosmic_sample >= COSMIC_CGC_SAMPLE_MIN:
        return True
    return False


def count_somatic_layers(layer1: bool, layer2: bool, layer3: bool, layer4: bool) -> int:
    return sum([layer1, layer2, layer3, layer4])


def assign_evidence_level(n_somatic_layers: int, layer4_cosmic: bool) -> str:
    if n_somatic_layers >= 3 and layer4_cosmic:
        return "high"
    if n_somatic_layers >= 2:
        return "medium"
    return "low"


def is_hotspot_match(gene: str, protein_change: str, hotspot_dict: dict) -> bool:
    if gene not in hotspot_dict or not protein_change:
        return False
    pc = protein_change[2:] if protein_change.startswith("p.") else protein_change
    for pattern in hotspot_dict[gene]:
        if not pc.startswith(pattern):
            continue
        remainder = pc[len(pattern):]
        if remainder == "" or not remainder[0].isdigit():
            return True
    return False


def is_truncating(consequence: str) -> bool:
    if not consequence:
        return False
    import re
    terms = re.split(r"[,&|]", consequence)
    return any(t.strip() in TRUNCATING_CONSEQUENCES for t in terms)


def classify_chip(
    gene: str,
    vaf: Optional[float],
    protein_change: str,
    consequence: str,
) -> str:
    if not gene:
        return "no_chip"

    is_tier1 = gene in CHIP_TIER1_GENES
    is_tier2 = gene in CHIP_TIER2_GENES
    if not is_tier1 and not is_tier2:
        return "no_chip"

    if is_tier1:
        if vaf is not None and vaf < CHIP_TIER1_VAF:
            return "chip_likely"
        if is_hotspot_match(gene, protein_change, CHIP_TIER1_HOTSPOTS):
            return "chip_likely"
        if gene in {"TET2", "ASXL1", "PPM1D"} and is_truncating(consequence):
            return "chip_likely"
        return "chip_possible"

    if gene == "TP53":
        is_hotspot = is_hotspot_match(gene, protein_change, CHIP_TIER2_HOTSPOTS)
        if is_hotspot and vaf is not None and vaf < CHIP_TP53_VAF:
            return "chip_likely"
    if gene == "KRAS":
        is_hotspot = is_hotspot_match(gene, protein_change, CHIP_TIER2_HOTSPOTS)
        if is_hotspot and vaf is not None and vaf < CHIP_KRAS_VAF:
            return "chip_likely"
    return "chip_possible"


def assign_slt_tier(
    posterior: Optional[float],
    evidence_level: str,
    chip_status: str,
) -> str:
    """Assign full-mode SLT tier; high evidence already encodes n_layers >= 3."""
    p = 0.0 if posterior is None else posterior
    cs = chip_status if chip_status else "no_chip"

    if p >= POSTERIOR_A and evidence_level in ("high", "medium") and cs != "chip_likely":
        return "SLT-A"

    if cs != "chip_likely":
        if p >= POSTERIOR_B:
            return "SLT-B"
        if evidence_level == "high":
            return "SLT-B"

    if p >= POSTERIOR_C:
        return "SLT-C"
    if evidence_level in ("high", "medium"):
        return "SLT-C"
    return "SLT-D"


def assign_annotation_only_tier(gnomad_state: str, layer4_cosmic: bool) -> tuple[str, str]:
    """Six-cell R1 annotation-only matrix.

    Only confidently common, COSMIC-negative variants route to SLT-D. All other
    annotation-only states are retained at capped SLT-C as cannot-rule-out cases.
    """
    if gnomad_state == GNOMAD_COMMON and not layer4_cosmic:
        return "SLT-D", "common_cosmic_negative"
    return "SLT-C", f"{gnomad_state}_{'cosmic_positive' if layer4_cosmic else 'cosmic_negative'}"


def classify_variant(row: dict, mode: str = "full") -> dict:
    if mode not in {"full", "annotation_only"}:
        raise ValueError(f"Unsupported SLT mode: {mode}")

    annotation_only = mode == "annotation_only"
    if not annotation_only and not any(field in row for field in POSTERIOR_FIELDS):
        raise ValueError(
            "Full SLT mode requires PureCN somatic probability column "
            "(POSTERIOR.SOMATIC or POSTERIOR_SOMATIC) to be present in the schema; "
            "use annotation_only mode when posterior support is unavailable."
        )

    posterior = None if annotation_only else safe_float(first_present(row, POSTERIOR_FIELDS))
    posterior_status = "not_applicable" if annotation_only else (
        "present" if posterior is not None else "absent"
    )
    popaf = None if annotation_only else safe_float(row.get("POPAF"))
    germq = None if annotation_only else safe_float(row.get("GERMQ"))
    gnomad_state, gnomad_reason = derive_gnomad_state(row)

    cosmic_confirmed = safe_int(row.get("COSMIC_CONFIRMED_SOMATIC"))
    cosmic_hotspot = safe_int(row.get("COSMIC_HOTSPOT"))
    cosmic_sample = safe_int(row.get("COSMIC_SAMPLE") or row.get("COSMIC_TOTAL_OCC"))
    is_cgc = str(row.get("CGC_GENE", "")).upper() in ("TRUE", "YES", "1")
    gene = row.get("Hugo_Symbol") or row.get("GENE") or row.get("gene") or ""
    vaf = safe_float(row.get("t_alt_freq") or row.get("AF") or row.get("VAF"))
    protein_change = row.get("HGVSp_Short") or row.get("Protein_Change") or row.get("AAChange") or ""
    consequence = row.get("Variant_Classification") or row.get("Consequence") or ""

    l1 = evaluate_layer1_popaf(popaf)
    l2 = evaluate_layer2_gnomad_state(gnomad_state)
    l3 = evaluate_layer3_germq(germq)
    l4 = evaluate_layer4_cosmic(cosmic_confirmed, cosmic_hotspot, is_cgc, cosmic_sample)
    n_layers = count_somatic_layers(l1, l2, l3, l4)
    ev_level = assign_evidence_level(n_layers, l4)
    chip = classify_chip(gene, vaf, protein_change, consequence)

    matrix_cell = ""
    if annotation_only:
        tier, matrix_cell = assign_annotation_only_tier(gnomad_state, l4)
    else:
        tier = assign_slt_tier(posterior, ev_level, chip)

    result = dict(row)
    result["slt_layer1_popaf"] = l1
    result["slt_gnomad_state"] = gnomad_state
    result["slt_gnomad_state_reason"] = gnomad_reason
    result["slt_layer2_gnomad"] = l2
    result["slt_layer3_germq"] = l3
    result["slt_layer4_cosmic"] = l4
    result["slt_n_somatic_layers"] = n_layers
    result["slt_evidence_level"] = ev_level
    result["slt_chip_status"] = chip
    result["slt_posterior_status"] = posterior_status
    result["slt_posterior_used"] = posterior if posterior is not None else "NA"
    result["slt_annotation_only_matrix_cell"] = matrix_cell
    result["slt_tier"] = tier
    result["slt_mode"] = mode
    result["slt_schema_version"] = "R1-callability-posterior-v2"
    result["slt_layer2_semantics"] = "gnomAD rare_callable only"
    return result


def process_file(
    input_path: str,
    output_path: str,
    mode: str = "full",
    allow_unevaluable_gnomad: bool = False,
):
    with open(input_path, "r") as fin:
        header_line = None
        for line in fin:
            if not line.startswith("#"):
                header_line = line
                break
        if not header_line:
            print("ERROR: No header found", file=sys.stderr)
            sys.exit(1)

        header_fields = header_line.strip().split("\t")
        reader = csv.DictReader(fin, fieldnames=header_fields, delimiter="\t")
        first_row = next(reader, None)
        if first_row is None:
            print("ERROR: No data rows", file=sys.stderr)
            sys.exit(1)
        raw_rows = [first_row] + list(reader)

        if mode == "full":
            if not any(field in header_fields for field in POSTERIOR_FIELDS):
                print(
                    "ERROR: Full SLT mode requires PureCN somatic probability column "
                    "(POSTERIOR.SOMATIC or POSTERIOR_SOMATIC) in the input schema.",
                    file=sys.stderr,
                )
                sys.exit(2)
            posterior_present = sum(
                safe_float(first_present(row, POSTERIOR_FIELDS)) is not None
                for row in raw_rows
            )
            if posterior_present == 0:
                print(
                    "ERROR: Full SLT mode requires at least one nonblank PureCN somatic "
                    "probability value in the cohort input. An all-blank posterior column "
                    "does not establish a PureCN-backed workflow; use annotation_only mode "
                    "or provide validated field mapping.",
                    file=sys.stderr,
                )
                sys.exit(2)
        else:
            posterior_present = 0

        try:
            classified_rows = [classify_variant(row, mode=mode) for row in raw_rows]
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(2)

    n = len(classified_rows)
    tier_counts = {}
    state_counts = {}
    reason_counts = {}
    for row in classified_rows:
        tier_counts[row["slt_tier"]] = tier_counts.get(row["slt_tier"], 0) + 1
        state_counts[row["slt_gnomad_state"]] = state_counts.get(row["slt_gnomad_state"], 0) + 1
        reason = row["slt_gnomad_state_reason"]
        reason_counts[reason] = reason_counts.get(reason, 0) + 1

    if mode == "full" and state_counts.get(GNOMAD_UNEVALUABLE, 0) and not allow_unevaluable_gnomad:
        print(
            "ERROR: full-mode input contains unevaluable gnomAD states. "
            "Regenerate gnomAD annotation with AN/callability or explicit gnomAD_state, "
            "or rerun with --allow-unevaluable-gnomad for a deliberately conservative analysis.",
            file=sys.stderr,
        )
        for state, count in sorted(state_counts.items()):
            print(f"  gnomAD state {state}: {count}", file=sys.stderr)
        for reason, count in sorted(reason_counts.items()):
            print(f"  gnomAD reason {reason}: {count}", file=sys.stderr)
        sys.exit(2)

    out_fields = list(classified_rows[0].keys())
    with open(output_path, "w", newline="") as fout:
        writer = csv.DictWriter(
            fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(classified_rows)

    print(f"SLT R1 Classification Complete [{mode}]", file=sys.stderr)
    print(f"  Total variants: {n}", file=sys.stderr)
    for tier in ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]:
        count = tier_counts.get(tier, 0)
        pct = 100 * count / n if n else 0
        print(f"  {tier}: {count} ({pct:.1f}%)", file=sys.stderr)
    for state in [GNOMAD_RARE_CALLABLE, GNOMAD_COMMON, GNOMAD_UNEVALUABLE]:
        count = state_counts.get(state, 0)
        pct = 100 * count / n if n else 0
        print(f"  gnomAD {state}: {count} ({pct:.1f}%)", file=sys.stderr)
    if mode == "full":
        pct = 100 * posterior_present / n if n else 0
        print(
            f"  PureCN posterior present: {posterior_present} ({pct:.1f}%)",
            file=sys.stderr,
        )
    print(f"  Output: {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="R1 SLT classifier with callability-aware gnomAD Layer 2",
    )
    parser.add_argument("--input", "-i", required=True, help="Input TSV/MAF file")
    parser.add_argument("--output", "-o", required=True, help="Output TSV with SLT tiers")
    parser.add_argument(
        "--mode",
        choices=["full", "annotation_only"],
        default="full",
        help="Classification mode",
    )
    parser.add_argument(
        "--degraded",
        action="store_true",
        help="Alias for --mode annotation_only, retained for submitted-script compatibility",
    )
    parser.add_argument(
        "--annotation-only",
        action="store_true",
        help="Alias for --mode annotation_only, retained for prior public-repository examples",
    )
    parser.add_argument(
        "--allow-unevaluable-gnomad",
        action="store_true",
        help="Allow full-mode output when gnomAD state is unevaluable; default is to fail fast",
    )
    args = parser.parse_args()

    mode = "annotation_only" if args.degraded or args.annotation_only else args.mode
    process_file(
        args.input,
        args.output,
        mode=mode,
        allow_unevaluable_gnomad=args.allow_unevaluable_gnomad,
    )


if __name__ == "__main__":
    main()
