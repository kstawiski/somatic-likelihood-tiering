#!/usr/bin/env python3
"""
Evaluate competitor caller VCF against SEQC2 truth set.

Computes sensitivity, precision, F1 for any VCF caller output
at various FILTER thresholds (PASS-only, all variants).
Also applies simple AF-based tiering for fair comparison.

Usage:
    python3 evaluate_competitor.py \
        --caller-vcf vardict/vardict.vcf.gz \
        --caller-name VarDict \
        --truth-snv   data/truth_snv_eval_regions.vcf.gz \
        --truth-indel data/truth_indel_eval_regions.vcf.gz \
        --eval-bed    data/evaluation_regions/final_evaluation_regions_with_cov.bed \
        --output      vardict/evaluation_summary.txt
"""

import argparse
import gzip
import bisect
import sys
from collections import defaultdict


def parse_vcf(path):
    """Parse VCF, return list of variant dicts."""
    variants = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            filt = parts[6]
            info = parts[7]

            # Parse AF from INFO or FORMAT
            af = None
            # Try INFO field first
            for field in info.split(";"):
                if field.startswith("AF="):
                    try:
                        af = float(field.split("=")[1].split(",")[0])
                    except (ValueError, IndexError):
                        pass

            # Try FORMAT field (for VarDict: AD field)
            if af is None and len(parts) >= 10:
                fmt_keys = parts[8].split(":")
                fmt_vals = parts[9].split(":")
                fmt_dict = dict(zip(fmt_keys, fmt_vals))
                if "AF" in fmt_dict:
                    try:
                        af = float(fmt_dict["AF"])
                    except ValueError:
                        pass
                elif "AD" in fmt_dict:
                    try:
                        ads = [int(x) for x in fmt_dict["AD"].split(",")]
                        if sum(ads) > 0:
                            af = ads[1] / sum(ads) if len(ads) > 1 else 0
                    except (ValueError, IndexError):
                        pass

            # Determine variant type
            vtype = "SNV" if len(ref) == 1 and len(alt) == 1 else "INDEL"

            # Handle multi-allelic (simplistic: take first alt)
            for a in alt.split(","):
                variants.append({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": a,
                    "filter": filt,
                    "af": af,
                    "vtype": "SNV" if len(ref) == 1 and len(a) == 1 else "INDEL",
                })

    return variants


def parse_truth_vcf(path):
    """Parse truth VCF, return set of (chrom, pos, ref, alt)."""
    truth = set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
            for a in alt.split(","):
                truth.add((chrom, pos, ref, a))
    return truth


def parse_bed(path):
    """Parse BED file, return dict of chrom -> sorted intervals."""
    regions = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            regions[parts[0]].append((int(parts[1]), int(parts[2])))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def in_region(chrom, pos, regions):
    """Check if position falls in any evaluation region."""
    if chrom not in regions:
        return False
    intervals = regions[chrom]
    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx >= 0 and intervals[idx][0] <= pos < intervals[idx][1]:
        return True
    if idx + 1 < len(intervals) and intervals[idx + 1][0] <= pos < intervals[idx + 1][1]:
        return True
    return False


def evaluate(caller_variants, truth_set, regions, filter_fn, label):
    """Evaluate caller variants matching filter_fn against truth."""
    filtered = [v for v in caller_variants if filter_fn(v)]

    # Filter to eval regions
    in_eval = [v for v in filtered if in_region(v["chrom"], v["pos"], regions)]

    truth_in_eval = {t for t in truth_set if in_region(t[0], t[1], regions)}
    n_truth = len(truth_in_eval)

    # Count TP/FP
    called_keys = set()
    tp = 0
    fp = 0
    for v in in_eval:
        key = (v["chrom"], v["pos"], v["ref"], v["alt"])
        if key in called_keys:
            continue  # Skip duplicate calls
        called_keys.add(key)
        if key in truth_in_eval:
            tp += 1
        else:
            fp += 1

    fn = n_truth - tp
    sens = tp / n_truth if n_truth > 0 else 0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0

    return {
        "label": label,
        "called": len(called_keys),
        "tp": tp, "fp": fp, "fn": fn,
        "sens": sens, "prec": prec, "f1": f1,
        "n_truth": n_truth,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--caller-vcf", required=True)
    parser.add_argument("--caller-name", required=True)
    parser.add_argument("--truth-snv", required=True)
    parser.add_argument("--truth-indel", required=True)
    parser.add_argument("--eval-bed", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    print(f"Evaluating {args.caller_name}...")

    print("Loading caller VCF...")
    caller_variants = parse_vcf(args.caller_vcf)
    print(f"  {len(caller_variants)} variants")

    print("Loading truth sets...")
    truth_snv = parse_truth_vcf(args.truth_snv)
    truth_indel = parse_truth_vcf(args.truth_indel)
    truth_all = truth_snv | truth_indel
    print(f"  Truth SNVs: {len(truth_snv)}, indels: {len(truth_indel)}")

    print("Loading evaluation regions...")
    regions = parse_bed(args.eval_bed)
    n_regions = sum(len(v) for v in regions.values())
    print(f"  {n_regions} regions")

    # Define filter functions
    filters = [
        ("All variants", lambda v: True),
        ("PASS only", lambda v: v["filter"] in ("PASS", ".")),
        ("SNV all", lambda v: v["vtype"] == "SNV"),
        ("SNV PASS", lambda v: v["vtype"] == "SNV" and v["filter"] in ("PASS", ".")),
    ]

    # AF-based tiers (approximate SLT comparison)
    af_tiers = [
        ("AF≥0.20", lambda v: v["af"] is not None and v["af"] >= 0.20),
        ("AF≥0.10", lambda v: v["af"] is not None and v["af"] >= 0.10),
        ("AF≥0.05", lambda v: v["af"] is not None and v["af"] >= 0.05),
        ("AF≥0.01", lambda v: v["af"] is not None and v["af"] >= 0.01),
    ]

    all_filters = filters + af_tiers

    results = []
    for label, fn in all_filters:
        r = evaluate(caller_variants, truth_all, regions, fn, label)
        results.append(r)

    # Write output
    with open(args.output, "w") as f:
        f.write("=" * 70 + "\n")
        f.write(f"{args.caller_name} — SEQC2 HCC1395 Benchmark Evaluation\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Total caller variants: {len(caller_variants)}\n")
        f.write(f"Truth SNVs: {len(truth_snv)}\n")
        f.write(f"Truth indels: {len(truth_indel)}\n\n")

        f.write("-" * 70 + "\n")
        f.write(f"{'Filter':<20} {'Called':>7} {'TP':>5} {'FP':>7} "
                f"{'FN':>5} {'Sens':>7} {'Prec':>7} {'F1':>7}\n")
        f.write("-" * 70 + "\n")

        for r in results:
            f.write(f"{r['label']:<20} {r['called']:>7} {r['tp']:>5} {r['fp']:>7} "
                    f"{r['fn']:>5} {r['sens']:>7.4f} {r['prec']:>7.4f} {r['f1']:>7.4f}\n")
        f.write("-" * 70 + "\n")

    # Also print to stdout
    with open(args.output) as f:
        print(f.read())

    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
