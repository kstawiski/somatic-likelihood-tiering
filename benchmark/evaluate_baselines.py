#!/usr/bin/env python3
"""
Evaluate PureCN-only and gnomAD-only baselines against SEQC2 truth set.

Uses the SLT classified TSV (which contains all required fields) and the
per-variant evaluation TSV (which has truth labels).

Protocol: 

Usage:
    python3 evaluate_baselines.py \
        --slt-tsv     slt/slt_classified.tsv \
        --eval-tsv    slt/evaluation_per_variant.tsv \
        --eval-bed    data/evaluation_regions/final_evaluation_regions_with_cov.bed \
        --eval-summary slt/evaluation_summary.txt \
        --output      results/baseline_evaluation.tsv
"""

import argparse
import bisect
import csv
import math
import sys
from collections import defaultdict


def safe_float(value, default=None):
    if value is None or value in ("", ".", "NA", "nan", "NaN"):
        return default
    try:
        f = float(value)
        return default if math.isnan(f) else f
    except (ValueError, TypeError):
        return default


def parse_bed(path):
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


def wilson_ci(successes, total, z=1.96):
    if total == 0:
        return (0.0, 0.0)
    p = successes / total
    denom = 1 + z**2 / total
    center = (p + z**2 / (2 * total)) / denom
    half_width = (z / denom) * math.sqrt(p * (1 - p) / total + z**2 / (4 * total**2))
    return (max(0.0, center - half_width), min(1.0, center + half_width))


def parse_eval_summary(path):
    """Parse evaluation_summary.txt for authoritative truth count."""
    n_truth = None
    n_truth_snv = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Truth in eval regions:"):
                n_truth = int(line.split(":")[1].strip())
            elif line.startswith("Truth SNVs in eval regions:"):
                n_truth_snv = int(line.split(":")[1].strip())
    return n_truth, n_truth_snv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slt-tsv", required=True)
    parser.add_argument("--eval-tsv", required=True)
    parser.add_argument("--eval-bed", required=True)
    parser.add_argument("--eval-summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load truth labels from evaluation per-variant TSV
    print("Loading truth labels...")
    truth_set = set()
    with open(args.eval_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["in_truth"] == "True":
                truth_set.add((row["chrom"], row["pos"], row["ref"], row["alt"]))
    print(f"  Truth variants in per-variant TSV: {len(truth_set)}")

    # Get authoritative truth count (includes uncalled variants)
    n_truth, n_truth_snv = parse_eval_summary(args.eval_summary)
    print(f"  Authoritative truth count: {n_truth} (SNVs: {n_truth_snv})")

    # Load evaluation regions
    regions = parse_bed(args.eval_bed)

    # Process SLT classified TSV
    print("Processing SLT classified TSV...")
    results = {
        "purecn_only_0.5": {"tp": 0, "fp": 0, "called": 0},
        "purecn_only_0.8": {"tp": 0, "fp": 0, "called": 0},
        "purecn_only_0.2": {"tp": 0, "fp": 0, "called": 0},
        "gnomad_only_0.0001": {"tp": 0, "fp": 0, "called": 0},
        "gnomad_only_0.001": {"tp": 0, "fp": 0, "called": 0},
    }

    n_in_eval = 0
    with open(args.slt_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("Chromosome", "")
            pos = row.get("Start_Position", "")
            ref = row.get("Reference_Allele", "")
            alt = row.get("Tumor_Seq_Allele2", "")

            if not chrom or not pos:
                continue

            if not in_region(chrom, int(pos), regions):
                continue
            n_in_eval += 1

            key = (chrom, pos, ref, alt)
            is_truth = key in truth_set

            # PureCN posterior
            posterior = safe_float(row.get("POSTERIOR.SOMATIC") or row.get("slt_posterior_used"))

            # gnomAD AF
            gnomad_af = safe_float(row.get("gnomAD_AF") or row.get("gnomAD_exome_AF"))

            # PureCN-only at 0.5
            if posterior is not None and posterior >= 0.5:
                results["purecn_only_0.5"]["called"] += 1
                if is_truth:
                    results["purecn_only_0.5"]["tp"] += 1
                else:
                    results["purecn_only_0.5"]["fp"] += 1

            # PureCN-only at 0.8
            if posterior is not None and posterior >= 0.8:
                results["purecn_only_0.8"]["called"] += 1
                if is_truth:
                    results["purecn_only_0.8"]["tp"] += 1
                else:
                    results["purecn_only_0.8"]["fp"] += 1

            # PureCN-only at 0.2
            if posterior is not None and posterior >= 0.2:
                results["purecn_only_0.2"]["called"] += 1
                if is_truth:
                    results["purecn_only_0.2"]["tp"] += 1
                else:
                    results["purecn_only_0.2"]["fp"] += 1

            # gnomAD-only at 0.0001
            if gnomad_af is None or gnomad_af < 0.0001:
                results["gnomad_only_0.0001"]["called"] += 1
                if is_truth:
                    results["gnomad_only_0.0001"]["tp"] += 1
                else:
                    results["gnomad_only_0.0001"]["fp"] += 1

            # gnomAD-only at 0.001
            if gnomad_af is None or gnomad_af < 0.001:
                results["gnomad_only_0.001"]["called"] += 1
                if is_truth:
                    results["gnomad_only_0.001"]["tp"] += 1
                else:
                    results["gnomad_only_0.001"]["fp"] += 1

    print(f"  Variants in eval regions: {n_in_eval}")

    # Compute metrics
    print("\n" + "=" * 80)
    print("Baseline Competitor Evaluation — SEQC2 HCC1395")
    print("=" * 80)

    header = f"{'Method':<25} {'Called':>7} {'TP':>5} {'FP':>7} {'FN':>5} {'Sens':>8} {'[95% CI]':>20} {'PPV':>8} {'[95% CI]':>20} {'F1':>7}"
    print(header)
    print("-" * len(header))

    output_rows = []
    for method, r in results.items():
        tp = r["tp"]
        fp = r["fp"]
        fn = n_truth - tp
        called = r["called"]
        sens = tp / n_truth if n_truth > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1 = 2 * sens * ppv / (sens + ppv) if (sens + ppv) > 0 else 0

        sens_lo, sens_hi = wilson_ci(tp, n_truth)
        ppv_lo, ppv_hi = wilson_ci(tp, tp + fp)

        sens_ci = f"[{sens_lo:.4f}, {sens_hi:.4f}]"
        ppv_ci = f"[{ppv_lo:.4f}, {ppv_hi:.4f}]"

        print(f"{method:<25} {called:>7} {tp:>5} {fp:>7} {fn:>5} {sens:>8.4f} {sens_ci:>20} {ppv:>8.4f} {ppv_ci:>20} {f1:>7.4f}")

        output_rows.append({
            "Method": method,
            "Called": called,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "Sensitivity": f"{sens:.4f}",
            "Sens_95CI_lo": f"{sens_lo:.4f}",
            "Sens_95CI_hi": f"{sens_hi:.4f}",
            "PPV": f"{ppv:.4f}",
            "PPV_95CI_lo": f"{ppv_lo:.4f}",
            "PPV_95CI_hi": f"{ppv_hi:.4f}",
            "F1": f"{f1:.4f}",
        })

    # Write TSV
    from pathlib import Path
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(output_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        for row in output_rows:
            writer.writerow(row)

    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
