#!/usr/bin/env python3
"""
Generate head-to-head comparison table for SEQC2 benchmark.

Compares all methods at their optimal/specified operating points:
- SLT tiers (cumulative and per-tier)
- PureCN-only baselines
- gnomAD-only baseline
- DeepSomatic PASS
- VarDict PASS (if available)

Output: TSV + formatted markdown table for manuscript

Usage:
    python3 head_to_head_table.py \
        --slt-metrics  results/seqc2_metrics_with_ci.tsv \
        --baseline     results/baseline_evaluation.tsv \
        --deepsomatic  results/deepsomatic_evaluation.txt \
        --vardict      results/vardict_evaluation.txt \
        --n-truth 455 \
        --output       results/benchmark_comparison.tsv
"""

import argparse
import csv
import math


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    hw = (z / denom) * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return (max(0.0, center - hw), min(1.0, center + hw))


def parse_eval_file(path, n_truth):
    """Parse evaluate_competitor.py output, return dict for PASS only."""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("PASS only"):
                parts = line.split()
                called = int(parts[2])
                tp = int(parts[3])
                fp = int(parts[4])
                fn = n_truth - tp
                sens = tp / n_truth
                ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
                f1 = 2 * sens * ppv / (sens + ppv) if (sens + ppv) > 0 else 0
                s_lo, s_hi = wilson_ci(tp, n_truth)
                p_lo, p_hi = wilson_ci(tp, tp + fp) if (tp + fp) > 0 else (0, 0)
                return {
                    "called": called, "tp": tp, "fp": fp, "fn": fn,
                    "sens": sens, "ppv": ppv, "f1": f1,
                    "sens_lo": s_lo, "sens_hi": s_hi,
                    "ppv_lo": p_lo, "ppv_hi": p_hi,
                }
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slt-metrics", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--deepsomatic", default=None)
    parser.add_argument("--vardict", default=None)
    parser.add_argument("--n-truth", type=int, default=455)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []

    # SLT cumulative tiers
    with open(args.slt_metrics) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Threshold"].startswith("≥SLT-"):
                tier = row["Threshold"]
                rows.append({
                    "Method": tier,
                    "Category": "SLT (this work)",
                    "Type": "Post-hoc classifier",
                    "Called": int(row["Called"]),
                    "TP": int(row["TP"]),
                    "FP": int(row["FP"]),
                    "FN": int(row["FN"]),
                    "Sensitivity": f"{float(row['Sensitivity']):.1%}",
                    "Sens_CI": f"[{float(row['Sens_95CI_lo']):.1%}--{float(row['Sens_95CI_hi']):.1%}]",
                    "PPV": f"{float(row['PPV']):.1%}",
                    "PPV_CI": f"[{float(row['PPV_95CI_lo']):.1%}--{float(row['PPV_95CI_hi']):.1%}]",
                    "F1": f"{float(row['F1']):.3f}",
                })

    # PureCN baselines
    label_map = {
        "purecn_only_0.8": "PureCN ≥0.8",
        "purecn_only_0.5": "PureCN ≥0.5",
        "purecn_only_0.2": "PureCN ≥0.2",
        "gnomad_only_0.0001": "gnomAD <1e-4",
    }
    with open(args.baseline) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            method = row["Method"]
            if method in label_map:
                rows.append({
                    "Method": label_map[method],
                    "Category": "Single-feature baseline",
                    "Type": "Post-hoc filter",
                    "Called": int(row["Called"]),
                    "TP": int(row["TP"]),
                    "FP": int(row["FP"]),
                    "FN": int(row["FN"]),
                    "Sensitivity": f"{float(row['Sensitivity']):.1%}",
                    "Sens_CI": f"[{float(row['Sens_95CI_lo']):.1%}--{float(row['Sens_95CI_hi']):.1%}]",
                    "PPV": f"{float(row['PPV']):.1%}",
                    "PPV_CI": f"[{float(row['PPV_95CI_lo']):.1%}--{float(row['PPV_95CI_hi']):.1%}]",
                    "F1": f"{float(row['F1']):.3f}",
                })

    # DeepSomatic
    if args.deepsomatic:
        d = parse_eval_file(args.deepsomatic, args.n_truth)
        if d:
            rows.append({
                "Method": "DeepSomatic PASS",
                "Category": "End-to-end caller",
                "Type": "Neural network",
                "Called": d["called"],
                "TP": d["tp"],
                "FP": d["fp"],
                "FN": d["fn"],
                "Sensitivity": f"{d['sens']:.1%}",
                "Sens_CI": f"[{d['sens_lo']:.1%}--{d['sens_hi']:.1%}]",
                "PPV": f"{d['ppv']:.1%}",
                "PPV_CI": f"[{d['ppv_lo']:.1%}--{d['ppv_hi']:.1%}]",
                "F1": f"{d['f1']:.3f}",
            })

    # VarDict
    if args.vardict:
        d = parse_eval_file(args.vardict, args.n_truth)
        if d:
            rows.append({
                "Method": "VarDict PASS",
                "Category": "End-to-end caller",
                "Type": "Heuristic caller",
                "Called": d["called"],
                "TP": d["tp"],
                "FP": d["fp"],
                "FN": d["fn"],
                "Sensitivity": f"{d['sens']:.1%}",
                "Sens_CI": f"[{d['sens_lo']:.1%}--{d['sens_hi']:.1%}]",
                "PPV": f"{d['ppv']:.1%}",
                "PPV_CI": f"[{d['ppv_lo']:.1%}--{d['ppv_hi']:.1%}]",
                "F1": f"{d['f1']:.3f}",
            })

    # Write TSV
    fieldnames = ["Method", "Category", "Type", "Called", "TP", "FP", "FN",
                  "Sensitivity", "Sens_CI", "PPV", "PPV_CI", "F1"]
    with open(args.output, "w") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote: {args.output}")
    print(f"  {len(rows)} methods compared")

    # Print markdown table
    print("\n### Markdown Table\n")
    print("| Method | Category | Called | TP | FP | FN | Sensitivity [95% CI] | PPV [95% CI] | F1 |")
    print("|--------|----------|--------|-----|------|------|-----|------|-----|")
    for r in rows:
        print(f"| {r['Method']} | {r['Category']} | {r['Called']} | {r['TP']} | "
              f"{r['FP']} | {r['FN']} | {r['Sensitivity']} {r['Sens_CI']} | "
              f"{r['PPV']} {r['PPV_CI']} | {r['F1']} |")


if __name__ == "__main__":
    main()
