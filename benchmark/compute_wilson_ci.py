#!/usr/bin/env python3
"""
Compute Wilson score confidence intervals for benchmark metrics.

Protocol: 
Adds 95% Wilson score CIs to sensitivity, PPV, and F1 metrics.

Usage:
    python3 compute_wilson_ci.py \
        --eval-summary slt/evaluation_summary.txt \
        --output       results/seqc2_metrics_with_ci.tsv
"""

import argparse
import math
import csv


def wilson_ci(successes: int, total: int, z: float = 1.96) -> tuple:
    """Wilson score confidence interval for a proportion.

    Returns (lower, upper) bounds.
    """
    if total == 0:
        return (0.0, 0.0)

    p = successes / total
    denom = 1 + z**2 / total
    center = (p + z**2 / (2 * total)) / denom
    half_width = (z / denom) * math.sqrt(p * (1 - p) / total + z**2 / (4 * total**2))

    lower = max(0.0, center - half_width)
    upper = min(1.0, center + half_width)
    return (lower, upper)


def parse_eval_summary(path: str) -> dict:
    """Parse evaluation_summary.txt for metrics."""
    summary = {"n_truth": None, "rows": []}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Truth in eval regions:"):
                summary["n_truth"] = int(line.split(":")[1].strip())
            elif line.startswith("Truth SNVs in eval regions:"):
                summary["n_truth_snv"] = int(line.split(":")[1].strip())
            # Parse metric rows
            # Format: Threshold  Called  TP  FP  FN  Sens  Prec  F1
            elif line.startswith(("≥SLT-", ">=SLT-", "SLT-", "SNV_")):
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        row = {
                            "threshold": parts[0],
                            "called": int(parts[1]),
                            "tp": int(parts[2]),
                            "fp": int(parts[3]),
                            "fn": parts[4],  # May be "N/A"
                            "sens": parts[5],
                            "prec": parts[6],
                            "f1": parts[7],
                        }
                        summary["rows"].append(row)
                    except (ValueError, IndexError):
                        pass
    return summary


def main():
    parser = argparse.ArgumentParser(description="Add Wilson CIs to SLT metrics")
    parser.add_argument("--eval-summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    summary = parse_eval_summary(args.eval_summary)
    n_truth = summary["n_truth"]

    print(f"Truth variants in eval regions: {n_truth}")
    print()

    header = [
        "Threshold", "Called", "TP", "FP", "FN",
        "Sensitivity", "Sens_95CI_lo", "Sens_95CI_hi",
        "PPV", "PPV_95CI_lo", "PPV_95CI_hi",
        "F1", "F1_note"
    ]

    rows_out = []

    print(f"{'Threshold':<18} {'Called':>7} {'TP':>5} {'FP':>5} {'FN':>5} "
          f"{'Sens':>8} {'[95% CI]':>18} {'PPV':>8} {'[95% CI]':>18} {'F1':>7}")
    print("-" * 105)

    for row in summary["rows"]:
        tp = row["tp"]
        fp = row["fp"]
        fn_str = row["fn"]
        called = row["called"]

        # Sensitivity CI: TP out of (TP + FN)
        if fn_str != "N/A":
            fn = int(fn_str)
            denom_sens = tp + fn
            sens = tp / denom_sens if denom_sens > 0 else 0
            sens_lo, sens_hi = wilson_ci(tp, denom_sens)
        else:
            fn = None
            sens = None
            sens_lo, sens_hi = None, None

        # PPV CI: TP out of (TP + FP)
        denom_ppv = tp + fp
        ppv = tp / denom_ppv if denom_ppv > 0 else 0
        ppv_lo, ppv_hi = wilson_ci(tp, denom_ppv)

        # F1 (no standard CI for F1; computed from point estimates)
        if sens is not None and ppv > 0:
            f1 = 2 * sens * ppv / (sens + ppv) if (sens + ppv) > 0 else 0
        else:
            f1 = None

        # Output
        out_row = {
            "Threshold": row["threshold"],
            "Called": called,
            "TP": tp,
            "FP": fp,
            "FN": fn if fn is not None else "N/A",
            "Sensitivity": f"{sens:.4f}" if sens is not None else "N/A",
            "Sens_95CI_lo": f"{sens_lo:.4f}" if sens_lo is not None else "N/A",
            "Sens_95CI_hi": f"{sens_hi:.4f}" if sens_hi is not None else "N/A",
            "PPV": f"{ppv:.4f}",
            "PPV_95CI_lo": f"{ppv_lo:.4f}",
            "PPV_95CI_hi": f"{ppv_hi:.4f}",
            "F1": f"{f1:.4f}" if f1 is not None else "N/A",
            "F1_note": "point estimate",
        }
        rows_out.append(out_row)

        # Print
        sens_str = f"{sens:.4f}" if sens is not None else "N/A"
        sens_ci = f"[{sens_lo:.4f}, {sens_hi:.4f}]" if sens_lo is not None else "N/A"
        ppv_ci = f"[{ppv_lo:.4f}, {ppv_hi:.4f}]"
        f1_str = f"{f1:.4f}" if f1 is not None else "N/A"
        fn_disp = str(fn) if fn is not None else "N/A"

        ppv_str = f"{ppv:.4f}"
        print(f"{row['threshold']:<18} {called:>7} {tp:>5} {fp:>5} {fn_disp:>5} "
              f"{sens_str:>8} {sens_ci:>18} {ppv_str:>8} {ppv_ci:>18} {f1_str:>7}")

    # Write TSV
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for r in rows_out:
            writer.writerow(r)

    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
