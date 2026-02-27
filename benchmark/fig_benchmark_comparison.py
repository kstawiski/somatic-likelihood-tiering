#!/usr/bin/env python3
"""
Benchmark Comparison Figure — SLT vs Competitor Methods

Creates a 2-panel publication-ready figure:
  (A) Sensitivity-PPV scatter comparing all methods
  (B) TP/FP bar chart for each method

Usage:
    python3 fig_benchmark_comparison.py \
        --slt-metrics  results/seqc2_metrics_with_ci.tsv \
        --baseline     results/baseline_evaluation.tsv \
        --output       figures/fig_benchmark_comparison.pdf
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    hw = (z / denom) * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return (max(0.0, center - hw), min(1.0, center + hw))


TIER_COLORS = {
    "SLT-A": "#1b9e77",
    "SLT-B": "#d95f02",
    "SLT-C": "#7570b3",
    "SLT-D": "#e7298a",
}

BASELINE_COLORS = {
    "PureCN >= 0.5": "#e6ab02",
    "PureCN >= 0.8": "#a6761d",
    "PureCN >= 0.2": "#66a61e",
    "gnomAD < 1e-4": "#666666",
}

CALLER_COLORS = {
    "DeepSomatic PASS": "#e41a1c",
    "VarDict PASS": "#377eb8",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slt-metrics", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--deepsomatic-eval", default=None,
                        help="DeepSomatic evaluation summary (evaluate_competitor.py output)")
    parser.add_argument("--n-truth", type=int, default=455,
                        help="Authoritative truth count for consistent CI calculation")
    parser.add_argument("--vardict-eval", default=None,
                        help="VarDict evaluation summary (evaluate_competitor.py output)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Parse SLT metrics
    slt_data = []
    with open(args.slt_metrics) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Sensitivity"] != "N/A" and row["Threshold"].startswith("≥SLT-"):
                slt_data.append({
                    "label": row["Threshold"].replace("≥", ">= "),
                    "sens": float(row["Sensitivity"]),
                    "ppv": float(row["PPV"]),
                    "tp": int(row["TP"]),
                    "fp": int(row["FP"]),
                    "fn": int(row["FN"]),
                    "called": int(row["Called"]),
                    "sens_lo": float(row["Sens_95CI_lo"]),
                    "sens_hi": float(row["Sens_95CI_hi"]),
                    "ppv_lo": float(row["PPV_95CI_lo"]),
                    "ppv_hi": float(row["PPV_95CI_hi"]),
                })

    # Parse baseline metrics
    baseline_data = []
    label_map = {
        "purecn_only_0.5": "PureCN >= 0.5",
        "purecn_only_0.8": "PureCN >= 0.8",
        "purecn_only_0.2": "PureCN >= 0.2",
        "gnomad_only_0.0001": "gnomAD < 1e-4",
    }
    with open(args.baseline) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            method = row["Method"]
            if method in label_map:
                baseline_data.append({
                    "label": label_map[method],
                    "sens": float(row["Sensitivity"]),
                    "ppv": float(row["PPV"]),
                    "tp": int(row["TP"]),
                    "fp": int(row["FP"]),
                    "fn": int(row["FN"]),
                    "called": int(row["Called"]),
                    "sens_lo": float(row["Sens_95CI_lo"]),
                    "sens_hi": float(row["Sens_95CI_hi"]),
                    "ppv_lo": float(row["PPV_95CI_lo"]),
                    "ppv_hi": float(row["PPV_95CI_hi"]),
                })

    # Parse DeepSomatic evaluation (if available)
    # Format: "PASS only                438   376      62    76  0.8319  0.8584  0.8449"
    # Columns: Filter(variable width)  Called  TP    FP    FN  Sens    Prec    F1
    caller_data = []
    if args.deepsomatic_eval:
        with open(args.deepsomatic_eval) as f:
            for line in f:
                line = line.strip()
                if line.startswith("PASS only"):
                    # "PASS only" is the filter name; numeric fields follow
                    parts = line.split()
                    # parts: ["PASS", "only", "438", "376", "62", "76", ...]
                    called = int(parts[2])
                    tp = int(parts[3])
                    fp = int(parts[4])
                    fn = args.n_truth - tp
                    sens = tp / args.n_truth
                    ppv = tp / (tp + fp)
                    s_lo, s_hi = wilson_ci(tp, args.n_truth)
                    p_lo, p_hi = wilson_ci(tp, tp + fp)
                    caller_data.append({
                        "label": "DeepSomatic PASS",
                        "sens": sens, "ppv": ppv,
                        "tp": tp, "fp": fp, "fn": fn,
                        "called": called,
                        "sens_lo": s_lo, "sens_hi": s_hi,
                        "ppv_lo": p_lo, "ppv_hi": p_hi,
                    })

    # Parse VarDict evaluation (if available)
    if args.vardict_eval:
        with open(args.vardict_eval) as f:
            for line in f:
                line = line.strip()
                if line.startswith("PASS only"):
                    parts = line.split()
                    called = int(parts[2])
                    tp = int(parts[3])
                    fp = int(parts[4])
                    fn = args.n_truth - tp
                    sens = tp / args.n_truth
                    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
                    s_lo, s_hi = wilson_ci(tp, args.n_truth)
                    p_lo, p_hi = wilson_ci(tp, tp + fp) if (tp + fp) > 0 else (0, 0)
                    caller_data.append({
                        "label": "VarDict PASS",
                        "sens": sens, "ppv": ppv,
                        "tp": tp, "fp": fp, "fn": fn,
                        "called": called,
                        "sens_lo": s_lo, "sens_hi": s_hi,
                        "ppv_lo": p_lo, "ppv_hi": p_hi,
                    })

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ─── Panel A: Sensitivity-PPV scatter ───────────────────────────────────
    ax = axes[0]

    # Plot SLT points with error bars
    slt_tiers = ["SLT-A", "SLT-B", "SLT-C", "SLT-D"]
    for i, d in enumerate(slt_data):
        tier = slt_tiers[i] if i < len(slt_tiers) else f"SLT-{i}"
        color = TIER_COLORS.get(tier, "#333333")

        xerr_lo = d["sens"] - d["sens_lo"]
        xerr_hi = d["sens_hi"] - d["sens"]
        yerr_lo = d["ppv"] - d["ppv_lo"]
        yerr_hi = d["ppv_hi"] - d["ppv"]

        ax.errorbar(d["sens"], d["ppv"],
                     xerr=[[xerr_lo], [xerr_hi]],
                     yerr=[[yerr_lo], [yerr_hi]],
                     fmt="o", color=color, markersize=10, markeredgecolor="black",
                     markeredgewidth=0.8, capsize=3, linewidth=1.2,
                     label=f">= {tier}")

    # Plot baseline points
    for d in baseline_data:
        color = BASELINE_COLORS.get(d["label"], "#999999")
        marker = "s" if "PureCN" in d["label"] else "D"

        xerr_lo = d["sens"] - d["sens_lo"]
        xerr_hi = d["sens_hi"] - d["sens"]
        yerr_lo = d["ppv"] - d["ppv_lo"]
        yerr_hi = d["ppv_hi"] - d["ppv"]

        ax.errorbar(d["sens"], d["ppv"],
                     xerr=[[xerr_lo], [xerr_hi]],
                     yerr=[[yerr_lo], [yerr_hi]],
                     fmt=marker, color=color, markersize=8, markeredgecolor="black",
                     markeredgewidth=0.8, capsize=3, linewidth=1.0,
                     label=d["label"])

    # Plot end-to-end caller points (star markers)
    for d in caller_data:
        color = CALLER_COLORS.get(d["label"], "#e41a1c")
        xerr_lo = d["sens"] - d["sens_lo"]
        xerr_hi = d["sens_hi"] - d["sens"]
        yerr_lo = d["ppv"] - d["ppv_lo"]
        yerr_hi = d["ppv_hi"] - d["ppv"]

        ax.errorbar(d["sens"], d["ppv"],
                     xerr=[[xerr_lo], [xerr_hi]],
                     yerr=[[yerr_lo], [yerr_hi]],
                     fmt="*", color=color, markersize=14, markeredgecolor="black",
                     markeredgewidth=0.8, capsize=3, linewidth=1.2,
                     label=d["label"])

    # Connect SLT points with line
    slt_sens = [d["sens"] for d in slt_data]
    slt_ppv = [d["ppv"] for d in slt_data]
    ax.plot(slt_sens, slt_ppv, "--", color="gray", linewidth=0.8, alpha=0.5, zorder=0)

    ax.set_xlabel("Sensitivity", fontsize=11)
    ax.set_ylabel("Positive Predictive Value", fontsize=11)
    ax.set_title("(A) Sensitivity vs PPV", fontsize=12, fontweight="bold")
    ax.legend(fontsize=7.5, loc="upper right", ncol=1, framealpha=0.9)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.2)

    # ─── Panel B: TP/FP bar comparison ──────────────────────────────────────
    ax = axes[1]

    # Combine all methods
    all_methods = []
    all_tp = []
    all_fp = []
    all_colors = []

    for d in slt_data:
        all_methods.append(d["label"])
        all_tp.append(d["tp"])
        all_fp.append(d["fp"])
        tier = d["label"].replace(">= ", "").replace(" ", "")
        all_colors.append(TIER_COLORS.get(tier, "#333333"))

    for d in baseline_data:
        all_methods.append(d["label"])
        all_tp.append(d["tp"])
        all_fp.append(d["fp"])
        all_colors.append(BASELINE_COLORS.get(d["label"], "#999999"))

    for d in caller_data:
        all_methods.append(d["label"])
        all_tp.append(d["tp"])
        all_fp.append(d["fp"])
        all_colors.append(CALLER_COLORS.get(d["label"], "#e41a1c"))

    y_pos = np.arange(len(all_methods))
    bar_height = 0.35

    # TP bars (positive direction)
    bars_tp = ax.barh(y_pos - bar_height/2, all_tp, bar_height,
                       color=all_colors, edgecolor="black", linewidth=0.5,
                       alpha=0.9, label="True Positives")

    # FP bars (lighter)
    bars_fp = ax.barh(y_pos + bar_height/2, all_fp, bar_height,
                       color=all_colors, edgecolor="black", linewidth=0.5,
                       alpha=0.3, label="False Positives", hatch="//")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_methods, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Variant Count", fontsize=11)
    ax.set_title("(B) TP and FP by Method", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xscale("log")
    ax.set_xlim(1, 10000)

    # Add count labels
    for i, (tp, fp) in enumerate(zip(all_tp, all_fp)):
        max_val = max(tp, fp)
        ax.text(max_val * 1.3, i, f"TP={tp}, FP={fp}",
                va="center", fontsize=6.5)

    fig.tight_layout(w_pad=2.0)

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outpath), dpi=300, bbox_inches="tight")
    fig.savefig(str(outpath.with_suffix(".png")), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")
    print(f"Saved: {outpath.with_suffix('.png')}")


if __name__ == "__main__":
    main()
