#!/usr/bin/env python3
"""
PureCN Component AUROC Evaluation — Protocol §2.4

Computes AUROC using PureCN POSTERIOR.SOMATIC as a continuous classifier
for somatic status. This evaluates the PureCN posterior alone (not SLT tiers).

Usage:
    python3 compute_purecn_auroc.py \
        --eval-tsv  slt/evaluation_per_variant.tsv \
        --slt-tsv   slt/slt_classified.tsv \
        --eval-bed  data/evaluation_regions/final_evaluation_regions_with_cov.bed \
        --output    results/purecn_auroc_results.tsv \
        --fig       figures/fig_purecn_roc.pdf
"""

import argparse
import bisect
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


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


def compute_auroc(labels, scores):
    """Compute AUROC from binary labels and continuous scores."""
    # Sort by score descending
    pairs = sorted(zip(scores, labels), reverse=True)

    tp = 0
    fp = 0
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos

    if n_pos == 0 or n_neg == 0:
        return 0.0, [], [], []

    tpr_list = [0.0]
    fpr_list = [0.0]
    thresholds = []
    prev_score = None

    for score, label in pairs:
        if score != prev_score:
            tpr_list.append(tp / n_pos)
            fpr_list.append(fp / n_neg)
            thresholds.append(score)
            prev_score = score

        if label == 1:
            tp += 1
        else:
            fp += 1

    tpr_list.append(tp / n_pos)
    fpr_list.append(fp / n_neg)

    # Trapezoidal AUC
    auroc = 0.0
    for i in range(1, len(fpr_list)):
        auroc += (fpr_list[i] - fpr_list[i-1]) * (tpr_list[i] + tpr_list[i-1]) / 2

    return auroc, fpr_list, tpr_list, thresholds


def delong_variance(labels, scores):
    """Estimate AUROC variance using DeLong's method (simplified).

    Returns standard error of AUROC estimate.
    """
    pos_scores = [s for s, l in zip(scores, labels) if l == 1]
    neg_scores = [s for s, l in zip(scores, labels) if l == 0]

    m = len(pos_scores)
    n = len(neg_scores)

    if m == 0 or n == 0:
        return 0.0

    # Placement values
    pos_placements = []
    for ps in pos_scores:
        count = sum(1 for ns in neg_scores if ps > ns) + 0.5 * sum(1 for ns in neg_scores if ps == ns)
        pos_placements.append(count / n)

    neg_placements = []
    for ns in neg_scores:
        count = sum(1 for ps in pos_scores if ps > ns) + 0.5 * sum(1 for ps in pos_scores if ps == ns)
        neg_placements.append(count / m)

    # Variance components
    s10 = np.var(pos_placements, ddof=1) if m > 1 else 0
    s01 = np.var(neg_placements, ddof=1) if n > 1 else 0

    var_auc = s10 / m + s01 / n
    se = np.sqrt(var_auc)
    return se


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eval-tsv", required=True)
    parser.add_argument("--slt-tsv", required=True)
    parser.add_argument("--eval-bed", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fig", required=True)
    args = parser.parse_args()

    # Load truth labels from eval TSV
    print("Loading evaluation data...")
    eval_truth = {}
    with open(args.eval_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            eval_truth[key] = row["in_truth"] == "True"

    # Load posteriors from SLT classified TSV
    print("Loading SLT classifications...")
    regions = parse_bed(args.eval_bed)

    labels = []
    scores = []
    n_with_posterior = 0
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

            posterior_str = row.get("POSTERIOR.SOMATIC", "")
            if not posterior_str or posterior_str in ("", "NA", "nan"):
                continue

            try:
                posterior = float(posterior_str)
            except ValueError:
                continue

            n_with_posterior += 1
            key = (chrom, pos, ref, alt)
            is_truth = eval_truth.get(key, False)

            labels.append(1 if is_truth else 0)
            scores.append(posterior)

    print(f"  Variants in eval regions: {n_in_eval}")
    print(f"  With PureCN posterior: {n_with_posterior}")
    print(f"  Truth positives: {sum(labels)}")
    print(f"  False positives: {len(labels) - sum(labels)}")

    # Compute AUROC
    auroc, fpr_list, tpr_list, thresholds = compute_auroc(labels, scores)
    se = delong_variance(labels, scores)
    ci_lo = max(0, auroc - 1.96 * se)
    ci_hi = min(1, auroc + 1.96 * se)

    print(f"\n  PureCN POSTERIOR.SOMATIC AUROC: {auroc:.4f} [95% CI: {ci_lo:.4f}–{ci_hi:.4f}]")
    print(f"  DeLong SE: {se:.4f}")

    # Find optimal threshold (Youden's J)
    best_j = -1
    best_thresh = 0
    for i, thresh in enumerate(thresholds):
        if i + 1 < len(tpr_list):
            j = tpr_list[i+1] - fpr_list[i+1]
            if j > best_j:
                best_j = j
                best_thresh = thresh
                best_tpr = tpr_list[i+1]
                best_fpr = fpr_list[i+1]

    print(f"  Optimal threshold (Youden's J): {best_thresh:.3f}")
    print(f"    TPR at optimal: {best_tpr:.4f}")
    print(f"    FPR at optimal: {best_fpr:.4f}")

    # Write output
    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        f.write("Metric\tValue\n")
        f.write(f"AUROC\t{auroc:.4f}\n")
        f.write(f"AUROC_95CI_lo\t{ci_lo:.4f}\n")
        f.write(f"AUROC_95CI_hi\t{ci_hi:.4f}\n")
        f.write(f"DeLong_SE\t{se:.4f}\n")
        f.write(f"N_variants\t{len(labels)}\n")
        f.write(f"N_positive\t{sum(labels)}\n")
        f.write(f"N_negative\t{len(labels) - sum(labels)}\n")
        f.write(f"Optimal_threshold\t{best_thresh:.3f}\n")
        f.write(f"Optimal_TPR\t{best_tpr:.4f}\n")
        f.write(f"Optimal_FPR\t{best_fpr:.4f}\n")

    print(f"\nSaved: {outpath}")

    # Plot ROC curve
    figpath = Path(args.fig)
    figpath.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(1, 1, figsize=(5.5, 5))

    ax.plot(fpr_list, tpr_list, color="#2166ac", linewidth=2,
            label=f"PureCN posterior (AUROC = {auroc:.3f})")
    ax.plot([0, 1], [0, 1], ":", color="gray", linewidth=1, alpha=0.5)

    # Mark SLT thresholds on ROC
    for thresh, label, color in [
        (0.8, "SLT-A (0.8)", "#1b9e77"),
        (0.5, "SLT-B (0.5)", "#d95f02"),
        (0.2, "SLT-C (0.2)", "#7570b3"),
    ]:
        # Find nearest point on ROC
        dists = [abs(t - thresh) for t in thresholds]
        if dists:
            idx = min(range(len(dists)), key=lambda i: dists[i])
            if idx + 1 < len(tpr_list):
                ax.scatter(fpr_list[idx+1], tpr_list[idx+1], c=color,
                           s=80, zorder=5, edgecolors="black", linewidth=0.8)
                ax.annotate(label,
                            (fpr_list[idx+1], tpr_list[idx+1]),
                            xytext=(fpr_list[idx+1] + 0.05, tpr_list[idx+1] - 0.05),
                            fontsize=8, fontweight="bold", color=color,
                            arrowprops=dict(arrowstyle="->", color=color, lw=0.8))

    ax.set_xlabel("False Positive Rate", fontsize=11)
    ax.set_ylabel("True Positive Rate", fontsize=11)
    ax.set_title("PureCN POSTERIOR.SOMATIC — ROC Curve", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    # Add text box with stats
    stats_text = (
        f"N = {len(labels)} variants\n"
        f"TP = {sum(labels)}, FP = {len(labels) - sum(labels)}\n"
        f"95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]"
    )
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
            fontsize=8, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    fig.tight_layout()
    fig.savefig(str(figpath), dpi=300, bbox_inches="tight")
    fig.savefig(str(figpath.with_suffix(".png")), dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {figpath}")
    print(f"Saved: {figpath.with_suffix('.png')}")


if __name__ == "__main__":
    main()
