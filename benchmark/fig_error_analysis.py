#!/usr/bin/env python3
"""
Error Analysis Figure — PureCN Posterior vs SLT Evidence Layer Rescue

Creates a 2-panel figure:
  (A) PureCN posterior distribution: truth vs non-truth (violin/histogram)
  (B) Evidence-layer rescue: Venn-like diagram showing how many TPs are
      captured by PureCN posterior vs evidence layers

Usage:
    python3 fig_error_analysis.py \
        --slt-classified  work/seqc2/slt/slt_classified.tsv \
        --evaluation      work/seqc2/slt/evaluation_per_variant.tsv \
        --output          figures/fig_error_analysis.pdf
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def is_true(val):
    return str(val).lower() in ("true", "1", "yes")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slt-classified", required=True)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load evaluation data
    eval_data = {}
    with open(args.evaluation) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
            eval_data[key] = row

    # Load SLT classified data
    truth_posts = []  # posteriors for truth variants
    nontruth_posts = []  # posteriors for non-truth variants
    # For rescue analysis
    tp_post_high = 0  # TP with posterior >= 0.5
    tp_evidence_only = 0  # TP with posterior < 0.5, tier <= SLT-C
    tp_both = 0  # TP with posterior >= 0.5 AND tier <= SLT-C (overlap)
    fn_count = 0  # FN at >= SLT-C
    total_tp = 0

    with open(args.slt_classified) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            parts = row["variant_key"].split(":")
            key = ":".join(parts[:4]) if len(parts) >= 4 else row["variant_key"]
            ev = eval_data.get(key, {})

            if not is_true(ev.get("in_eval_region", "False")):
                continue

            try:
                post = float(row.get("POSTERIOR.SOMATIC", "") or "0")
            except (ValueError, TypeError):
                post = 0.0

            tier = row.get("slt_tier", "SLT-D")
            is_tp = is_true(ev.get("in_truth", "False"))

            if is_tp:
                truth_posts.append(post)
                if tier in ("SLT-A", "SLT-B", "SLT-C"):
                    total_tp += 1
                    if post >= 0.5:
                        tp_post_high += 1
                    else:
                        tp_evidence_only += 1
                else:
                    fn_count += 1
            else:
                nontruth_posts.append(post)

    print(f"Truth variants in eval region: {len(truth_posts)}")
    print(f"Non-truth variants: {len(nontruth_posts)}")
    print(f"TP at >=SLT-C: {total_tp}")
    print(f"  Post>=0.5: {tp_post_high}")
    print(f"  Post<0.5 (evidence rescue): {tp_evidence_only}")
    print(f"FN (SLT-D truth): {fn_count}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ─── Panel A: Posterior distribution ───────────────────────────────
    ax = axes[0]

    bins = np.linspace(0, 1, 51)

    # Non-truth (germline/noise)
    ax.hist(nontruth_posts, bins=bins, alpha=0.6, color="#e41a1c",
            label=f"Non-truth (n={len(nontruth_posts)})", density=True)

    # Truth (somatic)
    ax.hist(truth_posts, bins=bins, alpha=0.7, color="#377eb8",
            label=f"Truth somatic (n={len(truth_posts)})", density=True)

    # Add SLT threshold lines
    for thresh, label in [(0.2, "SLT-C"), (0.5, "SLT-B"), (0.8, "SLT-A")]:
        ax.axvline(x=thresh, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.text(thresh + 0.01, ax.get_ylim()[1] * 0.95, label,
                fontsize=7, color="gray", va="top")

    ax.set_xlabel("PureCN POSTERIOR.SOMATIC", fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title("(A) PureCN Posterior Distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8, loc="upper center")
    ax.set_xlim(-0.02, 1.02)

    # ─── Panel B: Evidence-layer rescue stacked bar ────────────────────
    ax = axes[1]

    categories = ["Post ≥0.5\n(PureCN path)", "Post <0.5\n(evidence rescue)", "FN\n(SLT-D)"]
    values = [tp_post_high, tp_evidence_only, fn_count]
    colors = ["#377eb8", "#4daf4a", "#e41a1c"]

    bars = ax.bar(categories, values, color=colors, edgecolor="black", linewidth=0.8)

    # Add count and percentage labels
    total_truth = total_tp + fn_count
    for bar, val in zip(bars, values):
        pct = val / total_truth * 100 if total_truth > 0 else 0
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                f"{val}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=9)

    ax.set_ylabel("Number of Truth Variants", fontsize=11)
    ax.set_title("(B) Classification Pathway for True Somatics", fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(values) * 1.25)

    # Add horizontal line for total
    ax.axhline(y=total_tp, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.text(2.4, total_tp + 2, f"Total TP={total_tp}", fontsize=7, color="gray", va="bottom")

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
