#!/usr/bin/env python3
"""
Error Analysis — SLT False Positives and False Negatives

Analyzes:
  1. SLT-A/B False Positives: What non-somatic variants pass the high-confidence tiers?
  2. False Negatives (truth variants missed at ≥SLT-C): Why does SLT miss some true somatics?
  3. PureCN posterior failure modes: Where does the continuous classifier disagree with SLT tiers?

Usage:
    python3 error_analysis.py \
        --slt-classified  work/seqc2/slt/slt_classified.tsv \
        --evaluation      work/seqc2/slt/evaluation_per_variant.tsv \
        --output-dir      results/error_analysis/
"""

import argparse
import csv
import os
from collections import Counter, defaultdict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--slt-classified", required=True)
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load evaluation per-variant data
    eval_data = {}
    with open(args.evaluation) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
            eval_data[key] = row

    # Load SLT classified data
    slt_data = []
    with open(args.slt_classified) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            parts = row["variant_key"].split(":")
            if len(parts) >= 4:
                key = ":".join(parts[:4])
            else:
                key = row["variant_key"]
            ev = eval_data.get(key, {})
            row["in_truth"] = ev.get("in_truth", "False")
            row["in_eval_region"] = ev.get("in_eval_region", "False")
            row["_key"] = key
            slt_data.append(row)

    # Restrict to eval region variants (handle both "True"/"1" formats)
    def is_true(val):
        return str(val).lower() in ("true", "1", "yes")

    eval_variants = [r for r in slt_data if is_true(r["in_eval_region"])]

    print(f"Total SLT variants: {len(slt_data)}")
    print(f"In evaluation region: {len(eval_variants)}")

    # --- 1. False Positive Analysis ---
    # FP at ≥SLT-A: variants in eval region with tier A that are NOT in truth
    tier_order = {"SLT-A": 0, "SLT-B": 1, "SLT-C": 2, "SLT-D": 3}
    fp_a = [r for r in eval_variants
            if r["slt_tier"] in ("SLT-A",) and not is_true(r["in_truth"])]
    fp_ab = [r for r in eval_variants
             if r["slt_tier"] in ("SLT-A", "SLT-B") and not is_true(r["in_truth"])]
    fp_abc = [r for r in eval_variants
              if r["slt_tier"] in ("SLT-A", "SLT-B", "SLT-C") and not is_true(r["in_truth"])]

    # Analyze FP characteristics
    def analyze_fp(variants, label):
        lines = []
        lines.append(f"=== {label} False Positive Analysis (n={len(variants)}) ===\n")

        # Gene distribution
        genes = Counter(r.get("Hugo_Symbol", "unknown") for r in variants)
        lines.append(f"Top genes (top 20):")
        for gene, count in genes.most_common(20):
            lines.append(f"  {gene}: {count}")

        # Consequence distribution
        consequences = Counter(r.get("Variant_Classification", "unknown") for r in variants)
        lines.append(f"\nVariant consequences:")
        for cons, count in consequences.most_common():
            lines.append(f"  {cons}: {count}")

        # Posterior distribution
        posteriors = []
        for r in variants:
            try:
                posteriors.append(float(r.get("POSTERIOR.SOMATIC", "0") or "0"))
            except ValueError:
                posteriors.append(0.0)
        if posteriors:
            posteriors.sort()
            lines.append(f"\nPosterior somatic distribution:")
            lines.append(f"  Min: {min(posteriors):.4f}")
            lines.append(f"  Q1:  {posteriors[len(posteriors)//4]:.4f}")
            lines.append(f"  Median: {posteriors[len(posteriors)//2]:.4f}")
            lines.append(f"  Q3:  {posteriors[3*len(posteriors)//4]:.4f}")
            lines.append(f"  Max: {max(posteriors):.4f}")

        # Evidence level distribution
        evidence = Counter(r.get("slt_evidence_level", "unknown") for r in variants)
        lines.append(f"\nEvidence levels:")
        for ev, count in evidence.most_common():
            lines.append(f"  {ev}: {count}")

        # FILTER distribution
        filters = Counter(r.get("FILTER", "unknown") for r in variants)
        lines.append(f"\nMutect2 FILTER:")
        for filt, count in filters.most_common(10):
            lines.append(f"  {filt}: {count}")

        # CHIP status
        chip = Counter(r.get("slt_chip_status", "unknown") for r in variants)
        lines.append(f"\nCHIP status:")
        for c, count in chip.most_common():
            lines.append(f"  {c}: {count}")

        # gnomAD status
        gnomad_present = sum(1 for r in variants
                             if r.get("slt_layer2_gnomad", "0") == "0")
        lines.append(f"\ngnomAD status:")
        lines.append(f"  Layer fails (gnomAD AF >= 0.001): {gnomad_present}")
        lines.append(f"  Layer passes (rare/absent): {len(variants) - gnomad_present}")

        lines.append("")
        return "\n".join(lines)

    fp_report = []
    fp_report.append(analyze_fp(fp_a, "≥SLT-A"))
    fp_report.append(analyze_fp(fp_ab, "≥SLT-B (cumulative A+B)"))
    fp_report.append(analyze_fp(fp_abc, "≥SLT-C (cumulative A+B+C)"))

    fp_path = os.path.join(args.output_dir, "false_positive_analysis.txt")
    with open(fp_path, "w") as f:
        f.write("\n".join(fp_report))
    print(f"Wrote: {fp_path}")

    # Write detailed FP TSV for SLT-A
    fp_a_path = os.path.join(args.output_dir, "fp_slt_a_details.tsv")
    with open(fp_a_path, "w") as f:
        cols = ["_key", "Hugo_Symbol", "Variant_Classification", "FILTER",
                "POSTERIOR.SOMATIC", "slt_evidence_level", "slt_chip_status",
                "slt_layer1_popaf", "slt_layer2_gnomad", "slt_layer3_germq",
                "slt_layer4_cosmic", "t_alt_count", "t_depth"]
        f.write("\t".join(cols) + "\n")
        for r in fp_a:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"Wrote: {fp_a_path}")

    # --- 2. False Negative Analysis ---
    # Truth variants missed at cumulative ≥SLT-C
    tp_keys = set()
    for r in eval_variants:
        if is_true(r["in_truth"]) and r["slt_tier"] in ("SLT-A", "SLT-B", "SLT-C"):
            tp_keys.add(r["_key"])

    # Get all truth variant keys from evaluation data
    truth_keys = set()
    for key, ev in eval_data.items():
        if is_true(ev.get("in_truth")) and is_true(ev.get("in_eval_region")):
            truth_keys.add(key)

    # Missed truth = truth_keys not in tp_keys
    missed_keys = truth_keys - tp_keys

    # Try to find missed variants in SLT data (they may be SLT-D or not called at all)
    slt_by_key = {r["_key"]: r for r in slt_data}

    fn_report = []
    fn_report.append(f"=== False Negative Analysis (≥SLT-C threshold) ===")
    fn_report.append(f"Truth variants in eval region: {len(truth_keys)}")
    fn_report.append(f"True positives at ≥SLT-C: {len(tp_keys)}")
    fn_report.append(f"Missed (false negatives): {len(missed_keys)}")
    fn_report.append("")

    fn_called = []  # FN that were called but classified SLT-D
    fn_not_called = []  # FN not even in SLT output

    for key in sorted(missed_keys):
        if key in slt_by_key:
            fn_called.append(slt_by_key[key])
        else:
            fn_not_called.append(key)

    fn_report.append(f"Called by Mutect2 but SLT-D: {len(fn_called)}")
    fn_report.append(f"Not called by Mutect2: {len(fn_not_called)}")
    fn_report.append("")

    if fn_called:
        fn_report.append("--- FN called but SLT-D (missed by tier classification) ---")
        posteriors = []
        for r in fn_called:
            gene = r.get("Hugo_Symbol", "unknown")
            tier = r.get("slt_tier", "?")
            post = r.get("POSTERIOR.SOMATIC", "NA")
            ev = r.get("slt_evidence_level", "?")
            filt = r.get("FILTER", "?")
            fn_report.append(f"  {r['_key']}  gene={gene}  tier={tier}  "
                           f"post={post}  evidence={ev}  filter={filt}")
            try:
                posteriors.append(float(post or "0"))
            except ValueError:
                pass

        if posteriors:
            fn_report.append(f"\n  Posterior distribution for SLT-D truth variants:")
            posteriors.sort()
            fn_report.append(f"    Min: {min(posteriors):.4f}")
            fn_report.append(f"    Median: {posteriors[len(posteriors)//2]:.4f}")
            fn_report.append(f"    Max: {max(posteriors):.4f}")
            low_post = sum(1 for p in posteriors if p < 0.2)
            fn_report.append(f"    Posterior < 0.2: {low_post}/{len(posteriors)}")

        # Evidence level distribution
        ev_dist = Counter(r.get("slt_evidence_level", "?") for r in fn_called)
        fn_report.append(f"\n  Evidence levels:")
        for ev, count in ev_dist.most_common():
            fn_report.append(f"    {ev}: {count}")

    if fn_not_called:
        fn_report.append(f"\n--- FN not called by Mutect2 ({len(fn_not_called)} variants) ---")
        for key in fn_not_called[:20]:
            fn_report.append(f"  {key}")
        if len(fn_not_called) > 20:
            fn_report.append(f"  ... and {len(fn_not_called) - 20} more")

    fn_path = os.path.join(args.output_dir, "false_negative_analysis.txt")
    with open(fn_path, "w") as f:
        f.write("\n".join(fn_report))
    print(f"Wrote: {fn_path}")

    # --- 3. PureCN Posterior Failure Modes ---
    purecn_report = []
    purecn_report.append("=== PureCN Posterior Failure Mode Analysis ===\n")

    # Variants with high posterior but NOT in truth (FP driven by PureCN)
    high_post_fp = []
    for r in eval_variants:
        if not is_true(r["in_truth"]):
            try:
                post = float(r.get("POSTERIOR.SOMATIC", "0") or "0")
            except ValueError:
                post = 0.0
            if post >= 0.8:
                high_post_fp.append((post, r))

    purecn_report.append(f"High posterior (>=0.8) false positives: {len(high_post_fp)}")
    if high_post_fp:
        high_post_fp.sort(key=lambda x: -x[0])
        genes = Counter(r.get("Hugo_Symbol", "?") for _, r in high_post_fp)
        purecn_report.append(f"Top genes:")
        for gene, count in genes.most_common(10):
            purecn_report.append(f"  {gene}: {count}")

    # Truth variants with low posterior (PureCN disagrees)
    low_post_tp = []
    for r in eval_variants:
        if is_true(r["in_truth"]):
            try:
                post = float(r.get("POSTERIOR.SOMATIC", "0") or "0")
            except ValueError:
                post = 0.0
            if post < 0.2:
                low_post_tp.append((post, r))

    purecn_report.append(f"\nTrue somatics with low posterior (<0.2): {len(low_post_tp)}")
    if low_post_tp:
        low_post_tp.sort(key=lambda x: x[0])
        for post, r in low_post_tp[:10]:
            gene = r.get("Hugo_Symbol", "?")
            tier = r.get("slt_tier", "?")
            ev = r.get("slt_evidence_level", "?")
            purecn_report.append(f"  post={post:.4f}  gene={gene}  tier={tier}  evidence={ev}")

    # Variants with no PureCN posterior at all
    no_post = [r for r in eval_variants
               if not r.get("POSTERIOR.SOMATIC") or r.get("POSTERIOR.SOMATIC") == ""]
    tp_no_post = [r for r in no_post if is_true(r["in_truth"])]
    purecn_report.append(f"\nVariants with no PureCN posterior: {len(no_post)}")
    purecn_report.append(f"  Of which are true positives: {len(tp_no_post)}")

    # How many TPs are rescued by evidence layers alone?
    rescued_by_evidence = []
    for r in eval_variants:
        if is_true(r["in_truth"]):
            try:
                post = float(r.get("POSTERIOR.SOMATIC", "0") or "0")
            except ValueError:
                post = 0.0
            tier = r.get("slt_tier", "SLT-D")
            if post < 0.5 and tier in ("SLT-A", "SLT-B", "SLT-C"):
                rescued_by_evidence.append(r)

    purecn_report.append(f"\nTruth variants rescued by evidence layers (post<0.5, tier≤C): "
                        f"{len(rescued_by_evidence)}")
    if rescued_by_evidence:
        tier_dist = Counter(r["slt_tier"] for r in rescued_by_evidence)
        for t, c in sorted(tier_dist.items()):
            purecn_report.append(f"  {t}: {c}")

    purecn_path = os.path.join(args.output_dir, "purecn_failure_modes.txt")
    with open(purecn_path, "w") as f:
        f.write("\n".join(purecn_report))
    print(f"Wrote: {purecn_path}")

    # --- 4. Summary statistics for manuscript ---
    summary = []
    summary.append("=== Error Analysis Summary (for manuscript) ===\n")
    summary.append(f"FP at ≥SLT-A: {len(fp_a)} (of {len([r for r in eval_variants if r['slt_tier'] == 'SLT-A'])} called)")
    summary.append(f"FP at ≥SLT-B: {len(fp_ab)}")
    summary.append(f"FP at ≥SLT-C: {len(fp_abc)}")
    summary.append(f"FN at ≥SLT-C: {len(missed_keys)}")
    summary.append(f"  - Called but SLT-D: {len(fn_called)}")
    summary.append(f"  - Not called by Mutect2: {len(fn_not_called)}")
    summary.append(f"PureCN high-post FP (>=0.8): {len(high_post_fp)}")
    summary.append(f"True somatics low-post (<0.2): {len(low_post_tp)}")
    summary.append(f"TPs rescued by evidence layers: {len(rescued_by_evidence)}")

    summary_path = os.path.join(args.output_dir, "error_analysis_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary))
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
