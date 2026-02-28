#!/usr/bin/env python3
"""Convert a MAF file to SLT input format.

Reads a standard MAF (Mutation Annotation Format) file and produces a
tab-separated file compatible with slt_classify.py. Handles common MAF
column naming conventions from Funcotator, Oncotator, maf2maf, and
cBioPortal formats.

No external dependencies required (Python 3.7+ standard library only).

Usage:
    python maf_to_slt_input.py --input somatic.maf --output slt_input.tsv
    python maf_to_slt_input.py --input somatic.maf --output slt_input.tsv --purecn purecn_calls.csv
"""

import argparse
import csv
import sys
from typing import Dict, Optional


# Standard SLT input column names
SLT_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "HGVSp_Short",
    "t_alt_freq",
    "POPAF",
    "GERMQ",
    "gnomAD_AF",
    "COSMIC_CONFIRMED_SOMATIC",
    "COSMIC_HOTSPOT",
    "COSMIC_SAMPLE",
    "CGC_GENE",
    "POSTERIOR.SOMATIC",
    "Tumor_Sample_Barcode",
]


# Column name aliases: SLT name -> list of MAF alternatives
COLUMN_ALIASES = {
    "Hugo_Symbol": ["Hugo_Symbol", "gene", "GENE", "Gene", "SYMBOL"],
    "Chromosome": ["Chromosome", "chr", "Chrom", "CHROM", "#CHROM"],
    "Start_Position": ["Start_Position", "Start_position", "start", "POS", "Position"],
    "Reference_Allele": ["Reference_Allele", "ref", "REF", "Ref_Allele"],
    "Tumor_Seq_Allele2": ["Tumor_Seq_Allele2", "alt", "ALT", "Alt_Allele", "Tumor_Allele"],
    "Variant_Classification": ["Variant_Classification", "Consequence", "consequence",
                                "Variant_Type", "Effect"],
    "HGVSp_Short": ["HGVSp_Short", "Protein_Change", "protein_change", "AAChange",
                     "HGVSp", "amino_acid_change"],
    "t_alt_freq": ["t_alt_freq", "AF", "VAF", "tumor_f", "Tumor_Freq",
                   "allelic_fraction", "i_tumor_f"],
    "POPAF": ["POPAF", "popaf"],
    "GERMQ": ["GERMQ", "germq"],
    "gnomAD_AF": ["gnomAD_AF", "gnomAD_exome_AF", "gnomAD_genome_AF",
                   "gnomad_AF", "gnomad_exome_AF", "ExAC_AF", "exac_af"],
    "COSMIC_CONFIRMED_SOMATIC": ["COSMIC_CONFIRMED_SOMATIC", "cosmic_confirmed_somatic"],
    "COSMIC_HOTSPOT": ["COSMIC_HOTSPOT", "cosmic_hotspot"],
    "COSMIC_SAMPLE": ["COSMIC_SAMPLE", "COSMIC_TOTAL_OCC", "cosmic_total_occ",
                       "COSMIC_COUNT", "cosmic_count"],
    "CGC_GENE": ["CGC_GENE", "cgc_gene", "Cancer_Gene_Census"],
    "POSTERIOR.SOMATIC": ["POSTERIOR.SOMATIC", "POSTERIOR_SOMATIC",
                          "posterior_somatic", "posterior.somatic"],
    "Tumor_Sample_Barcode": ["Tumor_Sample_Barcode", "sample_id", "Sample",
                              "SAMPLE", "sample"],
}


def resolve_column(maf_header: list, slt_col: str) -> Optional[str]:
    """Find the MAF column name that maps to an SLT column."""
    aliases = COLUMN_ALIASES.get(slt_col, [slt_col])
    for alias in aliases:
        if alias in maf_header:
            return alias
    return None


def compute_vaf(row: dict, maf_header: list) -> Optional[str]:
    """Compute VAF from t_ref_count and t_alt_count if direct VAF column is missing."""
    ref_count_cols = ["t_ref_count", "n_ref_count_tumor", "ref_count"]
    alt_count_cols = ["t_alt_count", "n_alt_count_tumor", "alt_count"]

    ref_col = None
    alt_col = None
    for c in ref_count_cols:
        if c in maf_header:
            ref_col = c
            break
    for c in alt_count_cols:
        if c in maf_header:
            alt_col = c
            break

    if ref_col and alt_col:
        try:
            ref = int(row.get(ref_col, 0))
            alt = int(row.get(alt_col, 0))
            total = ref + alt
            if total > 0:
                return f"{alt / total:.6f}"
        except (ValueError, TypeError):
            pass
    return None


def parse_cosmic_from_existing_ids(row: dict, maf_header: list) -> dict:
    """Extract COSMIC occurrence counts from existing COSMIC ID fields if dedicated columns are missing."""
    # Some MAFs have COSMIC_ID or similar with format like "COSM1234|5"
    # This is a best-effort extraction
    cosmic_fields = {}

    # Check for COSMIC_TOTAL_OCC or similar
    for col in ["COSMIC_TOTAL_OCC", "cosmic_total_occ", "COSMIC_n_overlapping_mutations"]:
        if col in maf_header and row.get(col):
            val = row[col]
            if val not in ("", ".", "NA", "nan", "None"):
                cosmic_fields["COSMIC_SAMPLE"] = val
                break

    return cosmic_fields


def load_purecn(purecn_path: str) -> Dict[str, dict]:
    """Load PureCN output and index by chr:pos:ref:alt."""
    purecn_data = {}
    with open(purecn_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = row.get("chr", "")
            pos = row.get("start", "")
            ref = row.get("ref", "")
            alt = row.get("alt", "")
            key = f"{chrom}:{pos}:{ref}:{alt}"
            purecn_data[key] = {
                "POSTERIOR.SOMATIC": row.get("POSTERIOR.SOMATIC", ""),
                "POPAF": row.get("POPAF", ""),
                "GERMQ": row.get("GERMQ", ""),
            }
    return purecn_data


def main():
    parser = argparse.ArgumentParser(
        description="Convert MAF to SLT input format",
        epilog="Produces a TSV file ready for slt_classify.py",
    )
    parser.add_argument("--input", "-i", required=True, help="Input MAF file")
    parser.add_argument("--output", "-o", required=True, help="Output SLT-ready TSV")
    parser.add_argument("--purecn", help="Optional PureCN output CSV to merge posteriors/POPAF/GERMQ")
    args = parser.parse_args()

    # Load PureCN data if provided
    purecn_data = {}
    if args.purecn:
        print(f"Loading PureCN data from {args.purecn}...", file=sys.stderr)
        purecn_data = load_purecn(args.purecn)
        print(f"  {len(purecn_data)} variants loaded", file=sys.stderr)

    with open(args.input, "r") as fin:
        # Skip MAF comment lines (start with #)
        header_line = None
        for line in fin:
            if not line.startswith("#"):
                header_line = line.strip()
                break

        if not header_line:
            print("ERROR: No header found in MAF file", file=sys.stderr)
            sys.exit(1)

        maf_header = header_line.split("\t")

        # Build column mapping
        col_map = {}
        mapped = []
        unmapped = []
        for slt_col in SLT_COLUMNS:
            maf_col = resolve_column(maf_header, slt_col)
            if maf_col:
                col_map[slt_col] = maf_col
                mapped.append(f"  {slt_col} <- {maf_col}")
            else:
                unmapped.append(slt_col)

        print("Column mapping:", file=sys.stderr)
        for m in mapped:
            print(m, file=sys.stderr)
        if unmapped:
            print(f"  Unmapped (will be empty): {', '.join(unmapped)}", file=sys.stderr)

        # Check VAF availability
        vaf_col = col_map.get("t_alt_freq")
        vaf_from_counts = vaf_col is None
        if vaf_from_counts:
            print("  Note: No direct VAF column found; will compute from t_ref_count/t_alt_count if available",
                  file=sys.stderr)

        reader = csv.DictReader(fin, fieldnames=maf_header, delimiter="\t")

        with open(args.output, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=SLT_COLUMNS, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()

            n = 0
            n_purecn_matched = 0
            for row in reader:
                out = {}
                for slt_col in SLT_COLUMNS:
                    maf_col = col_map.get(slt_col)
                    if maf_col:
                        val = row.get(maf_col, "")
                        out[slt_col] = val if val not in (".", "NA", "None", "nan") else ""
                    else:
                        out[slt_col] = ""

                # Compute VAF from counts if needed
                if vaf_from_counts and not out.get("t_alt_freq"):
                    computed_vaf = compute_vaf(row, maf_header)
                    if computed_vaf:
                        out["t_alt_freq"] = computed_vaf

                # Extract COSMIC from existing fields if dedicated columns are empty
                if not out.get("COSMIC_SAMPLE"):
                    cosmic_extras = parse_cosmic_from_existing_ids(row, maf_header)
                    for k, v in cosmic_extras.items():
                        if k in SLT_COLUMNS:
                            out[k] = v

                # Merge PureCN data if available
                if purecn_data:
                    chrom = out.get("Chromosome", "")
                    pos = out.get("Start_Position", "")
                    ref = out.get("Reference_Allele", "")
                    alt = out.get("Tumor_Seq_Allele2", "")
                    key = f"{chrom}:{pos}:{ref}:{alt}"
                    if key in purecn_data:
                        pc = purecn_data[key]
                        for field in ["POSTERIOR.SOMATIC", "POPAF", "GERMQ"]:
                            if pc.get(field) and not out.get(field):
                                out[field] = pc[field]
                        n_purecn_matched += 1

                writer.writerow(out)
                n += 1

    print(f"\nConverted {n} variants to SLT input format", file=sys.stderr)
    if purecn_data:
        print(f"  PureCN matched: {n_purecn_matched}/{n}", file=sys.stderr)
    print(f"  Output: {args.output}", file=sys.stderr)
    print(f"\nNext step: python slt_classify.py --input {args.output} --output classified.tsv", file=sys.stderr)
    if not purecn_data and not col_map.get("POSTERIOR.SOMATIC"):
        print("  Tip: No PureCN data detected. Use --degraded mode:", file=sys.stderr)
        print(f"    python slt_classify.py --input {args.output} --output classified.tsv --degraded",
              file=sys.stderr)


if __name__ == "__main__":
    main()
