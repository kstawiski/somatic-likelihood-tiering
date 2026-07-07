# Somatic Likelihood Tiering (SLT)

An interpretable, open-source post-calling triage framework for review-prioritizing tumor-only whole-exome sequencing (WES) variants.

Current release: **v2.1.0**, aligned with the BIB R1 callability-aware SLT protocol.

## Overview

SLT ranks tumor-only WES variants into four review-priority tiers (SLT-A through SLT-D) using:

- **Four complementary evidence layers**: population allele frequency (POPAF), callability-aware gnomAD rarity, germline quality (GERMQ), and COSMIC recurrence
- **PureCN posterior somatic probabilities**: Bayesian posterior from copy number-aware classification
- **Integrated CHIP detection**: Clonal hematopoiesis of indeterminate potential flagging with curated gene lists

All thresholds are deterministic, convention-grounded, frozen, and fully auditable.

## Key Features

- **Zero dependencies**: Core classifier uses only the Python standard library (Python 3.7+)
- **Interpretable**: Every tier assignment is traceable to specific evidence layers
- **Graduated triage**: SLT-A/B high-priority SNV-calibrated queues, SLT-C conservative catchment, and SLT-D lowest-priority review
- **CHIP-aware**: Integrated clonal hematopoiesis detection prevents blood-derived variants from being promoted into the highest review-priority tiers
- **Caller-agnostic**: Works with Mutect2, VarDict, or any TSV/MAF-producing pipeline
- **Callability-aware gnomAD**: missing or allele-unmatched gnomAD evidence is `unevaluable`, not rarity-positive
- **Annotation-only mode**: Operates with gnomAD state + COSMIC when PureCN or BAM-level annotations are unavailable; only common/COSMIC-negative variants route to SLT-D

## Installation

### Standalone (no installation needed)

```bash
# Download and run directly - no pip install required
python3 slt_classify.py --input variants.tsv --output classified.tsv
```

### pip install

```bash
pip install .
# Then use the CLI:
slt-classify --input variants.tsv --output classified.tsv
```

### Development

```bash
pip install -e ".[test]"
pytest tests/ -v
```

## Quick Start

### From a MAF file (most common)

```bash
# 1. Convert MAF to SLT input format
python maf_to_slt_input.py --input somatic.maf --output slt_input.tsv

# 2. Run SLT (annotation-only mode if no PureCN)
python slt_classify.py --input slt_input.tsv --output classified.tsv --annotation-only

# Or with PureCN posteriors for full mode:
python maf_to_slt_input.py --input somatic.maf --purecn purecn_calls.csv --output slt_input.tsv
python slt_classify.py --input slt_input.tsv --output classified.tsv
```

### From a VCF file

```bash
# 1. Convert annotated VCF to SLT input format
python vcf_to_slt_input.py --input annotated.vcf.gz --output slt_input.tsv

# 2. Run SLT
python slt_classify.py --input slt_input.tsv --output classified.tsv
```

### Direct usage (pre-formatted TSV)

```bash
# Full mode (Mutect2 + PureCN)
python slt_classify.py --input annotated_variants.tsv --output slt_output.tsv

# Annotation-only mode (MAF-level, no PureCN)
python slt_classify.py --input variants.maf --output slt_output.tsv --annotation-only
```

### Python API

```python
from slt_classify import classify_variant

variant = {
    "POSTERIOR.SOMATIC": "0.95",
    "POPAF": "40",
    "gnomAD_AF": "0.00001",
    "gnomAD_AN": "20000",
    "GERMQ": "50",
    "COSMIC_CONFIRMED_SOMATIC": "20",
    "Hugo_Symbol": "FGFR3",
    "t_alt_freq": "0.30",
    "HGVSp_Short": "S249C",
    "Variant_Classification": "Missense_Mutation",
}

result = classify_variant(variant)
print(result["slt_tier"])           # SLT-A
print(result["slt_evidence_level"]) # high
print(result["slt_chip_status"])    # no_chip
```

## Preparing Input from Common Formats

SLT accepts a simple tab-separated input, but most variant callers produce VCF or MAF files. Two helper scripts are provided to convert these formats.

### From MAF files (`maf_to_slt_input.py`)

Converts standard MAF (Mutation Annotation Format) files to SLT input. Handles column naming conventions from Funcotator, Oncotator, maf2maf, and cBioPortal.

```bash
# Basic conversion (annotation-only mode — no PureCN)
python maf_to_slt_input.py --input somatic.maf --output slt_input.tsv

# With PureCN posteriors for full-mode SLT
python maf_to_slt_input.py --input somatic.maf --purecn purecn_calls.csv --output slt_input.tsv
```

**Features:**
- Auto-detects column names (Hugo_Symbol vs GENE vs gene, etc.)
- Computes VAF from t_ref_count/t_alt_count if no direct VAF column
- Merges PureCN posteriors, POPAF, and GERMQ by genomic coordinates
- Extracts COSMIC counts from existing annotation fields

### From VCF files (`vcf_to_slt_input.py`)

Converts annotated VCF files (Funcotator or VEP) to SLT input. The VCF must be annotated before conversion — see [Reference Files](#reference-files-for-annotation) below.

```bash
# From Funcotator-annotated VCF
python vcf_to_slt_input.py --input funcotator_annotated.vcf --output slt_input.tsv

# From VEP-annotated VCF
python vcf_to_slt_input.py --input vep_annotated.vcf.gz --output slt_input.tsv

# With PureCN posteriors
python vcf_to_slt_input.py --input annotated.vcf --purecn purecn_calls.csv --output slt_input.tsv

# Specify sample (for multi-sample VCFs)
python vcf_to_slt_input.py --input annotated.vcf --sample TUMOR --output slt_input.tsv
```

**Features:**
- Parses Funcotator FUNCOTATION and VEP CSQ INFO fields
- Extracts Mutect2 POPAF and GERMQ from INFO
- Computes VAF from FORMAT/AD fields
- Handles gzipped VCFs (.vcf.gz)

## Reference Files for Annotation

SLT itself has no reference file dependencies — it operates on pre-annotated TSV input. However, the upstream annotation pipeline requires several reference databases. Below is the recommended annotation workflow.

### Required for Full Mode

| Resource | Version | Purpose | Download |
|----------|---------|---------|----------|
| **Reference genome** | GRCh38 / hg38 | Alignment, variant calling | [GATK bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) |
| **Funcotator data sources** | v1.8 (hg38) | Gene annotation, gnomAD AF, COSMIC | `gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download` |
| **gnomAD** | v4.1 exomes | Population allele frequencies (Layer 2) | [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) |
| **COSMIC** | v103+ | Somatic recurrence (Layer 4) | [COSMIC](https://cancer.sanger.ac.uk/cosmic/download) (registration required) |
| **PureCN panel of normals** | Project-specific | Posterior somatic probability | Built from matched normals or unmatched panel |

### Deployment Guidelines

| Setting | Recommended Mode | Maximum Tier | Notes |
|---------|-----------------|--------------|-------|
| Process-matched NormalDB available | Full SLT | SLT-A | Best performance; requires BAM files + matched normals |
| No NormalDB / no BAM access | Annotation-only | SLT-C | Uses only gnomAD + COSMIC; `--annotation-only` flag |
| Matched normal tissue available | Paired T/N analysis | N/A | Preferred over any tumor-only approach |

### Recommended Annotation Pipeline

```bash
# 1. Variant calling (Mutect2 tumor-only)
gatk Mutect2 \
    -R reference.fa \
    -I tumor.bam \
    --germline-resource af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals pon.vcf.gz \
    -O raw.vcf.gz

# 2. Filter variants
gatk FilterMutectCalls \
    -R reference.fa \
    -V raw.vcf.gz \
    --contamination-table contamination.table \
    -O filtered.vcf.gz

# 3. Annotate with Funcotator (provides gene, consequence, gnomAD, COSMIC)
gatk Funcotator \
    -R reference.fa \
    -V filtered.vcf.gz \
    --ref-version hg38 \
    --data-sources-path funcotator_dataSources.v1.8.hg38.20230908s \
    -O annotated.vcf \
    --output-file-format VCF

# 4. Run PureCN (provides posterior somatic probability)
# See PureCN documentation: https://bioconductor.org/packages/PureCN/
Rscript PureCN.R \
    --sampleid SAMPLE \
    --tumor filtered.vcf.gz \
    --normaldb normalDB.rds \
    --intervals targets_intervals.txt \
    --genome hg38

# 5. Convert to SLT input and classify
python vcf_to_slt_input.py \
    --input annotated.vcf \
    --purecn PureCN_output/SAMPLE.csv \
    --output slt_input.tsv

python slt_classify.py --input slt_input.tsv --output slt_classified.tsv
```

### Annotation-Only Mode (Minimal Requirements)

For MAF-level reanalysis without BAM files, gnomAD evaluability/state and COSMIC annotations are needed:

| Resource | Purpose | SLT Layer |
|----------|---------|-----------|
| **gnomAD AF/AN or state** | Population frequency plus evaluability | Layer 2 |
| **COSMIC** | Somatic recurrence database | Layer 4 |

Layers 1 (POPAF) and 3 (GERMQ) require Mutect2 BAM-level output and are disabled in annotation-only mode. PureCN posterior is unavailable. Maximum achievable tier is SLT-C. Missing, allele-unmatched, or insufficient-AN gnomAD records are `unevaluable` and are not counted as rarity evidence.

### Alternative Annotation with VEP

```bash
# Ensembl VEP (alternative to Funcotator)
vep --input_file filtered.vcf.gz \
    --output_file annotated.vcf \
    --format vcf --vcf \
    --assembly GRCh38 \
    --cache --dir_cache /path/to/vep_cache \
    --fasta reference.fa \
    --everything --pick \
    --plugin COSMIC,/path/to/CosmicCodingMuts.vcf.gz \
    --custom gnomAD_exomes.vcf.gz,gnomADe,vcf,exact,0,AF
```

## Input Format

Tab-separated file with the following columns (flexible naming):

| Field | Alternative Names | Description |
|-------|-------------------|-------------|
| `POSTERIOR.SOMATIC` | `POSTERIOR_SOMATIC` | PureCN posterior somatic probability |
| `POPAF` | | Mutect2 negative log10 population allele frequency |
| `gnomAD_AF` | `gnomAD_exome_AF` | gnomAD exome allele frequency |
| `gnomAD_AN` | `gnomAD_exome_AN` | gnomAD allele number / evaluability support |
| `gnomAD_state` | `gnomad_state` | Optional explicit state: `rare_callable`, `common`, or `unevaluable` |
| `GERMQ` | | Mutect2 germline quality score |
| `COSMIC_CONFIRMED_SOMATIC` | | COSMIC confirmed somatic count |
| `COSMIC_HOTSPOT` | | COSMIC hotspot count |
| `COSMIC_SAMPLE` | `COSMIC_TOTAL_OCC` | COSMIC sample count |
| `CGC_GENE` | | Cancer Gene Census membership (TRUE/FALSE) |
| `Hugo_Symbol` | `GENE`, `gene` | Gene symbol |
| `t_alt_freq` | `AF`, `VAF` | Variant allele frequency |
| `HGVSp_Short` | `Protein_Change`, `AAChange` | Protein change |
| `Variant_Classification` | `Consequence` | Variant consequence |

Missing fields are handled conservatively. Full mode requires a PureCN posterior column in the schema and fails by default if gnomAD state is unevaluable; use `--allow-unevaluable-gnomad` only for a deliberately conservative analysis. Annotation-only mode treats unevaluable gnomAD as cannot-rule-out SLT-C retention, not as positive rarity evidence.

## Output

The classifier appends audit columns to each input row:

| Column | Type | Description |
|--------|------|-------------|
| `slt_layer1_popaf` | bool | POPAF >= 5.0 |
| `slt_gnomad_state` | str | `rare_callable`, `common`, or `unevaluable` |
| `slt_gnomad_state_reason` | str | reason for gnomAD state assignment |
| `slt_layer2_gnomad` | bool | true only for `rare_callable` |
| `slt_layer3_germq` | bool | GERMQ >= 30 |
| `slt_layer4_cosmic` | bool | COSMIC composite criterion |
| `slt_n_somatic_layers` | int | Number of passing layers (0-4) |
| `slt_evidence_level` | str | high / medium / low |
| `slt_chip_status` | str | chip_likely / chip_possible / no_chip |
| `slt_posterior_used` | float | PureCN posterior (or "NA") |
| `slt_tier` | str | SLT-A / SLT-B / SLT-C / SLT-D |
| `slt_mode` | str | full / annotation_only |

## Tier Definitions

### Review-Priority Cascade

Conditions are evaluated in order; the first match determines the tier:

| Tier | Meaning | Conditions |
|------|---------|------------|
| **SLT-A** | Highest review priority | Posterior >= 0.8 AND evidence in {high, medium} AND not chip_likely |
| **SLT-B** | Second review priority | Posterior >= 0.5, OR (high evidence AND >= 3 layers); chip_likely blocked |
| **SLT-C** | Conservative catchment tier | Posterior >= 0.2, OR evidence in {high, medium} |
| **SLT-D** | Lowest review priority | All remaining variants |

### Clinical Workflow Safety Note

SLT is a triage queue for tumor-only review. It does not replace matched-normal sequencing, orthogonal confirmation, germline-risk workflows, tumor board adjudication, or a validated somatic caller. Do not use an SLT tier alone to report somatic status, include a variant in tumor mutational burden, or select targeted therapy.

### Suggested Review Use

| Tier | Recommended Action |
|------|-------------------|
| **SLT-A** | Review first; require matched-normal, orthogonal, or expert adjudication before clinical reporting or TMB use |
| **SLT-B** | Review after SLT-A; prioritize clinically relevant genes for confirmation |
| **SLT-C** | Preserve as an extended cannot-rule-out review queue when clinically relevant |
| **SLT-D** | Lowest-priority queue; revisit only for specific clinical or research hypotheses |
| **CHIP-likely** | Treat as a safety flag; route through hematology/germline-aware review when clinically indicated |

### Evidence Levels

| Level | Criteria |
|-------|----------|
| **High** | >= 3 somatic-supporting layers AND COSMIC layer passes |
| **Medium** | >= 2 somatic-supporting layers |
| **Low** | < 2 somatic-supporting layers |

### CHIP Classification

Variants in CHIP-associated genes are evaluated for clonal hematopoiesis:

- **Tier 1 genes** (13 canonical CHIP drivers): DNMT3A, TET2, ASXL1, PPM1D, JAK2, SF3B1, SRSF2, U2AF1, IDH1, IDH2, ZBTB33, GNB1, CBL
- **Tier 2 genes** (28 extended): TP53, KRAS, NRAS, FLT3, KIT, NPM1, and others

`chip_likely` status blocks both SLT-A and SLT-B assignment, preventing CHIP-context variants from reaching the highest review-priority tiers. `chip_likely` variants are assigned to SLT-C (maximum) or SLT-D. Tier 2 genes (including TP53, KRAS, NRAS) receive `chip_possible` annotation but are NOT blocked from SLT-A/B.

### Annotation-Only Mode

When PureCN is unavailable (e.g., MAF-level reanalysis without BAM files):

- Layers 1 (POPAF) and 3 (GERMQ) are disabled
- PureCN posterior is set to null
- Maximum achievable tier is SLT-C
- SLT-A and SLT-B are unreachable
- Rare-callable/COSMIC-positive, rare-callable/COSMIC-negative, common/COSMIC-positive, unevaluable/COSMIC-positive, and unevaluable/COSMIC-negative profiles route to capped SLT-C
- Only common/COSMIC-negative profiles route to SLT-D

Use the `--annotation-only` flag (or `--degraded` for backward compatibility).

## Thresholds

All thresholds are convention-grounded, frozen, and deterministic:

| Parameter | Value | Source |
|-----------|-------|--------|
| Posterior SLT-A gate | >= 0.8 | PureCN recommended (Riester et al., 2016) |
| Posterior SLT-B gate | >= 0.5 | PureCN recommended |
| Posterior SLT-C gate | >= 0.2 | PureCN recommended |
| POPAF threshold | >= 5.0 | Mutect2 standard (AF <= 1e-5) |
| gnomAD AF threshold | < 0.001 with adequate AN or explicit state | ACMG/AMP BA1/BS1 aligned, callability-aware |
| GERMQ threshold | >= 30 | Mutect2 standard (<=0.1% germline probability) |
| COSMIC confirmed min | >= 5 | Conservative recurrence threshold |
| COSMIC hotspot min | >= 10 | Conservative recurrence threshold |
| CGC + COSMIC sample min | >= 2 | Conservative recurrence threshold |

## Benchmark Performance

### SEQC2 HCC1395 (primary benchmark)

Validated on the SEQC2 HCC1395 breast cancer truth set (455 variants in evaluation regions), processed in tumor-only mode:

| Threshold | Called | TP | Sensitivity | PPV | NNR |
|-----------|-------:|---:|------------:|----:|----:|
| >= SLT-A | 101 | 78 | 17.1% | 77.2% | 1.29 |
| >= SLT-B | 157 | 99 | 21.8% | 63.1% | 1.59 |
| >= SLT-C | 2,246 | 352 | 77.4% | 15.7% | 6.38 |
| >= SLT-D | 4,471 | 430 | 94.5% | 9.6% | 10.40 |

The R1 callability-aware correction narrowed SLT-C relative to the original release because missing gnomAD annotation no longer counts as rarity evidence. All tiers together recover the frozen Mutect2-detected truth ceiling (430/455).

### External validation (BostonGene cell lines)

| Sample | Tumor type | Truth frame | SLT-A/B recovered | SLT-A/B/C recovered | Paired tumor-normal Mutect2 PASS recovered |
|--------|-----------|------------:|-------------------:|---------------------:|------------------------------------------:|
| COLO829 | Melanoma | 445 | 37/445 | 125/445 | 264/445 |
| NCI-H1770 | NSCLC | 1,203 | 151/1,203 | 459/1,203 | 817/1,203 |

These rows position tumor-only SLT as post-calling triage on the local candidate surface, not as a replacement for paired tumor-normal PASS calling.

### Clinical validation (HdM-BLCA-1 bladder cancer)

Validated on 27 FFPE metastatic urothelial carcinoma patients (Boll et al. 2023, *Sci Rep*; EGA: EGAS00001007086). Reference standard: matched tumor-normal Mutect2 PASS variants (TLOD >= 20, VAF >= 5%). PureCN v2.12.0 with FFPE-adapted min.base.quality = 20.

| Threshold | Recall [patient-level 95% CI] | Concordance | NNR |
|-----------|---------------------|-----|-----|
| >= SLT-A | 18.2% [14.0-23.5%] | 6.8% | 14.6 |
| >= SLT-B | 31.8% [25.5-38.7%] | 4.8% | 20.9 |
| >= SLT-C | 72.3% [68.0-76.7%] | 0.65% | 152.8 |

Clinical rows are concordance/recall analyses versus partially dependent matched-normal references, not independent sensitivity validation. SLT reduces first-pass candidate counts but does not measure review time or create a tumor-board-ready reportable list by itself.

**FFPE deployment note:** Set PureCN `--min-base-quality 20` (or lower) for FFPE samples. SLT-A/B precision tiers require PureCN purity >= ~0.25.

## Tests

```bash
# Install test dependencies
pip install pytest

# Run tests
pytest tests/ -v
```

The test suite covers all evidence layers, CHIP classification, tier assignment cascade, annotation-only mode, boundary conditions, and integration tests.

## Repository Structure

```
somatic-likelihood-tiering/
├── slt_classify.py          # Core classifier (standalone, no dependencies)
├── maf_to_slt_input.py      # MAF → SLT input converter
├── vcf_to_slt_input.py      # VCF → SLT input converter (Funcotator/VEP)
├── tests/
│   └── test_slt_classify.py # Unit test suite (39+ tests)
├── examples/
│   ├── example_input.tsv    # Example input
│   └── example_output.tsv   # Example output
├── benchmark/               # Benchmark evaluation scripts
│   └── README.md
├── pyproject.toml           # Python packaging
├── LICENSE                  # MIT License
└── README.md                # This file
```


## License

MIT License. See [LICENSE](LICENSE) for details.
