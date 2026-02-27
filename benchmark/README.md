# Benchmark Scripts

Scripts used to reproduce the SEQC2 HCC1395 benchmark results reported in the manuscript. These scripts require additional dependencies (numpy, scikit-learn, matplotlib) and external data.

## Scripts

| Script | Dependencies | Description |
|--------|-------------|-------------|
| `evaluate_seqc2.py` | stdlib | Evaluate SLT tiers against SEQC2 truth set |
| `evaluate_competitor.py` | stdlib | Evaluate competitor methods (DeepSomatic, VarDict) |
| `evaluate_baselines.py` | stdlib | Evaluate PureCN-only and gnomAD-only baselines |
| `compute_wilson_ci.py` | stdlib | Add Wilson score confidence intervals |
| `compute_purecn_auroc.py` | numpy, sklearn | PureCN posterior AUROC analysis |
| `head_to_head_table.py` | stdlib | Generate benchmark comparison table |
| `error_analysis.py` | stdlib | False positive/negative analysis |
| `fig_slt_architecture.py` | matplotlib, numpy | SLT framework architecture figure |
| `fig_seqc2_benchmark.py` | matplotlib, numpy | Benchmark results figures |
| `fig_benchmark_comparison.py` | matplotlib, numpy | Method comparison figures |
| `fig_mc3_degraded.py` | matplotlib, numpy | MC3 degraded mode figure |
| `fig_error_analysis.py` | matplotlib, numpy | Error analysis figure |

## Installation

```bash
pip install -e ".[benchmark]"
```

## Data Requirements

These scripts require the SEQC2 HCC1395 truth set and evaluation regions, which are available from the SEQC2 consortium.
