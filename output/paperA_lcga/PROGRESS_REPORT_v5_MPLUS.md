# Paper A v5 Progress Report: Mplus Replication

**Date:** 2026-05-08  
**Scope:** Mplus 7 replication of the primary 9-item LCGA model.

## What Was Added

Mplus 7 was installed and available at:

- `C:/Program Files/Mplus/Mplus.exe`

Because Mplus 7 can be brittle with non-ASCII paths, the run directory was set to:

- `D:/mplus_kcyps_paperA`

Project copies of the inputs, outputs, parsed summaries, and report were saved to:

- `output/paperA_lcga/primary_9item/mplus/`

## Mplus Model Specification

The Mplus model was matched to the Python/R primary LCGA:

- 9-item depression composite, excluding E04 suicidal ideation.
- Analytic sample `N=2,452`, requiring at least 3 valid waves.
- Quadratic growth factors with time scores `-3,-2,-1,0,1,2,3`.
- LCGA constraints: `i@0; s@0; q@0;`.
- One residual variance constrained equal across all seven waves and classes.
- `STARTS = 1000 250`; K=1-5 all terminated normally and replicated the best loglikelihood.

## Fit Replication

The full Mplus fit table includes number of free parameters, log-likelihood, AIC, BIC, sample-size adjusted BIC (aBIC), and entropy.

| K | Free parameters | Mplus logL | Python logL | Mplus AIC | Python AIC | Mplus BIC | Python BIC | Mplus aBIC | Python aBIC | Entropy Mplus | Entropy Python |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4 | -14143.189 | -14143.188 | 28294.378 | 28294.376 | 28317.597 | 28317.595 | 28304.888 | 28304.886 | NA | 1.000 |
| 2 | 8 | -13055.844 | -13055.844 | 26127.689 | 26127.687 | 26174.126 | 26174.124 | 26148.708 | 26148.706 | .675 | .674 |
| 3 | 12 | -12746.201 | -12746.200 | 25516.402 | 25516.400 | 25586.058 | 25586.056 | 25547.931 | 25547.929 | .719 | .719 |
| 4 | 16 | -12551.601 | -12551.600 | 25135.202 | 25135.200 | 25228.076 | 25228.075 | 25177.240 | 25177.239 | .661 | .661 |
| 5 | 20 | -12485.230 | -12485.230 | 25010.461 | 25010.459 | 25126.554 | 25126.553 | 25063.009 | 25063.008 | .667 | .667 |

Conclusion: Mplus reproduces the Python EM and R lcmm fit indices to rounding precision.

## K=4 Label Alignment

Mplus raw classes differ only by arbitrary label order:

| Manuscript label | Mplus raw class | Python n | Mplus n | Difference |
|---|---:|---:|---:|---:|
| C1 Stable-Low | 2 | 615 | 616 | +1 |
| C2 Persistent-High | 4 | 245 | 245 | 0 |
| C3 Decreasing | 1 | 726 | 725 | -1 |
| C4 Late-Increasing | 3 | 866 | 866 | 0 |

After label alignment:

- Percent agreement = 99.96%.
- Cohen's kappa = .9994.

## Files

Generated scripts:

- `scripts/33_prepare_mplus_lcga_9item.py`
- `scripts/34_parse_mplus_lcga_9item.py`

Primary Mplus outputs:

- `primary_9item/mplus/MPLUS_REPLICATION_REPORT.md`
- `primary_9item/mplus/mplus_fit_summary_9item.csv`
- `primary_9item/mplus/mplus_vs_python_fit_summary_9item.csv`
- `primary_9item/mplus/mplus_class_profiles_K4_9item.csv`
- `primary_9item/mplus/python_vs_mplus_class_comparison_9item.csv`
- `primary_9item/mplus/lcga9_k1.out` through `lcga9_k5.out`

## Manuscript Update

`MANUSCRIPT_DRAFT.md` was updated from v4.2 to v4.3:

- Methods now states independent R lcmm and Mplus cross-validation.
- Results now reports Mplus K=4 BIC/logL replication and modal agreement.
- Section 3.2 now includes Mplus raw-label alignment.
- Sensitivity analyses now mention R lcmm and Mplus together.
- Output tree and script list now include the Mplus replication files.

## Bottom Line

The primary 9-item LCGA is now replicated across three independent implementations:

1. Python EM
2. R lcmm
3. Mplus 7

This substantially strengthens reviewer defense for JAD/JYA-level submission.
