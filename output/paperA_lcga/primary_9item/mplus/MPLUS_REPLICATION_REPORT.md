# Mplus 9-item LCGA Replication Report

**Run date:** 2026-05-08  
**Mplus version:** 7  
**Run directory:** `D:/mplus_kcyps_paperA`  
**Project copy:** `output/paperA_lcga/primary_9item/mplus`

## Model

- 9-item depression composite, suicidal ideation item E04 excluded.
- Analytic N = 2,452, >=3 valid depression waves.
- Quadratic LCGA with time scores -3,-2,-1,0,1,2,3.
- Growth-factor variances fixed to zero (`i@0; s@0; q@0;`).
- One residual variance constrained equal across waves and classes.
- K=1..5 estimated with `STARTS = 1000 250`; every solution terminated normally and replicated the best loglikelihood.

## Fit Replication

Mplus reproduces the Python EM/R lcmm likelihoods and information criteria to rounding precision.

|   K |   free_parameters |   logL_Mplus |   logL_Python |   AIC_Mplus |   AIC_Python |   BIC_Mplus |   BIC_Python |   aBIC_Mplus |   aBIC_Python |   entropy_Mplus |   entropy_Python | best_logL_replicated   |
|----:|------------------:|-------------:|--------------:|------------:|-------------:|------------:|-------------:|-------------:|--------------:|----------------:|-----------------:|:-----------------------|
|   1 |                 4 |     -14143.2 |      -14143.2 |     28294.4 |      28294.4 |     28317.6 |      28317.6 |      28304.9 |       28304.9 |         nan     |           1      | True                   |
|   2 |                 8 |     -13055.8 |      -13055.8 |     26127.7 |      26127.7 |     26174.1 |      26174.1 |      26148.7 |       26148.7 |           0.675 |           0.6745 | True                   |
|   3 |                12 |     -12746.2 |      -12746.2 |     25516.4 |      25516.4 |     25586.1 |      25586.1 |      25547.9 |       25547.9 |           0.719 |           0.7192 | True                   |
|   4 |                16 |     -12551.6 |      -12551.6 |     25135.2 |      25135.2 |     25228.1 |      25228.1 |      25177.2 |       25177.2 |           0.661 |           0.6607 | True                   |
|   5 |                20 |     -12485.2 |      -12485.2 |     25010.5 |      25010.5 |     25126.6 |      25126.6 |      25063   |       25063   |           0.667 |           0.6668 | True                   |

## K=4 Label Alignment

Mplus raw class order differs from the manuscript order, as expected. Alignment:

| manuscript_label   |   manuscript_class |   Mplus_raw_class |   Python_modal_n |   Mplus_modal_n |   difference_Mplus_minus_Python |
|:-------------------|-------------------:|------------------:|-----------------:|----------------:|--------------------------------:|
| C1 Stable-Low      |                  1 |                 2 |              615 |             616 |                               1 |
| C2 Persistent-High |                  2 |                 4 |              245 |             245 |                               0 |
| C3 Decreasing      |                  3 |                 1 |              726 |             725 |                              -1 |
| C4 Late-Increasing |                  4 |                 3 |              866 |             866 |                               0 |

Python vs Mplus modal class agreement after alignment:

- Percent agreement = 99.96%
- Cohen's kappa = 0.9994

The tiny mismatch reflects rounding/optimizer implementation differences, not a different class solution.

## K=4 Mplus Class Profiles

| manuscript_label   |   mplus_raw_class |     I |      S |      Q |   mean_w1 |   mean_w4 |   mean_w7 |
|:-------------------|------------------:|------:|-------:|-------:|----------:|----------:|----------:|
| C1 Stable-Low      |                 2 | 1.326 | -0.016 |  0.009 |     1.455 |     1.326 |     1.359 |
| C2 Persistent-High |                 4 | 2.614 | -0.043 | -0.039 |     2.392 |     2.614 |     2.134 |
| C3 Decreasing      |                 1 | 1.932 | -0.1   |  0.005 |     2.277 |     1.932 |     1.677 |
| C4 Late-Increasing |                 3 | 1.864 |  0.073 | -0.017 |     1.492 |     1.864 |     1.93  |

## Manuscript-ready wording

Independent re-estimation in Mplus 7 using the same 9-item quadratic LCGA specification returned fit indices indistinguishable from the Python EM and R lcmm solutions (Mplus K=4 BIC = 25,228.08; logL = -12,551.60). After label alignment, Mplus modal class counts differed from Python by at most 1 individuals across classes, supporting cross-platform reproducibility of the four-class solution.
