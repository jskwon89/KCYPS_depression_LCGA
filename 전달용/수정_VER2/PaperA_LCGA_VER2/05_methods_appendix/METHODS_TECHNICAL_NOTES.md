# Paper A — Methods Technical Notes

본문 §2.3 Methods에 대응. 작성 시 detail 부족하면 여기 참조.

---

## 1. LCGA fitting (Python EM)

**Spec**: Y_{it} = β_{0k} + β_{1k}(t-4) + β_{2k}(t-4)² + ε_{it}, ε_{it} ~ N(0, σ²)
- Time origin: wave 4 (mid-adolescence)
- Class-invariant residual variance σ²
- Within-class growth-parameter variance = 0 (Nagin GBTM)

**Algorithm**: EM with FIML for missing waves
- 8 random restarts (avoid local optima)
- Convergence tolerance 1e-5
- Max 300 iterations per restart

**Filter**: ≥3 valid waves of 9-item depression → N=2,452

**Class enumeration criteria**:
- BIC, aBIC, AIC (monotonic in K-1, K, K+1)
- Entropy ≥ 0.7 desirable
- Smallest class proportion ≥ 5% (Nagin 2005)
- BLRT (parametric bootstrap, K=2 vs K=1; B=30 due to compute)

**Files**:
- Script: `06_reproducibility/scripts/22_lcga_9item_primary.py`
- Output: `04_results_data/class_assignments_K4_9item.csv`, `class_profiles_K4_9item.csv`, `fit_summary_9item.csv`

## 2. R lcmm cross-validation

**Spec**: same polynomial-in-time mean trajectories; `random = ~ -1` (no random effects, equivalent to LCGA)
- 30 gridsearch starts
- Identical filter (≥3 valid waves)

**Cross-validation result**: K=4 BIC R lcmm = 25,228.08 vs Python EM = 25,228.07; modal-n agreement ±1 individual across all 4 classes after explicit label alignment.

**R class label mapping** (note for reviewers):
- R class 1 ↔ C1 Stable-Low
- R class 3 ↔ C2 Persistent-High
- R class 2 ↔ C3 Decreasing
- R class 4 ↔ C4 Late-Increasing

**Files**:
- Script: `06_reproducibility/scripts/26_lcmm_9item.R`, `31_lcmm_label_alignment.py`
- Output: `04_results_data/lcmm_R/` (fit_summary, posterior, class_pred, label_alignment)

## 3. BCH 3-step distal outcome correction

**Procedure** (Bolck-Croon-Hagenaars; Vermunt 2010):
1. Fit unconditional LCGA → posterior P(C=s | Y_i)
2. Build classification-error matrix D where D[s,t] = E_post[modal=t | true class=s]
3. Use D⁻¹ to weight individuals:
   P̂(Y_dist=1 | true class=s) = Σ_i (D⁻¹)[modal_i, s] · Y_dist_i / Σ_i (D⁻¹)[modal_i, s]

**Bootstrap**: 2000 nonparametric resamples, percentile 95% CI, **clipped to [0,1]**

**Why BCH**: When entropy is moderate (0.66 in our K=4 solution), modal-class logistic systematically underestimates class-distal associations. Diagonal of D (correct classification rate): 84.6%, 81.0%, 78.1%, 78.9%.

**Files**:
- Script: `06_reproducibility/scripts/23_bch_cox_9item.py`
- Output: `04_results_data/bch/` (proportions, ORs, classification matrix)

## 4. Cox time-to-event analysis

**Outcome**: time of first endorsement of strong SI (YPSY4E04 ≥ 3) across W1–W7, right-censored at last observed wave.

**Cox PH model**: class indicators (C1 reference) → average HR per class.

**PH assumption test**: Schoenfeld residuals via `lifelines.statistics.proportional_hazard_test`, `time_transform='rank'`.

### Results
- C2: p=.477 (PH OK)
- C3: **p<.001 (PH violated)**
- C4: **p<.001 (PH violated)**

**Time-varying coefficient model** (when PH violated):
- Long-form: each individual contributes (start, stop, event) intervals per wave
- Class-by-wave interaction terms: c2_x_t, c3_x_t, c4_x_t
- Per-wave HR = exp(β_main + w · β_interaction)

### Why cluster bootstrap, not lifelines `cluster_col=ID, robust=True`?

`lifelines.CoxPHFitter.variance_matrix_` is the **Hessian inverse (model-based covariance)**, NOT the sandwich/cluster-robust covariance, even when fitted with `robust=True`. The robust SEs are stored in `summary['se(coef)']` but the variance_matrix_ used for delta-method per-wave CI computation is naive. Therefore:

1. We use **ID-level cluster bootstrap** (B=300; resample unique IDs with replacement, refit Cox per resample)
2. Per-wave HR computed from each bootstrap fit
3. **95% percentile CI** + **two-sided percentile bootstrap p-value** (add-one rule, floor 1/301 ≈ .003)
4. Cross-checked against `CoxTimeVaryingFitter(id_col='ID')` for point estimates (identical)

**Files**:
- Schoenfeld test: `04_results_data/cox_PH_schoenfeld_test.csv`
- Naive long-form Cox (reference): `06_reproducibility/scripts/24_cox_time_varying_9item.py`
- Lifelines reported "robust" (Hessian inverse, kept as cautionary): `29_cox_time_varying_robust.py`
- ⭐ **Cluster bootstrap (PRIMARY)**: `30_cox_cluster_bootstrap.py` → `cox_HR_per_wave_with_CI_9item_BOOT.csv`

## 5. Bootstrap LRT (BLRT) for K-vs-K-1 class enumeration

**Procedure** (Nylund-Asparouhov-Muthen 2007; McLachlan 1987):
1. Fit K-1 model on observed data, record LL_(K-1)
2. Fit K model, record LL_K → observed LRT_obs = 2(LL_K - LL_{K-1})
3. For b=1..B: simulate Y* from K-1 fit, refit K-1 and K, record LRT_b
4. **Add-one rule**: p_BLRT = (1 + #{LRT_b ≥ LRT_obs}) / (B + 1)

**Floor for B=30**: 1/31 ≈ .032 (cannot report p<.032 with B=30 and zero exceedances).
**Floor for B=999**: 1/1000 = .001 exactly. **B≥1000** needed for p<.001.

### K=2 vs K=1 result
- observed LRT = 2,174.7
- B=30 null distribution: mean=3.74, SD=2.50, max=11.51
- All 30 bootstrap LRTs < observed (188× smaller than max)
- **p_BLRT = .032** (= 1/31)

K>2 vs K-1 BLRTs not run due to compute cost; relied on monotonic BIC drop, entropy ≥ 0.66, 5%-class-size rule.

**Files**:
- Script: `06_reproducibility/scripts/28_blrt_9item.py`
- Output: `04_results_data/blrt/blrt_summary_9item.csv`, `blrt_dist_K1vs2_9item.npy`

## 6. Longitudinal sampling weight sensitivity

**Weight**: KCYPS standardized longitudinal weight WEIGHTB2w7 (mean=1, applied at individual level)

**Procedure**: re-fit LCGA with weighted EM (each individual's likelihood multiplied by weight).

**Comparison**: Hungarian alignment of class labels (linear_sum_assignment), then Cohen's κ on aligned modal assignments.

### Result
- **Post-alignment percent agreement: 95.6%, Cohen κ = 0.939**
- Class structure robust to differential attrition

**Files**:
- Script: `25_weighted_sensitivity_9item.py`
- Output: `04_results_data/sensitivity_weighted/` (weighted profiles, alignment crosstab, kappa.txt)

## 7. Software and reproducibility

- **Python 3.11**: numpy, pandas, scipy, statsmodels, lifelines, scikit-learn
- **R 4.5.2**: lcmm 2.x, haven, dplyr, tidyr
- All random seeds fixed (`seed=42`) for EM restarts and bootstrap.
- Total runtime on a single CPU: LCGA ≈ 5 min; BLRT ≈ 6 min; cluster bootstrap ≈ 30 sec; R lcmm ≈ 30 min (gridsearch).
