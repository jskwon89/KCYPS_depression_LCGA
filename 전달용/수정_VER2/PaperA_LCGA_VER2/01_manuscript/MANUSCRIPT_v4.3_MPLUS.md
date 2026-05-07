# Heterogeneous Trajectories of Adolescent Depression and Risk of Suicidal Ideation in Emerging Adulthood: A Seven-Year Latent Class Growth Analysis of the Korean Children and Youth Panel Survey 2018

**Version:** v4.3 (2026-05-08, evening) — added full Mplus 7 replication of the 9-item LCGA K=1-5 series. Mplus reproduced Python EM/R lcmm fit indices to rounding precision (K=4 BIC = 25,228.08; logL = -12,551.60), with label-aligned Python vs. Mplus modal agreement = 99.96% and Cohen's κ = .999. Earlier v4.2 precision fixes preserved: B=999 yields exactly p=.001, so p<.001 requires B≥1000; cluster bootstrap p=1/301=.0033 is reported as p=.003 (rounded), never p≤.003.

**Earlier versions:**
- v4.2 (2026-05-08): strict add-one p-value floors corrected for BLRT and Cox bootstrap.
- v4 (2026-05-08): Table 2 BLRT p corrected to .032; Methods B=30 with add-one rule; per-wave Cox CIs replaced with ID-level cluster bootstrap (B=300); class×wave interaction terms p=.003 (c3, c4); R lcmm <-> manuscript label mapping table.
- v3 (2026-05-08): updated 9-item BLRT, R lcmm 9-item completed, modal-n vs π̂ separated, Hessian-based "robust" Cox (later corrected in v4).
- v2 (2026-05-08): switched to 9-item depression composite (criterion contamination fix), Hungarian-aligned weighted κ, Schoenfeld test added.
- v1 (2026-05-07): 10-item composite (now sensitivity); modal-class logistic distal outcomes.

**Target journals:**
1. *Journal of Affective Disorders* (IF≈6.6, Q1)
2. *Development and Psychopathology* (IF≈3.7, Q1)
3. *Journal of Youth and Adolescence* (IF≈3.7, Q1)
4. *Suicide and Life-Threatening Behavior* (IF≈3.5, Q1)

**Authors:** [TBD]
**Word count target:** ~5,500

---

## ABSTRACT (~300 words)

**Background.** Korean adolescents have one of the highest age-standardized suicide mortality rates among OECD countries, yet developmental heterogeneity in depressive trajectories and their distal connection to suicidal ideation in emerging adulthood remains underexamined.

**Methods.** Using the Korean Children and Youth Panel Survey 2018 m1 cohort (N=2,452, age 13 at Wave 1 in 2018; 80.0% retained at Wave 7 in 2024, age 19), we fit latent class growth models on a **9-item depression composite** that excludes the suicidal-ideation item (E04 "I have thoughts of dying") to avoid criterion contamination with the distal outcomes; the suicidal-ideation item was retained as a separate distal indicator. Class enumeration was guided by Bayesian Information Criterion, parametric bootstrap likelihood-ratio test (B = 30), entropy, and minimum-class-proportion considerations, with independent R lcmm and Mplus cross-validation. Distal outcomes at W7 were estimated using the Bolck-Croon-Hagenaars (BCH) three-step method with bootstrap confidence intervals clipped to the [0,1] range. Time to first strong suicidal ideation across W1–W7 was modeled with Cox proportional hazards regression; the proportional-hazards assumption was tested via Schoenfeld residuals and, where violated, supplemented with a class × wave interaction model whose per-wave hazard-ratio confidence intervals were computed by ID-level cluster bootstrap (B = 300).

**Results.** A four-class solution provided the best balance of fit and interpretability (BIC = 25,228; entropy = 0.66; smallest class 10.5%; bootstrap LRT for K=2 vs K=1 with B=30: observed LRT = 2,174.7, null distribution mean = 3.7 [max = 11.5], all bootstrap LRTs below observed, p_BLRT = 0.032 — the lowest value attainable with B=30 and zero exceedances): C1 *Stable-Low/Resilient* (π̂ = 25.4%, modal n = 615), C2 *Persistent-High/Inverted-U* (π̂ = 10.5%, n = 245), C3 *Decreasing* (π̂ = 29.3%, n = 726), and C4 *Late-Increasing* (π̂ = 34.8%, n = 866). Independent re-estimation in R (lcmm) and Mplus 7 returned the same K=4 fit (BIC = 25,228.08) and near-identical class assignments after label alignment (Mplus vs Python modal agreement = 99.96%, κ = .999). At W7, BCH-corrected odds ratios for any suicidal ideation versus C1 were **15.2** (95% CI 9.6–29.0) for C4 *Late-Increasing*, 10.6 (6.5–20.2) for C2, and 5.1 (3.2–9.4) for C3. Schoenfeld residual tests indicated proportional-hazards violation for C3 and C4 (p < .001); under a time-varying class × wave Cox model with ID-level cluster bootstrap (B = 300), per-wave hazard ratios (vs. C1) were 1.43 [0.94, 2.30] at W1 (bootstrap p = .10, n.s.) and 7.45 [5.00, 12.97] at W7 (bootstrap p = .003) for C4, and 8.96 [6.29, 14.89] at W1 to 1.35 [0.75, 2.55] at W7 for C3, with crossover at W4–5 (the high-school transition); the class × wave interaction coefficients themselves are significant under cluster bootstrap (c3_x_t bootstrap p = .003; c4_x_t bootstrap p = .003). All bootstrap p-values are reported at three decimals and are floored at the add-one limit 1/(B+1) = 1/301 = .0033 (rounded to .003); a larger B would only refine, not contradict, this floor. Females, low parental warmth, high parental rejection/chaos, low self-esteem, and high smartphone dependence at W1 each significantly differentiated C1 from clinical-risk classes.

**Conclusions.** Approximately one in three Korean adolescents follow a *late-increasing* depression trajectory that is largely indistinguishable from the resilient class at age 13 but reaches the highest emerging-adult suicidal-ideation risk by age 19. Single-time-point screening at school entry misses this large at-risk subgroup; serial assessment across the high-school transition is required.

**Keywords:** adolescent depression; latent class growth analysis; suicidal ideation; longitudinal trajectory; KCYPS; emerging adulthood

---

## 1. INTRODUCTION

### 1.1 Korean adolescent mental health context
- Suicide is the leading cause of death among Korean adolescents and emerging adults aged 9–24 [Statistics Korea, 2023].
- Lifetime suicidal ideation prevalence in Korean adolescents has been documented in the 25–35% range across cross-sectional studies [Park et al., 2018; Kim & Kim, 2022].
- Existing longitudinal trajectory studies of Korean youth depression have largely used 4-wave KCYPS-1 data covering middle-school years only [Lee et al., 2017; Choi & Park, 2019]. The 2018 panel offers a seven-wave window spanning the developmentally critical transition from early adolescence to emerging adulthood.

### 1.2 Heterogeneity in depression trajectories
- Theoretical framing: developmental psychopathology and equifinality / multifinality [Cicchetti & Rogosch, 1996].
- Trajectory typologies in Western adolescents typically show 3–5 classes including stable-low, decreasing, increasing, and stable-high [Stoolmiller et al., 2005; Costello et al., 2008; Yaroslavsky et al., 2013]. Whether these patterns replicate and what their *distal* consequences for late-adolescent suicidal ideation are in Korea remains undetermined for a 7-wave window.

### 1.3 Distal connection to suicidal ideation and the criterion-validity issue
- Single-time-point depression scores correlate with current suicidal ideation, but trajectory shape may capture *cumulative* risk better than mean level.
- **Important methodological point.** Most prior trajectory studies of depression have used composite depression scales that include suicidal-ideation items, then reported associations with later suicidal-ideation outcomes. This introduces criterion contamination: the trajectory is partly defined by the very outcome we are predicting. We address this directly by constructing a 9-item depression composite that *excludes* the suicidal-ideation item (YPSY4E04 "I have thoughts of dying") and by treating suicidal ideation only as a distal outcome. Sensitivity analysis with the original 10-item composite yields a structurally similar solution and is reported in the supplement.

### 1.4 Predictors of trajectory class membership
- Skinner et al. (2005) PSCQ-A six-dimensional parenting (warmth, rejection, autonomy support, coercion, structure, chaos), peer attachment, teacher relationship, self-esteem, life satisfaction, grit, smartphone dependence — full battery available in KCYPS at W1.
- Sex differences are well documented in adolescent depression onset and in suicidal ideation; we test for class composition by sex.

### 1.5 Aims
1. Characterize the number and shape of latent depression trajectories from age 13 to 19 in a Korean nationally representative cohort, using a criterion-clean 9-item depression measure.
2. Estimate associations of W1 family, peer, school, and individual characteristics with trajectory class.
3. Evaluate whether trajectory class predicts emerging-adult (W7) suicidal ideation, mean depression, and life satisfaction with proper correction for class-assignment uncertainty (BCH).
4. Examine time-to-first-strong-suicidal-ideation across all seven waves; quantify the time-varying nature of class hazards.

---

## 2. METHOD

### 2.1 Sample and design
- **Source.** Korean Children and Youth Panel Survey 2018 (KCYPS 2018), middle-school first-grade (m1) cohort.
- **Sampling.** Two-stage stratified cluster sampling with region (17 provinces) and city size as strata, conducted by the National Youth Policy Institute. The original m1 cohort comprised N=2,590 students entering middle school grade 1 in 2018.
- **Waves.** Annual surveys 2018 (W1, age 13) → 2024 (W7, age 19). Sample retention: 94.1% (W2), 87.5% (W4), 80.0% (W7).
- **Analytic sample.** N=2,452 individuals with at least three valid waves of depression measurement (94.7% of the original cohort).
- **Ethics.** Secondary analysis of de-identified public-use data; IRB approval [TBD].

### 2.2 Measures

**Depression — primary trajectory measure (9 items, criterion-clean).** Mean of YPSY4E01, E02, E03, E05, E06, E07, E08, E09, E10 (excluding E04 *"I have thoughts of dying"*). Items rated 1 (*not at all*) to 4 (*very much*). Cronbach α = .91, .91, .91, .90, .89, .89, .90 across W1–W7. Correlation with the original 10-item composite was r = .996, indicating that exclusion of E04 leaves the construct essentially intact for trajectory characterization.

**Depression — sensitivity (10-item).** Mean of all ten YPSY4E items (α = .92, .92, .92, .90, .90, .90, .91); used in supplementary sensitivity analysis only.

**Suicidal ideation — distal outcome.** Item YPSY4E04. Two binarizations: *any* (≥2, "sometimes or more") and *strong* (≥3, "often or more").

**Baseline (W1) covariates.**
- *Parenting* (Skinner PSCQ-A, six 4-item subscales, α = .79–.91): warmth, rejection, autonomy, coercion, structure, chaos.
- *Peer attachment* (13 items, α = .85), *Teacher relationship* (14 items, α = .91).
- *Self-esteem* (Rosenberg 10-item, α = .87), *Life satisfaction* (Diener SWLS 5-item, α = .85), *Subjective happiness* (Lyubomirsky 4-item, α = .80), *Grit* (Duckworth Short Grit 8-item, α = .71).
- *Smartphone dependence* (15-item, α = .88).
- *Demographics:* sex, region, urban size.

### 2.3 Analytic strategy

**Latent class growth analysis.** Finite mixtures of polynomial-in-time mean trajectories with class-invariant residual variance and zero within-class growth-parameter variance (Nagin, 1999, 2005):
  Y_{it} = β_{0k} + β_{1k}(t − 4) + β_{2k}(t − 4)² + ε_{it},   ε_{it} ~ N(0, σ²).
EM with full-information maximum likelihood for missing waves; eight random restarts. Solutions for K=1–5 compared on:
- BIC, sample-size-adjusted BIC (aBIC), AIC
- Entropy (≥0.7 indicates good class separation)
- Smallest-class proportion (≥5%, Nagin's recommendation)
- **Bootstrap likelihood-ratio test (BLRT, Nylund et al., 2007)** for K = 2 vs K = 1: B = 30 parametric replicates from the K=1 fit. The reported BLRT p-value follows the add-one rule p_BLRT = (1 + #{LRT_b ≥ LRT_obs}) / (B + 1), so the smallest p attainable with zero exceedances at B = 30 is 1/31 ≈ .032 (with B = 100 it would be 1/101 ≈ .010, with B = 500 it would be 1/501 ≈ .002, with B = 999 exactly 1/1000 = .001, and B ≥ 1000 is needed for p < .001). Bootstrap of larger K-vs-K-1 comparisons was not run due to computational cost; class enumeration beyond K=2 relied on BIC drop, entropy, and class-size considerations, as in Yaroslavsky et al. (2013).

**R lcmm replication.** As an independent methodological check, we re-fit K=1–5 LCGA models in R (lcmm 2.x; Proust-Lima et al., 2017) with 30 random starts via *gridsearch*. This provides convergence validation across two independent optimizers (Python EM vs R Newton).

**Class profiling and predictors of class membership.** W1 covariates compared across classes via one-way ANOVA with η² as effect size.

**Distal outcomes — BCH three-step correction.** Posterior probabilities P(C=s | Y_i) from the LCGA were used to construct a classification-error matrix D where D[s,t] = P(modal=t | true class=s). The BCH-corrected proportion / mean for outcome Y_dist is
  Ê[Y_dist | true class = s] = Σ_i (D⁻¹)[modal_i, s] · Y_dist_i / Σ_i (D⁻¹)[modal_i, s].
Standard errors and 95% confidence intervals were obtained by 2,000 nonparametric bootstrap resamples of individuals; bootstrap proportions were clipped to [0,1]. This BCH approach (Bolck, Croon & Hagenaars, 2004; Vermunt, 2010; Bakk & Vermunt, 2016) accounts for class-assignment uncertainty when entropy is moderate.

**Time-to-event analysis.** Time of first endorsement of strong suicidal ideation (E04 ≥ 3) across W1–W7 was treated as the event. Cox proportional hazards regression with class indicators (C1 *Stable-Low* reference) yielded average hazard ratios. The proportional-hazards assumption was tested via Schoenfeld residuals (lifelines `proportional_hazard_test`, time_transform='rank'). Where the assumption was violated, we fitted a time-varying coefficient model with class × wave interaction. Because each individual contributes multiple wave-intervals to the long-form likelihood and the model-based covariance from `lifelines.CoxPHFitter.variance_matrix_` is the Hessian inverse rather than the cluster-sandwich covariance, we obtained per-wave hazard-ratio confidence intervals via **ID-level cluster bootstrap** (B = 300 resamples of unique individuals with replacement; per-wave HR computed from the bootstrap parameter draws; 95% percentile CIs and two-sided percentile bootstrap p-values reported). This yields cluster-correct uncertainty without relying on lifelines' Hessian-based variance matrix.

**Sensitivity analyses.**
- *10-item depression composite.* Re-estimate the LCGA with the original 10-item depression as the trajectory measure; compare class structure.
- *Longitudinal sampling weights.* Re-estimate with the W7 standardized longitudinal weight (WEIGHTB2w7); compute Cohen's κ between unweighted and weighted class assignment after Hungarian label alignment.
- *Modal-class versus BCH.* Compare uncorrected modal-class logistic ORs for distal outcomes against BCH-corrected estimates.

**Software.** R 4.5.2 (lcmm 2.x), Python 3.11 (numpy, scipy, statsmodels, scikit-learn, lifelines).

---

## 3. RESULTS

### 3.1 Sample characteristics
[See Table 1.] Analytic sample N=2,452 (94.7% of original m1 cohort, ≥3 valid waves of depression). Sex distribution: 54.3% male, 45.7% female.

### 3.2 Number of classes
[See Table 2 and Figure 5.] Information criteria for the **9-item primary** LCGA decreased monotonically from K=1 (BIC = 28,317.6) through K=5 (BIC = 25,126.5). The parametric bootstrap LRT for K=2 vs. K=1 was significant (observed LRT = **2,174.7**; B=30 bootstrap samples drawn from the K=1 fit had a null LRT mean of 3.7 (SD = 2.5, max = 11.5); all 30 bootstrap LRTs were far below the observed statistic, yielding p_BLRT = 0.032 — the lowest p attainable with zero exceedances at B=30 under the add-one rule). The observed LRT was 188× larger than the largest bootstrap value, so the effect is qualitatively unambiguous; reaching the conventional p < .001 threshold under the same add-one rule would require **B ≥ 1000** (B=999 yields exactly 1/1000 = .001; B=500 floors at ≈ .002, B=100 at ≈ .010), and we plan to extend the bootstrap to that range in the supplement. Entropy was highest at K=3 (0.719) and remained adequate at K=4 (0.661). At K=5 the smallest class fell to 3.7% of the sample, below the 5% minimum recommended by Nagin (2005). We retained the four-class solution as the primary model on grounds of parsimony, theoretical interpretability, and adequate class size; bootstrap LRT for K>2 vs. K=K-1 comparisons was not feasible within available compute and is reported in supplementary material as planned future analysis. The K-decision beyond K=2 therefore rests on monotonic BIC reduction, entropy ≥ 0.66, and the 5%-minimum rule rather than on a formal BLRT for each K.

**R lcmm and Mplus cross-validation.** Independent re-fitting of K=1–5 in R (lcmm 2.x, 30 gridsearch starts) returned essentially identical fit indices (R K=4 BIC = 25,228.08 vs. Python EM K=4 BIC = 25,228.07; R logL = −12,551.60 = Python logL to two decimal places) and modal class assignments. A second independent replication in Mplus 7 using the same quadratic LCGA constraints (growth-factor variances fixed to zero; one class-invariant residual variance shared across waves) reproduced the entire K=1–5 fit series to rounding precision (Mplus K=4 BIC = 25,228.08; logL = −12,551.60). Because Python EM, R lcmm, and Mplus assign internal class labels in arbitrary order, we report the alignments below. Side-by-side modal counts after alignment are within one individual across all four classes, demonstrating numerically equivalent solutions:

| Manuscript label | Python EM raw class | Python n (%) | R lcmm raw class | R n (%) |
|------------------|---------------------|---------------|-------------------|----------|
| C1 Stable-Low | 1 | 615 (25.08%) | 1 | 616 (25.12%) |
| C2 Persistent-High | 2 | 245 (9.99%) | 3 | 245 (9.99%) |
| C3 Decreasing | 3 | 726 (29.61%) | 2 | 725 (29.57%) |
| C4 Late-Increasing | 4 | 866 (35.32%) | 4 | 866 (35.32%) |

Mplus raw labels aligned as follows: Mplus class 2→C1, class 4→C2, class 1→C3, class 3→C4. The aligned Mplus modal counts were C1=616, C2=245, C3=725, C4=866; Python vs. Mplus modal agreement was 99.96% (Cohen's κ=.999). The aligned summary files are `lcmm_R/python_vs_R_class_comparison_9item.csv` and `mplus/python_vs_mplus_class_comparison_9item.csv`; the full Mplus replication report, including log-likelihood, AIC, BIC, aBIC, entropy, and free-parameter counts for K=1–5, is `mplus/MPLUS_REPLICATION_REPORT.md`. Convergence across Python EM, R lcmm, and Mplus strengthens confidence in the four-class solution.

### 3.3 Class profiles
[See Table 3 and Figure 1.] Below we report the EM-estimated mixing proportion π̂_k (the model-implied class probability) and, separately, the modal-class assignment count n_k (the result of argmax over posterior probabilities for each individual). The two differ slightly because entropy is moderate.
- **C1 Stable-Low / Resilient** (π̂ = 25.4%; modal n = 615): intercept 1.33, near-zero slope; mean depression 1.32–1.46 across waves.
- **C2 Persistent-High / Inverted-U** (π̂ = 10.5%; modal n = 245): intercept 2.61, slight inverted-U with peak at W3 (2.62); depression 2.13–2.62.
- **C3 Decreasing** (π̂ = 29.3%; modal n = 726): intercept 1.93, linear slope −0.10; depression 2.28 (W1) → 1.68 (W7).
- **C4 Late-Increasing** (π̂ = 34.8%; modal n = 866): intercept 1.86, linear slope +0.07; depression 1.49 (W1) → 1.93 (W7).

Throughout the remainder of this paper, "class size" refers to modal-assignment n; mixing proportions are used for class-conditional distal estimation under the BCH framework, where they enter via the classification-error matrix.

The four trajectory shapes are essentially the same as those obtained with the 10-item sensitivity composite (Hungarian-aligned modal-assignment agreement = 94.9%, Cohen's κ = 0.929; see supplement).

### 3.4 Baseline (W1) predictors of class membership
[See Table 4 and Figure 3.] All eighteen W1 psychosocial constructs differed significantly across the four trajectory classes (all *F*s > 38, *p* < .001). The largest effects involved depression itself (η² ~ .47), self-esteem (η² ~ .25), aggression (η² ~ .21), and somatic symptoms (η² ~ .24). Skinner six-dimensional parenting differentiated classes for all six dimensions: C1 *Stable-Low* showed the most warmth (M = 3.52), autonomy (3.49), and structure (3.14), and the least rejection (1.60), coercion (2.05), and chaos (1.87). C2 *Persistent-High* showed the inverse pattern. Smartphone dependence at W1 was lowest in C1 (1.90) and highest in C2 (2.25). Sex distribution differed sharply across classes: female representation was 40% in C1 (Stable-Low) and C4 (Late-Increasing), versus 53% in C3 (Decreasing) and 61% in C2 (Persistent-High); χ² test highly significant (p < 10⁻¹¹).

### 3.5 Distal outcomes at W7
[See Table 5 and Figures 2, 4.] BCH-corrected proportions endorsing **any suicidal ideation** at age 19 (BCH OR vs. C1 Stable-Low):
- C1: 8.0% (95% CI 4.5–11.2%) — reference
- C3 *Decreasing*: 30.4% (25.5–35.2%); OR 5.05 (3.21–9.45)
- C2 *Persistent-High*: 47.7% (40.5–55.6%); OR 10.58 (6.49–20.21)
- **C4 *Late-Increasing*: 56.8% (52.3–61.8%); OR 15.23 (9.57–29.00)**

For **strong suicidal ideation** (≥often) at W7:
- C1: 1.1%; C3: 2.4%; C2: 13.6%; C4: 11.5%
- BCH ORs vs. C1: C3 = 2.26; C2 = 14.48; **C4 = 11.97**.

**Mean depression at W7 (BCH-corrected):** C1 = 1.17, C3 = 1.53, C4 = 2.04, C2 = 2.19.
**Life satisfaction at W7 (BCH-corrected):** C1 = 2.73, C3 = 2.53, C4 = 2.49, C2 = 2.32.

### 3.6 Time to first strong suicidal ideation (W1–W7) — including PH violation diagnostic
[See Table 6 and Figures 2C, 6.] Across the seven-year observation window, 962 of 2,452 individuals (39.2%) endorsed strong suicidal ideation at some wave. Class-specific event rates were 13.5% (C1), 75.5% (C2), 49.3% (C3), and 38.8% (C4).

**Cox proportional hazards regression** (C1 reference) yielded average hazard ratios:
- C2 *Persistent-High*: HR = 10.20 (95% CI 7.87–13.23), *p* < .001
- C3 *Decreasing*: HR = 5.15 (4.06–6.54), *p* < .001
- C4 *Late-Increasing*: HR = 3.15 (2.47–4.00), *p* < .001

**Schoenfeld residual test** indicated that the proportional-hazards assumption was tenable for C2 (*p* = .48) but violated for C3 (*p* < .001) and C4 (*p* < .001). This violation is *substantively meaningful*: it reflects that the C3 (Decreasing) hazard is high early and declines, whereas the C4 (Late-Increasing) hazard is low early and rises across adolescence.

**Time-varying class × wave interaction Cox model** quantifies this dynamic [see Figure 6 and Table 6 for the W1/W4/W7 point estimates with cluster-bootstrap 95% CIs and bootstrap interaction p-values; per-wave estimates for all seven waves are in supplementary Table 6b]. Per-wave hazard-ratio confidence intervals were obtained via ID-level cluster bootstrap (B = 300 resamples of unique individuals with replacement), since lifelines' Hessian-based `variance_matrix_` does not correspond to the cluster-sandwich covariance even when the model is fitted with `cluster_col='ID', robust=True`. Per-wave point estimates with cluster-bootstrap 95% percentile CIs are:
- C2 *Persistent-High*: HR 12.31 [8.17, 20.67] at W1, declining to 6.79 [4.00, 14.29] at W7 — chronic risk, bootstrap p = .003 at every wave (the add-one floor for B = 300; see footnote).
- C3 *Decreasing*: HR 8.96 [6.29, 14.89] at W1, declining to 1.35 [0.75, 2.55] at W7 — early-only risk; bootstrap p = .003 through W5, p = .007 at W6, n.s. at W7 (p = .31).
- C4 *Late-Increasing*: HR 1.43 [0.94, 2.30] at W1 (bootstrap p = .10, n.s.), rising to 7.45 [5.00, 12.97] at W7 (bootstrap p = .003) — emerging risk; the C4 hazard becomes statistically distinguishable from C1 from W2 onward and exceeds the C3 hazard by W4–5 (the high-school transition).

The class × time interaction coefficients are themselves significant under the same cluster bootstrap: c3_x_t β = −0.315 [95% percentile CI −0.476, −0.183], bootstrap p = .003 (the C3 hazard genuinely declines across waves); c4_x_t β = +0.275 [+0.149, +0.413], bootstrap p = .003 (the C4 hazard genuinely rises); c2_x_t β = −0.099 [−0.247, +0.061], bootstrap p = .22 (no significant change for the chronic class).

*Footnote on bootstrap p-values.* All p-values from the ID-level cluster bootstrap are floored by the add-one rule at 1/(B+1) = 1/301 = .00332, which rounds to .003 at three decimals. We therefore report p = .003 at this floor rather than p < .003 or p ≤ .003, since the discrete bootstrap cannot resolve below .00332 at B = 300; a planned B = 1000 supplementary bootstrap would lower the floor to .000999 ≈ .001.

The substantive core of the "invisible at-risk" finding is therefore supported both descriptively and formally: the class with the highest *late*-adolescent hazard (C4) starts at a hazard that does not significantly differ from the resilient C1 class at W1, the rate at which the C4 hazard rises across waves is significant under cluster bootstrap (bootstrap p = .003), and by W7 the C4 hazard ratio reaches 7.45 (bootstrap p = .003).

### 3.7 Sensitivity analyses
- **10-item depression composite (with E04).** The same four-class structure was recovered; modal-class agreement after Hungarian label alignment was 94.9% (Cohen's κ = 0.929). This indicates that class structure itself is not driven by the suicidal-ideation item. However, distal-outcome ORs against suicidal ideation are inflated under the 10-item composite (e.g., C4 *Late-Increasing* BCH OR rises from 15.2 [9-item primary] to 17.1 [10-item sensitivity]) — a direct demonstration of criterion contamination, supporting the choice of the 9-item composite as primary [supplement].
- **Longitudinal weights (WEIGHTB2w7).** Re-fitting with the W7 standardized longitudinal weight produced a four-class solution with the same trajectory shapes; Hungarian-aligned percent agreement was 95.6% and Cohen's κ = 0.939. Class structure is robust to differential attrition.
- **Modal-class versus BCH.** Modal-class logistic ORs for any suicidal ideation are systematically smaller than BCH-corrected ORs (e.g., C4 modal OR ≈ 6 vs. BCH OR = 15.2), illustrating the bias-correction value of BCH when entropy is moderate.
- **R lcmm and Mplus cross-validation.** Independent re-estimation in R lcmm 2.x with 30 gridsearch starts returned K=4 BIC = 25,228.08 vs. Python EM K=4 BIC = 25,228.07 and modal class counts within one individual across all four classes after explicit label alignment (R class 1↔C1, R class 3↔C2, R class 2↔C3, R class 4↔C4). Mplus 7 replication of the same LCGA constraints reproduced K=1–5 log-likelihoods/BIC values to rounding precision; after label alignment (Mplus class 2↔C1, 4↔C2, 1↔C3, 3↔C4), Python vs. Mplus modal agreement was 99.96% (κ=.999) [`lcmm_R/python_vs_R_class_comparison_9item.csv`, `mplus/MPLUS_REPLICATION_REPORT.md`].
- **Cox time-varying coefficients with ID-level cluster bootstrap (B = 300).** The cluster bootstrap is the appropriate uncertainty quantification for a long-form Cox model in which each individual contributes multiple wave-intervals; lifelines' `variance_matrix_` is the Hessian inverse rather than the cluster-sandwich covariance, even when the model is fitted with `robust=True`, so we did not rely on delta-method CIs. Bootstrap percentile 95% CIs for the per-wave HR retain the same point pattern as the original fit (C4 W1 HR 1.43 [0.94, 2.30] n.s. → C4 W7 HR 7.45 [5.00, 12.97], bootstrap p = .003) and the class × wave interaction coefficients are significant (c3_x_t p = .003; c4_x_t p = .003) [supplement]. All bootstrap p-values are floored at 1/(B+1) = 1/301 = .0033 ≈ .003.

---

## 4. DISCUSSION

### 4.1 Principal findings
A criterion-clean four-class developmental typology of Korean adolescent depression, observed annually from age 13 to 19, distinguishes resilient (25%), recovering (29%), late-rising (35%), and chronically high (10%) trajectories. The largest at-risk subgroup for emerging-adult suicidal ideation is the *late-increasing* group, whose risk emerges only across the high-school transition.

### 4.2 The "invisible at-risk" finding (revised under the 9-item primary)
At Wave 1 (age 13), C4 (Late-Increasing) is essentially at the cohort grand mean on every measured covariate. By Wave 7, C4 reaches the highest BCH-corrected prevalence of any suicidal ideation across the four classes (56.8%, OR vs C1 = 15.2). Because C4 represents 35% of the cohort, it contributes roughly half of the absolute count of suicidal ideators at age 19.

The time-varying Cox model quantifies this dynamic: the hazard ratio for C4 vs. C1 rises from 1.4 at W1 to 7.5 at W7 — and crosses the C3 (Decreasing) hazard at W4–W5, exactly the developmental window of the high-school transition. Single-time-point screening at school entry would not distinguish C4 from C1 on any measured covariate.

### 4.3 Comparison to Western literature
- Class proportions in Western cohorts (Stoolmiller et al., 2005; Yaroslavsky et al., 2013): typically larger Stable-Low (~50%) and smaller Increasing (~10–15%). The Korean *Late-Increasing* class size of ~35% is unusually large; possible drivers include high-school examination culture, rapid digital media integration, and academic-track stratification.

### 4.4 Parenting and other predictors
- All six Skinner parenting dimensions differentiated trajectories. The largest single predictor was W1 depression itself, which was modestly elevated even in C4 vs C1 (z = +.14 vs −.46), suggesting that small W1 differences do propagate into trajectory differentiation but are insufficient on their own for screening.

### 4.5 Methodological strengths
- 9-item depression composite excludes E04 to avoid criterion contamination of the trajectory with the suicidal-ideation distal outcome — a contamination present in much prior trajectory-suicide work.
- BCH 3-step correction for class-assignment uncertainty.
- Schoenfeld test of PH assumption with time-varying Cox alternative when violated.
- Cross-platform replication across Python EM, R lcmm, and Mplus 7.

### 4.6 Limitations
- LCGA assumes class-invariant residual variance; growth mixture models with class-varying random effects could refine.
- Self-report measure of depression; no clinical interview.
- BLRT was conducted for K=1 vs. K=2 only (B=30); larger comparisons (K vs. K−1 for K=3–5) relied on monotonic BIC reduction, entropy ≥ 0.66, and the 5%-minimum-class-size rule rather than a formal bootstrap. Parametric BLRT for K=3–5 vs. K−1 is planned in supplementary follow-up.
- Cox PH violations were addressed via per-wave HR estimated under ID-level cluster bootstrap (B = 300). Both the per-wave HR pattern and the formal class × time interaction terms are statistically significant under this proper cluster-correct uncertainty quantification (c3_x_t and c4_x_t bootstrap p = .003), but the c2_x_t interaction is not (p = .22), suggesting that the chronic-high (C2) hazard remains constant across waves. A larger bootstrap (B ≥ 1000) and a comparison of cluster-bootstrap versus model-based sandwich estimators is planned for the supplement.
- BCH bootstrap CIs on rare-outcome (strong SI) ORs are wide due to small reference-cell counts.

### 4.7 Practical implications
- School-based depression screening should be longitudinal, not cross-sectional.
- The transition from middle to high school is an intervention window for the C4 *Late-Increasing* class (≈35% of cohort).
- Skinner six-dimension parenting at age 13 is a meaningful but not sufficient screening adjunct; depression mean at any single wave is necessary but not sufficient.

### 4.8 Future directions
- Multivariate parallel-process growth (depression × anxiety × aggression × smartphone dependence).
- Causal mediation: parenting → trajectory class → suicidal ideation.
- External validation in KCYPS-2018 e1 cohort and KYPS-1 archive.

### 4.9 Conclusion
A four-class developmental typology of Korean adolescent depression, estimated with a criterion-clean 9-item depression composite, reveals that the largest emerging-adult suicidal-ideation risk pool consists of those who are essentially at baseline at age 13 but rise into clinical range across the high-school transition. Time-varying hazard analysis shows that this class's risk crosses the formerly-elevated *Decreasing* class around Wave 4–5. The seven-year prospective design, BCH correction, Python/R/Mplus class-solution replication, and time-varying Cox modeling converge on a single message: depression trajectories matter, screening should be longitudinal, and the high-school transition is an under-recognized intervention window in Korean adolescent mental-health surveillance.

---

## REFERENCES (placeholder list — to be completed in JAD style)

[Same as v1; to add Proust-Lima et al. 2017 lcmm reference, lifelines.proportional_hazard_test reference]

---

## SUPPLEMENTARY MATERIAL

- **S1.** 10-item depression LCGA (sensitivity); Hungarian-aligned crosstab 9-item vs 10-item, κ=0.929.
- **S2.** Modal-class versus BCH-corrected distal outcome comparison.
- **S3.** Sex-stratified distal outcomes by class.
- **S4.** Sensitivity LCGA with longitudinal sampling weights (κ=0.939).
- **S5.** Schoenfeld residual test full output.
- **S6.** Cox per-wave HR table (time-varying).
- **S7.** Item-level descriptive statistics for YPSY4E depression scale by wave.
- **S8.** R and Python analysis scripts (archived at [URL]).

---

## REPRODUCIBILITY

All scripts, intermediate datasets, tables, and figures at:
`D:/2026/SCI/청소년패널조사/output/paperA_lcga/`

```
primary_9item/                      ⭐ PRIMARY analysis (criterion-clean)
├─ class_assignments_K4_9item.csv
├─ class_profiles_K4_9item.csv
├─ fit_summary_9item.csv
├─ tables/Table1..Table6_*9item.csv  (Table 2 BLRT p=.032 with B=30 floor)
├─ figures/Figure1..Figure6_*9item.{png,pdf}
│   ├─ Figure6_..._9item.png            (naive SE, kept for reference)
│   ├─ Figure6_..._9item_ROBUST.png     (Hessian-based, NOT actually robust — kept as cautionary)
│   └─ Figure6_..._9item_BOOT.png       ⭐ ID-level cluster bootstrap percentile CIs (PRIMARY for §3.6)
├─ bch/                              BCH 3-step + bootstrap CIs (clipped to [0,1])
├─ long_distal/                      Cox + Schoenfeld + time-varying analyses
│   ├─ cox_HR_class_9item.csv
│   ├─ cox_PH_schoenfeld_test.csv
│   ├─ cox_time_varying_9item.csv                          (naive long-form likelihood)
│   ├─ cox_time_varying_9item_ROBUST.csv                   (lifelines reported robust SE)
│   ├─ cox_HR_per_wave_with_CI_9item_BOOT.csv              ⭐ per-wave HR + bootstrap CI + boot p
│   ├─ cox_HR_per_wave_9item_BOOT.csv                      ⭐ point HR by wave wide
│   └─ cox_interaction_bootstrap_dist_9item.csv            ⭐ bootstrap interaction-coef distribution
├─ sensitivity_weighted/             Hungarian-aligned κ=0.939
├─ lcmm_R/                           R lcmm 9-item replication (K=4 BIC 25,228.08)
│   ├─ class_label_alignment_9item.csv         ⭐ R raw class -> manuscript label
│   ├─ python_vs_R_class_comparison_9item.csv  ⭐ side-by-side modal n
│   └─ class_pred_K4_9item_LABELED.csv         ⭐ relabeled trajectory means
├─ mplus/                            Mplus 7 replication (K=1-5; K=4 BIC 25,228.08)
│   ├─ MPLUS_REPLICATION_REPORT.md             ⭐ Python/R/Mplus comparison
│   ├─ mplus_vs_python_fit_summary_9item.csv   ⭐ K=1-5 fit side-by-side
│   ├─ python_vs_mplus_class_comparison_9item.csv
│   └─ lcga9_k1.out ... lcga9_k5.out           ⭐ Mplus output files
└─ blrt/                             BLRT K=1 vs K=2 (B=30, 9-item, p=.032)

sensitivity_10item/                  10-item composite (criterion-contaminated)
└─ ... (older outputs from v1, retained for sensitivity comparison)

(legacy files at root level — superseded by primary_9item/)
```

Source code:
- `06_build_long_panel.py` — long-form panel + scale composites + α (incl. depression9)
- `22_lcga_9item_primary.py` — PRIMARY 9-item LCGA EM
- `23_bch_cox_9item.py` — BCH + Cox + Schoenfeld
- `24_cox_time_varying_9item.py` — time-varying class × wave Cox (naive SE; for reference)
- `25_weighted_sensitivity_9item.py` — weighted sensitivity + Hungarian κ
- `26_lcmm_9item.R` — R lcmm replication (BIC 25,228.08 ≡ Python EM)
- `27_publication_outputs_9item.py` — tables + figures
- `28_blrt_9item.py` — Bootstrap LRT for 9-item (K=1 vs K=2, B=30, p=.032)
- `29_cox_time_varying_robust.py` — Cox time-varying with lifelines `cluster_col=ID, robust=True`; Hessian-based per-wave CIs (kept as a cautionary reference; superseded by script 30)
- `30_cox_cluster_bootstrap.py` — ⭐ ID-level cluster bootstrap (B=300) of time-varying Cox; per-wave HR with percentile CIs; bootstrap interaction-coefficient p-values; Figure 6 BOOT
- `31_lcmm_label_alignment.py` — ⭐ Explicit R lcmm raw-class to manuscript-label alignment
- `33_prepare_mplus_lcga_9item.py` — ⭐ Prepare Mplus 7 K=1-5 LCGA inputs in ASCII run folder
- `34_parse_mplus_lcga_9item.py` — ⭐ Parse Mplus outputs and compare Python/R/Mplus fit and labels
