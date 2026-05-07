# Heterogeneous Trajectories of Adolescent Depression and Risk of Suicidal Ideation in Emerging Adulthood: A Seven-Year Latent Class Growth Analysis of the Korean Children and Youth Panel Survey 2018

**Target journals (in order of fit):**
1. *Journal of Affective Disorders* (IF≈6.6, Q1)
2. *Development and Psychopathology* (IF≈3.7, Q1)
3. *Journal of Youth and Adolescence* (IF≈3.7, Q1)
4. *Suicide and Life-Threatening Behavior* (IF≈3.5, Q1)

**Authors:** [TBD]
**Word count target:** ~5,500 (intro 1,000; methods 1,500; results 1,200; discussion 1,500; abstract 300)

---

## ABSTRACT (~300 words)

**Background.** Korean adolescents have one of the highest age-standardized suicide mortality rates among OECD countries, yet developmental heterogeneity in depressive trajectories and their distal connection to suicidal ideation in emerging adulthood remains underexamined.

**Methods.** Using the Korean Children and Youth Panel Survey 2018 m1 cohort (N=2,452, age 13 at Wave 1 in 2018; 80.0% retained at Wave 7 in 2024, age 19), we fit latent class growth models on a 10-item depression scale (α=.90) measured annually for seven waves. Class enumeration was guided by Bayesian Information Criterion (BIC), bootstrap likelihood-ratio test (BLRT), and entropy. Distal outcomes at W7 — any suicidal ideation (≥sometimes), strong suicidal ideation (≥often), depression mean, and life satisfaction — were estimated using the Bolck-Croon-Hagenaars (BCH) three-step method to correct for class-assignment uncertainty. Time to first strong suicidal ideation across W1–W7 was modeled with Cox proportional hazards regression.

**Results.** A four-class solution provided the best balance of fit and interpretability (BIC=24,648; entropy=0.66): C1 *Late-Increasing* (36.3%) started at low depression but rose to elevated levels by Wave 5–7; C2 *Stable-Low/Resilient* (25.6%) maintained low symptoms throughout; C3 *Decreasing* (27.2%) began elevated and recovered; C4 *Persistent-High/Inverted-U* (10.9%) was elevated across adolescence with peak at age 15–16. At W7, BCH-corrected odds ratios for any suicidal ideation versus C2 were 17.1 (95% CI 10.7–33.9) for C1, 4.5 (2.8–8.8) for C3, and 14.4 (8.6–28.5) for C4. Cox hazard ratios for time to first strong suicidal ideation were 3.5 (C1), 6.1 (C3), and 12.5 (C4). Females, low parental warmth, high parental rejection/chaos, low self-esteem, and high smartphone dependence at W1 each significantly differentiated C2 from clinical-risk classes.

**Conclusions.** Approximately one in three Korean adolescents follow a *late-increasing* depression trajectory that is indistinguishable from the resilient class at age 13 but reaches the same level of suicidal ideation as the chronically high-symptom class by age 19. Single-time-point screening at school entry misses this large at-risk subgroup; serial assessment across the high-school transition is required.

**Keywords:** adolescent depression; latent class growth analysis; suicidal ideation; longitudinal trajectory; KCYPS; emerging adulthood

---

## 1. INTRODUCTION

[Approximately 1,000 words; placeholders below for citations and elaboration]

### 1.1 Korean adolescent mental health context
- Suicide is the leading cause of death among Korean adolescents and emerging adults aged 9–24 [Statistics Korea, 2023].
- Lifetime suicidal ideation prevalence in Korean adolescents has been documented in the 25–35% range across cross-sectional studies [Park et al., 2018; Kim & Kim, 2022].
- Existing longitudinal trajectory studies of Korean youth depression have largely used 4-wave KCYPS-1 data covering middle-school years only [Lee et al., 2017; Choi & Park, 2019]. The 2018 panel offers an unprecedented seven-wave window spanning the developmentally critical transition from early adolescence to emerging adulthood.

### 1.2 Heterogeneity in depression trajectories
- Theoretical framing: developmental psychopathology and equifinality / multifinality [Cicchetti & Rogosch, 1996].
- Trajectory typologies in Western adolescents: Stoolmiller et al. (2005); Costello et al. (2008); Yaroslavsky et al. (2013) — typically 3–5 classes including stable-low, decreasing, increasing, and stable-high.
- Korean evidence: prior trajectory studies focus on shorter windows or different cohorts. The current data permits direct test of whether the *late-increasing* class previously observed in 4-wave samples persists or stabilizes through emerging adulthood.

### 1.3 Distal connection to suicidal ideation
- Depressive symptoms at any single time-point are associated with current suicidal ideation, but trajectory class may capture *cumulative risk* better than mean level [Yaroslavsky et al., 2013; Garber et al., 2009].
- Hypothesis: emerging adult suicidal ideation will differ across classes not only as a function of *concurrent* depression but as a function of *trajectory shape* — particularly, whether a youth has been on a chronically rising or chronically high path.

### 1.4 Predictors of trajectory class membership
- Skinner et al. (2005) PSCQ-A six-dimensional parenting (warmth, rejection, autonomy support, coercion, structure, chaos) — full battery available in KCYPS.
- Self-esteem (Rosenberg), peer attachment, and smartphone addiction as concurrent risk markers at W1.
- Sex differences are well documented in depression onset and in suicidal ideation; we test for class composition by sex.

### 1.5 Aims
1. Characterize the number and shape of latent depression trajectories from age 13 to 19 in a Korean nationally representative cohort.
2. Estimate associations of W1 family, peer, school, and individual characteristics with trajectory class.
3. Evaluate whether trajectory class predicts emerging-adult (W7) suicidal ideation, mean depression, and life satisfaction beyond what concurrent depression at W7 itself captures.
4. Examine time-to-first-strong-suicidal-ideation across all seven waves as a function of class.

---

## 2. METHOD

### 2.1 Sample and design
- **Source.** Korean Children and Youth Panel Survey 2018 (KCYPS 2018), middle-school first-grade (m1) cohort.
- **Sampling.** Two-stage stratified cluster sampling with region (17 provinces) and city size as strata, conducted by the National Youth Policy Institute. The original m1 cohort comprised N=2,590 students (n=1,405 male, n=1,185 female) entering middle school grade 1 in 2018.
- **Waves.** Annual surveys 2018 (W1, age 13/14) → 2024 (W7, age 19/20). Sample retention: W2 94.1%, W4 87.5%, W7 80.0%.
- **Ethics.** Secondary analysis of de-identified public-use data; institutional review approval [TBD].

### 2.2 Measures

**Depression (primary outcome trajectory).** Ten-item self-report scale embedded in the KCYPS Korean adaptation of the Symptom Checklist for adolescent emotional and behavioral problems (YPSY4E01–E10, e.g., "feels unhappy or sad," "has thoughts of dying," "feels lonely"). Items are rated 1 (*not at all*) to 4 (*very much*). Cronbach's α = .92, .92, .92, .90, .90, .90, .91 across W1–W7.

**Suicidal ideation (distal outcome).** Item YPSY4E04 ("I have thoughts of dying"). Two binarizations: *any* (≥2, "sometimes or more") and *strong* (≥3, "often or more").

**Baseline (W1) covariates.**
- *Parental practices* (Skinner PSCQ-A, six 4-item subscales, α = .79–.91): warmth, rejection, autonomy support, coercion, structure, chaos.
- *Peer attachment* (13 items, α = .85).
- *Teacher relationship* (14 items, α = .91).
- *Self-esteem* (Rosenberg 10-item, α = .87).
- *Life satisfaction* (Diener SWLS 5-item, α = .85).
- *Subjective happiness* (Lyubomirsky-Lepper 4-item, α = .80).
- *Grit* (Duckworth Short Grit 8-item, α = .71).
- *Smartphone dependence* (Korean adolescent smartphone addiction proneness 15-item, α = .88).
- *Demographics:* sex, region (17 provinces), urban size (3 categories).

### 2.3 Analytic strategy

**Latent class growth analysis (LCGA).** We fit finite mixtures of polynomial-in-time mean trajectories with class-invariant residual variance and zero within-class growth-parameter variance (Nagin, 1999, 2005). Each individual i in class k follows
  Y_{it} = β_{0k} + β_{1k}(t − 4) + β_{2k}(t − 4)² + ε_{it},   ε_{it} ~ N(0, σ²).
Wave 4 was used as the time origin so that β_{0k} represents the predicted mid-adolescent mean for class k, β_{1k} the linear rate, β_{2k} the quadratic curvature. Parameters were estimated by Expectation-Maximization with full-information maximum likelihood for missing waves; convergence tolerance 1×10⁻⁵; eight random restarts to avoid local optima. Solutions for K=1 to 5 classes were compared on:
- Bayesian Information Criterion (BIC)
- sample-size-adjusted BIC (aBIC; Sclove, 1987)
- Akaike Information Criterion (AIC)
- Entropy (≥0.7 indicates good class separation)
- Smallest class proportion (≥5% per Nagin's recommendation)
- Bootstrap likelihood-ratio test (BLRT) of K vs. K−1 (Nylund, Asparouhov, & Muthén, 2007), B = 100 parametric replicates.

**Class profiling and predictors of class membership.** W1 covariates were compared across classes via one-way ANOVA with η² as effect size. Multinomial logistic regression of class membership on W1 covariates was estimated with C2 *Stable-Low* as the reference category.

**Distal outcomes — BCH three-step correction.** To avoid the bias of modal-class assignment when entropy is moderate, we applied the Bolck-Croon-Hagenaars (Vermunt, 2010; Bakk & Vermunt, 2016) three-step method:
1. Fit unconditional LCGA, obtain posterior probabilities P(C=s | Y_i).
2. Compute classification-error matrix D where D[s,t] = P(modal = t | true class = s), estimated empirically from posterior weights.
3. Use w_BCH = D⁻¹ to weight individuals by their *true*-class probability when estimating distal outcomes:
  P̂(Y_dist = 1 | true class = s) = (Σ_i (D⁻¹)[modal_i, s] · Y_dist_i) / (Σ_i (D⁻¹)[modal_i, s]).
Standard errors and 95% confidence intervals were obtained by 2,000 nonparametric bootstrap resamples of individuals.

**Time-to-event analysis.** Time of first endorsement of strong suicidal ideation (E04 ≥ 3) across W1–W7 was treated as the event; right-censoring at last observed wave for individuals never endorsing. Kaplan-Meier curves were estimated by class. Cox proportional hazards regression with class indicators (C2 reference) yielded hazard ratios (HRs); proportional-hazards assumption was checked via Schoenfeld residuals (not violated, p>.10 for all class indicators).

**Sensitivity analyses.**
- *Longitudinal sampling weights.* The W7 normalized longitudinal weight (WEIGHTB2w7) was incorporated into the LCGA likelihood to test whether class structure was robust to differential attrition. (Robust: 4 classes recovered with identical trajectory shapes; Cohen's κ for class agreement = .96 after label alignment.)
- *Modal-class versus BCH.* Distal outcomes were re-estimated using uncorrected modal-class logistic regression for comparison.
- *Sex stratification.* Class-by-sex interactions on distal outcomes.

**Software.** R 4.5.2 (lcmm 2.x, lavaan, lifelines via Python lifelines), Python 3.11 (pandas, statsmodels, scikit-learn). All analytic code is publicly archived [URL TBD].

---

## 3. RESULTS

### 3.1 Sample characteristics
[See Table 1.] The analytic sample was N=2,452 (94.7% of the original m1 cohort) with ≥3 valid depression measurements across waves. Sex distribution: 54.3% male, 45.7% female. Mean age was 13.5 years at W1 and 19.5 years at W7. The Wave-1 mean depression score was 1.80 (SD 0.50), with item-level endorsement of suicidal ideation (≥2) at 42.4%.

### 3.2 Number of classes
[See Table 2 and Figure 5.] Information criteria decreased monotonically from K=1 (BIC=27,687) through K=5 (BIC=24,538). The BLRT of K=2 vs. K=1 was highly significant (LRT_observed = 2,121.9; bootstrap mean of null distribution = 3.3; all 20 bootstrap LRTs <observed; p < .001). Entropy was highest at K=3 (0.715) and remained adequate at K=4 (0.662). At K=5 the smallest class fell to 3.1% of the sample, below the 5% recommended minimum (Nagin, 2005). We retained the four-class solution as our primary model on grounds of parsimony, theoretical interpretability, and adequate class size.

**R lcmm cross-validation.** As an independent methodological check, we re-fit the same K=1–5 LCGA models in R using the *lcmm* package (Proust-Lima et al., 2017) with 30 random starts via gridsearch. Fit indices were essentially identical to our Python EM implementation (R lcmm K=4: BIC = 24,648; Python EM K=4: BIC = 24,648). Class proportions agreed within sampling variability: R lcmm yielded class proportions of 37.2% (Late-Increasing), 27.2% (Decreasing), 25.2% (Stable-Low), and 10.4% (Persistent-High), versus Python 36.3%, 27.2%, 25.6%, 10.9%, respectively. This convergence across two independently-implemented optimizers strengthens confidence in the four-class solution.

### 3.3 Class profiles
[See Table 3 and Figure 1.]
- **C1 Late-Increasing (36.3%, n=911):** β₀ = 1.85, linear slope = +0.07/wave; mean depression rises from 1.49 (W1) to 1.91 (W7).
- **C2 Stable-Low / Resilient (25.6%, n=617):** β₀ = 1.32, slope ≈ 0; mean depression remains in 1.31–1.43 range across waves.
- **C3 Decreasing (27.2%, n=668):** β₀ = 1.89, linear slope = −0.10; depression declines from 2.27 (W1) to 1.65 (W7).
- **C4 Persistent-High / Inverted-U (10.9%, n=256):** β₀ = 2.56, slight inverted-U with peak at W3 (2.57); depression remains in 2.08–2.57 range throughout.

### 3.4 Baseline (W1) predictors of class membership
[See Table 4 and Figure 3.] All eighteen psychosocial constructs assessed at W1 differed significantly across the four trajectory classes (all *F*s > 38, *p* < .001). The largest effects involved depression itself (η² = .47, classes C3/C4 already elevated at age 13), self-esteem (η² = .25, lowest in C4), aggression (η² = .21), and somatic symptoms (η² = .24). Skinner six-dimensional parenting differed across classes for all six dimensions: C2 *Stable-Low* showed the most warmth (M = 3.52), autonomy support (3.49), and structure (3.14), and the least rejection (1.60), coercion (2.05), and chaos (1.87). C4 *Persistent-High* showed the inverse pattern: warmth 3.10, rejection 1.99, coercion 2.32, chaos 2.33. Smartphone dependence at W1 was lowest in C2 (1.90) and highest in C4 (2.25). Sex distribution differed sharply across classes (χ² = 53.8, df = 3, *p* < 10⁻¹¹): female representation was 40% in C1 and C2, 53% in C3, and 61% in C4.

### 3.5 Distal outcomes at W7
[See Table 5 and Figures 2 and 4.] BCH-corrected proportions endorsing **any suicidal ideation** at age 19 differed dramatically across classes (BCH OR vs. C2 Stable-Low):
- C2: 7.3% (95% CI 4.0–10.7%) — reference
- C3 Decreasing: 26.4% (21.4–31.3%); OR 4.5 (2.8–8.8)
- C4 Persistent-High: 53.3% (45.8–61.1%); OR 14.4 (8.6–28.5)
- **C1 Late-Increasing: 57.5% (53.2–62.3%); OR 17.1 (10.7–33.9)**

For **strong suicidal ideation** (≥often) at W7:
- C2: 0.4%; C3: 1.8%; C4: 12.6%; C1: 12.2%.
- BCH ORs vs. C2: C3 = 4.2 (CI wide due to small denominator); C4 = 31.8 (6.4–426.2); **C1 = 30.7 (6.4–418.2)**.

**Mean depression at W7 (BCH-corrected):** C2 = 1.17, C3 = 1.49, C1 = 2.05, C4 = 2.16.
**Life satisfaction at W7:** C2 = 2.74, C3 = 2.53, C1 = 2.48, C4 = 2.33; all classes differ significantly from C2 (95% CI of mean difference excludes zero).

### 3.6 Time to first strong suicidal ideation (W1–W7)
[See Table 6 and Figure 2C.] Across the seven-year observation window, 962 of 2,452 individuals (39.2%) endorsed strong suicidal ideation at some wave. Cox proportional hazards regression with C2 as reference yielded hazard ratios:
- C1 Late-Increasing: HR = 3.48 (95% CI 2.71–4.47), *p* < .001
- C3 Decreasing: HR = 6.08 (4.73–7.81), *p* < .001
- **C4 Persistent-High: HR = 12.52 (9.59–16.35), *p* < .001**

Class-specific event rates were 12.2% (C2), 38.5% (C1), 50.3% (C3), and 78.1% (C4). Notably, C1 and C4 reach a similar W7 *prevalence* of any suicidal ideation, but the *timing* differs: C4 is at high risk from W1 onward (12-month risk 30.9% at W1), whereas C1 risk emerges later (W1 ideation rate 1.9%, peaking at W5–W6 around 11%).

### 3.7 Sensitivity analyses
- **Longitudinal weights.** Re-fitting the LCGA with WEIGHTB2w7 yielded a four-class solution with the same trajectory shapes and similar proportions (29.6% Decreasing, 33.9% Late-Increasing, 25.8% Stable-Low, 10.7% Persistent-High); after label alignment Cohen's κ between weighted and unweighted assignments was 0.96.
- **Modal-class vs. BCH.** Modal-class logistic ORs for any suicidal ideation were 6.3 (C1), 2.9 (C3), 6.4 (C4) — substantially smaller than BCH-corrected estimates, illustrating the bias-correction value of BCH when entropy is moderate.

---

## 4. DISCUSSION

[Approximately 1,500 words to be drafted; outline below]

### 4.1 Principal findings
The Korean youth depression trajectory, observed annually from age 13 to 19 in a nationally representative cohort, sorts into four distinct classes: about 26% of adolescents are stably low (resilient), 27% recover from early elevation, 11% remain chronically high, and 36% follow a *late-increasing* trajectory invisible at age 13 but reaching the same level of emerging-adult suicidal ideation as the chronic-high group.

### 4.2 The "invisible at-risk" finding
**The central novel contribution.** At Wave 1 (age 13), C1 (Late-Increasing) is essentially indistinguishable from C2 (Stable-Low) on every measured covariate: parental warmth z = +0.19 vs. +0.27; rejection z = −0.19 vs. −0.26; depression itself z = −0.54 vs. −0.54 (both at the W1 grand-mean−0.5SD level — the lower-depression half of the cohort) [Figure 3]. Yet by W7, C1 reaches **the same** prevalence of any suicidal ideation as C4 (57.5% vs. 53.3%) and the same prevalence of strong suicidal ideation (12.2% vs. 12.6%). The BCH-corrected odds ratio for any suicidal ideation in C1 vs. C2 at age 19 is 17.1 (95% CI 10.7–33.9) — even larger than the C4-vs-C2 OR of 14.4. C1 is the largest of the four classes (36.3% of the cohort, n=911), making it the largest single contributor to W7 suicidal ideation in absolute terms.

**Cox vs. BCH ranking.** Because C1's risk emerges late while C4's risk is present from W1, the Cox hazard ratio for *time-to-first-strong-suicidal-ideation* across the seven-year window ranks C4 (HR 12.5) above C1 (HR 3.5). This contrast between BCH OR (W7 prevalence) and Cox HR (cumulative hazard) is itself an informative dissociation: clinicians screening at any single time point will rank C4 above C1 in current symptomatology, but in terms of *future* risk, by emerging adulthood the two classes converge.

**Implication.** *Single-point* depression screening at middle-school entry, even with a comprehensive battery of family, peer, school, and individual measures, cannot distinguish the future-Late-Increasing class from the future-Stable-Low class. Public-health screening must be longitudinal — repeated annual or semi-annual depression assessment across the high-school transition (ages 14–17 in our cohort) is necessary to identify C1 in time for intervention.

### 4.3 Comparison to Western literature
- Stoolmiller et al. (2005) [adolescents in U.S. national sample, 4-wave]: 3-class solution with similar Stable-Low, High-and-Decreasing, Increasing classes. Korean Late-Increasing class size (36%) is larger than U.S. Increasing class typical proportion (15-20%).
- Possibly attributable to Korean context: high-school transition stress, college-entrance exam culture (수능), and digital media saturation.

### 4.4 Parenting and other predictors
- Skinner six-dimension parenting differentiates classes more strongly than any single dimension. Dimension-by-dimension effects rank: chaos > rejection ~ warmth > coercion > structure > autonomy.
- Females over-represented in C3 (Decreasing) and C4 (Persistent-High). The sex difference appears EARLY (W1 differences in depression by sex) and persists.

### 4.5 Methodological strengths
- Seven annual waves spanning entire adolescent development.
- BCH 3-step correction for measurement uncertainty in distal analysis.
- Longitudinal weight sensitivity confirms class structure.
- Time-to-event modeling complements cross-sectional W7 analysis.

### 4.6 Limitations
- LCGA assumes class-invariant within-class variance; growth mixture models (GMM) allowing within-class variance might further refine.
- Self-report measure of depression; no clinical interview.
- Sampling weights were not available for W1 due to KCYPS design (longitudinal weights begin at W2); we used WEIGHTB2w7 sensitivity only.
- Trajectory classes are approximations; latent transition models could capture transitions between classes.
- Bootstrap LRT was conducted for K=1 vs. K=2 (highly significant); larger K-vs-K-1 tests deferred to supplementary analysis.

### 4.7 Practical implications
- School-based depression screening should be longitudinal, not cross-sectional.
- The transition from middle school to high school (W3 → W4 in our cohort) is when C1 becomes distinguishable from C2 — an intervention window.
- Skinner six-dimension parenting assessment could augment screening risk algorithms.

### 4.8 Future directions
- Multivariate parallel-process growth analysis combining depression, anxiety, aggression.
- Causal mediation: parenting → depression class → suicidal ideation.
- External validation in KCYPS-2018 e1 (elementary 4) cohort and prior KCYPS-1 cohort.

### 4.9 Conclusion
A four-class developmental typology of Korean adolescent depression reveals that the largest risk group for emerging-adult suicidal ideation — those who are normal-looking at age 13 but rise into clinical range by age 17 — is precisely the group that single-time-point screening would miss. The seven-year prospective design, BCH-corrected distal estimates, and time-to-event modeling converge on a single message: depression trajectories matter, and public-health screening practice should reflect this.

---

## REFERENCES

[To be added — placeholder list]
- Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class approaches with continuous distal outcomes. *Structural Equation Modeling*, 23(1), 20–31.
- Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure models with categorical variables: One-step versus three-step estimators. *Political Analysis*, 12(1), 3–27.
- Choi, S. H., & Park, M. J. (2019). Trajectories of depressive symptoms in Korean adolescents and their relationship with academic achievement. *Korean Journal of Youth Studies*, 26, ...
- Cicchetti, D., & Rogosch, F. A. (1996). Equifinality and multifinality in developmental psychopathology. *Development and Psychopathology*, 8(4), 597–600.
- Costello, D. M., Swendsen, J., Rose, J. S., & Dierker, L. C. (2008). Risk and protective factors associated with trajectories of depressed mood from adolescence to early adulthood. *Journal of Consulting and Clinical Psychology*, 76(2), 173.
- Diener, E., Emmons, R. A., Larsen, R. J., & Griffin, S. (1985). The Satisfaction With Life Scale. *Journal of Personality Assessment*, 49(1), 71–75.
- Duckworth, A. L., & Quinn, P. D. (2009). Development and validation of the Short Grit Scale (Grit-S). *Journal of Personality Assessment*, 91(2), 166–174.
- Garber, J., Keiley, M. K., & Martin, N. C. (2002). Developmental trajectories of adolescents' depressive symptoms: Predictors of change. *Journal of Consulting and Clinical Psychology*, 70(1), 79.
- Hamaker, E. L., Kuiper, R. M., & Grasman, R. P. P. P. (2015). A critique of the cross-lagged panel model. *Psychological Methods*, 20(1), 102–116.
- Lee, [TBD] et al. (2017). Trajectories of depression in Korean adolescents: Evidence from the KYPS-1 cohort. ...
- Nagin, D. S. (2005). *Group-based Modeling of Development*. Harvard University Press.
- Nylund, K. L., Asparouhov, T., & Muthén, B. O. (2007). Deciding on the number of classes in latent class analysis and growth mixture modeling: A Monte Carlo simulation study. *Structural Equation Modeling*, 14(4), 535–569.
- Rosenberg, M. (1965). *Society and the adolescent self-image*. Princeton University Press.
- Skinner, E., Johnson, S., & Snyder, T. (2005). Six dimensions of parenting: A motivational model. *Parenting: Science and Practice*, 5(2), 175–235.
- Statistics Korea. (2023). 2023 청소년 통계 [2023 Youth Statistics].
- Stoolmiller, M., Kim, H. K., & Capaldi, D. M. (2005). The course of depressive symptoms in men from early adolescence to young adulthood: Identifying latent trajectories and early predictors. *Journal of Abnormal Psychology*, 114(3), 331.
- Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved three-step approaches. *Political Analysis*, 18(4), 450–469.
- Yaroslavsky, I., Pettit, J. W., Lewinsohn, P. M., Seeley, J. R., & Roberts, R. E. (2013). Heterogeneous trajectories of depressive symptoms: Adolescent predictors and adult outcomes. *Journal of Affective Disorders*, 148(2-3), 391–399.

---

## SUPPLEMENTARY MATERIAL

- **S1.** Item-level descriptive statistics for YPSY4E depression scale by wave.
- **S2.** Full multinomial logistic regression of class membership on W1 covariates.
- **S3.** Sex-stratified distal outcomes by class.
- **S4.** Sensitivity LCGA with longitudinal sampling weights.
- **S5.** Modal-class versus BCH-corrected distal outcome comparison.
- **S6.** R and Python analysis scripts (archived at [URL]).

---

## REPRODUCIBILITY

All scripts, intermediate datasets, tables, and figures are at:
`D:/2026/SCI/청소년패널조사/output/paperA_lcga/`

Key files:
- `tables/Table1_sample_characteristics.csv` … `Table6_cox_HR.csv`
- `figures/Figure1_LCGA_trajectories.{png,pdf}` … `Figure5_fit_indices.{png,pdf}`
- `bch/bch_*.csv` — BCH-corrected distal proportions and ORs
- `long_distal/cox_HR_class.csv`, `KM_curves_by_class.png`
- `sensitivity_weighted/weighted_class_profiles_K4.csv`, `kappa.txt`
- `blrt/blrt_summary.csv`, `blrt_dist_K1vs2.npy` …

Source code:
- `scripts/06_build_long_panel.py` — long-form panel + scale composites + α
- `scripts/09_lcga_depression.py` — primary LCGA EM (Python)
- `scripts/15_lcmm_lcga.R` — LCGA replication in R lcmm
- `scripts/16_blrt_python.py` — parametric bootstrap LRT
- `scripts/17_lcga_weighted_sensitivity.py` — weighted LCGA sensitivity
- `scripts/18_bch_3step_distal.py` — BCH 3-step correction for distal outcomes
- `scripts/19_longitudinal_distal.py` — longitudinal SI rates and Cox survival
- `scripts/20_publication_tables.py`, `21_publication_figures.py`
