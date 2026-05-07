# Paper A — Key Results Cheat Sheet (본문 인용 즉시 사용)

> 모든 수치는 `04_results_data/` raw CSV에서 검증 가능. Class 라벨: C1 Stable-Low, C2 Persistent-High, C3 Decreasing, C4 Late-Increasing.

---

## 1. 표본 (Table 1, 04_results_data/comparison_to_10item.txt)

- **N (analytic sample)** = 2,452 (≥3 valid waves of 9-item depression)
- N (original m1 cohort) = 2,590; W7 retention 80.0%
- 성별: 남자 54.3%, 여자 45.7%
- 나이: W1=13세 (Korean age 14, 2018), W7=19세 (2024)

## 2. 척도 신뢰도

- 9-item depression α: **W1 .914 → W7 .898** (모두 ≥ .89)
- 10-item depression (E04 포함, sensitivity) α: W1 .922 → W7 .905
- 9문항 vs 10문항 wave 1 r = **.996** (구성 거의 동일)

## 3. LCGA 모델 fit (Table 2)

| K | logL | BIC | aBIC | Entropy | Min class % |
|---|------|------|------|---------|-------------|
| 1 | -14143.19 | 28317.60 | 28304.89 | 1.000 | 100.0 |
| 2 | -13055.84 | 26174.12 | 26148.71 | 0.674 | 42.0 |
| 3 | -12746.20 | 25586.06 | 25547.93 | 0.719 | 13.5 |
| **4** | **-12551.60** | **25228.07** | **25177.24** | **0.661** | **10.5** |
| 5 | -12485.23 | 25126.55 | 25063.01 | 0.667 | **3.7 (<5%)** |

→ K=4 채택 (parsimony + 5% rule)

**BLRT K=2 vs K=1**: observed LRT = **2,174.7**, B=30 bootstrap null mean=3.7 (max 11.5), all 30 < observed, **p_BLRT = .032** (lowest attainable with B=30, add-one rule). 본문에 "p<.001"로 쓰지 말 것.

**R lcmm 재현**: K=4 BIC 25,228.08 ≈ Python EM 25,228.07 (소수점 일치).

## 4. Class profiles (Table 3)

| Class | π̂ (%) | Modal n | β₀ (intercept) | β₁ (linear) | β₂ (quadratic) | M W1 → W7 |
|-------|--------|---------|----------------|--------------|------------------|-------------|
| C1 Stable-Low | 25.4 | **615** | 1.326 | -0.016 | +0.009 | 1.46 → 1.36 |
| C2 Persistent-High | 10.5 | **245** | 2.613 | -0.043 | -0.039 | 2.39 → 2.13 |
| C3 Decreasing | 29.3 | **726** | 1.932 | -0.100 | +0.005 | 2.28 → 1.68 |
| C4 Late-Increasing | 34.8 | **866** | 1.864 | +0.073 | -0.017 | 1.49 → 1.93 |

총 modal n = 615+245+726+866 = 2,452 ✓

## 5. Baseline (W1) covariates (Table 4)

핵심 효과 크기 (η², ANOVA across 4 classes, 모두 p < .001):
- Depression W1 itself: η² = .47 (하지만 9문항 primary)
- Self-esteem: η² = .25
- Aggression: η² = .21
- Somatic symptoms: η² = .24
- Withdrawal: η² = .15
- Parent warmth: η² = .08
- Parent rejection: η² = .08
- Smartphone dependence: η² = .09

**Skinner 6차원 양육 (W1 평균)** by class:
| Dimension | C1 | C2 | C3 | C4 |
|-----------|------|------|------|------|
| Warmth | 3.52 | 3.10 | 3.18 | 3.48 |
| Rejection | 1.60 | 1.99 | 2.00 | 1.65 |
| Coercion | 2.05 | 2.32 | 2.36 | 2.08 |
| Chaos | 1.87 | 2.33 | 2.27 | 1.91 |
| Autonomy | 3.49 | 3.12 | 3.13 | 3.47 |
| Structure | 3.14 | 2.88 | 2.93 | 3.15 |

**성별 분포** (χ² = 53.8, df=3, p < 10⁻¹¹):
- C1 Stable-Low: 40% female
- C2 Persistent-High: 61% female
- C3 Decreasing: 53% female
- C4 Late-Increasing: 40% female

## 6. Distal outcomes at W7 — BCH 3-step (Table 5)

### Suicidal ideation any (≥2)
| Class | BCH proportion [95% CI] | OR vs C1 [95% CI] |
|-------|---------------------------|----------------------|
| C1 Stable-Low (ref) | 8.0% [4.5, 11.2] | 1.00 (ref) |
| C2 Persistent-High | 47.7% [40.5, 55.6] | 10.58 [6.49, 20.21] |
| C3 Decreasing | 30.4% [25.5, 35.2] | 5.05 [3.21, 9.45] |
| **C4 Late-Increasing** | **56.8% [52.3, 61.8]** | **15.23 [9.57, 29.00]** |

### Suicidal ideation strong (≥3)
| Class | BCH proportion [95% CI] | OR vs C1 [95% CI] |
|-------|---------------------------|----------------------|
| C1 | 1.1% [0.0, 2.6] | 1.00 (ref) |
| C2 | 13.6% [9.0, 18.8] | 14.48 [4.84, 149.14] |
| C3 | 2.4% [0.4, 4.6] | 2.26 [0.38, 23.60] |
| C4 | 11.5% [9.0, 14.4] | 11.97 [4.20, 121.43] |

### Mean depression W7 (BCH)
- C1: 1.17, C2: 2.19, C3: 1.53, C4: 2.04

### Life satisfaction W7 (BCH)
- C1: 2.73, C2: 2.32, C3: 2.53, C4: 2.49

## 7. Time-to-first strong SI — Cox + cluster bootstrap (Table 6, Figure 6)

### Schoenfeld PH test
- C2 vs C1: p = .477 (PH satisfied)
- C3 vs C1: **p < .001 (PH violated)**
- C4 vs C1: **p < .001 (PH violated)**

### Per-wave HR with cluster-bootstrap 95% CI (B=300)

| Wave | C2 Persistent-High | C3 Decreasing | C4 Late-Increasing |
|------|---------------------|------------------|------------------------|
| W1 | 12.31 [8.17, 20.67]*** | 8.96 [6.29, 14.89]*** | **1.43 [0.94, 2.30] (n.s., p=.10)** |
| W2 | 11.14 [8.29, 16.82]*** | 6.54 [4.90, 9.58]*** | 1.88 [1.38, 2.76]*** |
| W3 | 10.09 [7.86, 14.03]*** | 4.77 [3.69, 6.44]*** | 2.48 [1.97, 3.44]*** |
| W4 | 9.14 [7.01, 12.82]*** | 3.48 [2.71, 4.73]*** | 3.26 [2.61, 4.42]*** |
| W5 | 8.28 [5.89, 12.95]*** | 2.54 [1.84, 3.73]*** | 4.30 [3.38, 6.20]*** |
| W6 | 7.49 [4.97, 13.40]*** | 1.85 [1.16, 3.07]** | 5.66 [4.14, 8.81]*** |
| W7 | 6.79 [4.00, 14.29]*** | 1.35 [0.75, 2.55] (n.s., p=.31) | **7.45 [5.00, 12.97]*** |

`*** p = .003 (bootstrap floor 1/301 = .0033 ≈ .003); ** p = .007 (W6 c3 only).`

### Class × wave interaction (cluster bootstrap, B=300)
- c2_x_t: β = -0.099 [-0.247, +0.061], boot p = .22 (n.s. — chronic class HR ~ stable)
- **c3_x_t: β = -0.315 [-0.476, -0.183], boot p = .003** (significant decline)
- **c4_x_t: β = +0.275 [+0.149, +0.413], boot p = .003** (significant rise)

### Cumulative event rate (KM)
- C1: 13.5% experienced strong SI at any wave
- C2: 75.5%
- C3: 49.3%
- C4: 38.8%

## 8. Sensitivity analyses

- **10-item composite**: same 4-class structure; modal-class agreement after Hungarian alignment = 94.9% (κ = .929). Distal SI ORs slightly inflated under 10-item (C4 OR 17.1 vs 9-item 15.2) — direct evidence of criterion contamination.
- **Longitudinal weight (WEIGHTB2w7)**: same trajectory shapes; agreement 95.6%, **κ = .939**. Class structure robust to differential attrition.
- **Modal vs BCH**: modal-class logistic ORs systematically smaller (C4 modal OR ≈ 6 vs BCH OR 15.2) → BCH correction necessary.

---

## 본문 인용 시 즉시 사용할 핵심 한 줄

> "A four-class developmental typology of Korean adolescent depression (Stable-Low 25.4%, Persistent-High 10.5%, Decreasing 29.3%, Late-Increasing 34.8%) was identified across seven annual waves (age 13–19), with the Late-Increasing class — invisible from baseline covariates at age 13 (W1 hazard ratio 1.43 [0.94, 2.30] vs. resilient class) — reaching the highest emerging-adult suicidal-ideation risk by age 19 (BCH-corrected odds ratio 15.2 [9.6, 29.0]; W7 cluster-bootstrap hazard ratio 7.45 [5.00, 12.97], bootstrap p = .003)."
