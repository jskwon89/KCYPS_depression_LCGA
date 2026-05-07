# Figures — 사용 지침

## ⭐ 투고용 Primary: `R_ggplot2/` 사용

작성담당자께: 투고 시에는 **`R_ggplot2/` 폴더 안의 R/ggplot2 버전**을 사용해주세요. 모든 figure는 6종 (Figure 1-6) 모두 R로 생성됨.
- **600 dpi PNG** (rasterized, review/supplement용)
- **600 dpi TIFF (LZW compressed)** ← JAD/Elsevier 본문 figure 업로드 표준
- (PDF 버전은 작성자 요청에 따라 제거. 필요 시 R 스크립트에서 재생성 가능)

R 버전이 publication 표준에 가까운 이유:
1. ggplot2 그리드 시스템 + APA-friendly classic theme
2. Color-blind safe palette (Okabe-Ito 계열)
3. 600 dpi 고해상도 + vector PDF (zoom 무손실)
4. TIFF LZW (Elsevier 권장)
5. 한글 폰트 의존성 없음 (sans serif fallback)

## 폴더 구성

```
03_figures/
├─ README_FIGURES.md                              (이 문서)
├─ R_ggplot2/                                     ⭐⭐ 투고용 (R/ggplot2, 600 dpi)
│  ├─ Figure1_LCGA_trajectories_R.{png,tiff}        trajectories + bootstrap CI
│  ├─ Figure2_suicidal_ideation_R.{png,tiff}        SI rates + KM curve (3-panel)
│  ├─ Figure3_baseline_heatmap_R.{png,tiff}         W1 covariate z-score heatmap
│  ├─ Figure4_distal_forest_R.{png,tiff}            BCH ORs + Cox HRs forest
│  ├─ Figure5_fit_indices_R.{png,tiff}              BIC/aBIC/AIC + entropy
│  └─ Figure6_time_varying_HR_R.{png,tiff}          Cox cluster bootstrap HR
└─ (Python matplotlib 버전 6개 — 600 dpi 초안 검토용)
   Figure1..6_*_9item.png + Figure6_..._BOOT.png
```

## 도형 설명 (Figure caption 초안)

### Figure 1
> Latent class growth analysis of 9-item depression composite (criterion-clean; YPSY4E04 suicidal ideation excluded) across seven annual waves of the KCYPS-2018 m1 cohort (N=2,452). Solid lines indicate observed class means at each wave; shaded ribbons indicate 95% nonparametric bootstrap confidence intervals (B=500). The shaded clinical range (depression ≥ 2.0) is shown for reference.

### Figure 2
> Suicidal ideation by depression trajectory class. (A) Proportion endorsing any suicidal ideation (item YPSY4E04 ≥ 2 = "sometimes or more") at each wave. (B) Proportion endorsing strong suicidal ideation (≥ 3 = "often"). (C) Kaplan-Meier survival curves for time to first endorsement of strong suicidal ideation across W1–W7. Shaded bands in panel C indicate 95% pointwise log-log CIs.

### Figure 5
> LCGA model fit indices for K=1 through 5 latent classes on the 9-item primary depression composite. (A) Bayesian Information Criterion (BIC), sample-size-adjusted BIC (aBIC), and Akaike Information Criterion (AIC). (B) Classification entropy (≥ 0.7 indicates good class separation; horizontal dotted line). The red dashed line marks the K=4 solution selected as primary on the basis of BIC reduction, entropy, and the 5%-minimum-class-size rule (Nagin, 2005).

### Figure 6
> Time-varying class hazard ratios for first endorsement of strong suicidal ideation, derived from a Cox proportional hazards model with class × wave interaction terms. Reference: C1 Stable-Low. Shaded ribbons indicate 95% percentile confidence intervals from ID-level cluster bootstrap (B=300 resamples of unique individuals). Bootstrap p-values are floored at 1/(B+1) = 1/301 = .003. The C4 Late-Increasing class crosses the C3 Decreasing class around W4–W5 (the high-school transition).

## Re-generation

R 그래프 재생성:
```r
Rscript 06_reproducibility/scripts/33_publication_figures_R.R
```

요구 패키지: `ggplot2, dplyr, tidyr, readr, scales, patchwork, ggsci, cowplot, haven`

원자료 (.dta) 필요 (Figure 1 trajectory 평균 계산용).
