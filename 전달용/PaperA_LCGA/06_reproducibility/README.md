# Reproducibility — Paper A scripts

전체 분석은 아래 순서로 재실행 가능. 모든 random seed는 42로 고정.

## 환경
- **Python 3.11**: numpy, pandas, scipy, statsmodels, lifelines, scikit-learn, pyreadstat, matplotlib
- **R 4.5.2**: lcmm, haven, dplyr, tidyr
- 원본 데이터: `KCYPS2018m1Yw{1..7}.dta` (NYPI public-use 2025.12), 본 패키지에는 미포함

## 실행 순서

| 단계 | 스크립트 | 설명 | 출력 |
|------|----------|------|------|
| 0 | `06_build_long_panel.py` | 7파동 long-form 패널 + 18개 척도 점수 + α | `long_panel.parquet` |
| 1 | `22_lcga_9item_primary.py` | 9문항 LCGA K=1..5 EM | `class_assignments_K4_9item.csv`, `fit_summary_9item.csv` |
| 2 | `28_blrt_9item.py` | Bootstrap LRT K=1 vs K=2 (B=30) | `blrt/blrt_summary_9item.csv` |
| 3 | `26_lcmm_9item.R` | R lcmm cross-validation | `lcmm_R/lcmm_fit_summary_9item.csv` |
| 4 | `31_lcmm_label_alignment.py` | R↔manuscript class label 매핑 | `lcmm_R/python_vs_R_class_comparison_9item.csv` |
| 5 | `25_weighted_sensitivity_9item.py` | 종단 가중치 sensitivity | `sensitivity_weighted/kappa_aligned.txt` |
| 6 | `23_bch_cox_9item.py` | BCH 3-step + Cox + Schoenfeld | `bch/`, `cox_HR_class_9item.csv`, `cox_PH_schoenfeld_test.csv` |
| 7 | `24_cox_time_varying_9item.py` | Cox time-varying (naive SE; reference) | `cox_time_varying_9item.csv` |
| 8 | `29_cox_time_varying_robust.py` | lifelines `cluster_col=ID, robust=True` (caution: Hessian-based) | `cox_time_varying_9item_ROBUST.csv` |
| 9 | **`30_cox_cluster_bootstrap.py`** ⭐ | **ID-level cluster bootstrap (B=300) — PRIMARY for §3.6** | `cox_HR_per_wave_with_CI_9item_BOOT.csv`, `Figure6_..._BOOT.{png,pdf}` |
| 10 | `27_publication_outputs_9item.py` | Tables 1-5 + Figures 1-5 | `tables/`, `figures/` |
| 11 | `32_rebuild_table6.py` | Table 6 / 6b (cluster bootstrap 기반, note row 자동 추가) | `tables/Table6*.csv`, `Table6_notes.txt` |

## 의존성 그래프
```
0 (long panel)
├─ 1 (LCGA EM) ──┬─ 2 (BLRT)
│                ├─ 3 (R lcmm) ── 4 (label align)
│                ├─ 5 (weighted sensitivity)
│                └─ 6 (BCH + Cox + Schoenfeld) ──┬─ 7 (TV Cox naive)
│                                                  ├─ 8 (TV Cox lifelines robust)
│                                                  └─ 9 (cluster bootstrap) ── 11 (Table 6)
└─ 10 (Tables 1-5, Figures 1-5)
```

## 전체 실행 시간 (단일 CPU 기준)
- 0–11 풀 파이프라인: 약 50–60분
- 가장 무거운 단계: 3 R lcmm (gridsearch 30 starts × K=2..5) ≈ 30분
- 두 번째: 2 BLRT (B=30 parametric bootstrap) ≈ 6분

## 옵션 — 본격 BLRT/Bootstrap (게재 전)
- BLRT B=1000 → p<.001 보고 가능. 추정: ~3–4시간
- Cluster bootstrap B=1000 → floor 1/1001 ≈ .001. 추정: ~2분 (현재 B=300)

## 키 결과 일관성 점검 (디버깅용)

```python
import pandas as pd
df = pd.read_csv('class_assignments_K4_9item.csv')
modal = df['class_K4_9item'].value_counts().sort_index()
# Expected: {1: 615, 2: 245, 3: 726, 4: 866}
assert modal.tolist() == [615, 245, 726, 866], f'Modal counts changed: {modal.tolist()}'

fit = pd.read_csv('fit_summary_9item.csv')
# Expected K=4 BIC = 25228.07
assert abs(fit.loc[fit.K==4, 'BIC'].iloc[0] - 25228.07) < 0.1
```

```r
library(lcmm)
# Expected R lcmm K=4 BIC = 25228.08
```
