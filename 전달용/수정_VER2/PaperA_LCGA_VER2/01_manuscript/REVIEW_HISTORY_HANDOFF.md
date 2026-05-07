# KCYPS Paper A Review History and Handoff

**Date:** 2026-05-08  
**Project folder:** `D:/2026/SCI/청소년패널조사`  
**Primary output folder:** `output/paperA_lcga`  
**Current manuscript:** `output/paperA_lcga/MANUSCRIPT_DRAFT.md`  

이 문서는 Codex/Chat 기록이 프로젝트 화면에 보이지 않을 경우에도 Paper A 분석과 검토 이력을 이어받을 수 있도록 남기는 핸드오프 메모다.

## 1. 현재 논문 개요

**Paper A:** KCYPS 2018 중1 코호트 1-7차 자료를 이용한 청소년 우울 발달궤적 LCGA와 후기 청소년기/초기 성인기 자살사고 위험 분석.

핵심 메시지:

> 중1 시점에는 stable-low 집단과 거의 구별되지 않는 late-increasing 우울 집단이 고등학교 전환기 이후 급격히 위험군으로 부상하며, W7 자살사고 위험이 가장 높다. 따라서 단일 시점 선별보다 반복적 종단 선별이 필요하다.

주요 타깃 저널:

1. `Journal of Affective Disorders`
2. `Journal of Youth and Adolescence`
3. `Development and Psychopathology`
4. `Suicide and Life-Threatening Behavior`
5. 빠른 게재 대안: `Frontiers in Psychiatry`, `Behavioral Sciences`

## 2. 분석 설계

자료:

- KCYPS 2018 중1 코호트, Wave 1-7.
- 분석표본: `N=2,452`, 3개 이상 우울 파동 유효.
- 원 코호트: `N=2,590`.

주요 측정:

- Primary 우울척도: `YPSY4E01-E03, E05-E10` 9문항 평균.
- `YPSY4E04` 자살사고 문항은 우울 궤적에서 제외하고 distal outcome으로 분리.
- 이 조치는 criterion contamination 방지를 위한 핵심 수정사항이다.

주요 분석:

- LCGA, K=1-5 비교.
- 모형 선택 기준: BIC, aBIC, AIC, entropy, 최소 class 비율, K=1 vs K=2 BLRT.
- Python EM 결과를 R `lcmm`로 재현 검증.
- W7 distal outcome은 BCH 3-step으로 class assignment error 보정.
- W1-W7 strong suicidal ideation 최초 발생은 Cox 생존분석.
- PH 위반이 확인된 C3/C4는 class x wave time-varying Cox로 보완.
- time-varying Cox uncertainty는 ID-level cluster bootstrap `B=300`으로 산출.

## 3. 핵심 결과

4-class LCGA primary solution:

| Class | Label | Estimated proportion | Modal n |
|---|---|---:|---:|
| C1 | Stable-Low / Resilient | 25.4% | 615 |
| C2 | Persistent-High / Inverted-U | 10.5% | 245 |
| C3 | Decreasing | 29.3% | 726 |
| C4 | Late-Increasing | 34.8% | 866 |

W7 suicidal ideation, BCH-corrected OR vs C1:

| Class | OR |
|---|---:|
| C2 Persistent-High | 10.57 |
| C3 Decreasing | 5.05 |
| C4 Late-Increasing | 15.23 |

Time-varying Cox key pattern:

- C2 Persistent-High: W1 HR 12.31 -> W7 HR 6.79, chronic high risk.
- C3 Decreasing: W1 HR 8.96 -> W7 HR 1.35, early risk that dissipates.
- C4 Late-Increasing: W1 HR 1.43, n.s. -> W7 HR 7.45, late-emerging risk.
- C3 x time and C4 x time interactions are significant under cluster bootstrap.
- C4 becomes statistically distinct from C1 from W2 onward and crosses C3 around W4-W5.

## 4. Review Findings Closed

아래 항목들은 reviewer-style 검토에서 발견되었고 현재 산출물 기준으로 처리된 것으로 확인했다.

1. **Criterion contamination fixed**
   - 기존 10문항 우울척도는 `YPSY4E04` 자살사고 문항을 포함했다.
   - Primary 분석을 9문항 우울척도로 전면 교체했다.
   - 10문항 결과는 sensitivity로 내렸다.

2. **Weighted sensitivity kappa fixed**
   - 라벨 정렬 전 kappa가 잘못 보고되는 문제가 있었다.
   - Hungarian label alignment 후 Cohen kappa를 보고하도록 수정했다.

3. **BLRT reporting fixed**
   - 10-item BLRT 숫자 잔재를 제거하고 9-item BLRT로 정정했다.
   - K=1 vs K=2 observed LRT = 2,174.7.
   - B=30, zero exceedance, add-one p = 1/31 = .032.
   - K>2 BLRT는 수행하지 않았고, K=4 선택은 BIC 감소, entropy, 최소 class 비율, 해석가능성 기준으로 명시했다.

4. **BLRT add-one floor precision fixed**
   - B=999 gives exactly p=.001.
   - p<.001 requires B>=1000.
   - Methods, Results, version note가 이 기준으로 정리됐다.

5. **R lcmm replication completed**
   - R K=4 BIC = 25,228.08.
   - Python EM K=4 BIC = 25,228.07.
   - 라벨 정렬 후 modal count 차이는 class별 0-1명 수준.
   - R raw label mapping:
     - R class 1 -> C1 Stable-Low
     - R class 3 -> C2 Persistent-High
     - R class 2 -> C3 Decreasing
     - R class 4 -> C4 Late-Increasing

6. **Modal n vs estimated class proportion separated**
   - 본문에서 mixing proportion과 modal assignment count를 구분했다.

7. **Cox PH assumption handled**
   - Schoenfeld test 산출물을 추가했다.
   - C2는 PH assumption tenable.
   - C3/C4는 PH violation; time-varying Cox로 substantive pattern을 해석했다.

8. **Time-varying Cox uncertainty fixed**
   - 기존 `variance_matrix_` 기반 CI는 Hessian covariance라 robust CI가 아니었다.
   - ID-level cluster bootstrap `B=300`으로 per-wave HR CI와 interaction p를 재산출했다.

9. **Cox bootstrap p-value floor fixed**
   - B=300이면 add-one floor는 `1/301=.00332`.
   - 원고에서 `p<=.003` 표현을 제거하고 `p=.003` rounded floor로 통일했다.
   - footnote를 추가해 discrete bootstrap floor를 설명했다.

10. **Table 6 rebuilt**
    - 기존 Table 6은 average HR와 Schoenfeld summary 위주였다.
    - 현재 Table 6은 W1/W4/W7 per-wave HR, cluster-bootstrap CI, bootstrap interaction p를 포함한다.
    - Table 6b는 W1-W7 전체 per-wave HR를 supplementary table로 제공한다.

## 5. 주요 산출물 위치

Manuscript and reports:

- `output/paperA_lcga/MANUSCRIPT_DRAFT.md`
- `output/paperA_lcga/PROGRESS_REPORT_v1.md`
- `output/paperA_lcga/PROGRESS_REPORT_v2.md`
- `output/paperA_lcga/PROGRESS_REPORT_v3.md`
- `output/paperA_lcga/PROGRESS_REPORT_v4.md`
- `output/paperA_lcga/REVIEW_HISTORY_HANDOFF.md`

Primary 9-item outputs:

- `output/paperA_lcga/primary_9item/tables/Table1_sample_characteristics_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table2_model_fit_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table3_class_profiles_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table4_baseline_by_class_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table5_BCH_distal_outcomes_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table6_cox_HR_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table6b_cox_HR_per_wave_full_9item.csv`
- `output/paperA_lcga/primary_9item/tables/Table6_notes.txt`

Key technical outputs:

- `output/paperA_lcga/primary_9item/blrt/blrt_summary_9item.csv`
- `output/paperA_lcga/primary_9item/lcmm_R/lcmm_fit_summary_9item.csv`
- `output/paperA_lcga/primary_9item/lcmm_R/python_vs_R_class_comparison_9item.csv`
- `output/paperA_lcga/primary_9item/lcmm_R/class_label_alignment_9item.csv`
- `output/paperA_lcga/primary_9item/mplus/MPLUS_REPLICATION_REPORT.md`
- `output/paperA_lcga/primary_9item/mplus/mplus_vs_python_fit_summary_9item.csv`
- `output/paperA_lcga/primary_9item/mplus/python_vs_mplus_class_comparison_9item.csv`
- `output/paperA_lcga/primary_9item/long_distal/cox_HR_per_wave_with_CI_9item_BOOT.csv`
- `output/paperA_lcga/primary_9item/long_distal/cox_interaction_bootstrap_dist_9item.csv`
- `output/paperA_lcga/primary_9item/figures/Figure6_time_varying_HR_9item_BOOT.png`

Key scripts:

- `scripts/22_lcga_9item_primary.py`
- `scripts/23_bch_cox_9item.py`
- `scripts/25_weighted_sensitivity_9item.py`
- `scripts/26_lcmm_9item.R`
- `scripts/27_publication_outputs_9item.py`
- `scripts/28_blrt_9item.py`
- `scripts/30_cox_cluster_bootstrap.py`
- `scripts/31_lcmm_label_alignment.py`
- `scripts/32_rebuild_table6.py`
- `scripts/33_prepare_mplus_lcga_9item.py`
- `scripts/34_parse_mplus_lcga_9item.py`

## 5a. Mplus Replication Added

2026-05-08에 Mplus 7로 primary 9-item LCGA K=1-5를 재현했다.

Mplus model:

- 9-item depression composite, E04 excluded.
- Analytic N = 2,452.
- Quadratic LCGA with time scores -3,-2,-1,0,1,2,3.
- Growth-factor variances fixed to zero: `i@0; s@0; q@0;`.
- One residual variance constrained equal across all seven waves and classes.
- `STARTS = 1000 250`; K=1-5 모두 best loglikelihood replicated.

Fit replication:

| K | Mplus logL | Python logL | Mplus BIC | Python BIC | Entropy |
|---:|---:|---:|---:|---:|---:|
| 1 | -14143.189 | -14143.188 | 28317.597 | 28317.595 | NA |
| 2 | -13055.844 | -13055.844 | 26174.126 | 26174.124 | .675 |
| 3 | -12746.201 | -12746.200 | 25586.058 | 25586.056 | .719 |
| 4 | -12551.601 | -12551.600 | 25228.076 | 25228.075 | .661 |
| 5 | -12485.230 | -12485.230 | 25126.554 | 25126.553 | .667 |

Mplus K=4 raw label mapping:

- Mplus raw class 2 -> C1 Stable-Low
- Mplus raw class 4 -> C2 Persistent-High
- Mplus raw class 1 -> C3 Decreasing
- Mplus raw class 3 -> C4 Late-Increasing

Python vs Mplus modal class counts after alignment:

| Label | Python n | Mplus n | Difference |
|---|---:|---:|---:|
| C1 Stable-Low | 615 | 616 | +1 |
| C2 Persistent-High | 245 | 245 | 0 |
| C3 Decreasing | 726 | 725 | -1 |
| C4 Late-Increasing | 866 | 866 | 0 |

Agreement after alignment:

- Percent agreement = 99.96%.
- Cohen's kappa = .9994.

Interpretation:

- Mplus, R lcmm, and Python EM all recover the same four-class solution.
- This is now a strong cross-platform replication suitable for JAD/JYA reviewer defense.

## 6. 현재 판단

분석 패키지는 JAD 혹은 JYA에 투고해도 부끄럽지 않은 수준까지 왔다.

강점:

- 전국 패널 7파동.
- 우울 trajectory와 suicidal ideation의 임상적 연결.
- late-increasing invisible risk group이라는 명확한 narrative.
- criterion contamination, class uncertainty, PH violation, bootstrap p floor 등 reviewer 공격 지점을 선제적으로 처리.
- Python/R cross-platform replication.

남은 승부처:

- Discussion 완성도.
- clinical/public-health framing.
- JAD 스타일에 맞춘 concise abstract와 implications.
- limitation을 선제적으로 인정하는 균형감.

현재 게재 준비도:

- 분석/산출물 기준: 약 80%대 초중반.
- 원고 작성과 Discussion polish가 잘 되면 JAD 도전 가능.
- 빠른 게재가 더 중요하면 `Journal of Youth and Adolescence`, `Frontiers in Psychiatry`, `Behavioral Sciences`도 현실적 대안.

## 7. 다음 할 일

1. `scripts/32_rebuild_table6.py`가 현재 Table 6/6b note row 또는 `Table6_notes.txt`까지 재현하는지 최종 확인.
2. `MANUSCRIPT_DRAFT.md`의 Discussion을 1,500단어 내외로 정리.
3. Abstract를 JAD 스타일로 압축.
4. References 실제 DOI/형식 점검.
5. Tables/Figures를 journal formatting에 맞게 정리.
6. 최종 투고 전, 원고 전체에서 아래 표현 검색:
   - `10-item primary`
   - `p<=.003`
   - `p ≤ .003`
   - `B>=999`
   - `B ≥ 999`
   - `robust CI` without bootstrap clarification
   - `BLRT for K>2`
