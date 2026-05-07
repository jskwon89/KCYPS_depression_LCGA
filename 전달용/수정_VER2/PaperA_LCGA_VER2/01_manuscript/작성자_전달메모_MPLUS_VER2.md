# 작성자 전달메모: Mplus 이후 VER2 수정사항

생성일: 2026-05-08

## 한 줄 결론

기존 v4.2 원고의 핵심 결론은 유지됩니다. 다만 Mplus 7 재현분석을 추가하면서 방법, 결과, 보충표, 보충그림, 재현 스크립트를 보강했습니다. 작성 중인 원고에는 아래 항목만 합치면 됩니다.

## 반드시 반영할 수정사항

1. Methods에 Mplus cross-validation 추가
   - 9-item depression LCGA를 Python EM으로 primary 추정했고, R lcmm와 Mplus 7로 독립 재현했다는 문장을 넣습니다.
   - Mplus 모델은 quadratic LCGA, time scores = -3 to 3, class-specific growth means, growth factor variances fixed to zero, residual variance equal across waves/classes입니다.
   - Mplus K=1-5 모두 정상 종료했고 best loglikelihood가 replicated되었습니다.

2. Results 3.2에 Mplus 재현 결과 추가
   - K=4 Mplus BIC = 25,228.08
   - Python EM K=4 BIC = 25,228.07
   - 차이는 반올림 수준입니다.
   - Mplus raw class는 manuscript label과 순서가 달라서 label alignment가 필요합니다.

3. Mplus label alignment 명시
   - Mplus raw 2 -> C1 Stable-Low
   - Mplus raw 4 -> C2 Persistent-High
   - Mplus raw 1 -> C3 Decreasing
   - Mplus raw 3 -> C4 Late-Increasing
   - Label-aligned modal agreement = 99.96%, Cohen's kappa = .9994입니다.

4. AIC, BIC, aBIC 포함 여부
   - 모두 포함했습니다.
   - `02_tables/TableS_Mplus_fit_replication_9item.csv`에 K=1-5 free parameters, logL, AIC, BIC, aBIC, entropy, Mplus-Python 차이가 들어 있습니다.
   - 본문 Table 2는 primary model selection용이고, Mplus 전체 비교는 supplementary table로 두는 것을 권장합니다.

5. Table 6 최신본 유지
   - `02_tables/Table6_cox_HR_9item.csv`가 최신 cluster-bootstrap Cox per-wave HR 표입니다.
   - 예전 average HR/Schoenfeld 중심 Table 6로 돌아가면 안 됩니다.
   - bootstrap p-value floor note도 CSV와 `02_tables/Table6_notes.txt`에 포함했습니다.

6. Figure 보완
   - 기존 R Figure 1-6은 그대로 유지합니다.
   - 새 보충그림 `03_figures/Mplus_validation/FigureS1_Mplus_replication_9item.*`를 추가했습니다.
   - Figure S1은 Python vs Mplus BIC overlap, AIC/BIC/aBIC 차이, K=4 class count alignment를 보여줍니다.

## 원고에 넣을 수 있는 짧은 문장 예시

Methods:

> As a cross-platform validation, the primary 9-item LCGA solution was independently re-estimated in Mplus 7 using the same quadratic LCGA specification and K = 1-5 class enumeration.

Results:

> Mplus closely replicated the Python EM solution. For the selected four-class model, Mplus yielded BIC = 25,228.08 compared with Python EM BIC = 25,228.07, and label-aligned modal classification agreement was 99.96% (Cohen's kappa = .9994).

Supplementary material:

> Supplementary Figure S1 and Supplementary Tables S1-S2 present the Mplus replication, including AIC, BIC, sample-size adjusted BIC, entropy, and class-label alignment.

## 주요 수치

| 항목 | 결과 |
|---|---:|
| 분석표본 | N = 2,452 |
| Primary scale | 9-item depression, suicidal ideation item E04 excluded |
| Mplus K=4 BIC | 25,228.08 |
| Python EM K=4 BIC | 25,228.07 |
| Mplus K=4 AIC | 25,135.20 |
| Mplus K=4 aBIC | 25,177.24 |
| Mplus K=4 entropy | .661 |
| Modal agreement | 99.96% |
| Cohen's kappa | .9994 |

## 작성자가 사용할 파일 체크리스트

- 최신 원고: `01_manuscript/MANUSCRIPT_v4.3_MPLUS.md`
- Mplus fit 보충표: `02_tables/TableS_Mplus_fit_replication_9item.csv`
- Mplus class alignment 보충표: `02_tables/TableS_Mplus_class_alignment_9item.csv`
- Mplus 보충그림: `03_figures/Mplus_validation/FigureS1_Mplus_replication_9item.tiff`
- Mplus 전체 산출물: `04_results_data/mplus/MPLUS_REPLICATION_REPORT.md`
- 재현 스크립트: `06_reproducibility/scripts/33_prepare_mplus_lcga_9item.py`, `34_parse_mplus_lcga_9item.py`, `35_mplus_replication_figure.py`

## 표현상 주의점

- Mplus 결과는 "primary 결과 변경"이 아니라 "cross-platform replication 강화"로 쓰는 편이 좋습니다.
- K=5는 AIC/BIC가 더 낮지만, 최소 class 3.7% 수준과 해석가능성/재현성/기존 reviewer 대응 흐름을 고려해 K=4를 primary로 유지합니다.
- BLRT와 bootstrap p-value는 add-one floor 규칙을 유지해야 합니다. BLRT B=30이면 p=.032, Cox bootstrap B=300이면 floor가 1/301=.00332입니다.
