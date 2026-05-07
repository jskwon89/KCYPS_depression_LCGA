# Paper A — Korean Adolescent Depression Trajectories (KCYPS 2018) — 작성 담당자 인계 패키지

**최종 버전**: Manuscript v4.2 (2026-05-08)
**상태**: 통계 분석 완료, 본문 작성 단계
**1순위 타겟 저널**: *Journal of Affective Disorders* (IF≈6.6, Q1)

---

## 폴더 구조 안내

```
PaperA_LCGA/
├─ README.md                              (이 문서)
├─ 01_manuscript/                         매뉴스크립트 + version history
│  ├─ MANUSCRIPT_v4.2.md                  ⭐ 최신 본문 (작업 시작점)
│  ├─ HISTORY_v3.md                       이전 버전 변경 로그
│  ├─ HISTORY_v4.md                       v4 변경 로그
│  └─ PAPER_PORTFOLIO_RANKING.md          KCYPS 자료에서 가능한 논문 6개 우선순위
├─ 02_tables/                             출판 품질 표 (CSV, APA 형식 가까움)
│  ├─ Table1_sample_characteristics_9item.csv
│  ├─ Table2_model_fit_9item.csv
│  ├─ Table3_class_profiles_9item.csv
│  ├─ Table4_baseline_by_class_9item.csv
│  ├─ Table5_BCH_distal_outcomes_9item.csv
│  ├─ Table6_cox_HR_9item.csv             (cluster-bootstrap 기반)
│  └─ Table6b_cox_HR_per_wave_full_9item.csv (보충 supplement)
├─ 03_figures/                            출판 품질 그림
│  ├─ README_FIGURES.md                   ⭐ figure 사용 가이드 + caption 초안
│  ├─ R_ggplot2/                          ⭐⭐ 투고용 — R/ggplot2 버전 사용 권장
│  │  ├─ Figure1_LCGA_trajectories_R.{png,pdf,tiff}     (600 dpi PNG, vector PDF, 300 dpi TIFF)
│  │  ├─ Figure2_suicidal_ideation_R.{png,pdf,tiff}
│  │  ├─ Figure5_fit_indices_R.{png,pdf,tiff}
│  │  └─ Figure6_time_varying_HR_R.{png,pdf,tiff}
│  └─ (Python matplotlib 버전 — 초안 검토용)
│     Figure1_..._9item.{png,pdf}, Figure2..6 동일
├─ 04_results_data/                       모든 통계 결과 raw CSV
│  ├─ class_assignments_K4_9item.csv      개인별 class 할당 + posterior probs
│  ├─ class_profiles_K4_9item.csv         4개 class trajectory parameters
│  ├─ fit_summary_9item.csv               K=1..5 모델 fit indices
│  ├─ comparison_to_10item.txt            9문항 vs 10문항 sensitivity (κ=0.929)
│  ├─ cox_HR_class_9item.csv              Cox 평균 HR
│  ├─ cox_PH_schoenfeld_test.csv          PH 가정 검정
│  ├─ cox_HR_per_wave_with_CI_9item_BOOT.csv  per-wave HR + bootstrap CI + p
│  ├─ cox_HR_per_wave_9item_BOOT.csv      wide format
│  ├─ cox_interaction_bootstrap_dist_9item.csv  B=300 bootstrap distributions
│  ├─ suic_long.csv, survival_data.csv    7파동 자살사고 종단 자료
│  ├─ suic_any_by_class_wave.csv          파동별 자살사고 비율 (≥2)
│  ├─ suic_strong_by_class_wave.csv       파동별 강 자살사고 비율 (≥3)
│  ├─ bch/                                BCH 3-step distal outcome 보정
│  ├─ blrt/                               Bootstrap LRT (B=30, p=.032)
│  ├─ sensitivity_weighted/               longitudinal weight sensitivity (κ=0.939)
│  └─ lcmm_R/                             R lcmm 9-item 재현 (BIC 일치)
├─ 05_methods_appendix/                   방법 부록
│  ├─ KEY_RESULTS_CHEATSHEET.md           ⭐ 본문 인용용 핵심 수치 모음
│  ├─ METHODS_TECHNICAL_NOTES.md          ⭐ 방법 상세 (LCGA EM, BCH, Cox bootstrap)
│  ├─ KCYPS_DATA_NOTES.md                 ⭐ 자료/척도/문항 정의
│  └─ REFERENCES_TEMPLATE.md              레퍼런스 초안
└─ 06_reproducibility/
   └─ scripts/                            전체 재현 스크립트 (Python + R)
```

---

## 작업 우선순위 (시간 추정)

| 순위 | 항목 | 추정 시간 | 비고 |
|------|------|-----------|------|
| 1 | Discussion 1,500단어 본문 작성 | 3-5일 | outline은 v4.2에 있음 |
| 2 | Reference list (50-70개) JAD style | 2일 | 본문 인용 + 핵심 선행연구 |
| 3 | Abstract 다듬기 (300단어) | 0.5일 | 현재 길이 OK, 표현 다듬기 |
| 4 | IRB / author info / affiliations | 1일 | 기관 IRB 신청 또는 secondary data 면제 확인 |
| 5 | Cover letter 작성 | 0.5일 | JAD 양식 |
| 6 | (옵션) BLRT B=1000 supplement | 24-48h compute | p<.001 보고를 위해 |
| 7 | (옵션) Cluster bootstrap B=1000 | 1-2h compute | floor을 .001로 낮춤 |

**전체 게재 준비까지 추정**: 1-2주

---

## 핵심 메시지 (본문 narrative 골격)

> 한국 청소년 우울 trajectory는 **4개 잠재 계급**으로 분류된다 (LCGA 9문항 기준):
> 
> - **C1 Stable-Low/Resilient (25.4%)**: 항상 낮음, W7 자살사고 8.0%
> - **C2 Persistent-High (10.5%)**: 만성 고위험, W7 자살사고 47.7% (BCH OR 10.6)
> - **C3 Decreasing (29.3%)**: 초기 높았다가 회복, W7 자살사고 30.4% (BCH OR 5.1)
> - **C4 Late-Increasing (34.8%)** ⚠️: W1엔 정상이지만 점진 상승, **W7 자살사고 56.8% (BCH OR 15.2)** — C2와 동등 또는 그 이상
> 
> Cox 시간변화 hazard 분석에서 **C4 위험은 W1 HR 1.43 (n.s.) → W7 HR 7.45 (p=.003)** 로 emerge하며, **W4-5 (고등학교 전환기)에서 C3 (하향) class와 cross-over**.
> 
> **임상적 의의**: 중1 시점 단면 평가는 코호트의 1/3에 달하는 잠재 위험군 (C4)을 발견 못함. 고등학교 전환기 반복 평가가 필수.

---

## 방법 핵심 (논문 신뢰성 방어선)

1. **Criterion-clean 9문항 우울**: E04 자살사고 항목을 trajectory 척도에서 제외 → distal SI outcome과의 contamination 없음. 10문항은 sensitivity로 보존 (κ=0.929).
2. **R lcmm cross-validation**: Python EM K=4 BIC 25,228.07, R lcmm 25,228.08 (소수점 일치). modal n도 ±1명 일치.
3. **BCH 3-step distal**: classification error 보정. Bootstrap CI는 [0,1]로 clipped.
4. **Schoenfeld test로 PH 위반 진단** + Cox **time-varying class × wave 모형**으로 정량화.
5. **ID-level cluster bootstrap** (B=300) per-wave HR CI: lifelines `variance_matrix_`는 Hessian inverse로 cluster-sandwich가 아님 → bootstrap이 적절. p=.003은 add-one floor.
6. **종단 가중치 sensitivity** (WEIGHTB2w7): Hungarian-aligned κ=0.939, class 구조 동일.

---

## ⚠️ 본문 작성 시 주의사항

1. **Class 라벨**: C1=Stable-Low, C2=Persistent-High, C3=Decreasing, C4=Late-Increasing. *5개 모든 표/그림에서 통일됨.*
2. **n vs π̂**: "Class size"는 modal-assignment n (615/245/726/866). mixing proportion π̂ (25.4/10.5/29.3/34.8%)은 BCH distal에 사용.
3. **R lcmm raw class 매핑** (혹시 reviewer 질문 대비): R 1=C1, R 3=C2, R 2=C3, R 4=C4. 표 `04_results_data/lcmm_R/python_vs_R_class_comparison_9item.csv` 참조.
4. **BLRT p**: B=30, observed LRT 2,174.7, p=.032 (1/31 floor). 본문에서 "p<.001"로 쓰지 말 것 — add-one rule 위반.
5. **Cox bootstrap p**: B=300, floor 1/301=.0033 → .003 (반올림). p=.003은 floor 명시 footnote 필수. p≤.003 같은 표현 금지.
6. **숫자 인용 시**: `02_tables/`의 표 그대로 옮기지 말고, 표는 표대로 두고 본문은 narrative 형태로. 모든 수치는 `04_results_data/`에서 검증 가능.
7. **그림 — 투고용은 R 버전 사용 필수**: `03_figures/R_ggplot2/`의 PNG (600 dpi), PDF (vector, 폰트 임베드), TIFF (300 dpi LZW)을 사용해주세요. Python matplotlib 버전은 초안 검토용입니다. JAD/Elsevier는 보통 TIFF (LZW) 또는 vector PDF를 권장. 상세는 `03_figures/README_FIGURES.md`에 figure caption 초안과 함께 정리.

---

## 인계 시점 게재 준비도: ~80%

남은 20%는 본문 작성·레퍼런스·IRB·cover letter. 통계 측면에서 reviewer 가져갈 큰 결함은 처리됨 (criterion contamination, BLRT precision, Cox cluster correction, R cross-validation, label alignment 모두 해결).

질문/이슈 발생 시 `06_reproducibility/scripts/`의 모든 스크립트가 self-contained로 재실행 가능.
