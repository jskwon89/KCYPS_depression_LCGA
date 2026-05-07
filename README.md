# KCYPS 2018 — Korean Adolescent Depression Trajectories Analysis

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)

본 저장소는 **한국아동·청소년패널조사 2018 (KCYPS 2018) 중1 코호트 1~7차 (2018–2024)** 자료를 활용한 청소년 우울 trajectory + 자살사고 distal outcome 분석의 통계·재현 패키지입니다.

---

## 📚 패키지 구성

이 데이터셋에서 도출 가능한 SSCI 후보 논문 6편 중, **Paper A — 우울 trajectory LCGA + 자살사고 distal outcome** 가 최우선 진행 중입니다.

| Paper | DV | 방법 | 상태 |
|-------|-----|------|------|
| **A** | 9-item depression trajectory (W1–W7) | LCGA + BCH 3-step + Cox time-varying + cluster bootstrap | **v4.2 (작성 중)** |
| B | Cyber/offline delinquency | FE panel + GEE onset + ML | 별도 분석 패키지 (외부) |
| B' | Smartphone × Depression | RI-CLPM | pilot only |
| C | Suicidal ideation prediction | XGBoost + SHAP | pilot only |
| D | Symptom network | Graphical Lasso GGM | pilot only |

상세는 [`output/paperA_lcga/RANKING_FINAL.md`](output/paperA_lcga/RANKING_FINAL.md) 참조.

---

## ⭐ Paper A 핵심 결과

- **Sample**: N=2,452 (m1 cohort, ≥3 valid waves of 9-item depression)
- **4-class LCGA solution** (BIC 25,228; entropy 0.66; B=30 BLRT p=.032 — add-one floor):
  - **C1 Stable-Low** (25.4%) — resilient, W7 SI 8.0%
  - **C2 Persistent-High** (10.5%) — chronic, W7 SI 47.7% (BCH OR 10.6 vs C1)
  - **C3 Decreasing** (29.3%) — recovering, W7 SI 30.4% (BCH OR 5.1)
  - **C4 Late-Increasing** (34.8%) ⚠️ — invisible at W1, **W7 SI 56.8%** (BCH OR 15.2)
- **Cox cluster-bootstrap** (B=300): C4 W1 HR 1.43 (n.s., p=.10) → W7 HR 7.45 (p=.003); class-by-time interaction p=.003
- **R lcmm cross-validation**: K=4 BIC 25,228.08 ≡ Python EM 25,228.07 (modal n ±1 individual)

---

## 🗂 폴더 구조

```
청소년패널조사/
├─ README.md                              (이 문서)
├─ .gitignore                             (KCYPS .dta 데이터 제외)
├─ scripts/                               분석 스크립트 (Python + R)
│  ├─ 06_build_long_panel.py              long-form 패널 + 18개 척도 점수
│  ├─ 22_lcga_9item_primary.py            ⭐ 9문항 LCGA primary
│  ├─ 23_bch_cox_9item.py                 BCH 3-step + Cox + Schoenfeld
│  ├─ 25_weighted_sensitivity_9item.py    가중치 sensitivity
│  ├─ 26_lcmm_9item.R                     R lcmm 재현
│  ├─ 27_publication_outputs_9item.py     Tables 1-5 + Figures 1-5
│  ├─ 28_blrt_9item.py                    Bootstrap LRT
│  ├─ 30_cox_cluster_bootstrap.py         ⭐ ID-level cluster bootstrap
│  ├─ 31_lcmm_label_alignment.py          R lcmm ↔ manuscript label
│  └─ 32_rebuild_table6.py                Table 6 / 6b (bootstrap 기반)
│
├─ output/                                전체 분석 결과
│  ├─ paperA_lcga/                        ⭐ Paper A 메인
│  │  ├─ MANUSCRIPT_DRAFT.md              v4.2 매뉴스크립트
│  │  ├─ PROGRESS_REPORT_v{1,2,3,4}.md    버전별 변경 로그
│  │  ├─ RANKING_FINAL.md                 6편 portfolio 우선순위
│  │  └─ primary_9item/                   9문항 primary 결과 (criterion-clean)
│  │     ├─ class_assignments_K4_9item.csv
│  │     ├─ tables/Table1..6.csv
│  │     ├─ figures/Figure1..6.{png,pdf}
│  │     ├─ bch/, blrt/, long_distal/, sensitivity_weighted/, lcmm_R/
│  │     └─ comparison_to_10item.txt
│  ├─ paperB_riclpm/                      Smartphone × Depression RI-CLPM
│  ├─ paperC_ml/                          ML 자살사고 예측
│  ├─ paperD_network/                     증상 GGM
│  ├─ codebook/                           KCYPS 코드북 파싱 결과
│  └─ desc/                               기술통계, ICC, 신뢰도
│
└─ 전달용/                                 작성 담당자 인계 패키지
   └─ PaperA_LCGA/                        ⭐⭐ 큐레이트된 인계 폴더
      ├─ 01_manuscript/                   매뉴스크립트 + history
      ├─ 02_tables/                       publication tables
      ├─ 03_figures/                      publication figures (300dpi)
      ├─ 04_results_data/                 모든 raw 결과 CSV
      ├─ 05_methods_appendix/             ⭐ KEY_RESULTS_CHEATSHEET, METHODS, DATA, REFERENCES
      ├─ 06_reproducibility/              스크립트 + 실행 가이드
      └─ README.md                        인계 안내
```

---

## 🚀 재현 방법

### 환경
```bash
# Python
pip install numpy pandas scipy statsmodels scikit-learn lifelines pyreadstat matplotlib

# R (4.5+)
install.packages(c("lcmm", "haven", "dplyr", "tidyr"))
```

### 데이터 준비
KCYPS 2018 m1 코호트 .dta 파일 (W1–W7)이 필요합니다. 한국청소년정책연구원에서 신청 후 다운로드:
- https://www.nypi.re.kr/archive/contents/siteMain.do
- 본 저장소에는 자료 미포함 (NYPI 라이센스).

`KCYPS2018m1Yw{1..7}.dta`, `KCYPS2018m1Pw{1..7}.dta` 파일을 저장소 루트에 배치.

### 실행
```bash
# 1. 패널 구축 (5분)
python scripts/06_build_long_panel.py

# 2. Paper A 본 분석 (약 50분)
python scripts/22_lcga_9item_primary.py
python scripts/23_bch_cox_9item.py
python scripts/25_weighted_sensitivity_9item.py
python scripts/27_publication_outputs_9item.py
python scripts/28_blrt_9item.py
python scripts/30_cox_cluster_bootstrap.py
python scripts/31_lcmm_label_alignment.py
python scripts/32_rebuild_table6.py

# 3. R lcmm cross-validation (약 30분)
Rscript scripts/26_lcmm_9item.R
```

상세는 `전달용/PaperA_LCGA/06_reproducibility/README.md` 참조.

---

## 📝 인용

> Korean Children and Youth Panel Survey 2018 (KCYPS 2018), middle-school first-grade (m1) cohort, Waves 1–7 (2018–2024). National Youth Policy Institute, Republic of Korea.

논문이 게재되면 본 README에 citation 업데이트 예정.

---

## 📄 라이선스

- 분석 스크립트 / 결과 산출물: MIT License
- KCYPS 2018 원자료: NYPI public-use license (본 저장소에 미포함)

---

## 👤 연구진

- 분석: Jeong-Sung Kwon (jskwon89)
- 자료 출처: 한국청소년정책연구원 (NYPI)
- 통계 자문: TBD
- 작성 담당: TBD

---

## 📅 버전 이력

| 버전 | 날짜 | 주요 변경 |
|------|------|-----------|
| v4.2 | 2026-05-08 | BLRT add-one rule strict, Cox bootstrap p=.003 (floor noted), Table 6 reproducibility |
| v4.1 | 2026-05-08 | Table 6 rebuilt from cluster bootstrap, Table 6b for supplement |
| v4 | 2026-05-08 | ID-level cluster bootstrap (B=300), R lcmm ↔ manuscript label mapping |
| v3 | 2026-05-08 | 9-item BLRT, R lcmm 9-item complete, modal-n vs π̂ 분리 |
| v2 | 2026-05-08 | 9-item primary (criterion contamination fix), Hungarian-aligned weighted κ |
| v1 | 2026-05-07 | 10-item composite (now sensitivity), modal-class logistic distal |
