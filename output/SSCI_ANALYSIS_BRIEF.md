# KCYPS 2018 (m1 코호트) — SSCI 게재용 분석 Brief

> **참고**: 본 폴더에는 기존 세션에서 작성된 별도 분석 결과 `output/kcyps_ssci_report.md`가 있습니다 (사이버·현실 비행을 종속변수로 한 fixed-effects panel + GEE onset + K-means trajectory + ML). 본 brief의 분석은 **우울/자살사고/내재화** 측면을 다루므로 상호보완적이며, 두 분석을 합치면 하나의 데이터셋에서 **2~3개의 독립 SSCI 논문**이 가능합니다.

**데이터**: KCYPS 2018 중1 코호트 1~7차 (2018~2024, age 13→19)
**표본**: N=2,590 (W1) → 2,073 (W7), 7-wave completion 70.1%
**핵심 척도 신뢰도**: 우울 α=.90 (7차 평균), 자아존중감 α=.83, 부모온정 α=.88, 또래애착 α=.85

작업폴더: `D:/2026/SCI/청소년패널조사/`

---

## ⭐ Paper A — Latent Class Growth Analysis of Adolescent Depression with Suicidal Ideation as Distal Outcome

> **권장 1순위 — Journal of Affective Disorders (IF≈6.6)**

### 핵심 결과
**LCGA K=4 모델 채택** (BIC: 1→27,687 / 2→25,596 / 3→25,008 / **4→24,648** / 5→24,538; entropy=0.66)
- K=5는 한 클래스가 3.1%로 너무 작음 → K=4 채택 (Nagin 권장 5%+ 룰)

| 클래스 | 비율 | W1 우울 | W7 우울 | 패턴 | W7 자살사고 (≥2) | OR vs C2 |
|--------|------|---------|---------|------|---------------------|----------|
| **C1 Late-Increasing** | 36.3% | 1.49 | 1.91 | 점진 상승 | **49.6%** | **6.34** |
| **C2 Stable Low/Resilient** | 25.6% | 1.43 | 1.35 | 항상 낮음 | 13.5% | 1.00 (ref) |
| **C3 Decreasing** | 27.2% | 2.27 | 1.65 | 회복 | 31.7% | 2.94 |
| **C4 Persistent High/Inv-U** | 10.9% | 2.36 | 2.08 | 만성·역U형 | **50.4%** | **6.37** |

### 임상적·이론적 핵심 메시지
1. **"보이지 않는 위험군"의 발견**: C1 (36.3%) 청소년은 W1엔 정상 우울 수준이지만, 6년 추적 결과 W7 자살사고 위험이 만성 고위험 C4와 **동등** (OR 6.34 vs 6.37). 중1 시점의 단면 평가는 1/3에 달하는 잠재 위험군을 놓침.
2. **여성 과대 표상**: C4 (만성 고위험) 60.9% 여성, C3 (회복) 52.8% 여성 (χ²=53.8, p<10⁻¹¹)
3. **Skinner 6차원 양육이 모두 차별화**: C2(저위험) → C4(고위험) 비교에서 부모 온정 -.42, 거부 +.38, 혼돈 +.42 (모두 표준화 평균차이)
4. **자아존중감 W1**: C2=3.23 vs C4=2.59 (대형 효과)

### 게재 잠재 저널 (현실적)
- **Journal of Affective Disorders** (IF 6.6, Q1) — primary target. 한국 청소년 longitudinal trajectory + suicide 강력 매치
- **Development & Psychopathology** (IF 3.7, Q1) — secondary
- **Journal of Youth and Adolescence** (IF 3.7, Q1) — backup
- **Suicide and Life-Threatening Behavior** (IF 3.5)

### 산출 파일 위치
```
output/paperA_lcga/
  ├─ class_assignments_K4.csv          # 개인별 계급 할당 + posterior
  ├─ class_profiles_K4.csv             # 4개 계급의 성장 모수 (절편/선형/이차)
  ├─ baseline_by_class.csv             # W1 평균 by 계급 (anova_baseline_by_class.csv 포함)
  ├─ logit_suic_OR.csv                 # 자살사고 OR (C2 ref)
  ├─ suicidal_ideation_by_class.csv
  ├─ fit_summary.csv                   # K=1~5 모델 비교
  └─ lcga_K4_observed_means.png        # 4개 궤적 시각화 (영문 라벨)
```

### 추후 강화 작업
- [ ] **Mplus 본분석**: 본 파이썬 LCGA를 Mplus type=mixture (FIML+MAR + BCH 3-step) 으로 재현. 동일성 검정 + LMR-LRT/BLRT
- [ ] **Bootstrap LRT**: K=3 vs K=4 vs K=5 BLRT
- [ ] **Sensitivity**: 종단 가중치(WEIGHTB1) 적용, 결측 패턴 강건성 (Pattern-mixture)
- [ ] **구체적 기여**: COVID-19 영향 wave (W3=2020, W4=2021) 내 trajectory shift 검토

---

## ⭐ Paper B — RI-CLPM: Smartphone Addiction ↔ Depression (Within-Person Bidirectional Effects)

> **권장 2순위 — Computers in Human Behavior (IF≈9.0) 또는 Cyberpsychology (IF≈4.2)**

### 핵심 결과 (constrained, standardized RI-CLPM, N=2,286)
| 경로 | β_std (constrained) | p | β 범위 (free) |
|------|---------------------|---|----------------|
| Smartphone(t-1) → Smartphone(t) (a) | .220–.268 | <.001 | 0.14–0.32 |
| Depression(t-1) → Depression(t) (c) | .163–.225 | <.001 | 0.16–0.24 |
| Smartphone(t-1) → **Depression(t)** (d) | **.057–.071** | <.001 | -0.04 ~ **+0.13** |
| Depression(t-1) → Smartphone(t) (b) | .028–.032 | .004 | -0.03 ~ +0.09 |

- **X→Y가 Y→X의 약 2배** (d>b)
- **고3→대학 전이 시기(W6→W7)**가 가장 강한 X→Y 효과 (β=.10)
- Random Intercept covariance: 무의미 (-.001) → 두 변수의 trait-level 안정 차이는 독립

### 이론적 의의
- "Displacement hypothesis" (스마트폰이 정신건강 침해)에 대한 **within-person 인과 추론** 강화
- Korean adolescents의 emerging adulthood 진입 단계에서 디지털 사용→우울 효과가 강화됨
- 정책 함의: 고등학교 졸업 / 대학 입학 전후 시기 디지털 웰빙 개입 표적화

### 산출 파일
```
output/paperB_riclpm/
  ├─ riclpm_sm_dep_constrained_est.csv  # constrained model
  ├─ riclpm_sm_dep_free_est.csv         # 자유 추정 (시간이질성 검토)
  └─ riclpm_sm_dep_data.csv             # 표준화 데이터
```

### 보강 권장
- [ ] **측정 동일성**: SM(스마트폰 의존) 척도의 longitudinal MI 검정 후 latent factor RI-CLPM
- [ ] **다집단 분석**: 성별
- [ ] **공변량 추가**: parent_warmth/peer_attach 시변
- [ ] **lavaan/Mplus 재추정**: 적합도 (CFI/TLI/RMSEA/SRMR) — semopy에서는 calc_stats 오류

---

## Paper C — ML/SHAP Prediction of Late-Adolescent Risk

> **권장 3순위 — JMIR Mental Health (IF≈4.8) 또는 Journal of Adolescent Health**

### 5-fold CV AUC (W1 baseline → W7 outcomes)
| 결과변수 | LogL1 | RF | GB | LightGBM | XGBoost | Best |
|----------|-------|------|------|----------|---------|------|
| 자살사고(W7≥2) | .555 | **.592** | .564 | .583 | .577 | RF |
| 고우울(W7 ≥75pctile) | .557 | .557 | .567 | .550 | **.571** | XGBoost |
| **고위험궤적(C1+C4)** | .618 | **.652** | .642 | .634 | .644 | RF |

### 임상 메시지
- 6년 후 개별 자살사고는 W1 변수만으로 AUC=.59 (literature 0.55-0.70 범위 내)
- **궤적 멤버십 예측은 AUC=.65로 가장 높음** — 즉, "위험 vs 회복" 이분 분류는 1/3 더 정확
- 상위 SHAP 변수 (high depression W7): 지역(서울/대도시), 교사관계, 행복감, 부모강요, 스마트폰의존

### 산출 파일
```
output/paperC_ml/
  ├─ summary.csv               # 모델별 AUC/AP/Brier
  ├─ cv_results.csv            # fold-level
  ├─ oof_*.csv                 # Out-of-fold predictions for calibration
  ├─ shap_importance_*.csv     # SHAP 글로벌 중요도
  └─ shap_top15_*.png
```

### 보강 권장
- [ ] **Item-level features** 추가 (현재 scale-level만): YPSY4 36문항, YDLQ 30문항 → SHAP-아이템 단위
- [ ] **Causal Forest**: 양육 효과의 이질성 (HTE estimation)
- [ ] **Decision Curve Analysis**: 임상 효용성 평가
- [ ] **Calibration plots & DCA**

---

## Paper D — Dynamic Symptom Network Analysis (Internalizing/Externalizing 36 Items)

> **권장 — Journal of Affective Disorders 또는 Psychological Medicine (IF≈7.0)**

### 핵심 결과
- 36개 증상 GGM (graphical lasso) at W1, W4, W7 — 모두 density ~0.6
- **W7에서 가장 중심적 증상**: E03 "걱정이 많다", C08 "발열감", D02 "부끄럼", C06 "피곤"
- **가장 강한 bridge symptom**: B06 "이유 없이 운다" (외현화 분류이지만 우울 가교) bridge=0.65
- **연령에 따른 변화 (W1→W7)**:
  - **자살사고(E04) bridge strength +0.16** — emerging adulthood에서 자살사고가 외현화 증상과 더 강하게 연결
  - 주의집중 문제(A01, A07)의 strength + bridge 모두 증가
  - 신체증상(C08)의 중심성 증가

### 의의
- 한국 청소년에서 internalizing/externalizing의 communities 구조가 13 → 19세 동안 어떻게 변형되는가에 대한 첫 종단 evidence
- 자살사고가 단순 우울 군집 내가 아닌 **외현화 증상과의 cross-domain bridging**으로 진화
- 임상 함의: 후기 청소년의 자살 risk 평가 시 우울 점수만이 아닌 attention/aggression 패턴까지 함께 고려

### 산출 파일
```
output/paperD_network/
  ├─ centrality_w1.csv / w4.csv / w7.csv
  ├─ centrality_change_w1_to_w7.csv
  ├─ pcor_w1.npy / w4.npy / w7.npy   (36x36 partial correlation matrices)
  └─ partial_corr_heatmaps.png
```

### 보강 권장
- [ ] **Network Comparison Test (NCT)**: W1 vs W7 통계적 차이 검정
- [ ] **Time-varying GGM (mgm package, R)**: 더 부드러운 시간 변화
- [ ] **Bootstrap stability** (case-drop, edge-stability) — bootnet
- [ ] **R/qgraph 시각화**: 더 publication-ready

---

## 종합 권장 — 진행 우선순위

### 🥇 1순위: Paper A (LCGA + 자살사고 outcome)
**근거**: 
- 데이터 fit 가장 좋음 (BIC drop 명확, 4개 의미 있는 클래스)
- "Late-Increasing 36.3%" 의 발견은 임상적·정책적으로 강력한 메시지
- Korean adolescent + 7-wave + suicide ideation = 매우 publishable 조합
- IF 6+ 저널 게재 현실 가능

**다음 단계**:
1. Mplus로 LCGA + GMM 재추정 (BCH 3-step 정식 적용)
2. wave 가중치 적용
3. 가족 SES, 거주지역 보강 공변량 (W1 부모자료 P 파일에서)
4. 본문 서론·문헌 검토 시작 (한국 청소년 우울 trajectory 선행 KCYPS 4차 패널 논문 비교)

### 🥈 2순위: Paper B (RI-CLPM 스마트폰↔우울)
**근거**: 
- 결과 명확 (β X→Y > Y→X, p<.001)
- "스마트폰 → 우울" 인과 방향에 대한 within-person 증거 → top journal 매력적
- 구현 잘 되어 있음

**다음 단계**:
1. semopy의 fit indices 문제 → R lavaan으로 재실행 (CFI/TLI/RMSEA 보고 필수)
2. 측정동일성 검정
3. 다집단 분석 (성별)
4. 공변량 (양육, 또래) 통합 모형

### 🥉 3순위: Paper D (Network Analysis)
**근거**: 매우 신선한 각도. 자살사고-외현화 bridge 발견은 publishable.

### 4순위: Paper C (ML/SHAP)
**근거**: AUC 0.59-0.65로 modest. paper로는 가능하나 임상 실용성 강조 필요.

---

## 데이터 활용 가능성 추가 메모

### 분석 안 한 영역 (추후 탐색)
- **YPSY7 그릿** (8 items, 7 waves) — 학업 적응 매개 변수로 활용
- **YEDU3 교사 관계** (14 items) — 학교 폭력 가해/피해 보호 요인
- **YACT3 fan/celebrity 행동** — 대중문화 몰입과 정체성 발달
- **YFAM1** 부모와 시간 (32 items, weekday/weekend × 부모유형) — fatherhood 변수로 변환 가능
- **YFUR2** 진로 관련 부모 대화 (9 items) — 진로 발달 trajectory
- **부모 보고 (m1Pw1~7)** 미사용 — multi-informant 분석 가능 (parent-child agreement, cross-informant validation)

### 추가 분석 후보 (5번째 paper)
- **부모-자녀 일치도 (multi-informant)**: YFAM2 자녀 보고 vs 보호자 보고 (P file의 동일 척도) → ICC, parent-child concordance trajectory
- **Multivariate GMM**: 우울+자아존중감 parallel process trajectories
- **Dynamic SEM (DSEM)**: brms나 Stan으로 within-person 무선 효과
- **Causal Forest** (Athey & Wager): 양육 → 우울 변화의 HTE
- **Continuous-time SEM** (ctsem): 측정 간격 고려 모형

---

## 산출물 구조

```
D:/2026/SCI/청소년패널조사/
├─ KCYPS2018m1Yw1.dta ~ Yw7.dta   # 원본 (청소년)
├─ KCYPS2018m1Pw1.dta ~ Pw7.dta   # 보호자 (미사용)
├─ KCYPS2018m1Sw2.dta ~ Sw7.dta   # 형제자매 (미사용)
├─ scripts/                        # 모든 분석 스크립트
│   ├─ 01~05_explore_codebook.py
│   ├─ 06_build_long_panel.py     # 장기 패널 빌드 + 척도 점수
│   ├─ 07_pilot_descriptive.py    # 결측·기술통계
│   ├─ 08_icc_check.py            # ICC 수동
│   ├─ 09_lcga_depression.py      # LCGA 본분석 ⭐
│   ├─ 12_lcga_distal_outcomes.py # distal 결과
│   ├─ 10_riclpm_cyber_depression.py
│   ├─ 11_riclpm_smartphone_depression.py # RI-CLPM ⭐
│   ├─ 13_ml_predict_suicidal.py  # ML+SHAP
│   └─ 14_network_symptoms.py     # 증상 네트워크
├─ output/
│   ├─ codebook/                  # 코드북 파싱 (변수 dictionary)
│   ├─ desc/                      # long_panel.parquet, scale_alpha, ICC
│   ├─ paperA_lcga/               # ⭐ LCGA 결과
│   ├─ paperB_riclpm/             # ⭐ RI-CLPM 결과
│   ├─ paperC_ml/                 # ML+SHAP
│   ├─ paperD_network/            # GGM
│   ├─ ANALYSIS_DESIGN.md         # 분석 설계
│   └─ SSCI_ANALYSIS_BRIEF.md     # 본 문서
```

---

## 결론

본 KCYPS 2018 m1 코호트 7-wave 자료는 **현재 단계에서 즉시 SSCI 1편 (Paper A)을, 보강 후 추가 2-3편을 산출 가능**한 양·질을 갖추고 있습니다. 가장 명확한 1차 contribution은 **"Korean adolescent depression trajectories 13→19세, 4개 잠재 계급, late-increasing 36% 가 자살사고 위험 고위험과 동등"** 이라는 메시지입니다.

다음 회기에 Mplus 재현 + 본문 작성 시작 권장.
