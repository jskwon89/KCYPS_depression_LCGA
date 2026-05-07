# KCYPS 2018 (m1 cohort, 1차~7차, 2018-2024) — SSCI Analysis Design

## 데이터 요약
- **표본**: 중1 코호트 원패널 청소년 N=2,590 → 2,073 (7차, 80% retention)
- **시점**: 2018(중1) → 2024(고3 직후/대1)
- **자료**: 청소년 본인 응답 (m1Yw1~7)
- **추가**: 보호자 응답 (m1Pw1~7), 형제자매 (m1Sw2~7)

## 핵심 척도 (모두 7차에서 측정, 신뢰도 우수)

| 약자 | 척도 | 문항수 | α (1차/4차/7차) | 출처 |
|------|------|--------|-----------------|------|
| YPSY1 | 삶의 만족 (SWLS) | 5 | .85/.81/.75 | Diener et al. 1985 |
| YPSY2 | 주관적 행복 (SHS) | 4 | .80/.70/.58 | Lyubomirsky & Lepper 1999 |
| YPSY3 | 자아존중감 (Rosenberg SES) | 10 | .87/.83/.77 | Rosenberg 1965 |
| YPSY4A | 주의집중 문제 | 7 | .82/.84/.84 | KCYPS 자체 |
| YPSY4B | 공격성 | 6 | .84/.85/.85 | KCYPS 자체 |
| YPSY4C | 신체증상 | 8 | .87/.86/.87 | KCYPS 자체 |
| YPSY4D | 사회적 위축 | 5 | .87/.88/.88 | KCYPS 자체 |
| YPSY4E | **우울 (자살사고 포함)** | 10 | **.92/.90/.91** | 김광일 등 BDI 변형 |
| YPSY7 | 그릿 (Short Grit) | 8 | .71/.61/.55 | Duckworth & Quinn 2009 |
| YFAM2A~F | Skinner PSCQ-A 양육 (6하위×4) | 24 | .77~.91 | Skinner et al. 2005 |
| YEDU2 | 또래 관계 | 13 | .85 | KCYPS 자체 |
| YEDU3 | 교사 관계 | 14 | .91~.97 | KCYPS 자체 |
| YDLQ1 | 오프라인 비행 (15유형, 흡연/음주/폭력/도박 등) | 15 | n/a (count) | KCYPS |
| YDLQ2 | 사이버 비행 (스토킹/험담/신상털기 등) | 15 | n/a (count) | KCYPS |
| YMDA1C | 스마트폰 의존 | 15 | .85~.88 | 한국정보화진흥원 |

---

## 분석 설계 — 4개 SSCI 후보

### **Paper A (Primary, 우선 진행)** — Latent Class Growth Analysis of Adolescent Depression with Suicide Risk Outcomes (LCGA + GMM)

> **RQ**: 한국 청소년의 우울 궤적은 몇 개의 잠재 계급으로 분류되며, 1차(중1) 시점의 양육·또래·자아 요인이 어느 계급으로 진입할지를 예측하는가? 7차 시점의 자살사고·삶의만족·비행은 계급별로 어떻게 다른가?

- **표본**: N=2,590 (1차) → 7차까지 추적
- **DV (반복측정)**: 우울 (YPSY4E) 1~7차
- **방법**:
  1. **Step 1** — 무조건 잠재계급성장분석(LCGA, 2~5개 계급) → BIC/aBIC/entropy/LMR-LRT/BLRT
  2. **Step 2** — Growth Mixture Model (GMM, class-varying random effects)
  3. **Step 3** — 분류 Manual 3-step (BCH 가중치) → 1차 공변량(부모온정/거부, 또래애착, 자아존중감, 그릿, 성별)이 계급 멤버십에 미치는 영향 (다항 로지스틱)
  4. **Step 4** — 7차 distal outcomes: 자살사고(YPSY4E04 dichotomized), 삶의만족, 오프라인 비행 — 계급별 평균/비율 비교
  5. **Step 5** — 강건성: 가중치(WEIGHTB1) 적용, 결측 FIML, 성별 다집단

- **타겟 저널**: 
  - **Journal of Affective Disorders** (IF≈6.6, Q1) — primary
  - **Development & Psychopathology** (IF≈3.7) — secondary
  - **Journal of Youth and Adolescence** (IF≈3.7) — backup

- **기여**: (a) 한국 동일 코호트의 7년 우울 궤적, (b) early adolescent → emerging adult 전환의 계급 식별, (c) Skinner 6차원 양육이 계급에 미치는 차별적 영향, (d) 자살사고 distal outcome.

---

### **Paper B (Secondary)** — RI-CLPM of Cyberbullying Perpetration and Depression (Within-Person Bidirectional Effects)

> **RQ**: 사이버 비행 가해와 우울은 within-person 수준에서 양방향으로 영향을 주는가?

- **변수**: 사이버 가해 (YDLQ2, count) ↔ 우울 (YPSY4E) — 7파동
- **방법**:
  - **Random-Intercept Cross-Lagged Panel Model** (RI-CLPM; Hamaker, Kuiper & Grasman, 2015)
  - 측정동일성 검정 후 lagged paths의 점진적 단순화
  - 부모 양육·또래 애착을 시변·시불변 공변량
  - **확장**: Multiple-Group RI-CLPM (성별), Latent Growth Curve Mediation
- **타겟**: **Computers in Human Behavior** (IF≈9.0), Cyberpsychology Behavior & Social Networking (IF≈4.2)

---

### **Paper C (Machine Learning)** — Predicting Suicidal Ideation in Late Adolescence from Early-Wave Multi-Domain Features (XGBoost + SHAP)

> **RQ**: 중1(1차) 시점의 가족·학교·심리·미디어 다영역 특성으로 5~7차의 자살사고 발현을 어느 정도 예측할 수 있고, 어떤 요인이 가장 기여하는가?

- **DV**: 자살사고 (YPSY4E04 ≥2 또는 우울점수 상위 10%) 5~7차 중 한 번이라도 발생
- **predictors (≈100개)**: 1차의 30+ 척도 점수, 인구통계, 부모/가구 특성, 시간사용 등
- **방법**: 
  - 5-fold nested CV with hyperparameter tuning
  - XGBoost + LightGBM + Logistic ridge — 비교
  - SHAP global/local importance
  - Decision Curve Analysis (clinical utility)
  - Calibration curves
- **타겟**: **JMIR Mental Health** (IF≈4.8), Journal of Affective Disorders, Journal of Adolescent Health (IF≈5.2)

---

### **Paper D (Network Analysis)** — Dynamic Symptom Network of Internalizing/Externalizing Across Adolescence

> **RQ**: 우울/불안/공격성/위축/주의집중 36개 증상의 네트워크 구조는 청소년기 동안 어떻게 변하는가? 어떤 증상이 bridge symptom인가?

- **변수**: YPSY4 36문항 × 7차
- **방법**:
  - 차수별 GGM (qgraph/bootnet 식)
  - **Time-varying GGM** (mgm 패키지 또는 Python `gimme`)
  - Network Comparison Test (NCT)
  - Bridge centrality (between internalizing/externalizing communities)
- **타겟**: **Journal of Affective Disorders**, Psychological Medicine (IF≈7.0)

---

## 우선순위 결정

이 세션에서 1) 데이터/척도 검증 2) **Paper A (LCGA/GMM)** 3) **Paper B (RI-CLPM)** 4) **Paper C (ML+SHAP)** 까지 파일럿을 모두 시도. 결과 해석에 따라 1순위 결정.
