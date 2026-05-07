# KCYPS 2018 m1 코호트 — 6개 분석 통합 우선순위 Ranking

## 검증 단계에서 새로 발견한 사실 (먼저)

### ⚠️ 기존 사이버비행 분석의 W7 현실비행 prevalence 급등(58.5%)은 **construct contamination**

W1 → W7 **항목별** 현실비행 응답률(>=2):

| 항목 | W1 | W7 | delta |
|------|-----|-----|-------|
| **술** | **4.3%** | **52.4%** | **+48.1%p** |
| **담배** | 1.3% | 13.8% | +12.5%p |
| 무단결석 | 1.6% | 6.5% | +4.9%p |
| 가출 | 1.4% | 3.1% | +1.7%p |
| 놀림 | 8.9% | 1.8% | -7.0%p |
| 폭행 | 2.5% | 0.4% | -2.1%p |
| 협박 | 2.0% | 0.3% | -1.7%p |
| (12개 폭력·재산 항목 모두 감소)|  |  | |

**핵심**: W7=2024 시점 응답자 평균 19~20세 → **음주가 합법**, 흡연도 합법. "현실비행" 합성지수의 W7 폭등은 **불법 행동의 증가가 아니라 합법 행동의 정상적 채택**.

**파장**: 기존 리포트가 사이버비행 trajectory 모형에서 lag_offline_freq를 통제변수로 쓰는데, W7→? 전이에서 "직전 현실비행"의 의미가 W1→W2 전이에서와 **다른 잣대**가 됨. K-means 군집 결과의 W7 점은 이 contamination을 이미 반영.

다행히 **사이버비행 자체는 문항 수준에서 안정적**(욕설 31%→16%, 험담 12%→5% 자연 감소), 사이버비행을 결과로 보는 모형은 fixable.

---

## 6개 분석 통합 비교표

| # | 분석 | DV | 핵심 method | 가장 중요한 효과 | 표본 | publishable AS-IS? | 보완 시간 |
|---|------|-----|-------------|------------------|------|---------------------|-----------|
| **A** | **우울 LCGA + 자살사고 distal** | 우울 7파동 | LCGA K=4 EM | C1 36% **OR=6.34** vs C2 (자살사고) | N=2,452 | ⚠ Mplus 재현 필요 | **3-4주** |
| **B** | RI-CLPM 스마트폰↔우울 | 우울/스마트폰 | RI-CLPM 7-wave | X→Y β=.06 (p<.001) > Y→X β=.03 | N=2,286 | ⚠ lavaan 재현 필요 | 2-3주 |
| **C** | ML 자살사고/궤적 예측 | 위험궤적(C1+C4) | RF/XGBoost+SHAP | AUC=0.65 (궤적), 0.59 (자살사고) | N=2,406 | ✓ AUC modest | 1-2주 |
| **D** | 증상 GGM 네트워크 | 36 증상 × 3시점 | Graphical Lasso | 자살사고 bridge +0.16 (W1→W7) | N=2,073-2,590 | ⚠ NCT/qgraph 필요 | 2-3주 |
| **E (기존)** | 사이버비행 cascade | 사이버비행 7파동 | FE+GEE+ K-means+ML | onset OR: SM=1.11, agg=1.12, parent_reject=1.09 | N=2,428 | ⚠ 현실비행 contamination 처리 필요 | 2-3주 |
| **F (기존)** | 현실비행 cascade | 현실비행 | FE | lag_teacher β=.04, lag_LS β=-.04 | - | ❌ W7 음주 문제 | 4주+ |

---

## 우선순위 Ranking (유일한 1위 → 6위)

### 🥇 1위 — Paper A (우울 LCGA + 자살사고 distal)

**왜 1위인가 (3가지 이유)**:
1. **메시지 강도가 가장 높음**: "초기 정상이지만 점진 상승하는 36%가 만성 고위험군과 동등한 자살사고 OR을 가진다" — clinical screening 프로토콜에 직결 (한국 청소년 자살은 사망원인 1위, 정책 적합도 높음)
2. **모델 fit 깨끗함**: BIC 단조감소 + 4계급 비율 36/26/27/11 (균형), entropy 0.66 (수용 범위), 모든 baseline 변수가 ANOVA로 차별화 (p<.001)
3. **타겟 저널 적합도**: Journal of Affective Disorders (IF 6.6), Development & Psychopathology, JYA 모두 한국 longitudinal trajectory + suicide 출판 전례 다수

**현재 상태에서 부족한 것**:
- Mplus LCGA로 재현 (BCH 3-step + LMR-LRT/BLRT) — 1주
- 종단 가중치 적용 + sensitivity — 1주
- 본문 작성 + 선행연구 비교 (특히 KCYPS-1 4차 패널 우울 trajectory 논문들) — 2주

**현실 평가**: 본 분석을 Mplus + lavaan으로 재현 후 SSCI Q1(IF 5+) 게재 확률 **70-80%**.

---

### 🥈 2위 — Paper E (기존 사이버비행 cascade) — 단, 현실비행 처리 후

**왜 2위인가**:
1. 가장 distinct한 DV (사이버비행)로 1위와 **DV 중복 없음** → portfolio strategy에서 최적
2. **GEE onset 모델**은 이론적으로 정밀 (위험집합을 직전 비경험자로 제한하는 risk-set design)
3. **시간차 holdout ML**은 다른 longitudinal ML 논문 대비 차별점 (1-5차 학습 → 6→7차 검증)
4. AUC=0.63은 사회과학 ML 기준 **수용 범위** (Walsh et al. suicide prediction 0.55-0.70 비교)
5. CHB·Aggressive Behavior·Cyberpsychology 출판 전례 풍부

**현재 상태에서 부족한 것 (필수)**:
- **Offline 합성지수에서 음주·흡연 분리** — 본 분석 방향성은 그대로지만 lag_offline_freq를 (a) 음주·흡연 빼고 산출, (b) 19+ 합법 행동 대 미성년 전용 행동 분리
- **K-means 군집(91.8/6.9/1.3 비율)을 LCGA로 교체** — 1.3%(29명) 군집은 reviewer challenge 받기 쉬움. 본인이 직접 LCGA 다시 돌려서 의미 있는 4계급으로 재정의 권장
- FE 동적 모형의 Nickell bias 처리 (Arellano-Bond GMM 대안)

**현실 평가**: 위 처리 후 SSCI Q1(IF 4-6) 게재 확률 **60-70%**. CHB 매력도 매우 높음.

---

### 🥉 3위 — Paper B (RI-CLPM 스마트폰↔우울)

**왜 3위인가**:
- **결과 자체는 깨끗함**: X→Y β=.06 (p<.001) > Y→X β=.03 (p=.004), 7-wave 전 구간
- 그러나 **DV가 1위 (우울)와 겹침** → portfolio에서는 1위에 흡수하거나 mediator로 통합하는 게 살라미 슬라이싱 방지 유리
- 단독 논문으로 가면 1위 LCGA 게재 후 12개월 이상 시차 두고 진행 권장

**살라미 회피 전략**:
- Paper A의 distal outcomes 옆에 "스마트폰 의존 → 우울 → 자살사고" 매개 분석 통합
- 또는 별도 논문으로 출판하되 1위와 다른 저널 (Paper A를 JAD에, Paper B를 CHB에)

**현재 상태에서 부족한 것**:
- semopy의 fit indices 미산출 → R lavaan 재실행 (CFI/TLI/RMSEA/SRMR 보고 필수)
- 측정동일성 검정
- 다집단 분석 (성별)

**현실 평가**: Paper A 출간 후 별도 paper로 SSCI Q1-Q2 게재 확률 60%.

---

### 4위 — Paper D (증상 GGM 네트워크)

**왜 4위인가**:
- **각도가 매우 신선함**: 자살사고(E04)의 bridge strength가 W1 0.085 → W7 0.246 (+189%) — emerging adulthood에서 자살사고가 외현화 증상과 cross-domain 연결 강화. 이런 종단 발견은 별로 없음.
- 그러나 **방법적 완성도 부족**: NCT(network comparison test), bootstrap stability, qgraph 시각화 모두 미실행
- Python sklearn graphical lasso는 publication-ready로는 부족 → R bootnet/qgraph 재실행 필요

**현실 평가**: 본격 보강 시간 (2-3주) 들이면 Psychological Medicine·JAD에 도전 가능. 시간 대비 효과는 1-2위보다 낮음.

---

### 5위 — Paper C (ML/SHAP 예측)

**왜 5위인가**:
- AUC 0.59-0.65: **방어 가능한 수준**이지만 압도적이지 않음
- "위험 궤적(C1+C4) 예측" AUC=0.65는 최고치
- **단독 논문으로는 약함**, Paper A의 supplementary로 통합 권장

**현실 평가**: 단독 게재는 어렵고, Paper A 안의 한 섹션 또는 supplementary material로 활용이 가장 현실적.

---

### 6위 — 기존 현실비행 cascade

**왜 6위인가**:
- 단독 논문으로 가기엔 W7 음주 contamination 너무 큼
- "Korean adolescent의 현실비행 변화" 자체로는 메시지가 약함 (사이버비행만큼 hot topic 아님)
- Paper E(사이버비행)의 통제변수로만 활용 권장

---

## 권장 Portfolio (현실적 시나리오)

### 시나리오 A: 단일 1편 게재 우선 (안전)
**Paper A만 진행 → JAD 또는 JYA 투고**
- 3-4주 보강 후 투고
- 게재 가능성 70-80%
- 1편이지만 강력한 message

### 시나리오 B: 2편 portfolio (적극)
**Paper A (우울→자살사고) + Paper E (사이버비행 cascade)**
- DV가 다르므로 살라미 슬라이싱 우려 없음
- Paper A를 JAD에 1차 투고 → 동시에 Paper E 보강
- Paper E를 CHB나 Aggressive Behavior에 4-6주 후 투고
- 두 논문 모두 게재 확률 50-60% (둘 다 게재 시 큰 성과)

### 시나리오 C: 3편 portfolio (도전적)
**Paper A + Paper E + Paper B (또는 D)**
- 12-18개월 timeline
- 데이터셋 1개로 3편은 각 저널 reviewer가 의식 가능 → ‘본 논문은 KCYPS 자료를 활용한 ___ 시리즈의 일부’ 명시 필요
- Paper A 출간 후 Paper B/D 진행 (선후행 명확하게)

---

## 결론 — 한 문장 요약
> **1순위는 Paper A (우울 LCGA + 자살사고 distal)**, 2순위는 기존 사이버비행 cascade(현실비행 contamination 처리 후), 3순위는 Paper B(RI-CLPM)이며, **Paper A + 사이버비행 = 2편 portfolio**가 가장 현실적 게재 전략입니다.
