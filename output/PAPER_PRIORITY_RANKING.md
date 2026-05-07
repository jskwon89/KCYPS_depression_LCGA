# KCYPS 2018 중1 코호트 SSCI 논문 후보 우선순위

작성일: 2026-05-07

## 결론

최우선 투고 트랙은 **Paper A: 우울 궤적 LCGA + 자살사고 distal outcome**이다. 현재 결과의 메시지, 표본 적합성, 저널 적합도, 임상적 함의가 가장 강하다. 다만 Python EM 기반 LCGA 결과이므로, 원고 착수 전 **Mplus 또는 R tidyLPA/lcmm/flexmix 기반 재현과 BCH 3-step**을 반드시 거쳐야 한다.

두 번째 트랙은 **Paper B: 스마트폰 의존과 우울의 RI-CLPM**이다. CHB/Cyberpsychology 계열의 관심 주제이고 within-person 양방향 효과라는 기여가 분명하다. 단, 현재 semopy 산출물은 적합도 지표와 측정동일성 검정이 부족하므로 **R lavaan/Mplus 재추정**이 통과 조건이다.

세 번째는 **Paper D: 증상 네트워크**다. 참신성은 높지만 현재는 기술적 네트워크 변화 분석에 가까워, Network Comparison Test와 bootstrap stability 없이는 고위 SSCI에서 방어가 약하다.

기존 사이버비행 분석은 단독 논문으로는 4순위다. GEE 진입모형의 결과는 명확하지만 효과크기가 중간 이하이고, K-means 궤적은 투고용 잠재궤적모형으로 대체해야 한다. Paper B나 Paper A의 보조/확장 논문으로 묶는 것이 좋다.

Paper C ML/SHAP은 5순위다. AUC가 .59-.65라 “예측 논문”으로는 임상적 설득력이 약하다. 다른 논문의 보조 분석, 또는 calibration/decision curve/item-level SHAP까지 보강한 뒤 별도 투고를 고려한다.

## 종합 순위

| 순위 | 후보 | SSCI 가능성 | 현재 결과 강도 | 방법론 준비도 | 핵심 리스크 | 추천 타겟 |
|---:|---|---|---|---|---|---|
| 1 | Paper A. Depression trajectories LCGA + suicide outcome | 매우 높음 | 매우 강함 | 중간 | Python LCGA 파일럿, Mplus/BCH 필요 | Journal of Affective Disorders |
| 2 | Paper B. Smartphone dependence ↔ depression RI-CLPM | 높음 | 강함 | 중간 | fit indices, 측정동일성, lavaan 재현 필요 | Computers in Human Behavior, Cyberpsychology |
| 3 | Paper D. Symptom network across adolescence | 중간-높음 | 신선함 | 낮음-중간 | NCT/bootstrap stability 없음 | JAD, Psychological Medicine는 보강 후 |
| 4 | Cyber delinquency onset/developmental cascade | 중간 | 중간 | 중간 | 효과크기 작음, 궤적은 K-means 예비 | Youth & Society, Children and Youth Services Review |
| 5 | Paper C. ML/SHAP late-adolescent risk prediction | 중간 | 약-중간 | 중간 | AUC .59-.65, 임상 효용 약함 | JMIR Mental Health 보강 후 |

## 후보별 판단

### 1순위: Paper A

현재 핵심 결과:

- K=4 우울 궤적이 해석 가능하다.
- Late-Increasing 36.3%, Stable-Low 25.6%, Decreasing 27.2%, Persistent-High 10.9%.
- W7 자살사고율은 Stable-Low 13.5% 대비 Late-Increasing 49.6%, Persistent-High 50.4%.
- OR은 Late-Increasing 6.34, Persistent-High 6.37로 거의 동일하다.
- "초기 단면평가로는 포착되지 않는 late-increasing 위험군"이라는 메시지가 매우 강하다.

왜 1순위인가:

- 우울 궤적과 자살사고는 JAD/Development & Psychopathology/Journal of Youth and Adolescence에 명확히 맞는다.
- 7-wave, 13세에서 19세 전환, 한국 청소년이라는 데이터 기여가 선명하다.
- C1 late-increasing 계급이 전체 36.3%로 작지 않아 정책적 함의가 크다.
- Skinner 양육, 자아존중감, 또래/교사, 성별까지 이론적 설명축이 충분하다.

반드시 보완할 점:

- Mplus `TYPE=MIXTURE`로 K=2-5 LCGA/GMM 재추정.
- BLRT, LMR-LRT, entropy, average posterior probability, OCC 보고.
- K=5가 BIC는 더 낮으므로, K=4 선택 사유를 min class 3.1%, 해석가능성, parsimony로 엄격히 정리.
- 자살사고 문항이 우울척도 YPSY4E의 일부라 distal outcome과 trajectory indicator가 겹친다. 투고용에서는 우울척도에서 E04를 제외한 9문항 우울 궤적으로 재분석하고, E04를 distal outcome으로 두는 민감도 분석이 필요하다. 이 부분이 가장 중요한 심사 리스크다.
- 종단 가중치와 결측 패턴 검토.

권장 원고 제목:

**Hidden Risk in Plain Sight: Depression Trajectories from Early Adolescence to Emerging Adulthood and Suicidal Ideation in a Seven-Wave Korean Cohort**

### 2순위: Paper B

현재 핵심 결과:

- RI-CLPM constrained model에서 스마트폰 의존 → 우울의 표준화 효과가 약 .057-.071.
- 우울 → 스마트폰 의존은 약 .028-.032.
- 양방향이지만 스마트폰 의존 → 우울이 약 2배 강하다.
- free model에서는 일부 구간, 특히 W6→W7 전후 효과가 커지는 패턴이 보인다.

왜 2순위인가:

- 스마트폰/디지털 미디어와 청소년 우울은 CHB 계열에서 강한 주제다.
- 단순 CLPM이 아니라 RI-CLPM이라 between-person confounding을 줄였다는 방법론적 장점이 있다.
- Paper A와 결과 영역이 겹치지만, 저널 시장은 다르다.

반드시 보완할 점:

- semopy 결과만으로는 부족하다. R lavaan 또는 Mplus로 fit indices(CFI/TLI/RMSEA/SRMR)를 산출해야 한다.
- 스마트폰 의존과 우울의 longitudinal measurement invariance를 검정해야 한다.
- 우울척도에서 자살사고 문항 포함 여부를 명시하고 민감도 분석을 해야 한다.
- 성별 다집단 RI-CLPM을 넣으면 기여가 커진다.
- 공변량으로 SES, 부모양육, 또래/교사관계를 넣은 모형과 무공변량 모형을 비교한다.

권장 원고 제목:

**Does Smartphone Dependence Precede Depression? A Seven-Wave Random-Intercept Cross-Lagged Panel Study of Korean Adolescents**

### 3순위: Paper D

현재 핵심 결과:

- W1, W4, W7에서 36개 증상 GGM을 구성했다.
- W7에서 E03 걱정, C08 발열감, D02 부끄럼, C06 피곤 등이 중심성이 높다.
- E04 자살사고의 bridge strength가 W1 .085에서 W7 .246으로 증가했다.
- A07, C08, C03, A01, E04의 bridge 증가가 크다.

왜 3순위인가:

- 우울/주의/공격/신체/위축 증상의 network development는 신선하다.
- Paper A의 후속 또는 companion paper로 좋다.
- 자살사고가 후기 청소년기에 cross-domain bridge로 부상한다는 메시지는 매력적이다.

리스크:

- 현재 GraphicalLassoCV 기반이라 심리 네트워크 표준인 EBICglasso/qgraph/bootnet과 다르다.
- Network Comparison Test가 아직 없다.
- 안정성 계수(CS-coefficient), edge accuracy bootstrap이 없다.
- density가 약 .6으로 높아 네트워크가 너무 조밀할 수 있다.

보강 후 투고:

- R `bootnet`, `qgraph`, `NetworkComparisonTest`, `networktools`로 재분석.
- W1-W7 차이를 permutation 기반으로 검정.
- Paper A가 제출된 뒤 별도 원고로 분리하는 것이 안전하다.

### 4순위: Cyber Delinquency Developmental Cascade

현재 핵심 결과:

- 직전 파형 사이버비행 미경험자를 위험집합으로 둔 GEE onset model에서 다음 파형 사이버비행 진입 위험은 남학생, 스마트폰 의존, 공격성, 현실비행, 부모 거부가 설명한다.
- 예: 스마트폰 의존 OR=1.11, 공격성 OR=1.12, 현실비행 OR=1.06, 부모거부 OR=1.09, 여학생 OR=.57.
- ML holdout AUC는 .63 수준이다.

왜 4순위인가:

- 범죄학/청소년비행 분야에서는 충분히 논문화 가능하다.
- KCYPS의 사이버비행 15문항 7-wave 반복은 장점이다.
- 다만 효과크기가 크지 않고, 고정효과 모형의 결과가 덜 선명하다.

보완 방향:

- 사이버비행 onset/discontinuation/recurrent event의 이산시간 생존분석으로 재설계.
- 학교폭력 W6 모듈과 연결하거나, 스마트폰 의존과 공격성의 상호작용을 검토.
- K-means trajectory는 Mplus LCGA/zero-inflated Poisson GBTM으로 대체.
- Paper B와 통합하지 말고 범죄학 저널용 별도 논문으로 두는 편이 낫다.

권장 타겟:

- Youth & Society
- Journal of Interpersonal Violence
- Children and Youth Services Review
- Crime & Delinquency는 보강 후 가능

### 5순위: Paper C

현재 핵심 결과:

- W1 baseline feature만으로 W7 자살사고 AUC=.59, 고우울 AUC=.57, 고위험궤적 AUC=.65.
- 가장 좋은 결과도 clinical prediction paper로는 modest하다.

왜 5순위인가:

- ML/SHAP 자체는 최신성은 있지만 AUC가 높지 않다.
- "예측 가능하다"보다 "초기 단면변수만으로 6년 후 위험 예측은 제한적이다"에 가까운 결과다.
- 단독 SSCI보다는 Paper A의 부록/보조분석으로 가치가 크다.

보강 방향:

- baseline scale-level이 아니라 item-level, parent-report, time-use, school, region 변수를 모두 투입.
- W1만이 아니라 W1-W3 누적 정보로 W7 예측하는 landmark model로 바꾼다.
- calibration plot, decision curve analysis, net benefit를 제시한다.
- 예측 대상은 자살사고보다 high-risk trajectory membership이 낫다.

## 실행 로드맵

### 이번 주 바로 할 일

1. Paper A의 우울척도에서 자살사고 문항 E04를 제외한 9문항 우울 궤적을 재산출한다.
2. Paper A K=2-5를 Mplus로 재현한다.
3. K=4 선택 논리를 확정한다: BIC, entropy, class size, interpretability, posterior probability.
4. BCH/3-step으로 W7 자살사고, 삶의 만족, 비행, 자아존중감을 distal outcome으로 비교한다.
5. Paper A용 Figure 1 trajectory, Table 1 class profile, Table 2 predictors/distal outcomes를 완성한다.

### Paper A가 안정되면

1. Paper B를 lavaan RI-CLPM으로 재추정한다.
2. 성별 다집단과 측정동일성 검정을 추가한다.
3. Paper D는 R 네트워크 표준 패키지로 재분석한다.
4. Paper C는 단독 논문이 아니라 Paper A의 prediction extension으로 보류한다.

## 최종 판단

현재 상태에서 “가장 빨리, 가장 높은 저널에, 가장 명확한 스토리로” 갈 수 있는 원고는 **Paper A**다. Paper B는 방법론 재현만 통과하면 두 번째 독립 논문으로 충분하다. Paper D는 참신하지만 통계 검정 보강이 선행되어야 하고, Paper C와 사이버비행 패키지는 보조 또는 후속 논문으로 두는 것이 효율적이다.
