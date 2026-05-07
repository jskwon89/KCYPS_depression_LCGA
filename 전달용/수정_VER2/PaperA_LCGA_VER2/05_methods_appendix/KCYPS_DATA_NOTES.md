# Paper A — KCYPS 2018 Data Notes

자료 / 척도 / 문항 정의 정리. 본문 §2.1, §2.2 작성 시 참조.

---

## 1. KCYPS 2018 자료 개요

- **Sponsor**: 한국청소년정책연구원 (National Youth Policy Institute, NYPI)
- **Cohort**: m1 (중학교 1학년 cohort) — 2018년 중1 재학생
- **Sampling**: 2-stage stratified cluster sampling
  - Strata: 17 시·도 (regions) × 3 도시 규모 (urban size)
  - Stage 1: school cluster 추출
  - Stage 2: 학교 내 학생 추출
- **Original m1 cohort N**: 2,590
- **Waves**: 2018 (W1) ~ 2024 (W7), annual face-to-face/online survey
- **Sample retention**: W2 94.1%, W3 92.1%, W4 87.5%, W5 86.9%, W6 85.9%, W7 80.0%
- **Public-use data**: 청소년(Y), 보호자(P), 형제자매(S) 별도 .dta 파일 (2025.12판)

## 2. 분석에 사용된 척도

### 2.1 우울 (Primary outcome trajectory)
- **Variable prefix**: `YPSY4E` (YPSY4 = 정서·사회 적응 척도; E = 우울 하위)
- **9-item primary** (E04 자살사고 항목 제외):
  - YPSY4E01 "기운이 별로 없다"
  - YPSY4E02 "불행하다고 생각하거나 슬퍼하고 우울해한다"
  - YPSY4E03 "걱정이 많다"
  - ~~YPSY4E04 "죽고 싶은 생각이 든다"~~ ← distal outcome으로만 사용
  - YPSY4E05 "울기를 잘한다"
  - YPSY4E06 "어떤 일이 잘못 되었을 때 나 때문이라는 생각을 자주 한다"
  - YPSY4E07 "외롭다"
  - YPSY4E08 "모든 일에 관심과 흥미가 없다"
  - YPSY4E09 "장래가 희망적이지 않은 것 같다"
  - YPSY4E10 "모든 일이 힘들다"
- **Response scale**: 1=전혀 그렇지 않다, 2=그렇지 않은 편이다, 3=그런 편이다, 4=매우 그렇다
- **Score**: mean of 9 items if ≥5 items answered
- **Reliability (9-item, 7 waves)**: α = .914, .908, .908, .899, .891, .888, .898

### 2.2 자살사고 (Distal outcome)
- **Variable**: `YPSY4E04w{1..7}` "죽고 싶은 생각이 든다"
- **Binarizations**:
  - **any** (≥2): "sometimes or more"
  - **strong** (≥3): "often or more"

### 2.3 Skinner PSCQ-A 6차원 양육 (W1 covariate)
- **Prefix**: `YFAM2` (24 items, 6 dimensions × 4 items)
- 6 dimensions:
  - A: **Warmth** (온정) — 4 items
  - B: **Rejection** (거부) — 4 items
  - C: **Autonomy support** (자율성 지지) — 4 items
  - D: **Coercion** (강요) — 4 items
  - E: **Structure** (구조 제공) — 4 items
  - F: **Chaos** (혼돈) — 4 items
- **Source**: Skinner, Johnson & Snyder (2005) Parents as Social Context Questionnaire — Adolescent
- **Response**: 1-4 Likert
- **Reliability** range: α = .79–.91 across waves and dimensions

### 2.4 또래/교사/심리 (W1 covariates)
- **YEDU2 Peer attachment**: 13 items, α=.85
- **YEDU3 Teacher relationship**: 14 items, α=.91
- **YPSY3 Self-esteem**: Rosenberg 10-item, α=.87
- **YPSY1 Life satisfaction**: Diener SWLS 5-item, α=.85
- **YPSY2 Subjective happiness**: Lyubomirsky-Lepper 4-item, α=.80
- **YPSY7 Grit**: Duckworth Short Grit 8-item, α=.71

### 2.5 스마트폰 의존 (W1 covariate)
- **Prefix**: `YMDA1C01..C15` (15-item smartphone dependence scale)
- **Source**: 한국정보화진흥원 청소년 스마트폰 중독 자가진단 척도
- α=.88

### 2.6 기타 정서·행동 (YPSY4 5 subscales — internalizing/externalizing battery)
- A: **Attention problems** 7 items
- B: **Aggression** 6 items
- C: **Somatic symptoms** 8 items
- D: **Social withdrawal** 5 items
- E: **Depression** 10 items (E04 포함, 우리 primary는 E04 제외)

### 2.7 인구학·맥락
- **YGENDER**: 1=남, 2=여
- **YBRT1A** (birth year): 99% 2005년생 (W1=14세 한국나이)
- **ARA1A** (학교 시도, 17개)
- **ARA2A** (학교 도시규모: 1=대도시, 2=중소도시, 3=읍면지역)
- **WEIGHTB2w7**: W7 종단면 표준화 가중치 (sensitivity 분석)

## 3. 결측 코딩 정책

- KCYPS 시스템 결측: 음수 (-1, -9 등)
- "해당 없음": 99
- 분석에서: `(value < 1) | (value > 4)` → NaN (Likert 1-4 척도)
- Scale composite: ≥50% 문항 응답한 경우만 평균 계산
- Trajectory analysis filter: ≥3 valid waves of dep9

## 4. 윤리/IRB

- 연구진 측 IRB 면제 또는 secondary data exemption 신청 필요
- KCYPS 자체는 NYPI에서 IRB 승인 받음, 데이터 사용 등록 필요
- 본 연구는 deidentified public-use 자료의 secondary analysis

## 5. 자료 활용 시 주의 (작성담당자용)

1. **Wave label**: KCYPS는 "1차~7차" 또는 "Wave 1~7"로 표기. 영문 본문은 "Wave 1 (2018) through Wave 7 (2024)" 또는 "W1–W7" 통일.
2. **m1 vs e1**: 본 연구는 **m1 (middle-school grade 1 cohort)**만 사용. e1 (elementary 4 cohort)은 별도.
3. **본 패널의 cohort 연령**: W1 (중1)=대부분 만 13세 (Korean age 14, 2005년생), W7 (대1 또는 고3 직후)=만 19세 (Korean age 20).
4. **Korean age vs 만 나이**: 한국 나이로 보면 W1 14세 → W7 20세, 만 나이로는 13세 → 19세. 본문은 **만 나이 (international age)** 통일 권장.

## 6. 자료 인용 형식

> Korean Children and Youth Panel Survey 2018 (KCYPS 2018), middle-school first-grade (m1) cohort, Waves 1–7 (2018–2024). National Youth Policy Institute, Republic of Korea. Public-use file release 2025.12.

---

상세 코드북: `D:/2026/SCI/청소년패널조사/output/codebook/`
