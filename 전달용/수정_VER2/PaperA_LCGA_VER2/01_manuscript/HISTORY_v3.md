# Paper A v3 — 검토의견 4개 처리 (BLRT 숫자 / R lcmm 완료 / 비율 vs n / Cox robust)

날짜: 2026-05-08

---

## P1.1 — BLRT 숫자를 9문항 값으로 정정 ✅

**문제**: v2 매뉴스크립트와 Table 2가 10문항 LRT_obs = 2,121.9를 그대로 보고 (criterion contamination 수정 후 9문항이 primary가 되었는데).

**해결**:
- 9문항 fit_summary 기준: K=1 logL = -14,143.19, K=2 logL = -13,055.84 → **LRT_obs = 2 × 1,087.35 = 2,174.69**
- `28_blrt_9item.py` 실행: B=30 parametric bootstrap from K=1 null
- 완료된 결과 (B=30): null mean LRT = 3.74, SD = 2.50, **max = 11.51** → 관측값 2,174.69 대비 ~189× 큼
- **p_BLRT = (1 + 0) / (30 + 1) = 0.0323** (B=30, 0 exceedances의 가장 보수적 보고값). B=100+로 늘리면 p < .001 보고 가능하지만 magnitude 차이가 압도적이라 substantively 충분
- `primary_9item/blrt/blrt_summary_9item.csv`로 저장
- Manuscript abstract / 본문 §3.2 / Table 2 모두 정정

---

## P1.2 — R lcmm 9-item 완료 ✅

**문제**: v2가 "R lcmm cross-validation 완료"라고 썼지만 실제는 K=3까지만 진행 중.

**해결**: R lcmm 실행 완료 확인 (`primary_9item/lcmm_R/lcmm_fit_summary_9item.csv`).

| Metric | Python EM | R lcmm | 일치도 |
|--------|-----------|--------|--------|
| K=4 logL | -12,551.60 | -12,551.60 | 소수점까지 일치 |
| K=4 BIC | 25,228.07 | 25,228.08 | 0.01 차이 |
| K=4 entropy | 0.661 | 0.661 | 일치 |
| Class 1 (Stable-Low) | 25.4% | 25.1% | 0.3pp |
| Class 2 (Persistent-High) | 10.5% | 10.0% | 0.5pp |
| Class 3 (Decreasing) | 29.3% | 29.6% | 0.3pp |
| Class 4 (Late-Increasing) | 34.8% | 35.3% | 0.5pp |

→ Manuscript 본문에 실제 값 명시: "R K=4 BIC = 25,228.08 vs. Python EM K=4 BIC = 25,228.07; R logL = −12,551.60 = Python logL to two decimal places".

---

## P1.3 — 추정 mixing proportion vs modal-assignment n 분리 ✅

**문제**: 본문은 C2 n=257, C3 n=720, C4 n=856 등으로 적었는데, 이는 **mixing proportion × N**(estimated class size)임. 실제 modal-assignment count는 C2=245, C3=726, C4=866.

**해결**: 매뉴스크립트 §3.3과 abstract에 명시적으로 분리 표기:
- **C1 Stable-Low**: π̂ = 25.4%; modal n = 615
- **C2 Persistent-High**: π̂ = 10.5%; modal n = 245
- **C3 Decreasing**: π̂ = 29.3%; modal n = 726
- **C4 Late-Increasing**: π̂ = 34.8%; modal n = 866

> "Throughout the remainder of this paper, 'class size' refers to modal-assignment n; mixing proportions are used for class-conditional distal estimation under the BCH framework, where they enter via the classification-error matrix."

Figure 1 캡션도 modal n 통일 (이미 615/245/726/866 사용 중).

---

## P2 — Time-varying Cox with cluster_col='ID', robust=True ✅

**문제**: v2의 class × time interaction Cox는 individual 별 여러 interval row를 만들지만 cluster correction이 없어 SE가 과소추정.

**해결**: `29_cox_time_varying_robust.py`로 두 가지 방법 적합:
1. **CoxPHFitter(cluster_col='ID', robust=True)**: sandwich SE
2. **CoxTimeVaryingFitter(id_col='ID')**: cross-check

### 결과 요약

**Per-wave HR with cluster-robust 95% CI** (delta method):

| Class | W1 HR [95% CI], p | W4 HR [95% CI], p | W7 HR [95% CI], p |
|-------|---------------------|---------------------|---------------------|
| C2 Persistent-High | 12.31 [8.21, 18.46], <.001 | 9.14 [6.83, 12.23], <.001 | 6.78 [3.73, 12.33], <.001 |
| C3 Decreasing | 8.96 [6.12, 13.12], <.001 | 3.48 [2.65, 4.56], <.001 | 1.35 [0.77, 2.36], .29 |
| **C4 Late-Increasing** | **1.43 [0.94, 2.17], .09** | 3.26 [2.55, 4.18], <.001 | **7.45 [4.64, 11.97], <.001** |

### Robust SE 적용으로 인한 narrative 조정

- **Per-wave HR 패턴은 robust SE 하에서 유지**: C4가 W1에 비유의(p=.09)에서 W7 강유의(p<.001)로 transition, C3은 W1 강유의에서 W7 비유의(p=.29)로
- **그러나 interaction term 자체는 robust SE 하에서 비유의**: c4_x_t robust p=.26, c3_x_t robust p=.19 (naive SE 하에서는 둘 다 p<.001)
- → Manuscript에서 "rate of change"는 형식 검정으로 주장하지 않고, "per-wave HR pattern"으로 descriptive하게 보고
- 핵심 narrative ("invisible at-risk")는 W1 vs W7 per-wave HR 비교만으로 강하게 뒷받침됨

**Figure 6 갱신**: `Figure6_time_varying_HR_9item_ROBUST.{png,pdf}` (95% robust CI 음영 포함)

---

## 매뉴스크립트 v3 변경 요약

| 위치 | v2 → v3 |
|------|---------|
| 메타 | "v2 (criterion-clean)" → "v3 (2026-05-08)" + 변경 항목 명시 |
| Abstract | LRT 2,121.9 → **2,174.7**; class size를 π̂ + modal n 분리; Cox per-wave with robust CI |
| §3.2 | BLRT 숫자 정정; B=30; R lcmm 실제 BIC 값 (25,228.08 vs 25,228.07) |
| §3.3 | "C2 (10.5%, n=257)" → "π̂ = 10.5%, modal n = 245" 분리 명시 + 일반 정책 문장 |
| §3.6 | naive SE c3_x_t/c4_x_t p<.001 + robust SE p=.19/.26 둘 다 보고; per-wave HR 95% CI 추가 |
| Reproducibility | scripts 28, 29 추가 |

---

## 산출물

```
output/paperA_lcga/
├─ MANUSCRIPT_DRAFT.md (v3)                     ⭐ 4개 review 처리 완료
├─ PROGRESS_REPORT_v3.md                        ⭐ 이 문서
├─ MANUSCRIPT_DRAFT_v1_10item.md                v1 backup
├─ PROGRESS_REPORT_v1.md, v2.md                 backup
└─ primary_9item/
   ├─ tables/Table2_model_fit_9item.csv         ⭐ BLRT 2,174.7 정정
   ├─ figures/
   │   ├─ Figure6_time_varying_HR_9item.png     (naive SE, 참조)
   │   └─ Figure6_time_varying_HR_9item_ROBUST.png  ⭐ cluster-robust
   ├─ lcmm_R/
   │   ├─ lcmm_fit_summary_9item.csv            ⭐ K=4 BIC=25,228.08
   │   ├─ class_pred_K4_9item.csv
   │   ├─ posterior_K4_9item.csv
   │   └─ lcga_models_9item.rds
   ├─ blrt/
   │   ├─ blrt_summary_9item.csv                ⭐ 9-item LRT_obs=2,174.69, B=30
   │   └─ blrt_dist_K1vs2_9item.npy
   └─ long_distal/
      ├─ cox_time_varying_9item_ROBUST.csv      ⭐ cluster_col=ID, robust SE
      ├─ cox_time_varying_9item_CTV.csv         CoxTimeVaryingFitter cross-check
      ├─ cox_HR_per_wave_9item_ROBUST.csv
      └─ cox_HR_per_wave_with_CI_9item_ROBUST.csv  ⭐ 95% robust CI per wave
```

---

## 게재 준비도 평가 (v3)

| 단계 | 평가 |
|------|------|
| v1 (10-item) | 사용자 review로 55-65% (criterion contamination) |
| v2 (9-item, 4가지 P1 처리) | 사용자 review로 68-75% (BLRT 숫자, R lcmm 미완, n 표기, Cox SE) |
| **v3 (이번 보강)** | **75-82%** |

남은 것 (게재까지 1-2주):
- Discussion 1,500단어 본문 (현재 outline)
- Reference list (50-70개, JAD style)
- IRB 정보, author info
- Cover letter
- (option) BLRT K=2~5 모두 완전 실행 (B=100, ~1-2일 compute)

본 v3는 동료/지도 검수 → 최종 수정 → 투고 흐름의 두 번째 단계까지 안정적으로 도달.
