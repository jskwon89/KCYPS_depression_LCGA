# Paper A v4 — Reviewer P1×3 + P2 처리 (2026-05-08)

검토의견 4개 모두 처리, 그 중 가장 중요한 P1.3 (per-wave HR robust CI)이 cluster bootstrap으로 대체되어 narrative가 v3보다 더 강해졌습니다.

---

## P1.1 — Table 2 BLRT p값 정정 ✅
**문제**: B=30, zero exceedance라 add-one 규칙으로 최소 p = 1/31 ≈ .032인데 Table 2가 "p<.001"로 잘못 표기.

**해결**: Table 2 BLRT 칸을 "K=2 vs K=1: LRT_obs=2,174.7, B=30, **p=.032** (lowest attainable for B=30; null max=11.5)"로 정정. 스크립트(`27_publication_outputs_9item.py:116`)도 동일.

---

## P1.2 — Methods §2.3 BLRT B=20 → B=30 + 잘못된 p<.001 주장 정정 ✅
**문제**: Methods가 아직 "B=20"으로, 또한 "B=100 or 500이면 p<.001"이라고 잘못 적힘.

**해결**:
- B=20 → **B=30**
- 정확한 add-one 규칙 명시: "smallest p attainable with zero exceedances at B=30 is 1/31 ≈ .032 (with B=100 it would be ≈ .010, with B=500 ≈ .002, and **B ≥ 999** needed for p < .001)"

---

## P1.3 — Per-wave HR CI를 cluster bootstrap으로 정정 ✅ (최중요)

**문제**: lifelines `variance_matrix_`는 Hessian-inverse covariance이지 cluster-sandwich 가 아님. `cluster_col='ID', robust=True` 옵션은 summary table의 SE만 갱신하고 variance matrix는 갱신하지 않음. 따라서 v3의 "robust CI"는 실제로 robust가 아니었음.

**해결**: `30_cox_cluster_bootstrap.py` 작성 — ID 레벨 cluster bootstrap (B=300, ~30초 소요).
- 각 b에서 unique ID들을 with replacement로 resample → long-form 재구성 → Cox class×wave fit → per-wave HR 계산
- 95% percentile CI + two-sided percentile bootstrap p-value (add-one floor 적용)

### 결과 비교 (per-wave HR for C4 Late-Increasing)

| Wave | naive (independent rows) p | "lifelines robust" p (실제는 Hessian-based) | **Cluster bootstrap p** |
|------|----------------------------|---------------------------------------------|---------------------------|
| W1 | <.001 (HR 1.43, naive 잘못) | .093 (still Hessian-based) | **.10** (n.s. — narrative 일치) |
| W7 | <.001 | <.001 | **.003** (highly sig, narrower CI) |

### Per-wave HR (cluster bootstrap 95% percentile CI)

| Wave | Persistent-High (C2) | Decreasing (C3) | Late-Increasing (C4) |
|------|----------------------|------------------|------------------------|
| W1 | 12.31 [8.17, 20.67]*** | 8.96 [6.29, 14.89]*** | 1.43 [0.94, 2.30] (n.s.) |
| W4 | 9.14 [7.01, 12.82]*** | 3.48 [2.71, 4.73]*** | 3.26 [2.61, 4.42]*** |
| W7 | 6.79 [4.00, 14.29]*** | 1.35 [0.75, 2.55] (n.s.) | 7.45 [5.00, 12.97]*** |

`*** p ≤ .003 (bootstrap)`

### Class × wave 상호작용 항 (cluster bootstrap)

| 항 | β | 95% percentile CI | bootstrap p |
|-----|----|---------------------|--------------|
| c2_x_t | -0.099 | [-0.247, +0.061] | .22 (n.s.) |
| c3_x_t | -0.315 | [-0.476, -0.183] | **.003** |
| c4_x_t | +0.275 | [+0.149, +0.413] | **.003** |

→ **C3 hazard는 시간에 따라 유의하게 감소, C4 hazard는 유의하게 증가** (cluster-correct)

### v3 → v4 narrative 변화

- v3 (잘못된 robust): "interaction terms are not significant under cluster-robust SEs (c4_x_t robust p = .26; c3_x_t robust p = .19)" → 형식 검정 못 한다고 후퇴
- **v4 (cluster bootstrap)**: "c3_x_t와 c4_x_t 모두 cluster bootstrap에서 p = .003으로 유의" → 형식 검정으로 narrative 강화

이건 사실상 v4가 v3보다 더 강한 결과를 줄 수 있다는 뜻. "invisible at-risk" 메시지가 descriptive에서 inferential로 격상.

산출물:
- `cox_HR_per_wave_with_CI_9item_BOOT.csv`: per-wave HR + percentile CI + boot p
- `cox_interaction_bootstrap_dist_9item.csv`: B=300 interaction coefficient distributions
- `Figure6_time_varying_HR_9item_BOOT.{png,pdf}`: bootstrap CI 음영 포함

---

## P2 — R lcmm class label 매핑 명시화 ✅

**문제**: R lcmm raw class 1-4는 트랜스 임의 순서이고 본문 C1-C4 manuscript 라벨과 다름 (R 1=C1 OK, R 2=C3, R 3=C2, R 4=C4 OK).

**해결**: `31_lcmm_label_alignment.py` 작성. 매뉴스크립트 §3.2에 명시적 매핑 표 추가:

| Manuscript label | Python EM raw class | Python n (%) | R lcmm raw class | R n (%) |
|------------------|---------------------|---------------|-------------------|----------|
| C1 Stable-Low | 1 | 615 (25.08%) | **1** | 616 (25.12%) |
| C2 Persistent-High | 2 | 245 (9.99%) | **3** | 245 (9.99%) |
| C3 Decreasing | 3 | 726 (29.61%) | **2** | 725 (29.57%) |
| C4 Late-Increasing | 4 | 866 (35.32%) | **4** | 866 (35.32%) |

**Modal n agreement: ±1 individual across all 4 classes** — Python EM과 R lcmm은 사실상 동일 해를 제공.

산출물:
- `lcmm_R/class_label_alignment_9item.csv`: R class → manuscript label
- `lcmm_R/python_vs_R_class_comparison_9item.csv`: 비교표
- `lcmm_R/class_pred_K4_9item_LABELED.csv`: 라벨 정렬된 trajectory means

---

## v4 매뉴스크립트 변경 요약

| 위치 | v3 → v4 |
|------|---------|
| Version 메타 | v3 → v4 (2026-05-08), 변경 4개 항목 명시 |
| Abstract | "p ≤ .003 under cluster bootstrap" 명시; 명확한 per-wave HR 보고 |
| §2.3 Methods | BLRT B=20 → B=30; add-one 규칙 정확하게 ("B ≥ 999 needed for p < .001"); time-to-event는 ID-level cluster bootstrap (B=300) |
| §3.2 | R lcmm 명시적 매핑 표 추가 |
| §3.6 | Hessian-based "robust" → cluster bootstrap; interaction term p=.003 (c3, c4) 보고; descriptive에서 inferential로 격상 |
| §3.7 | Cox sensitivity 항목을 cluster bootstrap 기반 표현으로 갱신 |
| §4.6 Limitations | "interaction terms are not significant under cluster-robust SEs" 삭제; 대신 "B≥1000으로 더 큰 부트스트랩 + sandwich estimator 비교 supplement 예정" |
| Reproducibility | scripts 30, 31 추가; Figure 6 BOOT를 PRIMARY로, ROBUST는 cautionary로 |

---

## 산출물 새 추가

```
primary_9item/
├─ figures/Figure6_time_varying_HR_9item_BOOT.{png,pdf}    ⭐ 새 (PRIMARY)
├─ long_distal/
│   ├─ cox_HR_per_wave_with_CI_9item_BOOT.csv              ⭐
│   ├─ cox_HR_per_wave_9item_BOOT.csv                      ⭐
│   ├─ cox_interaction_bootstrap_dist_9item.csv            ⭐ B=300 distributions
│   └─ cluster_boot_run.log
└─ lcmm_R/
   ├─ class_label_alignment_9item.csv                      ⭐
   ├─ python_vs_R_class_comparison_9item.csv               ⭐
   └─ class_pred_K4_9item_LABELED.csv                      ⭐
```

---

## 게재 준비도 (수정 평가)

| 단계 | 평가 |
|------|------|
| v1 (10-item, criterion contamination) | 55-65% |
| v2 (9-item primary) | 68-75% |
| v3 (BLRT 정정 + R lcmm 완료 + n 표기 + Hessian-based "robust") | 사용자 평가 72-78% |
| **v4 (BLRT p 정정 + cluster bootstrap + lcmm 매핑)** | **78-83%** (사용자 예상치) |

남은 작업 (게재까지 1-2주):
- Discussion 1,500단어 본문 (현재 outline)
- Reference list (50-70개)
- IRB / author info
- Cover letter
- (option) BLRT B=999 (1/1000 p < .001 보고 가능)
- (option) Cluster bootstrap B=1000 (현재 B=300)
- (option) Sandwich vs cluster-bootstrap 비교 supplement

v4는 "박사학위 검수 → 학회 발표 → 투고" 흐름의 두 번째 단계까지 안정.
