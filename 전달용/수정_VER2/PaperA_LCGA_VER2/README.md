# PaperA_LCGA 수정_VER2 안내

생성일: 2026-05-08

이 폴더는 기존 전달본 `PaperA_LCGA`를 그대로 복사한 뒤, Mplus 재현분석 이후 수정 및 보완된 파일만 추가/교체한 작성자 전달용 VER2입니다. 기존 전달본은 건드리지 않았습니다.

## 작성 시 우선 사용할 파일

1. 원고 최신본: `01_manuscript/MANUSCRIPT_v4.3_MPLUS.md`
2. 작성자 메모: `01_manuscript/작성자_전달메모_MPLUS_VER2.md`
3. Mplus 보강 진행보고: `01_manuscript/PROGRESS_REPORT_v5_MPLUS.md`
4. 전체 review history: `01_manuscript/REVIEW_HISTORY_HANDOFF.md`
5. Mplus 보충그림: `03_figures/Mplus_validation/FigureS1_Mplus_replication_9item.*`

## VER2에서 바뀐 핵심

- Mplus 7로 9-item primary LCGA K=1-5를 독립 재현했습니다.
- Mplus K=4 BIC = 25,228.08로 Python EM BIC = 25,228.07과 반올림 수준에서 일치했습니다.
- K=4 modal class assignment도 label alignment 후 99.96% 일치했고 Cohen's kappa = .9994였습니다.
- Table 2의 AIC, BIC, aBIC, entropy가 Mplus 보충표와 대조 가능하도록 정리되었습니다.
- Table 6은 cluster-bootstrap Cox per-wave HR 표가 최신본입니다. 예전 average HR/Schoenfeld 중심 표를 사용하지 마세요.
- Figure 1-6의 기존 R 그림은 유지하고, Mplus 재현성 검증용 Figure S1을 보충그림으로 추가했습니다.

## 새로 추가한 주요 파일

- `02_tables/TableS_Mplus_fit_replication_9item.csv`: K=1-5 Mplus vs Python fit indices 비교
- `02_tables/TableS_Mplus_class_alignment_9item.csv`: Python vs Mplus K=4 class label alignment 및 modal n 비교
- `03_figures/Mplus_validation/FigureS1_Mplus_replication_9item.png`
- `03_figures/Mplus_validation/FigureS1_Mplus_replication_9item.tiff`
- `03_figures/Mplus_validation/FigureS1_Mplus_replication_9item.pdf`
- `04_results_data/mplus/`: Mplus input/output/class probability/final report 전체
- `06_reproducibility/scripts/33_prepare_mplus_lcga_9item.py`
- `06_reproducibility/scripts/34_parse_mplus_lcga_9item.py`
- `06_reproducibility/scripts/35_mplus_replication_figure.py`

## 주의

Mplus는 새로운 primary 분석으로 제시하기보다, Python EM 및 R lcmm 결과를 교차검증하는 cross-platform replication으로 쓰는 것이 가장 안정적입니다. 원고에서는 "Python EM primary, R lcmm and Mplus cross-validation" 흐름을 유지하는 편이 좋습니다.
