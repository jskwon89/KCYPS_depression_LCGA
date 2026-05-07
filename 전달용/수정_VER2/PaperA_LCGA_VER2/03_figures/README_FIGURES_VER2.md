# Figures VER2 안내

기존 R 기반 Figure 1-6은 유지했습니다. Mplus 분석은 본문의 핵심 그림을 교체하기보다, 재현성 검증을 보여주는 보충 Figure S1로 추가하는 것이 가장 안전합니다.

## Main figures

- Figure 1-6: 기존 `03_figures/R_ggplot2/` 및 상위 폴더 그림을 그대로 사용
- Figure 6는 cluster-bootstrap Cox 결과가 반영된 BOOT 버전을 primary로 사용

## Added supplementary figure

위치: `03_figures/Mplus_validation/`

- `FigureS1_Mplus_replication_9item.png`
- `FigureS1_Mplus_replication_9item.tiff`
- `FigureS1_Mplus_replication_9item.pdf`

Figure S1 구성:

1. Python EM과 Mplus 7의 BIC curve overlap
2. AIC, BIC, aBIC의 절대 차이
3. K=4 label-aligned modal class counts 비교

권장 caption:

> Supplementary Figure S1. Cross-platform replication of the primary 9-item LCGA solution in Mplus 7. Panel A compares BIC values across K = 1-5 between Python EM and Mplus. Panel B shows absolute differences in AIC, BIC, and sample-size adjusted BIC. Panel C compares modal class counts after label alignment. Mplus replicated the selected K = 4 solution to rounding precision, with 99.96% modal agreement and Cohen's kappa = .9994.
