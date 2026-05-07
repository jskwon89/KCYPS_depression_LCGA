# ==============================================================================
# Paper A — Publication figures in R (ggplot2)
# Mirrors the Python matplotlib figures but with publication-grade ggplot2 styling.
# Use these versions for journal submission.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(scales); library(patchwork); library(ggsci); library(cowplot)
})

ROOT <- 'D:/2026/SCI/청소년패널조사'
P9   <- file.path(ROOT, 'output', 'paperA_lcga', 'primary_9item')
FIGS <- file.path(P9, 'figures_R')
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

# Class metadata
class_levels <- c('C1 Stable-Low', 'C2 Persistent-High',
                  'C3 Decreasing', 'C4 Late-Increasing')
class_colors <- c('C1 Stable-Low'      = '#2ca02c',
                  'C2 Persistent-High' = '#d62728',
                  'C3 Decreasing'      = '#1f77b4',
                  'C4 Late-Increasing' = '#ff7f0e')

# Publication theme (JAD/APA-friendly)
pub_theme <- theme_classic(base_size = 11, base_family = 'sans') +
  theme(
    plot.title       = element_text(face = 'bold', size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 10, color = 'grey30'),
    axis.title       = element_text(size = 12, face = 'bold'),
    axis.text        = element_text(size = 10, color = 'black'),
    axis.line        = element_line(size = 0.5),
    axis.ticks       = element_line(size = 0.4),
    panel.grid.major = element_line(color = 'grey92', size = 0.3),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 10, face = 'bold'),
    legend.text      = element_text(size = 9),
    legend.position  = 'top',
    legend.background = element_blank(),
    legend.key       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = 'bold', size = 11),
    plot.caption     = element_text(size = 8, color = 'grey40', hjust = 0)
  )

save_pub <- function(p, name, w = 8, h = 5.5) {
  # 600 dpi for both raster formats; vector PDF for unlimited zoom
  ggsave(file.path(FIGS, paste0(name, '.png')), p, width = w, height = h,
         dpi = 600, bg = 'white')
  ggsave(file.path(FIGS, paste0(name, '.pdf')), p, width = w, height = h,
         device = cairo_pdf, bg = 'white')
  ggsave(file.path(FIGS, paste0(name, '.tiff')), p, width = w, height = h,
         dpi = 600, bg = 'white', compression = 'lzw')
}

# ============================================================
# Figure 1 — 4-class depression trajectories (cluster-bootstrap CI envelopes)
# ============================================================
class_assign <- read_csv(file.path(P9, 'class_assignments_K4_9item.csv'),
                          show_col_types = FALSE) %>%
  rename(class_K4 = class_K4_9item)

# Long format depression data per ID per wave
# Build trajectory means + bootstrap CI per wave per class
# Use saved long_panel data via Python parquet — but R can't easily read parquet.
# Instead recompute means directly from raw .dta files (or use the per-class
# predicted means from class_profiles).
profiles <- read_csv(file.path(P9, 'class_profiles_K4_9item.csv'),
                      show_col_types = FALSE)
profiles_long <- profiles %>%
  select(class, proportion, starts_with('mean_w')) %>%
  pivot_longer(cols = starts_with('mean_w'),
               names_to = 'wave', values_to = 'depression') %>%
  mutate(wave = as.integer(gsub('mean_w', '', wave))) %>%
  mutate(class_label = factor(case_when(
    class == 1 ~ 'C1 Stable-Low',
    class == 2 ~ 'C2 Persistent-High',
    class == 3 ~ 'C3 Decreasing',
    class == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels))

# Bootstrap CI per class per wave: read suic_long for class assignments and
# pull observed depression means from desc panel via R's readstat
# For simplicity, compute means + 95% CI at each wave using haven::read_dta
suppressPackageStartupMessages(library(haven))
items_9 <- c(1,2,3,5,6,7,8,9,10)
long_dep <- list()
for (w in 1:7) {
  d <- read_dta(file.path(ROOT, sprintf('KCYPS2018m1Yw%d.dta', w)))
  cols <- sprintf('YPSY4E%02dw%d', items_9, w)
  m <- as.data.frame(d[, cols])
  m[m < 1 | m > 4] <- NA
  n_obs <- rowSums(!is.na(m))
  dep <- ifelse(n_obs >= 5, rowMeans(m, na.rm = TRUE), NA_real_)
  long_dep[[w]] <- data.frame(ID = d$ID, wave = w, dep = dep)
}
ldep <- do.call(rbind, long_dep) %>%
  filter(!is.na(dep)) %>%
  inner_join(class_assign %>% select(ID, class_K4), by = 'ID') %>%
  mutate(class_label = factor(case_when(
    class_K4 == 1 ~ 'C1 Stable-Low',
    class_K4 == 2 ~ 'C2 Persistent-High',
    class_K4 == 3 ~ 'C3 Decreasing',
    class_K4 == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels))

# Bootstrap CI
set.seed(42)
B <- 500
boot_means <- ldep %>%
  group_by(class_label, wave) %>%
  summarise(mean_obs = mean(dep), n = n(), .groups = 'drop')

boot_results <- list()
for (k in class_levels) {
  for (w in 1:7) {
    sub <- ldep %>% filter(class_label == k, wave == w) %>% pull(dep)
    if (length(sub) > 1) {
      boot_vals <- replicate(B, mean(sample(sub, length(sub), replace = TRUE)))
      boot_results[[paste(k, w)]] <- data.frame(
        class_label = k, wave = w,
        mean_boot = mean(boot_vals),
        ci_lo = quantile(boot_vals, 0.025),
        ci_hi = quantile(boot_vals, 0.975)
      )
    }
  }
}
ci_df <- bind_rows(boot_results) %>%
  mutate(class_label = factor(class_label, levels = class_levels))

# Class size labels
n_class <- class_assign %>% count(class_K4) %>%
  mutate(class_label = factor(case_when(
    class_K4 == 1 ~ 'C1 Stable-Low',
    class_K4 == 2 ~ 'C2 Persistent-High',
    class_K4 == 3 ~ 'C3 Decreasing',
    class_K4 == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels)) %>%
  left_join(profiles %>% select(class, proportion) %>%
              mutate(class_label = factor(case_when(
                class == 1 ~ 'C1 Stable-Low',
                class == 2 ~ 'C2 Persistent-High',
                class == 3 ~ 'C3 Decreasing',
                class == 4 ~ 'C4 Late-Increasing'),
                levels = class_levels)) %>%
              select(class_label, proportion), by = 'class_label')

class_legend <- n_class %>%
  mutate(label = sprintf('%s (%.1f%%, n=%d)',
                         class_label, proportion*100, n)) %>%
  pull(label, name = class_label)

p1 <- ci_df %>%
  inner_join(boot_means %>% select(class_label, wave, mean_obs),
              by = c('class_label', 'wave')) %>%
  ggplot(aes(x = wave, y = mean_obs, color = class_label, fill = class_label)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2.0, ymax = 3.0,
           alpha = 0.05, fill = '#d62728') +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, shape = 21, color = 'white', stroke = 1.2,
              aes(fill = class_label)) +
  scale_color_manual(values = class_colors, labels = class_legend) +
  scale_fill_manual(values = class_colors, labels = class_legend) +
  scale_x_continuous(breaks = 1:7,
                      labels = c('W1\n(13y)', 'W2\n(14y)', 'W3\n(15y)',
                                 'W4\n(16y)', 'W5\n(17y)', 'W6\n(18y)',
                                 'W7\n(19y)')) +
  scale_y_continuous(limits = c(1.0, 3.0), breaks = seq(1.0, 3.0, 0.25)) +
  labs(title = 'Figure 1. LCGA of Korean adolescent depression trajectories',
       subtitle = 'KCYPS-2018 m1 cohort, N=2,452, 9-item criterion-clean depression composite',
       x = 'Wave (W1=2018 Grade 7 — W7=2024 ≈ 1 yr post Grade 12)',
       y = 'Depression score (mean of 9 items, 1–4)',
       color = 'Latent class', fill = 'Latent class',
       caption = paste('Bands: 95% nonparametric bootstrap CI (B = 500).',
                       'Shaded clinical range: depression score ≥ 2.0.', sep='\n')) +
  pub_theme +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         fill  = guide_legend(nrow = 2, byrow = TRUE))

save_pub(p1, 'Figure1_LCGA_trajectories_R', w = 9, h = 6)
cat('Figure 1 (R) saved.\n')

# ============================================================
# Figure 2 — Suicidal ideation by class across waves + KM
# ============================================================
suic <- read_csv(file.path(P9, 'long_distal', 'suic_long.csv'),
                  show_col_types = FALSE) %>%
  filter(!is.na(class_K4)) %>%
  mutate(class_label = factor(case_when(
    class_K4 == 1 ~ 'C1 Stable-Low',
    class_K4 == 2 ~ 'C2 Persistent-High',
    class_K4 == 3 ~ 'C3 Decreasing',
    class_K4 == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels))

rates_any <- suic %>% group_by(class_label, wave) %>%
  summarise(rate = mean(suic_any, na.rm = TRUE), n = n(), .groups = 'drop')
rates_strong <- suic %>% group_by(class_label, wave) %>%
  summarise(rate = mean(suic_strong, na.rm = TRUE), .groups = 'drop')

p2a <- ggplot(rates_any, aes(x = wave, y = rate, color = class_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, shape = 21, fill = 'white', stroke = 1.2) +
  scale_color_manual(values = class_colors) +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(limits = c(0, 0.85), labels = percent_format(accuracy = 1)) +
  labs(title = 'A. Any suicidal ideation (≥2 = "sometimes")',
       x = 'Wave', y = 'Endorsement rate', color = NULL) +
  pub_theme + theme(legend.position = 'none')

p2b <- ggplot(rates_strong, aes(x = wave, y = rate, color = class_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, shape = 21, fill = 'white', stroke = 1.2) +
  scale_color_manual(values = class_colors) +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(limits = c(0, 0.40), labels = percent_format(accuracy = 1)) +
  labs(title = 'B. Strong suicidal ideation (≥3 = "often")',
       x = 'Wave', y = 'Endorsement rate', color = NULL) +
  pub_theme + theme(legend.position = 'none')

# KM curve via custom plot
surv <- read_csv(file.path(P9, 'long_distal', 'survival_data.csv'),
                  show_col_types = FALSE) %>%
  mutate(class_label = factor(case_when(
    class_K4 == 1 ~ 'C1 Stable-Low',
    class_K4 == 2 ~ 'C2 Persistent-High',
    class_K4 == 3 ~ 'C3 Decreasing',
    class_K4 == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels))

km_data <- surv %>%
  group_by(class_label) %>%
  arrange(time) %>%
  mutate(at_risk = n() - cumsum(lag(event, default = 0))) %>%
  group_by(class_label, time) %>%
  summarise(n_event = sum(event), n_at_risk = first(at_risk), .groups = 'drop') %>%
  group_by(class_label) %>%
  mutate(surv = cumprod(1 - n_event / n_at_risk),
         se = sqrt(surv^2 * cumsum(n_event / (n_at_risk * (n_at_risk - n_event)))),
         ci_lo = pmax(0, surv - 1.96 * se),
         ci_hi = pmin(1, surv + 1.96 * se)) %>%
  ungroup()

# Add wave 0 baseline for proper step plot
km_baseline <- km_data %>%
  group_by(class_label) %>% slice(1) %>%
  mutate(time = 0, surv = 1, se = 0, ci_lo = 1, ci_hi = 1) %>%
  ungroup()
km_full <- bind_rows(km_baseline, km_data)

p2c <- ggplot(km_full, aes(x = time, y = surv, color = class_label, fill = class_label)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, color = NA) +
  geom_step(linewidth = 1.0) +
  scale_color_manual(values = class_colors) +
  scale_fill_manual(values = class_colors) +
  scale_x_continuous(breaks = 0:7, limits = c(0, 7)) +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format(accuracy = 1)) +
  labs(title = 'C. Kaplan-Meier: time to first strong SI',
       x = 'Wave', y = 'P(free of strong SI)', color = NULL, fill = NULL) +
  pub_theme + theme(legend.position = 'none')

p2 <- (p2a + p2b + p2c) +
  plot_annotation(
    title = 'Figure 2. Suicidal ideation by depression trajectory class',
    subtitle = 'Across-wave endorsement rates and time-to-event analysis',
    theme = theme(plot.title    = element_text(face = 'bold', size = 14),
                  plot.subtitle = element_text(size = 11, color = 'grey30'))
  ) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')
# Add a shared bottom legend
p2_legend <- p2a + theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1))
legend <- cowplot::get_legend(p2_legend)
p2_final <- cowplot::plot_grid(
  cowplot::plot_grid(p2a, p2b, p2c, ncol = 3, align = 'h'),
  legend, ncol = 1, rel_heights = c(1, 0.08))

ggsave(file.path(FIGS, 'Figure2_suicidal_ideation_R.png'),
       p2_final, width = 16, height = 6, dpi = 600, bg = 'white')
ggsave(file.path(FIGS, 'Figure2_suicidal_ideation_R.pdf'),
       p2_final, width = 16, height = 6, device = cairo_pdf, bg = 'white')
ggsave(file.path(FIGS, 'Figure2_suicidal_ideation_R.tiff'),
       p2_final, width = 16, height = 6, dpi = 600, bg = 'white',
       compression = 'lzw')
cat('Figure 2 (R) saved.\n')

# ============================================================
# Figure 6 — Time-varying class HR with cluster bootstrap CI
# ============================================================
boot_hr <- read_csv(file.path(P9, 'long_distal',
                                'cox_HR_per_wave_with_CI_9item_BOOT.csv'),
                     show_col_types = FALSE) %>%
  mutate(class_label = factor(case_when(
    class == 'Persistent-High' ~ 'C2 Persistent-High',
    class == 'Decreasing'      ~ 'C3 Decreasing',
    class == 'Late-Increasing' ~ 'C4 Late-Increasing'),
    levels = class_levels[2:4]))

p6 <- ggplot(boot_hr, aes(x = wave, y = HR, color = class_label, fill = class_label)) +
  geom_ribbon(aes(ymin = HR_lo, ymax = HR_hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, shape = 21, color = 'white', stroke = 1.2,
             aes(fill = class_label)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50') +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20),
                limits = c(0.7, 25)) +
  scale_color_manual(values = class_colors[2:4]) +
  scale_fill_manual(values = class_colors[2:4]) +
  scale_x_continuous(breaks = 1:7) +
  labs(title = 'Figure 6. Time-varying class hazard ratios for first strong SI',
       subtitle = paste0(
         'Cox class × wave model with ID-level cluster bootstrap (B=300); ',
         'reference: C1 Stable-Low'),
       x = 'Wave (1–7)',
       y = 'Hazard ratio (log scale, 95% cluster-bootstrap CI)',
       color = NULL, fill = NULL,
       caption = paste('p-values floored at 1/(B+1) = 1/301 = .003 (B = 300).',
                       'A planned B = 1000 supplementary bootstrap would lower the floor to .001.',
                       sep='\n')) +
  pub_theme +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

save_pub(p6, 'Figure6_time_varying_HR_R', w = 9, h = 6)
cat('Figure 6 (R) saved.\n')

# ============================================================
# Figure 5 — Information criteria + entropy
# ============================================================
fit <- read_csv(file.path(P9, 'fit_summary_9item.csv'),
                 show_col_types = FALSE)
fit_long <- fit %>%
  select(K, BIC, aBIC, AIC) %>%
  pivot_longer(cols = c(BIC, aBIC, AIC), names_to = 'criterion', values_to = 'value')

p5a <- ggplot(fit_long, aes(x = K, y = value, color = criterion, shape = criterion)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 4, linetype = 'dashed', color = '#d62728') +
  annotate('text', x = 4.05, y = max(fit_long$value)*0.98,
            label = 'K = 4 selected', color = '#d62728', size = 3.5, hjust = 0) +
  scale_color_manual(values = c('AIC' = '#2ca02c', 'BIC' = '#1f77b4', 'aBIC' = '#ff7f0e')) +
  scale_x_continuous(breaks = 1:5) +
  labs(title = 'A. Information criteria',
        x = 'Number of latent classes (K)', y = 'IC value', color = NULL, shape = NULL) +
  pub_theme

p5b <- ggplot(fit, aes(x = K, y = entropy)) +
  geom_line(linewidth = 1, color = '#d62728') +
  geom_point(size = 3.5, color = '#d62728') +
  geom_hline(yintercept = 0.7, linetype = 'dotted', color = 'grey50') +
  geom_vline(xintercept = 4, linetype = 'dashed', color = '#d62728') +
  annotate('text', x = 1.1, y = 0.71, label = 'Entropy ≥ 0.7 (good)',
            size = 3, color = 'grey50', hjust = 0) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(limits = c(0.6, 1.05)) +
  labs(title = 'B. Classification entropy',
       x = 'Number of latent classes (K)', y = 'Entropy') +
  pub_theme

p5 <- (p5a + p5b) +
  plot_annotation(
    title = 'Figure 5. LCGA model fit indices for K = 1 to 5 (9-item primary)',
    theme = theme(plot.title = element_text(face = 'bold', size = 14)))
ggsave(file.path(FIGS, 'Figure5_fit_indices_R.png'), p5,
       width = 12, height = 5, dpi = 600, bg = 'white')
ggsave(file.path(FIGS, 'Figure5_fit_indices_R.pdf'), p5,
       width = 12, height = 5, device = cairo_pdf, bg = 'white')
ggsave(file.path(FIGS, 'Figure5_fit_indices_R.tiff'), p5,
       width = 12, height = 5, dpi = 600, bg = 'white', compression = 'lzw')
cat('Figure 5 (R) saved.\n')

# ============================================================
# Figure 3 — Baseline (W1) covariate z-score heatmap by class
# ============================================================
covars <- c('parent_warmth', 'parent_reject', 'parent_coerce', 'parent_chaos',
            'parent_autonomy', 'parent_structure', 'peer_attach', 'teacher_rel',
            'self_esteem', 'life_sat', 'happy', 'grit', 'depression9', 'somatic',
            'aggression', 'withdrawal', 'attention', 'sm_addiction')
covar_labels <- c(
  'parent_warmth'   = 'Parent: Warmth',
  'parent_reject'   = 'Parent: Rejection',
  'parent_coerce'   = 'Parent: Coercion',
  'parent_chaos'    = 'Parent: Chaos',
  'parent_autonomy' = 'Parent: Autonomy',
  'parent_structure'= 'Parent: Structure',
  'peer_attach'     = 'Peer attachment',
  'teacher_rel'     = 'Teacher relationship',
  'self_esteem'     = 'Self-esteem',
  'life_sat'        = 'Life satisfaction',
  'happy'           = 'Subjective happiness',
  'grit'            = 'Grit',
  'depression9'     = 'Depression (9-item) W1',
  'somatic'         = 'Somatic symptoms',
  'aggression'      = 'Aggression',
  'withdrawal'      = 'Social withdrawal',
  'attention'       = 'Attention problems',
  'sm_addiction'    = 'Smartphone dependence')

# Read W1 data with covariates from Python long_panel via .dta directly
df1 <- read_dta(file.path(ROOT, 'KCYPS2018m1Yw1.dta'))
# Compute scale composites for W1 (mirrors 06_build_long_panel.py)
score_scale <- function(d, items, reverse = NULL, lo = 1, hi = 4, mode = NULL) {
  m <- as.data.frame(d[, items])
  m[m < 1 | m > hi | m > 90] <- NA
  if (!is.null(reverse) && length(reverse) > 0) {
    for (i in reverse) m[, i] <- (lo + hi) - m[, i]
  }
  n_obs <- rowSums(!is.na(m))
  ifelse(n_obs >= max(2, ncol(m) %/% 2), rowMeans(m, na.rm = TRUE), NA_real_)
}

w1_score <- data.frame(ID = df1$ID)
w1_score$parent_warmth <- score_scale(df1, sprintf('YFAM2A%02dw1', 1:4))
w1_score$parent_reject <- score_scale(df1, sprintf('YFAM2B%02dw1', 1:4))
w1_score$parent_autonomy <- score_scale(df1, sprintf('YFAM2C%02dw1', 1:4))
w1_score$parent_coerce <- score_scale(df1, sprintf('YFAM2D%02dw1', 1:4))
w1_score$parent_structure <- score_scale(df1, sprintf('YFAM2E%02dw1', 1:4))
w1_score$parent_chaos <- score_scale(df1, sprintf('YFAM2F%02dw1', 1:4))
w1_score$peer_attach <- score_scale(df1, sprintf('YEDU2A%02dw1', 1:13),
                                      reverse = 9:13)
w1_score$teacher_rel <- score_scale(df1, sprintf('YEDU3A%02dw1', 1:14))
w1_score$self_esteem <- score_scale(df1, sprintf('YPSY3A%02dw1', 1:10),
                                      reverse = c(2, 5, 6, 8, 9))
w1_score$life_sat <- score_scale(df1, sprintf('YPSY1A0%dw1', 1:5))
w1_score$happy <- score_scale(df1, sprintf('YPSY2A0%dw1', 1:4),
                                reverse = 4, lo = 1, hi = 7)
w1_score$grit <- score_scale(df1, sprintf('YPSY7A0%dw1', 1:8),
                              reverse = c(1, 3, 5, 6), lo = 1, hi = 5)
w1_score$depression9 <- score_scale(df1, sprintf('YPSY4E%02dw1', items_9))
w1_score$somatic <- score_scale(df1, sprintf('YPSY4C%02dw1', 1:8))
w1_score$aggression <- score_scale(df1, sprintf('YPSY4B%02dw1', 1:6))
w1_score$withdrawal <- score_scale(df1, sprintf('YPSY4D%02dw1', 1:5))
w1_score$attention <- score_scale(df1, sprintf('YPSY4A%02dw1', 1:7))
w1_score$sm_addiction <- score_scale(df1, sprintf('YMDA1C%02dw1', 1:15),
                                       reverse = c(5, 10, 15))

w1m <- w1_score %>%
  inner_join(class_assign %>% select(ID, class_K4), by = 'ID') %>%
  mutate(class_label = factor(case_when(
    class_K4 == 1 ~ 'C1 Stable-Low',
    class_K4 == 2 ~ 'C2 Persistent-High',
    class_K4 == 3 ~ 'C3 Decreasing',
    class_K4 == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels))

# Compute z-score per covariate, then class means
z_long <- w1m %>%
  mutate(across(all_of(covars), ~ scale(.)[, 1])) %>%
  pivot_longer(all_of(covars), names_to = 'covar', values_to = 'z') %>%
  group_by(class_label, covar) %>%
  summarise(z_mean = mean(z, na.rm = TRUE), .groups = 'drop') %>%
  mutate(covar_label = factor(covar_labels[covar], levels = rev(covar_labels)))

p3 <- ggplot(z_long, aes(x = class_label, y = covar_label, fill = z_mean)) +
  geom_tile(color = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%+.2f', z_mean),
                color = abs(z_mean) > 0.5),
            size = 3.2, fontface = 'bold') +
  scale_fill_gradient2(low = '#1f4d7a', mid = 'white', high = '#a3211a',
                        midpoint = 0, limits = c(-1, 1),
                        oob = scales::squish,
                        name = 'z-score') +
  scale_color_manual(values = c(`TRUE` = 'white', `FALSE` = 'black'),
                      guide = 'none') +
  scale_x_discrete(labels = c('C1\nStable-Low', 'C2\nPersistent-High',
                                'C3\nDecreasing', 'C4\nLate-Increasing')) +
  labs(title = 'Figure 3. Baseline (W1) z-score profile by trajectory class',
       subtitle = 'Cell value = (class mean − grand mean) / SD',
       x = NULL, y = NULL) +
  pub_theme +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        axis.line = element_blank(),
        axis.ticks = element_blank())

save_pub(p3, 'Figure3_baseline_heatmap_R', w = 8, h = 9)
cat('Figure 3 (R) saved.\n')

# ============================================================
# Figure 4 — BCH-corrected ORs and Cox HRs forest plot
# ============================================================
or_any <- read_csv(file.path(P9, 'bch', 'bch_suic_any_OR.csv'),
                    show_col_types = FALSE) %>%
  mutate(group = 'Suicidal ideation any (BCH OR, W7)')
or_strong <- read_csv(file.path(P9, 'bch', 'bch_suic_strong_OR.csv'),
                       show_col_types = FALSE) %>%
  mutate(group = 'Suicidal ideation strong (BCH OR, W7)')
cox_hr <- read_csv(file.path(P9, 'long_distal', 'cox_HR_class_9item.csv'),
                    show_col_types = FALSE) %>%
  rename(label = `...1`) %>%
  mutate(class = case_when(
    grepl('c2', label) ~ 2L, grepl('c3', label) ~ 3L, grepl('c4', label) ~ 4L),
    OR_vs_ref = HR, CI_lo = HR_lower, CI_hi = HR_upper,
    group = 'Time to strong SI (Cox HR, average)') %>%
  select(class, OR_vs_ref, CI_lo, CI_hi, group)

# Combine
forest_df <- bind_rows(or_any, or_strong, cox_hr) %>%
  mutate(class_label = factor(case_when(
    class == 2 ~ 'C2 Persistent-High',
    class == 3 ~ 'C3 Decreasing',
    class == 4 ~ 'C4 Late-Increasing'),
    levels = class_levels[2:4]),
    label = paste0(group, ': ', class_label, ' vs C1'),
    group_short = factor(group, levels = c(
      'Suicidal ideation any (BCH OR, W7)',
      'Suicidal ideation strong (BCH OR, W7)',
      'Time to strong SI (Cox HR, average)'))) %>%
  arrange(group_short, class_label) %>%
  mutate(label = factor(label, levels = rev(label)))

p4 <- ggplot(forest_df, aes(x = OR_vs_ref, y = label, color = class_label)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey50') +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.5, shape = 21, fill = 'white', stroke = 1.5) +
  scale_color_manual(values = class_colors[2:4]) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 500),
                labels = c('0.5', '1', '2', '5', '10', '20', '50', '100', '500')) +
  facet_grid(group_short ~ ., scales = 'free_y', space = 'free_y',
              switch = 'y') +
  labs(title = 'Figure 4. BCH-corrected odds ratios and Cox hazard ratios',
       subtitle = 'Reference: C1 Stable-Low; effects on log scale',
       x = 'Effect size (OR or HR), log scale',
       y = NULL, color = NULL,
       caption = 'BCH bootstrap CI is 2.5/97.5 percentile of B=2,000 nonparametric resamples.') +
  pub_theme +
  theme(legend.position = 'top',
        strip.placement = 'outside',
        strip.text.y.left = element_text(angle = 0, hjust = 0, face = 'bold'),
        panel.spacing = unit(0.4, 'lines'))

save_pub(p4, 'Figure4_distal_forest_R', w = 11, h = 7)
cat('Figure 4 (R) saved.\n')

cat('\nAll R/ggplot2 figures saved to:', FIGS, '\n')
cat('Formats: PNG 600 dpi, vector PDF (Cairo), TIFF 600 dpi LZW\n')
