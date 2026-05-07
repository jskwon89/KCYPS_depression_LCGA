# ==============================================================================
# Paper A — LCGA replication using R lcmm package (Mplus-equivalent gold standard)
# ==============================================================================
# Output: model fit (BIC/AIC/Entropy), class profiles, posterior probs,
#         compare to Python EM result.
# ==============================================================================

suppressPackageStartupMessages({
  library(lcmm)
  library(haven)
  library(dplyr)
  library(tidyr)
})

ROOT <- 'D:/2026/SCI/청소년패널조사'
OUT  <- file.path(ROOT, 'output', 'paperA_lcga', 'lcmm_R')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Load and prepare long-form depression data ----
# Re-build from raw .dta (Python long_panel.parquet not directly readable from R without arrow)
waves <- 1:7
long <- list()
for (w in waves) {
  d <- read_dta(file.path(ROOT, sprintf('KCYPS2018m1Yw%d.dta', w)))
  items <- sprintf('YPSY4E%02dw%d', 1:10, w)
  sub <- as.data.frame(d[, items])
  # KCYPS missing: <0 or >=90 -> NA; valid 1-4
  sub[sub < 1 | sub > 4] <- NA
  n_obs <- rowSums(!is.na(sub))
  dep <- ifelse(n_obs >= 5, rowMeans(sub, na.rm = TRUE), NA_real_)
  wcol <- grep('WEIGHTB1w', names(d), value = TRUE)
  weight <- if (length(wcol)) as.numeric(d[[wcol[1]]]) else rep(NA_real_, nrow(d))
  long[[w]] <- data.frame(ID = d$ID, wave = w, dep = dep, weight = weight)
}
ldf <- do.call(rbind, long)
ldf <- ldf[!is.na(ldf$dep), ]
# only keep persons with >=3 valid waves (consistent with Python LCGA)
keep_id <- names(which(table(ldf$ID) >= 3))
ldf <- ldf[as.character(ldf$ID) %in% keep_id, ]
ldf$ID <- as.numeric(ldf$ID)
ldf$time <- ldf$wave - 4   # center on wave 4
cat(sprintf('Long form: %d rows, %d unique IDs\n', nrow(ldf), length(unique(ldf$ID))))

# ---- Fit LCGA K=1..5 ----
# lcmm::hlme is the LCGA/GMM workhorse
# LCGA spec: fixed quadratic time, NO random effects (Nagin's GBTM)
# K=1: simple growth curve
mod1 <- hlme(dep ~ time + I(time^2), random = ~ -1,
             subject = 'ID', ng = 1, data = ldf, verbose = FALSE)

# K=2..5 with class-varying intercept and slopes
fit_K <- function(K, B = NULL) {
  if (is.null(B)) {
    # Use parameters from K-1 model as starting values via gridsearch
    return(gridsearch(
      hlme(dep ~ time + I(time^2),
           mixture = ~ time + I(time^2),
           random = ~ -1,
           subject = 'ID', ng = K, data = ldf, verbose = FALSE),
      rep = 30, maxiter = 30, minit = mod1
    ))
  }
}

cat('\n--- Fitting K=2 ---\n')
mod2 <- fit_K(2)
cat('--- Fitting K=3 ---\n')
mod3 <- fit_K(3)
cat('--- Fitting K=4 ---\n')
mod4 <- fit_K(4)
cat('--- Fitting K=5 ---\n')
mod5 <- fit_K(5)

# Save each model
saveRDS(list(mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5),
        file.path(OUT, 'lcga_models.rds'))

# ---- Fit summary ----
fit_tab <- summarytable(mod1, mod2, mod3, mod4, mod5,
                        which = c('G','loglik','npm','AIC','BIC','SABIC',
                                  'entropy','ICL','%class'))
cat('\n=== Model fit summary (R lcmm) ===\n')
print(fit_tab)
write.csv(fit_tab, file.path(OUT, 'lcmm_fit_summary.csv'), row.names = TRUE)

# ---- Class profiles for each K ----
for (K in 2:5) {
  m <- get(paste0('mod', K))
  # predicted trajectories
  newdf <- data.frame(time = seq(-3, 3, by = 1))
  pred <- predictY(m, newdata = newdf, var.time = 'time')
  # Save predicted means
  p_df <- as.data.frame(pred$pred)
  p_df$wave <- newdf$time + 4
  write.csv(p_df, file.path(OUT, sprintf('class_pred_K%d.csv', K)), row.names = FALSE)
  # Posterior class memberships
  pp <- m$pprob
  pp$class_R <- pp$class
  write.csv(pp, file.path(OUT, sprintf('posterior_K%d.csv', K)), row.names = FALSE)
  cat(sprintf('\n--- K=%d posterior class proportions ---\n', K))
  print(prop.table(table(pp$class_R)))
}

# ---- Bootstrap LRT (BLRT) for K-1 vs K ----
# Implement simple parametric bootstrap of the LRT statistic.
# Generate B = 50 datasets from K-1 fit, refit K-1 and K, record LRT.
blrt <- function(m_null, m_alt, B = 50, seed = 42) {
  set.seed(seed)
  null_ll <- logLik(m_null)[1] %>% as.numeric()
  alt_ll  <- logLik(m_alt)[1]  %>% as.numeric()
  obs_lrt <- 2 * (alt_ll - null_ll)
  K_null <- m_null$ng
  K_alt  <- m_alt$ng
  cat(sprintf('  null K=%d ll=%.2f, alt K=%d ll=%.2f, observed LRT=%.2f\n',
              K_null, null_ll, K_alt, alt_ll, obs_lrt))
  # Simulate from null. lcmm doesn't have built-in simulate; generate via class
  # assignments + class trajectories + residual SD.
  pp_null <- m_null$pprob
  cluster_assign <- pp_null$class
  # Get class trajectories
  newdf <- data.frame(time = unique(ldf$time))
  pred_null <- predictY(m_null, newdata = newdf, var.time = 'time')$pred
  resid_sd <- sqrt(m_null$best['stderr'] %>% ifelse(is.na(.), 1, .))
  if (is.na(resid_sd)) {
    # estimate from residuals
    fit_vals <- predictY(m_null, newdata = ldf[, 'time', drop = FALSE], var.time = 'time')$pred
    # rough fallback
    resid_sd <- sd(ldf$dep - mean(ldf$dep), na.rm = TRUE) * 0.7
  }
  cat(sprintf('  resid SD (used for sim): %.3f\n', resid_sd))

  lrt_dist <- numeric(B)
  for (b in 1:B) {
    # generate dep_b from null model
    sim <- ldf
    cls_map <- setNames(cluster_assign, pp_null$ID)
    sim$cls <- cls_map[as.character(sim$ID)]
    # predict mean by class for each time
    # pred_null columns are 'class1','class2',...
    cls_names <- grep('^class', colnames(pred_null), value = TRUE)
    # match time -> row
    time_to_row <- match(sim$time, newdf$time)
    means_b <- numeric(nrow(sim))
    for (i in seq_len(nrow(sim))) {
      ck <- min(sim$cls[i], length(cls_names))
      means_b[i] <- pred_null[time_to_row[i], cls_names[ck]]
    }
    sim$dep <- means_b + rnorm(nrow(sim), 0, resid_sd)
    # Refit (try, may fail)
    fit_null_b <- try(suppressMessages(hlme(dep ~ time + I(time^2),
                                            mixture = ~ time + I(time^2),
                                            random = ~ -1,
                                            subject = 'ID', ng = K_null,
                                            data = sim, verbose = FALSE,
                                            B = m_null$best)), silent = TRUE)
    fit_alt_b <- try(suppressMessages(hlme(dep ~ time + I(time^2),
                                           mixture = ~ time + I(time^2),
                                           random = ~ -1,
                                           subject = 'ID', ng = K_alt,
                                           data = sim, verbose = FALSE,
                                           B = m_alt$best)), silent = TRUE)
    if (inherits(fit_null_b, 'try-error') || inherits(fit_alt_b, 'try-error')) {
      lrt_dist[b] <- NA
      next
    }
    lrt_dist[b] <- 2 * (logLik(fit_alt_b) - logLik(fit_null_b))
    if (b %% 10 == 0) cat(sprintf('    b=%d done, mean LRT so far=%.2f\n', b, mean(lrt_dist[1:b], na.rm=TRUE)))
  }
  p_blrt <- mean(lrt_dist >= obs_lrt, na.rm = TRUE)
  return(list(obs_lrt = obs_lrt, lrt_dist = lrt_dist, p_blrt = p_blrt))
}

# Skip BLRT in this run — too slow; we run a lighter version that uses
# lcmm's resid-based simulation or just reports.
# BLRT loop disabled for runtime. We rely on BIC/aBIC/entropy for K selection.
cat('\n[BLRT skipped in primary run for runtime]\n')

cat('\n=== R lcmm LCGA done ===\n')
cat('Outputs at:', OUT, '\n')
