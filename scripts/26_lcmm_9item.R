# R lcmm replication of 9-item primary LCGA
suppressPackageStartupMessages({
  library(lcmm); library(haven); library(dplyr); library(tidyr)
})

ROOT <- 'D:/2026/SCI/청소년패널조사'
OUT  <- file.path(ROOT, 'output', 'paperA_lcga', 'primary_9item', 'lcmm_R')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Build 9-item depression long-form
items_9_idx <- c(1,2,3,5,6,7,8,9,10)
long <- list()
for (w in 1:7) {
  d <- read_dta(file.path(ROOT, sprintf('KCYPS2018m1Yw%d.dta', w)))
  items <- sprintf('YPSY4E%02dw%d', items_9_idx, w)
  sub <- as.data.frame(d[, items])
  sub[sub < 1 | sub > 4] <- NA
  n_obs <- rowSums(!is.na(sub))
  dep9 <- ifelse(n_obs >= 5, rowMeans(sub, na.rm = TRUE), NA_real_)
  long[[w]] <- data.frame(ID = d$ID, wave = w, dep9 = dep9)
}
ldf <- do.call(rbind, long)
ldf <- ldf[!is.na(ldf$dep9), ]
keep_id <- names(which(table(ldf$ID) >= 3))
ldf <- ldf[as.character(ldf$ID) %in% keep_id, ]
ldf$ID <- as.numeric(ldf$ID)
ldf$time <- ldf$wave - 4
cat(sprintf('9-item long form: %d rows, %d unique IDs\n', nrow(ldf), length(unique(ldf$ID))))

# K=1
mod1 <- hlme(dep9 ~ time + I(time^2), random = ~ -1,
             subject = 'ID', ng = 1, data = ldf, verbose = FALSE)

# K=2..5 with gridsearch
fit_K <- function(K) {
  gridsearch(
    hlme(dep9 ~ time + I(time^2),
         mixture = ~ time + I(time^2),
         random = ~ -1,
         subject = 'ID', ng = K, data = ldf, verbose = FALSE),
    rep = 30, maxiter = 30, minit = mod1
  )
}

cat('Fitting K=2..5 ...\n')
mod2 <- fit_K(2); cat('K=2 done\n')
mod3 <- fit_K(3); cat('K=3 done\n')
mod4 <- fit_K(4); cat('K=4 done\n')
mod5 <- fit_K(5); cat('K=5 done\n')

saveRDS(list(mod1, mod2, mod3, mod4, mod5), file.path(OUT, 'lcga_models_9item.rds'))

# Fit summary
fit_tab <- summarytable(mod1, mod2, mod3, mod4, mod5,
                        which = c('G','loglik','npm','AIC','BIC','SABIC',
                                  'entropy','%class'))
write.csv(fit_tab, file.path(OUT, 'lcmm_fit_summary_9item.csv'))
cat('\n=== R lcmm fit summary (9-item) ===\n'); print(fit_tab)

# Class predictions for K=4
newdf <- data.frame(time = -3:3)
pred4 <- predictY(mod4, newdata = newdf, var.time = 'time')
p4_df <- as.data.frame(pred4$pred); p4_df$wave <- newdf$time + 4
write.csv(p4_df, file.path(OUT, 'class_pred_K4_9item.csv'), row.names = FALSE)

pp4 <- mod4$pprob
write.csv(pp4, file.path(OUT, 'posterior_K4_9item.csv'), row.names = FALSE)

cat('\n=== R lcmm K=4 9-item proportions ===\n')
print(prop.table(table(pp4$class)))
cat('\n=== R lcmm K=4 9-item predicted means by class and wave ===\n')
print(p4_df)

cat('\nDone — R lcmm 9-item outputs in', OUT, '\n')
