options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- c('lcmm','NetworkComparisonTest','mixtools')
miss <- need[!need %in% rownames(installed.packages())]
if (length(miss)) {
  install.packages(miss, dependencies = TRUE, quiet = TRUE)
}
inst <- rownames(installed.packages())
cat('After install — present:', need[need %in% inst], sep='\n')
cat('\nMissing:', need[!need %in% inst], sep='\n')
