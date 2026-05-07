need <- c('ggplot2','dplyr','tidyr','readr','scales','patchwork','ggsci','cowplot')
inst <- rownames(installed.packages())
cat('Have:\n'); cat(need[need %in% inst], sep='\n')
cat('\nMissing:\n'); cat(need[!need %in% inst], sep='\n')
