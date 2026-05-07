need <- c('lcmm','flexmix','poLCA','lavaan','semTools','MplusAutomation','haven',
          'dplyr','tidyr','survey','tidyLPA','bootnet','qgraph',
          'NetworkComparisonTest','mgm','mclust','mixtools','psych','semtree','rpart')
inst <- rownames(installed.packages())
have <- need[need %in% inst]
miss <- need[!need %in% inst]
cat('Installed:\n'); cat(have, sep='\n'); cat('\n\nMissing:\n'); cat(miss, sep='\n')
