#!/usr/bin/env Rscript

library(fst)

df <- read.fst('06augmented_subset.fst')
#dffilt <- df[df$totbits >= 16384 & df$totbits <= 524288,]
dffilt <- df
write.fst(dffilt, '07filtered_subset.fst')
