#!/usr/bin/env Rscript

library(fst)

df <- read.fst('06augmented_subset-large.fst')
dffilt <- df[df$totbits >= 16384 & df$totbits <= 524288,]
write.fst(dffilt, '07filtered_subset-large.fst')
