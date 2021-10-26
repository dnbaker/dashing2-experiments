#!/usr/bin/env Rscript

library(fst)

df <- read.fst('02augmented-large.fst')
dffilt <- df[df$totbits >= 16384 & df$totbits <= 524288,]
write.fst(dffilt, '03filtered-large.fst')
