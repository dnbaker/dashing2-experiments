#!/usr/bin/env Rscript

library(fst)

df <- read.fst('02augmented.fst')
#dffilt <- df[df$totbits >= 16384 & df$totbits <= 524288,]
dffilt <- df
write.fst(dffilt, '03filtered.fst')
