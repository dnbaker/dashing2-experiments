#!/usr/bin/env Rscript

library(fst)
library(dplyr)

columns=c('totbits', 'k', 'ANI', 'ani_diff', 'type')
dfsse <- read.fst('07filtered_subset-large.fst', columns=columns) %>%
  mutate(anibin=round(ANI/100.0, digits=2)) %>%
  dplyr::group_by(k, totbits, anibin, type) %>%
  dplyr::summarise(sse=sum((ani_diff)**2), se=sum(ani_diff))

write.fst(dfsse, '08sse_subset-large.fst')
