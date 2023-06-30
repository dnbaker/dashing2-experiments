#!/usr/bin/env Rscript

library(fst)
library(dplyr)

columns=c('sketchsize', 'k', 'ANI', 'ani_diff', 'type')
dfsse <- read.fst('02augmented.fst', columns=columns) %>%
  mutate(anibin=round(ANI/100.0, digits=2)) %>%
  filter(ANI >= 80.0) %>%
  dplyr::group_by(k, sketchsize, anibin, type) %>%
  dplyr::summarise(sse=sum((ani_diff)**2), se=sum(ani_diff))

write.fst(dfsse, '04sse.fst')
