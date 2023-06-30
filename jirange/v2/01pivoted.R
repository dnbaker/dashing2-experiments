#!/usr/bin/env Rscript

library(fst)
library(tidyr)

df.wide <- read.fst('00wide.fst')
df.wide$JItrue <- df.wide$JI
df.wide$WJItrue <- df.wide$WJI
cols <- c(Mash, Dash1,
          BD8, BD4, BD2, BD1,
          SS8, SS4, SS2, SS1,
          FSS8, FSS4, FSS2, FSS1,
          JItrue, WJItrue,
          PMH8Exact, 'PMH8-50000000', 'PMH8-500000',
          PMH4Exact, 'PMH4-50000000', 'PMH4-500000',
          PMH2Exact, 'PMH2-50000000', 'PMH2-500000',
          PMH1Exact, 'PMH1-50000000', 'PMH1-500000')
df <- df.wide %>% pivot_longer(
    cols=cols, names_to = "type", values_to = "jest")
write.fst(df, '01pivoted.fst')
