#!/usr/bin/env Rscript

library(fst)
library(tidyr)

df.wide <- read.fst('00wide.fst')
df.wide$JItrue <- df.wide$JI
df.wide$WJItrue <- df.wide$WJI
cols <- c('Mash', 'Dash1',
          'BD8', 'BD1',
          'SS8', 'SS1',
          'FSS8', 'FSS1',
          'JItrue', 'WJItrue',
          'SM',
          'PMH8Exact', 'PMH8-50000000',
          'PMH1Exact', 'PMH1-50000000')
df <- df.wide %>% pivot_longer(
    cols=all_of(cols), names_to = "type", values_to = "jest")
write.fst(df, '01pivoted.fst')
