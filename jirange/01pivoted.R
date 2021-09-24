#!/usr/bin/env Rscript

library(fst)
library(tidyr)

df.wide <- read.fst('00wide.fst')
df <- df.wide %>% pivot_longer(cols=c(Mash, Dash1,
                                 BD8, BD4, BD2, BD1, BDN,
                                 SS8, SS4, SS2, SS1, SSN,
                                 FSS8, FSS4, FSS2, FSS1, FSSN,
                                 MH8, MH4, MH2, MH1, MHN,
                                 FMH8, FMH4, FMH2, FMH1, FMHN,
                                 PMH8Exact, 'PMH8-50000000', 'PMH8-500000',
                                 PMH4Exact, 'PMH4-50000000', 'PMH4-500000',
                                 PMH2Exact, 'PMH2-50000000', 'PMH2-500000',
                                 PMH1Exact, 'PMH1-50000000', 'PMH1-500000',
                                 PMHNExact, 'PMHN-50000000', 'PMHN-500000',
                                 BMH8Exact, 'BMH8-50000000', 'BMH8-500000',
                                 BMH4Exact, 'BMH4-50000000', 'BMH4-500000',
                                 BMH2Exact, 'BMH2-50000000', 'BMH2-500000',
                                 BMH1Exact, 'BMH1-50000000', 'BMH1-500000',
                                 BMHNExact, 'BMHN-50000000', 'BMHN-500000'),
                          names_to = "type", values_to = "jest")
write.fst(df, '01pivoted.fst')
