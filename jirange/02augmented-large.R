#!/usr/bin/env Rscript

library(fst)
library(dplyr)

df <- read.fst('01pivoted-large.fst')
df$k <- df$K
df <- df[df$ANI > 0,]
df$totbits <- df$sketchsize
# some must be corrected (Mash for k > 32)
df$jest <- ifelse(df$k > 32 & (df$type == 'Mash'), -1., df$jest)
# some are negative
df$ani_est <- ifelse(df$jest > 0, 1 + 1/df$k * log(2*df$jest/(1+df$jest)), 0)
df$ani_est <- pmax(df$ani_est, 0.0)
stopifnot(all(!is.nan(df$ani_est)))
stopifnot(all(df$ani_est >= 0.0))
df$ji_diff <- df$jest - df$JI
df$wji_diff <- df$jest - df$WJI
stopifnot(all(0 < df$ANI))
stopifnot(all(df$ANI <= 100.0))

# r ani_diff_and_badness
df$ani_diff <- df$ani_est - (df$ANI/100.0)
df$se <- df$ani_diff ** 2
df$ani_badness <- factor(ifelse(df$jest < 0, 'J est = -1',
                          ifelse(df$jest > 1, 'J est = Inf',
                            ifelse(df$se < 0.08, 'SE<.08',
                              ifelse(df$se < 0.4, '.08<=SE<.4',
                                ifelse(df$se < 0.8, '.4<=SE<.8',
                                    '.8<=SE<=1'))))))

#r add_totbits
flat8 <- c('BD8', 'SS8', 'FSS8', 'MH8', 'FMH8', 'Mash')
flat4 <- c('BD4', 'SS4', 'FSS4', 'MH4', 'FMH4')
flat2 <- c('BD2', 'SS2', 'FSS2', 'MH2', 'FMH2')
flat1 <- c('BD1', 'SS1', 'FSS1', 'MH1', 'FMH1', 'Dash1', 'JItrue', 'WJItrue')
flatn <- c('BDN', 'FSSN', 'MHN', 'FMHN', 'SSN')
type_subsets <- list('n'=flatn, '1'=flat1, '2'=flat2, '4'=flat4, '8'=flat8)
r64 <- df$type %in% c(flat8, 'PMH8Exact', 'PMH8-50000000', 'PMH8-500000', 'BMH8Exact', 'BMH8-50000000', 'BMH8-500000')
r32 <- df$type %in% c(flat4, 'PMH4Exact', 'PMH4-50000000', 'PMH4-500000', 'BMH4Exact', 'BMH4-50000000', 'BMH4-500000')
r16 <- df$type %in% c(flat2, 'PMH2Exact', 'PMH2-50000000', 'PMH2-500000', 'BMH2Exact', 'BMH2-50000000', 'BMH2-500000')
r8  <- df$type %in% c(flat1, 'PMH1Exact', 'PMH1-50000000', 'PMH1-500000', 'BMH1Exact', 'BMH1-50000000', 'BMH1-500000')
r4  <- df$type %in% c(flatn, 'PMHNExact', 'PMHN-50000000', 'PMHN-500000', 'BMHNExact', 'BMHN-50000000', 'BMHN-500000')
df$flat <- ifelse(df$type %in% c(flat8, flat4, flat2, flat1, flatn, 'JItrue'), 1, 0)

df$totbits[r64] <- df$totbits[r64] * 64
df$totbits[r32] <- df$totbits[r32] * 32
df$totbits[r16] <- df$totbits[r16] * 16
df$totbits[r8] <- df$totbits[r8] * 8
df$totbits[r4] <- df$totbits[r4] * 4
stopifnot(all(df$totbits != df$sketchsize))

#r low_and_high_hi
print('Invalid low JI')
print('==============')
print('Histogram of how often each sketch type has low JI:')
table(df$type[df$jest < 0.0])
print('Histogram of values of low JIs:')
table(df$jest[df$jest < 0.0])
df <- df[df$jest >= 0.0,]

print('Invalid high JI')
print('===============')
print('Histogram of how often each sketch type has high JI:')
table(df$type[df$jest > 1.0])
print('Histogram of values of high JIs:')
table(df$jest[df$jest > 1.0])
df <- df[df$jest <= 1.0,]

stopifnot(all(0 <= df$jest))
stopifnot(all(df$jest <= 1.0))
stopifnot(all(0 <= df$ani_est))
stopifnot(all(df$ani_est <= 1.0))

df$ani_quartile <- ntile(df$ANI, 4)

write.fst(df, '02augmented-large.fst')
