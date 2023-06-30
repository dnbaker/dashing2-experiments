#!/usr/bin/env Rscript

library(fst)
library(dplyr)

df <- read.fst('01pivoted.fst')
df$k <- df$K
df <- df[df$ANI >= 89,]

print('Invalid high JI')
print('===============')
print('Histogram of how often each sketch type has high JI:')
table(df$type[df$jest > 1.0])
print('Histogram of values of high JIs:')
table(df$jest[df$jest > 1.0])
df$jest <- pmin(df$jest, 1.0)

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

#r low_and_high_hi
print('Invalid low JI')
print('==============')
print('Histogram of how often each sketch type has low JI:')
table(df$type[df$jest < 0.0])
print('Histogram of values of low JIs:')
table(df$jest[df$jest < 0.0])
df <- df[df$jest >= 0.0,]

stopifnot(all(0 <= df$jest))
stopifnot(all(df$jest <= 1.0))
stopifnot(all(0 <= df$ani_est))
stopifnot(all(df$ani_est <= 1.0))

df$ani_quartile <- ntile(df$ANI, 4)

write.fst(df, '06augmented_subset.fst')
