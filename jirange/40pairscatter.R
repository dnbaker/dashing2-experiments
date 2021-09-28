#!/usr/bin/env Rscript

library(fst)
library(tidyr)
library(dplyr)
library(ggplot2)

columns=c('totbits', 'k', 'ji_diff', 'ani_diff', 'type')
df <- read.fst('03filtered.fst', columns=columns)
df <- df %>% filter(!(type %in% c('BDN', 'BD1', 'BD2', 'BD4', 'BD8')))
stopifnot(length(unique(table(df$type))) == 1)

hist_type <- function(type1, type2, diff_col='ji_diff', est='JI') {
  columns=c('totbits', 'k', 'type', diff_col)
  df <- read.fst('03filtered.fst', columns=columns)
  d1 <- df[[diff_col]][df$type == type1]
  d2 <- df[[diff_col]][df$type == type2]
  k <- df$k[df$type == type1]
  bits <- df$totbits[df$type == type1]
  mx <- max(abs(c(d1, d2)))
  stopifnot(length(d1) == length(d2))
  newdf <- data.frame(diff1=d1, diff2=d2, k=k, bits=bits)
  ggplot(newdf, aes(x=d1, y=d2)) +
    scale_x_continuous(limits=c(-mx, mx)) +
    scale_y_continuous(limits=c(-mx, mx)) +
    facet_grid(bits ~ k) +
    labs(x=paste(type1, 'True', est, '- Est', est),
         y=paste(type2, 'True', est, '- Est', est)) +
    geom_point(col=rgb(0,0,.7,alpha=.2)) + geom_rug(col=rgb(0,0,0.7,alpha=.15)) +
       theme_bw()
}

pdf(file='40pairscatter_ss1_mash.pdf', width=7, height=14)
print(hist_type('SS1', 'Mash'))
dev.off()

pdf(file='40pairscatter_ss1_dash1.pdf', width=7, height=14)
print(hist_type('SS1', 'Dash1'))
dev.off()

pdf(file='40pairscatter_fss1_mash.pdf', width=7, height=14)
print(hist_type('FSS1', 'Mash'))
dev.off()

pdf(file='40pairscatter_fss1_dash1.pdf', width=7, height=14)
print(hist_type('FSS1', 'Dash1'))
dev.off()

pdf(file='40pairscatter_fss1_ss1.pdf', width=7, height=14)
print(hist_type('FSS1', 'SS1'))
dev.off()
