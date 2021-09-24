#!/usr/bin/env Rscript

library(fst)
library(tidyr)
library(dplyr)
library(ggplot2)

error_boxplot <- function(df, bits, k, column='ji_diff', est='JI',
                          exclude_types=c(), include_types=c())
{
  dffilt <- df %>% filter((totbits == !!bits) & (k == !!k))
  if(length(include_types) > 0) {
      dffilt <- dffilt %>% filter(type %in% include_types)
  } else if(length(exclude_types) > 0) {
      dffilt <- dffilt %>% filter(!(type %in% exclude_types))
  }
  stopifnot(nrow(dffilt) > 0)
  dffilt %>% ggplot(aes(x=type, color=type, y=get(column))) +
    labs(x='Sketch type', y=paste0('Est ', est, ' - True ', est)) +
    ggtitle(paste0('Bytes = ', bits/8, ', k = ', k)) +
    theme_bw()
}

incl_types_ji <- c('Dash1', 'Mash', 'FSSN', 'FSS1', 'SS1')
incl_types_ani <- c('Dash1', 'Mash', 'FSSN', 'FSS1', 'SS1',
                    'PMH8-50000000', 'PMH1-50000000',
                    'BMH8-50000000', 'BMH1-50000000')

columns=c('totbits', 'k', 'ji_diff', 'wji_diff', 'ani_diff', 'type', 'flat')
df <- read.fst('03filtered.fst', columns=columns)

for(bits in unique(df$totbits)) {
  for(k in unique(df$k)) {
    pdf(file=paste0('30boxplot_ji_', k, '_', bits, '.pdf'))
    print(error_boxplot(df, bits, k, include_types=incl_types_ji) + geom_boxplot())
    dev.off()

    pdf(file=paste0('30violin_ji_', k, '_', bits, '.pdf'))
    print(error_boxplot(df, bits, k, include_types=incl_types_ji) + geom_violin())
    dev.off()

    pdf(file=paste0('30boxplot_ani_', k, '_', bits, '.pdf'))
    print(error_boxplot(df, bits, k, include_types=incl_types_ani, column='ani_diff', est='ANI') + geom_boxplot())
    dev.off()

    pdf(file=paste0('30violin_ani_', k, '_', bits, '.pdf'))
    print(error_boxplot(df, bits, k, include_types=incl_types_ani, column='ani_diff', est='ANI') + geom_violin())
    dev.off()
  }
}
