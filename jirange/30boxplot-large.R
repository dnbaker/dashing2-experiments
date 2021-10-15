#!/usr/bin/env Rscript

# TODO: subset

library(fst)
library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)

for(do_subset in c(F, T)) {


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
      title_str <- if (do_subset) {'ANI >= 87%'} else {'All ANI'}
      dffilt %>% ggplot(aes(x=type, color=type, y=get(column))) +
        labs(x='Sketch type', y=paste0('Est ', est, ' - True ', est)) +
        ggtitle(paste0(title_str, ', bytes = ', bits/8, ', k = ', k)) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme_bw()
    }

    if(do_subset) {
        incl_types_ji <- c('Dash1', 'Mash', 'BD1', 'BD8', 'FSS1', 'SS1')
        incl_types_ani <- c('Dash1', 'Mash', 'BD1', 'BD8', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                            'PMH1-50000000', 'PMH1-500000',
                            'BMH1-50000000', 'BMH1-500000')
        incl_types_ani_lots <- c('Dash1', 'Mash', 'BD1', 'BD8', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                                 'PMH8-50000000', 'PMH4-50000000', 'PMH2-50000000', 'PMH1-50000000', 'PMHN-50000000',
                                 'PMH8-500000', 'PMH4-500000', 'PMH2-500000', 'PMH1-500000', 'PMHN-500000',
                                 'BMH8-50000000', 'BMH4-50000000', 'BMH2-50000000', 'BMH1-50000000', 'BMHN-50000000',
                                 'BMH8-500000', 'BMH4-500000', 'BMH2-500000', 'BMH1-500000', 'BMHN-500000')
    } else {
        incl_types_ji <- c('Dash1', 'Mash', 'FSS1', 'SS1')
        incl_types_ani <- c('Dash1', 'Mash', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                            'PMH1-50000000', 'PMH1-500000',
                            'BMH1-50000000', 'BMH1-500000')
        incl_types_ani_lots <- c('Dash1', 'Mash', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                                 'PMH8-50000000', 'PMH4-50000000', 'PMH2-50000000', 'PMH1-50000000', 'PMHN-50000000',
                                 'PMH8-500000', 'PMH4-500000', 'PMH2-500000', 'PMH1-500000', 'PMHN-500000',
                                 'BMH8-50000000', 'BMH4-50000000', 'BMH2-50000000', 'BMH1-50000000', 'BMHN-50000000',
                                 'BMH8-500000', 'BMH4-500000', 'BMH2-500000', 'BMH1-500000', 'BMHN-500000')
    }

    fn <- if(do_subset) { '07filtered_subset-large.fst' } else { '03filtered-large.fst' }

    columns=c('totbits', 'k', 'ji_diff', 'wji_diff', 'ani_diff', 'type', 'flat')
    df <- read.fst(fn, columns=columns)

    suffix <- if(do_subset) { '_subset.pdf' } else { '.pdf' }

    for(bits in unique(df$totbits)) {
      for(k in unique(df$k)) {
        df.tmp <- df %>% mutate(type = fct_relevel(type, incl_types_ji))

        pdf(file=paste0('30boxplot_ji_', k, '_', bits, suffix))
        print(error_boxplot(df.tmp, bits, k, include_types=incl_types_ji) + geom_boxplot())
        dev.off()

        pdf(file=paste0('30violin_ji_', k, '_', bits, suffix))
        print(error_boxplot(df.tmp, bits, k, include_types=incl_types_ji) + geom_violin())
        dev.off()

        df.tmp <- df %>% mutate(type = fct_relevel(type, incl_types_ani))

        #pdf(file=paste0('30boxplot_ani_', k, '_', bits, suffix))
        #print(error_boxplot(df.tmp, bits, k, include_types=incl_types_ani, column='ani_diff', est='ANI') + geom_boxplot())
        #dev.off()

        pdf(file=paste0('30violin_ani_', k, '_', bits, suffix))
        print(error_boxplot(df.tmp, bits, k, include_types=incl_types_ani, column='ani_diff', est='ANI') + geom_violin())
        dev.off()
      }
    }

}
