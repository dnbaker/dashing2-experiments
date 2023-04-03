#!/usr/bin/env Rscript

library(fst)
library(tidyr)
library(dplyr)

for(do_subset in c(F, T)) {

    fn <- if(do_subset) { '07filtered_subset-large.fst' } else { '03filtered-large.fst' }

    se_df <- function(ani=F, bins=10,
                      exclude_types=c(),
                      include_types=c())
    {
        columns=c('totbits', 'k', 'ANI', 'JI', 'ani_diff', 'ji_diff', 'type')
        dfsse <- read.fst(fn, columns=columns) %>%
          mutate(mest=(if(!!ani) { ANI/100.0 } else { JI })) %>%
          mutate(bin=as.integer(bins * (mest*0.99999))/bins) %>%
          mutate(mdiff=(if(!!ani) { ani_diff } else { ji_diff })) %>%
          dplyr::group_by(k, totbits, bin, type) %>%
          dplyr::summarise(sse=sum(mdiff**2), se=sum(mdiff), sae=sum(abs(mdiff)), n=dplyr::n())
        if(length(include_types) > 0) {
            dfsse <- dfsse %>% filter(type %in% include_types)
        } else if(length(exclude_types) > 0) {
            dfsse <- dfsse %>% filter((!type %in% exclude_types))
        }
        dfsse
    }

    options("width"=200)

    if(do_subset) {
        incl_types_ji <- c('Dash1', 'Mash', 'BD1', 'BD2', 'BD4', 'BD8', 'FSS1', 'SS1', 'FSS8', 'SS8')
        incl_types_ani <- c('Dash1', 'Mash', 'BD1', 'BD2', 'BD4', 'BD8', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                            'PMH1-50000000', 'PMH1-500000',
                            'BMH1-50000000', 'BMH1-500000')
        incl_types_ani_lots <- c('Dash1', 'Mash', 'BD1', 'BD2', 'BD4', 'BD8', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                                 'PMH8-50000000', 'PMH4-50000000', 'PMH2-50000000', 'PMH1-50000000', 'PMHN-50000000',
                                 'PMH8-500000', 'PMH4-500000', 'PMH2-500000', 'PMH1-500000', 'PMHN-500000',
                                 'BMH8-50000000', 'BMH4-50000000', 'BMH2-50000000', 'BMH1-50000000', 'BMHN-50000000',
                                 'BMH8-500000', 'BMH4-500000', 'BMH2-500000', 'BMH1-500000', 'BMHN-500000')
    } else {
        incl_types_ji <- c('Dash1', 'Mash', 'FSS1', 'FSS8', 'SS1', 'SS8')
        incl_types_ani <- c('Dash1', 'Mash', 'FSSN', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                            'PMH1-50000000', 'PMH1-500000',
                            'BMH1-50000000', 'BMH1-500000', 'FSS8', 'SS8')
        incl_types_ani_lots <- c('Dash1', 'Mash', 'FSSN', 'FSS1', 'SS1', 'JItrue', 'WJItrue',
                                 'PMH8-50000000', 'PMH4-50000000', 'PMH2-50000000', 'PMH1-50000000', 'PMHN-50000000',
                                 'PMH8-500000', 'PMH4-500000', 'PMH2-500000', 'PMH1-500000', 'PMHN-500000',
                                 'BMH8-50000000', 'BMH4-50000000', 'BMH2-50000000', 'BMH1-50000000', 'BMHN-50000000',
                                 'BMH8-500000', 'BMH4-500000', 'BMH2-500000', 'BMH1-500000', 'BMHN-500000')
    }

    suffix <- if(do_subset) { '_subset.tbl' } else { '.tbl' }

    tb <- se_df(ani=T, bins=1, include_types=incl_types_ani) %>%
      select(-se, -sae, -bin) %>%
      pivot_wider(names_from='type', values_from='sse')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ani_sse', suffix))

    tb <- se_df(ani=T, bins=1, include_types=incl_types_ani) %>%
      select(-sse, -sae, -bin) %>%
      pivot_wider(names_from='type', values_from='se')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ani_se', suffix))

    tb <- se_df(ani=T, bins=1, include_types=incl_types_ani) %>%
      select(-se, -sse, -bin) %>%
      pivot_wider(names_from='type', values_from='sae')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ani_sae', suffix))

    tb <- se_df(ani=F, bins=1, include_types=incl_types_ji) %>%
      select(-se, -sae, -bin) %>%
      pivot_wider(names_from='type', values_from='sse')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ji_sse', suffix))

    tb <- se_df(ani=F, bins=1, include_types=incl_types_ji) %>%
      select(-sse, -sae, -bin) %>%
      pivot_wider(names_from='type', values_from='se')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ji_se', suffix))

    tb <- se_df(ani=F, bins=1, include_types=incl_types_ji) %>%
      select(-se, -sse, -bin) %>%
      pivot_wider(names_from='type', values_from='sae')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ji_sae', suffix))

    tb <- se_df(ani=T, bins=20, include_types=incl_types_ani) %>%
      select(-se, -sae) %>%
      pivot_wider(names_from='type', values_from='sse')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ani_bins_sse', suffix))

    tb <- se_df(ani=T, bins=20, include_types=incl_types_ani) %>%
      select(-sse, -sae) %>%
      pivot_wider(names_from='type', values_from='se')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ani_bins_se', suffix))

    tb <- se_df(ani=T, bins=20, include_types=incl_types_ani) %>%
      select(-se, -sse) %>%
      pivot_wider(names_from='type', values_from='sae')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ani_bins_sae', suffix))

    tb <- se_df(ani=F, bins=20, include_types=incl_types_ji) %>%
      select(-se, -sae) %>%
      pivot_wider(names_from='type', values_from='sse')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ji_bins_sse', suffix))

    tb <- se_df(ani=F, bins=20, include_types=incl_types_ji) %>%
      select(-sse, -sae) %>%
      pivot_wider(names_from='type', values_from='se')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ji_bins_se', suffix))

    tb <- se_df(ani=F, bins=20, include_types=incl_types_ji) %>%
      select(-se, -sse) %>%
      pivot_wider(names_from='type', values_from='sae')
    capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_large_ji_bins_sae', suffix))
}
