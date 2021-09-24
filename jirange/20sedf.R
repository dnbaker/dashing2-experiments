#!/usr/bin/env Rscript

library(fst)
library(tidyr)
library(dplyr)

se_df <- function(flat=T, ani=F,
                  bins=10,
                  squared=T,
                  exclude_types=c(),
                  include_types=c())
{
    columns=c('totbits', 'k', 'ANI', 'JI', 'ani_diff', 'ji_diff', 'type', 'flat')
    dfsse <- read.fst('03filtered.fst', columns=columns) %>%
      filter(flat == !!flat) %>%
      mutate(mest=(if(!!ani) { ANI/100.0 } else {1-(if(!!flat) {JI} else {WJI})})) %>%
      mutate(bin=as.integer(bins * (mest*0.99999))/bins) %>%
      mutate(mdiff=(if(!!ani) { ani_diff } else { ji_diff })) %>%
      dplyr::group_by(k, totbits, bin, type, flat) %>%
      dplyr::summarise(sse=sum(mdiff**2), se=sum(mdiff), n=dplyr::n())
    if(length(include_types) > 0) {
        dfssetmp <- dfsse %>% filter(type %in% include_types)
    } else if(length(exclude_types) > 0) {
        dfssetmp <- dfsse %>% filter((!type %in% exclude_types))
    }
    dfssetmp
}

options("width"=200)

incl_types <- c('Dash1', 'Mash', 'FSSN', 'FSS1', 'SSN', 'SS1')

tb <- se_df(ani=T, bins=1, include_types=incl_types, squared=T) %>%
  select(-se, -flat, -bin) %>%
  pivot_wider(names_from='type', values_from='sse')
capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_ani_sse.tbl'))

tb <- se_df(ani=T, bins=1, include_types=incl_types, squared=F) %>%
  select(-sse, -flat, -bin) %>%
  pivot_wider(names_from='type', values_from='se')
capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_ani_se.tbl'))

tb <- se_df(ani=T, bins=20, include_types=incl_types, squared=T) %>%
  select(-se, -flat) %>%
  pivot_wider(names_from='type', values_from='sse')
capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_ani_bins_sse.tbl'))

tb <- se_df(ani=T, bins=20, include_types=incl_types, squared=F) %>%
  select(-sse, -flat) %>%
  pivot_wider(names_from='type', values_from='se')
capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_ani_bins_se.tbl'))

tb <- se_df(ani=F, bins=20, include_types=incl_types, squared=T) %>%
  select(-se, -flat) %>%
  pivot_wider(names_from='type', values_from='sse')
capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_ji_bins_sse.tbl'))

tb <- se_df(ani=F, bins=20, include_types=incl_types, squared=F) %>%
  select(-sse, -flat) %>%
  pivot_wider(names_from='type', values_from='se')
capture.output(print(tb, n=nrow(tb), width=200), file=paste0('20sedf_ji_bins_se.tbl'))
