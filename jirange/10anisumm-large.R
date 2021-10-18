#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(fst)

dfsse <- read.fst('04sse-large.fst')

plot_ani_sum <- function(bits=bits, k=k, max_error=Inf, squared=T, types=NULL, ani_min=0.7) {
    if(squared) {
      bad_types <- dfsse %>% dplyr::group_by(type) %>%
        dplyr::summarise(max_error=max(sse)) %>% filter(max_error > !!max_error)
    } else {
      bad_types <- dfsse %>% dplyr::group_by(type) %>%
        dplyr::summarise(max_error=max(se)) %>% filter(max_error > !!max_error)
    }

    dfssetmp <- dfsse %>% dplyr::filter(totbits == !!bits & k == !!k & !(type %in% bad_types) & anibin >= ani_min)
    if(!is.null(types)) {
        dfssetmp <- dfssetmp %>% filter(type %in% types)
    }
    if(squared) {
        pl <- ggplot(dfssetmp, aes(x=anibin, color=type, y=sse))
    } else {
        pl <- ggplot(dfssetmp, aes(x=anibin, color=type, y=se))
    }
    pl + geom_point(position=position_dodge(width = 0.01)) +
      labs(title=paste('Sketch total of', as.character(bits), 'bits, k =', k),
           y=paste0(ifelse(squared, 'SSE', 'SE'), '(Est ANI - True ANI)')) +
      theme_bw()
}

pdf(file='10largeanisumm_16384_21.pdf')
plot_ani_sum(bits=16384, k=21, squared=F)
dev.off()

pdf(file='10largeanisumm_65536_21.pdf')
plot_ani_sum(bits=65536, k=21, squared=F)
dev.off()

pdf(file='10largeanisumm_262144_21.pdf')
plot_ani_sum(bits=262144, k=21, squared=F)
dev.off()

pdf(file='10largeanisumm_1048576_21.pdf')
plot_ani_sum(bits=1048576, k=21, squared=F)
dev.off()

pdf(file='10largeanisumm_16384_31.pdf')
plot_ani_sum(bits=16384, k=31, squared=F)
dev.off()

pdf(file='10largeanisumm_65536_31.pdf')
plot_ani_sum(bits=65536, k=31, squared=F)
dev.off()

pdf(file='10largeanisumm_262144_31.pdf')
plot_ani_sum(bits=262144, k=31, squared=F)
dev.off()

pdf(file='10largeanisumm_1048576_31.pdf')
plot_ani_sum(bits=1048576, k=31, squared=F)
dev.off()


pdf(file='10largeanisumm_16384_71.pdf')
plot_ani_sum(bits=16384, k=71, squared=F)
dev.off()

pdf(file='10largeanisumm_65536_71.pdf')
plot_ani_sum(bits=65536, k=71, squared=F)
dev.off()

pdf(file='10largeanisumm_262144_71.pdf')
plot_ani_sum(bits=262144, k=71, squared=F)
dev.off()

pdf(file='10largeanisumm_1048576_71.pdf')
plot_ani_sum(bits=1048576, k=71, squared=F)
dev.off()

pdf(file='10largeanisumm_squared_16384_21.pdf')
plot_ani_sum(bits=16384, k=21, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_65536_21.pdf')
plot_ani_sum(bits=65536, k=21, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_262144_21.pdf')
plot_ani_sum(bits=262144, k=21, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_1048576_21.pdf')
plot_ani_sum(bits=1048576, k=21, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_16384_31.pdf')
plot_ani_sum(bits=16384, k=31, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_65536_31.pdf')
plot_ani_sum(bits=65536, k=31, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_262144_31.pdf')
plot_ani_sum(bits=262144, k=31, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_1048576_31.pdf')
plot_ani_sum(bits=1048576, k=31, squared=T)
dev.off()


pdf(file='10largeanisumm_squared_16384_71.pdf')
plot_ani_sum(bits=16384, k=71, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_65536_71.pdf')
plot_ani_sum(bits=65536, k=71, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_262144_71.pdf')
plot_ani_sum(bits=262144, k=71, squared=T)
dev.off()

pdf(file='10largeanisumm_squared_1048576_71.pdf')
plot_ani_sum(bits=1048576, k=71, squared=T)
dev.off()
