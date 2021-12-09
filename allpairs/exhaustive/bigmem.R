#!/usr/bin/env Rscript

library(ggplot2)
library(RColorBrewer)

fn <- 'bigmem-archfun-st.tbl'
dfs <- list()
do_stampede <- F
big_pmh <- T

dfs[['bigmem']] <- read.table(fn, comment.char='#', sep='\t')
columns <-  c('method', 'k', 'nreg', 'regsize', 'nthread',
                  'sketchtime', 'distancetime')
colnames(dfs[['bigmem']]) <- columns
dfs[['bigmem']]$totsize <- dfs[['bigmem']]$nreg * dfs[['bigmem']]$regsize

if(do_stampede) {
    dfs[['stampede2']] <- read.table('stampede2.tbl', comment.char='#')
    colnames(dfs[['stampede2']]) <- columns
    dfs[['stampede2']]$totsize <- dfs[['stampede2']]$nreg * dfs[['stampede2']]$regsize
}

speedplot <- function(df, k, totsize, regsize) {
    dftmp <- df[df$k == k & df$totsize == totsize & !grepl('txt', df$method),]
    #dftmp <- dftmp[!grepl('-2500000$', dftmp$method),]
    dftmp$method <- gsub('-bin-[1248]-500000$', '-500k', dftmp$method)
    dftmp$method <- gsub('-bin-[1248]-2500000$', '-2.5M', dftmp$method)
    dftmp$method <- gsub('-[1248]$', '', dftmp$method)
    dftmp$method <- gsub('-0.5$', '', dftmp$method)
    dftmp$method <- gsub('D2FSS-bin', 'D2-full', dftmp$method)
    dftmp$method <- gsub('D2OP-bin', 'D2', dftmp$method)
    dftmp$method <- gsub('^PMH', 'D2W', dftmp$method)
    colours <- c("Mash" = "green", "D2" = "red",
                 "D2W-500k" = "brown", "D2W-2.5M" = "black",
                 "Bindash" = "blue", "D2-full" = "purple")
    ggplot(dftmp,
        aes(x=sketchtime, y=distancetime,
            color=method, shape=factor(regsize))) +
        scale_color_manual(values = colours) +
        geom_point(size=3) + theme_bw() +
        labs(x='Sketch time (seconds)', y='Distance time (seconds)',
             title=paste0('Sketch size = ', totsize, ' bits')) +
        expand_limits(x=0, y=0)
}

for(machine in c('bigmem')) {
    for(k in c(31)) {
        for(totsz in c(8192, 16384, 32768)) {
              pdf(file=paste0(machine, '_k', k, '_sz', totsz, '.pdf'),
                width=8, height=4.5)
              print(speedplot(dfs[[machine]], k, totsz, regsz))
              dev.off()
        }
    }
}
