#!/usr/bin/env Rscript

library(data.table)
library(fst)

fns <- commandArgs(trailingOnly = TRUE)
getcmd <- function(name) {
    return(ifelse(endsWith(name, '.xz'), 'xz -dc', ifelse(endsWith(name, '.gz'), 'gzip -dc', 'cat')))
}
cmd <- getcmd(fns[1])
df.wide <- fread(cmd = paste(cmd, fns[1]), sep='\t', header=T)
if(length(fns) > 1) {
  for(fn in fns[2:length(fns)]) {
    cmd <- getcmd(fn)
    df.wide <- rbind(df.wide, fread(cmd = paste(cmd, fn), sep='\t', header=T))
  }
}
write.fst(df.wide, '00wide-large.fst')

print('Number of rows where JI == 0:')
print(sum(df.wide$JI == 0))

print('Number of rows where ANI == 0:')
print(sum(df.wide$ANI == 0))

print('Number of rows where jest < 0:')
print(sum(df.wide$jest < 0))

print('Number of rows where jest > 1:')
print(sum(df.wide$jest > 1.0))
