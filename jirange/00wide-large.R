#!/usr/bin/env Rscript

library(data.table)
library(fst)

fns <- c('outfiles.cat.tbl.xz')
cmd <- 'xz -d -c'
df.wide <- fread(cmd = paste(cmd, fns[1]), sep='\t', header=T)
write.fst(df.wide, '00wide-large.fst')

print('Number of rows where JI == 0:')
print(sum(df.wide$JI == 0))

print('Number of rows where ANI == 0:')
print(sum(df.wide$ANI == 0))

print('Number of rows where jest < 0:')
print(sum(df.wide$jest < 0))

print('Number of rows where jest > 1:')
print(sum(df.wide$jest > 1.0))
