#!/usr/bin/env Rscript

library(data.table)
library(fst)

fns <- c('wanghashbatch/20.2021-09-15-20:53:26.256-10.tbl.xz',
         'wanghashbatch/27.2021-09-15-20:53:36.256-10.tbl.xz')
cmd <- '/opt/miniconda3/bin/xz -d -c'
df.wide <- fread(cmd = paste(cmd, fns[1]), sep='\t', header=T)
for(fn in fns[2:length(fns)]) {
  df.wide <- rbind(df.wide, fread(cmd = paste(cmd, fn), sep='\t', header=T))
}
write.fst(df.wide, '00wide.fst')

print('Number of rows where JI == 0:')
print(sum(df.wide$JI == 0))

print('Number of rows where ANI == 0:')
print(sum(df.wide$ANI == 0))

print('Number of rows where jest < 0:')
print(sum(df.wide$jest < 0))

print('Number of rows where jest > 1:')
print(sum(df.wide$jest > 1.0))
