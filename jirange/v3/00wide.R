#!/usr/bin/env Rscript

library(data.table)
library(fst)

fns <- c('20230325_k21_008kbit_out.txt',
         '20230325_k21_032kbit_out.txt',
         '20230325_k21_128kbit_out.txt',
         '20230325_k31_008kbit_out.txt',
         '20230325_k31_032kbit_out.txt',
         '20230325_k31_128kbit_out.txt')

df.wide <- fread(fns[1], sep='\t', header=T)
if(length(fns) > 1) {
    for(fn in fns[2:length(fns)]) {
      df.wide <- rbind(df.wide, fread(fn, sep='\t', header=T))
    }
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
