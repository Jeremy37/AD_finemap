#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
pCol <- as.integer(args[2])
header <- args[3]

if (is.na(header) | as.logical(header) == F) {
  header = F
} else {
  header = T
}
df <- read.table(infile, header=header, row.names=NULL)

df$fdr <- p.adjust(df[,pCol], method="fdr")

write.table(df, file="", sep="\t", col.names=header, row.names=F, quote=F)
