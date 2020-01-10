#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)
options(stringsAsFactors = F)

option_list <- list(
  make_option(c("--file"), type="character", default=NULL,
              help="Input file, which must include column name header"),
  make_option(c("--groupby"), type="character", default=NULL,
              help="Column to use for grouping"),
  make_option(c("--aggregate"), type="character", default=NULL,
              help="Column in the file to aggregate"),
  make_option(c("--sep"), type="character", default=",",
              help="Character to use to separate aggregated values")
)

opt <- parse_args(OptionParser(option_list=option_list))
df = readr::read_tsv(opt$file)
# This is to maintain the same order of values in the grouped summary as
# the order of unique values in the input
df[,opt$groupby] = factor(sapply(df[,opt$groupby], as.character), levels=unique(sapply(df[,opt$groupby], as.character)))

summaryCall = sprintf("paste0(%s, collapse=%s)", opt$aggregate, "opt$sep")
df2 = df %>%
  dplyr::group_by_(.dots = list(opt$groupby)) %>%
  dplyr::summarise_(.dots = setNames(summaryCall, opt$aggregate))

write.table(df2, file="", col.names=T, row.names=F, sep="\t", quote=F, na="")
