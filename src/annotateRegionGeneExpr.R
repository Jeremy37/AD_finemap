#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
genesFile = args[1]
exprFile = args[2]

trimGeneID = function(geneID) {
  gsub("\\..*", "", geneID)
}

genes.df = readr::read_tsv(genesFile) %>%
  dplyr::mutate(geneID = trimGeneID(geneID))

expr_in = readr::read_tsv(exprFile)
expr.mat = expr_in %>% column_to_rownames("ensgene") %>% select(-symbol)
# View(expr_in %>% filter(symbol %in% c("TSPAN14", "SPRED2", "NCK2", "CCDC6")) %>% select(-c(1,2)) %>% t())
# View(expr_in %>% filter(symbol %in% c("FCER1G")) %>% select(-c(1,2)) %>% t())
# Determine the percentile of microglia expression relative to GTEx tissues
# And do the same for the average brain expression relative to other tissues
getExprQuantile = function(x) {
  ranks = rank(x, ties.method = "average")
  ranks[is.na(x) | sum(is.na(x)) > 10] = NA
  ranks / sum(!is.na(x))
}
geneExprQuantiles = apply(expr.mat, MARGIN = 1, FUN=getExprQuantile) %>% t()
microglia_quantile = geneExprQuantiles[, grepl("primary_microglia", colnames(geneExprQuantiles))]
microglia_quantile[is.na(microglia_quantile)] = 0

brainExprAvg = expr.mat %>% select(contains("brain")) %>% rowMeans(na.rm=T)

expr.nobrain = expr.mat %>% select(-contains("brain"))
brainExprAvg[sapply(brainExprAvg, is.nan)] = NA
expr.brainAvg = cbind(brainExprAvg, expr.nobrain)
brain_quantile = apply(expr.brainAvg, MARGIN = 1, FUN=function(x) {ranks=rank(x)/(ncol(expr.nobrain)+1); ranks[is.na(x)] = NA; ranks}) %>% t() %>% .[,1]
brain_quantile[is.na(brain_quantile)] = 0

expr.to_annotate = cbind(microglia_quantile, brain_quantile, expr.mat)

# Round values to 3 decimal places
roundToString = function(x) {
  if (is.na(x)) { return("") }
  sprintf("%.3g", round(x, digits = 3))
}
expr.short = cbind(expr_in %>% select(ensgene), apply(expr.to_annotate, MARGIN=c(1, 2), FUN=roundToString))

genes.df = genes.df %>% dplyr::left_join(expr.short, by=c("geneID" = "ensgene"))

write.table(genes.df, file="", sep="\t", quote=F, na="", row.names=F, col.names=T)
