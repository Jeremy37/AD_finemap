#!/usr/bin/env Rscript
library(tidyverse)
library(pheatmap)

ad_dir = "/Users/jeremys/work/opentargets/AD_finemap"
exprFileSingleCell = file.path(ad_dir, "expression/brain_cell_type_tpm.tsv.gz")
exprFileBulk = file.path(ad_dir, "expression/tissues.combined.tpm.tsv.gz")

# Converts expression to a uniform quantile within each tissue or cell type
getExprNorm = function(expr.df) {
  # The first two columns of the expression table should be gene_id and symbol,
  # and the remaining columns are any kind of gene expression quantification.
  expr.df = expr.df %>% select(gene_id, symbol, everything())
  expr.mat = expr.df[, -(1:2)]
  expr.mat[is.na(expr.mat)] = 0
  N = nrow(expr.mat)
  expr.norm = apply(expr.mat, MARGIN=2, FUN=function(x) rank(x, ties.method="average") / N)
  cbind(expr.df[, 1:2], expr.norm)
}

expr.sc.df = read_tsv(exprFileSingleCell) %>%
  rename(symbol = hgnc_symbol) %>%
  select(-Exclude, -VLMC) # Remove "exclude" cell type, and VLMC which has too few cells
expr.sc.norm = getExprNorm(expr.sc.df)
colnames(expr.sc.norm) = gsub(" |/", "_", colnames(expr.sc.norm))
write_tsv(expr.sc.norm, path=file.path(ad_dir, "expression/expr.sc.brain.norm.tsv.gz"))

expr.bulk.df = read_tsv(exprFileBulk) %>% rename(gene_id = ensgene)
colnames(expr.bulk.df) = gsub(" - ", "_", colnames(expr.bulk.df))
colnames(expr.bulk.df) = gsub(" \\(", "_", colnames(expr.bulk.df))
colnames(expr.bulk.df) = gsub(" |/|\\+|\\-", "_", colnames(expr.bulk.df))
colnames(expr.bulk.df) = gsub("\\)", "", colnames(expr.bulk.df))
expr.bulk.norm = getExprNorm(expr.bulk.df)

# Select a subset of the GTEx brain regions so that the overall expression matrix
# isn't heavily weighted towards brain. We also select a few of our cell line
# models to know how well they are enriched.
gtex_cols = colnames(expr.bulk.df %>% select(Adipose_Subcutaneous:Whole_Blood))
gtex_selected_cols = colnames(expr.bulk.df %>% select(Adipose_Subcutaneous:Whole_Blood, -starts_with("Brain"), Brain_Cerebellum, Brain_Cortex, Brain_Hippocampus, Brain_Substantia_nigra))
expr.bulk.gtex = expr.bulk.df %>% select(gene_id:symbol, iNeuron, neuron, primary_microglia, ipsc_microglia, one_of(gtex_selected_cols))
write_tsv(expr.bulk.gtex, path=file.path(ad_dir, "expression/expr.bulk.gtex.selected.tpm.tsv.gz"))

N = ncol(expr.bulk.gtex) - 2
expr.bulk.gtex.norm = getExprNorm(expr.bulk.gtex %>% select(gene_id:symbol, iNeuron, neuron, primary_microglia, ipsc_microglia, one_of(gtex_selected_cols)))
write_tsv(expr.bulk.gtex.norm, path=file.path(ad_dir, "expression/expr.bulk.gtex.norm.tsv.gz"))

expr.bulk.relative = t(apply(expr.bulk.gtex.norm %>% select(-gene_id:-symbol), MARGIN=1, FUN=function(x) rank(x) / N * 100))
expr.bulk.relative = cbind(expr.bulk.gtex.norm %>% select(gene_id:symbol), expr.bulk.relative)
write_tsv(expr.bulk.relative, path=file.path(ad_dir, "expression/expr.bulk.gtex.relative.tsv.gz"))

# Select the eQTL catalogue datasets, as well as a few of our cell lines.
expr.bulk.eqtl_catalogue = expr.bulk.df %>% select(gene_id:symbol, iNeuron, neuron, primary_microglia, ipsc_microglia, everything(), -one_of(gtex_cols))
write_tsv(expr.bulk.eqtl_catalogue, path=file.path(ad_dir, "expression/expr.bulk.eqtl_catalogue.tpm.tsv.gz"))
expr.bulk.eqtl_catalogue.norm = expr.bulk.norm %>% select(gene_id:symbol, iNeuron, neuron, primary_microglia, ipsc_microglia, everything(), -one_of(gtex_cols))
write_tsv(expr.bulk.eqtl_catalogue.norm, path=file.path(ad_dir, "expression/expr.bulk.eqtl_catalogue.norm.tsv.gz"))

N = ncol(expr.bulk.eqtl_catalogue.norm) - 2
expr.bulk.relative = t(apply(expr.bulk.eqtl_catalogue.norm %>% select(-gene_id:-symbol), MARGIN=1, FUN=function(x) rank(x) / N * 100))
expr.bulk.relative = cbind(expr.bulk.norm %>% select(gene_id:symbol), expr.bulk.relative)
write_tsv(expr.bulk.relative, path=file.path(ad_dir, "expression/expr.bulk.eqtl_catalogue.relative.tsv.gz"))


N = ncol(expr.sc.norm) - 2
expr.sc.relative = t(apply(expr.sc.norm %>% select(-gene_id:-symbol), MARGIN=1, FUN=function(x) rank(x) / N * 100))
expr.sc.relative = cbind(expr.sc.norm %>% select(gene_id:symbol), expr.sc.relative)
write_tsv(expr.sc.relative, path=file.path(ad_dir, "expression/expr.sc.brain.relative.tsv.gz"))

