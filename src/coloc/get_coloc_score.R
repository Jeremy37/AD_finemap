#!/usr/bin/env Rscript
# Determines a summary coloc score per locus and gene, aggregating over QTL datasets
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
coloc_file = args[1]
dataset_groups_file = args[2]
category_weights_file = args[3]
output_root = args[4]
#root = "AD_finemap"
# coloc_file = file.path(root, "coloc/output/coloc.AD.meta.eqtl_sqtl_colocs.txt")
# dataset_groups_file = file.path(root, "coloc/dataset_groups.tsv")
# category_weights_file = file.path(root, "coloc/coloc_group_weights.tsv")
# output_root = file.path(root, "coloc/coloc.AD.meta")

colocs.df = readr::read_tsv(coloc_file) %>%
  filter(!is.na(locus_name))
# Should have columns locus_name, dataset_short, ensembl_id, geneSymbol, H4

# Make sure there is only one entry per locus : gene : dataset.
# In the case of sQTLs there can be more than one, so we take the
# one with the highest H4.
colocs.df = colocs.df %>% group_by(locus_name, dataset_short, ensembl_id) %>%
  arrange(desc(H4)) %>%
  filter(!duplicated(paste(locus_name, signal, dataset_short, ensembl_id, sep=".")))

# Table of weights for combining colocs
dataset_groups = readr::read_tsv(dataset_groups_file)
category_weights = readr::read_tsv(category_weights_file)

colocs.df = colocs.df %>%
  left_join(dataset_groups, by="dataset_short") %>%
  left_join(category_weights %>% select(category, category_min), by="category")

# For computing the score, filter out colocs which are below the min for that group.
# But we retain these in the main data.frame for display later.
colocs.toscore.df = colocs.df %>% filter(H4 >= category_min)

combineColocScores = function(H4vec) {
  H4vec = H4vec[order(H4vec, decreasing = T)]
  score = 0
  for (i in 1:length(H4vec)) {
    score = score + (1 - score) * H4vec[i] / i
  }
  score
}

combineScoreGroups = function(scoreVec) {
  scoreVec = scoreVec[order(scoreVec, decreasing = T)]
  score = 0
  for (i in 1:length(scoreVec)) {
    score = score + (1 - score) * scoreVec[i]
  }
  score
}

# The first scoring method I tried, which does 3 steps:
# * within each category, combine H4 values as a harmonic sum
# * multiply these scores by the weight for the category
# * combine scores across categories (without harmonic sum weighting)
getColocScoreV1 = function(gene_colocs.df, groupvars) {
  category_scores.df = gene_colocs.df %>% group_by(category) %>%
    summarise(score = combineColocScores(H4))
  category_scores.df = category_scores.df %>%
    left_join(category_weights %>% select(category, weight), by="category") %>%
    mutate(weighted_score = score * weight)
  combined_score = combineScoreGroups(category_scores.df$weighted_score)
  # Strangely, this group_map function works differently on my home mac vs. on the cluster
  # despite having the same version of tidyverse, etc. On my home mac, it doesn't allow me
  # to include the groupvars in the data.frame, whereas on the cluster I must include it.
  #data.frame(colocScore = combined_score)
  data.frame(groupvars, colocScore = combined_score)
}

# The second scoring method I tested attempts to reduce further the impact of
# datasets with many QTL subsets, e.g. the 15 mostly T-cell subsets of Schmiedel et al.
# It does this by taking the simple max H4 for each gene within each dataset.
# These values are then combined as before - harmonic sum by category,
# and standard combining method across categories.
getColocScoreV2 = function(gene_colocs.df, groupvars) {
  dataset_group_scores.df = gene_colocs.df %>% group_by(dataset_group) %>%
    summarise(category = first(category),
              maxH4 = max(H4))
  category_scores.df = dataset_group_scores.df %>% group_by(category) %>%
    summarise(score = combineColocScores(maxH4))
  category_scores.df = category_scores.df %>%
    left_join(category_weights %>% select(category, weight), by="category") %>%
    mutate(weighted_score = score * weight)
  combined_score = combineScoreGroups(category_scores.df$weighted_score)
  data.frame(groupvars, colocScore = combined_score)
}

# The simplest scoring method, this approach does two things:
# * first multiply each H4 score by a weight based on its category
# * combine all scores together (arranged from highest to lowest) as a harmonic sum
getColocScoreV3 = function(gene_colocs.df, groupvars) {
  gene_colocs.df = gene_colocs.df %>%
    left_join(category_weights %>% select(category, weight), by="category") %>%
    mutate(weightedH4 = H4 * weight)
  combined_score = combineColocScores(gene_colocs.df$weightedH4)
  data.frame(groupvars, colocScore = combined_score)
}


scores.df = colocs.toscore.df %>% group_by(locus_name, signal, ensembl_id) %>%
  group_map(.f = ~ getColocScoreV3(.x, .y), keep=T)
scores.df = bind_rows(scores.df)

scores.df = scores.df %>%
  arrange(locus_name, signal, desc(colocScore))

colocs.spread.df = colocs.df %>%
  select(locus_name, signal, dataset_short, H4, ensembl_id, geneSymbol) %>%
  spread(key = dataset_short, value = H4) %>%
  left_join(scores.df, by=c("locus_name", "signal", "ensembl_id")) %>%
  select(locus_name, signal, ensembl_id, geneSymbol, colocScore, everything()) %>%
  arrange(locus_name, signal, desc(colocScore))

write.table(colocs.spread.df, file=paste0(output_root, ".colocScores.per_signal.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")
write.table(scores.df, file=paste0(output_root, ".colocScores.per_signal.tidy.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")

scores.perGene.df = scores.df %>%
  group_by(locus_name, ensembl_id) %>%
  summarise(colocScore = max(colocScore)) %>%
  arrange(desc(colocScore))

colocs.perGene.spread.df = colocs.df %>%
  select(locus_name, signal, dataset_short, H4, ensembl_id, geneSymbol) %>%
  group_by(locus_name, ensembl_id, geneSymbol, dataset_short) %>%
  summarise(H4 = max(H4)) %>%
  ungroup() %>%
  spread(key = dataset_short, value = H4) %>%
  left_join(scores.perGene.df, by=c("locus_name", "ensembl_id"))

colocs.perGene.spread.df = colocs.perGene.spread.df %>%
  arrange(desc(colocScore)) %>%
  mutate(scoreRank = 1:nrow(colocs.perGene.spread.df)) %>%
  arrange(locus_name, desc(colocScore)) %>%
  select(locus_name, ensembl_id, geneSymbol, colocScore, scoreRank, everything())

write.table(colocs.perGene.spread.df, file=paste0(output_root, ".colocScores.per_gene.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")
write.table(scores.perGene.df, file=paste0(output_root, ".colocScores.per_gene.tidy.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")

write.table(colocs.perGene.spread.df %>% filter(colocScore > 0.5),
            file=paste0(output_root, ".colocScores.per_gene.gt_0.5.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")


###############################################################################
# Get coloc scores per gene per cell type

getColocCelltypeScore = function(gene_colocs.df, groupvars) {
  group_scores.df = gene_colocs.df %>% group_by(celltype_group) %>%
    summarise(cellTypeScore = combineColocScores(H4))
  data.frame(groupvars, group_scores.df)
}

celltype_groups = unique(colocs.toscore.df$celltype_group)

scores.celltype.df = colocs.toscore.df %>% group_by(locus_name, signal, ensembl_id, geneSymbol) %>%
  group_map(.f = ~ getColocCelltypeScore(.x, .y), keep=T)
scores.celltype.df = bind_rows(scores.celltype.df)

scores.celltype.perGene.df = bind_rows(scores.celltype.df) %>%
  group_by(locus_name, ensembl_id, geneSymbol, celltype_group) %>%
  summarise(cellTypeScore = max(cellTypeScore)) %>% # Get the max across signals
  arrange(locus_name, desc(cellTypeScore))

write.table(scores.celltype.perGene.df, file=paste0(output_root, ".cellTypeScores.per_gene.tidy.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")

scores.celltype.perGene.spread.df = scores.celltype.perGene.df %>%
  spread(key = celltype_group, value = cellTypeScore) %>%
  left_join(scores.perGene.df, by=c("locus_name", "ensembl_id")) %>%
  select(locus_name, ensembl_id, geneSymbol, colocScore, everything()) %>%
  arrange(locus_name, desc(colocScore))
write.table(scores.celltype.perGene.spread.df, file=paste0(output_root, ".cellTypeScores.per_gene.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")

