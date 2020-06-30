#!/usr/bin/env Rscript
library(tidyverse)
library(annotables)

options(stringsAsFactors = F)

args <- commandArgs(trailingOnly = TRUE)
lit_genes_file = args[1]
nearest_genes_file = args[2]
coloc_file = args[3]
network_genes_file = args[4]

grch38_mainchrs = grch38 %>% filter(grepl("^([\\d+]|X|Y|MT)", chr, perl=T))

# setwd("/Users/jeremys/work/opentargets/AD_finemap/")
# lit_genes_file="genes/literature_genes.tsv"
# nearest_genes_file="genes/AD.loci.nearest_gene.pc.tsv"
# coloc_file="coloc/output/coloc.AD.meta.merged_cond.H4_gt_0.5.all.tsv"
# network_genes_file = "network/network_genes.txt"

lit.df = readr::read_tsv(lit_genes_file) %>%
  mutate(literature_gene = 1)

nearest.df = readr::read_tsv(nearest_genes_file) %>%
  mutate(gene_id = gsub("\\..*", "", gene_id)) %>%
  select(locus, gene_id, symbol, dist) %>%
  mutate(nearest_gene = 1)

coloc.df = readr::read_tsv(coloc_file) %>%
  select(gene_id = ensembl_id, symbol = geneSymbol, locus = locus_name, coloc_H4 = H4) %>%
  filter(coloc_H4 >= 0.8) %>%
  mutate(coloc_gene = 1,
         locus = gsub("MIR142", "TSPOAP1", locus)) # Hack due to historical changing of locus name

network_genes.df = readr::read_tsv(network_genes_file) %>%
  mutate(in_network = 1)

seed_genes.df = lit.df %>%
  full_join(nearest.df, by=c("locus", "gene_id", "symbol")) %>%
  full_join(coloc.df, by=c("locus", "gene_id", "symbol")) %>%
  left_join(network_genes.df, by="gene_id") %>%
  left_join(grch38_mainchrs %>% select(ensgene, biotype) %>% unique(), by=c("gene_id"="ensgene")) %>%
  group_by(gene_id) %>%
  mutate(in_network = if_else(is.na(in_network) | !in_network, 0, 1),
         gene_sum = sum(nearest_gene, coloc_gene, literature_gene, na.rm = T) * in_network) %>%
  group_by(locus) %>%
  mutate(locus_sum = sum(gene_sum)) %>%
  group_by(gene_id) %>%
  mutate(gene_weight = gene_sum / locus_sum) %>%
  select(locus, gene_id, symbol, biotype, nearest_gene, coloc_gene, literature_gene, in_network, gene_weight, locus_sum, dist, coloc_H4, publication, note)

cat(readr::format_tsv(seed_genes.df, na="", col_names = T))
