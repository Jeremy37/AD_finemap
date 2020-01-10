#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
regionsFile = args[1]
genesFileNearest = args[2]
genesFile1 = args[3]
genesFile2 = args[4]

regions.df = readr::read_tsv(regionsFile, col_names = T, col_types = cols(Chr="c", pos="c", p="c")) %>%
  rename(snp=lead_SNP)

# Merge together the 200 kb and 1 Mb window for genes near AD regions
genesNearest.df = readr::read_tsv(genesFileNearest) %>% dplyr::select(snp, nearest_gene = symbol)
genes1.df = readr::read_tsv(genesFile1) %>% dplyr::rename(snp = leadSNP, gene_within_100kb = symbol)
genes2.df = readr::read_tsv(genesFile2) %>% dplyr::rename(snp = leadSNP, gene_within_500kb = symbol)

region_genes.df = regions.df %>%
  dplyr::left_join(genesNearest.df, by="snp") %>%
  dplyr::left_join(genes1.df, by="snp") %>%
  dplyr::left_join(genes2.df, by="snp")

write.table(region_genes.df, file="", sep="\t", quote=F, na="", row.names=F, col.names=T)
