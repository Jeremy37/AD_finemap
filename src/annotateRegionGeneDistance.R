#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(annotables)

option_list <- list(
  make_option(c("--loci"), type="character", default=NULL,
              help="Input file with list of independent associated loci (produced by earlier pipeline steps)"),
  make_option(c("--genes"), type="character", default=NULL,
              help="Genes overlapping a window around AD loci")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt = list(loci = "AD.loci.tsv", genes = "genes/AD.loci.1Mb_window.gene_overlaps.all.tsv")
trimGeneID = function(geneID) {
  gsub("\\..*", "", geneID)
}

regions.df = readr::read_tsv(opt$loci, col_names = T, col_types = cols(.default = col_character(), start="i", stop="i", lead_p="d"))
genes.df = readr::read_tsv(opt$genes, col_names=T) %>%
  mutate(geneID = trimGeneID(geneID))

grch38_noDupEnsg = grch38 %>% filter(!duplicated(ensgene))
genes.df = genes.df %>%
  left_join(grch38_noDupEnsg %>% select(ensgene, biotype), by=c("geneID" = "ensgene"))

getLeadSNPs = function(i) {
  positions = as.integer(strsplit(regions.df$pos[i], ",", fixed = T)[[1]])
  snps = strsplit(regions.df$SNPs[i], ",", fixed = T)[[1]]
  pvals = as.numeric(strsplit(regions.df$p[i], ",", fixed = T)[[1]])
  data.frame(locus_name = regions.df$locus_name[i], snp = snps, chr = regions.df$Chr[i], lead_pos = positions, lead_p = pvals)
}
leadSNPs.df = bind_rows(lapply(1:nrow(regions.df), getLeadSNPs))
  
getDist = function(pos, region_start, region_end) {
  if (pos < region_start) {
    pos - region_start
  } else if (pos > region_end) {
    pos - region_end
  } else {
    0
  }
}

getDists = function(pos, region_start, region_end) {
  sapply(1:length(pos), FUN=function(i) getDist(pos[i], region_start[i], region_end[i]))
}

# Add region info to genes
genes.df2 = genes.df %>%
  left_join(leadSNPs.df %>% select(locus = locus_name, lead_pos, snp, lead_p), by="locus") %>%
  group_by(leadSNP, geneID) %>%
  mutate(dist = getDists(lead_pos, gene_start, gene_end)) %>%
  group_by(leadSNP, geneID) %>%
  summarise(minInd = which.min(dist),
            locus = first(locus),
            symbol = first(symbol),
            chr = first(chr),
            gene_start = first(gene_start),
            gene_end = first(gene_end),
            lead_pos = lead_pos[minInd],
            snp = snp[minInd],
            lead_p = lead_p[minInd],
            dist = dist[which.min(abs(dist))],
            biotype = first(biotype)) %>%
  ungroup() %>%
  select(leadSNP, locus, geneID, symbol, biotype, chr, gene_start, gene_end, lead_pos, snp, lead_p, dist) %>%
  mutate(lead_neglog10p = -log10(lead_p)) %>%
  mutate(chrInt = as.integer(gsub("chr", "", chr))) %>%
  arrange(chrInt, locus, dist) %>%
  select(-chrInt)

write.table(genes.df2, file="", sep="\t", quote=F, na="", row.names=F, col.names=T)
