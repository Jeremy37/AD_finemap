#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(annotables)

grch38_mainchrs = grch38 %>% filter(grepl("^([\\d+]|X|Y|MT)", chr, perl=T))

coloc_main.df = readr::read_tsv("coloc/output/coloc.AD.meta.eqtl_sqtl_colocs.txt", col_names = T)
coloc_cond1.df = readr::read_tsv("coloc/output/coloc.AD.meta.cond_1.eqtl_sqtl_colocs.txt", col_names = T) %>%
  mutate(locus_name = paste0(locus_name, ".cond_1"))
coloc_cond2.df = readr::read_tsv("coloc/output/coloc.AD.meta.cond_2.eqtl_sqtl_colocs.txt", col_names = T) %>%
  mutate(locus_name = paste0(locus_name, ".cond_2"))

coloc.df = bind_rows(coloc_main.df, coloc_cond1.df, coloc_cond2.df) %>%
  mutate(gwas_neglog10_pval = -log10(gwas_pval)) %>%
  select(locus_name, feature, ensembl_id, geneSymbol, H3=PP.H3, H4, gwas_neglog10_pval, qtl_pval, gwas_pval, qtl_lead, gwas_lead, chr, gwas_pos, qtl_pos, dataset, dataset_short)

# Separate into lists, and keep only the top coloc per gene
coloc.all.df = coloc.df %>%
  filter(!grepl("mac_sqtl", dataset_short),
         !grepl("mac_tx", dataset_short),
         H4 > 0.5) %>%
  left_join(grch38_mainchrs %>% select(ensgene, symbol) %>% unique(), by=c("geneSymbol"="symbol")) %>%
  mutate(ensembl_id = ifelse(!is.na(ensembl_id), ensembl_id, ensgene)) %>%
  filter(!is.na(ensembl_id)) %>%
  select(-ensgene)

getTopColocForEachGene = function(df) {
  df %>% group_by(geneSymbol) %>%
    arrange(desc(H4)) %>%
    summarise(locus_name = first(locus_name),
              feature = first(feature),
              ensembl_id = first(ensembl_id),
              H3 = first(H3),
              H4 = first(H4),
              gwas_neglog10_pval = first(gwas_neglog10_pval),
              qtl_pval = first(qtl_pval),
              gwas_pval = first(gwas_pval),
              qtl_lead = first(qtl_lead),
              gwas_lead = first(gwas_lead),
              chr = first(chr),
              gwas_pos = first(gwas_pos),
              qtl_pos = first(qtl_pos),
              dataset = first(dataset),
              dataset_short = first(dataset_short) )
}


coloc.all.flt.df = coloc.all.df %>% group_by(geneSymbol) %>%
  getTopColocForEachGene()

coloc.relevant.flt.df = coloc.all.df %>%
  filter(grepl("microglia|mono_eQTL|gtex.whole_blood|brain|gtex.spleen|xQTL_eQTL|gtex.brain|gtex.nerve|mac_eqtl", dataset_short, perl=T, ignore.case = T)) %>%
  getTopColocForEachGene()

coloc.irrrelevant.flt.df = coloc.all.df %>%
  filter(!grepl("microglia|mono_eQTL|gtex.whole_blood|brain|gtex.spleen|xQTL_eQTL|gtex.brain|gtex.nerve|mac_eqtl", dataset_short, perl=T, ignore.case = T)) %>%
  getTopColocForEachGene()

write.table(coloc.all.flt.df, "coloc/output/coloc.AD.meta.merged_cond.H4_gt_0.5.all.tsv",
            sep="\t", row.names=F, col.names=T, quote=F)

write.table(coloc.relevant.flt.df, "coloc/output/coloc.AD.meta.merged_cond.H4_gt_0.5.relevant.tsv",
            sep="\t", row.names=F, col.names=T, quote=F)

write.table(coloc.irrrelevant.flt.df, "coloc/output/coloc.AD.meta.merged_cond.H4_gt_0.5.irrelevant.tsv",
            sep="\t", row.names=F, col.names=T, quote=F)
