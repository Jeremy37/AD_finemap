#!/usr/bin/env Rscript
library(tidyverse)

setwd("/Users/jeremys/work/opentargets/AD_finemap/replication")

#ad_loci.df = read_tsv("AD_indep_snp_loci.txt")
# ad_meta.df = read_tsv("AD.IGAP1_GWAX_exclude_firsts_v5.meta.AD_indep_snps.txt")
# grace.df = read_tsv("GRACE_StageI.AD_indep_snps.txt")
# finngen_wide.df = read_tsv("finngen_r3_G6_AD_WIDE.AD_indep_snps.txt")
# finngen_lo.df = read_tsv("finngen_r3_AD_LO.AD_indep_snps.txt")

ad_loci.df = read_tsv("AD_lead_snp_loci.txt")
ad_meta.df = read_tsv("AD.IGAP1_GWAX_exclude_firsts_v5.meta.AD_lead_snps.txt")
grace.df = read_tsv("GRACE_StageI.AD_lead_snps.txt")
finngen_wide.df = read_tsv("finngen_r3_G6_AD_WIDE.AD_lead_snps.txt")
finngen_lo.df = read_tsv("finngen_r3_AD_LO.AD_lead_snps.txt")

global_meta = read_tsv("global_meta/global_meta.assoc_loci.tsv.gz") %>%
  select(CHR, pos_hg19, pos_hg38, SNP, ref, alt, locus,
         global_p, global_beta, global_se, META_P, META_BETA, META_SE, FREQ,
         grace_beta, grace_SE, grace_P, finngen_adwide_beta = finngen_beta, finngen_adwide_SE = finngen_SE, finngen_adwide_P = finngen_P)

merged.df = ad_loci.df %>%
  left_join(global_meta) %>%
  left_join(finngen_lo.df %>% select(SNP=rsids, finngen_adLO_beta=beta, finngen_adLO_se=sebeta, finngen_adLO_P=pval))


merged.df = merged.df %>%
  mutate(Grace_concordant = if_else(sign(grace_beta) == sign(META_BETA), "T", ""),
         Finngen_adwide_concordant = if_else(sign(finngen_adwide_beta) == sign(META_BETA), "T", ""),
         Finngen_adLO_concordant = if_else(sign(finngen_adLO_beta) == sign(META_BETA), "T", "") ) %>%
  select(CHR, pos_hg19, SNP, ref, alt, locus, global_p, global_beta, global_se, META_BETA, META_SE, META_P, FREQ, starts_with("grace"), starts_with("finngen_adwide"), starts_with("finngen_adLO"))

write_tsv(merged.df, path="AD_replication_merged.lead_snps.tsv", na="")
