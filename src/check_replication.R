#!/usr/bin/env Rscript
library(tidyverse)

ad_loci.df = read_tsv("AD_indep_snp_loci.txt")
ad_meta.df = read_tsv("AD.IGAP1_GWAX_exclude_firsts_v5.meta.AD_indep_snps.txt")
grace.df = read_tsv("GRACE_StageI.AD_indep_snps.txt")
finngen_wide.df = read_tsv("finngen_r3_G6_AD_WIDE.AD_indep_snps.txt")
finngen_lo.df = read_tsv("finngen_r3_AD_LO.AD_indep_snps.txt")

ad_meta.df = read_tsv("AD.IGAP1_GWAX_exclude_firsts_v5.meta.AD_lead_snps.txt")
grace.df = read_tsv("GRACE_StageI.AD_lead_snps.txt")
finngen_wide.df = read_tsv("finngen_r3_G6_AD_WIDE.AD_lead_snps.txt")
finngen_lo.df = read_tsv("finngen_r3_AD_LO.AD_lead_snps.txt")


merged.df = ad_meta.df %>%
  select(CHR, BP, SNP, EFF_ALLELE=A1, A2, META_BETA, META_SE, META_P, FREQ) %>%
  left_join(ad_loci.df) %>%
  left_join(grace.df %>% select(SNP=rsID, Grace_eff_allele=Effect_Allele, NonEffect_Allele, Grace_freq=FREQ_Effect_Allele, Grace_OR=OR, Grace_SE=SE, Grace_P=P)) %>%
  left_join(finngen_wide.df %>% select(SNP=rsids, ref, Finngen_eff_allele=alt, Finngen_maf=maf, Finngen_adwide_beta=beta, Finngen_adwide_se=sebeta, Finngen_adwide_P=pval)) %>%
  left_join(finngen_lo.df %>% select(SNP=rsids, Finngen_adLO_beta=beta, Finngen_adLO_se=sebeta, Finngen_adLO_P=pval)) %>%
  select(CHR, BP, SNP, EFF_ALLELE, A2, locus_name, everything())

merged.df = merged.df %>%
  mutate(risk_allele = if_else(META_BETA > 0, EFF_ALLELE, A2),
         Grace_risk_allele = if_else(Grace_OR > 1, Grace_eff_allele, NonEffect_Allele),
         Finngen_adwide_risk_allele = if_else(Finngen_adwide_beta > 0, Finngen_eff_allele, ref),
         Finngen_adLO_risk_allele = if_else(Finngen_adLO_beta > 0, Finngen_eff_allele, ref),
         Grace_concordant = if_else(Grace_risk_allele == risk_allele, "T", ""),
         Finngen_adwide_concordant = if_else(Finngen_adwide_risk_allele == risk_allele, "T", ""),
         Finngen_adLO_concordant = if_else(Finngen_adLO_risk_allele == risk_allele, "T", "") ) %>%
  select(CHR, BP, SNP, EFF_ALLELE, A2, locus_name, META_BETA, META_SE, META_P, FREQ, starts_with("Grace"), Finngen_eff_allele, Finngen_maf, starts_with("Finngen_adwide"), starts_with("Finngen_adLO"))

write_tsv(merged.df, path="AD_replication_merged.lead_snps.tsv", na="")
