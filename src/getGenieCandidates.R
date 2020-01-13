#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
leadSnpsFile = args[1]
annotatedSnpsFile = args[2]
rpkmsFile = args[3]
outputRoot = args[4]

# setwd("/Users/jeremys/work/opentargets/AD_finemap/")
# annotatedSnpsFile = "/Users/jeremys/work/opentargets/AD_finemap/annotated/AD.IGAP1_GWAX_v3.annotated.selected.probable.tsv"
# rpkmsFile = "/Users/jeremys/work/opentargets/reference/tissueRPKM/tissues.selected.rpkm_average.txt.gz"
# outputRoot = "/Users/jeremys/work/opentargets/AD_finemap/annotated/AD.credible_sets.annotated.genie"

isTranscribedConsequence = function(consequence) {
  grepl("3_prime_UTR_variant|intron_variant|synonymous_variant|non_coding_transcript_exon_variant|5_prime_UTR_variant|missense_variant|splice_region_variant|stop_retained_variant|NMD_transcript_variant|stop_gained", consequence)
}

# Not used currently
#leadSnps.df = readr::read_tsv(leadSnpsFile)

snps.df = readr::read_tsv(annotatedSnpsFile, col_names = T,
                          col_types = cols(SYMBOL="c", EXON="c", INTRON="c", cDNA_position="c", CDS_position="c", Protein_position="c",
                                           PUBMED="c", MCAP_score="n", MetaLR_score="n", CADD_PHRED="n",
                                           fathmm_MKL_coding_score="n", integrated_fitCons_score="n", finemap_prob3="n")) %>%
  dplyr::filter(finemap_prob_nc > 0.001) %>%
  dplyr::mutate(is_transcribed_variant = isTranscribedConsequence(Consequence),
                transcr_prob = finemap_prob_nc * is_transcribed_variant,
                cum_transcr_prob = ave(transcr_prob, lead_chrpos, FUN=cumsum)) %>%
  dplyr::select(lead_chrpos, locus_name, snp, rsids, chr, pos_hg37, pos_hg38, Eff_allele,	A2, ref, alt,
                META_P, META_BETA, finemap_prob=finemap_prob_nc, cum_transcr_prob, everything())

# A table with sums of probability per locus for transcribed variants
locus.df = snps.df %>%
  dplyr::group_by(lead_chrpos, locus_name) %>%
  dplyr::summarise(total_transcr_prob = max(cum_transcr_prob),
                   total_transcr_prob_gt_0.01 = max(cum_transcr_prob * (finemap_prob > 0.01)),
                   max_transcr_snp_prob = ifelse(any(is_transcribed_variant), max(finemap_prob[is_transcribed_variant]), NA),
                   nSnps_to_0.5_prob = sum(is_transcribed_variant & cum_transcr_prob <= 0.5) + 1,
                   nSnps_to_0.8_prob = sum(is_transcribed_variant & cum_transcr_prob <= 0.8) + 1,
                   nSnps_to_0.9_prob = sum(is_transcribed_variant & cum_transcr_prob <= 0.9) + 1,
                   nSnps_gt_0.01_prob = sum(is_transcribed_variant & finemap_prob >= 0.01)) %>%
  dplyr::mutate(nSnps_to_0.5_prob = ifelse(total_transcr_prob < 0.5, NA, nSnps_to_0.5_prob),
                nSnps_to_0.8_prob = ifelse(total_transcr_prob < 0.8, NA, nSnps_to_0.8_prob),
                nSnps_to_0.9_prob = ifelse(total_transcr_prob < 0.9, NA, nSnps_to_0.9_prob)) %>%
  ungroup() %>%
  dplyr::arrange(-total_transcr_prob_gt_0.01, -nSnps_to_0.9_prob, -nSnps_to_0.8_prob, nSnps_to_0.5_prob)

gene.df = snps.df %>%
  dplyr::group_by(lead_chrpos, SYMBOL) %>%
  dplyr::summarise(gene_transcr_prob_gt_0.01 = sum(transcr_prob[transcr_prob > 0.01])) %>%
  ungroup()
  
# Add summary from the locus table so that we can order SNPs by loci of most interest
# (those with high total causal probability covered by transcribed variants)
snps.df = snps.df %>%
  dplyr::left_join(locus.df %>% select(-locus_name), by="lead_chrpos") %>%
  dplyr::left_join(gene.df, by=c("lead_chrpos", "SYMBOL")) %>%
  dplyr::select(lead_chrpos:cum_transcr_prob, total_transcr_prob:nSnps_gt_0.01_prob, gene_transcr_prob_gt_0.01, everything()) %>%
  dplyr::arrange(-total_transcr_prob_gt_0.01, nSnps_to_0.9_prob, nSnps_to_0.8_prob, nSnps_to_0.5_prob, desc(gene_transcr_prob_gt_0.01))


write.table(locus.df, file = paste0(outputRoot, ".locus.tsv"),
            sep="\t", quote=F, row.names=F, col.names=T, na="")

write.table(snps.df, file = paste0(outputRoot, ".snps.tsv"),
            sep="\t", quote=F, row.names=F, col.names=T, na="")

write.table(snps.df %>% filter(is_transcribed_variant), file = paste0(outputRoot, ".snps.transcribed.tsv"),
            sep="\t", quote=F, row.names=F, col.names=T, na="")

# Make a more detailed table that splits out sums by gene the SNPs are found in,
# rather than just by locus.
snps.df = readr::read_tsv(annotatedSnpsFile, col_names = T, col_types = cols(finemap_prob3="n")) %>%
  dplyr::filter(finemap_prob_nc > 0.001) %>%
  dplyr::mutate(is_transcribed_variant = isTranscribedConsequence(Consequence),
                transcr_prob = finemap_prob_nc * is_transcribed_variant,
                cum_transcr_prob = ave(transcr_prob, lead_chrpos, FUN=cumsum)) %>%
  dplyr::select(lead_chrpos, locus_name, snp, rsids, chr, pos_hg37, pos_hg38, Eff_allele,	A2, ref, alt,
                META_P, META_BETA, finemap_prob=finemap_prob_nc, cum_transcr_prob, everything())

locus.gene.df = snps.df %>%
  dplyr::group_by(lead_chrpos, locus_name, Gene, SYMBOL) %>%
  dplyr::summarise(gene_transcr_prob = sum(transcr_prob),
                   gene_transcr_prob_gt_0.01 = sum(transcr_prob[transcr_prob > 0.01]),
                   max_transcr_snp_prob = ifelse(any(is_transcribed_variant), max(finemap_prob[is_transcribed_variant]), NA),
                   nSnps_to_0.5_prob = sum(is_transcribed_variant & cum_transcr_prob <= 0.5) + 1,
                   nSnps_to_0.8_prob = sum(is_transcribed_variant & cum_transcr_prob <= 0.8) + 1,
                   nSnps_to_0.9_prob = sum(is_transcribed_variant & cum_transcr_prob <= 0.9) + 1,
                   nSnps_gt_0.01_prob = sum(is_transcribed_variant & finemap_prob >= 0.01)) %>%
  dplyr::mutate(nSnps_to_0.5_prob = ifelse(gene_transcr_prob < 0.5, NA, nSnps_to_0.5_prob),
                nSnps_to_0.8_prob = ifelse(gene_transcr_prob < 0.8, NA, nSnps_to_0.8_prob),
                nSnps_to_0.9_prob = ifelse(gene_transcr_prob < 0.9, NA, nSnps_to_0.9_prob)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(Gene))

locus.gene.df = locus.gene.df %>%
  dplyr::left_join(locus.df %>% dplyr::select(lead_chrpos, locus_transcr_prob = total_transcr_prob, locus_transcr_prob_gt_0.01 = total_transcr_prob_gt_0.01), by="lead_chrpos") %>%
  dplyr::select(lead_chrpos, locus_name, locus_transcr_prob, locus_transcr_prob_gt_0.01, everything()) %>%
  dplyr::arrange(-locus_transcr_prob_gt_0.01, -gene_transcr_prob_gt_0.01)

trimGeneID = function(geneID) {
  gsub("\\..*", "", geneID)
}
rpkm.df = readr::read_tsv(rpkmsFile)

locus.gene.df = locus.gene.df %>% dplyr::left_join(rpkm.df, by=c("Gene" = "gene_id"))

write.table(locus.gene.df, file = paste0(outputRoot, ".locus.gene.tsv"),
            sep="\t", quote=F, row.names=F, col.names=T, na="")
