#!/usr/bin/env Rscript
library(tidyverse)
options(stringsAsFactors = F)
#theme_set(theme_bw())

snp_meta = function(betas, ses) {
  n = length(betas)
  w = 1.0 / ses^2
  se = sqrt(1.0 / sum(w))
  b = sum(betas * w) / sum(w)
  z = b / se
  p = pnorm(abs(z), lower.tail = F) * 2
  q = sum( w * (betas - b)^2 )
  if (q != 0) {
    i2 = (q-1) / q
  } else {
    q = 0.0
  }
  if (i2 < 0) {
    i2 = 0
  }
  p_het = pchisq(q, df = n-1, lower.tail = F)
  return( list(b = b,
               se = se,
               p = p,
               i2 = i2,
               p_het = p_het) )
}

#snp_meta(c(1.4,1.3,1.7), c(0.4,0.3,0.2))


dir = "/Users/jeremys/work/opentargets/AD_finemap/"
library(Rsamtools)
library(liftOver)
library(rtracklayer)

get_header = function(fpath) {
  gzf = gzfile(fpath)
  header = readLines(gzf, 1)
  close(gzf)
  header
}

# Get the header for each GWAS to enable reading in Tabix results
ad_meta_header = get_header(file.path(dir, "summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz"))
ad_meta_cols = str_split(ad_meta_header, "\\t")[[1]]

grace_header = get_header(file.path(dir, "replication/GRACE_StageI.txt.bgz"))
grace_cols = str_split(grace_header, "\\t")[[1]]

finngen_header = get_header(file.path(dir, "replication/finngen_r3_G6_AD_WIDE.bgz"))
finngen_header = gsub("#", "", finngen_header)
finngen_cols = str_split(finngen_header, "\\t")[[1]]


do_three_study_meta = function(range) {
  ad_meta_sumstats = scanTabix(file = file.path(dir, "summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz"), param = range)
  ad_meta = read_tsv(ad_meta_sumstats[[1]], col_names = ad_meta_cols) %>%
    dplyr::rename(Effect_Allele = A1, NonEffect_Allele = A2)

  grace_sumstats = scanTabix(file = file.path(dir, "replication/GRACE_StageI.txt.bgz"), param = range)
  grace = read_tsv(grace_sumstats[[1]], col_names = grace_cols) %>%
    dplyr::rename(grace_OR = OR, grace_SE = SE, grace_P = P, grace_Effect_Allele = Effect_Allele, grace_NonEffect_Allele = NonEffect_Allele, grace_freq = FREQ_Effect_Allele) %>%
    mutate(grace_beta = log(grace_OR))

  finngen_sumstats = scanTabix(file = file.path(dir, "replication/finngen_r3_G6_AD_WIDE.bgz"), param = range)
  finngen = read_tsv(finngen_sumstats[[1]], col_names = finngen_cols) %>%
    dplyr::rename(CHR = chrom, pos_hg38 = pos, finngen_beta = beta, finngen_P = pval, finngen_SE = sebeta,
                  finngen_maf = maf, finngen_maf_cases = maf_cases, finngen_maf_controls = maf_controls)

  # LiftOver FinnGen back to GRCh37 coordinates
  path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  chain = import.chain(path)

  finngen_hg38 = GRanges(seqnames = finngen$CHR,
                         ranges = IRanges(start = finngen$pos_hg38, end = finngen$pos_hg38, width = 1,
                                          index = 1:nrow(finngen)))
  seqlevelsStyle(finngen_hg38) = "UCSC"
  finngen_hg19_gr = liftOver(finngen_hg38, chain)
  #class(finngen_hg19)
  #finngen_hg19_x = unlist(finngen_hg19)
  finngen_hg19_pos = as_tibble(as.data.frame(finngen_hg19_gr))
  # Remove rows that are duplicated, i.e. when an hg38 variant maps to two positions
  # in hg19
  finngen_hg19_pos.nodup = finngen_hg19_pos %>%
    filter(!duplicated(index)) %>%
    dplyr::select(seqnames, pos_hg19 = start, index)

  # Join the hg19 coords to the original FinnGen data
  finngen$index = 1:nrow(finngen)
  finngen_hg19 = finngen %>%
    inner_join(finngen_hg19_pos.nodup, by="index") %>%
    mutate(pos_id = paste(CHR, pos_hg19, ref, alt, sep = "_")) %>%
    dplyr::select(-seqnames, -index) %>%
    filter(!duplicated(pos_hg19))


  ###############################################################################
  # We have all the GWAS in hg19 coords. Now we need to join them by chr:pos:ref:alt.
  # For AD_meta and Gr@ace we have only effect and non-effect alleles, not ref/alt.
  # So we will first match with FinnGen to see if chr:pos:non-eff:eff matches
  # with chr:pos:ref:alt. We'll keep those matches, and then add back any that
  # match chr:pos:eff:non-eff, after we reverse the beta.

  # Join AD meta with FinnGen
  # To join two GWAS together, we need to ensure that effect directions match.
  # The alt allele is effect allele in FinnGen. We first join cases where ref/alt
  # match non-effect/effect allele in AD meta (or Grace). These should have the same
  # direction for beta. Then we join the others, and reverse the effect direction
  # so that ref/alt match non-eff/eff.
  finngen_join = finngen_hg19 %>% dplyr::select(-CHR, -rsids)
  ad_meta_harmonised_1 = ad_meta %>%
    mutate(pos_id = paste(CHR, BP, NonEffect_Allele, Effect_Allele, sep = "_")) %>%
    inner_join(finngen_join, by="pos_id")
  ad_meta_harmonised_2 = ad_meta %>%
    mutate(pos_id = paste(CHR, BP, Effect_Allele, NonEffect_Allele, sep = "_")) %>%
    inner_join(finngen_join, by="pos_id")

  ad_meta_harmonised_tmp = bind_rows(
    ad_meta_harmonised_1,
    ad_meta_harmonised_2 %>% mutate(Eff_Allele_tmp = Effect_Allele,
                                    Effect_Allele = NonEffect_Allele,
                                    NonEffect_Allele = Eff_Allele_tmp,
                                    META_BETA = -META_BETA,
                                    GWAS_BETA = -GWAS_BETA,
                                    GWAX_UKBB_BETA = -GWAX_UKBB_BETA) %>%
      dplyr::select(-Eff_Allele_tmp)
  )
  ad_meta_leftover = ad_meta %>%
    mutate(pos_id = paste(CHR, BP, NonEffect_Allele, Effect_Allele, sep = "_"),
           pos_id_rev = paste(CHR, BP, Effect_Allele, NonEffect_Allele, sep = "_")) %>%
    filter(!(pos_id %in% ad_meta_harmonised_1$pos_id) & !(pos_id_rev %in% ad_meta_harmonised_2$pos_id)) %>%
    dplyr::select(-pos_id_rev) %>%
    left_join(finngen_join, by="pos_id")
  ad_meta_harmonised = bind_rows(ad_meta_harmonised_tmp, ad_meta_leftover)
  if (nrow(ad_meta_harmonised) != nrow(ad_meta)) {
    warning(sprintf("Different number of rows in ad_meta_harmonised (%d) and ad_meta (%d)", nrow(ad_meta_harmonised), nrow(ad_meta)))
  }

  # We do the same to join with Grace.
  all_harmonised_1 = grace %>%
    mutate(pos_id = paste(CHR, BP, grace_NonEffect_Allele, grace_Effect_Allele, sep = "_")) %>%
    filter(!duplicated(pos_id)) %>%
    dplyr::select(-CHR, -BP, -rsID) %>%
    inner_join(ad_meta_harmonised, by="pos_id")
  all_harmonised_2 = grace %>%
    mutate(pos_id = paste(CHR, BP, grace_Effect_Allele, grace_NonEffect_Allele, sep = "_")) %>%
    filter(!duplicated(pos_id)) %>%
    dplyr::select(-CHR, -BP, -rsID) %>%
    inner_join(ad_meta_harmonised, by="pos_id")

  all_harmonised_tmp = bind_rows(
    all_harmonised_1,
    all_harmonised_2 %>% mutate(Eff_Allele_tmp = grace_Effect_Allele,
                                grace_Effect_Allele = grace_NonEffect_Allele,
                                grace_NonEffect_Allele = Eff_Allele_tmp,
                                grace_beta = -grace_beta) %>%
      dplyr::select(-Eff_Allele_tmp)
  )
  all_harmonised_leftover = ad_meta_harmonised %>%
    mutate(pos_id_rev = paste(CHR, BP, Effect_Allele, NonEffect_Allele, sep = "_")) %>%
    filter(!(pos_id %in% all_harmonised_tmp$pos_id) & !(pos_id_rev %in% all_harmonised_tmp$pos_id)) %>%
    dplyr::select(-pos_id_rev) %>%
    left_join(grace %>%
                mutate(pos_id = paste(CHR, BP, grace_NonEffect_Allele, grace_Effect_Allele, sep = "_")) %>%
                dplyr::select(-CHR, -BP, -rsID), by="pos_id")
  all_harmonised = bind_rows(all_harmonised_tmp, all_harmonised_leftover)
  print(sprintf("Removing %d SNPs with invalid beta in the AD meta-analysis", sum(is.na(all_harmonised$META_BETA))))
  all_harmonised = all_harmonised %>% filter(!is.na(META_BETA))

  missing_finngen = sum(is.na(all_harmonised$finngen_beta))
  missing_grace = sum(is.na(all_harmonised$grace_beta))
  missing_both = sum(is.na(all_harmonised$finngen_beta) & is.na(all_harmonised$grace_beta))
  print(sprintf("%d SNPs total, %d missing in FinnGen, %d missing in Grace, %d missing in both",
                nrow(all_harmonised), missing_finngen, missing_grace, missing_both))

  # Check that the betas are generally all in the same direction.

  # ggplot(all_harmonised %>% filter(META_P < 1e-8),
  #        aes(x=META_BETA, y=finngen_beta)) +
  #   geom_point(alpha=0.5)
  #
  # ggplot(all_harmonised %>% filter(META_P < 1e-8),
  #        aes(x=grace_beta, y=finngen_beta)) +
  #   geom_point(alpha=0.5)
  # They are, especially for the best-powered SNPs... with a few exceptions

  do_ukb_grace_finngen_meta = function(i, df) {
    if (is.na(df$finngen_beta[i]) && is.na(df$grace_beta[i])) {
      return( list(b = df$META_BETA[i],
                   se = df$META_SE[i],
                   p = df$META_P[i],
                   i2 = df$I2[i],
                   p_het = df$HET_P[i]) )
    }
    if (is.na(df$finngen_beta[i])) {
      beta = c(df$META_BETA[i], df$grace_beta[i])
      se = c(df$META_SE[i], df$grace_SE[i])
    } else if (is.na(df$grace_beta[i])) {
      beta = c(df$META_BETA[i], df$finngen_beta[i])
      se = c(df$META_SE[i], df$finngen_SE[i])
    } else {
      beta = c(df$META_BETA[i], df$finngen_beta[i], df$grace_beta[i])
      se = c(df$META_SE[i], df$finngen_SE[i], df$grace_SE[i])
    }
    return(snp_meta(betas = beta, ses = se))
  }

  do_ukb_finngen_meta = function(i, df) {
    if (is.na(df$finngen_beta[i])) {
      return( list(b = df$META_BETA[i],
                   se = df$META_SE[i],
                   p = df$META_P[i],
                   i2 = df$I2[i],
                   p_het = df$HET_P[i]) )
    }
    beta = c(df$META_BETA[i], df$finngen_beta[i])
    se = c(df$META_SE[i], df$finngen_SE[i])
    return(snp_meta(betas = beta, ses = se))
  }


  #test_harmonised = all_harmonised[1:10,]
  #meta_res = lapply(1:nrow(test_harmonised), FUN = do_ukb_grace_finngen_meta, test_harmonised)
  global_meta_res = lapply(1:nrow(all_harmonised), FUN = do_ukb_grace_finngen_meta, all_harmonised)
  # Make the output table, renaming with "global_" for the meta-analysis of all
  # studies: Kunkle, UKB, Grace, FinnGen meta-analysis.
  global_meta_df = bind_rows(global_meta_res) %>%
    dplyr::rename(global_beta = b,
                  global_se = se,
                  global_p = p,
                  global_i2 = i2,
                  global_p_het = p_het)

  kfu_meta_res = lapply(1:nrow(all_harmonised), FUN = do_ukb_finngen_meta, all_harmonised)
  # Make the output table, renaming with "global_" for the meta-analysis of all
  # studies: Kunkle, UKB, Grace, FinnGen meta-analysis.
  kfu_meta_df = bind_rows(kfu_meta_res) %>%
    dplyr::rename(kfu_beta = b,
                  kfu_se = se,
                  kfu_p = p,
                  kfu_i2 = i2,
                  kfu_p_het = p_het)


  global_df = bind_cols(all_harmonised, global_meta_df, kfu_meta_df) %>%
    dplyr::select(pos_id, CHR, pos_hg19=BP, pos_hg38, SNP, ref, alt,
                  global_beta, global_se, global_p, global_i2, global_p_het,
                  kfu_beta, kfu_se, kfu_p, kfu_i2, kfu_p_het,
                  META_BETA, META_SE, META_P, DIRECT, I2, HET_P, FREQ, INFO,
                  Kunkle_BETA = GWAS_BETA, Kunkle_SE = GWAS_SE, Kunkle_P = GWAS_P,
                  GWAX_UKBB_BETA, GWAX_UKBB_SE, GWAX_UKBB_P,
                  finngen_beta, finngen_SE, finngen_P, finngen_maf, finngen_maf_cases, finngen_maf_controls, nearest_genes,
                  grace_beta, grace_SE, grace_P, grace_OR, grace_freq) %>%
    arrange(pos_hg19)
  return(global_df)
}


chroms = c(1:22)
for (chr in chroms) {
  print(sprintf("Doing chr%d", chr))
  range = GRanges(seqnames = chr, ranges = IRanges(start = 1, end = 3e8))
  chrom_global_df = do_three_study_meta(range)
  chr_fname = file.path(dir, sprintf("replication/global_meta/chr%d.global_meta.tsv.gz", chr))
  write_tsv(chrom_global_df, path = chr_fname)
}


#ggplot(chrom_global_df, aes(x=-log10(global_p_het))) + geom_histogram()

