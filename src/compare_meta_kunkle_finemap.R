library(tidyverse)

setwd("/Users/jeremys/work/opentargets/AD_finemap/")
theme_set(theme_bw())

#meta_annotated = read_tsv("annotated/AD.meta.annotated.selected.probable.paintor.tsv") %>%
#  dplyr::select(locus_name, rsids)

locusFile = "AD.loci.exceptAPOE.tsv"
loci.df = readr::read_tsv(locusFile) %>%
  mutate(lead_chrpos = paste(Chr, lead_pos, sep = "_"))

# Load FINEMAP probabilities
getFinemapSnpID = function(snpid_a1_a2) {
  snpid_parts = str_split(snpid_a1_a2, "_", Inf)
  sapply(snpid_parts, FUN = function(s) paste(s[1:(length(s)-2)], collapse="_"))
}

# Load FINEMAP results for the meta-analysis
finemapFile1 = "finemap/output/AD.meta.finemap.ncausal_1.snp"
finemapFile2 = "finemap/output/AD.meta.finemap.ncausal_2.snp"
finemapFile3 = "finemap/output/AD.meta.finemap.ncausal_3.snp"

meta.ncausal1.df = readr::read_tsv(finemapFile1) %>%
  select(snp=rsid, meta_prob1=prob, everything()) %>%
  mutate(snp=getFinemapSnpID(snp))

meta.ncausal2.df = readr::read_tsv(finemapFile2) %>%
  select(snp=rsid, meta_prob2=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

meta.ncausal3.df = readr::read_tsv(finemapFile3) %>%
  select(snp=rsid, meta_prob3=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

meta_finemap = meta.ncausal1.df %>%
  left_join(meta.ncausal2.df, by="snp") %>%
  left_join(meta.ncausal3.df, by="snp")

# getFinemapProb = function(df) {
#   nCausal = df$locus_nSnps[1]
#   if (nCausal == 1) {
#     df$finemap_prob1
#   } else if (nCausal == 2) {
#     df$finemap_prob2
#   } else if (nCausal == 3) {
#     df$finemap_prob3
#   } else if (nCausal == 4) {
#     df$finemap_prob4
#   } else {
#     rep(NA, nrow(df))
#   }
# }
# meta_finemap$finemap_prob_nc = unlist(by(meta_finemap, meta_finemap$lead_chrpos, getFinemapProb))


# Load FINEMAP results for Kunkle et al
finemapFile1 = "finemap_kunkle/output/AD.meta.finemap.ncausal_1.snp"
finemapFile2 = "finemap_kunkle/output/AD.meta.finemap.ncausal_2.snp"
finemapFile3 = "finemap_kunkle/output/AD.meta.finemap.ncausal_3.snp"

kunkle.ncausal1.df = readr::read_tsv(finemapFile1) %>%
  select(snp=rsid, kunkle_prob1=prob, everything()) %>%
  mutate(snp=getFinemapSnpID(snp))

kunkle.ncausal2.df = readr::read_tsv(finemapFile2) %>%
  select(snp=rsid, kunkle_prob2=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

kunkle.ncausal3.df = readr::read_tsv(finemapFile3) %>%
  select(snp=rsid, kunkle_prob3=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

kunkle_finemap = kunkle.ncausal1.df %>%
  left_join(kunkle.ncausal2.df, by="snp") %>%
  left_join(kunkle.ncausal3.df, by="snp")

merged_finemap = meta_finemap %>%
  inner_join(kunkle_finemap %>% select(snp, kunkle_prob1, kunkle_prob2, kunkle_prob3),
             by = "snp")

merged_finemap = merged_finemap %>%
  mutate(meta_prob1_cat = if_else(meta_prob1 < 0.0001, "< 0.0001", "> 0.0001"),
         kunkle_prob1_cat = if_else(kunkle_prob1 < 0.0001, "< 0.0001", "> 0.0001"),
         meta_prob2_cat = if_else(meta_prob2 < 0.0001, "< 0.0001", "> 0.0001"),
         kunkle_prob2_cat = if_else(kunkle_prob2 < 0.0001, "< 0.0001", "> 0.0001"))

xtabs(~ meta_prob1_cat + kunkle_prob1_cat, data = merged_finemap)
xtabs(~ meta_prob2_cat + kunkle_prob2_cat, data = merged_finemap)

# Annotate which loci SNPs are at
merged_finemap$locus = NA
for (i in 1:nrow(loci.df)) {
  chr = loci.df$Chr[i]
  start = loci.df$lead_pos[i] - 1e6
  end = loci.df$lead_pos[i] + 1e6
  name = loci.df$locus_name[i]
  merged_finemap$locus[merged_finemap$chromosome == loci.df$Chr[i] & between(merged_finemap$position, start, end)] = name
}

merged_finemap = merged_finemap %>%
  left_join(loci.df %>% select(locus_name, lead_p, n_snps), by=c("locus" = "locus_name"))

merged_finemap = merged_finemap %>%
  rowwise() %>%
  mutate(meta_prob_nc = if_else(n_snps == 1, meta_prob1, if_else(n_snps == 2, meta_prob2, meta_prob3)),
         kunkle_prob_nc = if_else(n_snps == 1, kunkle_prob1, if_else(n_snps == 2, kunkle_prob2, kunkle_prob3)))

# Plot SNP probabilities for Kunkle vs. the Meta for each locus

merged_finemap_flt = merged_finemap %>%
  filter(meta_prob1 > 0.001 | kunkle_prob1 > 0.001)
merged_finemap_flt2 = merged_finemap %>%
  filter(meta_prob2 > 0.001 | kunkle_prob2 > 0.001)
merged_finemap_flt_nc = merged_finemap %>%
  filter(meta_prob_nc > 0.001 | kunkle_prob_nc > 0.001)

# For ncausal == 1
p = ggplot(merged_finemap_flt %>%
             mutate(meta_prob1 = if_else(meta_prob1 < 1e-3, 1e-3, meta_prob1),
                    kunkle_prob1 = if_else(kunkle_prob1 < 1e-3, 1e-3, kunkle_prob1)),
           aes(x=meta_prob1, y=kunkle_prob1)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ locus, ncol = 6) +
  theme_bw(13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  xlab("Meta-analysis FINEMAP probability") + ylab("Kunkle FINEMAP prabability")
pdf("replication/finemap_probabiity.by_locus.ncausal_1.pdf", width=10, height=10)
print(p)
dev.off()

# For ncausal == 2
p = ggplot(merged_finemap_flt2 %>%
             mutate(meta_prob2 = if_else(meta_prob2 < 1e-3, 1e-3, meta_prob2),
                    kunkle_prob2 = if_else(kunkle_prob2 < 1e-3, 1e-3, kunkle_prob2)),
           aes(x=meta_prob2, y=kunkle_prob2)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ locus, ncol = 6) +
  theme_bw(13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  xlab("Meta-analysis FINEMAP probability") + ylab("Kunkle FINEMAP prabability")
pdf("replication/finemap_probabiity.by_locus.ncausal_2.pdf", width=10, height=10)
print(p)
dev.off()

# For ncausal == GCTA ncausal from meta-analysis
p = ggplot(merged_finemap_flt_nc %>%
             mutate(meta_prob_nc = if_else(meta_prob_nc < 1e-3, 1e-3, meta_prob_nc),
                    kunkle_prob_nc = if_else(kunkle_prob_nc < 1e-3, 1e-3, kunkle_prob_nc)),
           aes(x=meta_prob_nc, y=kunkle_prob_nc)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ locus, ncol = 6) +
  theme_bw(13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  xlab("Meta-analysis FINEMAP probability") + ylab("Kunkle FINEMAP prabability")
pdf("replication/finemap_probabiity.by_locus.ncausal_gcta_nc.pdf", width=10, height=10)
print(p)
dev.off()


# Read in FINEMAP credible set files to compare overlap between the two
meta_cred_all = NULL
for (i in 1:nrow(loci.df)) {
  locus = paste(loci.df$Chr[i], loci.df$lead_pos[i], sep="_")
  locus_name = loci.df$locus_name[i]
  if (loci.df$n_snps[i] == 1) {
    meta_cred_in = read_delim(sprintf("finemap/out_ncausal_1/%s.IGAP1_GWAX.1.cred", locus), delim=" ", col_names = T)
    colnames(meta_cred_in)[2:3] = c("snp", "prob")
    meta_cred_all = bind_rows(meta_cred_all, cbind(locus_name = locus_name, meta_cred_in[, 2:3]))
  } else if (loci.df$n_snps[i] == 2) {
    meta_cred_in = read_delim(sprintf("finemap/out_ncausal_2/%s.IGAP1_GWAX.2.cred", locus), delim=" ", col_names = T)
    colnames(meta_cred_in)[2:5] = c("snp", "prob", "snp", "prob")
    meta_cred_all = bind_rows(meta_cred_all, cbind(locus_name = locus_name, bind_rows(meta_cred_in[, 2:3], meta_cred_in[, 4:5])))
  } else if (loci.df$n_snps[i] == 3) {
    meta_cred_in = read_delim(sprintf("finemap/out_ncausal_3/%s.IGAP1_GWAX.3.cred", locus), delim=" ", col_names = T)
    colnames(meta_cred_in)[2:7] = c("snp", "prob", "snp", "prob", "snp", "prob")
    meta_cred_all = bind_rows(meta_cred_all, cbind(locus_name = locus_name, bind_rows(meta_cred_in[, 2:3], meta_cred_in[, 4:5], meta_cred_in[, 6:7])))
  }
}


kunkle_cred_all = NULL
for (i in 1:nrow(loci.df)) {
  locus = paste(loci.df$Chr[i], loci.df$lead_pos[i], sep="_")
  locus_name = loci.df$locus_name[i]
  if (loci.df$n_snps[i] == 1) {
    kunkle_cred_in = read_delim(sprintf("finemap_kunkle/out_ncausal_1/%s.Kunkle.1.cred", locus), delim=" ", col_names = T)
  } else if (loci.df$n_snps[i] == 2) {
    kunkle_cred_in = read_delim(sprintf("finemap_kunkle/out_ncausal_2/%s.Kunkle.2.cred", locus), delim=" ", col_names = T)
  } else if (loci.df$n_snps[i] == 3) {
    kunkle_cred_in = read_delim(sprintf("finemap_kunkle/out_ncausal_3/%s.Kunkle.3.cred", locus), delim=" ", col_names = T)
  }
  print(sprintf("%s, N=%d", locus_name, nrow(kunkle_cred_in)))
  colnames(kunkle_cred_in)[2:3] = c("snp", "prob")
  kunkle_cred_all = bind_rows(kunkle_cred_all, cbind(locus_name = locus_name, kunkle_cred_in[, 2:3]))
  if (ncol(kunkle_cred_in) > 4) {
    colnames(kunkle_cred_in)[4:5] = c("snp", "prob")
    kunkle_cred_all = bind_rows(kunkle_cred_all, cbind(locus_name = locus_name, kunkle_cred_in[, 4:5]))
  }
  if (ncol(kunkle_cred_in) > 6) {
    colnames(kunkle_cred_in)[6:7] = c("snp", "prob")
    meta_cred_all = bind_rows(meta_cred_all, cbind(locus_name = locus_name, meta_cred_in[, 6:7]))
  }
}


meta_cred_all = meta_cred_all %>%
  filter(!is.na(snp)) %>%
  mutate(in_kunkle_credset = (snp %in% kunkle_cred_all$snp))
meta_cred_summary = meta_cred_all %>%
  group_by(locus_name) %>%
  summarise(meta_cred_set_size = n(),
            meta_cred_set_overlap = sum(in_kunkle_credset),
            meta_cred_set_prob_overlap = sum(prob[in_kunkle_credset], na.rm=T))

kunkle_cred_all = kunkle_cred_all %>%
  filter(!is.na(snp)) %>%
  mutate(in_meta_credset = (snp %in% meta_cred_all$snp))
kunkle_cred_summary = kunkle_cred_all %>%
  group_by(locus_name) %>%
  summarise(kunkle_cred_set_size = n(),
            kunkle_cred_set_overlap = sum(in_meta_credset),
            kunkle_cred_set_prob_overlap = sum(prob[in_meta_credset], na.rm=T))

credset_summary = meta_cred_summary %>%
  left_join(kunkle_cred_summary, by="locus_name")

write_tsv(credset_summary, path = "replication/credset_kunkle_comparison.tsv")


