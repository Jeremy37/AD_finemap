library(tidyverse)

setwd("/Users/jeremys/work/opentargets/AD_finemap/")
theme_set(theme_bw())

# Read in just a few columns from the full dataset. This will still crank
# for a while before finishing.
meta = read_tsv("replication/global_meta/global_meta.assoc_loci.tsv.gz")


locus_summary = meta %>%
  group_by(locus) %>%
  arrange(global_p) %>%
  summarise(SNP = first(SNP),
            global_p = min(global_p, na.rm = T),
            kfu_p = min(kfu_p, na.rm = T),
            Kunkle_P = min(Kunkle_P, na.rm = T),
            META_P = min(META_P, na.rm = T),
            GWAX_UKBB_P = min(GWAX_UKBB_P, na.rm = T),
            finngen_P = min(finngen_P, na.rm = T))

write_tsv(locus_summary, path = "replication/locus_summary_gobal_gwsig.tsv")

plot_df = locus_summary %>%
  dplyr::select(locus, global=global_p, discovery=META_P) %>%
  mutate(log_p_diff = log10(discovery) - log10(global)) %>%
  gather(key=meta_analysis, value=p_val, -locus, -log_p_diff) %>%
  mutate(neg_log_p = -log10(p_val)) %>%
  filter(locus != "APOE")

p = ggplot(plot_df, aes(x=neg_log_p, y=fct_reorder(locus, neg_log_p), fill=meta_analysis)) +
  geom_bar(stat="identity", position="dodge", width=0.75) +
  xlab("-log10(min p value)") + ylab("Locus") +
  geom_vline(xintercept=-log10(5e-8), col="blue", alpha=0.6, lty=5) +
  scale_fill_manual(values = c("global"="cornflowerblue", "discovery"="grey80"))

pdf("replication/global_meta_vs_discovery.pdf", width=7, height=7)
p
dev.off()



###############################################################################
## OLD CODE

meta = read_tsv(gzfile("summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz"),
                col_types = cols_only(SNP = "c", GWAS_BETA = "i", GWAX_UKBB_BETA = "c", META_P = "d", kfu_p = "d", META_P = "d"))

orig_loci = read_tsv("AD.loci.tsv") %>%
  rowwise() %>%
  mutate(snp_id = str_split(lead_SNP, "_")[[1]][1])

meta_orig_loci = orig_loci %>%
  left_join(meta, by=c("snp_id" = "SNP")) %>%
  dplyr::select(snp_id, global_p, kfu_p, META_P, locus_name, lead_freq, everything())

write_tsv(meta_orig_loci, "replication/global_meta/meta.orig_loci.tsv")

ad_annotated = read_tsv("annotated/AD.meta.annotated.selected.probable.paintor.tsv") %>%
  dplyr::select(locus_name, rsids)
#, chr, pos_hg37)

meta_gwsig = meta %>% filter(global_p < 5e-8) %>%
  left_join(ad_annotated %>% dplyr::select(rsids, locus_name), by=c("SNP"="rsids")) %>%
  mutate(pos_diff = pos_hg19 - lag(pos_hg19))
meta_gwsig$pos_diff[meta_gwsig$CHR != lag(meta_gwsig$CHR)] = NA
write_tsv(meta_gwsig, "replication/global_meta/global_meta.p_lt_5e-8.tsv")

meta_gwsig = meta %>% filter(kfu_p < 5e-8) %>%
  left_join(ad_annotated %>% dplyr::select(rsids, locus_name), by=c("SNP"="rsids")) %>%
  mutate(pos_diff = pos_hg19 - lag(pos_hg19))
meta_gwsig$pos_diff[meta_gwsig$CHR != lag(meta_gwsig$CHR)] = NA
write_tsv(meta_gwsig, "replication/global_meta/global_meta.kfu_p_lt_5e-8.tsv")

meta_gwsig_loci = meta_gwsig %>%
  group_by(CHR, locus_name) %>%
  arrange(global_p) %>%
  summarise(SNP = first(SNP),
            global_p = min(global_p),
            kfu_p = min(kfu_p),
            META_P = min(META_P))
write_tsv(meta_gwsig_loci, "replication/global_meta/global_meta.kfu.gwsig_loci.tsv")

meta_sig = meta %>% filter(global_p < 5e-7) %>%
  left_join(ad_annotated %>% dplyr::select(rsids, locus_name), by=c("SNP"="rsids")) %>%
  mutate(pos_diff = pos_hg19 - lag(pos_hg19))
meta_sig$pos_diff[meta_gwsig$CHR != lag(meta_gwsig$CHR)] = NA
write_tsv(meta_sig, "replication/global_meta/global_meta.p_lt_5e-7.tsv")


# Look at change in -log10(p) relative to our original meta-analysis
meta_all_loci = meta %>%
  inner_join(ad_annotated %>% dplyr::select(rsids, locus_name), by=c("SNP"="rsids")) %>%
  mutate(pos_diff = pos_hg19 - lag(pos_hg19))
meta_all_loci$pos_diff[meta_all_loci$CHR != lag(meta_all_loci$CHR)] = NA
meta_all_loci = meta_all_loci %>%
  group_by(CHR, locus_name) %>%
  arrange(global_p) %>%
  summarise(SNP = first(SNP),
            global_p = min(global_p),
            kfu_p = min(kfu_p),
            META_P = min(META_P))


meta_all_loci = meta_all_loci %>%
  mutate(diff_global_meta = -log10(global_p) + log10(META_P),
         diff_kfu_meta = -log10(kfu_p) + log10(META_P))
ggplot(meta_all_loci %>% filter(!is.na(locus_name)), aes(x=diff_global_meta, y=fct_rev(locus_name))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggplot(meta_all_loci %>% filter(!is.na(locus_name)), aes(x=diff_kfu_meta, y=fct_rev(locus_name))) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))




