library(tidyverse)

snps.df = read_tsv("/Users/jeremys/work/opentargets/AD_finemap/finemap/input/all.snplist.txt", col_names = F)
excl.df = read_tsv("/Users/jeremys/work/opentargets/AD_finemap/finemap/input/all_excluded_snps.txt")

excl.df$exclusion_reason2 = excl.df$exclusion_reason
excl.df$exclusion_reason2[grepl("SNP has freq", excl.df$exclusion_reason2)] = "SNP freq"
excl.df$exclusion_reason2[grepl("SNP info", excl.df$exclusion_reason2)] = "SNP info"
excl.df$exclusion_reason2[grepl("het_p value", excl.df$exclusion_reason2)] = "het_p value"

tbl = table(excl.df$exclusion_reason2)
tbl

n_total_variants = nrow(snps.df) + nrow(excl.df)

sprintf("Excluded for SNP freq: %d of %d (%.3f%%)", tbl["SNP freq"], n_total_variants, 100 * tbl["SNP freq"] / n_total_variants)
sprintf("Excluded for SNP info: %d of %d (%.3f%%)", tbl["SNP info"], n_total_variants, 100 * tbl["SNP info"] / n_total_variants)
sprintf("Excluded for SNP het_p: %d of %d (%.3f%%)", tbl["het_p value"], n_total_variants, 100 * tbl["het_p value"] / n_total_variants)

sprintf("All exclusions: %d of %d (%.3f%%)", nrow(excl.df), n_total_variants, 100 * nrow(excl.df) / n_total_variants)

