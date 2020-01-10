library(tidyverse)
library(optparse)
library(pheatmap)

option_list <- list(
  make_option(c("--file"), type="character", default=NULL,
              help="Main file of annotations"),
  make_option(c("--paintor_file"), type="character", default=NULL,
              help="File of paintor output merged together for all loci")
)
opt <- parse_args(OptionParser(option_list=option_list))
# opt <- list(file = "annotated/AD.meta.annotated.selected.probable.tsv", paintor_file = "paintor_cred/out_nc2/model3/all_loci.results")

annot.df = readr::read_tsv(opt$file, col_names = T, col_types = cols(.default = col_character(), finemap_prob_nc="d", gcta_maxCondProb="d"))

paintor.df = readr::read_tsv(opt$paintor_file, col_names = T, col_types = cols(.default = col_character(), Posterior_Prob="d"))
paintor.df$locus = gsub(".IGAP1_GWAX", "", paintor.df$locus, fixed = T)

# We have to match on SNP ID without alleles, since in the annot.df we don't
# have the alleles in the same order as in the SNP ID in the PAINTOR FILE.
getSnp = function(id) {
  s = str_split(id, "_")[[1]]
  toIndex = max(1, length(s) - 2)
  paste(s[1:toIndex], collapse="_")
}
paintor.df$snpkey = sapply(paintor.df$rsid, getSnp)
#sum(paintor.df$snpkey %in% annot.df$snp)

annot.df2 = annot.df %>%
  left_join(paintor.df %>% select(locus, paintor_pp=Posterior_Prob, snpkey), by=c("lead_chrpos"="locus", "snp"="snpkey")) %>%
  group_by(lead_chrpos, snp) %>%
  mutate(mean_prob = mean(c(finemap_prob_nc, gcta_maxCondProb, paintor_pp), na.rm=T)) %>%
  select(lead_chrpos:locus_nSnps, gcta_pCond, gcta_bCond, gcta_bCondSE, gcta_cond_snps, gcta_maxCondProb, mean_prob, paintor_pp, everything())
#sum(is.na(annot.df2$Posterior_Prob))
# Only the SNPs at locus 4_11027619	(CLNK) don't have a paintor
# probability - those were excluded because the model didn't work
# when the locus was present (for unknown reasons).

write.table(annot.df2, file="", quote=F, sep="\t", col.names=T, row.names=F, na="")
