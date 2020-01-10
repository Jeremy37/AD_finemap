#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
leadSnpsFile = args[1]
allSnpsFile = args[2]
spliceAIfname = args[3]
outputFname = args[4]

leadSnps.df = readr::read_tsv(leadSnpsFile)
snps.df = readr::read_tsv(allSnpsFile)

# Processing the spliceAI datasets is slow, so cache the file
# in case we need to redo the SNP annotation
if (file.exists("spliceai.tmp.tsv.gz")) {
  spliceai.df = readr::read_tsv("spliceai.tmp.tsv.gz")
  dim(spliceai.df)
} else {
  spliceaiColnames = c("chr", "pos", "ref", "alt", "symbol", "strand", "type", "dist", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
  spliceai.dflist = list()
  for (i in 1:nrow(leadSnps.df)) {
    chr = leadSnps.df$Chr[i]
    minpos = max(1, leadSnps.df$start[i] - 1e6)
    maxpos = leadSnps.df$stop[i] + 1e6
    cmd = sprintf("tabix %s %s:%d-%d", spliceAIfname, chr, minpos, maxpos, i)
    print(cmd)
    spliceAI_data_str = system(cmd, intern = T)
    if (length(spliceAI_data_str) > 1) {
      spliceAI_data_str_single = paste0(spliceAI_data_str, collapse = "\n")
      spliceai.dflist = c(spliceai.dflist, list(readr::read_tsv(spliceAI_data_str_single, col_names = spliceaiColnames)))
    } else {
      stop(sprintf("Error: no data received from tabix command to get SpliceAI data for %s:%d-%d", chr, minpos, maxpos))
    }
  }
  
  spliceai.df1 = bind_rows(spliceai.dflist)
  
  # Make a second copy of the spliceAI table, but swapping ref/alt alleles,
  # as well as AG/AL, and DG/DL columns. This is because in our annotated SNP file
  # we have only an "effect allele", but don't know which is ref/alt. When we then
  # join with this table, we can retrieve the spliceai score for the relevant change.
  spliceai.df2 = spliceai.df1 %>%
    dplyr::rename(ref2 = ref, ref = alt, DS_AG2 = DS_AG, DS_AG = DS_AL, DS_DG2 = DS_DG, DS_DG = DS_DL) %>%
    dplyr::rename(alt = ref2, DS_AL = DS_AG2, DS_DL = DS_DG2)
  
  spliceai.df = bind_rows(spliceai.df1, spliceai.df2)
  spliceai.df = spliceai.df %>%
    mutate(max_DS = pmax(DS_AG, DS_AL, DS_DG, DS_DL))
  
  gzf = gzfile("spliceai.tmp.tsv.gz", "w")
  write.table(spliceai.df, file=gzf, sep="\t", quote=F, col.names=T, row.names=F, na="")
  close(gzf)
}

print("Joining SpliceAI table...")
snps.df = snps.df %>%
  select(chr=CHR, pos=BP, SNP, Eff_allele=A1, other_allele=A2) %>%
  left_join(spliceai.df %>% select(everything(), -symbol, -dist, spliceai_strand = strand, spliceai_type = type),
            by=c("chr", "pos", "Eff_allele"="alt", "other_allele"="ref"))

write.table(snps.df, file=outputFname, sep="\t", quote=F, col.names=T, row.names=F, na="")
