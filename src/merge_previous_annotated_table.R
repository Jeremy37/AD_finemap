library(tidyverse)
options(stringsAsFactors = F)


annot.df = readr::read_tsv("annotated/AD.meta.annotated.selected.probable.paintor.tsv", col_types = cols(.default = col_character()))

old_annot.df = readr::read_tsv("annotated/AD.meta_v5.annotated.selected.probable.with_notes.tsv", col_types = cols(.default = col_character()), guess_max = 100000)

which(!colnames(annot.df) %in% colnames(old_annot.df))
colnames(annot.df)[!colnames(annot.df) %in% colnames(old_annot.df)]

annot.df2 = annot.df %>%
  left_join(old_annot.df %>% select(snp, `Top candidate variant`, NOTES), by="snp") %>%
  select(lead_chrpos:finemap_prob3, `Top candidate variant`, NOTES, everything())

write.table(annot.df2, file="annotated/AD.meta.annotated.selected.probable.paintor.with_notes.tsv", quote=F, sep="\t", col.names=T, row.names=F, na="")

