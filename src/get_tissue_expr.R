#!/usr/bin/env Rscript
# Gathers expression data from multiple datasets to create a large table of
# median TPMs per gene & dataset.
library(tidyverse)

root = ".."
outdir = "./expression"
combined.df = NULL

getTPMs = function(counts, gene_lengths) {
  rpK = apply(counts, MARGIN=2, FUN=function(x) x / gene_lengths) * 1e3
  tpm.mat = apply(rpK, MARGIN=2, FUN=function(x) x / sum(x)) * 1e6
  tpm.mat
}

###########################################################################
# Fiona / Gaffney microglia
microglia.counts.df = readr::read_tsv(file.path(root, "datasets/gaffney_microglia/RNA", "microglia_counts.gz"))
microglia_gene_lengths = microglia.counts.df$Length

microglia.counts.mat = microglia.counts.df %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()
# Calculate TPM
microglia.sample.tpms = getTPMs(microglia.counts.mat, microglia_gene_lengths)

microglia.avg.tpms = data.frame(ensgene_full = rownames(microglia.sample.tpms),
                                ensgene = rownames(microglia.sample.tpms),
                                primary_microglia = apply(microglia.sample.tpms, MARGIN=1, FUN=median))


###########################################################################
# Erica's iPSC-derived microglia
ips_mic.meta.df = readr::read_tsv(file.path(root, "experiment/RNA/CCDC6_clones_new", "CCDC6_clones_all.meta.tsv")) %>%
  filter(cell_type == "microglia", condition == "WT")
ips_mic.counts.df = readr::read_tsv(file.path(root, "experiment/RNA/CCDC6_clones_new", "CCDC6_clones_all.counts.tsv.gz")) %>%
  select(gene_id, length, one_of(ips_mic.meta.df$sample_name))

ips_mic_gene_lengths = ips_mic.counts.df$length

ips_mic.counts.mat = ips_mic.counts.df %>%
  select(-length) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

ips_mic.sample.tpms = getTPMs(ips_mic.counts.mat, ips_mic_gene_lengths)

ips_mic.avg.tpms = data.frame(ensgene_full = rownames(ips_mic.sample.tpms),
                              ensgene = rownames(ips_mic.sample.tpms),
                              ipsc_microglia = apply(ips_mic.sample.tpms, MARGIN=1, FUN=mean)) %>%
  mutate(ensgene = gsub("\\.[\\d]+", "", ensgene_full, perl=T))
  


###########################################################################
# GTEx - provides a table of TPMs

gtex.tpm_file = file.path(root, "reference/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
gtex.tpm.df = readr::read_tsv(gtex.tpm_file, skip = 2) %>%
  dplyr::rename(ensgene_full = Name, symbol = Description) %>%
  group_by(ensgene_full) %>%
  mutate(ensgene = gsub("\\.[\\d]+", "", ensgene_full, perl=T)) %>%
  select(ensgene_full, ensgene, symbol, everything()) %>%
  ungroup()


###########################################################################
# ipsNeurons iPSC-derived cells

ipsneurons.counts.df = readr::read_tsv(file.path(root, "ipsneurons/GRCh38/RNA/analysis/ipsneurons.counts.txt.gz")) %>%
  dplyr::rename(ensgene_full = gene_id)
ipsneurons.gene_lengths = ipsneurons.counts.df$length
ipsneurons.counts.mat = ipsneurons.counts.df %>% column_to_rownames(var = "ensgene_full") %>% select(-length) %>% as.matrix()
ipsneurons.sample.tpms = getTPMs(ipsneurons.counts.mat, ipsneurons.gene_lengths) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensgene_full")

ipsneurons.meta = readr::read_tsv(file.path(root, "ipsneurons/GRCh38/RNA/meta/ipsneurons.rna.metadata.txt")) %>%
  filter(sampleID %in% colnames(ipsneurons.sample.tpms))
celltypes = unique(ipsneurons.meta$celltype)
ipsneurons.tpms.df = data.frame(ensgene_full = ipsneurons.sample.tpms$ensgene_full)
for (celltype in celltypes) {
  sampleIDs = ipsneurons.meta$sampleID[ipsneurons.meta$celltype == celltype]
  group.tpm = ipsneurons.sample.tpms %>% dplyr::select(one_of(sampleIDs))
  ipsneurons.tpms.df = cbind(ipsneurons.tpms.df, rowMeans(as.matrix(group.tpm)))
}
colnames(ipsneurons.tpms.df) = c("ensgene_full", celltypes)
ipsneurons.tpms.df = ipsneurons.tpms.df%>%
  mutate(ensgene = gsub("\\.[\\d]+", "", ensgene_full, perl=T)) %>%
  select(ensgene_full, ensgene, everything())


###########################################################################
# eQTL catalogue provides tables of TPMs
catalogue_names = c("Alasoo_2018", "BLUEPRINT_PE", "BLUEPRINT_SE", "BrainSeq", "GENCORD", "GEUVADIS", "HipSci", "Lepik_2017", "Nedelec_2016", "Quach_2016", "ROSMAP", "Schmiedel_2018", "Schwartzentruber_2018", "TwinsUK", "van_de_Bunt_2015")
dflist = list()
for (dataset in catalogue_names) {
  fname = paste0(dataset, "_median_tpm.tsv.gz")
  df = readr::read_tsv(file.path(root, "datasets/eqtl_catalogue/median_tpm", fname))
  dflist = c(dflist, list(df))
}

eqtl_catalogue.df = bind_rows(dflist) %>%
  rename(ensgene = phenotype_id)

eqtl_catalogue.tpms = eqtl_catalogue.df %>%
  mutate(qtl_study = paste(study, qtl_group)) %>%
  select(-study, -qtl_group) %>%
  spread(key = qtl_study, value = median_tpm)



###########################################################################
# Combine them all together
combined.df = ipsneurons.tpms.df %>% select(-ensgene_full) %>%
  dplyr::full_join(microglia.avg.tpms %>% select(-ensgene_full), by="ensgene") %>%
  dplyr::full_join(ips_mic.avg.tpms %>% select(-ensgene_full), by="ensgene") %>%
  dplyr::full_join(eqtl_catalogue.tpms, by="ensgene") %>%
  dplyr::full_join(gtex.tpm.df %>% select(-ensgene_full, -symbol), by="ensgene")

# Now annotate the full combined table with HGNC gene symbols (and their previous symbols)
library(annotables)
grch38_nodup = grch38 %>%
  select(ensgene, symbol) %>% unique()
combined.df = combined.df %>%
  left_join(grch38_nodup, by="ensgene") %>%
  select(ensgene, symbol, everything())

gzf = gzfile(file.path(outdir, "tissues.combined.tpm.tsv.gz"), "w")
write.table(combined.df, file=gzf, row.names=F, col.names=T, quote=F, sep="\t", na="")
close(gzf)

subset.df = combined.df %>%
  select(ensgene, symbol, `HipSci iPSC`, iNeuron, NPC, neuron, primary_microglia, ipsc_microglia,
         `BLUEPRINT monocyte`, `Alasoo_2018 macrophage_naive`, `Alasoo_2018 macrophage_IFNg`,
         ROSMAP_brain = `ROSMAP brain_naive`, `BrainSeq brain`, GTEx_hippocampus=`Brain - Hippocampus`,
         `Brain - Cerebellum`, `Brain - Cortex`, `Whole Blood`)
gzf = gzfile(file.path(outdir, "tissues.selected.tpm.tsv.gz"), "w")
write.table(subset.df, file=gzf, row.names=F, col.names=T, quote=F, sep="\t", na="")
close(gzf)


# Round values to 3 decimal places
roundToString = function(x) {
  if (is.na(x)) { return("") }
  sprintf("%.3g", round(x, digits = 3))
}
subset.mat = subset.df %>% select(-ensgene, -symbol)
subset.short = apply(subset.mat, MARGIN=c(1, 2), FUN=roundToString)
subset.short.df = cbind(subset.df %>% select(ensgene, symbol),
                        subset.short)
gzf = gzfile(file.path(outdir, "tissues.selected.tpm.short.tsv.gz"), "w")
write.table(subset.short.df, file=gzf, row.names=F, col.names=T, quote=F, sep="\t", na="")
close(gzf)

library(annotables)
subset.short.pc.df = subset.short.df %>%
  left_join(grch38 %>% select(ensgene, biotype)) %>%
  filter(biotype == "protein_coding") %>%
  select(-biotype)
gzf = gzfile(file.path(outdir, "tissues.selected.tpm.protein_coding.short.tsv.gz"), "w")
write.table(subset.short.pc.df, file=gzf, row.names=F, col.names=T, quote=F, sep="\t", na="")
close(gzf)

