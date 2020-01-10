library(tidyverse)
library(annotables)
options(stringsAsFactors = T)

ad_dir = "/Users/jeremys/work/opentargets/AD_finemap"
dataset_groups_file = file.path(ad_dir, "coloc/dataset_groups.tsv")
setwd(ad_dir)

theme_set(theme_bw(8))

coloc_main.df = readr::read_tsv("coloc/output/coloc.AD.meta.eqtl_sqtl_colocs.txt", col_names = T) %>% filter(!is.na(locus_name))

dataset_groups = readr::read_tsv(dataset_groups_file) %>%
  mutate(dataset_short = gsub("Schmiedel_mono_CD16_naive", "Schmiedel_CD16", dataset_short),
         dataset_short = gsub("_mono", "", dataset_short))

# Shorten some dataset labels further
coloc_main.df = coloc_main.df %>%
  mutate(dataset_short = gsub("Schmiedel_mono_CD16_naive", "Schmiedel_naive", dataset_short),
         dataset_short = gsub("_mono", "", dataset_short)) %>%
  left_join(dataset_groups %>% select(dataset_short, celltype_group), by="dataset_short")


colocScores.df = readr::read_tsv("coloc/coloc.AD.meta.colocScores.per_gene.tidy.tsv") %>%
  arrange(desc(colocScore)) %>%
  mutate(colocScoreRank = 1:n())

getSymbol = function(s) {
  df = data.frame(gene = s) %>% left_join(grch38 %>% select(ensgene, symbol), by=c("gene"="ensgene"))
  df$gene[grepl("ENSG", df$gene)] = df$symbol[grepl("ENSG", df$gene)]
  df$gene
}

coloc.all.df = coloc_main.df %>%
  select(locus_name, signal, feature, ensembl_id, geneSymbol, H3=PP.H3, H4, qtl_pval, gwas_pval, qtl_lead, gwas_lead, chr, qtl_pos, dataset, dataset_short, celltype_group) %>%
  mutate(geneSymbol = getSymbol(geneSymbol))

locus_gwas_pos = coloc_main.df %>%
  select(locus_name, gwas_pos) %>%
  group_by(locus_name) %>%
  summarise(gwas_pos = first(gwas_pos))
coloc.all.df = coloc.all.df %>%
  left_join(locus_gwas_pos, by="locus_name")

# Determine the max H4 score across signals, and focus on eQTLs only (not hQTLs)
coloc.eqtl_sqtl.df = coloc.all.df %>%
  filter(!grepl("k27", dataset_short), !grepl("k4me", dataset_short)) %>%
  arrange(desc(H4)) %>%
  filter(!duplicated(paste(dataset_short, geneSymbol)))

coloc.eqtl_sqtl.df = coloc.eqtl_sqtl.df %>%
  group_by(geneSymbol) %>%
  mutate(max_H4 = max(H4),
         num_H4_gt_0.8 = sum(H4 > 0.8, na.rm=T)) %>%
  left_join(colocScores.df, by=c("locus_name", "ensembl_id"))

coloc.gtex.df = coloc.eqtl_sqtl.df %>%
  filter(grepl("gtex", dataset_short)) %>%
  group_by(geneSymbol) %>%
  mutate(max_H4 = max(H4),
         num_H4_gt_0.8 = sum(H4 > 0.8, na.rm=T))

#displayDatasets = c("microglia", "brain_meta", "ROSMAP_brain", "brainseq_brain", "xQTL_eQTL", "blueprint", "Fairfax_LPS24", "Fairfax_LPS2", "Fairfax_IFN24", "Schmiedel_CD16_naive", "Schmiedel_naive", "Quach_Pam3CSK4", "Quach_LPS", "CEDAR", "Quach_R848", "Quach_naive", "Quach_IAV")
displayDatasets = c("brain_meta", "ROSMAP_brain", "brainseq_brain", "xQTL_eQTL", "microglia", "blueprint", "Fairfax_LPS24", "Fairfax_IFN24", "Schmiedel_naive", "Quach_LPS", "Quach_IAV", "CEDAR", "Other")
allDatasets = unique(c(displayDatasets, coloc.eqtl_sqtl.df$dataset_short, "other"))
coloc.eqtl_sqtl.df$dataset_short = factor(as.character(coloc.eqtl_sqtl.df$dataset_short), levels=allDatasets)

coloc.sel.df = coloc.eqtl_sqtl.df %>%
  filter(dataset_short %in% displayDatasets) %>%
  ungroup()
length(unique(coloc.sel.df$geneSymbol))
#filter(max_H4 > 0.5) %>%

# Select the top N genes. When we spread the data out by gene we can select
# on row_number at that point.
nGenesToShow = 40
coloc.sel.filled.df = coloc.sel.df %>%
  select(locus_name, geneSymbol, H4, chr, gwas_pos, dataset_short, colocScore, num_H4_gt_0.8) %>%
  spread(dataset_short, H4, fill = NA) %>%
  arrange(desc(colocScore)) %>% ungroup() %>%
  filter(row_number() <= nGenesToShow) %>%
  gather("dataset_short", "H4", one_of(displayDatasets), na.rm=F) %>%
  mutate(dataset_short = factor(as.character(dataset_short), levels=displayDatasets))
length(unique(coloc.sel.filled.df$geneSymbol))

getStringNoNA = function(x) {
  if (is.na(x)) {
    NA
  } else {
    sprintf("%.2g", round(x, digits = 2))
  }
}
# Add some labels with the H4 probability, formatted with 2 digits
coloc.sel.filled.df = coloc.sel.filled.df %>%
  group_by(geneSymbol, dataset_short) %>%
  mutate(label = getStringNoNA(H4)) %>%
  ungroup() %>%
  arrange(chr, gwas_pos, as.character(geneSymbol)) %>%
  mutate(locus_name = factor(as.character(locus_name), levels = unique(as.character(locus_name)))) %>%
  arrange(as.integer(locus_name)) %>%
  mutate(locus_grp = as.factor(as.integer(locus_name) %% 2)) %>%
  mutate(geneSymbol = factor(as.character(geneSymbol), levels = unique(as.character(geneSymbol))))

# Add back in the celltype group, which we had to remove when we spread/gathered
coloc.sel.filled.df = coloc.sel.filled.df %>%
  left_join(dataset_groups %>% select(dataset_short, celltype_group), by="dataset_short") %>%
  mutate(dataset_short = factor(as.character(dataset_short), levels=allDatasets))

dataset_colors = c("brain" = "orange", "microglia" = "#00BA38", "monocyte" = "#619CFF", "other" = "#F8766D")

p1 = ggplot(coloc.sel.filled.df, aes(x=dataset_short, y=fct_reorder(geneSymbol, -as.integer(geneSymbol)))) +
  geom_raster(aes(fill = locus_grp), alpha=0.3) +
  geom_point(mapping = aes(size=H4^2, col=celltype_group, alpha=(0.5 + (H4-0.2)/3))) +
  geom_text(mapping = aes(label=label), size=1.8, nudge_x = 0.4) +
  scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
  scale_color_manual(values = dataset_colors, guide = F) +
  scale_size_continuous(range = c(0.5, 4), guide=F) +
  scale_alpha_continuous(range = c(0.4, 1), guide=F) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
p1

# Make a plot similar to the above, but which only shows the number of colocs
# across all loci
coloc.summary.plot.df = coloc.sel.filled.df %>%
  select(locus_name, geneSymbol, colocScore, chr, gwas_pos, num_H4_gt_0.8) %>%
  unique() %>%
  ungroup() %>%
  arrange(chr, gwas_pos, as.character(geneSymbol)) %>%
  mutate(locus_name = factor(as.character(locus_name), levels = unique(as.character(locus_name)))) %>%
  arrange(as.integer(locus_name)) %>%
  mutate(locus_grp = as.factor(as.integer(locus_name) %% 2)) %>%
  mutate(geneSymbol = factor(as.character(geneSymbol), levels = unique(as.character(geneSymbol)))) %>%
  group_by(geneSymbol) %>%
  mutate(numH4size = min(5, sqrt(num_H4_gt_0.8 + 1)),
         colocScoreLabel = sprintf("%.2g", colocScore),
         numH4label = as.character(num_H4_gt_0.8))

p2 = ggplot(coloc.summary.plot.df, aes(x="Num colocs", y=fct_reorder(geneSymbol, -as.integer(geneSymbol)))) +
  geom_raster(aes(fill = locus_grp), alpha=0.3) +
  geom_point(mapping = aes(size=numH4size, col="Num colocs"), alpha=0.8) +
  scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
  scale_color_manual(values = c("Num colocs" = "#F8766D"), guide = F) +
  scale_size_continuous(range = c(0.5, 4), guide=F) +
  geom_text(mapping = aes(label=numH4label), size=2.1, nudge_x = 0.25) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
p2

# A bubble plot with just the coloc score per locus. It doesn't look interesting,
# since we've already selected just the top colocs, so we're not including it.
p3 = ggplot(coloc.summary.plot.df, aes(x="Coloc Score", y=fct_reorder(geneSymbol, as.integer(locus_name)))) +
  geom_raster(aes(fill = locus_grp), alpha=0.3) +
  geom_point(mapping = aes(size=colocScore, col="Coloc Score")) +
  scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
  scale_color_manual(values = c("Coloc Score" = "cornflowerblue"), guide = F) +
  scale_size_continuous(range = c(1, 5), guide=F) +
  geom_text(mapping = aes(label=colocScoreLabel), size=2.6, nudge_x = 0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

# coloc.gtex.plot.df = coloc.gtex.plot.df %>%
#   select(locus_name, geneSymbol, num_H4_gt_0.8) %>%
#   mutate(dataset_short = "GTEx") %>% unique() %>%
#   filter(geneSymbol %in% coloc.sel3.df$geneSymbol) %>%
#   mutate(label = as.character(num_H4_gt_0.8)) %>%
#   ungroup() %>%
#   mutate(locus_name = as.factor(locus_name)) %>%
#   arrange(as.integer(locus_name)) %>%
#   mutate(locus_grp = as.factor(as.integer(locus_name) %% 2))
# coloc.gtex.plot.df$num_H4_gt_0.8[coloc.gtex.plot.df$num_H4_gt_0.8 > 10] = 10

pdf(file = file.path(ad_dir, "plots", "coloc.summary.mainSignals.pdf"), width = 6.8, height=5)
cowplot::plot_grid(plotlist = list(p1, p2), align = "hv", rel_widths = c(6, 1.5))
dev.off()

png(filename = file.path(ad_dir, "plots", "coloc.summary.mainSignals.png"), width = 6.8, height=5, res = 300, unit="in")
cowplot::plot_grid(plotlist = list(p1, p2), align = "hv", rel_widths = c(6, 1.5))
dev.off()

########################################################################
# Make a similar plot, but showing each conditional signal independently

coloc.eqtl_sqtl.df = coloc.all.df %>%
  filter(!grepl("k27", dataset_short), !grepl("k4me", dataset_short)) %>%
  arrange(desc(H4)) %>%
  group_by(geneSymbol, signal) %>%
  mutate(geneSignal = ifelse(grepl("cond", signal), paste0(geneSymbol, ".", signal), geneSymbol)) %>%
  filter(!duplicated(paste(dataset_short, geneSignal)))

coloc.eqtl_sqtl.signals.df = coloc.eqtl_sqtl.df %>%
  group_by(geneSignal) %>%
  mutate(max_H4 = max(H4),
         num_H4_gt_0.8 = sum(H4 > 0.8, na.rm=T)) %>%
  left_join(colocScores.df, by=c("locus_name", "geneSymbol"))

coloc.sel.signals.df = coloc.eqtl_sqtl.signals.df %>%
  filter(dataset_short %in% displayDatasets) %>%
  ungroup()
length(unique(coloc.sel.signals.df$geneSignal))

conditional_loci = c("TREM2", "BIN1", "PTK2B-CLU", "EPHA1", "ACE", "ABCA7", "ADAM10", "PLCG2", "NCK2", "APP-ADAMTS1")
colocBubblePlotPerSignal = function(coloc.sel.signals.df, nGenesToShow, loci) {
  # Select the top N genes. When we spread the data out by gene we can select
  # on row_number at that point.
  coloc.sel.filled.df = coloc.sel.signals.df %>%
    filter(locus_name %in% loci) %>%
    select(locus_name, geneSignal, H4, dataset_short, colocScore, num_H4_gt_0.8) %>%
    spread(dataset_short, H4, fill = NA) %>%
    arrange(desc(colocScore)) %>%
    filter(row_number() <= nGenesToShow) %>%
    gather("dataset_short", "H4", one_of(displayDatasets), na.rm=F) %>%
    mutate(dataset_short = factor(as.character(dataset_short), levels=displayDatasets))
  
  # Add some labels with the H4 probability, formatted with 2 digits
  coloc.sel.filled.df = coloc.sel.filled.df %>%
    group_by(geneSignal, dataset_short) %>%
    mutate(label = getStringNoNA(H4)) %>%
    ungroup() %>%
    mutate(locus_name = as.factor(locus_name)) %>%
    arrange(as.integer(locus_name)) %>%
    mutate(locus_grp = as.factor(as.integer(locus_name) %% 2))
  
  p1 = ggplot(coloc.sel.filled.df, aes(x=dataset_short, y=fct_reorder(geneSignal, as.integer(locus_name)))) +
    geom_raster(aes(fill = locus_grp), alpha=0.3) +
    geom_point(mapping = aes(size=H4^2, col=dataset_short, alpha=(0.5 + (H4-0.2)/3))) +
    geom_text(mapping = aes(label=label), size=1.8, nudge_x = 0.4) +
    scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
    scale_color_discrete(guide = F) +
    scale_size_continuous(guide=F) +
    scale_alpha_continuous(guide=F) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  p1
}

pdf(file = file.path(ad_dir, "plots", "coloc.summary.cond_signals.pdf"), width = 6, height=6)
for (locus in conditional_loci) {
  print(colocBubblePlotPerSignal(coloc.sel.signals.df, 100, locus))
}
dev.off()



########################################################################
# Make a plot which includes mostly main gene signals, but conditional
# signals for selected loci

nGenesToShow = 32

selected_conditional_loci = c("BIN1", "PTK2B-CLU", "EPHA1")
selected_conditional_signals = c("BIN1.cond_1", "BIN1.cond_2", "CLU.cond_1", "CLU.cond_2", "EPHA1-AS1.cond_1", "EPHA1-AS1.cond_2", "ZYX.cond_1", "ZYX.cond_2", "PTK2B.cond_1", "PTK2B.cond_2")

coloc.sel.filled.df = coloc.sel.df %>%
  select(locus_name, geneSymbol, chr, gwas_pos, H4, dataset_short, colocScore, num_H4_gt_0.8) %>%
  spread(dataset_short, H4, fill = NA) %>%
  arrange(desc(colocScore)) %>% ungroup() %>%
  filter(!locus_name %in% selected_conditional_loci) %>%
  filter(row_number() <= nGenesToShow) %>%
  gather("dataset_short", "H4", one_of(displayDatasets), na.rm=F) %>%
  mutate(dataset_short = factor(as.character(dataset_short), levels=displayDatasets),
         geneSignal = geneSymbol)
length(unique(coloc.sel.filled.df$geneSymbol))

coloc.sel.filled.df = coloc.sel.filled.df %>%
  group_by(geneSymbol, dataset_short) %>%
  mutate(label = getStringNoNA(H4))

coloc.sel.filled.signals.df = coloc.sel.signals.df %>%
  select(locus_name, geneSignal, chr, gwas_pos, H4, dataset_short, colocScore, num_H4_gt_0.8) %>%
  spread(dataset_short, H4, fill = NA) %>%
  arrange(desc(colocScore)) %>%
  filter(row_number() <= nGenesToShow) %>%
  gather("dataset_short", "H4", one_of(displayDatasets), na.rm=F) %>%
  mutate(dataset_short = factor(as.character(dataset_short), levels=displayDatasets))

# Subset to our few signals/genes of interest
coloc.sel.filled.signals.df = coloc.sel.filled.signals.df %>%
  group_by(geneSignal) %>%
  filter(geneSignal %in% selected_conditional_signals)
#filter(any(sapply(conditional_locus_genes, function(gene) {grepl(gene, geneSignal)})))

coloc.sel.filled.signals.df = coloc.sel.filled.signals.df %>%
  group_by(geneSignal, dataset_short) %>%
  mutate(label = getStringNoNA(H4))


coloc.sel.plot.df = bind_rows(coloc.sel.filled.df, coloc.sel.filled.signals.df) %>%
  ungroup() %>%
  arrange(chr, gwas_pos, as.character(geneSignal)) %>%
  mutate(locus_name = factor(as.character(locus_name), levels = unique(as.character(locus_name)))) %>%
  arrange(as.integer(locus_name)) %>%
  mutate(locus_grp = as.factor(as.integer(locus_name) %% 2)) %>%
  mutate(geneSignal = gsub(".cond_1", " sig1", geneSignal),
         geneSignal = gsub(".cond_2", " sig2", geneSignal),
         geneSignal = factor(as.character(geneSignal), levels = unique(as.character(geneSignal))))

# Add back in the celltype group, which we had to remove when we spread/gathered
coloc.sel.plot.df = coloc.sel.plot.df %>%
  left_join(dataset_groups %>% select(dataset_short, celltype_group), by="dataset_short") %>%
  mutate(dataset_short = factor(as.character(dataset_short), levels=allDatasets))


p1 = ggplot(coloc.sel.plot.df, aes(x=dataset_short, y=fct_reorder(geneSignal, -as.integer(geneSignal)))) +
  geom_raster(aes(fill = locus_grp), alpha=0.3) +
  geom_point(mapping = aes(size=H4^2, col=celltype_group, alpha=(0.5 + (H4-0.2)/3))) +
  geom_text(mapping = aes(label=label), size=1.8, nudge_x = 0.4) +
  scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
  scale_color_manual(values = dataset_colors, guide = F) +
  scale_size_continuous(range = c(0.5, 4), guide=F) +
  scale_alpha_continuous(guide=F) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
p1

coloc.summary.plot.df = coloc.sel.plot.df %>%
  select(locus_name, geneSignal, colocScore, num_H4_gt_0.8) %>%
  unique() %>%
  mutate(locus_name = as.factor(locus_name)) %>%
  arrange(as.integer(locus_name)) %>%
  mutate(locus_grp = as.factor(as.integer(locus_name) %% 2)) %>%
  group_by(geneSignal) %>%
  mutate(numH4size = min(5, sqrt(num_H4_gt_0.8 + 1)),
         colocScoreLabel = sprintf("%.2g", colocScore),
         numH4label = as.character(num_H4_gt_0.8))

p2 = ggplot(coloc.summary.plot.df, aes(x="Num colocs", y=fct_reorder(geneSignal, -as.integer(geneSignal)))) +
  geom_raster(aes(fill = locus_grp), alpha=0.3) +
  geom_point(mapping = aes(size=numH4size, col="Num colocs"), alpha=0.8) +
  scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
  scale_color_manual(values = c("Num colocs" = "#F8766D"), guide = F) +
  scale_size_continuous(range = c(0.5, 4), guide=F) +
  geom_text(mapping = aes(label=numH4label), size=2.1, nudge_x = 0.35) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
p2

pdf(file = file.path(ad_dir, "plots", "coloc.summary.selected_signals.pdf"), width = 6.8, height=5)
cowplot::plot_grid(plotlist = list(p1, p2), align = "hv", rel_widths = c(6, 1.5))
dev.off()



########################################################################
# Plot top coloc scores aggregated by cell type group
scores.celltype.perGene.df = readr::read_tsv("coloc/coloc.AD.meta.cellTypeScores.per_gene.tidy.tsv")
write.table(scores.perGene.df, file=paste0(output_root, ".colocScores.per_gene.tidy.tsv"), quote=F, row.names=F, col.names=T, sep="\t", na="")

nGenesToShow = 40
celltype.plot.df = scores.celltype.perGene.df %>%
  left_join(colocScores.df, by=c("locus_name", "geneSymbol")) %>%
  select(locus_name, geneSymbol, colocScore, everything()) %>%
  spread(celltype_group, cellTypeScore, fill = NA) %>%
  arrange(desc(colocScore)) %>% ungroup() %>%
  filter(row_number() <= nGenesToShow) %>%
  gather("celltype_group", "cellTypeScore", one_of(celltype_groups), na.rm=F) %>%
  mutate(locus_name = as.factor(locus_name)) %>%
  arrange(as.integer(locus_name)) %>%
  mutate(locus_grp = as.factor(as.integer(locus_name) %% 2))
length(unique(celltype.plot.df$geneSymbol))

p.celltypes = ggplot(celltype.plot.df, aes(x=celltype_group, y=fct_reorder(geneSymbol, as.integer(locus_name)))) +
  geom_raster(aes(fill = locus_grp), alpha=0.3) +
  geom_point(mapping = aes(size=-log(1.01-cellTypeScore), col=celltype_group, alpha=(0.5 + (cellTypeScore-0.2)/3))) +
  scale_fill_manual(values = c("1"="white", "0"="grey"), guide = F) +
  scale_color_discrete(guide = F) +
  scale_size_continuous(guide=F) +
  scale_alpha_continuous(guide=F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
p.celltypes

pdf(file = file.path(ad_dir, "plots", "coloc.summary.mainSignals.cellTypeScores.pdf"), width = 6.8, height=5)
p.celltypes
dev.off()


########################################################################
# Correlation of colocs across datasets

colocs.spread.df = readr::read_tsv("coloc/coloc.AD.meta.colocScores.per_gene.gt_0.5.tsv") %>% select(-contains("CEDAR"), -contains("Fairfax"))
dataset.cor = cor(colocs.spread.df %>% select(-locus_name, -geneSymbol, -colocScore), use="pairwise.complete.obs")
library(pheatmap)
pheatmap(dataset.cor)

########################################################################
# Number of colocs

length(unique(coloc.eqtl_sqtl.df$geneSymbol))

nrow(coloc.eqtl_sqtl.df %>% filter(H4 > 0.8) %>% group_by(locus_name, geneSymbol) %>% summarise(H4 = max(H4)))
nrow(coloc.eqtl_sqtl.df %>% filter(H4 > 0.5) %>% group_by(locus_name, geneSymbol) %>% summarise(H4 = max(H4)))

# Get the top coloc per locus across selected datasets
coloc.sel.filled.df = coloc.sel.df %>%
  select(locus_name, geneSymbol, H4, chr, gwas_pos, dataset_short, colocScore, num_H4_gt_0.8) %>%
  spread(dataset_short, H4, fill = NA) %>%
  arrange(desc(colocScore)) %>% ungroup()

# coloc.sel.filled.gather = coloc.sel.filled.df %>%
#   gather("dataset_short", "H4", one_of(displayDatasets), na.rm=F) %>%
#   mutate(dataset_short = factor(as.character(dataset_short), levels=displayDatasets)) %>%
#   filter(!is.na)

coloc.sel.maxH4 = coloc.sel.df %>%
  group_by(locus_name, geneSymbol) %>%
  summarise(maxH4 = max(H4),
            dataset = dataset[which.max(H4)])

coloc.sel.top_per_locus.df = coloc.sel.filled.df %>%
  group_by(locus_name) %>%
  summarise(colocScore = max(colocScore),
            geneSymbol = geneSymbol[which.max(colocScore)]) %>%
  left_join(coloc.sel.maxH4, by=c("locus_name", "geneSymbol"))


# Get the top coloc per locus across all datasets
coloc.filled.df = coloc.eqtl_sqtl.df %>%
  select(locus_name, geneSymbol, H4, chr, gwas_pos, dataset_short, colocScore, num_H4_gt_0.8) %>%
  spread(dataset_short, H4, fill = NA) %>%
  arrange(desc(colocScore)) %>% ungroup() %>%
  filter(!is.na(colocScore))

coloc.maxH4.perGene = coloc.eqtl_sqtl.df %>%
  group_by(locus_name, geneSymbol) %>%
  summarise(maxH4 = max(H4),
            dataset = dataset[which.max(H4)])

coloc.maxH4.perLocus = coloc.maxH4.perGene %>%
  group_by(locus_name) %>%
  summarise(geneSymbol = geneSymbol[which.max(maxH4)],
            dataset = dataset[which.max(maxH4)],
            maxH4 = max(maxH4))

coloc.maxColocScore.perLocus = coloc.eqtl_sqtl.df %>%
  group_by(locus_name) %>%
  summarise(geneSymbol = geneSymbol[which.max(colocScore)],
            colocScore = max(colocScore, na.rm=T))

coloc.maxH4.perLocus.gtexOnly = coloc.eqtl_sqtl.df %>%
  filter(grepl("gtex", dataset)) %>%
  group_by(locus_name, geneSymbol) %>%
  summarise(maxH4 = max(H4),
            dataset = dataset[which.max(H4)]) %>%
  group_by(locus_name) %>%
  summarise(geneSymbol = geneSymbol[which.max(maxH4)],
            dataset = dataset[which.max(maxH4)],
            maxH4 = max(maxH4))

coloc.maxH4.perLocus$isMaxColocScore = coloc.maxH4.perLocus$geneSymbol %in% coloc.maxColocScore.perLocus$geneSymbol
coloc.maxColocScore.perLocus$isMaxH4 = coloc.maxColocScore.perLocus$geneSymbol %in% coloc.maxH4.perLocus$geneSymbol
coloc.maxH4.perLocus.gtexOnly$isMaxColocScore = coloc.maxH4.perLocus.gtexOnly$geneSymbol %in% coloc.maxColocScore.perLocus$geneSymbol


coloc.all.top_per_locus.df = coloc.filled.df %>%
  group_by(locus_name) %>%
  summarise(colocScore = max(colocScore),
            geneSymbol = geneSymbol[which.max(colocScore)]) %>%
  left_join(coloc.maxH4, by=c("locus_name", "geneSymbol"))

coloc.all.top_H4per_locus.df = coloc.filled.df %>%
  left_join(coloc.maxH4, by=c("locus_name", "geneSymbol")) %>%
  group_by(locus_name) %>%
  summarise(geneMaxH4 = maxH4[which.max(maxH4)],
            dataset_maxColoc = dataset[which.max(colocScore)],
            dataset_maxH4 = dataset[which.max(maxH4)],
            colocScore = colocScore[which.max(maxH4)],
            geneSymbol_maxColoc = geneSymbol[which.max(colocScore)],
            geneSymbol_maxH4 = geneSymbol[which.max(maxH4)]) %>%
  select(locus_name, colocScore, geneMaxH4, geneSymbol_maxColoc, geneSymbol_maxH4, dataset_maxColoc, dataset_maxH4)





###############################################################################
# Get max coloc score or max H4 per locus, across different data subsets

colocs.perGene.df = coloc.eqtl_sqtl.df %>%
  group_by(locus_name, ensembl_id, geneSymbol) %>%
  summarise(maxInd = which.max(H4),
            signal = signal[maxInd],
            maxH4 = max(H4),
            dataset = dataset[maxInd],
            dataset_short = dataset_short[maxInd]) %>%
  select(-maxInd)

colocs.perGene.gtex.df = coloc.eqtl_sqtl.df %>%
  filter(grepl("gtex", dataset)) %>%
  group_by(locus_name, ensembl_id, geneSymbol) %>%
  summarise(maxInd = which.max(H4),
            gtex_signal = signal[maxInd],
            gtex_maxH4 = max(H4),
            gtex_dataset = dataset[maxInd],
            gtex_dataset_short = dataset_short[maxInd]) %>%
  select(-maxInd)

coloc.summary.perGene = colocs.perGene.df %>%
  left_join(colocs.perGene.gtex.df, by=c("locus_name", "ensembl_id", "geneSymbol")) %>%
  left_join(scores.perGene.df, by=c("locus_name", "ensembl_id")) %>%
  select(locus_name, ensembl_id, geneSymbol, colocScore, everything()) %>%
  arrange(desc(colocScore)) %>% ungroup() %>%
  mutate(colocScore_rank = row_number()) %>%
  arrange(desc(maxH4)) %>%
  mutate(maxH4_rank = row_number()) %>%
  arrange(desc(gtex_maxH4)) %>%
  mutate(gtex_maxH4_rank = row_number()) %>%
  arrange(desc(colocScore)) %>%
  filter(!is.na(ensembl_id))

# Add in merged gene evidence score
mergedEvidence.df = readr::read_tsv(file.path(root, "genes/AD.loci.1Mb_window.expressed_genes.all.geneScores.tsv"))
coloc.summary.perGene = coloc.summary.perGene %>%
  left_join(mergedEvidence.df %>%
              select(ensembl_id = geneID, geneSymbol = symbol, totalScore, colocScore, networkScore, geneDistScore, codingScore, exprScore) %>%
              mutate(totalScoreNoColoc = totalScore - colocScore) %>% select(-colocScore), by = "ensembl_id") %>%
  arrange(desc(totalScoreNoColoc)) %>%
  filter(!is.na(colocScore)) %>%
  mutate(totalScoreNoColoc_rank = row_number())

write_tsv(coloc.summary.perGene, path=file.path(ad_dir, "coloc", "colocSummary.perGene.tsv"), na="")


colocs.perLocus.df = coloc.eqtl_sqtl.df %>%
  group_by(locus_name) %>%
  summarise(maxInd = which.max(H4),
            ensembl_id = ensembl_id[maxInd],
            geneSymbol = geneSymbol[maxInd],
            signal = signal[maxInd],
            maxH4 = max(H4),
            dataset = dataset[maxInd],
            dataset_short = dataset_short[maxInd]) %>%
  select(-maxInd)

colocs.perLocus.gtex.df = colocs.df %>%
  filter(grepl("gtex", dataset)) %>%
  group_by(locus_name) %>%
  summarise(maxInd = which.max(H4),
            gtex_ensembl_id = ensembl_id[maxInd],
            gtex_geneSymbol = geneSymbol[maxInd],
            gtex_signal = signal[maxInd],
            gtex_maxH4 = max(H4),
            gtex_dataset = dataset[maxInd],
            gtex_dataset_short = dataset_short[maxInd]) %>%
  select(-maxInd)

library(annotables)
grch38.nodup = grch38 %>% filter(!duplicated(symbol))

maxColocScore.perLocus = scores.perGene.df %>%
  left_join(grch38.nodup %>% select(ensgene, geneSymbol = symbol), by=c("ensembl_id"="ensgene")) %>%
  group_by(locus_name) %>%
  summarise(maxInd = which.max(colocScore),
            colocScore = max(colocScore),
            colocScore_ensembl_id = ensembl_id[maxInd],
            colocScore_geneSymbol = geneSymbol[maxInd]) %>%
  select(-maxInd)


coloc.summary.perLocus = colocs.perLocus.df %>%
  left_join(colocs.perLocus.gtex.df, by=c("locus_name")) %>%
  left_join(maxColocScore.perLocus, by=c("locus_name")) %>%
  select(locus_name, colocScore_ensembl_id, colocScore_geneSymbol, colocScore, everything())

write_tsv(coloc.summary.perLocus, path=file.path(ad_dir, "coloc", "colocSummary.perLocus.tsv"), na="")


# Is colocScore any better than maxH4 or gtex_maxH4 at selecting genes
# prioritized by all non-coloc lines of evidence?
coloc.summary.perGene = coloc.summary.perGene %>% filter(!is.na(totalScoreNoColoc))

ggplot(coloc.summary.perGene, aes(x=colocScore, y=maxH4)) + geom_point() + theme_bw() + geom_smooth()

cor(coloc.summary.perGene$colocScore, coloc.summary.perGene$totalScoreNoColoc, method="spearman", use="pairwise.complete")
cor(coloc.summary.perGene$colocScore, coloc.summary.perGene$networkScore, method="spearman", use="pairwise.complete")

ggplot(coloc.summary.perGene, aes(x=colocScore, y=totalScoreNoColoc)) + geom_point() + theme_bw() + geom_smooth()
ggplot(coloc.summary.perGene, aes(x=colocScore_rank, y=totalScoreNoColoc_rank)) + geom_point() + theme_bw() + geom_smooth()
summary(lm(totalScoreNoColoc~colocScore, data=coloc.summary.perGene))
summary(lm(totalScoreNoColoc_rank~colocScore_rank, data=coloc.summary.perGene))

ggplot(coloc.summary.perGene, aes(x=maxH4, y=totalScoreNoColoc)) + geom_point() + theme_bw() + geom_smooth()
ggplot(coloc.summary.perGene, aes(x=maxH4_rank, y=totalScoreNoColoc_rank)) + geom_point() + theme_bw() + geom_smooth()
summary(lm(totalScoreNoColoc~maxH4, data=coloc.summary.perGene))
summary(lm(totalScoreNoColoc_rank~maxH4_rank, data=coloc.summary.perGene))

ggplot(coloc.summary.perGene, aes(x=gtex_maxH4, y=totalScoreNoColoc)) + geom_point() + theme_bw() + geom_smooth()
ggplot(coloc.summary.perGene, aes(x=gtex_maxH4_rank, y=totalScoreNoColoc_rank)) + geom_point() + theme_bw() + geom_smooth()
summary(lm(totalScoreNoColoc~gtex_maxH4, data=coloc.summary.perGene))
summary(lm(totalScoreNoColoc_rank~gtex_maxH4_rank, data=coloc.summary.perGene))

ggplot(coloc.summary.perGene, aes(x=colocScore, y=networkScore)) + geom_point() + theme_bw() + geom_smooth()
summary(lm(networkScore~colocScore, data=coloc.summary.perGene))
summary(lm(networkScore~colocScore_rank, data=coloc.summary.perGene))

ggplot(coloc.summary.perGene, aes(x=maxH4, y=networkScore)) + geom_point() + theme_bw() + geom_smooth()
summary(lm(networkScore~maxH4, data=coloc.summary.perGene))
summary(lm(networkScore~maxH4_rank, data=coloc.summary.perGene))

ggplot(coloc.summary.perGene, aes(x=gtex_maxH4, y=networkScore)) + geom_point() + theme_bw() + geom_smooth()
summary(lm(networkScore~gtex_maxH4, data=coloc.summary.perGene))
summary(lm(networkScore~gtex_maxH4_rank, data=coloc.summary.perGene))

scoreDensity.df = coloc.summary.perGene %>%
  select(totalScoreNoColoc_rank, colocScore, maxH4, gtex_maxH4) %>%
  gather("scoreType", "score", colocScore:gtex_maxH4)
ggplot(scoreDensity.df, aes(x=score, fill=scoreType)) + geom_density(alpha=0.5) + theme_bw()


summary(lm(totalScoreNoColoc_rank <= 40 ~ gtex_maxH4, data=coloc.summary.perGene))
summary(lm(totalScoreNoColoc_rank <= 40 ~ maxH4, data=coloc.summary.perGene))
summary(lm(totalScoreNoColoc_rank <= 40 ~ colocScore, data=coloc.summary.perGene))
ggplot(scoreDensity.df, aes(x=(totalScoreNoColoc_rank <= 40), y=score, fill=scoreType)) + geom_boxplot(outlier.shape=NA) + theme_bw()

ggplot(scoreDensity.df %>% filter(scoreType == "colocScore"), aes(x=(totalScoreNoColoc_rank <= 40), y=score, fill=scoreType)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(aes(group=(totalScoreNoColoc_rank <= 40)), alpha=0.5, width=0.5) + theme_bw() + ggtitle("colocScore")

ggplot(scoreDensity.df %>% filter(scoreType == "maxH4"), aes(x=(totalScoreNoColoc_rank <= 40), y=score, fill=scoreType)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(aes(group=(totalScoreNoColoc_rank <= 40)), alpha=0.5, width=0.5) + theme_bw() + ggtitle("maxH4")

ggplot(scoreDensity.df %>% filter(scoreType == "gtex_maxH4"), aes(x=(totalScoreNoColoc_rank <= 40), y=score, fill=scoreType)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(aes(group=(totalScoreNoColoc_rank <= 40)), alpha=0.5, width=0.5) + theme_bw() + ggtitle("gtex_maxH4")

