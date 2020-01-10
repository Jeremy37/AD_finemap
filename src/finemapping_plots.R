library(tidyverse)
library(cowplot)
# library(GenomicRanges)
library("GenomicFeatures")
library(rtracklayer)
library(wiggleplotr)
library(biomaRt)
library(annotables)
library("EnsDb.Hsapiens.v86")
options(stringsAsFactors = F)

setwd("/Users/jeremys/work/opentargets/AD_finemap/")

# Just get names of genes near TSPAN14

region_start_hg38 = 80125000
region_end_hg38 = 80740000

hg19_diff = 82280137 - 80520381
region_start_hg19 = region_start_hg38 + hg19_diff
region_end_hg19 = region_end_hg38 + hg19_diff

# Get transcripts and exons for plotting
#listEnsemblArchives()
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "sep2019.archive.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                        "external_gene_name", "strand", 
                        "gene_biotype", "transcript_biotype")
data = getBM(attributes = selected_attributes, mart = ensembl_dataset)
head(data)
transcript_metadata = dplyr::rename(data, 
                                    transcript_id = ensembl_transcript_id, 
                                    gene_id = ensembl_gene_id, 
                                    gene_name = external_gene_name)


txdb_file_hg38 = "/Users/jeremys/work/opentargets/reference/hsapiens_gene_ensembl.txdb.rds"
txdb_file_hg19 = "/Users/jeremys/work/opentargets/reference/hsapiens_gene_ensembl.grch37.txdb.rds"
if (file.exists(txdb_file_hg19)) {
  txdb_hg38 = loadDb(txdb_file_hg38)
  txdb_hg19 = loadDb(txdb_file_hg19)
} else {
  txdb_hg19 = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                                  dataset = "hsapiens_gene_ensembl",
                                  host="grch37.ensembl.org")
  saveDb(txdb_hg19, txdb_file_hg19)
  
  txdb_hg38 = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                                  dataset = "hsapiens_gene_ensembl",
                                  host="sep2019.archive.ensembl.org")
  saveDb(txdb_hg38, txdb_file_hg38)
}

exons_hg19 = exonsBy(txdb_hg19, by = "tx", use.names = TRUE)
cdss_hg19 = cdsBy(txdb_hg19, by = "tx", use.names = TRUE)
exons_hg38 = exonsBy(txdb_hg38, by = "tx", use.names = TRUE)
cdss_hg38 = cdsBy(txdb_hg38, by = "tx", use.names = TRUE)


# Plot genes in the extended GWAS region including TSPAN14
tspan_region.df = grch38 %>%
  dplyr::filter(chr == 10,
                biotype == "protein_coding",
                (start - region_start_hg38 > 0 & start - region_end_hg38 < 0) | (end - region_start_hg38 > 0 & end - region_end_hg38 < 0))

tx_granges = transcripts(EnsDb.Hsapiens.v86, filter = GeneNameFilter(tspan_region.df$symbol))
tspan14_tx_df = as.data.frame(mcols(tx_granges)) %>% group_by(gene_id) %>% dplyr::summarise(tx_id = dplyr::first(tx_id))

p.tspan14_labels = plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = tspan_region.df$symbol, 
                             transcript_ids = tspan14_tx_df$tx_id, rescale_introns=F, transcript_label = T) +
  coord_cartesian(xlim=c(region_start_hg38, region_end_hg38))

p.tspan14_no_labels = plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = tspan_region.df$symbol, 
                                            transcript_ids = tspan14_tx_df$tx_id, rescale_introns=F, transcript_label = F) +
  coord_cartesian(xlim=c(region_start_hg38, region_end_hg38))

# Read in GWAS summary stats, and convert coords to hg38
assoc <- read.table("summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.col_subset.tsv.bgz",
                    header = T, colClasses = c("numeric","numeric","character","numeric"))
tspan_pvals = assoc %>%
  dplyr::filter(CHR == 10, BP > region_start_hg19, BP < region_end_hg19) %>%
  mutate(BP = BP - hg19_diff) %>%
  dplyr::rename(rsid = SNP, pos = BP, p_nominal = META_P)
# Get R2 with lead SNP, and fix SNP labels, which start off with _A1_A2 at the end
tspan_R2 = read_tsv("plots/tspan14_region_r2.rs1878036.tsv")
tspan_R2 = tspan_R2 %>%
  mutate(row = 1:nrow(tspan_R2)) %>%
  group_by(row) %>%
  mutate(rsid = str_split(rsid, "_")[[1]][1]) %>%
  ungroup() %>% dplyr::select(-row)

tspan_pvals = tspan_pvals %>%
  left_join(tspan_R2, by="rsid") %>%
  mutate(track_id = "GWAS",
         R2 = abs(R2))

r2levels = c("0.8", "0.6", "0.4", "0.2", "0")
getR2Groups = function(r2) {
  vals = rep("0", length(r2))
  vals[r2 > 0.2] = "0.2"
  vals[r2 > 0.4] = "0.4"
  vals[r2 > 0.6] = "0.6"
  vals[r2 > 0.8] = "0.8"
  factor(vals, levels = r2levels)
}

plotManhattan = function(pval_df) {
  pval_df = pval_df %>%
    mutate(R2group = getR2Groups(R2))
  pval_df_r2_0 = pval_df %>% dplyr::filter(R2 < 0.2)
  pval_df_r2_0.2 = pval_df %>% dplyr::filter(R2 >= 0.2)
  pval_df_r2_0.4 = pval_df %>% dplyr::filter(R2 >= 0.4)
  pval_df_r2_0.6 = pval_df %>% dplyr::filter(R2 >= 0.6)
  pval_df_r2_0.8 = pval_df %>% dplyr::filter(R2 >= 0.8)
  r2_colours = c("0" = "darkblue", "0.2" = "cyan3", "0.4" = "green3", "0.6" = "orange1", "0.8" = "red1")
  r2_df = data.frame(x=c(0, 0, 0, 0, 0), y=c(1, 1, 1, 1, 1), R2group = factor(r2levels, levels = r2levels))
  ggplot(mapping = aes(x = pos, y = -log10(p_nominal), colour = R2group, fill=R2group)) +
    facet_grid(track_id ~ .) +
    geom_point(data = pval_df_r2_0, alpha=0.8, size=1) + 
    geom_point(data = pval_df_r2_0.2, alpha=0.8, size=1) + 
    geom_point(data = pval_df_r2_0.4, alpha=0.8, size=1) + 
    geom_point(data = pval_df_r2_0.6, alpha=0.8, size=1) + 
    geom_point(data = pval_df_r2_0.8, alpha=0.8, size=1) + 
    theme_classic() + 
    ylab(expression(paste("-",log[10], "(p)"))) +
    geom_bar(data=r2_df, mapping=aes(x=x, y=y, fill=R2group, colour=R2group), stat="identity", position="dodge") +
    scale_colour_manual(values = r2_colours, name="R2") + scale_fill_manual(values = r2_colours, name="R2") +
    theme(legend.justification = c(1, 1), legend.position = c(.98, 1.0),
          legend.key.size = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"))
}


p.tspan14_manhattan = plotManhattan(tspan_pvals) +
  scale_x_continuous(limits = c(region_start_hg38, region_end_hg38), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8))

# Make a manhattan plot with gene track below
# Plot is in hg38 coords
cowplot::plot_grid(p.tspan14_manhattan, p.tspan14_labels, ncol = 1, align = "v", rel_heights = c(1,1))

pdf(file = "plots/tspan14_gwas_region.pdf", width=8, height=3)
cowplot::plot_grid(p.tspan14_manhattan, p.tspan14_labels, ncol = 1, align = "v", rel_heights = c(1,1))
cowplot::plot_grid(p.tspan14_manhattan, p.tspan14_no_labels, ncol = 1, align = "v", rel_heights = c(1,1))
cowplot::plot_grid(p.tspan14_manhattan + theme(legend.position = "none"), p.tspan14_no_labels, ncol = 1, align = "v", rel_heights = c(1,1))
dev.off()


# Make a coverage track for microglia... but these are in GRCh37 coords
microglia_bwdir = "/Users/jeremys/work/opentargets/datasets/glass_microglia/ATAC/bw"
track_data_mg = data.frame(sample_id = c("Kolf2_iPS_ATAC_1", "Kolf2_iPS_ATAC_2", "Kolf2_iPS_ATAC_3"),
                           bigWig = c(file.path(microglia_bwdir, "SRR5955079.bw"),
                                      file.path(microglia_bwdir, "SRR5955080.bw"),
                                      file.path(microglia_bwdir, "SRR5955081.bw"),
                                      file.path(microglia_bwdir, "SRR5955084.bw"),
                                      file.path(microglia_bwdir, "SRR5955085.bw"),
                                      file.path(microglia_bwdir, "SRR5955087.bw"),
                                      file.path(microglia_bwdir, "SRR5955090.bw"),
                                      file.path(microglia_bwdir, "SRR5955091.bw"),
                                      file.path(microglia_bwdir, "SRR5955092.bw")),
                           track_id = c("mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg"),
                           colour_group = c("mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg"),
                           scaling_factor = c(1, 58/81, 93/81, 138/81, 119/81, 133/81, 130/81, 47/81, 54/81))
# Include two microglia ATAC tracks which show a peak near rs1870138
track_data_mg_tspan14 = data.frame(sample_id = c("SRR5955090", "SRR5955091"),
                           bigWig = as.character(c(file.path(microglia_bwdir, "SRR5955090.bw"),
                                                   file.path(microglia_bwdir, "SRR5955091.bw"))),
                           track_id = c("mg", "mg"),
                           colour_group = c("mg", "mg"),
                           scaling_factor = c(130/81, 47/81))

readCoverageFromBigWig <- function(bigwig_path, gene_range){
  #Read coverage over a region from a bigWig file
  sel = rtracklayer::BigWigSelection(gene_range)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  GenomeInfoDb::seqlevels(coverage_ranges) = S4Vectors::as.vector.Rle(GenomicRanges::seqnames(gene_range), mode = "character")
  coverage_rle = GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(GenomicRanges::start(gene_range)):(GenomicRanges::end(gene_range))] #Keep the region of interest
}

track_data = track_data_mg_tspan14
#Read coverage tracks from BigWig file
sample_list = as.list(track_data$bigWig)
names(sample_list) = track_data$sample_id

tpsan14_start_hg38 = 80488000
tpsan14_end_hg38 = 80533123
tpsan14_start_hg19 = tpsan14_start_hg38 + hg19_diff
tpsan14_end_hg19 = tpsan14_end_hg38 + hg19_diff
#tpsan14_start_hg19 = 82269000
#tpsan14_end_hg19 = 82271000

tspan_range = data.frame(chrom="10", start = tpsan14_start_hg19, stop = tpsan14_end_hg19, strand=c("+"))
tspan_gr = makeGRangesFromDataFrame(tspan_range)
coverage_list = lapply(sample_list, readCoverageFromBigWig, tspan_gr)

coverageToDF = function(coverage, start, stop) {
  coverage = S4Vectors::as.vector.Rle(coverage, mode = "double")
  pos = seq(start, stop)
  assertthat::assert_that(assertthat::are_equal(length(pos), length(coverage)))
  new_coverage = dplyr::data_frame(pos = pos, coverage = coverage)
  return(new_coverage)
}
coverage_list = lapply(coverage_list, function(covg) {coverageToDF(covg, start=tpsan14_start_hg19, stop=tpsan14_end_hg19)})

coverage_df = purrr::map_df(coverage_list, identity, .id = "sample_id") %>% 
  as.data.frame() %>%
  dplyr::mutate_(.dots = stats::setNames(list(~as.character(sample_id)), c("sample_id")) ) #Convert factor to character
coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
  dplyr::mutate_(.dots = stats::setNames(list(~coverage/scaling_factor), c("coverage")) ) #Normalize by library size

subsampleCoverage = function(coverage, fraction = 0.1) {
  binsize = 1 / fraction
  minpos = min(coverage$pos, na.rm=T)
  maxpos = max(coverage$pos, na.rm=T)
  coverage %>%
    mutate(bins = as.integer((pos - minpos) / binsize)) %>%
    group_by(sample_id, bigWig, track_id, colour_group, bins) %>%
    summarise(coverage = mean(coverage, na.rm=T), pos = as.integer(mean(pos)), scaling_factor = dplyr::first(scaling_factor)) %>%
    ungroup() %>%
    mutate(bins = pos)
}
coverage_df2 = subsampleCoverage(coverage_df, fraction = 0.005) %>%
  group_by(colour_group, bins) %>%
  summarise(coverage = mean(coverage / scaling_factor, na.rm=T), pos = as.integer(mean(pos))) %>%
  ungroup()

fill_palette = c("blue")
alpha=0.8
tspan14.coverage_plot = ggplot(coverage_df2, aes(bins, coverage, group = colour_group, alpha = alpha)) + 
  geom_blank() +
  theme_classic() +
  geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") +
  geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity") +
  facet_grid(colour_group~.) +
  scale_x_continuous(limits = c(tpsan14_start_hg19, tpsan14_end_hg19), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = fill_palette, guide = F) +
  scale_fill_manual(values = fill_palette, guide = F) +
  ylab("FPM")  
#tspan14.coverage_plot #+ ylim(c(0, 10))

# Get TSPAN14 exons to plot
#tspan14_transcript_ids = transcript_metadata %>% dplyr::filter(gene_name == "TSPAN14") %>% .$transcript_id
tspan14_transcript_ids = "ENST00000429989" # Choose the single longest transcript

tspan14_exons = exons_hg38[intersect(names(exons_hg38), tspan14_transcript_ids)]
tspan14_cdss = cdss_hg38[intersect(names(cdss_hg38), tspan14_transcript_ids)]
#exons_hg19_tspan14 = transcript_metadata %>% dplyr::filter(gene_name == "TSPAN14") %>% dplyr::select(transcript_id) %>%
#  inner_join(as.data.frame(exons_hg19), by=c("transcript_id" = "group_name"))

tspan14.transcript_plot = plotTranscripts(tspan14_exons, tspan14_cdss, transcript_metadata, transcript_label = F, rescale_introns = F) +
  coord_cartesian(xlim = c(tpsan14_start_hg38, tpsan14_end_hg38))

# Read our annotated table to get causal probabilities for SNPs in the region
annot.df = readr::read_tsv("annotated/AD.meta.annotated.selected.probable.paintor.tsv",
                           col_types = cols_only(rsids="c", chr="i", pos_hg37="i", pos_hg38="i", META_P="d", mean_prob="d", finemap_prob_nc="d", paintor_pp="d"))
annot.tspan14.df = annot.df %>% dplyr::filter(chr == 10, pos_hg38 >= tpsan14_start_hg38, pos_hg38 <= tpsan14_end_hg38)

annot.tspan14.lines.df = annot.tspan14.df %>% dplyr::select(pos_hg38, mean_prob)
annot.tspan14.lines.df = bind_rows(annot.tspan14.lines.df,
                                   annot.tspan14.lines.df %>% mutate(mean_prob = 0))
tspan14.prob_plot = ggplot(annot.tspan14.df %>% mutate(track_id = "Fine-mapping"), aes(x=pos_hg38, y=mean_prob*100)) +
  geom_point(col="blue") +
  geom_line(data = annot.tspan14.lines.df, mapping = aes(x=pos_hg38, y=mean_prob*100, group=pos_hg38), col="darkblue") +
  theme_classic() +
  scale_x_continuous(limits = c(tpsan14_start_hg38, tpsan14_end_hg38), expand = c(0,0)) + ylab("mean prob") +
  facet_grid(track_id ~ .)

cowplot::plot_grid(tspan14.prob_plot, tspan14.transcript_plot, tspan14.coverage_plot, ncol = 1, align = "v", rel_heights = c(1,1,1))

pdf(file = "plots/tspan14_finemap.pdf", width=8, height=3.1)
cowplot::plot_grid(tspan14.prob_plot, tspan14.transcript_plot, tspan14.coverage_plot, ncol = 1, align = "v", rel_heights = c(1,1,0.9))
dev.off()


###############################################################################
# APH1B

region_start_hg38 = 63008000
region_end_hg38 = 63608000
hg19_diff = 63569902 - 63277703
region_start_hg19 = region_start_hg38 + hg19_diff
region_end_hg19 = region_end_hg38 + hg19_diff


# Plot genes in the extended GWAS region including APH1B
aph1b_region.df = grch38 %>%
  dplyr::filter(chr == 15,
                biotype == "protein_coding",
                (start - region_start_hg38 > 0 & start - region_end_hg38 < 0) | (end - region_start_hg38 > 0 & end - region_end_hg38 < 0),
                !duplicated(symbol))

tx_granges = transcripts(EnsDb.Hsapiens.v86, filter = GeneNameFilter(aph1b_region.df$symbol))
tx_df = as.data.frame(mcols(tx_granges)) %>% group_by(gene_id) %>% dplyr::summarise(tx_id = dplyr::first(tx_id))

p1.labels = plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = aph1b_region.df$symbol,
                                         transcript_ids = tx_df$tx_id, rescale_introns=F, transcript_label = T) +
  coord_cartesian(xlim=c(region_start_hg38, region_end_hg38))
#p1.transcripts

p1.no_labels = plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = aph1b_region.df$symbol,
                                            transcript_ids = tx_df$tx_id, rescale_introns=F, transcript_label = F) +
  coord_cartesian(xlim=c(region_start_hg38, region_end_hg38))


# Read in GWAS summary stats, and convert coords to hg38
aph1b_pvals = assoc %>%
  dplyr::filter(CHR == 15, BP > region_start_hg19, BP < region_end_hg19) %>%
  mutate(BP = BP - hg19_diff) %>%
  dplyr::rename(rsid = SNP, pos = BP, p_nominal = META_P)
# Get R2 with lead SNP, and fix SNP labels, which start off with _A1_A2 at the end
aph1b_R2 = read_tsv("plots/aph1b_region_r2.rs117618017.tsv")
aph1b_R2 = aph1b_R2 %>%
  mutate(row = 1:nrow(aph1b_R2)) %>%
  group_by(row) %>%
  mutate(rsid = str_split(rsid, "_")[[1]][1]) %>%
  ungroup() %>% dplyr::select(-row)

aph1b_pvals = aph1b_pvals %>%
  left_join(aph1b_R2, by="rsid") %>%
  mutate(track_id = "GWAS",
         R2 = abs(R2))

p.aph1b_manhattan = plotManhattan(aph1b_pvals) +
  scale_x_continuous(limits = c(region_start_hg38, region_end_hg38), expand = c(0,0))

# Make a manhattan plot with gene track below
# Plot is in hg38 coords
cowplot::plot_grid(p.aph1b_manhattan, p1.labels, ncol = 1, align = "v", rel_heights = c(1,1))
pdf(file = "plots/aph1b_gwas_region.pdf", width=4.6, height=3)
cowplot::plot_grid(p.aph1b_manhattan, p1.labels, ncol = 1, align = "v", rel_heights = c(1,1))
cowplot::plot_grid(p.aph1b_manhattan + theme(legend.position = "none"), p1.no_labels, ncol = 1, align = "v", rel_heights = c(1,1))
dev.off()

track_data = track_data_mg
#Read coverage tracks from BigWig file
sample_list = as.list(track_data$bigWig)
names(sample_list) = track_data$sample_id

aph1b_start_hg38 = 63276000
aph1b_end_hg38 = 63306000
aph1b_start_hg19 = aph1b_start_hg38 + hg19_diff
aph1b_end_hg19 = aph1b_end_hg38 + hg19_diff

aph1b_range = data.frame(chrom="15", start = aph1b_start_hg19, stop = aph1b_end_hg19, strand=c("+"))
aph1b_gr = makeGRangesFromDataFrame(aph1b_range)
coverage_list = lapply(sample_list, readCoverageFromBigWig, aph1b_gr)

coverage_list = lapply(coverage_list, function(covg) {coverageToDF(covg, start=aph1b_start_hg19, stop=aph1b_end_hg19)})

coverage_df = purrr::map_df(coverage_list, identity, .id = "sample_id") %>% 
  as.data.frame() %>%
  dplyr::mutate_(.dots = stats::setNames(list(~as.character(sample_id)), c("sample_id")) ) #Convert factor to character
coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
  dplyr::mutate_(.dots = stats::setNames(list(~coverage/scaling_factor), c("coverage")) ) #Normalize by library size

coverage_df2 = subsampleCoverage(coverage_df, fraction = 0.01) %>%
  group_by(colour_group, bins) %>%
  summarise(coverage = mean(coverage / scaling_factor, na.rm=T), pos = as.integer(mean(pos))) %>%
  ungroup()

fill_palette = c("blue")
alpha=0.8
aph1b.coverage_plot = ggplot(coverage_df2, aes(bins, coverage, group = colour_group, alpha = alpha)) + 
  geom_blank() +
  theme_classic() +
  geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") +
  geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity") +
  facet_grid(colour_group~.) +
  scale_x_continuous(limits = c(aph1b_start_hg19, aph1b_end_hg19), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = fill_palette, guide = F) +
  scale_fill_manual(values = fill_palette, guide = F) +
  ylab("FPM")  
#aph1b.coverage_plot

# Get exons to plot
#aph1b_transcript_ids = transcript_metadata %>% dplyr::filter(gene_name == "APH1B") %>% .$transcript_id
aph1b_transcript_ids = c("ENST00000380343", "ENST00000560353") # Choose the single main transcript

aph1b_exons = exons_hg38[intersect(names(exons_hg38), aph1b_transcript_ids)]
aph1b_cdss = cdss_hg38[intersect(names(cdss_hg38), aph1b_transcript_ids)]
# exons_hg38_aph1b = transcript_metadata %>% dplyr::filter(gene_name == "APH1B") %>% dplyr::select(transcript_id) %>%
#  inner_join(as.data.frame(exons_hg38), by=c("transcript_id" = "group_name"))

# plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = "APH1B", 
#                              transcript_ids = transcript_metadata %>% dplyr::filter(gene_name == "APH1B") %>% .$transcript_id, rescale_introns=F, transcript_label = T)
aph1b.transcript_plot = plotTranscripts(aph1b_exons, aph1b_cdss, transcript_metadata, transcript_label = F, rescale_introns = F) +
  coord_cartesian(xlim = c(aph1b_start_hg38, aph1b_end_hg38))
# Read our annotated table to get causal probabilities for SNPs in the region
# annot.df = readr::read_tsv("annotated/AD.meta.annotated.selected.probable.paintor.tsv",
#                            col_types = cols_only(rsids="c", chr="i", pos_hg37="i", pos_hg38="i", META_P="d", mean_prob="d", finemap_prob_nc="d", paintor_pp="d"))
annot.aph1b.df = annot.df %>% dplyr::filter(chr == 15, pos_hg38 >= aph1b_start_hg38, pos_hg38 <= aph1b_end_hg38)

annot.aph1b.lines.df = annot.aph1b.df %>% dplyr::select(pos_hg38, mean_prob)
annot.aph1b.lines.df = bind_rows(annot.aph1b.lines.df,
                                   annot.aph1b.lines.df %>% mutate(mean_prob = 0))
aph1b.prob_plot = ggplot(annot.aph1b.df %>% mutate(track_id = "Fine-mapping"), aes(x=pos_hg38, y=mean_prob*100)) +
  geom_point(col="blue") +
  geom_line(data = annot.aph1b.lines.df, mapping = aes(x=pos_hg38, y=mean_prob*100, group=pos_hg38), col="darkblue") +
  theme_classic() +
  scale_x_continuous(limits = c(aph1b_start_hg38, aph1b_end_hg38), expand = c(0,0), breaks = c(63280000, 63290000, 63300000)) + ylab("mean prob") +
  facet_grid(track_id ~ .)

cowplot::plot_grid(aph1b.prob_plot, aph1b.transcript_plot, aph1b.coverage_plot, ncol = 1, align = "v", rel_heights = c(1,1,1))

pdf(file = "plots/aph1b_finemap.pdf", width=4.6, height=2.3)
cowplot::plot_grid(aph1b.prob_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                   aph1b.coverage_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
                   aph1b.transcript_plot + scale_x_continuous(breaks = c(63280000, 63290000, 63300000), expand = c(0,0)) + theme(axis.title.x = element_blank()),
                   ncol = 1, align = "v", rel_heights = c(1,1,1))
dev.off()



###############################################################################
# CASS4

region_start_hg38 = 56224944
region_end_hg38 = 56224944 + 400000
hg19_diff = 54998544 - 56423488
region_start_hg19 = region_start_hg38 + hg19_diff
region_end_hg19 = region_end_hg38 + hg19_diff


# Plot genes in the extended GWAS region including CASS4
region.df = grch38 %>%
  dplyr::filter(chr == 20,
                biotype == "protein_coding",
                (start - region_start_hg38 > 0 & start - region_end_hg38 < 0) | (end - region_start_hg38 > 0 & end - region_end_hg38 < 0),
                !duplicated(symbol))

tx_granges = transcripts(EnsDb.Hsapiens.v86, filter = GeneNameFilter(region.df$symbol))
tx_df = as.data.frame(mcols(tx_granges)) %>% group_by(gene_id) %>% dplyr::summarise(tx_id = dplyr::first(tx_id))

p1.labels = plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = region.df$symbol,
                                         transcript_ids = tx_df$tx_id, rescale_introns=F, transcript_label = T) +
  coord_cartesian(xlim=c(region_start_hg38, region_end_hg38))
#p1.labels

p1.no_labels = plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = region.df$symbol,
                                            transcript_ids = tx_df$tx_id, rescale_introns=F, transcript_label = F) +
  coord_cartesian(xlim=c(region_start_hg38, region_end_hg38))
#p1.no_labels

# Read in GWAS summary stats, and convert coords to hg38
pvals = assoc %>%
  dplyr::filter(CHR == 20, BP > region_start_hg19, BP < region_end_hg19) %>%
  mutate(BP = BP - hg19_diff) %>%
  dplyr::rename(rsid = SNP, pos = BP, p_nominal = META_P)
# Get R2 with lead SNP, and fix SNP labels, which start off with _A1_A2 at the end
cass4_R2 = read_tsv("plots/cass4_region_r2.rs6014724.tsv")
cass4_R2 = cass4_R2 %>%
  mutate(row = 1:nrow(cass4_R2)) %>%
  group_by(row) %>%
  mutate(rsid = str_split(rsid, "_")[[1]][1]) %>%
  ungroup() %>% dplyr::select(-row)

pvals = pvals %>%
  left_join(cass4_R2, by="rsid") %>%
  mutate(track_id = "GWAS",
         R2 = abs(R2))

p.manhattan = plotManhattan(pvals) +
  scale_x_continuous(limits = c(region_start_hg38, region_end_hg38), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))

# Make a manhattan plot with gene track below
# Plot is in hg38 coords
cowplot::plot_grid(p.manhattan, p1.labels, ncol = 1, align = "v", rel_heights = c(1,1))
pdf(file = "plots/cass4_gwas_region.pdf", width=4.6, height=3)
#cowplot::plot_grid(p.manhattan, p1.labels, ncol = 1, align = "v", rel_heights = c(1,1))
cowplot::plot_grid(p.manhattan + theme(legend.position = "none"), p1.no_labels, ncol = 1, align = "v", rel_heights = c(1,1))
dev.off()

track_data = track_data_mg
#Read coverage tracks from BigWig file
sample_list = as.list(track_data$bigWig)
names(sample_list) = track_data$sample_id

cass4_start_hg38 = 56406000
cass4_end_hg38 = 56460000
#56,428,409
cass4_start_hg19 = cass4_start_hg38 + hg19_diff
cass4_end_hg19 = cass4_end_hg38 + hg19_diff

cass4_range = data.frame(chrom="20", start = cass4_start_hg19, stop = cass4_end_hg19, strand=c("+"))
cass4_gr = makeGRangesFromDataFrame(cass4_range)
coverage_list = lapply(sample_list, readCoverageFromBigWig, cass4_gr)

coverage_list = lapply(coverage_list, function(covg) {coverageToDF(covg, start=cass4_start_hg19, stop=cass4_end_hg19)})

coverage_df = purrr::map_df(coverage_list, identity, .id = "sample_id") %>% 
  as.data.frame() %>%
  dplyr::mutate_(.dots = stats::setNames(list(~as.character(sample_id)), c("sample_id")) ) #Convert factor to character
coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
  dplyr::mutate_(.dots = stats::setNames(list(~coverage/scaling_factor), c("coverage")) ) #Normalize by library size

coverage_df2 = subsampleCoverage(coverage_df, fraction = 0.01) %>%
  group_by(colour_group, bins) %>%
  summarise(coverage = mean(coverage / scaling_factor, na.rm=T), pos = as.integer(mean(pos))) %>%
  ungroup()

fill_palette = c("blue")
alpha=0.8
cass4.coverage_plot = ggplot(coverage_df2, aes(bins, coverage, group = colour_group, alpha = alpha)) + 
  geom_blank() +
  theme_classic() +
  geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") +
  geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity") +
  facet_grid(colour_group~.) +
  scale_x_continuous(limits = c(cass4_start_hg19, cass4_end_hg19), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = fill_palette, guide = F) +
  scale_fill_manual(values = fill_palette, guide = F) +
  ylab("FPM")  
#cass4.coverage_plot

# Get exons to plot
#transcript_ids = transcript_metadata %>% dplyr::filter(gene_name == "CASS4") %>% .$transcript_id
# plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = "CASS4",
#                              transcript_ids = transcript_metadata %>% dplyr::filter(gene_name == "CASS4") %>% .$transcript_id, rescale_introns=F, transcript_label = T)

transcript_ids = c("ENST00000360314") # Choose the single main transcript

cass4_exons = exons_hg38[intersect(names(exons_hg38), transcript_ids)]
cass4_cdss = cdss_hg38[intersect(names(cdss_hg38), transcript_ids)]
# exons_hg38_cass4 = transcript_metadata %>% dplyr::filter(gene_name == "CASS4") %>% dplyr::select(transcript_id) %>%
#  inner_join(as.data.frame(exons_hg38), by=c("transcript_id" = "group_name"))

cass4.transcript_plot = plotTranscripts(cass4_exons, cass4_cdss, transcript_metadata, transcript_label = F, rescale_introns = F) +
  scale_x_continuous(breaks = c(56410000, 56430000, 56450000), expand = c(0,0)) +
  coord_cartesian(xlim = c(cass4_start_hg38, cass4_end_hg38))
# Read our annotated table to get causal probabilities for SNPs in the region
# annot.df = readr::read_tsv("annotated/AD.meta.annotated.selected.probable.paintor.tsv",
#                            col_types = cols_only(rsids="c", chr="i", pos_hg37="i", pos_hg38="i", META_P="d", mean_prob="d", finemap_prob_nc="d", paintor_pp="d"))
annot.cass4.df = annot.df %>% dplyr::filter(chr == 20, pos_hg38 >= cass4_start_hg38, pos_hg38 <= cass4_end_hg38)

annot.cass4.lines.df = annot.cass4.df %>% dplyr::select(pos_hg38, mean_prob)
annot.cass4.lines.df = bind_rows(annot.cass4.lines.df,
                                 annot.cass4.lines.df %>% mutate(mean_prob = 0))
cass4.prob_plot = ggplot(annot.cass4.df %>% mutate(track_id = "Fine-mapping"), aes(x=pos_hg38, y=mean_prob*100)) +
  geom_point(col="blue") +
  geom_line(data = annot.cass4.lines.df, mapping = aes(x=pos_hg38, y=mean_prob*100, group=pos_hg38), col="darkblue") +
  theme_classic() +
  scale_x_continuous(limits = c(cass4_start_hg38, cass4_end_hg38), breaks = c(56410000, 56430000, 56450000), expand = c(0,0)) + ylab("mean prob") +
  facet_grid(track_id ~ .)

cowplot::plot_grid(cass4.prob_plot, cass4.coverage_plot, cass4.transcript_plot, ncol = 1, align = "v", rel_heights = c(1,1,1))

pdf(file = "plots/cass4_finemap.pdf", width=4.6, height=2.3)
cowplot::plot_grid(cass4.prob_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                   cass4.coverage_plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
                   cass4.transcript_plot + theme(axis.title.x = element_blank()),
                   ncol = 1, align = "v", rel_heights = c(1,1,1))
dev.off()
