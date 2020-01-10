library(optparse)
library(tidyverse)
options(stringsAsFactors = F)

# Read loci from input file, and write input file of loci for PAINTOR
option_list <- list(
  make_option(c("--locus_file"), type="character", default=NULL,
              help="Input file with list of independent associated loci (produced by earlier pipeline steps)"),
  make_option(c("--annotated"), type="character", default=NULL,
              help="File with all annotations for all SNPs"),
  make_option(c("--indir"), type="character", default=NULL,
              help="Directory for input files"),
  make_option(c("--outdir"), type="character", default=NULL,
              help="Directory for output files")
)

opt <- parse_args(OptionParser(option_list=option_list))
# opt <- list(locus_file = "AD.loci.exceptAPOE_CLNK.tsv",
#             annotated = "annotated/AD.meta.annotated.all.tsv.gz",
#             indir = "paintor/input",
#             outdir = "paintor/input")

loci.df = readr::read_tsv(opt$locus_file)
loci.df$locus = paste0(loci.df$Chr, "_", loci.df$lead_pos, ".IGAP1_GWAX")

annot.df = readr::read_tsv(opt$annotated, 
                           col_types = cols_only(snp="c", Consequence="c", spliceai_max_DS="d", DeepSEA_funsig="d",
                                                 phastCons="d", phyloP="d", GERP_RS="d",
                                                 captureSeqExonDist="i", captureSeqSpliceDist="i", CADD_PHRED="d",
                                                 iPSC="i", microglia="i", ipsMacrophage="i", NPC="i", Neuron="i", DNase="i",
                                                 `Brain DNase`="i", `Blood & Immune DNase`="i", `Fantom Enh`="i",
                                                 `RoadmapEnh`="i", `Brain Enh`="i", `Blood & Immune Enh`="i")) %>%
  rename(Brain_DNase=`Brain DNase`, Blood_Immune_DNase=`Blood & Immune DNase`, Fantom_Enh=`Fantom Enh`,
         Roadmap_Enh=`RoadmapEnh`, Brain_Enh=`Brain Enh`, Blood_Immune_Enh=`Blood & Immune Enh`)

# Decide on some thresholds to use for turning our annotations into binary
# annotations
checkThresholds = function(annot.df) {
  tbl = table(annot.df$Consequence) %>% as.data.frame() %>% arrange(desc(Freq))
  tbl %>% ggplot(aes(x=fct_reorder(Var1, desc(Freq)), y=Freq)) + geom_bar(stat="identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  sum(is.na(annot.df$spliceai_max_DS))
  quantile(annot.df$spliceai_max_DS, na.rm = T, probs = seq(0.0, 1, 0.01))
  
  sum(is.na(annot.df$DeepSEA_funsig))
  quantile(annot.df$DeepSEA_funsig, na.rm = T, probs = seq(0.0, 1, 0.01))
  quantile(-log10(annot.df$DeepSEA_funsig), na.rm = T, probs = seq(0.0, 1, 0.01))
  
  sum(is.na(annot.df$captureSeqSpliceDist))
  quantile(annot.df$captureSeqSpliceDist, na.rm = T, probs = seq(0.0, 0.1, 0.001))
  
  sum(is.na(annot.df$CADD_PHRED))
  quantile(annot.df$CADD_PHRED, na.rm = T, probs = seq(0.0, 1, 0.01))
  
  sum(is.na(annot.df$phastCons))
  quantile(annot.df$phastCons, na.rm = T, probs = seq(0.0, 1, 0.01))
  
  sum(is.na(annot.df$phyloP))
  quantile(annot.df$phyloP, na.rm = T, probs = seq(0.0, 1, 0.01))
  
  sum(is.na(annot.df$GERP_RS))
  quantile(annot.df$GERP_RS, na.rm = T, probs = seq(0.0, 1, 0.01))
  
  sum(annot.df$iPSC)
  sum(annot.df$microglia)
  sum(annot.df$ipsMacrophage)
  sum(annot.df$NPC)
  sum(annot.df$Neuron)
  sum(annot.df$DNase)
  sum(annot.df$Brain_DNase)
  sum(annot.df$Blood_Immune_DNase)
  sum(annot.df$Fantom_Enh)
}

annot.df = annot.df %>% mutate(spliceai_max_DS = if_else(is.na(spliceai_max_DS), 0, spliceai_max_DS),
                               DeepSEA_funsig = if_else(is.na(DeepSEA_funsig), 1, DeepSEA_funsig),
                               CADD_PHRED = if_else(is.na(CADD_PHRED), 0, CADD_PHRED))
annot.df$spliceai_max_DS[is.na(annot.df$spliceai_max_DS)] = 0
annot.df$DeepSEA_funsig[is.na(annot.df$DeepSEA_funsig)] = 0

annot.df$intronic = grepl("intron_variant", annot.df$Consequence)
annot.df$upstream_gene = grepl("upstream_gene_variant", annot.df$Consequence)
annot.df$downstream_gene = grepl("downstream_gene_variant", annot.df$Consequence)
annot.df$regulatory_region = grepl("regulatory_region_variant", annot.df$Consequence)
annot.df$non_coding_transcript_exon = grepl("non_coding_transcript_exon_variant", annot.df$Consequence)
annot.df$three_prime_UTR = grepl("3_prime_UTR_variant", annot.df$Consequence)
annot.df$missense = grepl("missense_variant", annot.df$Consequence)
annot.df$synonymous = grepl("synonymous_variant", annot.df$Consequence)
annot.df$five_prime_UTR = grepl("5_prime_UTR_variant", annot.df$Consequence)
annot.df$TF_binding_site = grepl("TF_binding_site_variant", annot.df$Consequence)
annot.df$stop_gained = grepl("stop_gained", annot.df$Consequence)
annot.df$frameshift = grepl("frameshift_variant", annot.df$Consequence)
annot.df$inframe_insertion = grepl("inframe_insertion", annot.df$Consequence)
annot.df$inframe_deletion = grepl("inframe_deletion", annot.df$Consequence)
annot.df$start_lost = grepl("start_lost", annot.df$Consequence)

annot.df = annot.df %>%
  mutate(coding_nonsyn = if_else(missense | stop_gained | frameshift | inframe_insertion | inframe_deletion | start_lost, T, F),
         UTR = if_else(five_prime_UTR | three_prime_UTR, T, F),
         exonic = if_else(coding_nonsyn | synonymous | UTR, T, F),
         regulatory = if_else(regulatory_region | TF_binding_site, T, F),
         gene_proximal = if_else(upstream_gene | downstream_gene, T, F))

annot.binary.df = annot.df %>%
  select(snp, intronic, upstream_gene, downstream_gene, regulatory_region,
         non_coding_transcript_exon, three_prime_UTR, missense, synonymous,
         five_prime_UTR, TF_binding_site, stop_gained, frameshift,
         inframe_insertion, inframe_deletion, start_lost, coding_nonsyn,
         UTR, exonic, regulatory, gene_proximal,
         phastCons, phyloP, GERP_RS, microglia, ipsMacrophage,
         DNase, Brain_DNase, Blood_Immune_DNase, Roadmap_Enh,
         Brain_Enh, Blood_Immune_Enh, Fantom_Enh) %>%
  mutate(phastCons_gt_0.1 = phastCons > 0.1,
         phastCons_gt_0.5 = phastCons > 0.5,
         phastCons_gt_0.95 = phastCons > 0.95,
         phyloP_gt_0.5 = phyloP > 0.5,
         phyloP_gt_1 = phyloP > 1,
         phyloP_gt_2 = phyloP > 2,
         GERP_gt_1 = GERP_RS > 1,
         GERP_gt_2 = GERP_RS > 2,
         GERP_gt_3 = GERP_RS > 3,
         spliceai_gt_0 = annot.df$spliceai_max_DS > 0,
         spliceai_gt_0.01 = annot.df$spliceai_max_DS > 0.01,
         spliceai_gt_0.1 = annot.df$spliceai_max_DS > 0.1,
         deepSEA_lt_0.1 = annot.df$DeepSEA_funsig < 0.1,
         deepSEA_lt_0.05 = annot.df$DeepSEA_funsig < 0.05,
         deepSEA_lt_0.01 = annot.df$DeepSEA_funsig < 0.01,
         captureSeqSpliceDist_lt_10 = annot.df$captureSeqSpliceDist < 10,
         CADD_PHRED_gt_5 = annot.df$CADD_PHRED > 5,
         CADD_PHRED_gt_10 = annot.df$CADD_PHRED > 10,
         CADD_PHRED_gt_20 = annot.df$CADD_PHRED > 20,
         microglia_or_macrophage_atac = microglia | ipsMacrophage,
         DNase_gteq10 = DNase >= 10,
         DNase = DNase > 0,
         Brain_DNase = Brain_DNase > 0,
         Blood_Immune_DNase = Blood_Immune_DNase > 0,
         Roadmap_Enh_gteq10 = Roadmap_Enh >= 10,
         Roadmap_Enh = Roadmap_Enh > 0,
         Brain_Enh = Brain_Enh > 0,
         Blood_Immune_Enh = Blood_Immune_Enh > 0) %>%
  select(-phastCons, -phyloP, -GERP_RS)


getLocusDF = function(locus.df) {
  getSnpID = function(id) {
    strs = str_split(id, "_")[[1]]
    if (length(strs) > 3) {
      strs[1] = paste(strs[1:3], collapse = "_")
    }
    strs[1]
  }
  locus.df$snp = sapply(locus.df$rsid, getSnpID)
  locus.df
}
# For each locus, get the list of SNPs from the paintor locus input file
# and join to this our selected annotations from our full annotations file.
# To make binary annotations, we threshold on certain arbitrarily chosen values.
# Write out the annotations file.
for (locus in loci.df$locus) {
  locus.df = readr::read_delim(file.path(opt$indir, locus), delim = " ")
  locus.df = getLocusDF(locus.df) %>% mutate_if(is.logical, as.integer)
  locus.annot.df = locus.df %>% select(snp) %>%
    left_join(annot.binary.df, by = "snp") %>%
    mutate_if(is.logical, as.integer) %>%
    select(-snp)
  if (any(is.na(locus.annot.df))) {
    stop("Error: found NA value when preparing annotations file. This shouldn't happen as we need to match every input SNP to annotations.")
  }
  fpath = file.path(opt$outdir, paste0(locus, ".annotations"))
  write.table(locus.annot.df, file=fpath, col.names = T, row.names=F, sep=" ", quote=F, na="")
}

