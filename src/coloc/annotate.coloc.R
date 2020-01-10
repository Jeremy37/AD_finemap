#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)
library(annotables)
options(stringsAsFactors = F)

args <- commandArgs(trailingOnly = TRUE)
root = args[1]
prefix = args[2]
annotated_gwas = args[3]
outputroot = args[4]

#root = "AD_finemap"
#annotated_gwas = file.path(root, "annotated", "AD.meta.annotated.selected.probable.tsv")
#outputroot = file.path(root, "annotated", "AD.meta.annotated")
#output_file = paste0(outputroot, ".colocs.txt")

#prefix = "coloc.AD.meta"
# outputroot = file.path(root, "annotated", "AD.meta.annotated")
# output_file = paste0(outputroot, ".colocs.txt")

# prefix = "coloc.AD.meta.cond_1"
# outputroot = file.path(root, "annotated", "AD.meta.cond_1.annotated")
# output_file = paste0(outputroot, ".colocs.txt")

# prefix = "coloc.AD.meta.cond_2"
# outputroot = file.path(root, "annotated", "AD.meta.cond_2.annotated")
# output_file = paste0(outputroot, ".colocs.txt")


getAllColocs = function(gwas_prefix) {
  ###############################################################################
  # eQTL catalogue
  catalogue_datasets = readr::read_tsv(file.path(root, "coloc/qtl_data/eqtl_catalogue_datasets.tsv"))
  
  eqtl_catalogue_qtl_list = list()
  for (i in 1:nrow(catalogue_datasets)) {
    dataset = catalogue_datasets$dataset[i]
    dataset_short = catalogue_datasets$dataset_short[i]
    qtl.df = readr::read_tsv(file.path(root, "coloc/output", sprintf("%s.%s.5e+05.txt", gwas_prefix, dataset))) %>%
      dplyr::mutate(dataset = dataset, dataset_short = dataset_short, ensembl_id = feature, geneSymbol = NA)
    eqtl_catalogue_qtl_list = c(eqtl_catalogue_qtl_list, list(qtl.df))
  }
  coloc.eqtl_catalogue.all.df = bind_rows(eqtl_catalogue_qtl_list)
  
  # Some of the eQTL catalogue studies were based on microarray, and the feature
  # is a probe ID rather than gene ID. Here we update these to the gene ID, and
  # just take the top scoring probe per gene per dataset
  probe_map = readr::read_tsv(file.path(root, "coloc/qtl_data/eqtl_catalogue/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv"))
  coloc.eqtl_catalogue.all.df = coloc.eqtl_catalogue.all.df %>%
    left_join(probe_map %>% select(phenotype_id, array_gene_id = gene_id), by=c("feature" = "phenotype_id"))
  # For any feature with "ILMN" in the name, replace the feature with the gene ID
  coloc.eqtl_catalogue.all.df$feature[grepl("ILMN", coloc.eqtl_catalogue.all.df$feature)] = coloc.eqtl_catalogue.all.df$array_gene_id[grepl("ILMN", coloc.eqtl_catalogue.all.df$feature)]
  coloc.eqtl_catalogue.all.df$ensembl_id = coloc.eqtl_catalogue.all.df$feature
  
  coloc.eqtl_catalogue.all.df = coloc.eqtl_catalogue.all.df %>%
    arrange(desc(PP.H4)) %>%
    filter(!duplicated(paste(feature, dataset_short)))
  
  ###############################################################################
  # Blueprint
  coloc.mono.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.mono_gene_eQTL.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint monocyte eQTL", dataset_short = "mono_eQTL",
                  ensembl_id = feature, geneSymbol = NA)
  
  coloc.neut.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.neut_gene_eQTL.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint neutrophil eQTL", dataset_short = "neut_eQTL",
                  ensembl_id = feature, geneSymbol = NA)
  
  coloc.tcel.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.tcel_gene_eQTL.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint t-cell eQTL", dataset_short = "tcel_eQTL",
                  ensembl_id = feature, geneSymbol = NA)
  
  coloc.mono.h3k27ac.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.mono_K27AC.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint monocyte H3K27ac", dataset_short = "mono_h3k27ac",
                  ensembl_id = NA, feature = paste0("peak_", feature), geneSymbol = feature)
  
  coloc.mono.h3k4me1.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.mono_K4ME1.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint monocyte H3K4me1", dataset_short = "mono_h3k4me1",
                  ensembl_id = NA, feature = paste0("peak_", feature), geneSymbol = feature)
  
  coloc.neut.h3k27ac.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.neut_K27AC.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint neutrophil H3K27ac", dataset_short = "neut_h3k27ac",
                  ensembl_id = NA, feature = paste0("peak_", feature), geneSymbol = feature)
  
  coloc.neut.h3k4me1.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.neut_K4ME1.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint neutrophil H3K4me1", dataset_short = "neut_h3k4me1",
                  ensembl_id = NA, feature = paste0("peak_", feature), geneSymbol = feature)
  
  coloc.tcel.h3k27ac.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.tcel_K27AC.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint t-cell H3K27ac", dataset_short = "tcel_h3k27ac",
                  ensembl_id = NA, feature = paste0("peak_", feature), geneSymbol = feature)
  
  coloc.tcel.h3k4me1.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".blueprint.tcel_K4ME1.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Blueprint t-cell H3K4me1", dataset_short = "tcel_h3k4me1",
                  ensembl_id = NA, feature = paste0("peak_", feature), geneSymbol = feature)
  
  
  ###############################################################################
  # xQTL
  # xQTL colocs have the gene symbol in the EnsemblID column
  coloc.xqtl.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".xQTL_eQTL.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Brain ROSMAP eQTL", dataset_short = "xQTL_eQTL",
                  ensembl_id = NA, geneSymbol = feature)
  
  # When getting haQTL and mQTLs, we want to provide some indication as to which
  # peak / CpG a QTL refers to. We need to extract this from another file, in this
  # case the *.minp file for each peak.
  coloc.xqtl.haqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".xQTL_haQTL.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Brain ROSMAP H3K27ac", dataset_short = "xQTL_h3k27ac", ensembl_id = NA, geneSymbol = NA)
  xqtl.haqtl.signals.df =  readr::read_tsv(file.path(root, "coloc", "qtl_data", "xQTL", "xQTL_haQTL.qtl_signals.txt.gz"))
  coloc.xqtl.haqtl.df2 = coloc.xqtl.haqtl.df %>%
    dplyr::left_join(xqtl.haqtl.signals.df %>% dplyr::select(feature, featureChromosome, featurePositionStart), by="feature") %>%
    mutate(feature = paste(feature, featureChromosome, featurePositionStart, sep = "_")) %>%
    select(-featureChromosome, -featurePositionStart)
  rm(xqtl.haqtl.signals.df)
  
  # coloc.xqtl.mqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".xQTL_mQTL.5e+05.txt"))) %>%
  #   dplyr::mutate(dataset = "Brain ROSMAP DNA methylation", dataset_short = "xQTL_meth", ensembl_id = NA)
  # xqtl.mqtl.signals.df =  readr::read_tsv(file.path(root, "coloc", "qtl_data", "xQTL", "xQTL_mQTL.qtl_signals.txt.gz"))
  # xqtl.mqtl.signals.df = xqtl.mqtl.signals.df %>% dplyr::mutate(geneSymbol = paste(feature, featureChromosome, featurePositionStart, sep = "_"))
  # coloc.xqtl.mqtl.df = coloc.xqtl.mqtl.df %>%
  #   dplyr::left_join(xqtl.mqtl.signals.df %>% dplyr::select(feature, geneSymbol), by="feature")
  # rm(xqtl.mqtl.signals.df)
  
  
  ###############################################################################
  # Microglia
  coloc.microglia.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".microglia.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Microglia", dataset_short = "microglia",
                  ensembl_id = feature, geneSymbol = NA)
  
  
  ###############################################################################
  # Brain cortex QTL meta-analysis (1433 samples)
  coloc.brain_meta.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", paste0(gwas_prefix, ".Cortex_Meta.cis_eQTL.5e+05.txt"))) %>%
    dplyr::mutate(dataset = "Brain eQTL meta-analysis", dataset_short = "brain_meta",
                  ensembl_id = feature, geneSymbol = NA)
  
  
  ###############################################################################
  # GTEx
  gtex.tissues = readr::read_tsv(file.path(root, "coloc", "qtl_data", "GTEx_v8", "gtex_tissues.table.tsv"))
  coloc.gtex.eqtl.df = data.frame()
  for (i in 1:nrow(gtex.tissues)) {
    tissue = gtex.tissues[i,]$tissue
    tissue_short = gtex.tissues[i,]$tissue_short
    dataset = paste0("gtex.", tissue_short)
    coloc.tissue.eqtl.df = readr::read_tsv(file.path(root, "coloc", "output", "GTEx_eqtl", paste0(gwas_prefix, ".", tissue, "_eqtl.5e+05.txt"))) %>%
      dplyr::mutate(dataset = paste0("gtex.", tissue), dataset_short = paste0("gtex.", tissue_short),
                    ensembl_id = feature, geneSymbol = NA)
    coloc.gtex.eqtl.df = rbind(coloc.gtex.eqtl.df, coloc.tissue.eqtl.df)
  }
  
  ###############################################################################
  # Make a big table of all the colocs
  
  # Write out a smaller version of the table, including each region-coloc as one row.
  colocs.no_gtex.df = bind_rows(coloc.eqtl_catalogue.all.df,
                                coloc.mono.h3k27ac.df,
                                coloc.mono.h3k4me1.df,
                                coloc.neut.h3k27ac.df,
                                coloc.neut.h3k4me1.df,
                                coloc.tcel.h3k27ac.df,
                                coloc.tcel.h3k4me1.df,
                                coloc.xqtl.eqtl.df,
                                coloc.xqtl.haqtl.df,
                                coloc.microglia.eqtl.df,
                                coloc.brain_meta.eqtl.df)

  colocs.df = bind_rows(colocs.no_gtex.df, coloc.gtex.eqtl.df)
  
  #colocs.no_gtex.df = colocs.no_gtex.df %>% dplyr::rename(gwas_snp = gwas_lead)
  colocs.df = colocs.df %>% dplyr::rename(gwas_snp = gwas_lead)
  
  return(colocs.df)
}


coloc_main.df = getAllColocs("coloc.AD.meta") %>% mutate(signal = "all")
coloc_cond1.df = getAllColocs("coloc.AD.meta.cond_1") %>% mutate(signal = "cond_1")
coloc_cond2.df = getAllColocs("coloc.AD.meta.cond_2") %>% mutate(signal = "cond_2")

colocs.df = bind_rows(coloc_main.df, coloc_cond1.df, coloc_cond2.df)

colocs.df$gwas_signal_rsid = gsub("_.*", "", colocs.df$gwas_signal_rsid)


###############################################################################
# Add in gene symbol where it isn't already present. Since there are sometimes
# more than one ensgene ID for gene symbol, we just take one of these to avoid
# duplication.
grch38_noDupEnsg = grch38 %>% filter(!duplicated(ensgene))
colocs.df = colocs.df %>%
  left_join(grch38_noDupEnsg %>% select(ensgene = ensgene, grch38_symbol = symbol), by=c("ensembl_id" = "ensgene")) %>%
  group_by(feature) %>%
  mutate(geneSymbol = ifelse(is.na(geneSymbol), grch38_symbol, geneSymbol)) %>%
  select(-grch38_symbol)


###############################################################################
# Add in ensgene ID where it isn't already present (xQTL dataset mainly). Since
# there is sometimes more than one gene symbol for an ensgene ID for gene symbol,
# we just take one of these to avoid duplication.
grch38_noDupSymbols = grch38 %>% filter(!duplicated(symbol))
colocs.df = colocs.df %>%
  left_join(grch38_noDupSymbols %>% select(grch38_ensgene = ensgene, grch38_symbol = symbol), by=c("geneSymbol" = "grch38_symbol")) %>%
  group_by(feature) %>%
  mutate(ensembl_id = ifelse(is.na(ensembl_id) & !is.na(grch38_ensgene), grch38_ensgene, ensembl_id)) %>%
  select(-grch38_ensgene)


###############################################################################
# Identify the GWAS signal that each coloc correspond with
# First read in the GWAS fine mapping table
annotated.df = readr::read_tsv(annotated_gwas, guess_max = 10000) %>%
  dplyr::rename(enhPromScores = geneScores, enhPromScoreDetails = geneScoreDetails)

annotated.df$gwas_lead = NA
for (locus in unique(annotated.df$lead_chrpos)) {
  sig.df = annotated.df %>% filter(lead_chrpos == locus) %>%
    dplyr::arrange(META_P)
  annotated.df[annotated.df$lead_chrpos == locus,]$gwas_lead = sig.df[1,]$snp
}
signals.df = annotated.df %>% 
  dplyr::mutate(chrpos = paste0(chr, ".", pos_hg37)) %>%
  dplyr::select(lead_chrpos, locus_name, snp, chrpos, gwas_lead)

# First replace the coloc GWAS SNPs with the rsID from the annotated table,
# matching on chr:pos. This is necessary because the coloc code changes the
# GWAS lead SNP IDs to those from the QTL variant info, and if the QTL variant
# info used a different ID (e.g. chr17_40051621 rather than rs...) then the
# IDs won't match.
colocs.df$chrpos = paste0(colocs.df$chr, ".", colocs.df$gwas_pos)
colocs.df %<>% dplyr::left_join(signals.df, by="chrpos")
colocs.df %<>% dplyr::select(-chrpos, -snp)

# Remove any duplicate entries with the same gwas_lead and qtl_lead SNP,
# and arrange in order of decreasing H4 probability
colocs.df %<>% dplyr::arrange(qtl_pval) %>%
  dplyr::mutate(sortkey = paste(dataset_short, signal, gwas_lead, qtl_lead, feature, sep=".")) %>%
  dplyr::filter(!duplicated(sortkey), !is.na(PP.H4)) %>%
  dplyr::rename(H4 = PP.H4) %>%
  dplyr::arrange(-H4) %>%
  dplyr::select(-gwas_lead) %>%
  dplyr::select(feature:qtl_lead, gwas_lead=gwas_signal_rsid, everything()) %>%
  group_by(sortkey) %>% # Make sure "feature" has ensgene ID whenever available
  mutate(feature = if_else(grepl("ENSG", ensembl_id) & !grepl("ENSG", feature), ensembl_id, feature)) %>%
  ungroup() %>% dplyr::select(-sortkey)
  # dplyr::filter(gwas_pval < 1e-7) %>%
  # dplyr::filter(qtl_pval < 1e-3)

write.table(colocs.df, file = file.path(root, "coloc", "output", paste0(prefix, ".all_colocs.txt")),
            sep="\t", row.names=F, col.names=T, quote=F, na="")

colocs.eqtl.df = colocs.df %>% dplyr::filter(!grepl("k27ac|meth|k4me1", dataset_short, ignore.case = T))
write.table(colocs.eqtl.df, file = file.path(root, "coloc", "output", paste0(prefix, ".eqtl_sqtl_colocs.txt")),
            sep="\t", row.names=F, col.names=T, quote=F, na="")

write.table(colocs.df %>%
              filter(!is.na(locus_name)) %>% # Exclude APOE since colocs aren't reliable
              select(locus_name, signal, feature, ensembl_id, geneSymbol, dataset, nsnps, contains("PP."), H4, qtl_pval, gwas_pval, qtl_lead, gwas_lead, chr, qtl_pos_hg19=qtl_pos, gwas_pos_hg19=gwas_pos, note) %>%
              group_by(locus_name) %>% mutate(minpos = min(qtl_pos_hg19)) %>%
              arrange(chr, minpos, locus_name, desc(H4)) %>%
              ungroup() %>% select(-minpos),
            file = file.path(root, "coloc", "output", paste0(prefix, ".all_colocs.supp_table.txt")),
            sep="\t", row.names=F, col.names=T, quote=F, na="")


colocs.top.df  = colocs.df %>%
  dplyr::filter(H4 > 0.5)
# AD.colocs.top.df %<>% dplyr::filter(dataset_short != "xQTL_meth")

# Summarize colocs for a locus across all QTL datasets
getColocsPerLocus = function(coloc.df, includeDataset = F) {
  coloc.df %<>% dplyr::arrange(-H4)
  
  if (includeDataset) {
    coloc.df$colocStr = sprintf("%s,%s,%.2g", coloc.df$geneSymbol, coloc.df$dataset_short, coloc.df$H4)
  } else {
    coloc.df$colocStr = sprintf("%s,%.2g", coloc.df$geneSymbol, coloc.df$H4)
  }
  topColocs = coloc.df %>% group_by(signal, gwas_lead) %>% summarise() %>% mutate(maxH4 = NA, colocStr = "")
  #topColocs = data.frame(gwas_lead = unique(coloc.df$gwas_lead), maxH4=NA, colocStr="")
  for (i in 1:nrow(topColocs)) {
    coloc.locus.df = coloc.df %>% filter(gwas_lead == topColocs$gwas_lead[i])
    topColocs[i,]$maxH4 = NA
    if (nrow(coloc.locus.df) > 0) {
      topColocs[i,]$maxH4 = max(coloc.locus.df$H4, na.rm=T)
    }
    topColocs[i,]$colocStr = paste(coloc.locus.df$colocStr, collapse=" / ")
  }
  topColocs %>% ungroup()
}

topColocs.df = getColocsPerLocus(colocs.top.df, includeDataset = T)

annotated.df %<>% dplyr::left_join(topColocs.df, by=c("snp" = "gwas_lead")) %>%
  dplyr::rename(topColocs = colocStr)

# Write out the full table, with the top colocs annotated only in the rows of GWAS
# lead SNPs
write.table(annotated.df %>% dplyr::select(-gwas_lead) %>% dplyr::arrange(chr, lead_chrpos, desc(finemap_prob_nc)),
            file=paste0(outputroot, ".colocs.txt"),
            sep="\t", row.names=F, col.names=T, quote=F, na="")



# Write out a more detailed table including all colocs
coloc.summary.df = NULL
for (curDataset in unique(colocs.df$dataset)) {
  topColocs.df = getColocsPerLocus(colocs.df %>% dplyr::filter(dataset == curDataset))
  topColocs.df$dataset = curDataset
  coloc.summary.df = bind_rows(coloc.summary.df, topColocs.df)
}

dataset_levels = c(unique(coloc.summary.df$dataset)[!grepl("gtex", unique(coloc.summary.df$dataset))],
                   unique(coloc.summary.df$dataset)[grepl("gtex", unique(coloc.summary.df$dataset))])
coloc.summary.df$dataset = factor(as.character(coloc.summary.df$dataset), levels=dataset_levels)

# Add in the APOE locus
snps.df = rbind(annotated.df %>% dplyr::select(lead_chrpos, locus_name, snp, chr, pos_hg37),
                data.frame(lead_chrpos = "19_45365447", locus_name = "APOE", snp = "rs146275714", chr = 19, pos_hg37 = "45365447"))

coloc.summary.df = coloc.summary.df %>%
  na.omit() %>%
  dplyr::left_join(snps.df, by=c("gwas_lead" = "snp")) %>%
  dplyr::arrange(chr, pos_hg37, signal, desc(maxH4), dataset) %>%
  dplyr::rename(colocDetails = colocStr) %>%
  dplyr::select(lead_chrpos, locus_name, signal, gwas_lead, dataset, maxH4, colocDetails)

write.table(coloc.summary.df, file=paste0(outputroot, ".coloc_details.txt"),
            sep="\t", row.names=F, col.names=T, quote=F, na="")

