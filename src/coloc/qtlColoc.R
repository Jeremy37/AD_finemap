#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(coloc)
library(GenomicRanges)
library(Rsamtools)
library(locuscomparer)

gene_id_map = NULL

#Parse command-line options
option_list <- list(
  make_option(c("--qtl_name"), type="character", default=NULL,
              help="Name of QTLs, used for output."),
  make_option(c("--qtl_nominal"), type="character", default=NULL,
              help="Path to the QTL nominal p values file"),
  make_option(c("--qtl_variant_info"), type="character", default=NULL,
              help="Path to a file with QTL SNP information, i.e. rsID, chr, pos, MAF, which will replace similar columns in the QTL summary files. This is useful for switching between GRCh37/38 coords."),
  make_option(c("--qtl_signals"), type="character", default=NULL,
              help="Path to the QTL min p value per gene file"),
  make_option(c("--qtl_samplesize"), type="integer", default=NULL,
              help="Integer sample size of the QTL study"),
  make_option(c("--qtl_list"), type="character", default=NULL,
              help="Path to a file listing QTL dataset to use for coloc, with each line giving qtl_name, qtl_nominal file, qtl_signals file, and qtl_samplesize"),
  make_option(c("--gwas_name"), type="character", default=NULL,
              help="Name of the GWAS trait"),
  make_option(c("--gwas_file"), type="character", default=NULL,
              help="Path to GWAS summary stats directory."),
  make_option(c("--gwas_signals"), type="character", default=NULL,
              help="Threshold of GWAS p value to use locus for colocalisation"),
  make_option(c("--window"), type="numeric", default=NULL,
              help="Size of the cis window around both the QTL and GWAS lead SNP where summary statistics should be used."),
  make_option(c("--overlap_dist"), type="numeric", default="5e5",
              help="Max distance between QTL lead SNP and GWAS lead SNP."),
  make_option(c("--plot_style"), type="character", default=NULL,
              help="Type of plots to generate as a PDF output - either 'overlay' or 'line' or 'all'"),
  make_option(c("--plot_H4_threshold"), type="numeric", default=0,
              help="Minimum H4 result from coloc for a plot to be generated"),
  make_option(c("--plot_p_threshold"), type="numeric", default=1,
              help="In plots, points with p greater than this threshold in both GWAS and QTL datasets will not be shown."),
  make_option(c("--plot_snp_threshold"), type="numeric", default=1000,
              help="In plots, show at most this many top SNPs for both GWAS and QTL datasets."),
  make_option(c("--outdir"), type="character", default=NULL,
              help="Path to the output directory."),
  make_option(c("--match_snps_by_chrpos"), type="logical", default=FALSE,
              help="Whether to match GWAS and QTL SNPs using chr:pos rather than rsid"),
  make_option(c("--p1"), type="numeric", default=1e-4,
              help="Prior probability of SNP affecting expression"),
  make_option(c("--p2"), type="numeric", default=1e-4,
              help="Prior probability of SNP affecting GWAS trait"),
  make_option(c("--p12"), type="numeric", default=1e-5,
              help="Prior probability of SNP affecting both traits"),
  make_option(c("--gene_id_map"), type="character", default=NULL,
              help="Path to file for mapping from ensembl IDs to HGNC symbols. File should have columns ensembl_gene_id,symbol_type,symbol, where symbol_type is either 'symbol' (the main HGNC symbol) or alias1, alias2, prev_symbol1, etc.")
)

importQTLSignals = function(qtlSignalsPath) {
  # This reads a table of significant QTLs, all of which will be tested for coloc
  # with any nearby gwas signal. SNP chr and pos will be mapped using the QTL
  # variant info passed in
  table = readr::read_tsv(qtlSignalsPath, col_names = T, guess_max = 1e5)
  expectedCols = c("feature", "chr", "pos", "rsid")
  if (!hasCols(table, expectedCols)) {
    stop(sprintf("In table %s, expected columns with header names %s.",
                 qtlSignalsPath, paste(expectedCols, collapse = ", ")))
  }
  return(table %>% dplyr::select(feature, chr, pos, rsid))
}

importVariantInfo = function(variantInfoPath) {
  table = readr::read_tsv(variantInfoPath, col_names = T, col_types = "cicd")
  expectedCols = c("chr", "pos", "rsid", "MAF")
  if (!hasCols(table, expectedCols)) {
    stop(sprintf("In table %s, expected columns with header names %s.",
                 variantInfoPath, paste(expectedCols, collapse = ", ")))
  }
  table$chr = as.character(table$chr)
  return(table %>% dplyr::select(chr, pos, rsid, MAF))
}

importGWASSignals = function(gwasSignalsPath) {
  table = readr::read_tsv(gwasSignalsPath)
  expectedCols = c("chr", "pos", "rsid")
  if (!hasCols(table, expectedCols)) {
    stop(sprintf("In table %s, expected columns with header names %s.",
                 gwasSignalsPath, paste(expectedCols, collapse = ", ")))
  }
  table$chr = as.character(table$chr)
  return(table %>% dplyr::select(chr, pos, rsid))
}


tabixFetchQTLs = function(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes) {
  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "feature"))
  qtl_colnames = c("feature", "chr", "pos", "rsid", "p_nominal")
  qtl_coltypes = c("ccicd")
  result = list()
  for (i in seq_along(phenotype_ranges)) {
    selected_feature = phenotype_ranges[i]$feature
    #print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], 
                                     col_names = qtl_colnames, col_types = qtl_coltypes)[[1]] %>%
      dplyr::filter(feature == selected_feature)
    
    #Add additional columns
    result[[selected_feature]] = tabix_table
  }
  return(result[[1]])
}


main = function() {
  opt <- parse_args(OptionParser(option_list=option_list))
  # opt = list(qtl_name = "Cortex_Meta.cis_eQTL", qtl_samplesize = 1433, gwas_name = "AD.meta",
  #            qtl_nominal = file.path(root, "qtl_data/brain_meta/Cortex_Meta.cis_eQTL.nominals.txt.gz"),
  #            qtl_variant_info = file.path(root, "qtl_data/brain_meta/Cortex_Meta.cis_eQTL.variant_info.txt.gz"),
  #            qtl_signals = file.path(root, "qtl_data/brain_meta/Cortex_Meta.cis_eQTL.qtl_signals.FDR_0.05.forColoc.txt"),
  #            gwas_file = file.path(root, "GWAS/Alzheimers_disease_Liu_2019_UKBB_GWAX_meta_v3.sorted.txt.gz"),
  #            gwas_signals = file.path(root, "AD.signals.for_coloc.txt"),
  #            window = 2e5, overlap_dist = 5e5,
  #            plot_style = "all", plot_H4_threshold = 0, plot_p_threshold = 0.1, plot_snp_threshold = 1000,
  #            p1 = 1e-4, p2 = 1e-4, p12 = 1e-5,
  #            match_snps_by_chrpos = F,
  #            outdir = file.path(root, "output"),
  #            gene_id_map = "/Users/jeremys/work/opentargets/reference/hgnc.ensembl.map.txt")
  print(opt)
  
  if (!is.null(opt$plot_style)) {
    if (!(grepl("overlay", opt$plot_style) | grepl("line", opt$plot_style) | grepl("locuscompare", opt$plot_style) | grepl("all", opt$plot_style))) {
      stop(sprintf("Error: unrecognized value for plot_style parameter: %s", opt$plot_style))
    }
  }
  
  gwas_signals = importGWASSignals(opt$gwas_signals) %>%
    dplyr::transmute(chr, gwas_pos = pos, gwas_rsid=rsid) %>%
    dplyr::filter(!is.na(gwas_pos))
  
  qtl_signals = importQTLSignals(opt$qtl_signals) %>%
    dplyr::rename(qtl_pos_orig = pos) %>% dplyr::select(-chr)
  
  qtl_variant_info = importVariantInfo(opt$qtl_variant_info)
  new_qtl_signals = qtl_signals %>%
    dplyr::left_join(qtl_variant_info, by="rsid") %>%
    dplyr::select(-MAF) %>%
    dplyr::rename(qtl_pos = pos, qtl_rsid = rsid)
  if (sum(is.na(new_qtl_signals$qtl_pos)) > 0) {
    warning(sprintf("WARNING: %d of %d QTL signals did not have positions mapped in the QTL variant info", sum(is.na(new_qtl_signals$qtl_pos)), nrow(qtl_signals)))
    new_qtl_signals = new_qtl_signals %>% dplyr::filter(!is.na(qtl_pos))
  }
  qtl_signals = new_qtl_signals
  
  colocPairs = getColocPairs(qtl_signals, gwas_signals, opt$overlap_dist)
  
  # Load table mapping ENSG ID to HGNC symbol, so we can use symbol in plots,
  # which are produced in the colocMolecularQTLs calls
  if (!is.null(opt$gene_id_map)) {
    gene_id_map = readr::read_tsv(opt$gene_id_map) %>%
      dplyr::filter(symbol_type == "symbol")
  }
  
  ###################### Run the coloc for each QTL-GWAS pair! ##################
  write(sprintf("Testing coloc for %d QTL-GWAS signal pairs.", nrow(colocPairs)), stderr())
  #colocPairs = colocPairs[1:2,]
  colocs_list = vector("list", nrow(colocPairs))
  for (i in 1:nrow(colocPairs)) {
      colocs_list[[i]] = colocMolecularQTLs(colocPairs[i,],
                             qtl_summary_path = opt$qtl_nominal,
                             gwas_summary_path = opt$gwas_file,
                             qtl_variant_info = qtl_variant_info,
                             N_qtl = opt$qtl_samplesize,
                             cis_dist = opt$window,
                             p1 = opt$p1, p2 = opt$p2, p12 = opt$p12,
                             match_snps_by_chrpos = opt$match_snps_by_chrpos,
                             plot_style = opt$plot_style, plot_H4_threshold = opt$plot_H4_threshold,
                             plot_p_threshold = opt$plot_p_threshold, plot_snp_threshold = opt$plot_snp_threshold,
                             gene_id_map = gene_id_map)
  }
  
  coloc_hits = bind_cols(colocPairs, bind_rows(lapply(colocs_list, FUN = function(x) x$summary)))
  
  coloc_hits$feature = gsub("\\.[\\d]+", "", coloc_hits$feature, perl=T)
  
  coloc_hits = coloc_hits %>%
    dplyr::rename(PP.H0=PP.H0.abf, PP.H1=PP.H1.abf, PP.H2=PP.H2.abf, PP.H3=PP.H3.abf, PP.H4=PP.H4.abf) %>%
    dplyr::select(feature, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval,
                  qtl_lead, gwas_lead, gwas_signal_rsid = gwas_rsid, chr, qtl_pos, gwas_pos, note)
  
  coloc_output = file.path(opt$outdir, paste("coloc", opt$gwas_name, opt$qtl_name, opt$overlap_dist, "txt", sep = "."))
  write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (!is.null(opt$plot_style)) {
    coloc_plots = lapply(colocs_list, FUN = function(x) x$plots)
    pdf(file.path(opt$outdir, sprintf("coloc.%s.%s.%s.plots.pdf", opt$gwas_name, opt$qtl_name, opt$overlap_dist)),
        width=8, height=5.5)
    for (colocPlot in coloc_plots) {
      if (!is.null(colocPlot)) {
        print(colocPlot)
      }
    }
    dev.off()
  }
}


###############################################################################
# colocTools

hasCols = function(table, cols) {
  all(sapply(cols, function(x) {x %in% colnames(table)}))
}

#' Determine the set of QTL-GWAS signal pairs to test for colocalisation, based
#' on an overlap window.
#'
#' @param qtl_signals List of data frames with QTL lead pvalues. Each data frame must contain
#' gene_id, snp_id and either p_val or FDR columns, and should not contain other columns.
#' @param gwas_signals Prefix of the GWAS summarystats file
#' @param overlap_dist Max distance between GWAS and QTL variants.
#'
#' @return List of data.frames with phenotype_ids and snp_ids to be tested with coloc.
#' @export
getColocPairs <- function(qtl_signals, gwas_signals, overlap_dist = 5e5) {
  coloc_list = qtl_signals %>%
    dplyr::left_join(gwas_signals, by="chr") %>%
    dplyr::mutate(qtl_gwas_distance = abs(gwas_pos - qtl_pos)) %>%
    dplyr::filter(qtl_gwas_distance < overlap_dist)
  coloc_list
}

colocMolecularQTLs <- function(coloc_pair, qtl_summary_path, gwas_summary_path,
                               qtl_variant_info = NULL, N_qtl = 84, cis_dist = 1e5,
                               p1 = 1e-04, p2 = 1e-04, p12 = 1e-05,
                               match_snps_by_chrpos = F,
                               plot_style = NULL, plot_H4_threshold = 0,
                               plot_p_threshold = 1, plot_snp_threshold = 1e5,
                               gene_id_map = NULL) {
  
  expectedCols = c("feature", "chr", "qtl_pos_orig", "gwas_pos", "qtl_rsid", "gwas_rsid")
  assertthat::assert_that(hasCols(coloc_pair, expectedCols))
  assertthat::assert_that(nrow(coloc_pair) == 1)
  
  
  #Print for debugging
  write(sprintf("\nTesting coloc for GWAS locus %s, %s:%d and QTL: %s, %s:%d, %s",
                coloc_pair$gwas_rsid, coloc_pair$chr, coloc_pair$gwas_pos,
                coloc_pair$feature, coloc_pair$chr, coloc_pair$qtl_pos, coloc_pair$qtl_rsid),
        stderr())
  
  coloc_pair = as.list(coloc_pair[1,])
  minPos = min(coloc_pair$qtl_pos, coloc_pair$gwas_pos)
  maxPos = max(coloc_pair$qtl_pos, coloc_pair$gwas_pos)
  diffPos = coloc_pair$gwas_pos - coloc_pair$qtl_pos
  
  qtl_summaries = NULL
  gwas_summaries = NULL
  colocPlots = NULL
  result = tryCatch({
    # Make GRanges objects with a distance of cis_dist upstream and downstream
    # of the GWAS and QTL variants. The QTL range should be relative to the coords
    # of the original QTL signal, before coords were replaced
    qtlStart = coloc_pair$qtl_pos_orig - cis_dist
    qtlEnd = coloc_pair$qtl_pos_orig + cis_dist
    if (diffPos > 0) {
      qtlEnd = qtlEnd + diffPos
    } else {
      qtlStart = qtlStart + diffPos
    }
    qtlRange = GRanges(feature = coloc_pair$feature, seqnames = coloc_pair$chr, ranges = IRanges(start = qtlStart, end = qtlEnd), strand = "*")
    gwasRange = GRanges(seqnames = coloc_pair$chr, ranges = IRanges(start = minPos - cis_dist, end = maxPos + cis_dist), strand = "*")
    
    if (!is.null(gene_id_map) & grepl("^ENSG", coloc_pair$feature)) {
      #coloc_pair$feature = gsub("\\.[\\d]+", "", coloc_pair$feature, perl=T)
      coloc_pair$feature = ifelse(grepl("ENSG", coloc_pair$feature), gsub("\\.[\\d]+", "", coloc_pair$feature, perl=T), coloc_pair$feature)
      coloc_pair$feature = sprintf("%s / %s", gene_id_map[gene_id_map$ensembl_gene_id == coloc_pair$feature,][1,]$symbol, coloc_pair$feature)
    }
    
    write(sprintf("QTL range: %s:%s:%d-%d\tGWAS range: %s:%d-%d",
                  coloc_pair$feature, coloc_pair$chr, qtlStart, qtlEnd,
                  coloc_pair$chr, minPos - cis_dist, maxPos + cis_dist),
          stderr())
    
    # Fetch QTL summary stats. The code calling colocMolecularQTLs needs to define
    # the function tabixFetchQTLs, which will handle fetching QTLs from potentially
    # differently formatted files.
    qtl_summaries = tabixFetchQTLs(qtlRange, qtl_summary_path)
    if (!is.null(qtl_variant_info)) {
      qtl_summaries = summaryReplaceCoordinates(qtl_summaries, qtl_variant_info)
    }
    
    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWASSummary(gwasRange, gwas_summary_path)[[1]]
    #if (match_snps_by_chrpos) {
    #  gwas_summaries = summaryReplaceIDs(gwas_summaries, qtl_variant_info)
    #}
    qtl_min_index = which.min(qtl_summaries$p_nominal)
    gwas_min_index = which.min(gwas_summaries$p_nominal)
    
    # Get DF with SNPs from both GWAS and QTL
    merged.df = getMergedDF(qtl_summaries, gwas_summaries, match_snps_by_chrpos)
    merged.flt.df = merged.df %>% filter(!is.na(qtl_rsid), !is.na(gwas_rsid))
    write(sprintf("%d SNPs in common, of %d QTL SNPs and %d GWAS SNPs", nrow(merged.flt.df), nrow(qtl_summaries), nrow(gwas_summaries)), stderr())
    if (nrow(merged.df) < 2) {
      stop("Too few SNPs in common between QTL and GWAS to run colocalisation.")
    }
    qtl_lead_rsid = qtl_summaries[qtl_min_index, ]$rsid
    gwas_lead_rsid = gwas_summaries[gwas_min_index, ]$rsid
    note = ""
    if (!qtl_lead_rsid %in% merged.flt.df$qtl_rsid) {
      note = sprintf("Warning: QTL lead SNP %s is not in set of SNPs in common", qtl_lead_rsid)
      write(note, stderr())
    }
    if (!gwas_lead_rsid %in% merged.flt.df$gwas_rsid) {
      note2 = sprintf("Warning: GWAS lead SNP %s is not in set of SNPs in common", gwas_lead_rsid)
      write(note2, stderr())
      note = paste(note, note2, sep = "; ")
    }
    
    #Perform coloc analysis
    coloc_res = colocQtlGWAS(merged.flt.df, N_qtl = N_qtl, p1, p2, p12)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_summaries[qtl_min_index,]$p_nominal, gwas_pval = gwas_summaries[gwas_min_index,]$p_nominal,
                    qtl_lead = qtl_summaries[qtl_min_index,]$rsid, gwas_lead = gwas_summaries[gwas_min_index,]$rsid,
                    note = note)
    # Old version - which used minp vals to get lead SNPs rather than using the ones passed in
    # coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
    #   dplyr::mutate(qtl_pval = qtl_summaries[qtl_min_index,]$p_nominal, gwas_pval = gwas_summaries[gwas_min_index,]$p_nominal,
    #                 qtl_lead = qtl_summaries[qtl_min_index,]$rsid, gwas_lead = gwas_summaries[gwas_min_index,]$rsid,
    #                 note = coloc_res$note)
    
    if (!is.null(plot_style) & coloc_summary$PP.H4.abf >= plot_H4_threshold) {
      plot.df = getPlotDF(qtl_summaries, gwas_summaries, coloc_pair, plot_p_threshold, plot_snp_threshold, match_snps_by_chrpos = match_snps_by_chrpos)
      if (grepl("overlay|all", plot_style)) {
        colocPlots = c(colocPlots, list(colocOverlayPlot(plot.df, coloc_pair, coloc_summary)))
      }
      if (grepl("line|all", plot_style)) {
        colocPlots = c(colocPlots, list(colocLinesPlot(plot.df, coloc_pair, coloc_summary)))
      }
      if (grepl("locuscompare|all", plot_style)) {
        plotTitle = sprintf("gwas %s, qtl %s: %s", coloc_pair$gwas_rsid, coloc_pair$qtl_rsid, coloc_pair$feature)
        colocPlots = c(colocPlots, list(locusComparePlot(plot.df, plotTitle)))
      }
    }
    
    result = list(summary = coloc_summary, data = list(qtl = qtl_summaries, gwas = gwas_summaries), plots = colocPlots)
  }, error = function(err) {
    print(paste("ERROR:", err))
    result = list(summary = data.frame(nsnps=NA, PP.H0.abf=NA, PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA, note=NA, qtl_pval=NA, gwas_pval=NA, qtl_lead=NA, gwas_lead=NA),
                  data = list(qtl = qtl_summaries, gwas = gwas_summaries), plots = colocPlots)
  }
  )
  return(result)
}


#' Test colocalisation between molecular QTL and GWAS summary stats
#'
#' @param qtl QTL summary stats (p_nominal, MAF, beta, snp_id)
#' @param gwas GWAS summary stats(beta, se, MAF, log_OR)
#' @param N_qtl Sample size of the QTL mapping study
#'
#' @return coloc.abf result object
#' @export
colocQtlGWAS <- function(merged, N_qtl, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05, case_control_proportion = 0.2) {
  #Check that QTL df has all correct names
  assertthat::assert_that(assertthat::has_name(merged, "qtl_rsid"))
  #assertthat::assert_that(assertthat::has_name(merged, "beta"))
  assertthat::assert_that(assertthat::has_name(merged, "qtl_MAF"))
  assertthat::assert_that(assertthat::has_name(merged, "qtl_p_nominal"))
  
  #Check that GWAS df has all correct names
  assertthat::assert_that(assertthat::has_name(merged, "gwas_beta"))
  assertthat::assert_that(assertthat::has_name(merged, "gwas_se"))
  assertthat::assert_that(assertthat::has_name(merged, "gwas_rsid"))
  assertthat::assert_that(assertthat::has_name(merged, "gwas_MAF"))
  assertthat::assert_that(assertthat::has_name(merged, "gwas_log_OR"))
  assertthat::assert_that(assertthat::has_name(merged, "gwas_p_nominal"))
  
  #Count NAs for log_OR and beta
  log_OR_NA_count = length(which(is.na(merged$gwas_log_OR)))
  beta_NA_count = length(which(is.na(merged$gwas_beta)))
  
  #Remove GWAS SNPs with NA std error
  merged = dplyr::filter(merged, !is.na(gwas_se))
  
  
  # Make QTL dataset object
  df1 = list(pvalues = merged$qtl_p_nominal, 
             N = N_qtl, 
             MAF = merged$qtl_MAF, 
             type = "quant", 
             beta = NA,
             snp = merged$gwas_rsid)

  #If beta is not specified in the GWAS then use log_OR
  if(beta_NA_count <= log_OR_NA_count){
    assertthat::assert_that(assertthat::has_name(merged, "gwas_MAF"))
    coloc_res = coloc::coloc.abf(dataset1 = df1,
                                 dataset2 = list(beta = merged$gwas_beta, 
                                                 varbeta = merged$gwas_se^2, 
                                                 type = "cc",
                                                 s = case_control_proportion,
                                                 snp = merged$gwas_rsid, 
                                                 MAF = merged$gwas_MAF),
                                 p1 = p1, p2 = p2, p12 = p12)
  } else{
    coloc_res = coloc::coloc.abf(dataset1 = df1,
                                 dataset2 = list(beta = merged$gwas_log_OR, 
                                                 varbeta = merged$gwas_se^2, 
                                                 type = "cc", 
                                                 s = case_control_proportion,
                                                 snp = merged$gwas_rsid),
                                 p1 = p1, p2 = p2, p12 = p12)
  }
  return(coloc_res)
}


#' Import a specific region from a tabix-indexed GWAS summary stats file
tabixFetchGWASSummary <- function(granges, summary_path) {
  gwas_col_names = c("rsid", "chr", "pos", "effect_allele", "MAF", 
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  #  gwas_col_names = c("snp_id", "chr", "pos", "MAF", "p_nominal", "beta", "se")
  #  gwas_col_types = c("ccidddd")
  gwas_pvalues = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names, col_types = gwas_col_types)
  return(gwas_pvalues)
}


#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList 
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...) {
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}


summaryReplaceCoordinates <- function(summary_df, variant_information) {
  #Remove MAF if it is present
  if ("MAF" %in% colnames(summary_df)) {
    summary_df = dplyr::select(summary_df, -chr, -pos, -MAF) %>%
      dplyr::left_join(variant_information, by = "rsid") %>%
      dplyr::filter(!is.na(pos))
  } else {
    summary_df = dplyr::select(summary_df, -chr, -pos) %>%
      dplyr::left_join(variant_information, by = "rsid") %>%
      dplyr::filter(!is.na(pos))
  }
  return(summary_df)
}

summaryReplaceIDs <- function(summary_df, variant_information) {
  #Remove MAF if it is present
  summary_df = summary_df %>%
      dplyr::left_join(variant_information %>% dplyr::select(pos, qtl_rsid = rsid), by = "pos")
  summary_df$rsid[!is.na(summary_df$qtl_rsid)] = summary_df$qtl_rsid[!is.na(summary_df$qtl_rsid)]
  return(summary_df)
}

getMergedDF <- function(qtl, gwas, match_snps_by_chrpos) {
  if (match_snps_by_chrpos) {
    merged = qtl %>% select(chr, pos, qtl_rsid = rsid, qtl_p_nominal = p_nominal, qtl_MAF = MAF, everything()) %>%
      dplyr::full_join(gwas %>% dplyr::select(chr, pos, gwas_rsid = rsid, gwas_p_nominal = p_nominal, gwas_MAF = MAF, gwas_beta = beta, gwas_se = se, gwas_log_OR = log_OR, everything()),
                        by=c("chr", "pos"))
  } else {
    merged = qtl %>% select(chr, pos, qtl_rsid, qtl_p_nominal = p_nominal, qtl_MAF = MAF, everything()) %>%
      dplyr::full_join(gwas %>% dplyr::select(chr, pos, rsid, gwas_p_nominal = p_nominal, gwas_MAF = MAF, gwas_beta = beta, gwas_se = se, gwas_log_OR = log_OR, everything()) %>%
                         dplyr::mutate(gwas_rsid = rsid),
                        by=c("qtl_rsid"="rsid"))
  }
  merged
}

getPlotDF <- function(qtl_summaries, gwas_summaries, coloc_pair, plot_p_threshold = 1, plot_snp_threshold = 1e5, match_snps_by_chrpos) {
  qtl.df = qtl_summaries %>% dplyr::select(one_of(c("rsid", "pos", "p_nominal"))) %>% dplyr::arrange(p_nominal)
  gwas.df = gwas_summaries %>% dplyr::select(one_of(c("rsid", "pos", "p_nominal"))) %>% dplyr::arrange(p_nominal)
  plot.df = bind_rows(qtl.df %>% dplyr::mutate(study = "QTL"),
                      gwas.df %>% dplyr::mutate(study = "GWAS"))
  
  # Determine which variants are present in both studies
  qtl.df = qtl.df %>% dplyr::mutate(keep_qtl = 1)
  gwas.df = gwas.df %>% dplyr::mutate(keep_gwas = 1)
  if (plot_snp_threshold < nrow(qtl.df)) {
    qtl.df[(plot_snp_threshold+1):nrow(qtl.df),]$keep_qtl = 0
  }
  if (plot_snp_threshold < nrow(gwas.df)) {
    gwas.df[(plot_snp_threshold+1):nrow(gwas.df),]$keep_gwas = 0
  }
  if (match_snps_by_chrpos) {
    sharing.df = dplyr::full_join(qtl.df %>% dplyr::rename(qtl_pos = pos, qtl_p = p_nominal, qtl_rsid = rsid),
                                  gwas.df %>% dplyr::rename(gwas_p = p_nominal) %>%
                                    dplyr::mutate(gwas_pos = pos),
                                  by=c("qtl_pos"="pos")) %>%
      mutate(rsid = if_else(is.na(rsid), qtl_rsid, rsid),
             pos = if_else(is.na(qtl_pos), gwas_pos, qtl_pos))
  } else {
    sharing.df = dplyr::full_join(qtl.df %>% dplyr::rename(qtl_pos = pos, qtl_p = p_nominal),
                                  gwas.df %>% dplyr::rename(gwas_pos = pos, gwas_p = p_nominal),
                                  by="rsid")
  }
  sharing.df$variant = "Shared"
  sharing.df$variant[is.na(sharing.df$gwas_pos) & !is.na(sharing.df$qtl_pos)] = "Not shared"
  sharing.df$variant[is.na(sharing.df$qtl_pos) & !is.na(sharing.df$gwas_pos)] = "Not shared"
  sharing.df$min_p = sapply(1:nrow(sharing.df), FUN=function(i) min(sharing.df[i,]$qtl_p, sharing.df[i,]$gwas_p, na.rm = T))
  sharing.df$keep_snp = sapply(1:nrow(sharing.df), FUN=function(i) (sharing.df[i,]$keep_qtl | sharing.df[i,]$keep_gwas))
  sharing.df$keep_snp[is.na(sharing.df$keep_snp)] = FALSE

  plot.df$label = NA
  plot.df$label[plot.df$study == "GWAS" & plot.df$rsid == coloc_pair$gwas_rsid] = coloc_pair$gwas_rsid
  plot.df$label[plot.df$study == "QTL" & plot.df$rsid == coloc_pair$qtl_rsid] = coloc_pair$qtl_rsid
  
  if (match_snps_by_chrpos) {
    plot.df = plot.df %>% dplyr::left_join(sharing.df %>% dplyr::select(pos, variant, min_p, keep_snp, shared_rsid = rsid), by="pos") %>%
      mutate(rsid = if_else(!is.na(shared_rsid), shared_rsid, rsid)) %>%
      select(-shared_rsid)
  } else {
    plot.df = plot.df %>% dplyr::left_join(sharing.df %>% dplyr::select(rsid, variant, min_p, keep_snp), by="rsid")
  }
  
  if (plot_p_threshold < 1) {
    plot.df$keep_snp = plot.df$keep_snp & (plot.df$min_p < plot_p_threshold)
  }
  plot.df
}

colocOverlayPlot <- function(plot.df, coloc_pair, coloc_summary) {
  nPoints = nrow(plot.df)
  plot.df = plot.df %>% dplyr::filter(keep_snp)
  nShown = nrow(plot.df)
  
  xmax = max(plot.df$pos)
  xrange = xmax - min(plot.df$pos)
  H4ratio = coloc_summary$PP.H4.abf / (coloc_summary$PP.H3.abf + coloc_summary$PP.H4.abf)
  colocLabel = sprintf("H4 = %.3f\nH3 = %.3f\nH4/(H3+H4) = %.3f\n%d points shown of %d",
                       coloc_summary$PP.H4.abf, coloc_summary$PP.H3.abf, H4ratio, nShown, nPoints)
  plotTitle = sprintf("gwas %s, qtl %s: %s", coloc_pair$gwas_rsid, coloc_pair$qtl_rsid, coloc_pair$feature)
  
  colocPlot = ggplot(plot.df, aes(x=pos, y=-log10(p_nominal), color=study, shape=variant)) +
    geom_point(alpha=0.7) + ggtitle(plotTitle) +
    theme_bw() + scale_shape_manual(values=c(`Not shared`=4, Shared=19)) +
    geom_text(aes(label=label), color="black", size=3, hjust="left", nudge_x=(xrange/100)) +
    annotate("text", max(plot.df$pos), max(-log10(plot.df$p_nominal), na.rm=T), label=colocLabel, hjust=1, vjust=1, size=3)
  return(colocPlot)
}

colocLinesPlot <- function(plot.df, coloc_pair, coloc_summary) {
  nPoints = nrow(plot.df)
  plot.df = plot.df %>% dplyr::filter(keep_snp)
  nShown = nrow(plot.df)
  
  plotTitle = sprintf("gwas %s, qtl %s: %s", coloc_pair$gwas_rsid, coloc_pair$qtl_rsid, coloc_pair$feature)
  H4ratio = coloc_summary$PP.H4.abf / (coloc_summary$PP.H3.abf + coloc_summary$PP.H4.abf)
  colocLabel = sprintf("H4 = %.3f\nH3 = %.3f\nH4/(H3+H4) = %.3f\n%d points shown of %d",
                       coloc_summary$PP.H4.abf, coloc_summary$PP.H3.abf, H4ratio, nShown, nPoints)
  
  qtlYMax  = min(300, max(-log10(plot.df %>% dplyr::filter(study == "QTL") %>% .$p_nominal), na.rm = T))
  gwasYMax = min(300, max(-log10(plot.df %>% dplyr::filter(study == "GWAS") %>% .$p_nominal), na.rm = T))
  
  plot.df$log10p = -log10(plot.df$p_nominal)
  plot.df[plot.df$study == "QTL",]$log10p = plot.df[plot.df$study == "QTL",]$log10p * gwasYMax / qtlYMax
  
  colocPlot = ggplot(plot.df, aes(x=study, y=log10p, color=rsid, shape=variant, group=rsid)) +
    geom_jitter(data=plot.df, mapping=aes(x=study, y=log10p, group=study, color=rsid, shape=variant), alpha=0.7, width=0.1) +
    geom_line(alpha=0.3) +
    theme_bw() + ggtitle(plotTitle) + scale_shape_manual(values=c(`Not shared`=4, Shared=19)) +
    theme(legend.position = "None") +
    scale_y_continuous(name = "GWAS:  -log10(p)",
                       sec.axis = sec_axis(~ . * qtlYMax / gwasYMax, name = "QTL:  -log10(p)"), limits = c(0, gwasYMax)) +
    annotate("text", "QTL", max(plot.df$log10p, na.rm=T), label=colocLabel, hjust=0, vjust=1, size=3)
  
  return(colocPlot)
}

locusComparePlot = function(plot.df, plotTitle) {
  plot.shared.df = plot.df %>% filter(variant == "Shared")
  write(sprintf("Subsetting locusCompare plot to %d of %d total SNPs. There were %d QTL-only SNPs and %d GWAS-only SNPs.",
                as.integer(sum(plot.df$keep_snp)/2), as.integer(nrow(plot.df)/2),
                plot.df %>% filter(study == "QTL", variant != "Shared") %>% nrow(),
                plot.df %>% filter(study == "GWAS", variant != "Shared") %>% nrow()), stderr())
  leadSnpID = plot.shared.df %>% arrange(study, p_nominal) %>% .[1,"rsid"]
  write(sprintf("Lead SNP ID = %s", leadSnpID), stderr())
  stopifnot(leadSnpID %in% plot.df[plot.df$keep_snp,]$rsid)
  p = locuscompare(in_fn1 = plot.shared.df %>% filter(keep_snp, study == "GWAS") %>% select(rsid, pval = p_nominal),
                   in_fn2 = plot.shared.df %>% filter(keep_snp, study == "QTL") %>% select(rsid, pval = p_nominal),
                   snp = leadSnpID, title = 'GWAS', title2 = 'QTL') +
    theme(plot.title = element_text(size=12)) +
    ggtitle(plotTitle)
  p
}


###############################################################################

main()

