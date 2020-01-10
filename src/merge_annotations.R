#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
options(stringsAsFactors = F)

source("ppaFunctions.R")

args <- commandArgs(trailingOnly = TRUE)
locusFile = args[1]
gwasFile = args[2]
grch38File = args[3]
dbSNPFile = args[4]
vepFile = args[5]
vepImpactFile = args[6]
overlapFile = args[7]
jemeFile = args[8]
spliceaiFile = args[9]
phastConsFile = args[10]
phyloPFile = args[11]
GERPFile = args[12]
deepSEAFile = args[13]
captureSeqClosestExonFile = args[14]
captureSeqClosestSpliceFile = args[15]
gctaCondFile = args[16]
finemapFile1 = args[17]
finemapFile2 = args[18]
finemapFile3 = args[19]
outputBase = args[20]

# locusFile = "AD.loci.tsv"
# gwasFile = "summary_stats/AD.meta.assoc_loci.gz"
# grch38File = "annotated/AD.meta.assoc_loci.GRCh38.bed"
# dbSNPFile = "annotated/dbsnp_matching_table.tsv"
# vepFile = "annotated/AD.meta.vep_output.tsv.gz"
# vepImpactFile = "reference/vep_impact_table.tsv"
# overlapFile = "annotated/AD.meta.overlaps.txt"
# jemeFile = "annotated/AD.meta.overlaps.JEME.tsv"
# spliceaiFile = "annotated/AD.meta.spliceAI.tsv"
# phastConsFile = "annotated/AD.meta.phastCons100way.bed"
# phyloPFile = "annotated/AD.meta.phyloP100way.bed"
# GERPFile = "annotated/AD.meta.GERP_RS.bed"
# deepSEAFile = "annotated/DeepSEA/AD_meta_v5.DeepSEA.all.funsig"
# captureSeqClosestExonFile = "overlaps/AD.meta.RNACaptureSeq.closestExon.tsv.gz"
# captureSeqClosestSpliceFile = "overlaps/AD.meta.RNACaptureSeq.closestSplice.tsv.gz"
# gctaCondFile = "gcta/output_1e-5/cond/merged_loci.cond.out.flt.tsv"
# finemapFile1 = "finemap/output/AD.meta.finemap.ncausal_1.snp"
# finemapFile2 = "finemap/output/AD.meta.finemap.ncausal_2.snp"
# finemapFile3 = "finemap/output/AD.meta.finemap.ncausal_3.snp"
# outputBase = "annotated/AD.annotated"

collapseCol = function(df, dupcol, collapsecol, newcol) {
  dt = data.table(df)
  dtnew <- dt[!duplicated(df[,dupcol]), ]
  temp = dt[, paste(eval(as.name(collapsecol)), collapse=","), by=dupcol]
  dtnew[, newcol] <- temp$V1  
  as_tibble(dtnew)
}

###############################################################################
# Load files

loci.df = readr::read_tsv(locusFile) %>%
  mutate(lead_chrpos = paste(Chr, lead_pos, sep = "_"))

gwas.df = readr::read_tsv(gwasFile) %>%
  rename(chr=CHR, pos_hg37=BP, snp=SNP, Eff_allele=A1)

# Because the +/- 1Mb window around a couple of loci overlap, there
# are some duplicated SNPs in some of the files. To avoid having these
# get amplified with each left_join, we make sure that the relevant data
# in annotation files is unique.
grch38.df = readr::read_tsv(grch38File, col_names = c("chr", "pos_hg38", "dummy", "snp"), col_types = "ciic") %>%
  unique() %>%
  select(-chr, -dummy)

dbSNP.df = readr::read_tsv(dbSNPFile, col_names = T, col_types = cols(chrom="c")) %>%
  select(gwas_snpid, rsid, ref, alt) %>%
  unique()

# We initially allow multiple entries per GWAS SNP ID so that we can capture the
# most severe VEP consequence, but later remove duplicates.
vep.df = readr::read_tsv(vepFile, col_names = T, na = "-", guess_max = 100000,
                         col_types = cols(.default = col_character(), SYMBOL="c", EXON="c", INTRON="c", cDNA_position="c", CDS_position="c", Protein_position="c",
                                          PUBMED="c", CADD_phred="n", `GERP++_RS`="n", `M-CAP_score`="n", MetaLR_score="n",
                                          `fathmm-MKL_coding_score`="n", integrated_fitCons_score="n", phastCons100way_vertebrate="n",
                                          phastCons20way_mammalian="n", phyloP100way_vertebrate="n", phyloP20way_mammalian="n",
                                          CADD_PHRED="n", CADD_RAW="n")) %>%
  select(-CADD_phred) # There are 2 CADD columns; this is the CADD score covering only exons; we want the whole-genome one

overlap.df = readr::read_tsv(overlapFile) %>%
  rename(snp = SNP) %>%
  unique()

jeme.df = readr::read_tsv(jemeFile) %>%
  unique() %>%
  select(snp, geneScores, geneScoreDetails)


vep.impact.df = readr::read_tsv(vepImpactFile)
vep.impact.df$consequence_priority = 1:nrow(vep.impact.df)

# Some of the consequences include multiple terms, e.g.
# "missense_variant,splice_region_variant,NMD_transcript_variant"
# We want to use the consequence priority that corresponds with the most
# severe of these
getWorstConsequencePriority = function(consequence_str) {
  terms.df = data.frame(term = strsplit(consequence_str, ",", fixed=T)[[1]])
  terms.df = terms.df %>% dplyr::left_join(vep.impact.df, by="term")
  min(terms.df$consequence_priority)
}
vep.consequence.df = data.frame(term = unique(vep.df$Consequence))
vep.consequence.df$consequence_priority = sapply(vep.consequence.df$term, getWorstConsequencePriority)

# Get conservation annotations
phastCons.df = readr::read_tsv(phastConsFile, col_names = c("chr", "start", "end", "snp", "phastCons")) %>%
  select(snp, phastCons)
phyloP.df = readr::read_tsv(phyloPFile, col_names = c("chr", "start", "end", "snp", "phyloP")) %>%
  select(snp, phyloP)
GERP.df = readr::read_tsv(GERPFile, col_names = c("chr", "start", "end", "snp", "GERP_RS")) %>%
  select(snp, GERP_RS)

# There are some SNPs with more than one spliceAI entry. I'm not sure why this is
# in every case, but sometimes it appears that SpliceAI gives two scores depending
# on transcript annotations, such as when a SNP could be in an exon or not. In these
# cases we'll select the annotation that has a higher max_DS score.
spliceai.df = readr::read_tsv(spliceaiFile) %>%
  select(snp = SNP, max_DS, spliceai_strand, spliceai_type, starts_with("DS"), starts_with("DP")) %>%
  arrange(desc(max_DS)) %>%
  filter(!duplicated(snp))


deepSEA.df = readr::read_csv(deepSEAFile) %>%
  select(snp = name, DeepSEA_funsig = `Functional significance score`) %>%
  arrange(DeepSEA_funsig) %>%
  filter(!duplicated(snp))

captureSeq.closestExon.df = readr::read_tsv(captureSeqClosestExonFile) %>%
  select(snp, captureSeqExon = feature, captureSeqExonDist = featureDist)

captureSeq.closestSplice.df = readr::read_tsv(captureSeqClosestSpliceFile) %>%
  select(snp, captureSeqSplice = feature, captureSeqSpliceDist = featureDist)


# Read in SNP probabilities from GCTA conditioning on independent lead SNPs at each locus
gctaCond.df = readr::read_tsv(gctaCondFile) %>%
  rename(cond_lead_snp = cond_snp, pCond = pC, bCond = bC, bCond_se = bC_se, snp=SNP)

# SNP IDs in this file are e.g. rs1234_A_C, but we don't want the alleles there.
# (However, for indels we have rs1234_AT_A_AT_A, and we want the first pair of alleles.)
fixSnpID = function(id) {
  vals = str_split(id, "_")[[1]]
  if (length(vals) > 3) {
    newID = paste(vals[1:(length(vals)-2)], collapse="_")
  } else {
    newID = vals[1]
  }
  newID
}
gctaCond.df$snp = sapply(gctaCond.df$snp, FUN=fixSnpID)

# The cond_lead_snp column indicates the lead SNP after conditioning, but we want
# to add a column that indicates the other SNPs that were conditioned on. We
# can get this by removing the cond_lead_snp from the list of lead SNPs at the locus.
updateCondSnp = function(SNPstr, cond_lead_snp) {
  snps = str_split(SNPstr, pattern = ",")
  cond_snps = lapply(1:length(snps), function(i) setdiff(snps[[i]], cond_lead_snp[i]))
  sapply(cond_snps, function(x) paste(x, collapse = ","))
}

gctaCond.df = gctaCond.df %>%
  left_join(loci.df %>% select(locus=locus_name, n_snps, SNPs), by = "locus") %>%
  mutate(cond_snps = updateCondSnp(SNPs, cond_lead_snp))

# Get a finemapping probability for each conditional signal, using the WTCCC
# simple Bayesian method.
gctaCond.df$signal = paste(gctaCond.df$locus, gctaCond.df$cond_lead_snp)
#gctaCond.df$gcta_condProb = getPPAsFromPval(gctaCond.df %>% rename(F = freq) %>% arrange(signal),
#                                           segmentCol = "signal", pCol = "pCond", defaultN = 10000)
gctaCond.df$gcta_condProb = getPPAsFromBetaSE(gctaCond.df %>% rename(F = freq) %>% arrange(signal),
                                            segmentCol = "signal", betaCol = "bCond", seCol = "bCond_se")

# For each locus, collapse the conditional probabilities for each SNP to
# select that with the max probability
gctaCond.simple.df = gctaCond.df %>%
  group_by(locus, snp) %>%
  summarise(gcta_maxCondProb = max(gcta_condProb),
            gcta_pCond = pCond[which.max(gcta_condProb)],
            gcta_bCond = bCond[which.max(gcta_condProb)],
            gcta_bCondSE = bCond_se[which.max(gcta_condProb)],
            gcta_cond_snps = cond_snps[which.max(gcta_condProb)]) %>%
  ungroup()


# Load FINEMAP probabilities
getFinemapSnpID = function(snpid_a1_a2) {
  snpid_parts = str_split(snpid_a1_a2, "_", Inf)
  sapply(snpid_parts, FUN = function(s) paste(s[1:(length(s)-2)], collapse="_"))
}
finemap.ncausal1.df = readr::read_tsv(finemapFile1) %>%
  select(snp=rsid, finemap_prob1=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

finemap.ncausal2.df = readr::read_tsv(finemapFile2) %>%
  select(snp=rsid, finemap_prob2=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

finemap.ncausal3.df = readr::read_tsv(finemapFile3) %>%
  select(snp=rsid, finemap_prob3=prob) %>%
  mutate(snp=getFinemapSnpID(snp))

# finemap.ncausal4.df = readr::read_tsv(finemapFile4) %>%
#   select(snp=rsid, finemap_prob4=prob) %>%
#   mutate(snp=getFinemapSnpID(snp))


###############################################################################
# Join everything together!

# Start by joining the GWAS with dbSNP and VEP. Arrange by descending severity of
# VEP consequence, and remove duplicate entries per GWAS SNP.
ann.df = gwas.df %>%
  left_join(grch38.df, by="snp") %>%
  left_join(dbSNP.df, by=c("snp" = "gwas_snpid")) %>%
  left_join(vep.df, by=c("rsid" = "variant")) %>%
  left_join(vep.consequence.df, by=c("Consequence" = "term")) %>%
  arrange(consequence_priority) %>%
  filter(!duplicated(snp)) %>%
  select(-consequence_priority)

# Next, collapse duplicated dbSNP IDs and join again, so that we have all
# possibly relevant dbSNP IDs
dbSNP.collapsed.df = dbSNP.df %>%
  select(gwas_snpid, rsid) %>%
  collapseCol(dupcol = "gwas_snpid", collapsecol = "rsid", newcol = "rsids") %>%
  select(-rsid)

ann.df = ann.df %>%
  select(-rsid) %>%
  left_join(dbSNP.collapsed.df, by=c("snp"="gwas_snpid")) %>%
  left_join(overlap.df, by="snp") %>%
  left_join(jeme.df, by="snp") %>%
  left_join(phastCons.df, by="snp") %>%
  left_join(phyloP.df, by="snp") %>%
  left_join(GERP.df, by="snp") %>%
  left_join(spliceai.df, by="snp") %>%
  left_join(deepSEA.df, by="snp") %>%
  left_join(captureSeq.closestExon.df, by="snp") %>%
  left_join(captureSeq.closestSplice.df, by="snp")

# Add locus annotations
window = 630000
ann.df$lead_chrpos = ""
for (i in 1:nrow(loci.df)) {
  chr = loci.df[i,]$Chr
  start = loci.df[i,]$start - window
  end = loci.df[i,]$stop + window
  in_locus = (ann.df$chr == chr & start <= ann.df$pos_hg37 & ann.df$pos_hg37 <= end)
  ann.df[in_locus, ]$lead_chrpos = loci.df[i,]$lead_chrpos
}

ann.df = ann.df %>%
  left_join(loci.df %>% select(locus_name, locus_nSnps = n_snps, lead_chrpos, lead_p), by="lead_chrpos") %>%
  filter(!is.na(locus_nSnps)) %>%
  left_join(gctaCond.simple.df %>% select(-locus), by="snp") %>%
  left_join(finemap.ncausal1.df, by="snp") %>%
  left_join(finemap.ncausal2.df, by="snp") %>%
  left_join(finemap.ncausal3.df, by="snp") %>%
  arrange(lead_chrpos)

# Make a column that has the Finemap probabilities corresponding to the
# number of SNPs selected by GCTA for the locus
getFinemapProb = function(df) {
  nCausal = df$locus_nSnps[1]
  if (nCausal == 1) {
    df$finemap_prob1
  } else if (nCausal == 2) {
    df$finemap_prob2
  } else if (nCausal == 3) {
    df$finemap_prob3
  } else if (nCausal == 4) {
    df$finemap_prob4
  } else {
    rep(NA, nrow(df))
  }
}
ann.df$finemap_prob_nc = unlist(by(ann.df, ann.df$lead_chrpos, getFinemapProb))

ann.df = ann.df %>%
  arrange(chr, lead_chrpos, desc(finemap_prob_nc)) %>%
  select(lead_chrpos, locus_name, lead_p, chr, pos_hg37, pos_hg38, Eff_allele, A2, ref, alt, snp, rsids, freq=FREQ,
         META_P, META_BETA, GWAS_P, GWAS_BETA, GWAX_UKBB_P, GWAX_UKBB_BETA, DIRECT, I2, HET_P, FREQ, INFO,
         locus_nSnps, gcta_maxCondProb, starts_with("gcta"), finemap_prob_nc, starts_with("finemap"), Consequence, SYMBOL, Gene, BIOTYPE, PUBMED,
         phastCons, phyloP, GERP_RS, spliceai_max_DS=max_DS, DeepSEA_funsig,
         captureSeqExon, captureSeqExonDist, captureSeqSplice, captureSeqSpliceDist,
         MCAP_score=`M-CAP_score`, MetaLR_score, fathmm_MKL_coding_score=`fathmm-MKL_coding_score`,
         CADD_PHRED, integrated_confidence_value, integrated_fitCons_score, iPSC:ipsSensNeuron, DNase:`Blood & Immune Enh`,
         starts_with("geneScore"), everything())

gzf = gzfile(paste0(outputBase, ".all.tsv.gz"), "w")
write.table(ann.df, file = gzf, sep = "\t", quote = F, col.names=T, row.names=F, na = "")
close(gzf)

# Remove a bunch of columns that aren't interesting. Some of these are because dbNSFP only annotates
# coding variation, which for conservation (GERP, phylop, phastcons) isn't very interesting.
ann.df.sel = ann.df %>%
  select(-Location, -IMPACT, -CADD_RAW, -PHENO, -SOMATIC, -starts_with("gnomAD"), -SYMBOL_SOURCE, -HGNC_ID,
         -Existing_variation, -APPRIS, -SIFT, -PolyPhen, -CLIN_SIG, -FLAGS,
         -AA_AF, -EA_AF, -integrated_confidence_value, -`GERP++_RS`, -phastCons100way_vertebrate,
         -phastCons20way_mammalian, -phyloP100way_vertebrate, -phyloP20way_mammalian) %>%
  rename(distance_to_transcript=DISTANCE, feature_strand=STRAND)

gzf = gzfile(paste0(outputBase, ".selected.tsv.gz"), "w")
write.table(ann.df.sel, file = gzf, sep = "\t", quote = F, col.names=T, row.names=F, na = "")
close(gzf)

# Write a table with only the "probable" SNPs, including those with
# finemap probability > 0.0001 based on the expected number of causals,
# or with probability > 0.01 with other numbers of causals. But in all
# loci, select at least the top 20 variants per locus.
fname = paste0(outputBase, ".selected.probable.tsv")
ann.df.sel.probable = ann.df.sel %>%
  mutate(select_for_causal_prob = (gcta_maxCondProb > 0.001 | finemap_prob_nc > 0.0001 | finemap_prob1 > 0.01 | finemap_prob2 > 0.01 | finemap_prob3 > 0.01)) %>%
  arrange(chr, lead_chrpos, desc(select_for_causal_prob), desc(finemap_prob_nc)) %>%
  group_by(lead_chrpos) %>%
  mutate(index = row_number()) %>%
  filter((index <= 20 & locus_name != "APOE") | select_for_causal_prob) %>%
  ungroup() %>%
  select(-index, -select_for_causal_prob)
write.table(ann.df.sel.probable, file = fname, sep = "\t", quote = F, col.names=T, row.names=F, na = "")
