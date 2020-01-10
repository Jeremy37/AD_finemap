#!/usr/bin/env Rscript
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
jemeOverlapDir = args[1]
jemeAnnotDir = args[2]
snpBedFile = args[3]
outname = args[4]


if (!file.exists(file.path(jemeOverlapDir, "overlaps.annotated.JEME.txt"))) {
  # Read the metadata file for JEME Roadmap samples
  roadmap.meta = readr::read_csv(file.path(jemeAnnotDir, "roadmap.description.csv"))
  
  # We read in each file of enhancer-promoter connections and scores which overlapped
  # with a SNP in our set. We'll store a data.frame for the associated gene scores.
  files <- list.files(path=jemeOverlapDir, pattern="*.txt", full.names=T, recursive=FALSE)
  numFiles = length(files)
  
  linkscore.df = data.frame(key=character(0))
  x = lapply(files, function(fname) {
    id = strsplit(basename(fname), "\\.")[[1]][1]
    EID = roadmap.meta[match(id, roadmap.meta$id), ]$sampleName
    eid.df <- readr::read_tsv(fname, col_names = c("chr", "pos", "gene", "score"))
    eid.df$key = paste(eid.df$chr, eid.df$pos, eid.df$gene, sep="_")
    # Important when full_joining multiple tables that you don't have duplicate
    # values, otherwise these combinations grown exponentially
    eid.df = eid.df[!duplicated(eid.df$key),]
    eid.df2 = eid.df %>% dplyr::rename(!!EID := "score") %>% dplyr::select(-chr, -pos, -gene)
    linkscore.df <<- linkscore.df %>% dplyr::full_join(eid.df2, by="key")
  })
  rm(x)
  
  getTopScores = function(df, i, ntop=0) {
    names = colnames(df)[!is.na(df[i,])]
    valsdf = data.frame(eid=names, vals=as.numeric(df[i,names])) %>% arrange(-vals)
    if (ntop > 0) {
      ntop = min(ntop, nrow(valsdf)-1)
      valsdf = valsdf[1:ntop,]
    }
    valsString = paste(valsdf$eid, valsdf$vals, sep=",", collapse="/")
    valsString
  }
  
  linkscore.df$chr = sapply(linkscore.df$key, function(x) {strsplit(x, "_", fixed=T)[[1]][1]})
  linkscore.df$pos = sapply(linkscore.df$key, function(x) {strsplit(x, "_", fixed=T)[[1]][2]})
  linkscore.df$geneName = sapply(linkscore.df$key, function(x) {name = strsplit(x, "_", fixed=T)[[1]][3]; strsplit(name, "$", fixed=T)[[1]][2]})
  
  linkscore.df$topScore = apply(linkscore.df[,2:(numFiles+1)], 1, function(x) sum(max(x, na.rm=T)))
  
  linkscore.df$topScoreStr = NA
  for (i in 1:nrow(linkscore.df)) {
    linkscore.df[i,]$topScoreStr = getTopScores(linkscore.df[,2:(numFiles+1)], i, 1)
  }
  
  linkscore.df$numScores = apply(linkscore.df[,2:(numFiles+1)], 1, function(x) sum(!is.na(x)))
  
  numTopScores = 5
  linkscore.df$topScoreDetails = NA
  for (i in 1:nrow(linkscore.df)) {
    linkscore.df[i,]$topScoreDetails = getTopScores(linkscore.df[,2:(numFiles+1)], i, numTopScores)
  }
  
  snplinks.df = linkscore.df %>% distinct(chr, pos)
  snplinks.df$geneScores = NA
  for (i in 1:nrow(snplinks.df)) {
    snp.df = linkscore.df %>% dplyr::filter(chr == snplinks.df[i,]$chr, pos == snplinks.df[i,]$pos) %>% arrange(-topScore)
    snplinks.df$geneScores[i] = paste(snp.df$topScore, snp.df$geneName, sep=",", collapse=" / ")
  }
  
  snplinks.df$geneScoreDetails = NA
  for (i in 1:nrow(snplinks.df)) {
    snp.df = linkscore.df %>% dplyr::filter(chr == snplinks.df[i,]$chr, pos == snplinks.df[i,]$pos) %>% arrange(-topScore)
    snplinks.df$geneScoreDetails[i] = paste(snp.df$geneName, snp.df$topScoreDetails, sep=":", collapse=" ## ")
  }
  
  
  write.table(snplinks.df, file=file.path(jemeOverlapDir, "overlaps.annotated.JEME.txt"),
              sep="\t", quote=F, col.names=T, row.names=F)
} else {
  snplinks.df = readr::read_tsv(file.path(jemeOverlapDir, "overlaps.annotated.JEME.txt"))  
}
  
#input.df = readr::read_tsv(snpFile, col_types = cols(PUBMED = "c"))
input.df = readr::read_tsv(snpBedFile, col_types = "ciic", col_names = c("chr", "dummy", "pos", "snp"))
input.df$chrpos = paste(input.df$chr, input.df$pos, sep="_")

snplinks.df$chr = gsub("chr", "", snplinks.df$chr)
snplinks.df$chrpos = paste(snplinks.df$chr, snplinks.df$pos, sep="_")
output.df = input.df %>%
  dplyr::left_join(snplinks.df %>% dplyr::select(-chr, -pos), by="chrpos") %>%
  dplyr::select(-chrpos)

write.table(output.df, file=file.path(outname), sep="\t", quote=F, col.names=T, row.names=F, na = "")

