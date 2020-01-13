library(motifbreakR)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP142.GRCh37")
library(SNPlocs.Hsapiens.dbSNP142.GRCh37) # dbSNP142 in hg19

#BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
#library(SNPlocs.Hsapiens.dbSNP151.GRCh38) # dbSNP151 in hg38

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg19)
#library(BSgenome.Hsapiens.UCSC.hg38)
data(motifbreakR_motif)
data(hocomoco)

#setwd("/Users/jeremys/work/opentargets/AD_finemap/")
setwd("/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/AD_finemap/")

ad.df = readr::read_tsv("annotated/AD.IGAP1_GWAX_v3.annotated.selected.probable.paintor.tsv")
ad.df$rsids[is.na(ad.df$rsids)] = ad.df$snp[is.na(ad.df$rsids)]
rsids = ad.df$rsids[1:10]

#rsids = c("rs4727455", "rs6733839", "rs11218343", "rs12444183", "rs7810606")

variants <- snps.from.rsid(rsid = rsids,
                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19)

data(hocomoco)
data(hocomoco)
breakr.results <- motifbreakR(snpList = variants, filterp = TRUE,
                       pwmList = motifbreakR_motif,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

# breakr.rs4727455 <- breakr.results[names(breakr.results) %in% "rs4727455"]
# breakr.rs6733839 <- breakr.results[names(breakr.results) %in% "rs6733839"]
# breakr.rs11218343 <- breakr.results[names(breakr.results) %in% "rs11218343"]
# breakr.rs12444183 <- breakr.results[names(breakr.results) %in% "rs12444183"]
# breakr.rs7810606 <- breakr.results[names(breakr.results) %in% "rs7810606"]
# 
# breakr.rs4727455 <- calculatePvalue(breakr.rs4727455)

results.df = elementMetadata(breakr.results)
results.df$rsid = names(breakr.results)
results.df = results.df %>% select(rsid, everything)

write.table(results.df, "motifbreakr.results.tsv", quote=F, sep="\t", row.names=F, col.names=T, na="")

# pdf("motifbreakr.results.rs4727455.pdf", width=8, height=8)
# plotMB(results = results, rsid = "rs4727455", effect = "strong")
# dev.off()

pdf("motifbreakr.results.5snps.pdf", width=8, height=8)
for (rsid in rsids) {
  plotMB(results = results, rsid = rsid, effect = "strong")
}
dev.off()



# pca.snps.file <- system.file("extdata", "pca.enhancer.snps", package = "motifbreakR")
# pca.snps <- as.character(read.table(pca.snps.file)[,1])
# snps.mb <- snps.from.rsid(rsid = pca.snps[1:5],
#                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
#                           search.genome = BSgenome.Hsapiens.UCSC.hg19)
# 
# 
# results <- motifbreakR(snpList = snps.mb[1:5], filterp = TRUE,
#                        pwmList = hocomoco,
#                        threshold = 1e-4,
#                        method = "ic",
#                        bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                        BPPARAM = BiocParallel::bpparam())
# 
# 
