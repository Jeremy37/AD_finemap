library(tidyverse)

hitsnps <- readr::read_tsv("AD.loci.tsv") %>%
  group_by(lead_SNP) %>%
  mutate(lead_rsid = str_split(lead_SNP, "_")[[1]][1]) %>%
  ungroup()

# list of new genes to highlight
recent_genes = c("ADAMTS4","CLNK","ECHDC3","SPPL2A","ADAM10","VKORC1","PLCG2","SCIMP","ACE","APP-ADAMTS1","APH1B")
#recent_genes = c()
new_genes = c("CCDC6","TSPAN14","NCK2","SPRED2")
new_subthreshold = c("IKZF1","TMEM163","TSPOAP1")

assoc <- read.table("summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.col_subset.tsv.bgz",
                    header = T,
                    colClasses = c("numeric","numeric","character","numeric"))

chroms <- seq(1,22)
# get chromosome positions
chrom_starts <- c()
chrom_stops <- c()
for (i in seq(1,22)) {
  chrom_starts <- append(chrom_starts, min(assoc$BP[which(assoc$CHR == i)]))
  chrom_stops <- append(chrom_stops, max(assoc$BP[which(assoc$CHR == i)]))
}
chrom_adjust <- c(0)
for (i in seq(2,22)) {
  adj <- chrom_stops[i-1] + chrom_adjust[i-1]
  chrom_adjust <- append(chrom_adjust,adj)
}

assoc_orig = assoc

getWindowSnps = function(assoc, chr, pos, windowsize) {
  assoc %>% filter(CHR == chr, abs(BP - pos) < windowsize)
}

assoc_nearhits = bind_rows(lapply(1:nrow(hitsnps), function(i) getWindowSnps(assoc_orig, hitsnps$Chr[i], hitsnps$lead_pos[i], 50000)))

assoc_sample = assoc_orig %>% filter(META_P >= 0.0001 & META_P < 0.1)
# Sample a fixed number of SNPs with probability of selection
# depending on the p value
select_indexes = sample(1:nrow(assoc_sample), size = 8e4, replace = F, 1 / sqrt(assoc_sample$META_P))
assoc_sample = assoc_sample[select_indexes,]

assoc_sample2 = assoc_orig %>% filter(META_P >= 0.1)
select_indexes2 = sample(1:nrow(assoc_sample2), size = 2e4, replace = F)
assoc_sample2 = assoc_sample2[select_indexes2,]

#quantile(assoc_sample$META_P, probs = seq(0,1,0.01), na.rm = T)

assoc = bind_rows(assoc_orig %>% filter(META_P < 0.0001),
                  assoc_nearhits,
                  assoc_sample,
                  assoc_sample2 ) %>%
  filter(!duplicated(SNP))

quantile(assoc$META_P, probs = seq(0,1,0.01), na.rm = T)


# colors
plot_cols <- rep(c("#578196","#005b7f"),12)
#annot_col <- "#2573ba"

pdf("AD_meta_v5_annotated_manhattan.pdf",
    width = 9,
    height = 4, pointsize = 9)

doplot = function() {
  annot_col_new <- "blue"
  annot_col_recent <- "black"
  annot_col_subthreshold <- "darkred"
  plot(c(min(assoc$BP[which(assoc$CHR == 1)]),
         max(chrom_adjust) + max(assoc$BP[which(assoc$CHR == 22)])),
       c(0,35),
       type = "n",
       xlab = "Chromosome",
       ylab = expression(-log[10](italic(p))),
       axes = F)
  par(cex = 0.75)
  
  
  #for (i in seq(1,22)) {
  for (i in seq(1,22)) {
    print(i)
    chrom_assoc <- assoc[which(assoc$CHR == i),]
    # chrom_assoc <- chrom_assoc[which(chrom_assoc$META_P < 0.005),]
    # get list of SNPs to highlight
    highlight_snps_new <- c()
    highlight_snps_recent <- c()
    highlight_snps_subthreshold <- c()
    hitsnps_chrom <- hitsnps[hitsnps$Chr == i,]
    if (nrow(hitsnps_chrom) > 0) {
      for (j in 1:nrow(hitsnps_chrom)) {
        if (hitsnps_chrom$locus_name[j] %in% new_genes) {
          hit_pos <- hitsnps_chrom$lead_pos[j]
          highlight_i <- intersect(which(chrom_assoc$BP > (hit_pos - 100000)), which(chrom_assoc$BP < (hit_pos + 100000)))
          highlight_snps_new <- c(highlight_snps_new, chrom_assoc$SNP[highlight_i])
        }
        if (hitsnps_chrom$locus_name[j] %in% recent_genes) {
          hit_pos <- hitsnps_chrom$lead_pos[j]
          highlight_i <- intersect(which(chrom_assoc$BP > (hit_pos - 100000)), which(chrom_assoc$BP < (hit_pos + 100000)))
          highlight_snps_recent <- c(highlight_snps_recent, chrom_assoc$SNP[highlight_i])
        }
        if (hitsnps_chrom$locus_name[j] %in% new_subthreshold) {
          hit_pos <- hitsnps_chrom$lead_pos[j]
          highlight_i <- intersect(which(chrom_assoc$BP > (hit_pos - 100000)), which(chrom_assoc$BP < (hit_pos + 100000)))
          highlight_snps_subthreshold <- c(highlight_snps_subthreshold, chrom_assoc$SNP[highlight_i])
        }
      }
    }
    new_snp_ind = which(chrom_assoc$SNP %in% highlight_snps_new)
    recent_snp_ind = which(chrom_assoc$SNP %in% highlight_snps_recent)
    subthreshold_snp_ind = which(chrom_assoc$SNP %in% highlight_snps_subthreshold)
    highlight_snp_ind = c(new_snp_ind, recent_snp_ind, subthreshold_snp_ind)
    if (length(highlight_snp_ind) > 0) {
      points(chrom_assoc$BP[-highlight_snp_ind]+chrom_adjust[i],
             -log10(chrom_assoc$META_P[-highlight_snp_ind]),
             pch = 20,
             col = plot_cols[i])
      points(chrom_assoc$BP[recent_snp_ind]+chrom_adjust[i],
             -log10(chrom_assoc$META_P[recent_snp_ind]),
             col = annot_col_recent,
             pch = 20)
      points(chrom_assoc$BP[subthreshold_snp_ind]+chrom_adjust[i],
             -log10(chrom_assoc$META_P[subthreshold_snp_ind]),
             col = annot_col_subthreshold,
             pch = 20)
      points(chrom_assoc$BP[new_snp_ind]+chrom_adjust[i],
             -log10(chrom_assoc$META_P[new_snp_ind]),
             col = annot_col_new,
             pch = 20)
    } else {
      points(chrom_assoc$BP+chrom_adjust[i],
             -log10(chrom_assoc$META_P),
             pch = 20,
             col = plot_cols[i])
      
    }
  }
  
  # adjust positions of some locus names to avoid overlaps
  #adjust_pos <- rep(0,nrow(hitsnps))
  adjust_pos <- c(0,0, # ADAMTS4, CR1
                  -20000000, # SPRED2,
                  -15000000, 0, # NCK2, BIN1
                  25000000, # TMEM163
                  0,0,0,0,0,0,0,0,0, # INPP5D, CLNK, HLA, TREM2, CD2AP, IKZF1, PILRA, EPHA1, PTK2B-CLU
                  -20000000, # ECHDC3
                  -15000000, # CCDC6
                  10000000, # TSPAN14
                  0,0,0,0,0,0, # SPI1, MS4A4A, PICALM, SORL1, FERMT2, SLC24A4
                  -60000000, # SPPL2A
                  -30000000, # ADAM10
                  0, # APH1B
                  -10000000, # VKORC1
                  -12500000, # PLCG2
                  12500000, # SCIMP
                  10000000, # MIR142
                  40000000, # ACE
                  0,0,0,0,0 # ABCA7, APOE, CD33, CASS4, APP
  ) * 1.5
  
  assoc_lead = assoc %>% inner_join(hitsnps %>% select(lead_rsid, lead_p), by=c("SNP"="lead_rsid"))
  
  # add gene names
  for (i in 1:nrow(hitsnps)) {
    gene <- hitsnps$locus_name[i]
    if (gene %in% new_genes | gene %in% recent_genes | gene %in% new_subthreshold) {
      lead_snp <- hitsnps$lead_rsid[i]
      chrom <- hitsnps$Chr[i]
      bp <- hitsnps$lead_pos[i] + chrom_adjust[chrom]
      p <- assoc_lead$META_P[which(assoc_lead$SNP == lead_snp)]
      bp_adj <- bp + adjust_pos[i]
      #annot_col = "#005b7f"
      annot_col = "black"
      if (gene %in% new_genes) {
        annot_col = "blue"
      }
      if (gene %in% new_subthreshold) {
        annot_col = "brown"
      }
      yheight = 20
      text_cex = 1
      if (gene %in% new_genes) {
        text_cex = 1.1
      }
      text(bp_adj, yheight,
           gene,
           font = 3,
           srt = 90,
           col = annot_col,
           adj = 0,
           cex = text_cex)
      lines(c(bp,bp),
            c(-log10(p)+0.5,15),
            col = annot_col)
      lines(c(bp,bp_adj),
            c(15,yheight - 0.5),
            col = annot_col)
    }
  }
  
  axis(2)
  abline(h = -log10(5e-8),
         lty = "dashed",
         col = "red")
  
  # chromosome labels
  uneven_pos <- c()
  even_pos <- c()
  for (i in 1:22) {
    if (i == 22) {
      chrom_mid <- (chrom_adjust[i] + max(assoc$BP[which(assoc$CHR == i)]) + chrom_adjust[i])/2
    } else {
      chrom_mid <- (chrom_adjust[i] + chrom_adjust[i+1])/2
    }
    if (i %% 2 == 1) {
      uneven_pos <- append(uneven_pos,chrom_mid)
    }
    if (i %% 2 == 0) {
      even_pos <- append(even_pos,chrom_mid)
    }
  }
  
  axis(1,at = uneven_pos,
       labels = seq(1,21,length.out = 11),
       tick = F,
       line = -0.5)
  axis(1,at = even_pos,
       labels = seq(2,22,length.out = 11),
       tick = F)
  
  axis(1,at = append(chrom_adjust,chrom_adjust[22] + max(assoc$BP[which(assoc$CHR == 22)])),
       labels = rep("",23))
}
doplot()

dev.off()

png("AD_meta_v5_annotated_manhattan.png",
    width = 9,
    height = 4,
    units = "in",
    res = 300)
doplot()

dev.off()
