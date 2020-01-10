library(tidyverse)
library(optparse)
library(pheatmap)

option_list <- list(
  make_option(c("--file"), type="character", default=NULL,
              help="Input file of annotations"),
  make_option(c("--minFraction"), type="double", default=0.001,
              help="Minimum fraction of SNPs with nonzero annotation value for the annotation to be included"),
  make_option(c("--out"), type="character", default="paintor.annotations",
              help="Name for PDF output file")
)
opt <- parse_args(OptionParser(option_list=option_list))
# opt <- list(file = "paintor/all_annotations.txt", minFraction = 0.001, out = "paintor/paintor_anotation_correlations")

ann.df = readr::read_delim(opt$file, delim=" ")
ann.sums = colSums(ann.df)
#View(ann.sums)

pdf(opt$out, width=8, height=7)

ann.summary.df = data.frame(annotation = names(ann.sums), 
                            n = ann.sums,
                            frac = ann.sums / nrow(ann.df),
                            remove = (ann.sums / nrow(ann.df) < opt$minFraction))
ggplot(ann.summary.df, aes(x=fct_reorder(annotation, n, .desc=T), y=n, fill=remove)) + geom_bar(stat="identity") +
  theme_bw(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  scale_y_log10() + xlab("Annotation") + ylab("Number of SNPs") +
  scale_fill_manual(values=c("TRUE" = "red", "FALSE" = "grey30"))

ggplot(ann.summary.df, aes(x=fct_reorder(annotation, n, .desc=T), y=frac, fill=remove)) + geom_bar(stat="identity") +
  theme_bw(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  xlab("Annotation") + ylab("Fraction of SNPs") +
  scale_fill_manual(values=c("TRUE" = "red", "FALSE" = "grey30"))


# Remove annotations which cover too few variants
colsToRemove = (colSums(ann.df) / nrow(ann.df) < opt$minFraction)
if (sum(colsToRemove) > 0) {
  print("Removing columns with low annotation fraction. Annotation sums:")
  print(ann.sums[colsToRemove])
  ann.df = ann.df %>% select_if(function(x) sum(x)/length(x) >= opt$minFraction)
}

ann.cor.mat = ann.df %>%
  as.matrix %>%
  cor %>%
  as.data.frame

pheatmap(ann.cor.mat, cluster_rows = T, cluster_cols = T, main = "Annotation correlations", fontsize = 8)
# Manually examine the heatmap. Within sets of correlated annotations,
# check which annotations have the highest model probability from a run
# of PAINTOR with the single annotation.

# View the highly correlated annotation pairs
ann.cor.df = ann.cor.mat %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1) %>% filter(!is.na(value), var1 != var2)
#View(ann.cor.df)

dev.off()
