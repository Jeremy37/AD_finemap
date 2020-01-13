library(optparse)
library(tidyverse)
options(stringsAsFactors = F)

# Read loci from input file, and write input file of loci for PAINTOR
option_list <- list(
  make_option(c("--locus_file"), type="character", default=NULL,
              help="Input file with list of independent associated loci (produced by earlier pipeline steps)")
)
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- list(locus_file = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/AD_finemap/AD.1e-5_loci.tsv", outdir = "finemap/input/")

loci.df = readr::read_tsv(opt$locus_file)
loci.df$locus = paste0(loci.df$Chr, "_", loci.df$lead_pos, ".IGAP1_GWAX")

write.table(loci.df %>% select(locus), file="", col.names = F, row.names=F, sep="\t", quote=F)
