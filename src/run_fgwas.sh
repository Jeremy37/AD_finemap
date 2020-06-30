ROOT=/path/to/AD_finemap
SRC=$ROOT/src
GWAS_NAME=AD.meta
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen2/jeremys
ROOT=$JS/AD_finemap
SRC=$ROOT/src
GWAS_NAME=AD.meta

cd $ROOT

###############################################################################
# Convert AD_meta summary stats to fgwas input format
mkdir fgwas
mkdir fgwas/input
mkdir fgwas/input/tmp
mkdir fgwas/out
# fgwas input columns:
# SNPID	CHR	POS	Z	SE	F	NCASE	NCONTROL
#
# Num cases: 21982+(54939+898)/4 = 35941
# Num controls: 41944+(355900)/4 = 130919
#(echo -e "SNPID\tCHR\tPOS\tZ\tSE\tF\tNCASE\tNCONTROL"; \
# tabix summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz 1:60000000-70000000 | awk 'BEGIN{OFS="\t"}{print $3,$1,$2,$12/$13,$13,$18,35941,130919}') \
# | sort -k2,2n -k3,3n -k4,4nr | awk '!seen[$3]++' \
# | bgzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.tsv.bgz
(echo -e "SNPID\tCHR\tPOS\tZ\tSE\tF\tNCASE\tNCONTROL"; \
 (zcat summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{if ($13 > 0 && $13 != "nan") {print $3,$1,$2,$12/$13,$13,$18,35941,130919}}' \
 | sort -k2,2 -k3,3n -k4,4nr | awk '!seen[$3]++') ) \
 | bgzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.tsv.bgz

# Convert to bed file for intersection with annotations
#zcat fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{print $2,$3,$3}' | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.bed.gz
#FGWAS_BASE_TSV=$ROOT/fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.tsv.bgz
#FGWAS_BASE_BED=$ROOT/fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.bed.gz

# Exclude the APOE locus so that it doesn't dominate the results due to its very strong association.
(echo -e "SNPID\tCHR\tPOS\tZ\tSE\tF\tNCASE\tNCONTROL"; \
 (zcat summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' \
 | awk 'BEGIN{OFS="\t"}{if ($13 > 0 && $13 != "nan" && ($1 != "19" || $2 < 44000000 || $2 > 47000000)) {print $3,$1,$2,$12/$13,$13,$18,35941,130919}}' \
 | sort -k2,2 -k3,3n -k4,4nr | awk '!seen[$3]++') ) \
 | bgzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.noAPOE.fgwas.tsv.bgz

# Convert to bed file for intersection with annotations
zcat fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.noAPOE.fgwas.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{print $2,$3,$3}' | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.noAPOE.bed.gz
FGWAS_BASE_TSV=$ROOT/fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.noAPOE.fgwas.tsv.bgz
FGWAS_BASE_BED=$ROOT/fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.noAPOE.bed.gz


# Make annotation for whether a SNP is nearest to a gene with network score > Xth %ile of page.rank
# First prepare a bed file with network scores of genes. This has been done above in
# the network_analysis.Rmd script.
cat network/ad_network.all.bed | awk '$1 !~ /X|Y|MT|HG/' | sort -k 1,1 -k2,2n > fgwas/input/tmp/ad_network.all.bed
# Get the closest gene to each SNP. The awk line '!seen[$1"\t"$2]++' resolves SNPs
# nearest to more than one gene by just picking the first match. The -d option
# produces an extra line with the distance to the gene. We use this below to make
# another annotation which is just overlapping a gene, rather than nearest.
bedtools closest -a $FGWAS_BASE_BED -b fgwas/input/tmp/ad_network.all.bed -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.network.nearest.tsv.gz

# Check the input SNP file and annotation file have matching positions
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'

# Make the annotated input file
(paste <(gzhead 1 $FGWAS_BASE_TSV) \
       <(echo -e "network_50_60_nearest\tnetwork_60_70_nearest\tnetwork_70_80_nearest\tnetwork_80_90_nearest\tnetwork_90_95_nearest\tnetwork_gt_80_nearest\tnetwork_gt_90_nearest\tnetwork_gt_95_nearest\tnetwork_gt_80_overlapping\tnetwork_gt_90_overlapping\tgene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 50 && $8 < 60) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 60 && $8 < 70) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 70 && $8 < 80) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 80 && $8 < 90) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 90 && $8 < 95) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 80) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 90) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 95) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 80 && $9 == 0) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($8 >= 90 && $9 == 0) {print "1"} else {print "0"}}' ) \
       <(zcat fgwas/input/tmp/AD.snps.network.nearest.tsv.gz | awk '{if ($9 == 0) {print "1"} else {print "0"}}' ) ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz
 gzhead fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.head.tsv

###############################################################################
# Prepare gene expression annotations. This uses the gene expression files calculated
# earlier for single-cell brain expression and bulk RNA-seq from the QTL datasets.
# We want to use only genes present in the network so that enrichments are comparable.
$SRC/get_gene_expr_for_fgwas.Rmd

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.gtex.relative.gt_80.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_80.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.brain.relative.gt_80.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_80.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.eqtl_catalogue.relative.gt_80.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_80.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.sc.brain.relative.gt_80.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_80.tsv.gz

# Check the input SNP file and annotation file have matching positions
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_80.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_80.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_80.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_80.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'

# Make the annotated input files
# bulk expression (gtex) > 80th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.gtex.relative.gt_80.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_80.tsv.gz | cut -f 8-56) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_80.tsv.gz | awk '{if ($57 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.gt_80.tsv.gz

# bulk expression (brain) > 80th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.brain.relative.gt_80.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_80.tsv.gz | cut -f 8-20) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_80.tsv.gz | awk '{if ($21 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.brain.gt_80.tsv.gz

# bulk expression (eqtl_catalogue) > 80th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.eqtl_catalogue.relative.gt_80.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_80.tsv.gz | cut -f 8-56) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_80.tsv.gz | awk '{if ($57 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.eqtl_catalogue.gt_80.tsv.gz

# SC expression (brain) > 80th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.sc.brain.relative.gt_80.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_80.tsv.gz | cut -f 8-25) \
       <(zcat fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_80.tsv.gz | awk '{if ($26 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_80.tsv.gz


bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.gtex.relative.gt_90.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_90.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.brain.relative.gt_90.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_90.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.eqtl_catalogue.relative.gt_90.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_90.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.sc.brain.relative.gt_90.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_90.tsv.gz

# bulk expression (gtex) > 90th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.gtex.relative.gt_90.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_90.tsv.gz | cut -f 8-56) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.relative.nearest.gt_90.tsv.gz | awk '{if ($57 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.gt_90.tsv.gz

# bulk expression (brain) > 90th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.brain.relative.gt_90.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_90.tsv.gz | cut -f 8-20) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.brain.relative.nearest.gt_90.tsv.gz | awk '{if ($21 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.brain.gt_90.tsv.gz

# bulk expression (eqtl_catalogue) > 90th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.eqtl_catalogue.relative.gt_90.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_90.tsv.gz | cut -f 8-56) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.eqtl_catalogue.relative.nearest.gt_90.tsv.gz | awk '{if ($57 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.eqtl_catalogue.gt_90.tsv.gz

# SC expression (brain) > 90th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.sc.brain.relative.gt_90.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_90.tsv.gz | cut -f 8-25) \
       <(zcat fgwas/input/tmp/AD.snps.sc.brain.relative.nearest.gt_90.tsv.gz | awk '{if ($26 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_90.tsv.gz


################################################################################
# Do the same but for just "normalised" (within cell type) expression, not relative
# (normalised across cell types) expression
bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.gtex.norm.gt_80.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_80.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.sc.brain.norm.gt_80.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_80.tsv.gz

# Check the input SNP file and annotation file have matching positions
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_80.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'
paste <(zcat $FGWAS_BASE_TSV | cut -f 3 | sed '1d') <(zcat fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_80.tsv.gz | cut -f 2) | awk '{print $1, $2}' | awk '{if ($1 != $2) {print "Non-matching lines!"}}'

# Make the annotated input files
# bulk expression (gtex) > 80th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.gtex.relative.gt_80.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_80.tsv.gz | cut -f 8-56) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_80.tsv.gz | awk '{if ($57 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.norm.gt_80.tsv.gz

# SC expression (brain) > 80th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.sc.brain.relative.gt_80.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_80.tsv.gz | cut -f 8-25) \
       <(zcat fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_80.tsv.gz | awk '{if ($26 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.norm.gt_80.tsv.gz


bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.bulk.gtex.norm.gt_90.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_90.tsv.gz

bedtools closest -a $FGWAS_BASE_BED -b fgwas/expr.sc.brain.norm.gt_90.bed.gz -d \
  | awk '!seen[$1"\t"$2]++' | gzip > fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_90.tsv.gz

# bulk expression (gtex) > 90th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.bulk.gtex.relative.gt_90.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_90.tsv.gz | cut -f 8-56) \
       <(zcat fgwas/input/tmp/AD.snps.bulk.gtex.norm.nearest.gt_90.tsv.gz | awk '{if ($57 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.norm.gt_90.tsv.gz

# SC expression (brain) > 90th pctile
(paste <(gzhead 1 $FGWAS_BASE_TSV) fgwas/expr.sc.brain.relative.gt_90.colnames.tsv <(echo "gene_overlapping");
 paste <(zcat $FGWAS_BASE_TSV | sed '1d') <(zcat fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_90.tsv.gz | cut -f 8-25) \
       <(zcat fgwas/input/tmp/AD.snps.sc.brain.norm.nearest.gt_90.tsv.gz | awk '{if ($26 == 0) {print "1"} else {print "0"}}') ) \
 | gzip > fgwas/input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.norm.gt_90.tsv.gz




#######################################
# Run fgwas to get enrichments
cd fgwas
ANN=network_gt_80_nearest
submitJobs.py --MEM 5000 -j fgwas.$ANN -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN -o out/fgwas.$ANN"

ANN=network_gt_90_nearest
submitJobs.py --MEM 5000 -j fgwas.$ANN -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN -o out/fgwas.$ANN"

ANN=network_gt_80_overlapping
submitJobs.py --MEM 5000 -j fgwas.$ANN -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN -o out/fgwas.$ANN"

ANN=network_gt_90_overlapping
submitJobs.py --MEM 5000 -j fgwas.$ANN -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN -o out/fgwas.$ANN"

ANN=network_gt_80_overlapping
submitJobs.py --MEM 5000 -j fgwas.with_gene_overlap.$ANN -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN+gene_overlapping -o out/fgwas.with_gene_overlap.$ANN"

ANN=network_gt_90_overlapping
submitJobs.py --MEM 5000 -j fgwas.with_gene_overlap.$ANN -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN+gene_overlapping -o out/fgwas.with_gene_overlap.$ANN"

ANNOTS=(network_50_60_nearest network_60_70_nearest network_70_80_nearest network_80_90_nearest network_90_95_nearest network_gt_80_nearest network_gt_90_nearest network_gt_95_nearest network_gt_80_overlapping network_gt_90_overlapping gene_overlapping)
for ANN in "${ANNOTS[@]}"; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.$ANN -q yesterday --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w $ANN -o out/fgwas.single.$ANN"
done
submitJobs.py --MEM 5000 -j fgwas.network.combined -q yesterday \
    -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.network.tsv.gz -cc -w network_50_60_nearest+network_60_70_nearest+network_70_80_nearest+network_80_90_nearest+network_90_95_nearest+network_gt_95_nearest -o out/fgwas.network.combined"



cd fgwas
mkdir out/bulk_gtex
mkdir out/bulk_brain
mkdir out/bulk_eqtl_catalogue
mkdir out/sc_brain
mkdir out/bulk_gtex_norm
mkdir out/sc_brain_norm

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.gt_80.tsv.gz -cc -w $ANN -o out/bulk_gtex/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.gtex.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.brain.gt_80.tsv.gz -cc -w $ANN -o out/bulk_brain/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.brain.relative.gt_80.colnames.tsv | tr '\t' '\n')
grep "Successfully" FarmOut/fgwas.nearest.single*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code" FarmOut/fgwas.nearest.single*.txt | wc -l

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.eqtl_catalogue.gt_80.tsv.gz -cc -w $ANN -o out/bulk_eqtl_catalogue/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.eqtl_catalogue.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_80.tsv.gz -cc -w $ANN -o out/sc_brain/fgwas.nearest.single.$ANN"
done < <(cat expr.sc.brain.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.norm.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.norm.gt_80.tsv.gz -cc -w $ANN -o out/bulk_gtex_norm/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.gtex.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.norm.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.norm.gt_80.tsv.gz -cc -w $ANN -o out/sc_brain_norm/fgwas.nearest.single.$ANN"
done < <(cat expr.sc.brain.relative.gt_80.colnames.tsv | tr '\t' '\n')

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.nearest_gt_80.single.enrichments.tsv
cat out/bulk_gtex/fgwas.nearest.single*gt_80.params | grep "gt_80" | tr ' ' '\t' >> out/expr.bulk.gtex.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.brain.nearest_gt_80.single.enrichments.tsv
cat out/bulk_brain/fgwas.nearest.single*gt_80.params | grep "gt_80" | tr ' ' '\t' >> out/expr.bulk.brain.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.eqtl_catalogue.nearest_gt_80.single.enrichments.tsv
cat out/bulk_eqtl_catalogue/fgwas.nearest.single*gt_80.params | grep "gt_80" | tr ' ' '\t' >> out/expr.bulk.eqtl_catalogue.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.nearest_gt_80.single.enrichments.tsv
cat out/sc_brain/fgwas.nearest.single*gt_80.params | grep "gt_80" | tr ' ' '\t' >> out/expr.sc_brain.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.norm.nearest_gt_80.single.enrichments.tsv
cat out/bulk_gtex_norm/fgwas.nearest.single*gt_80.params | grep "gt_80" | tr ' ' '\t' >> out/expr.bulk.gtex.norm.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.norm.nearest_gt_80.single.enrichments.tsv
cat out/sc_brain_norm/fgwas.nearest.single*gt_80.params | grep "gt_80" | tr ' ' '\t' >> out/expr.sc_brain.norm.nearest_gt_80.single.enrichments.tsv


while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.gt_90.tsv.gz -cc -w $ANN -o out/bulk_gtex/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.gtex.relative.gt_90.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.brain.gt_90.tsv.gz -cc -w $ANN -o out/bulk_brain/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.brain.relative.gt_90.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.eqtl_catalogue.gt_90.tsv.gz -cc -w $ANN -o out/bulk_eqtl_catalogue/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.eqtl_catalogue.relative.gt_90.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_90.tsv.gz -cc -w $ANN -o out/sc_brain/fgwas.nearest.single.$ANN"
done < <(cat expr.sc.brain.relative.gt_90.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.norm.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.gtex.norm.gt_90.tsv.gz -cc -w $ANN -o out/bulk_gtex_norm/fgwas.nearest.single.$ANN"
done < <(cat expr.bulk.gtex.relative.gt_90.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.single.norm.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.norm.gt_90.tsv.gz -cc -w $ANN -o out/sc_brain_norm/fgwas.nearest.single.$ANN"
done < <(cat expr.sc.brain.relative.gt_90.colnames.tsv | tr '\t' '\n')

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.nearest_gt_90.single.enrichments.tsv
cat out/bulk_gtex/fgwas.nearest.single*gt_90.params | grep "gt_90" | tr ' ' '\t' >> out/expr.bulk.gtex.nearest_gt_90.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.brain.nearest_gt_90.single.enrichments.tsv
cat out/bulk_brain/fgwas.nearest.single*gt_90.params | grep "gt_90" | tr ' ' '\t' >> out/expr.bulk.brain.nearest_gt_90.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.eqtl_catalogue.nearest_gt_90.single.enrichments.tsv
cat out/bulk_eqtl_catalogue/fgwas.nearest.single*gt_90.params | grep "gt_90" | tr ' ' '\t' >> out/expr.bulk.eqtl_catalogue.nearest_gt_90.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.nearest_gt_90.single.enrichments.tsv
cat out/sc_brain/fgwas.nearest.single*gt_90.params | grep "gt_90" | tr ' ' '\t' >> out/expr.sc_brain.nearest_gt_90.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.norm.nearest_gt_90.single.enrichments.tsv
cat out/bulk_gtex_norm/fgwas.nearest.single*gt_90.params | grep "gt_90" | tr ' ' '\t' >> out/expr.bulk.gtex.norm.nearest_gt_90.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.norm.nearest_gt_90.single.enrichments.tsv
cat out/sc_brain_norm/fgwas.nearest.single*gt_90.params | grep "gt_90" | tr ' ' '\t' >> out/expr.sc_brain.norm.nearest_gt_90.single.enrichments.tsv

################################
# Do the same, but including microglia in each model (to condition on it)
# Can also do this with fgwas -w primary_microglia_gt_80+$ANN, which is less
# likely to fail getting the confidence intervals (I don't know why)
mkdir out/bulk_gtex_cond
mkdir out/bulk_brain_cond
mkdir out/sc_brain_cond

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.cond_mic.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.all.gt_80.tsv.gz -cc -w primary_microglia_gt_80 -cond $ANN -o out/bulk_gtex_cond/fgwas.nearest.cond_mic.$ANN"
done < <(cat expr.bulk.gtex.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.cond_mic.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_80.tsv.gz -cc -w Microglia_gt_80 -cond $ANN -o out/sc_brain_cond/fgwas.nearest.cond_mic.$ANN"
done < <(cat expr.sc.brain.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.cond_mic.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.all.gt_90.tsv.gz -cc -w primary_microglia_gt_90 -cond $ANN -o out/bulk_gtex_cond/fgwas.nearest.cond_mic.$ANN"
done < <(cat expr.bulk.gtex.relative.gt_90.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.nearest.cond_mic.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_90.tsv.gz -cc -w Microglia_gt_90 -cond $ANN -o out/sc_brain_cond/fgwas.nearest.cond_mic.$ANN"
done < <(cat expr.sc.brain.relative.gt_90.colnames.tsv | tr '\t' '\n')

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.nearest_gt_80.cond_mic.enrichments.tsv
cat out/bulk_gtex_cond/fgwas.nearest.cond_mic*gt_80.params | grep "gt_80" | grep -v "primary_microglia" | tr ' ' '\t' >> out/expr.bulk.gtex.nearest_gt_80.cond_mic.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.nearest_gt_80.cond_mic.enrichments.tsv
cat out/sc_brain_cond/fgwas.nearest.cond_mic*gt_80.params | grep "gt_80" | grep -v "Microglia" | tr ' ' '\t' >> out/expr.sc_brain.nearest_gt_80.cond_mic.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.nearest_gt_90.cond_mic.enrichments.tsv
cat out/bulk_gtex_cond/fgwas.nearest.cond_mic*gt_90.params | grep "gt_90" | grep -v "primary_microglia" | tr ' ' '\t' >> out/expr.bulk.gtex.nearest_gt_90.cond_mic.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.nearest_gt_90.cond_mic.enrichments.tsv
cat out/sc_brain_cond/fgwas.nearest.cond_mic*gt_90.params | grep "gt_90" | grep -v "Microglia" | tr ' ' '\t' >> out/expr.sc_brain.nearest_gt_90.cond_mic.enrichments.tsv

### When run with fgwas -w primary_microglia_gt_80+$ANN
echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.microglia_plus_ann.nearest_gt_80.single.enrichments.tsv
cat out/bulk_gtex_cond/fgwas.nearest.single.mic*gt_80.params | grep "gt_80" | grep -v "primary_microglia" | tr ' ' '\t' >> out/expr.bulk.gtex.microglia_plus_ann.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.microglia_plus_ann.nearest_gt_80.single.enrichments.tsv
cat out/sc_brain_cond/fgwas.nearest.single.mic*gt_80.params | grep "gt_80" | grep -v "Microglia" | tr ' ' '\t' >> out/expr.sc_brain.microglia_plus_ann.nearest_gt_80.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.bulk.gtex.microglia_plus_ann.nearest_gt_90.single.enrichments.tsv
cat out/bulk_gtex_cond/fgwas.nearest.single.mic*gt_90.params | grep "gt_90" | grep -v "primary_microglia" | tr ' ' '\t' >> out/expr.bulk.gtex.microglia_plus_ann.nearest_gt_90.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > out/expr.sc_brain.microglia_plus_ann.nearest_gt_90.single.enrichments.tsv
cat out/sc_brain_cond/fgwas.nearest.single.mic*gt_90.params | grep "gt_90" | grep -v "Microglia" | tr ' ' '\t' >> out/expr.sc_brain.microglia_plus_ann.nearest_gt_90.single.enrichments.tsv





while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.overlapping.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.all.gt_80.tsv.gz -cc -w $ANN -o out/bulk_gtex_cond/fgwas.overlapping.single.$ANN"
done < <(cat expr.bulk.all.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.overlapping.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.bulk.brain.gt_80.tsv.gz -cc -w $ANN -o out/bulk_brain_cond/fgwas.overlapping.single.$ANN"
done < <(cat expr.bulk.brain.relative.gt_80.colnames.tsv | tr '\t' '\n')

while read ANN; do
  echo $ANN
  submitJobs.py --MEM 5000 -j fgwas.overlapping.single.$ANN -q normal --nostdin \
      -c "fgwas -i input/AD.IGAP1_GWAX_exclude_firsts_v5.fgwas.ann.sc.brain.gt_80.tsv.gz -cc -w $ANN -o out/sc_brain_cond/fgwas.overlapping.single.$ANN"
done < <(cat expr.sc.brain.relative.gt_80.colnames.tsv | tr '\t' '\n')

echo -e "parameter\tCI_lo\testimate\tCI_hi" > bulk.gtex.overlapping.single.enrichments.tsv
cat out/bulk_gtex_cond/fgwas.overlapping.single*.params | grep "gt_80" | tr ' ' '\t' >> bulk.gtex.overlapping.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > bulk.brain.overlapping.single.enrichments.tsv
cat out/bulk_brain_cond/fgwas.overlapping.single*.params | grep "gt_80" | tr ' ' '\t' >> bulk.brain.overlapping.single.enrichments.tsv

echo -e "parameter\tCI_lo\testimate\tCI_hi" > sc_brain.overlapping.single.enrichments.tsv
cat out/sc_brain_cond/fgwas.overlapping.single*.params | grep "gt_80" | tr ' ' '\t' >> sc_brain.overlapping.single.enrichments.tsv



mathys_S6_celltype_subpopulation_signature_genes.tsv