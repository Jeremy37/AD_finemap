#!/bin/bash
QTLDIR_IN=$1
DATASET=$2
OUT_ROOT=$3

SRC=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/AD_finemap/src

echo "QTLDIR_IN = $QTLDIR_IN"
echo "DATASET = $DATASET"
echo "OUT_ROOT = $OUT_ROOT"

QTL_BASE=$OUT_ROOT/$DATASET
FNAME_QTL_IN=$QTLDIR_IN/$DATASET.nominal.sorted.txt.gz

# First sort by gene, to get an input file for getting lead SNPs per gene
# We need an efficient way to do this, since sorting the whole file can be VERY slow
# We do it by chromosome, then sort those files by gene. We also only select SNPs
# within 500 kb of the gene (rather than 1 Mb), which speeds this up and gives us
# more QTLs at FDR 5%.
zcat $FNAME_QTL_IN | awk '$7 < 500000 && $7 > -500000' | awk -v FN=$FNAME_QTL_IN '{ print > FN"."$2".tmp.txt" }'
for chrfile in $FNAME_QTL_IN.*.tmp.txt; do
  echo $chrfile
  # Sort file in place by gene ID
  sort -k 1,1 -o $chrfile $chrfile
done

FNAME_QTL=$QTLDIR_IN/$DATASET.nominal.500k.gene_sorted.txt.gz
# Merge chr files together
cat $FNAME_QTL_IN.*.tmp.txt | gzip > $FNAME_QTL
rm $FNAME_QTL_IN.*.tmp.txt

# Get lead SNP per gene for each chr file, and correct for the number of tests per gene
(echo -e "feature\tfeature_chr\tfeature_start\tfeature_end\tstrand\tntested\tdist\trsid\tchr\tpos\tsnp_pos_end\tp_value\tslope\tis_top_variant\tntest\tp_bonf"; \
  zcat $FNAME_QTL | getLeadSnps.pl -f stdin --genecol 1 --pcol 12 \
   | awk -F"\t" 'BEGIN{OFS="\t"}{ pCor = $12*$16; if (pCor > 1) pCor = 1; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,pCor }') \
   > $QTL_BASE.qtl_signals.tsv.tmp

# Next, add FDR column to the "signals" file
Rscript $SRC/addFDRCol.R $QTL_BASE.qtl_signals.tsv.tmp 16 T > $QTL_BASE.qtl_signals.tsv
rm $QTL_BASE.qtl_signals.tsv.tmp

(head -n 1 $QTL_BASE.qtl_signals.tsv; cat $QTL_BASE.qtl_signals.tsv | awk '$17 <= .05') > $QTL_BASE.qtl_signals.FDR_5.tsv
sed '1d' $QTL_BASE.qtl_signals.FDR_5.tsv | cut -f 1 > $QTL_BASE.qtl_egenes.FDR_5.tsv

# Subset the nominal p values file to those genes with FDR < 0.05
hashJoin.pl --hashFile $QTL_BASE.qtl_egenes.FDR_5.tsv --scanFile $FNAME_QTL --colHashFile 1 --colScanFile 1 \
  | awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$9,$10,$8,$12}' \
  | sort -k2,2 -k3,3n \
  | bgzip > $QTL_BASE.nominals.egenes_FDR_5.tsv.gz

tabix -s 2 -b 3 -e 3 $QTL_BASE.nominals.egenes_FDR_5.tsv.gz

