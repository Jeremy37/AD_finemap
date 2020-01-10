#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SRC=$JS/AD_finemap/src

FNAME_QTL=$1
OUT=$2

# GTEx v8 summary stats have the format: (and are in GRCh38 coords)
#,phenotype_id,variant_id,tss_distance,maf,ma_samples,ma_count,pval_nominal,slope,slope_se
#0,ENSG00000015171.19,chr10_11501_C_A_b38,-122964,0.09210526,21,21,0.7301513913589728,-0.04181598,0.12085175105293121


# GTEx haven't provided a set of "eGenes" at FDR 5%, so we determine this.
# First get the min P value per gene, then do Bonferroni correction on nTests
(echo -e "feature\tchr\tpos\trsid\tp_value\tp_bonf"; \
 zcat $FNAME_QTL | sed '1d' | awk 'BEGIN {FS=",";OFS="\t"}{if(-500000 < $4 && $4 < 500000) print $2,$3,$4,$5,$6,$7,$8,$9,$10}' \
  | getLeadSnps.pl -f stdin --genecol 1 --pcol 7 \
  | awk 'BEGIN{OFS="\t"}{ pCor = $7*$11; if (pCor > 1) pCor = 1; print $1,$2,$7,pCor }' \
  | perl -ane '@S=split(/_/, $F[1]); $chr=$S[0]; $chr =~ s/chr//; $id=join("_", @S[0..3]); print join("\t", $F[0], $chr, $S[1], $id, $F[2], $F[3])."\n";' ) \
  | gzip > $OUT.qtl_signals.txt.gz
