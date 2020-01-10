#!/bin/bash
SRC=$JS/AD_finemap/src

FNAME_QTL=$1
FNAME_MAF=$2
OUT=$3

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests should be in a particular format:
# geneid chr pos snp_id p.value beta Bonferroni.p.value FDR MAF std.error_of_beta

zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $2,$3,$1 }' \
  | sort -k1,1 -k2,2n | uniq \
  | bgzip > $OUT.variant_info.partial.txt.gz

Rscript $SRC/coloc/xQTL_add_MAF.R $OUT.variant_info.partial.txt.gz $FNAME_MAF $OUT.variant_info.txt.gz

# xQTL datasets haven't computed FDR, nor corrected p values for number of tests.
# We get the min P value per gene, then do Bonferroni correction on nTests
(echo -e "feature\tchr\tpos\trsid\tp_value\tfeatureChromosome\tfeaturePositionStart\tSpearmanRho\tp_bonf"; \
 zcat $FNAME_QTL | sed '1d' | getLeadSnps.pl -f stdin --genecol 4 --pcol 8 \
  | awk -F"\t" 'BEGIN{OFS="\t"}{ pCor = $8*$10; if (pCor > 1) pCor = 1; print $4,$2,$3,$1,$8,$5,$6,$7,pCor }') \
  | gzip > $OUT.qtl_signals.txt.tmp.gz

# Finally, compute FDR using the Bonferroni-corrected p values per gene.
Rscript $SRC/addFDRCol.R $OUT.qtl_signals.txt.tmp.gz 9 T | gzip > $OUT.qtl_signals.txt.gz
rm $OUT.qtl_signals.txt.tmp.gz

(gzhead 1 $OUT.qtl_signals.txt.gz; zcat $OUT.qtl_signals.txt.gz | awk '$10 <= .05') > $OUT.qtl_signals.FDR_0.05.txt


zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $4,$2,$3,$1,$8 }' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.nominals.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.nominals.txt.gz

