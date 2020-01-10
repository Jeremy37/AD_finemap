#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SRC=$JS/AD_finemap/src

FNAME_QTL=$1
OUT=$2

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests should be in a particular format:
# geneid chr pos snp_id p.value beta Bonferroni.p.value FDR MAF std.error_of_beta

(echo -e "chr\tpos\trsid\tMAF";
 zcat $FNAME_QTL | sed '1d' | awk -F"\t" 'BEGIN{OFS="\t"}{af=$13; if (af > 0.5){af = 1-af}; print $1,$2,$3,af}' | sort -k1,1 -k2,2n | uniq) \
 | gzip > $OUT.variant_info.txt.gz

# They haven't computed FDR, nor corrected p values for number of tests.
# We get the min P value per gene, then do Bonferroni correction on nTests
(echo -e "feature\tchr\tpos\trsid\tpvalue\tbeta\texpressionIncreasingAllele\tgeneBiotype\tgeneStartPosition\tgeneEndPosition\tpBonf"; \
 zcat $FNAME_QTL | sed '1d' | sort -k5,5 -k2,2n | getLeadSnps.pl -f stdin --genecol 5 --pcol 8 \
  | awk -F"\t" 'BEGIN{OFS="\t"}{ pCor = $8*$20; if (pCor > 1) pCor = 1; print $5,$1,$2,$3,$8,$10,$14,$16,$17,$18,pCor }') \
  | gzip > $OUT.qtl_signals.txt.tmp.gz

# Finally, compute FDR using the Bonferroni-corrected p values per gene.
Rscript $SRC/utils/addFDRCol.R $OUT.qtl_signals.txt.tmp.gz 11 T | gzip > $OUT.qtl_signals.txt.gz
#rm $OUT.qtl_signals.txt.tmp.gz

(gzhead 1 $OUT.qtl_signals.txt.gz; zcat $OUT.qtl_signals.txt.gz | awk '$12 <= .05 && $4 ~ /^rs/') > $OUT.qtl_signals.FDR_0.05.txt
cut -f 1-4 $OUT.qtl_signals.FDR_0.05.txt > $OUT.qtl_signals.FDR_0.05.forColoc.txt


zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $5,$1,$2,$3,$8 }' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.nominals.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.nominals.txt.gz

