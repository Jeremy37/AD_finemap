#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SRC=$JS/AD_finemap/src

FNAME_QTL=$1
TISSUE_EGENES=$2
OUT=$3

# GTEx v8 summary stats have the format: (and are in GRCh38 coords)
#,phenotype_id,variant_id,tss_distance,maf,ma_samples,ma_count,pval_nominal,slope,slope_se
#0,ENSG00000015171.19,chr10_11501_C_A_b38,-122964,0.09210526,21,21,0.7301513913589728,-0.04181598,0.12085175105293121


# Get the file of nominal p values per gene
# Do this only for eGenes to keep the files reasonably small
zcat $FNAME_QTL | sed '1d' \
  | grep -f $TISSUE_EGENES \
  | perl -F/,/ -ane '@S=split(/_/, $F[2]); $chr=$S[0]; $chr =~ s/chr//; $id=join("_", @S[0..3]); print join("\t", $F[1], $chr, $S[1], $id, $F[7])."\n";' \
  | awk '!x[$0]++' | sort -k2,2 -k3,3n \
  | bgzip > $OUT.eGenes.FDR_0.05.nominals.txt.gz


# For colocalisation analyses we need a file of lead SNPs per eGene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests.
(echo -e "chr\tpos\trsid\tref\talt\tMAF"; \
 zcat $FNAME_QTL | sed '1d' \
  | grep -f $TISSUE_EGENES \
  | perl -F/,/ -ane '@S=split(/_/, $F[2]); $chr=$S[0]; $chr =~ s/chr//; $id=join("_", @S[0..3]); print join("\t", $chr, $S[1], $id, $S[2], $S[3], $F[4])."\n";' \
  | awk '!x[$0]++' | sort -k1,1 -k2,2n ) \
  | gzip > $OUT.eGenes.FDR_0.05.variants.vcf.gz

