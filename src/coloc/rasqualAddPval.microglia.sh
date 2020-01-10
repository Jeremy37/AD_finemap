#!/bin/bash
SRC=$1
RASQUAL_IN=$2
OUT=$3

# The awk line below removes multi-allelic SNPs. Rasqual tests each combination
# of SNP alleles (e.g. A,C vs G, and A,G vs C, and A vs G,C) and it's not clear
# which test to use in colocalisastion.
python $SRC/rasqualAddPval.py --rasqualfile $RASQUAL_IN \
    | awk '$5 !~ /,/ && $6 !~ /,/' \
    | perl -ane 'print join("\t", $F[0], @F[2..3], $F[1], $F[25])."\n";' \
    | sort -k2,2 -k3,3n | uniq \
    | bgzip > $OUT

tabix -s 2 -b 3 -e 3 $OUT
