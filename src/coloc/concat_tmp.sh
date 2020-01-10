#!/bin/bash

DIR=$1
dat=$2

(while read -r -a myfile; do \
     cat $DIR/$myfile; \
 done < <(ls $DIR) ) | gzip > /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/AD_finemap/coloc/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.gene_sorted2.txt.gz
