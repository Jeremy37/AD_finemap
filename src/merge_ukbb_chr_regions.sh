#!/bin/bash

chrs=( 1 2 4 6 7 8 10 11 14 15 16 17 19 20 21 )
for chr in "${chrs[@]}"; do
  # Write a file list of the regions on this chr
  for filepath in gcta/input/ukbb_sample_regions/ukbb_sample.$chr.*.bed; do
    filepath_no_ext="${filepath%.*}"
    echo $filepath_no_ext >> gcta/input/ukbb_sample.$chr.file_list.txt
  done
  plink --merge-list gcta/input/ukbb_sample.$chr.file_list.txt --make-bed --out gcta/input/ukbb_sample.$chr.merged
  
  #rm gcta/input/ukbb_sample.$chr.file_list.txt
  echo -e "===========================================================\n\n\n"
done

