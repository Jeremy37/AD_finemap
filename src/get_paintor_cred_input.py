#!/usr/bin/env python
import argparse
import os
import gzip
import gc
import sys
from scipy import stats

parser = argparse.ArgumentParser(description="Produce input files needed for PAINTOR - subset to 'probable' SNPs")
parser.add_argument("--locus_file", type=str, required=True, metavar='FILE', help="File listing loci to prepare for fine-mapping") # From get_gcta_output_loci.py
parser.add_argument("--gctadir", type=str, required=True, metavar='DIR', help="Directory for plink .bim files")
parser.add_argument("--outdir", type=str, required=True, metavar='DIR', help="Directory for files")
args = parser.parse_args()

# This script is used to get LD for the set of SNPs that we're using for PAINTOR
# fine-mapping
for i in open(args.locus_file, 'r'):
    if i.split()[0] == "Chr":
        continue
    
    [chrom,lead_snp,lead_pos,snps,poss,start,stop,n,lead_freq,lead_beta,lead_p,freqs,betas,ps,dummy] = i.split()
    locus = chrom + "_" + lead_pos
    sys.stderr.write(locus + "\n")
    locus_fname = os.path.join(args.outdir, locus + ".IGAP1_GWAX")
    
    # snplist file for LD file
    snplist_fname = locus_fname + ".snplist"
    os.system("sed '1d' {} | cut -d ' ' -f 1 > {}".format(locus_fname, snplist_fname))
    ## IMPORTANT NOTE: It is essential that the order of SNPs in the LD file extracted here
    # matches that in the paintor locus file. If not, then PAINTOR will run
    # based on incorrect LD data. NOTE that when data are extracted with plink, it doesn't
    # output the SNPs in the order they are found in the --extract file, but in the order they
    # are in the dataset itself (.bim file). So there is a great potential for error here.
    os.system("plink --bfile {}/ukbb_sample.{}.merged --extract {} --r square --out {}".format(args.gctadir, chrom, snplist_fname, locus_fname))
    # PAINTOR wants space-delimiters
    os.system("sed -i 's/\t/ /g' " + locus_fname + ".ld")
    
    # Check that there are no differences between the input snplist and the order
    # of SNPs extracted by plink
    os.system("plink --bfile {}/ukbb_sample.{}.merged --extract {} --recode vcf --out {}".format(args.gctadir, chrom, snplist_fname, locus_fname))
    os.system("cat {}.vcf | grep -v '^#' | cut -f 3 > {}.vcf.ids".format(locus_fname, locus_fname))
    print "Diff:"
    os.system("diff {}.snplist {}.vcf.ids".format(locus_fname, locus_fname))


os.system("rm {0}/*.nosex {0}/*.vcf".format(args.outdir))

# Also need to replace "nan" values in .ld files (SNPs that are monomorphic in 1KG)
# sed won't do because diagonals are also nan
# os.system("sed -i 's/nan/0/g' finemap_in/" + locus + ".IGAP1.ld")
# method below is cumbersome, but maybe ok for now
for i in open(args.locus_file, 'r'):
    if i.split()[0] == "Chr":
        continue
    locus = i.split()[0] + "_" + i.split()[2]
    #sys.stderr.write(locus + "\n")
    locus_fname = os.path.join(args.outdir, locus + ".IGAP1_GWAX")
    write_ld = open("temp.ld",'wa')
    ld = [i.split() for i in open(locus_fname + ".ld", 'r')]
    for i_ in range(len(ld)):
        for j_ in range(len(ld)):
            rsq = ld[i_][j_]
            if i_ == j_:
                rsq = "1.0"
            else:
                if rsq == "nan":
                    rsq = "0"
            print >>write_ld, rsq,
        print >>write_ld, ""
    write_ld.close()
    os.system("mv temp.ld " + locus_fname + ".ld")
    #gc.collect()
