#!/usr/bin/env python
import argparse
import os
import gzip
import sys

parser = argparse.ArgumentParser(description="For each locus with multiple signals, prepare input files and run GCTA --cojo-cond")
parser.add_argument("--locus_file", type=str, required=True, metavar='FILE', help="File listing locus details from previous GCTA runs") # From get_gcta_output_loci.py
parser.add_argument("--indir", type=str, default="gcta/input", metavar='DIR', help="Directory for gcta input files")
parser.add_argument("--window", type=int, default=5e5, metavar='INT', help="Half of the window size around lead SNP")
parser.add_argument("--outdir", type=str, default="gcta/output/cond", metavar='DIR', help="Directory for gcta output files")
args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
if not os.path.exists(os.path.join(args.indir, "cond")):
    os.makedirs(os.path.join(args.indir, "cond"))

# While doing conditional analyses, add to a file listing all independent associations,
# to be used later for coloc with conditioned p values.
cond_sig_file = open(os.path.join(args.outdir, "AD.meta.cond.indep_signals.tsv"),'w')
print >>cond_sig_file, "chr\tpos\trsid\tlocus_name\tp"

for i in open(args.locus_file, 'r'):
    if i.split()[0] == "Chr":
        continue
    
    [chrom,lead_snp,lead_pos,snps,poss,start,stop,nSignals,lead_freq,lead_beta,lead_p,freqs,betas,ps,locus_name] = i.split()
    if int(nSignals) <= 1:
        continue
        
    start = int(start) - args.window
    stop = int(stop) + args.window
    locus_chrpos = chrom + "_" + lead_pos
    sys.stderr.write("\n" + locus_name + "\n")
    
    bfile_fname = os.path.join(args.indir, "ukbb_sample.{}.merged".format(chrom))
    cojofile_fname = os.path.join(args.indir, "AD.IGAP1_GWAX.{}.ma".format(chrom))
    
    snpsList = snps.split(",")
    posList = poss.split(",")
    pList = ps.split(",")
    # For each independent SNP association, run GCTA cojo-cond
    for snpIndex in range(len(snpsList)):
        effect_snpid = snpsList[snpIndex]
        snplist_fname = os.path.join(args.indir, "cond", "{}.cond.snplist.{}".format(locus_name, effect_snpid))
        output_fname = os.path.join(args.outdir, "{}.cond.out.{}".format(locus_name, effect_snpid))
        
        # Write a snplist file with all other independent SNPs except this one
        snplist_file = open(snplist_fname, 'w')
        with open(snplist_fname, 'w') as snplist_file:
            for j in range(len(snpsList)):
                if j != snpIndex:
                    print >>snplist_file, snpsList[j]
        cmd = "gcta64 --bfile {} --chr {} --cojo-file {} --cojo-cond {} --out {}".format(bfile_fname, chrom, cojofile_fname, snplist_fname, output_fname)
        sys.stderr.write(cmd)
        os.system(cmd)
        
        # Filter the output file to only include the region of interest
        filtered_fname = os.path.join(args.outdir, "{}.cond.out.{}.flt.tsv".format(locus_name, effect_snpid))
        cmd = "(head -n 1 {}.cma.cojo; cat {}.cma.cojo | awk '$3 > {} && $3 < {} && $13 !~ /NA|pC/' | sort -k13,13g) > {}".format(output_fname, output_fname, start, stop, filtered_fname)
        sys.stderr.write(cmd)
        os.system(cmd)
        
        # Write to conditional independent signals file, for later use with coloc with conditioned p values.
        cond_locusname = locus_name + "_" + effect_snpid
        print >>cond_sig_file, "{}\t{}\t{}\t{}\t{}".format(chrom, posList[snpIndex], effect_snpid, cond_locusname, pList[snpIndex])
        
cond_sig_file.close()
