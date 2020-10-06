#!/usr/bin/env python
import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Read loci from input file, and run finemap for each locus")
parser.add_argument("--input_loci", type=argparse.FileType('rt'), required=True, metavar='FILE', help="Input file with list of independent associated loci (produced by earlier pipeline steps)")
parser.add_argument("--input_dir", type=str, default="finemap/input/", metavar='FILE', help="Input file with list of independent associated loci (produced by earlier pipeline steps)")
parser.add_argument("--output_dir", type=str, default="finemap/output/", metavar='FILE', help="Input file with list of independent associated loci (produced by earlier pipeline steps)")
parser.add_argument("--ncausal", type=int, default=0, metavar='FILE', help="Input file with list of independent associated loci (produced by earlier pipeline steps)")
parser.add_argument("--method", type=str, default="sss", metavar='STR', help="Either 'sss' or 'cond' to choose the FINEMAP subprogram to run (default sss)")
args = parser.parse_args()

if not os.path.exists(args.input_dir):
    os.makedirs(args.input_dir)
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

indir = args.input_dir
if indir[len(indir)-1] != "/":
	indir = indir + "/"

outdir = args.output_dir
if outdir[len(outdir)-1] != "/":
	outdir = outdir + "/"

n_samples = str(21982+41944)

for i in args.input_loci:
    if i.split()[0] == "Chr":
        continue
    chrom = i.split()[0]
    locus = chrom + "_" + i.split()[2]
    # Skip APOE
    if locus == "19_45365447":
        continue
    nSignalsGCTA = i.split()[7]
    ncausal = nSignalsGCTA
    if args.ncausal > 0:
    	ncausal = args.ncausal
    ncausal = str(ncausal)
    write_in = open(indir + locus + ".Kunkle." + ncausal + ".in", 'a')
    print("z;ld;snp;config;cred;log;n_samples", file=write_in)
    print(indir + locus + ".Kunkle.z;" + indir + locus + ".Kunkle.ld;" + outdir + locus + ".Kunkle." + ncausal + ".snp;" + outdir + locus + ".Kunkle." + ncausal + ".config;" + outdir + locus + ".Kunkle." + ncausal + ".cred;" + outdir + locus + ".Kunkle." + ncausal + ".log;" + n_samples, file=write_in)
    write_in.close()
    
    cmd = "finemap --{} --in-files {}{}.Kunkle.{}.in --corr-config 0.9 --log --n-causal-snps {}".format(args.method, indir, locus, ncausal, ncausal)
    sys.stderr.write(cmd + "\n")
    os.system(cmd)
    #break
