#!/usr/bin/env python
##
import argparse
import sys
import os
import os.path
import gzip
import datetime
import numpy as np
import subprocess

parser = argparse.ArgumentParser(description="Given a set of SNPs and a list of bed files, uses bedtools to get overlap of each SNP with all of the bed files, and outputs a table as well as a summary file.")
parser.add_argument("--snpbed", required=True, type=argparse.FileType('r'), metavar='FILE', help="BED file of SNP positions to intersect")
parser.add_argument("--bedfilelist", required=True, type=argparse.FileType('r'), metavar='FILE', help="File with list of bed file paths to intersect with the SNP file")
parser.add_argument("--output", required=True, type=str, metavar='FILE', help="Base path for output")
parser.add_argument("--norun", action='store_true', help="don't run external commands (e.g. bedtools)")
parser.add_argument("--debug", action='store_true', help="print additional debugging information")
parser.add_argument('--verbose', '-v', action='count')
args = parser.parse_args()


def main():
	bedfilelist = [l.strip() for l in args.bedfilelist.readlines()]
	if args.verbose > 1:
		sys.stderr.write(str(datetime.datetime.now())+"\n")
	if args.verbose:
		sys.stderr.write("Num bed files: " + str(len(bedfilelist)) + "\n")

	if len(bedfilelist) < 1:
		die("No lines in --bedfilelist file: {0}".format(args.bedfilelist.name))
		
	num_snps = sum(1 for line in args.snpbed)
	args.snpbed.seek(0)
	if args.verbose:
		sys.stderr.write("Num SNPs: " + str(num_snps) + "\n")
	
	intersect = np.zeros((num_snps, len(bedfilelist)), dtype=np.int)
	
	bedNames = []
	for i,linestr in enumerate(bedfilelist):
		#for linestr in args.bedfilelist:
		lineVals = linestr.strip().split('\t')
		fdesc = lineVals[0]
		fpath = lineVals[1]
		bedNames.append(fdesc)
		fname = os.path.basename(fpath)
		#cmd = "bedtools intersect -c -a {0} -b {1} > {2}.intersect.{3}".format(args.snpbed.name, fpath, args.output, fname)
		#output = docheckoutput(cmd).strip()

		cmd = "bedtools intersect -c -a {0} -b {1} | cut -f 5".format(args.snpbed.name, fpath)
		output = docheckoutput(cmd).strip()
		intersect[:,i] = np.int32(output.split('\n'))
		#print(overlaps)
	
	sums = np.sum(intersect, axis=1)
	bedOverlapStrs = [",".join(getBedNames(intersect[i,:], bedNames)) for i in range(num_snps)]
	outFileName = args.output + ".overlap.summary.txt"
	with open(outFileName, 'w') as f:
		f.write("\t".join(["chr", "start", "end", "overlap_count", "overlap_names"]) + "\n")
		i = 0
		for linestr in args.snpbed:
			lineVals = linestr.strip().split('\t')
			f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(lineVals[0], lineVals[1], lineVals[2], sums[i], bedOverlapStrs[i]))
			i = i + 1
	
	outFileName = args.output + ".overlaps.txt"
	with open(outFileName, 'w') as f:
		f.write("\t".join(["chr", "start", "end"] + bedNames) + "\n")
		i = 0
		args.snpbed.seek(0)
		for linestr in args.snpbed:
			lineVals = linestr.strip().split('\t')
			f.write("{0}\t{1}\t{2}\t{3}\n".format(lineVals[0], lineVals[1], lineVals[2], "\t".join([str(x) for x in intersect[i,:]])))
			i = i + 1
	
	#print(sums)
	#print(bedNames)
	#print(bedOverlapStrs)
	
	if args.verbose:
		if args.verbose > 1: sys.stderr.write(str(datetime.datetime.now())+"\n")
		sys.stderr.write("\nDone!\n")


def getBedNames(overlaps, bedNames):
	return [bedNames[i] for i,overlap in enumerate(overlaps) if overlap > 0]

def docheckoutput(cmd):
	if args.verbose:
		sys.stderr.write(cmd + "\n")
	retval = ""
	if (not args.norun):
		retval = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	return retval

def openg(f, mode=None):
	if type(f) is file:
		if isgzfile(f.name):
			return gzip.GzipFile(fileobj=f)
		return f
	elif type(f) is str:
		if isgzfile(f):
			if mode is None:
				mode = 'rb'
			return gzip.open(f, mode)
		if mode is None:
			mode = 'r'
		return open(f, mode)
	else:
		die("openg: unrecognized type for parameter 'f': {0}".format(type(f)))


main()
