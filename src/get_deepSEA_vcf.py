#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="For each locus with multiple signals, prepare input files and run GCTA --cojo-cond")
parser.add_argument("--dbSNPfile", type=file, required=True, metavar='FILE', help="dbsnp_matching_table.tsv output by get_variant_rsids_in_dbsnp.py")
args = parser.parse_args()

for lineStr in args.dbSNPfile:
	[gwas_snpid,chrom,pos_str,a1,a2,rsid,dbsnp_chrpos,ref,alt,matched_by] = lineStr.split("\t")
	matched = False
	if ref and alt:
		#If the variant was in dbSNP, we get the ref/alt alleles from there
		# We need to check that we don't give multiple alt alleles, which will
		# be present in dbsNP.
		alts = alt.split(",")
		if len(alts) == 1:
			matched = True
		elif len(alts) > 1:
			for a in alts:
				if a == a2 or a == a1:
					alt = a
					matched = True
					break
	if matched:
		print "\t".join([chrom, pos_str, rsid, ref, alt])
	else:
		print "\t".join([chrom, pos_str, rsid, a1, a2])
		print "\t".join([chrom, pos_str, rsid, a2, a1])
