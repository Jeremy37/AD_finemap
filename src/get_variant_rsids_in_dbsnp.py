#!/usr/bin/env python
# Use dbSNP file to match SNPs by chr:pos_A1_A2 to get proper rsIDs for as many
# SNPs as we can. Write a table matching alleles in chr:pos_A1_A2 to rsID.
import argparse
import gzip
import os
import sys

parser = argparse.ArgumentParser(description="For variants with an ID in the form chrom:pos_A1_A2, identifies matching chrom, pos, ref, alt variants in dbSNP vcf")
parser.add_argument("--locus_file", type=str, required=True, metavar='FILE', help="File listing loci to prepare for fine-mapping") # From get_gcta_output_loci.py
parser.add_argument("--gwas_file", type=str, required=True, metavar='FILE', help="GWAS files with Z scores")
parser.add_argument("--window", type=int, default=5e5, metavar='INT', help="Half of the window size around lead SNP")
parser.add_argument("--dbsnp", type=str, required=True, metavar='FILE', help="Path to dbSNP VCF file")
parser.add_argument("--outdir", type=str, required=True, metavar='DIR', help="Directory for files")
args = parser.parse_args()

# Open a single file that will have a table matching SNP IDs for all loci
snp_file = open(os.path.join(args.outdir, "dbsnp_matching_table.tsv"), 'w')
print >>snp_file, "\t".join(["gwas_snpid", "chrom", "pos_str", "a1", "a2", "rsid", "dbsnp_chrpos", "ref", "alt", "matched_by"])
missing_file = open(os.path.join(args.outdir, "dbsnp_missing_table.tsv"), 'w')
print >>missing_file, "rsid\tmeta_pval\treason"

for i in open(args.locus_file, 'r'):
    if i.split()[0] == "Chr":
        continue
    
    [chrom,lead_snp,lead_pos,snps,poss,start,stop,n,lead_freq,lead_beta,lead_p,freqs,betas,ps,locus_gene] = i.split()
    start = int(start) - args.window
    stop = int(stop) + args.window
    locus = chrom + "_" + lead_pos
    sys.stderr.write(locus + "\n")
    locus_fname = os.path.join(args.outdir, locus + ".IGAP1_GWAX")
    
    gwas_locus_tmp_file = locus_fname + ".gwas.bgz"
    os.system("tabix {} {}:{}-{} | bgzip > {}".format(args.gwas_file, chrom, start, stop, gwas_locus_tmp_file))
    dbsnp_locus_tmp_file = locus_fname + ".dbsnp.bgz"
    os.system("tabix {} {}:{}-{} | bgzip > {}".format(args.dbsnp, chrom, start, stop, dbsnp_locus_tmp_file))
    
    # Read the dbSNP locus file into a dictionary
    dbsnp_chrpos_dict = {}
    dbsnp_rsid_dict = {}
    with gzip.open(dbsnp_locus_tmp_file) as f:
        for line in f:
            vals = line.split()
            snp_chrpos = "{}:{}_{}_{}".format(vals[0], vals[1], vals[3], vals[4])
            dbsnp_chrpos_dict[snp_chrpos] = [ vals[2], vals[3], vals[4] ]
            dbsnp_rsid_dict[vals[2]] = [ snp_chrpos, vals[3], vals[4] ]
    os.remove(dbsnp_locus_tmp_file)
    
    num_snps_written = 0
    with gzip.open(gwas_locus_tmp_file, 'r') as f:
        for line in f:
            vals = line.split()
            pos_str = vals[1]
            pos = int(vals[1])
            if pos >= start and pos <= stop:
                a1 = vals[3]
                a2 = vals[4]
                gwas_snpid = vals[2]
                snp_chrpos = "{}:{}_{}_{}".format(chrom, pos_str, a1, a2)
                p = float(vals[13])
                dbsnp_chrpos = ""
                rsid = ""
                matched_by = ""
                ref = ""
                alt = ""
                # If the snpID is already an rsID, then just write it and continue
                if gwas_snpid[0:2] == "rs":
                    rsid = gwas_snpid
                    if rsid in dbsnp_rsid_dict:
                        vals = dbsnp_rsid_dict[rsid]
                        dbsnp_chrpos = vals[0]
                        ref = vals[1]
                        alt = vals[2]
                        matched_by = "rsID in dbSNP"
                    else:
                        matched_by = "not matched"
                        print >>missing_file, "{}\t{}\trsID not in dbSNP".format(rsid, p)
                        # Don't skip this SNP, since VEP may be able to find it
                        #continue
                elif gwas_snpid in dbsnp_chrpos_dict:
                    dbsnp_chrpos = gwas_snpid
                    vals = dbsnp_chrpos_dict[gwas_snpid]
                    rsid = vals[0]
                    ref = vals[1]
                    alt = vals[2]
                    matched_by = "chr:pos_a1_a2 in dbSNP"
                else:
                    # See if we can match alleles after adusting e.g. TA/A to T/-
                    found = False
                    if len(a1) == 1 and len(a2) > 1:
                        a1 = "-"
                        a2 = a2[1:]
                    elif len(a2) == 1 and len(a1) > 1:
                        a2 = "-"
                        a1 = a1[1:]
                    if a1 == "-" or a2 == "-":
                        chrpos_testid = "{}:{}_{}_{}".format(chrom, pos_str, a1, a2)
                        if chrpos_testid in dbsnp_chrpos_dict:
                            dbsnp_chrpos = chrpos_testid
                            vals = dbsnp_chrpos_dict[chrpos_testid]
                            rsid = vals[0]
                            ref = vals[1]
                            alt = vals[2]
                            found = True
                            matched_by = "Adjusted chr:pos_a1_a2 in dbSNP"
                    if not found:
                        print >>missing_file, "{}\t{}\tchr:pos_a1_a2 not in dbSNP".format(gwas_snpid, p)
                        continue
                print >>snp_file, "\t".join([gwas_snpid, chrom, pos_str, a1, a2, rsid, dbsnp_chrpos, ref, alt, matched_by])



