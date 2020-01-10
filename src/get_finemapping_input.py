#!/usr/bin/env python
import argparse
import os
import gzip
import gc
import sys
import math
from scipy import stats

parser = argparse.ArgumentParser(description="Produce input files needed for FINEMAP (and PAINTOR)")
parser.add_argument("--locus_file", type=str, required=True, metavar='FILE', help="File listing loci to prepare for fine-mapping") # From get_gcta_output_loci.py
parser.add_argument("--gwas_file", type=str, required=True, metavar='FILE', help="GWAS files with Z scores")
parser.add_argument("--window", type=int, default=5e5, metavar='INT', help="Half of the window size around lead SNP")
parser.add_argument("--min_freq", type=float, default=0.0, metavar='FLOAT', help="Minimum allele frequency for SNPs to include [default 0]")
parser.add_argument("--min_info", type=float, default=0.0, metavar='FLOAT', help="Minimum INFO score for SNPs to include [default 0]")
parser.add_argument("--hetp_threshold", type=float, default=0.0, metavar='FLOAT', help="Minimum HET_P score for SNPs to include [default 0]")
parser.add_argument("--outdir", type=str, required=True, metavar='DIR', help="Directory for files")
args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# We read in 3 kinds of files:
# - The GWAS meta-analysis
# - The UK Biobank plink files (*.bim/bed, etc)
# - A file with UK Biobank allele frequencies
# Note that the plink files were made with SNP IDs in the format rs1234_A1_A2,
# whereas the other two files just have rsID, and separately the alleles. We need
# to make sure that the alleles match up so that we get the right allele freqs.

for i in open(args.locus_file, 'r'):
    if i.split()[0] == "Chr":
        continue
    
    [chrom,lead_snp,lead_pos,snps,poss,start,stop,n,lead_freq,lead_beta,lead_p,freqs,betas,ps,dummy] = i.split()
    start = int(start) - args.window
    stop = int(stop) + args.window
    locus = chrom + "_" + lead_pos
    sys.stderr.write(locus + "\n")
    locus_fname = os.path.join(args.outdir, locus + ".IGAP1_GWAX")
    
    a1_dict = {j.split()[1]:[j.split()[4],j.split()[5]] for j in open("gcta/input/ukbb_sample." + chrom + ".merged.bim") if int(j.split()[3]) >= start and int(j.split()[3]) <= stop}
    sys.stderr.write("Num items in a1_dict: {}\n".format(len(a1_dict)))
    nFlipped = 0
    ukbb_freq_fname = "summary_stats/ukbb_frequencies/ukbb.british." + chrom + ".frq.gz"
    sys.stderr.write("Loading allele frequency dictionary: {}\n".format(ukbb_freq_fname))
    freq_dict = {}
    for j in gzip.open(ukbb_freq_fname):
        line = j.split()
        if len(line) < 2:
            continue
        # Create SNP ID as rsID_A1_A2
        snp = "{}_{}_{}".format(line[1], line[2], line[3])
        if snp in a1_dict:
            freq = float(line[-2])
            allele = line[-4]
            if snp in freq_dict:
                # If the SNP is already in our dict, there is a duplicate. Add / overwrite
                # this one only if it has a higher allele frequency
                if freq <= freq_dict[snp]:
                    continue
            freq_dict[snp] = freq
        else:
            # See if reversed SNP alleles are in a1_dict
            snp = "{}_{}_{}".format(line[1], line[3], line[2])
            if snp in a1_dict:
                #print "Reversed SNP alleles {} found in a1_dict".format(snp)
                freq = 1 - float(line[-2])
                allele = line[-4]
                if snp in freq_dict:
                    # If the SNP is already in our dict, there is a duplicate. Add / overwrite
                    # this one only if it has a higher allele frequency
                    if abs(0.5 - freq) >= abs(0.5 - freq_dict[snp]):
                        continue
                freq_dict[snp] = freq
                nFlipped = nFlipped + 1
            
    sys.stderr.write("{} SNPs loaded in UKBB allele frequency dictionary\n".format(len(freq_dict)))
    sys.stderr.write("{} SNP alleles flipped between .bim file and UKBB .frq.gz file\n".format(nFlipped))
    sys.stderr.write("Num items in freq_dict: {}\n".format(len(freq_dict)))
    
    # order matters
    finemap_z_file = open(locus_fname + ".z", 'w')
    print >>finemap_z_file, "rsid chromosome position allele1 allele2 maf beta se"
    
    paintor_locus_file = open(locus_fname + ".paintor", 'w')
    print >>paintor_locus_file, "rsid chromosome position allele1 allele2 maf beta se Z"
    
    # snplist file for LD file
    snplist_fname = locus_fname + ".snplist"
    write_snplist = open(snplist_fname, 'w')
    write_missing = open(locus_fname + ".excluded_snps", 'w')
    locus_tmp_file = locus_fname + ".gwas.bgz"
    os.system("tabix {} {}:{}-{} | bgzip > {}".format(args.gwas_file, chrom, start, stop, locus_tmp_file))
    num_snps_written = 0
    for j in gzip.open(locus_tmp_file, 'r'):
        vals = j.split()
        pos = int(vals[1])
        if pos >= start and pos <= stop:
            p = float(vals[13])
            beta = float(vals[11])
            se = vals[12]
            z = stats.norm.isf(p/2)
            a1 = vals[3]
            a2 = vals[4]
            snpid = vals[2]
            het_p = float(vals[16])
            info = vals[18]
            snp = "{}_{}_{}".format(snpid, a1, a2)
            
            
            if math.isnan(p):
                print >>write_missing, snp + "\tp is nan"
                continue
            if het_p < args.hetp_threshold:
                print >>write_missing, snp + "\thet_p value {} below threshold {}".format(het_p, args.hetp_threshold)
                continue
            if float(info) < args.min_info:
                print >>write_missing, snp + "\tSNP info ({}) is less than cutoff of {}".format(info, args.min_info)
                continue
            if not snp in a1_dict:
                # Try the reverse order of alleles
                snp = "{}_{}_{}".format(snpid, a2, a1)
                if not snp in a1_dict:
                    print >>write_missing, snp + "\tnot in a1_dict"
                    continue
            
            # Check the alleles in the a1 dict. The SNPs are named e.g. rs1234_A_C, but we can't
            # trust the allele order coded in the SNP name, because it may be different in this BIM
            # file from when it was determined in the full dataset. Have to be VERY careful about
            # this!!
            if not a1 in a1_dict[snp]:
                print >>write_missing, snp + "\tfails allele match"
                continue
            bim_a1 = a1_dict[snp][0]
            if a2 == bim_a1:
                # SNP is there with reversed alleles, so we flip beta, as well as the alleles
                beta = -beta
                tmp = a1
                a1 = a2
                a2 = tmp
            if not snp in freq_dict:
                print >>write_missing, snp + "\tnot in freq_dict"
                continue
            freq = freq_dict[snp]
            if beta < 0:
                z = z * -1.0
            if freq > 0.5:
                freq = 1 - freq
            if freq < args.min_freq:
                print >>write_missing, snp + "\tSNP has freq <= {}".format(args.min_freq)
                continue
            print >>finemap_z_file,     "{} {} {} {} {} {} {} {}".format(snp, chrom, pos, a1, a2, freq, beta, se)
            print >>paintor_locus_file, "{} {} {} {} {} {} {} {} {}".format(snp, chrom, pos, a1, a2, freq, beta, se, z)
            print >>write_snplist, snp
            num_snps_written = num_snps_written + 1
    #os.remove(locus_tmp_file)
    finemap_z_file.close()
    write_snplist.close()
    write_missing.close()
    #gc.collect()
    sys.stderr.write("{} SNPs written to .snplist and .z files".format(num_snps_written))
    ## IMPORTANT NOTE: It is essential that the order of SNPs in the LD file extracted here
    # matches that in the finemap Z file (which we've output in the order the SNPs were found
    # in the GWAS meta-analysis file - chr:pos sorted order). If not, then FINEMAP will run
    # based on incorrect LD data. NOTE that when data are extracted with plink, it doesn't
    # output the SNPs in the order they are found in the --extract file, but in the order they
    # are in the dataset itself (.bim file). So there is a great potential for error here.
    os.system("plink --bfile gcta/input/ukbb_sample.{}.merged --extract {} --r square --out {}".format(chrom, snplist_fname, locus_fname))
    # Finemap seems to prefer space-delimiters
    os.system("sed -i 's/\t/ /g' " + locus_fname + ".ld")

# Also need to replace "nan" values in .ld files (SNPs that are monomorphic in 1KG)
# sed won't do because diagonals are also nan
# os.system("sed -i 's/nan/0/g' finemap_in/" + locus + ".IGAP1.ld")
# method below is cumbersome, but maybe ok for now
for i in open(args.locus_file, 'r'):
    if i.split()[0] == "Chr":
        continue
    locus = i.split()[0] + "_" + i.split()[2]
    sys.stderr.write(locus + "\n")
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
