#!/usr/bin/env python
import gzip
import os
import subprocess
import sys

if not os.path.exists("gcta/input/"):
    os.mkdir("gcta/input")

MIN_FREQ = 0.001
# Get allele frequencies from the UKBB proxy GWAX
freq_dict = {i.split()[0]:[i.split()[3],i.split()[4],i.split()[5]] for i in gzip.open("summary_stats/AD.proxy_exclude_firsts.af.tsv.gz", 'rt')}

# Format of AD.IGAP1_GWAX_exclude_firsts_v5.meta.chr{}.bgz:
#CHR	BP	SNP	A1	A2	GWAS_BETA	GWAS_SE	GWAS_P	GWAX_UKBB_BETA	GWAX_UKBB_SE	GWAX_UKBB_P	META_BETA	META_SE	META_P	DIRECT	I2	HET_P	FREQ	INFO
#1	662622	rs61769339	A	G	-0.1000	0.0457	0.02869	0.0403900106458	0.0226624987665	8.9E-02	0.012680341449	0.0203031660791	0.5322664744291312	-++	0.867977473714	0.005920095008868728	0.110178	0.777266
#1	693625	rs190214723	T	C	-0.0163	0.0685	0.812	-0.00180805029748	0.0439661246245	9.1E-01	-0.00603629002041	0.0370004463748	0.8704074146487125	--+	0	0.8586891674636307	0.950775	0.438968

# Numbers are: Kunkle cases + control + (GB proxy cases + GB cases)/4 + GB controls/4
#n = str(21982+41944+(54939+898)/4+(355900)/4)
# Kunkle only: number of cases + controls
n = str(21982+41944)
for do_chrom in range(1,23)[::-1]:
    do_chrom = str(do_chrom)
    print(do_chrom)
    #freq_dict = {i.split()[1]:[i.split()[2],i.split()[3],i.split()[4]] for i in gzip.open("summary_stats/ukbb_frequencies/ukbb.british." + do_chrom + ".frq.gz", 'rt')}
    
    # Read a dictionary of SNP IDs from the Plink file for each chromosome
    plink_file = "gcta/input/ukbb_sample.{}.merged.bim".format(do_chrom)
    if not os.path.isfile(plink_file):
        continue
    ukbb_id_dict = {i.split()[1]:i.split() for i in open(plink_file)}
    
    chr_file = "summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.chr{}.bgz".format(do_chrom)
    os.system("tabix summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz " + "{}:1-500000000".format(do_chrom) + " | bgzip > {}".format(chr_file))
    
    # Main file used by GCTA
    write_ma = open("gcta/input/AD.Kunkle." + str(do_chrom) + ".ma",'wt')
    print("SNP A1 A2 freq b se p N", file=write_ma)
    # A file where we save SNPs that are excluded
    write_excl = open("gcta/input/AD.Kunkle.excluded." + str(do_chrom) + ".txt",'wt')
    print("chr SNP A1 A2 freq b se p N exclusion", file=write_excl)
    for line in gzip.open(chr_file,'rt'):
        vals = line.split()
        if vals[0] != do_chrom:
            sys.exit("Unexpected chromosome!")
            continue
        snpid = vals[2]
        a1 = vals[3]
        a2 = vals[4]
        meta_beta = vals[11]
        meta_se = vals[12]
        meta_p = vals[13]
        kunkle_beta = vals[5]
        kunkle_se = vals[6]
        kunkle_p = vals[7]
        gwax_beta = vals[8]
        het_p = vals[16]
        freq_meta = vals[17]
        info = vals[18]
        
        if kunkle_p == "nan":
            print(do_chrom + " " + snpid + " " + a1 + " " + a2 + " " + freq_meta + " " + kunkle_beta + " " + kunkle_se + " " + kunkle_p + " " + n + " p_is_nan", file=write_excl)
            continue
        
        if float(info) < 0.85:
            print(do_chrom + " " + snpid + " " + a1 + " " + a2 + " " + freq_meta + " " + kunkle_beta + " " + kunkle_se + " " + kunkle_p + " " + n + " low_INFO_lt_0.85", file=write_excl)
            continue
        
        if not snpid in freq_dict:
            print(do_chrom + " " + snpid + " " + a1 + " " + a2 + " " + freq_meta + " " + kunkle_beta + " " + kunkle_se + " " + kunkle_p + " " + n + " freq_not_found_in_ukbb", file=write_excl)
            continue
        
        if float(het_p) < 0.001:
            print(do_chrom + " " + snpid + " " + a1 + " " + a2 + " " + freq_meta + " " + kunkle_beta + " " + kunkle_se + " " + kunkle_p + " " + n + " low_het_p", file=write_excl)
            continue
        
        newSnpID = "{}_{}_{}".format(snpid, a1, a2)
        if not newSnpID in ukbb_id_dict:
            # Try the reverse order of alleles
            newSnpID = "{}_{}_{}".format(snpid, a2, a1)
            if not newSnpID in ukbb_id_dict:
                print(newSnpID + "\tnot in ukbb_id_dict", file=write_excl)
                continue
        
        bim_a1 = ukbb_id_dict[newSnpID][4]
        if a2 == bim_a1:
            kunkle_beta = str(-float(kunkle_beta))
            tmp = a1
            a1 = a2
            a2 = tmp
        
        # Get the effect allele frequency from the dictionary
        # If the alleles are flipped with respect to the v3 meta-analysis, then
        # invert the allele frequency.
        freq = freq_dict[snpid][2]
        if a1 == freq_dict[snpid][1] or a2 == freq_dict[snpid][0]:
            freq = str(1-float(freq))
        
        if (float(freq) < MIN_FREQ or (1-float(freq)) < MIN_FREQ):
            print(newSnpID + "\tfreq ({0}) below MIN_FREQ ({1})".format(freq, MIN_FREQ), file=write_excl)
            continue
        
        print(newSnpID + " " + a1 + " " + a2 + " " + freq + " " + kunkle_beta + " " + kunkle_se + " " + kunkle_p + " " + n, file=write_ma)
    #os.remove(chr_file)
    write_ma.close()
