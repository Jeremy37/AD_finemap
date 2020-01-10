#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="Read GCTA output files and summarise lead hits into one file")
parser.add_argument("--gcta_dir", type=str, required=True, metavar='DIR', help="Directory with GCTA output")
parser.add_argument("--pthresh", type=str, default="1e-5", metavar='FILE', help="P value threshold for independent SNPs")
parser.add_argument("--window", type=int, default=500000, metavar='INT', help="Window in bp for calling two SNPs as being in the same locus")
parser.add_argument("--out", type=str, required=True, metavar='FILE', help="Name of output file")
args = parser.parse_args()

pthreshStr = args.pthresh
pthresh = float(args.pthresh)


write_hits = open(args.out,'w')
print >>write_hits, "Chr	lead_SNP	lead_pos	SNPs	pos	start	stop	n_snps	lead_freq	lead_beta	lead_p	freq	beta	p"
for chrom in range(1,23):
    chrom = str(chrom)
    gcta_filepath = os.path.join(args.gcta_dir, "AD.IGAP1_GWAX." + chrom + ".jma.cojo")
    if not os.path.exists(gcta_filepath):
    	gcta_filepath = os.path.join(args.gcta_dir, "AD.IGAP1_GWAX." + pthreshStr + "." + chrom + ".jma.cojo")
    	if not os.path.exists(gcta_filepath):
        	continue
    # 1. Collapse variants in the same locus
    cojo = [i.split() for i in open(gcta_filepath, 'r') if i.split()[0] != "Chr" and float(i.split()[7]) < pthresh]
    indep_pos = [[cojo[0]]]
    for i in range(1,len(cojo)):
        pos = cojo[i][2]
        add_index = "NA"
        for j in range(len(indep_pos)):
            for k in indep_pos[j]:
                if abs(int(pos) - int(k[2])) < args.window:
                    add_index = j
        if add_index == "NA":
            indep_pos.append([cojo[i]])
        else:
            indep_pos[add_index].append(cojo[i])

    # 2. Print out list of hits
    for i in indep_pos:
        # region with min p value is lead snp
        [lead_snp,lead_pos,lead_beta,lead_p,lead_freq] = ["NA","NA","NA","NA","NA"]
        min_p = 1.0
        start = 3e9
        stop = 0
        for j in i:
            p = float(j[7])
            if p < min_p:
                [lead_snp,lead_pos,lead_beta,lead_p,lead_freq] = [j[1],j[2],j[5],j[7],j[4]]
                min_p = p
            pos = int(j[2])
            if pos < start:
                start = pos
            if pos > stop:
                stop = pos
        snps = ",".join([j[1] for j in i])
        poss = ",".join([j[2] for j in i])
        ps = ",".join([j[7] for j in i])
        betas = ",".join([j[5] for j in i])
        freqs = ",".join([j[4] for j in i])
        print >>write_hits, "\t".join([chrom,lead_snp,lead_pos,snps,poss,str(start),str(stop),str(len(i)),lead_freq,lead_beta,lead_p,freqs,betas,ps])
write_hits.close()
