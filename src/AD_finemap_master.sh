###############################################################################
# This script contains all essential processing steps for the AD fine-mapping
# paper following the GWAS-GWAX meta-analysis, beginning with summary stats from
# the meta-analysis, and bed/bim/fam from 10k UK Biobank individuals for the AD
# regions.
# We use GCTA to identify number of causal SNPs, annotate SNPs and loci, and do
# eQTL colocalisations. We run FINEMAP and PAINTOR
# This script is not meant to be run directly; rather, steps can be run on the
# command line one at a time, and this file is an overview of the steps.

ROOT=/path/to/AD_finemap
SRC=$ROOT/src
GWAS_NAME=AD.meta
#SWDIR=define here
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen2/jeremys
ROOT=$JS/AD_finemap
SRC=$ROOT/src
GWAS_NAME=AD.meta

cd $ROOT

###############################################################################
# Summary stats
# Download summary stat files
md() {
	mkdir -p "$1" && cd "$1"
}
md summary_stats
wget https://zenodo.org/record/3531493/files/AD.proxy_exclude_firsts.bgen.stats.gz
wget https://zenodo.org/record/3531493/files/AD.IGAP1_GWAX_exclude_firsts_v5.meta.gz
zcat AD.IGAP1_GWAX_exclude_firsts_v5.meta.gz | tr ' ' '\t' | sort -k1,1n -k2,2n | bgzip > AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz
tabix -S 1 -s 1 -b 2 -e 2 AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz
rm AD.IGAP1_GWAX_exclude_firsts_v5.meta.gz

# Make a summary stats file formatted for GWAS catalog. GWAS catalog validator
# only likes rsIDs, so we have to set our indel IDs to NA.
(echo -e "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error\tSNP_ID\tGWAS_BETA\tGWAS_SE\tGWAS_P\tGWAX_UKBB_BETA\tGWAX_UKBB_SE\tGWAX_UKBB_P\tDIRECT\tI2\tHET_P\tINFO";
 zcat AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{rsid=$3; pval=$14; if ($3 !~ /^rs/) { rsid="NA" } if ($14 !~ /nan/) { print rsid,$14,$1,$2,$4,$5,$18,$12,$13,$3,$6,$7,$8,$9,$10,$11,$15,$16,$17,$19} }') \
 | head -n 200 > AD.IGAP1_GWAX_exclude_firsts_v5.meta.gwas_catalog.tsv.head

(echo -e "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error\tSNP_ID\tGWAS_BETA\tGWAS_SE\tGWAS_P\tGWAX_UKBB_BETA\tGWAX_UKBB_SE\tGWAX_UKBB_P\tDIRECT\tI2\tHET_P\tINFO";
 zcat AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{rsid=$3; pval=$14; if ($3 !~ /^rs/) { rsid="NA" } if ($14 !~ /nan/) { print rsid,$14,$1,$2,$4,$5,$18,$12,$13,$3,$6,$7,$8,$9,$10,$11,$15,$16,$17,$19} }') \
 | gzip > AD.IGAP1_GWAX_exclude_firsts_v5.meta.gwas_catalog.tsv.gz
md5sum AD.IGAP1_GWAX_exclude_firsts_v5.meta.gwas_catalog.tsv.gz > AD.IGAP1_GWAX_exclude_firsts_v5.meta.gwas_catalog.tsv.gz.md5sum

#pip3 install ss-validate
ss-validate -f AD.IGAP1_GWAX_exclude_firsts_v5.meta.gwas_catalog.tsv.gz --logfile AD.IGAP1_GWAX_exclude_firsts_v5.meta.gwas_catalog.tsv.validate.log



###############################################################################
# Note that a few UK Biobank files were used but are not provided for download.
# These are:
# gcta/input/ukbb_sample_regions/ukbb_sample.<region>.bim/bed/fam - PLINK files with genotype counts for 10k randomly sampled UKB individuals used for LD
# summary_stats/ukbb_frequencies/ukbb.british.<chr>.frq.gz	- allele freqs for the 10k randomly sampled UKB individuals

###############################################################################
# Use GCTA-COJO to identify significant independent loci
# Not that this was already done genome-wide by Jimmy, which identified the
# regions I started with.
# gcta/input/ukbb_sample_regions needs to be a folder with plink files for UKBB
# sample data around independent AD signals

cd $ROOT
mkdir gcta/output_1e-5

# Merge together plink bed files on the same chromosome.
# NOTE: these UK biobank files are not provided here
submitJobs.py --MEM 7000 -j merge_ukbb_chr_regions -q yesterday \
  -c "bash $SRC/merge_ukbb_chr_regions.sh"

# Get just the allele frequencies from the UKBB proxy GWAS.
(zcat $ROOT/summary_stats/AD.proxy_exclude_firsts.bgen.stats.gz | cut -f 1,2,3,5,6,7 | head -n 1;
 zcat $ROOT/summary_stats/AD.proxy_exclude_firsts.bgen.stats.gz | cut -f 1,2,3,5,6,7 | sort -k6,6rg) | gzip > summary_stats/AD.proxy_exclude_firsts.af.tsv.gz

# Here we separate IGAP1+GWAX summary statistics by chromosome, match UK Biobank allele
# frequencies, and write GCTA input into gcta/input/AD.IGAP1_GWAX.[chr].ma
submitJobs.py --MEM 9000 -j get_gcta_input -q yesterday \
  -c "python $SRC/get_gcta_input.py"

# Run GCTA-COJO on each chromosome. Use  P<1×10−5  p-value cutoff, since we can
# later subset to only genome-wide significant loci.
for filepath in gcta/input/*.bed; do
  echo "$filepath"
  filename=$(basename "$filepath") # filename like ukbb_sample.1.156156033_166156033.bed
  filepath_no_ext="${filepath%.*}"
  chr=`echo $filename | perl -ne '@Ar=split(/\./); print $Ar[1];'`
  
  submitJobs.py --MEM 2000 -q normal -j gcta.1e-5.$chr \
    -c "gcta64 \
  --bfile $filepath_no_ext \
  --chr $chr \
  --cojo-file gcta/input/AD.IGAP1_GWAX.$chr.ma \
  --cojo-slct \
  --cojo-p 1e-5 \
  --out gcta/output_1e-5/AD.IGAP1_GWAX.$chr"
done

# For APOE and HLA, do a special run with a threshold of 5e-8 to avoid overfitting
# the strong signal based on imperfect LD match
mkdir gcta/output_5e-8
chr_paths=( gcta/input/ukbb_sample.6.merged.bed gcta/input/ukbb_sample.19.merged.bed )
for filepath in "${chr_paths[@]}"; do
  echo "$filepath"
  filename=$(basename "$filepath") # filename like ukbb_sample.1.156156033_166156033.bed
  filepath_no_ext="${filepath%.*}"
  chr=`echo $filename | perl -ne '@Ar=split(/\./); print $Ar[1];'`
  
  submitJobs.py --MEM 2000 -q normal -j gcta.5e-8.$chr \
    -c "gcta64 \
  --bfile $filepath_no_ext \
  --chr $chr \
  --cojo-file gcta/input/AD.IGAP1_GWAX.$chr.ma \
  --cojo-slct \
  --cojo-p 5e-8 \
  --out gcta/output_5e-8/AD.IGAP1_GWAX.$chr"
done

grep "Successfully" FarmOut/gcta*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code" FarmOut/gcta*.txt | wc -l

# This python script reads the output from GCTA (per chromosome gcta_out/AD.IGAP1_GWAX.[chr].jma.cojo files) and
# creates list of independent loci. If distance between two SNPs is < 500kb, they will be collapsed into a single
# locus. Remove SNPs with marginal P>5×10−8. Output file (in this example, "AD.IGAP1_GWAX.indep.hits") columns are
# "Chr lead_SNP lead_pos SNPs pos start stop n_snps lead_freq lead_beta lead_p freq beta p".
python $SRC/get_gcta_output_loci.py --gcta_dir gcta/output_1e-5 --pthresh 5e-8 --out AD.IGAP1_GWAX.gwsig_indep.hits
# Use a distance of 620 kb so that APP-ADAMTS1 is considered 1 locus, as variants
# are in partial LD
python $SRC/get_gcta_output_loci.py --gcta_dir gcta/output_1e-5 --pthresh 1e-5 --window 620000 --out AD.IGAP1_GWAX.1e-5_indep.hits
python $SRC/get_gcta_output_loci.py --gcta_dir gcta/output_5e-8 --pthresh 5e-8 --out AD.IGAP1_GWAX.APOE_HLA.hits

cp AD.IGAP1_GWAX.1e-5_indep.hits AD.loci.tsv
cp AD.IGAP1_GWAX.1e-5_indep.hits AD.loci.1e-5.tsv
## NOTE: I manually edited these files to add in locus names based on
# nearby genes, and to incorporate the update to the HLA, APOE loci.
# I removed the loci with lead p < 5e-8 from the AD.loci.tsv file. It's important
# that this starts with info from the "1e-5" file, except for these loci, since
# that will include secondary signals. 
# Because APOE p values are so small (reported as zero), we set the lead SNP to
# the one with the lowest p value (rs429358_T_C) in the independent studies
# (Kunkle and GWAX), which is the same SNP in both.

# Extract SNP IDs for all loci with p < 1e-5, then look these up in the meta-analysis
# with IGAP stage2 to see if any are genome-wide significant.
sed '1d' AD.IGAP1_GWAX.1e-5_indep.hits | cut -f 2 | cut -f 1 -d "_" > AD.IGAP1_GWAX.1e-5_indep.snp_ids.txt
(zcat summary_stats/AD.IGAP2_GWAX_exclude_firsts_v5.meta.tsv.bgz | head -n 1;
 zcat summary_stats/AD.IGAP2_GWAX_exclude_firsts_v5.meta.tsv.bgz | grep -wFf AD.IGAP1_GWAX.1e-5_indep.snp_ids.txt) > AD.IGAP2_GWAX_exclude_firsts_v5.meta.1e-5_snps.tsv
# The output supports the TSPOAP1/MIR142 locus and slightly weakens CD33. All other
# loci were genome-wide signficant already, and our borderline loci aren't included
# in IGAP2 so give no new information.


###############################################################################
# Do GCTA conditional analysis for regions with multiple signals

mkdir gcta/output_1e-5/cond

# Conditional analysis for all loci except APOE, where the signal is so strong
# that our LD ref panel size isn't sufficient for fine-mapping.
cat AD.loci.tsv | grep -v "APOE" > AD.loci.exceptAPOE.tsv

submitJobs.py --MEM 2000 -q yesterday -j gcta_cojo_cond \
    -c "python $SRC/run_gcta_conditional.py --locus_file AD.loci.exceptAPOE.tsv --outdir gcta/output_1e-5/cond"

# Append together the output files from conditional analyses
paste <(echo -e "locus\tcond_snp") <(head -n 1 gcta/output_1e-5/cond/ACE.cond.out.rs4311_T_C.flt.tsv) > gcta/output_1e-5/cond/merged_loci.cond.out.flt.tsv
for filepath in gcta/output_1e-5/cond/*.cond.out.*.flt.tsv; do
  echo "$filepath"
  filename=$(basename "$filepath") # filename like ACE.cond.out.rs3730025_A_G.flt.tsv
  filepath_no_ext="${filepath%.*}"
  locus=`echo $filename | perl -ne '@Ar=split(/\./); print $Ar[0];'`
  cond_snp=`echo $filename | perl -ne '@Ar=split(/\./); print $Ar[3];'`
  sed '1d' $filepath | perl -sne 'print $locus."\t".$cond."\t".$_' -- -locus=$locus -cond=$cond_snp >> gcta/output_1e-5/cond/merged_loci.cond.out.flt.tsv
done


###############################################################################
# FINEMAP

mkdir $ROOT/finemap
mkdir $ROOT/finemap/input

(head -n 1 AD.loci.tsv; cat AD.loci.tsv | grep "HLA") > AD.loci.HLA.tsv
cat AD.loci.tsv | grep -v "HLA" > AD.loci.exceptHLA.tsv

(head -n 1 AD.loci.tsv; cat AD.loci.tsv | grep "TREM2") > AD.loci.TREM2.tsv

# Prepare input files for running FINEMAP (and PAINTOR). This takes a couple of hours.
# We only include SNPs with MAF > 0.002, since this is at the limit of our accuracy
# for LD estimation (20 individuals in 10,000).
submitJobs.py --MEM 10000 -n 5 -j get_finemapping_input.exceptHLA -q yesterday \
  -c "python $SRC/get_finemapping_input.py --locus_file AD.loci.exceptHLA.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.002 --min_info 0.85 --hetp_threshold 0.001 --outdir finemap/input"
submitJobs.py --MEM 19000 -n 5 -j get_finemapping_input.HLA -q yesterday \
  -c "python $SRC/get_finemapping_input.py --locus_file AD.loci.HLA.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.002 --min_info 0.85 --hetp_threshold 0.001 --outdir finemap/input"
submitJobs.py --MEM 10000 -n 5 -j get_finemapping_input.TREM2 -q yesterday \
  -c "python $SRC/get_finemapping_input.py --locus_file AD.loci.TREM2.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.0019 --min_info 0.85 --hetp_threshold 0.001 --outdir finemap/input"

# Running this mainly to get each of the meta-analysis locus summary stat files
cat AD.IGAP1_GWAX.1e-5_indep.hits | awk '$11 > 4e-7' > AD.loci.not_gwsig.tsv
submitJobs.py --MEM 28000 -n 5 -j get_finemapping_input.1e-5 -q yesterday \
  -c "python $SRC/get_finemapping_input.py --locus_file AD.loci.not_gwsig.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.002 --min_info 0.855 --hetp_threshold 0.001 --outdir finemap/input_1e-5"

# Determine how many SNPs were excluded by different criteria. First paste together
# all of the *.excluded_snps files and *.snplist files that were output.
cat finemap/input/*.excluded_snps > finemap/input/all_excluded_snps.txt
cat finemap/input/*.snplist > finemap/input/all.snplist.txt
Rscript $SRC/count_snp_exclusions.R

# The above step also extracts the meta-analysis SNPs at each associated locus,
# and now we merge these into a single file for later use.
(zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | head -n 1;
 zcat finemap/input/*.IGAP1_GWAX.gwas.bgz | sort -k 1,1n -k2,2n) | bgzip > summary_stats/$GWAS_NAME.assoc_loci.gz
tabix -s 1 -b 2 -e 2 -S 1 summary_stats/$GWAS_NAME.assoc_loci.gz

# Run FINEMAP statistical fine-mapping
submitJobs.py --MEM 2000 -j run_finemap_ncausal_1 -q yesterday \
  -c "python $SRC/run_finemap.py --input_loci AD.loci.tsv --input_dir finemap/input/ --output_dir finemap/out_ncausal_1  --ncausal 1"

submitJobs.py --MEM 2500 -j run_finemap_ncausal_2 -q yesterday \
  -c "python $SRC/run_finemap.py --input_loci AD.loci.tsv --input_dir finemap/input/ --output_dir finemap/out_ncausal_2  --ncausal 2"

submitJobs.py --MEM 2000 -j run_finemap_ncausal_3 -q yesterday \
  -c "python $SRC/run_finemap.py --input_loci AD.loci.TREM2.tsv --input_dir finemap/input/ --output_dir finemap/out_ncausal_3  --ncausal 3"

# Merge finemap locus results into one file
mkdir finemap/output
(head -n 1 finemap/out_ncausal_1/1_161155392.IGAP1_GWAX.1.snp; cat finemap/out_ncausal_1/*.IGAP1_GWAX.1.snp | grep -v "^index") | tr ' ' '\t' > finemap/output/$GWAS_NAME.finemap.ncausal_1.snp
(head -n 1 finemap/out_ncausal_2/1_161155392.IGAP1_GWAX.2.snp; cat finemap/out_ncausal_2/*.IGAP1_GWAX.2.snp | grep -v "^index") | tr ' ' '\t' > finemap/output/$GWAS_NAME.finemap.ncausal_2.snp
(head -n 1 finemap/out_ncausal_3/6_40942196.IGAP1_GWAX.3.snp; cat finemap/out_ncausal_3/*.IGAP1_GWAX.3.snp | grep -v "^index") | tr ' ' '\t' > finemap/output/$GWAS_NAME.finemap.ncausal_3.snp


###############################################################################
# Prepare files for LocusZoom plotting of the associated regions

echo -e "snp\tchr\tstart\tend\tflank\trun\tm2zargs" > AD.loci.1e-5.forLocusZoom.tsv
sed '1d' AD.loci.1e-5.tsv | perl -ane '($snp,$b,$c) = split(/_/, $F[1]); print sprintf("%s\t%s\t%d\t%d\tNA\tyes\ttitle=\"%s\"\n", $snp, $F[0], $F[5]-400000, $F[6]+400000, $F[1]);' >> AD.loci.1e-5.forLocusZoom.tsv

(echo -e "MarkerName\tP.value";
 zcat summary_stats/$GWAS_NAME.assoc_loci.gz | sed '1d' \
 | perl -ane '$marker=$F[2]; if ($marker !~ /^rs/) { $marker=$F[0].":".$F[1];} print $marker."\t".$F[13]."\n";') \
 | gzip > summary_stats/$GWAS_NAME.assoc_loci.forLocusZoom.gz

# Also make a file which has GCTA conditional p values for each independent association
echo -e "snp\tchr\tstart\tend\tflank\trun\tm2zargs" > AD.loci.conditional_1.forLocusZoom.tsv
sed '1d' AD.loci.tsv | awk '$8 == 2' | perl -ane '($snp,$b,$c) = split(/_/, $F[1]); print sprintf("%s\t%s\t%d\t%d\tNA\tyes\ttitle=\"%s\"\n", $snp, $F[0], $F[5]-400000, $F[6]+400000, $F[1]);' >> AD.loci.conditional_1.forLocusZoom.tsv
cp AD.loci.conditional_1.forLocusZoom.tsv AD.loci.conditional_2.forLocusZoom.tsv
# Edit this file to specify the secondary SNP at each association

(echo -e "MarkerName\tP.value";
 (sed '1d' gcta/output_1e-5/cond/NCK2.cond.out.rs143080277_T_C.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/BIN1.cond.out.rs6733839_C_T.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/EPHA1.cond.out.rs12703526_G_T.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/PTK2B-CLU.cond.out.rs867230_C_A.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/ADAM10.cond.out.rs442495_T_C.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/ACE.cond.out.rs4311_T_C.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/ABCA7.cond.out.rs12151021_A_G.flt.tsv;
  sed '1d' gcta/output_1e-5/cond/APP-ADAMTS1.cond.out.rs4817090_T_C.flt.tsv) | sort -k1,1n -k 3,3n | python $SRC/get_gcta_cond_p_for_locuszoom.py;) \
 > gcta/output_1e-5/cond/$GWAS_NAME.conditional_1.forLocusZoom.tsv
 
 gzip -c gcta/output_1e-5/cond/$GWAS_NAME.conditional_1.forLocusZoom.tsv > gcta/output_1e-5/cond/$GWAS_NAME.conditional_1.forLocusZoom.tsv.gz
 
 (echo -e "MarkerName\tP.value";
  (sed '1d' gcta/output_1e-5/cond/NCK2.cond.out.rs116038905_T_C.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/BIN1.cond.out.rs7584040_C_T.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/EPHA1.cond.out.rs10265814_C_G.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/PTK2B-CLU.cond.out.rs73223431_C_T.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/ADAM10.cond.out.rs4775044_T_G.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/ACE.cond.out.rs3730025_A_G.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/ABCA7.cond.out.rs4147918_A_G.flt.tsv;
   sed '1d' gcta/output_1e-5/cond/APP-ADAMTS1.cond.out.rs2830489_C_T.flt.tsv) | sort -k1,1n -k 3,3n | python $SRC/get_gcta_cond_p_for_locuszoom.py;) \
 > gcta/output_1e-5/cond/$GWAS_NAME.conditional_2.forLocusZoom.tsv

 gzip -c gcta/output_1e-5/cond/$GWAS_NAME.conditional_2.forLocusZoom.tsv > gcta/output_1e-5/cond/$GWAS_NAME.conditional_2.forLocusZoom.tsv.gz


################################################################################
# Identify genes in regions of association
mkdir genes
mkdir reference

# Get latest Gencode annotation (in GRCh37 coords)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz -O $ROOT/reference/gencode.v29lift37.annotation.gtf.gz
zcat $ROOT/reference/gencode.v29lift37.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /gene/ and $F[8] =~ /gene_id "(ENSG[^"]+).*gene_name "([^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1,$2)."\n"; }' \
    > $ROOT/reference/gencode.v29lift37.genes.bed
zcat $ROOT/reference/gencode.v29lift37.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /gene/ and $F[8] =~ /gene_type "protein_coding"/ and $F[8] =~ /gene_id "(ENSG[^"]+).*gene_name "([^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1,$2)."\n"; }' \
    > $ROOT/reference/gencode.v29lift37.genes.pc.bed

# Get the closest gene to each lead SNP, as well as all genes within 100 kb or 500 kb.
# Gencode annotations are 1-based. BED format is half-open, so start coordinate
# is the start of the feature, and end coordinate is not included in the feature.
sed '1d' AD.loci.tsv | awk 'BEGIN{OFS="\t"}{print "chr"$1,$3,$3+1,$2,$15}' > genes/AD.loci.leadSNP.chr.bed
sed '1d' AD.loci.tsv | awk 'BEGIN{OFS="\t"}{print "chr"$1,$6-100000,$7+100000,$2,$15}' > genes/AD.loci.200kb_window.chr.bed
sed '1d' AD.loci.tsv | awk 'BEGIN{OFS="\t"}{print "chr"$1,$6-500000,$7+500000,$2,$15}' > genes/AD.loci.1Mb_window.chr.bed

# Get closest gene, and annotate lead SNP file with nearest gene (protein coding, or any)
echo -e "chr\tpos\tsnp\tlocus\tstrand\tgene_id\tsymbol\tdist" > genes/AD.loci.nearest_gene.pc.tsv
bedtools closest -a genes/AD.loci.leadSNP.chr.bed -b $ROOT/reference/gencode.v29lift37.genes.pc.bed -d -t first -g $ROOT/reference/GRCh37.genome.chr.txt \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5,$9,$10,$11,$12}' >> genes/AD.loci.nearest_gene.pc.tsv

echo -e "chr\tpos\tsnp\tlocus\tstrand\tgene_id\tsymbol\tdist" > genes/AD.loci.nearest_gene.all.tsv
bedtools closest -a genes/AD.loci.leadSNP.chr.bed -b $ROOT/reference/gencode.v29lift37.genes.bed -d -t first -g $ROOT/reference/GRCh37.genome.chr.txt \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5,$9,$10,$11,$12}' >> genes/AD.loci.nearest_gene.all.tsv
    
    
# Get genes within a window of lead SNPs
echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > genes/AD.loci.200kb_window.gene_overlaps.pc.tsv
bedtools intersect -a genes/AD.loci.200kb_window.chr.bed -b $ROOT/reference/gencode.v29lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' >> genes/AD.loci.200kb_window.gene_overlaps.pc.tsv

echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > genes/AD.loci.1Mb_window.gene_overlaps.pc.tsv
bedtools intersect -a genes/AD.loci.1Mb_window.chr.bed -b $ROOT/reference/gencode.v29lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' >> genes/AD.loci.1Mb_window.gene_overlaps.pc.tsv

echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > genes/AD.loci.200kb_window.gene_overlaps.all.tsv
bedtools intersect -a genes/AD.loci.200kb_window.chr.bed -b $ROOT/reference/gencode.v29lift37.genes.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' >> genes/AD.loci.200kb_window.gene_overlaps.all.tsv

echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > genes/AD.loci.1Mb_window.gene_overlaps.all.tsv
bedtools intersect -a genes/AD.loci.1Mb_window.chr.bed -b $ROOT/reference/gencode.v29lift37.genes.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' >> genes/AD.loci.1Mb_window.gene_overlaps.all.tsv

Rscript $SRC/annotateRegionGeneDistance.R --loci AD.loci.tsv --genes genes/AD.loci.1Mb_window.gene_overlaps.pc.tsv > genes/AD.loci.1Mb_window.genedist.pc.tsv
(head -n 1 genes/AD.loci.1Mb_window.genedist.pc.tsv; sed '1d' genes/AD.loci.1Mb_window.genedist.pc.tsv | awk '$12 <= 100000 && $12 >= -100000') > genes/AD.loci.200kb_window.genedist.pc.tsv

Rscript $SRC/annotateRegionGeneDistance.R --loci AD.loci.tsv --genes genes/AD.loci.1Mb_window.gene_overlaps.all.tsv > genes/AD.loci.1Mb_window.genedist.all.tsv
(head -n 1 genes/AD.loci.1Mb_window.genedist.all.tsv; sed '1d' genes/AD.loci.1Mb_window.genedist.all.tsv | awk '$12 <= 100000 && $12 >= -100000') > genes/AD.loci.200kb_window.genedist.all.tsv


Rscript $SRC/aggregateColumn.R --file genes/AD.loci.200kb_window.gene_overlaps.pc.tsv --groupby leadSNP --aggregate symbol > genes/AD.loci.200kb_window.gene_overlaps.aggregated.pc.tsv
Rscript $SRC/aggregateColumn.R --file genes/AD.loci.1Mb_window.gene_overlaps.pc.tsv --groupby leadSNP --aggregate symbol > genes/AD.loci.1Mb_window.gene_overlaps.aggregated.pc.tsv

Rscript $SRC/aggregateColumn.R --file genes/AD.loci.200kb_window.gene_overlaps.all.tsv --groupby leadSNP --aggregate symbol > genes/AD.loci.200kb_window.gene_overlaps.aggregated.all.tsv
Rscript $SRC/aggregateColumn.R --file genes/AD.loci.1Mb_window.gene_overlaps.all.tsv --groupby leadSNP --aggregate symbol > genes/AD.loci.1Mb_window.gene_overlaps.aggregated.all.tsv

# Add nearby gene names to lead SNPs file
Rscript $SRC/annotateRegionGenes.R AD.loci.tsv \
                                   genes/AD.loci.nearest_gene.pc.tsv \
                                   genes/AD.loci.200kb_window.gene_overlaps.aggregated.pc.tsv \
                                   genes/AD.loci.1Mb_window.gene_overlaps.aggregated.pc.tsv \
                                   > genes/AD.loci.gene_overlaps.pc.tsv

Rscript $SRC/annotateRegionGenes.R AD.loci.tsv \
                                   genes/AD.loci.nearest_gene.all.tsv \
                                   genes/AD.loci.200kb_window.gene_overlaps.aggregated.all.tsv \
                                   genes/AD.loci.1Mb_window.gene_overlaps.aggregated.all.tsv \
                                   > genes/AD.loci.gene_overlaps.all.tsv


################################################################################
# Annotate gene expression across tissues / cell types

mkdir expression

# Get a single tab-separated file with TPM expression of genes across multiple tissues
# NOTE that these datasets are not all provided
# Produces files reference/tissues.combined.tpm.tsv.gz and reference/tissues.selected.tpm.tsv.gz
Rscript $SRC/get_tissue_expr.R

# Prepare single-cell data from Allen Brain institute, based on paper
# Hodge, R.D., et al. (2019). "Conserved cell types with divergent features in human versus mouse cortex." Nature 573:61-68.
# First put ~4 Gb file "transcrip.tome" and sample_annotations.csv into expression directory.
cd expression
submitJobs.py --MEM 16000 -j get_brain_singlecell_test -q yesterday \
  -c "Rscript $SRC/get_allen_brain_cortex_data.R"

# Convert the expression TPM files to expression percentile for each gene relative
# to all other tissues / cell types.
Rscript $SRC/get_relative_expr.R


Rscript $SRC/annotateRegionGeneExpr.R genes/AD.loci.1Mb_window.genedist.pc.tsv \
                                      reference/tissues.selected.tpm.tsv.gz \
                                      > genes/AD.loci.1Mb_window.gene_overlaps.pc.tpms.tmp.tsv

Rscript $SRC/annotateRegionGeneExpr.R genes/AD.loci.1Mb_window.genedist.all.tsv \
                                      reference/tissues.selected.tpm.tsv.gz \
                                      > genes/AD.loci.1Mb_window.gene_overlaps.all.tpms.tmp.tsv

# Subset to GTEx tissues plus microglia
zcat reference/tissues.combined.tpm.tsv.gz | cut -f1,2,6,52-105 | gzip > reference/tissues.microglia_and_gtex.tpm.tsv.gz
Rscript $SRC/annotateRegionGeneExpr.R genes/AD.loci.1Mb_window.genedist.pc.tsv \
                                      reference/tissues.microglia_and_gtex.tpm.tsv.gz \
                                      > genes/AD.loci.1Mb_window.gene_overlaps.pc.microglia_and_gtex.tpms.tsv

Rscript $SRC/annotateRegionGeneExpr.R genes/AD.loci.1Mb_window.genedist.all.tsv \
                                      reference/tissues.microglia_and_gtex.tpm.tsv.gz \
                                      > genes/AD.loci.1Mb_window.gene_overlaps.all.microglia_and_gtex.tpms.tsv

# Put together a file which includes the percentile columns for microglia and brain
# relative to GTEx, but which includes only gene expression for selected datasets.
paste <(cut -f 1-15 genes/AD.loci.1Mb_window.gene_overlaps.all.microglia_and_gtex.tpms.tsv) <(cut -f 1-15 --complement genes/AD.loci.1Mb_window.gene_overlaps.all.tpms.tmp.tsv) > genes/AD.loci.1Mb_window.gene_overlaps.all.selected_tpms.tsv
paste <(cut -f 1-15 genes/AD.loci.1Mb_window.gene_overlaps.pc.microglia_and_gtex.tpms.tsv) <(cut -f 1-15 --complement genes/AD.loci.1Mb_window.gene_overlaps.pc.tpms.tmp.tsv) > genes/AD.loci.1Mb_window.gene_overlaps.pc.selected_tpms.tsv


###############################################################################
# Annotate variants
# Start with VEP. For annotation it's best if we get rsIDs for all indels where
# this is possible. Since some of our meta-analysis indels don't have rsIDs, but
# instead the form 17:5016253_CA_C. We look these up in dbSNP.

mkdir $ROOT/annotated

# Download dbSNP file
mkdir $ROOT/reference/dbSNP
cd $ROOT/reference/dbSNP
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz.tbi
DBSNP_FILE=$ROOT/reference/dbSNP/All_20170710.vcf.gz

# Download conservation tracks
cd $ROOT/reference/
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/hg19.100way.phastCons.bw
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw

cd $ROOT

# Try to match indel alleles in dbSNP
# Produces output file annotated/dbsnp_matching_table.tsv
submitJobs.py --MEM 1000 -j get_variant_rsids_in_dbsnp -q yesterday \
   -c "python $SRC/get_variant_rsids_in_dbsnp.py --locus_file AD.loci.tsv --gwas_file summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --dbsnp $DBSNP_FILE --outdir annotated"
rm annotated/*.gwas.bgz

sed '1d' annotated/dbsnp_matching_table.tsv | cut -f 6 > annotated/dbsnp_ids.for_vep.txt
# Pass the file  file to VEP - run manually in web browser. Add options:
# - HGVS; ESP allele freqs; gnomAD allele freqs; include flagged variants;
#   exon/intron numbers; dbNSFP score - MetaLR score, M-CAP score, CADD_Phred,
#   fathmm-MKL_coding_score, integrated_fitcons_score, integrated_confidence_value,
#   GERP++RS, phyloP100way_vert, phylop20way_mammal, phastCons100way_vert,
#   phastCons20way_mammal; and CADD

# Save as AD.meta.vep_output.tsv
# Rename the #Uploaded_variation column to variant
# NOTE: HAVE TO CHANGE OUTPUT FROM WINDOWS CR-LF line endings to Unix LF
gzip annotated/$GWAS_NAME.vep_output.tsv

# Use CrossMap to get GRCh38 coords for SNPs. Output to *.GRCh38.bed (works on farm3)
# May need to set path to include CrossMap
CHAIN_FILE=/path/to/GRCh37_to_GRCh38.chain.gz
zcat summary_stats/$GWAS_NAME.assoc_loci.gz | sed '1d' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3}' > $GWAS_NAME.assoc_loci.bed
CrossMap.py bed $CHAIN_FILE $GWAS_NAME.assoc_loci.bed annotated/$GWAS_NAME.assoc_loci.GRCh38.bed
rm $GWAS_NAME.assoc_loci.bed

# Overlap SNPs with genome annotations - Roadmap enhancer states, cell type
# ATAC profiles, enhancer-promoter links, etc.
mkdir overlaps
# Prepare lists of bed files in overlaps/ used by the script below (takes 1.5 hrs)
submitJobs.py --MEM 5000 -j do_genome_overlaps -q yesterday \
    -c "bash $SRC/do_genome_overlaps.sh summary_stats/$GWAS_NAME.assoc_loci.gz $GWAS_NAME"

# First download: http://yiplab.cse.cuhk.edu.hk/jeme/encoderoadmap_lasso.zip
# These data are annotated in the fine-mapping table, but aren't used in any other steps
# This download is available from the paper web site: http://yiplab.cse.cuhk.edu.hk/jeme/
# The file roadmap.description.csv is made from the table at that site.
bash $SRC/overlap.JEME.sh $ROOT/reference/JEME/encoderoadmap_lasso overlaps/$GWAS_NAME.overlap_input.chr.bed
# Takes ~12 hrs
submitJobs.py --MEM 10000 -j AD.annotate.JEME -q yesterday \
    -c "Rscript $SRC/merge_JEME_overlaps.R overlaps/JEME $ROOT/reference/JEME overlaps/$GWAS_NAME.overlap_input.bed $ROOT/annotated/$GWAS_NAME.overlaps.JEME.tsv"

# Annotate with SpliceAI (takes 20 min)
submitJobs.py --MEM 10000 -j AD.annotate.spliceai -q yesterday \
    -c "Rscript $SRC/annotate_spliceai.R AD.loci.tsv summary_stats/$GWAS_NAME.assoc_loci.gz $ROOT/reference/spliceai.merge.tsv.bgz $ROOT/annotated/$GWAS_NAME.spliceAI.tsv"

# Annotate with distance to nearest RNA CaptureSeq exon
(echo -e "chr\tpos\tend\tsnp\tfeatureChr\tfeatureStart\tfeatureEnd\tfeature\tfeatureDist";
 bedtools closest -d -t first -a <(sort -k1,1 -k2,2n overlaps/$GWAS_NAME.overlap_input.bed) -b <(sort -k1,1 -k2,2n reference/GSE118158_filtered_hybrid_transcriptome.exons.bed) ) | gzip > overlaps/$GWAS_NAME.RNACaptureSeq.closestExon.tsv.gz

(echo -e "chr\tpos\tend\tsnp\tfeatureChr\tfeatureStart\tfeatureEnd\tfeature\tfeatureDist";
bedtools closest -d -t first -a <(sort -k1,1 -k2,2n overlaps/$GWAS_NAME.overlap_input.bed) -b <(sort -k1,1 -k2,2n reference/GSE118158_filtered_hybrid_transcriptome.spliceBoundaries.bed) ) | gzip > overlaps/$GWAS_NAME.RNACaptureSeq.closestSplice.tsv.gz

gunzip -c overlaps/$GWAS_NAME.RNACaptureSeq.closestExon.tsv.gz > overlaps/$GWAS_NAME.RNACaptureSeq.closestExon.tsv
gunzip -c overlaps/$GWAS_NAME.RNACaptureSeq.closestSplice.tsv.gz > overlaps/$GWAS_NAME.RNACaptureSeq.closestSplice.tsv

# Get annotation values for phastCons, phyloP, GERP
bigWigAverageOverBed reference/hg19.100way.phastCons.bw overlaps/$GWAS_NAME.overlap_input.chr.bed annotated/$GWAS_NAME.phastCons100way.txt -bedOut=annotated/$GWAS_NAME.phastCons100way.bed
bigWigAverageOverBed reference/hg19.100way.phyloP100way.bw overlaps/$GWAS_NAME.overlap_input.chr.bed annotated/$GWAS_NAME.phyloP100way.txt -bedOut=annotated/$GWAS_NAME.phyloP100way.bed
bigWigAverageOverBed reference/All_hg19_RS.bw overlaps/$GWAS_NAME.overlap_input.chr.bed annotated/$GWAS_NAME.GERP_RS.txt -bedOut=annotated/$GWAS_NAME.GERP_RS.bed

# Get variant format for DeepSEA predictions. If the variant was in dbSNP, we get
# the ref/alt alleles from there; otherwise, since I can't be sure which is ref vs. alt,
# I try both configurations. One will be discarded by DeepSEA.
python $SRC/get_deepSEA_vcf.py --dbSNPfile annotated/dbsnp_matching_table.tsv > annotated/AD.meta.variants_for_deepSEA.vcf

# Split into 3 folds since DeepSEA limits the number of variants at once.
sed '1d' annotated/AD.meta.variants_for_deepSEA.vcf | head -n 57500 > annotated/AD.meta.variants_for_deepSEA.1.vcf
sed '1d' annotated/AD.meta.variants_for_deepSEA.vcf | tail -n +57501 | head -n 57500 > annotated/AD.meta.variants_for_deepSEA.2.vcf
sed '1d' annotated/AD.meta.variants_for_deepSEA.vcf | tail -n +115001 > annotated/AD.meta.variants_for_deepSEA.3.vcf


md $ROOT/annotated/DeepSEA
wget http://deepsea.princeton.edu/media/output/7101b819-fd00-4ae7-9308-70a61fd0942e/jobs.zip
mv jobs.zip AD_meta_v5.DeepSEA.1.zip
wget http://deepsea.princeton.edu/media/output/ad872980-e31c-4262-97b8-225b08552314/jobs.zip
mv jobs.zip AD_meta_v5.DeepSEA.2.zip
wget http://deepsea.princeton.edu/media/output/4fd3cb86-786c-4e79-9a5b-f30a01c60480/jobs.zip
mv jobs.zip AD_meta_v5.DeepSEA.3.zip
unzip AD_meta_v5.DeepSEA.1.zip infile.vcf.out.funsig
mv infile.vcf.out.funsig AD_meta_v5.DeepSEA.1.funsig
unzip AD_meta_v5.DeepSEA.2.zip infile.vcf.out.funsig
mv infile.vcf.out.funsig AD_meta_v5.DeepSEA.2.funsig
unzip AD_meta_v5.DeepSEA.3.zip infile.vcf.out.funsig
mv infile.vcf.out.funsig AD_meta_v5.DeepSEA.3.funsig
(cat AD_meta_v5.DeepSEA.1.funsig; sed '1d' AD_meta_v5.DeepSEA.2.funsig; sed '1d' AD_meta_v5.DeepSEA.3.funsig) > AD_meta_v5.DeepSEA.all.funsig
cd $ROOT

# Merge together all our various annotations with our GWAS summary stats into a single file
submitJobs.py --MEM 5000 -j AD.merge_annotations -q yesterday \
    -c "Rscript $SRC/merge_annotations.R AD.loci.exceptAPOE.tsv \
                                         summary_stats/$GWAS_NAME.assoc_loci.gz \
                                         annotated/$GWAS_NAME.assoc_loci.GRCh38.bed \
                                         annotated/dbsnp_matching_table.tsv \
                                         annotated/$GWAS_NAME.vep_output.tsv.gz \
                                         reference/vep_impact_table.tsv \
                                         annotated/$GWAS_NAME.overlaps.txt \
                                         annotated/$GWAS_NAME.overlaps.JEME.tsv \
                                         annotated/$GWAS_NAME.spliceAI.tsv \
                                         annotated/$GWAS_NAME.phastCons100way.bed \
                                         annotated/$GWAS_NAME.phyloP100way.bed \
                                         annotated/$GWAS_NAME.GERP_RS.bed \
                                         annotated/DeepSEA/AD_meta_v5.DeepSEA.all.funsig \
                                         overlaps/$GWAS_NAME.RNACaptureSeq.closestExon.tsv \
                                         overlaps/$GWAS_NAME.RNACaptureSeq.closestSplice.tsv \
                                         gcta/output_1e-5/cond/merged_loci.cond.out.flt.tsv \
                                         finemap/output/$GWAS_NAME.finemap.ncausal_1.snp \
                                         finemap/output/$GWAS_NAME.finemap.ncausal_2.snp \
                                         finemap/output/$GWAS_NAME.finemap.ncausal_3.snp \
                                         annotated/$GWAS_NAME.annotated"
# Produces annotated/AD.meta.annotated.selected.probable.tsv, the key fine-mapping file



###############################################################################
# PAINTOR: 
# We already got most Paintor input when running get_finemapping_input.py above.
# But they need to all be in the same paintor input directory.
md paintor/input
cd $ROOT
for f in $ROOT/finemap/input/*.ld; do
  echo "$f"
  filename=$(basename "$f")
  cp $f paintor/input/$filename
  echo $filename
done
cp finemap/input/*.paintor paintor/input/
for f in paintor/input/*.paintor; do
  filepath_no_ext="${f%.*}"
  echo "mv $f $filepath_no_ext"
  mv $f $filepath_no_ext
done

# For unknown reasons, when the CLNK locus is included then PAINTOR always gives
# an error, so we exclude it.
cat AD.loci.exceptAPOE.tsv | grep -v "CLNK" > AD.loci.exceptAPOE_CLNK.tsv
Rscript $SRC/get_paintor_input.R --locus_file AD.loci.exceptAPOE_CLNK.tsv > paintor/input/paintor.loci

submitJobs.py --MEM 5000 -j get_paintor_annot -q yesterday \
  -c "Rscript $SRC/get_paintor_annotations.R --locus_file AD.loci.exceptAPOE_CLNK.tsv --annotated annotated/$GWAS_NAME.annotated.all.tsv.gz --indir paintor/input/ --outdir paintor/input/"

# Check correlation between annotations. To do this, first make a single file with all SNPs x annotations.
head -n 1 paintor/input/1_161155392.IGAP1_GWAX.annotations > paintor/all_annotations.txt
for f in paintor/input/*.annotations; do
  sed '1d' $f >> paintor/all_annotations.txt
done
# Check correlations in R
Rscript $SRC/check_paintor_annotation_correlations.R --file paintor/all_annotations.txt --minFraction 0.001 --out paintor/paintor_annotation_correlations.pdf


mkdir paintor/out_single_annotations
ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site stop_gained frameshift inframe_insertion inframe_deletion start_lost coding_nonsyn UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.01 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Blood_Immune_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
for ann in "${ANNOTATIONS[@]}"; do
  if [ ! -d "paintor/out_single_annotations/$ann" ]; then
    mkdir -p paintor/out_single_annotations/$ann
  fi
  submitJobs.py --MEM 25000 -j run_paintor_nc1.$ann -q normal \
    -c "bash $SRC/run_paintor.sh paintor/input/ paintor/out_single_annotations/$ann 1 $ann"
done

cat paintor/out_single_annotations/*/*.Values > paintor/out_single_annotations/out_single_annotations.enrichments.txt

# Running PAINTOR with all SNPs across loci is prohibitively slow, and also does
# not account for multiple causal variants. Therefore, we run only on SNPs with
# probability from FINEMAP > 0.0001. This is essentially solving a slightly different
# problem: what are causal SNP enrichments *among those with prob > 0.0001?

mkdir paintor_cred
rsync -rptgoD paintor/input/* paintor_cred/input/
for f in paintor_cred/input/*.IGAP1_GWAX; do
  echo "mv $f $f.in"
  mv $f $f.in
done
submitJobs.py --MEM 500 -j get_paintor_annot.credible -q yesterday \
  -c "Rscript $SRC/get_paintor_annotations.credible.R --locus_file AD.loci.exceptAPOE_CLNK.tsv --annotated annotated/$GWAS_NAME.annotated.selected.probable.tsv --indir paintor_cred/input/ --outdir paintor_cred/input/"

submitJobs.py --MEM 1500 -n 5 -j get_paintor_cred_ld -q yesterday \
  -c "python $SRC/get_paintor_cred_input.py --locus_file AD.loci.exceptAPOE_CLNK.tsv --gctadir gcta/input --outdir paintor_cred/input"

# Check correlation between annotations. To do this, first make a single file
# with all SNPs x annotations.
head -n 1 paintor_cred/input/1_161155392.IGAP1_GWAX.annotations > paintor_cred/all_annotations.txt
for f in paintor_cred/input/*.annotations; do
  sed '1d' $f >> paintor_cred/all_annotations.txt
done

# Check correlations in R
Rscript $SRC/check_paintor_annotation_correlations.R --file paintor_cred/all_annotations.txt --minFraction 0.001 --out paintor_cred/paintor_annotation_correlations.pdf

# Run PAINTOR with each single annotation in the model
# First using number of causals = 1
mkdir paintor_cred/out_single_annotations_nc1
#ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site stop_gained frameshift inframe_insertion inframe_deletion start_lost coding_nonsyn UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.01 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia ipsMacrophage microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Blood_Immune_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
# Exclude the annotations that cover too small a fraction of SNPs
ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site coding_nonsyn UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.01 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia ipsMacrophage microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Blood_Immune_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
for ann in "${ANNOTATIONS[@]}"; do
  if [ ! -d "paintor_cred/out_single_annotations_nc1/$ann" ]; then
    mkdir -p paintor_cred/out_single_annotations_nc1/$ann
  fi
  submitJobs.py --MEM 500 -j run_paintor_cred_nc1.$ann -q normal \
    -c "bash $SRC/run_paintor.sh paintor_cred/input/ paintor_cred/out_single_annotations_nc1/$ann 1 $ann"
done
grep "Successfully" FarmOut/run_paintor_cred_nc1*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/run_paintor_cred_nc1*.txt | wc -l

# Run PAINTOR again using number of causals = 2
mkdir paintor_cred/out_single_annotations_nc2
ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site coding_nonsyn UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.01 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia ipsMacrophage microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Blood_Immune_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
for ann in "${ANNOTATIONS[@]}"; do
  if [ ! -d "paintor_cred/out_single_annotations_nc2/$ann" ]; then
    mkdir -p paintor_cred/out_single_annotations_nc2/$ann
  fi
  submitJobs.py --MEM 500 -j run_paintor_cred_nc2.$ann -q normal \
    -c "bash $SRC/run_paintor.sh paintor_cred/input/ paintor_cred/out_single_annotations_nc2/$ann 2 $ann"
done
grep "Successfully" FarmOut/run_paintor_cred_nc2*txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/run_paintor_cred_nc2*txt | wc -l

cat paintor_cred/out_single_annotations_nc1/*/*.Values > paintor_cred/out_single_annotations_nc1/out_single_annotations_nc1.enrichments.txt
cat paintor_cred/out_single_annotations_nc2/*/*.Values > paintor_cred/out_single_annotations_nc2/out_single_annotations_nc2.enrichments.txt
grep "" $ROOT/paintor_cred/out_single_annotations_nc1/*/Log.BayesFactor | sed 's/\.\///g' | sed 's/\/Log.BayesFactor:/\t/g' | sed 's/.*\///g' > $ROOT/paintor_cred/out_single_annotations_nc1/out_single_annotations_nc1.LogBayesFactor.txt
grep "" $ROOT/paintor_cred/out_single_annotations_nc2/*/Log.BayesFactor | sed 's/\.\///g' | sed 's/\/Log.BayesFactor:/\t/g' | sed 's/.*\///g' > $ROOT/paintor_cred/out_single_annotations_nc2/out_single_annotations_nc2.LogBayesFactor.txt


################# MODEL 2
# Select the Blood_Immune_DNase annotation, and run PAINTOR again with this
# annotation plus each other annotation individually.
cd $ROOT
mkdir paintor_cred/out_nc2/model2
ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site coding_nonsyn UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.01 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia ipsMacrophage microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
#ANNOTATIONS=( deepSEA_lt_0.05 )
for ann in "${ANNOTATIONS[@]}"; do
  ANNOTS="Blood_Immune_DNase,$ann"
  submitJobs.py --MEM 500 -j run_paintor_cred_model2.$ann -q normal \
    -c "bash $SRC/run_paintor.sh paintor_cred/input/ paintor_cred/out_nc2/model2/$ann 2 $ANNOTS"
done
grep "Successfully" FarmOut/run_paintor_cred_model2*txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/run_paintor_cred_model2*txt | wc -l
cat $ROOT/paintor_cred/out_nc2/model2/*/*.Values > $ROOT/paintor_cred/out_nc2/model2/model2.enrichments.txt
grep "" $ROOT/paintor_cred/out_nc2/model2/*/Log.BayesFactor | sed 's/\.\///g' | sed 's/\/Log.BayesFactor:/\t/g' | sed 's/.*\///g' > $ROOT/paintor_cred/out_nc2/model2/model2.LogBayesFactor.txt


################# MODEL 3
# Select the coding_nonsyn annotation, and run PAINTOR again with this
# annotation plus each other annotation individually.
mkdir paintor_cred/out_nc2/model3
ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.01 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia ipsMacrophage microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
#ANNOTATIONS=( DNase )
for ann in "${ANNOTATIONS[@]}"; do
  ANNOTS="Blood_Immune_DNase,coding_nonsyn,$ann"
  submitJobs.py --MEM 500 -j run_paintor_cred_model3.$ann -q normal \
    -c "bash $SRC/run_paintor.sh paintor_cred/input/ paintor_cred/out_nc2/model3/$ann 2 $ANNOTS"
done
grep "Successfully" FarmOut/run_paintor_cred_model3*txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/run_paintor_cred_model3*txt | wc -l
cat $ROOT/paintor_cred/out_nc2/model3/*/*.Values > $ROOT/paintor_cred/out_nc2/model3/model3.enrichments.txt
grep "" $ROOT/paintor_cred/out_nc2/model3/*/Log.BayesFactor | sed 's/\.\///g' | sed 's/\/Log.BayesFactor:/\t/g' | sed 's/.*\///g' > $ROOT/paintor_cred/out_nc2/model3/model3.LogBayesFactor.txt


################# MODEL 4
# Select the spliceai_gt_0.01 annotation, and run PAINTOR again with this
# annotation plus each other annotation individually.
mkdir paintor_cred/out_nc2/model4
ANNOTATIONS=( intronic upstream_gene downstream_gene regulatory_region non_coding_transcript_exon three_prime_UTR missense synonymous five_prime_UTR TF_binding_site UTR exonic regulatory gene_proximal phastCons_gt_0.1 phastCons_gt_0.5 phastCons_gt_0.95 phyloP_gt_0.5 phyloP_gt_1 phyloP_gt_2 GERP_gt_1 GERP_gt_2 GERP_gt_3 spliceai_gt_0 spliceai_gt_0.1 deepSEA_lt_0.1 deepSEA_lt_0.05 deepSEA_lt_0.01 captureSeqSpliceDist_lt_10 CADD_PHRED_gt_5 CADD_PHRED_gt_10 CADD_PHRED_gt_20 microglia ipsMacrophage microglia_or_macrophage_atac DNase DNase_gteq10 Brain_DNase Roadmap_Enh Roadmap_Enh_gteq10 Brain_Enh Blood_Immune_Enh )
#ANNOTATIONS=( upstream_gene )
for ann in "${ANNOTATIONS[@]}"; do
  ANNOTS="Blood_Immune_DNase,coding_nonsyn,spliceai_gt_0.01,$ann"
  submitJobs.py --MEM 500 -j run_paintor_cred_model4.$ann -q normal \
    -c "bash $SRC/run_paintor.sh paintor_cred/input/ paintor_cred/out_nc2/model4/$ann 2 $ANNOTS"
done
grep "Successfully" FarmOut/run_paintor_cred_model4*txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/run_paintor_cred_model4*txt | wc -l
cat $ROOT/paintor_cred/out_nc2/model4/*/*.Values > $ROOT/paintor_cred/out_nc2/model4/model4.enrichments.txt
grep "" $ROOT/paintor_cred/out_nc2/model4/*/Log.BayesFactor | sed 's/\.\///g' | sed 's/\/Log.BayesFactor:/\t/g' | sed 's/.*\///g' > $ROOT/paintor_cred/out_nc2/model4/model4.LogBayesFactor.txt


# STOP selecting annotations, as all of the top ones from model 4 have unreasonable enrichments
# So the model outputs that we want to use are in:
# paintor_cred/out_nc2/model3/spliceai_gt_0.01


# Put together PAINTOR outputs for all loci
paste <(echo "locus") <(head -n 1 paintor_cred/out_nc2/model3/spliceai_gt_0.01/1_161155392.IGAP1_GWAX.results) > paintor_cred/out_nc2/model3/all_loci.results
for f in paintor_cred/out_nc2/model3/spliceai_gt_0.01/*.IGAP1_GWAX.results; do
  #echo $f
  locusname=`echo $f | perl -ne '@A=split(/\//); @B=split(/\./, $A[$#A]); print join(".", @B[0..1])."\n";'`
  echo $locusname
  sed '1d' $f | perl -sne 'print $locus." ".$_' -- -locus=$locusname >> paintor_cred/out_nc2/model3/all_loci.results
done
sed -i 's/ /\t/g' paintor_cred/out_nc2/model3/all_loci.results


# Add the PAINTOR probabilities to the main annotated file
Rscript $SRC/add_paintor_probabilities.R --file annotated/$GWAS_NAME.annotated.selected.probable.tsv --paintor_file paintor_cred/out_nc2/model3/all_loci.results > annotated/$GWAS_NAME.annotated.selected.probable.paintor.tsv
# Check that we have only changed the file by adding one column
diff <(cut -f 1-24,30-107 annotated/$GWAS_NAME.annotated.selected.probable.tsv) <(cut -f 1-24,32-109 annotated/$GWAS_NAME.annotated.selected.probable.paintor.tsv) | head

# When redoing the annotation, add back in the handmade notes we had added
Rscript $SRC/merge_previous_annotated_table.R


################################################################################
# COLOC: colocalisation of associations with molecular QTLs using summary stats

# Convert AD summary stats into format needed for our coloc scripts
mkdir coloc/
mkdir coloc/GWAS
F=$ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz
#gzhead 1 $F
#CHR	BP	SNP	A1	A2	GWAS_BETA	GWAS_SE	GWAS_P	GWAX_BRITISH_BETA	GWAX_BRITISH_SE	GWAX_BRITISH_P	GWAX_EU_BETA	GWAX_EU_SE	GWAX_EU_P	META_BETA	META_SE	META_P	DIRECT	I2	HET_P	FREQ	INFO
# Prepare GWAS file in input format expected by the coloc scripts
time (echo -e "SNP\tCHR\tBP\tA1\tMAF\tMETA_P\tMETA_BETA\tOR\tlog_OR\tMETA_SE\tz_score\ttrait\tPMID\tused_file"; 
 zcat $F | sed '1d' \
 | perl -ane 'print join("\t", $F[2],$F[0],$F[1],$F[3],"",$F[13],$F[11],"","",$F[12],"","AD","","")."\n";' \
 | sort -k2,2 -k3,3n) \
 | bgzip > coloc/GWAS/Alzheimers_disease_Liu_2019_UKBB_GWAX_meta_v5.sorted.txt.gz

tabix -f -s 2 -b 3 -e 3 -S 1 coloc/GWAS/Alzheimers_disease_Liu_2019_UKBB_GWAX_meta_v5.sorted.txt.gz

echo -e "chr\tpos\trsid\tlocus_name\tp" > coloc/AD.signals.for_coloc.txt
sed '1d' AD.loci.tsv | perl -ane 'print join("\t", $F[0], $F[2], $F[1], $F[14], $F[10])."\n";' >> coloc/AD.signals.for_coloc.txt


# Also prepare a file of conditioned p values representing the "first" independent
# association at each locus that has multiple causal variants.
GCTA_OUT=gcta/output_1e-5/cond/
primary_cond_signal_files=( $GCTA_OUT/ABCA7.cond.out.rs12151021_A_G.flt.tsv \
                            $GCTA_OUT/ACE.cond.out.rs4311_T_C.flt.tsv \
                            $GCTA_OUT/ADAM10.cond.out.rs442495_T_C.flt.tsv \
                            $GCTA_OUT/APP-ADAMTS1.cond.out.rs4817090_T_C.flt.tsv \
                            $GCTA_OUT/BIN1.cond.out.rs6733839_C_T.flt.tsv \
                            $GCTA_OUT/EPHA1.cond.out.rs12703526_G_T.flt.tsv \
                            $GCTA_OUT/NCK2.cond.out.rs143080277_T_C.flt.tsv \
                            $GCTA_OUT/PTK2B-CLU.cond.out.rs867230_C_A.flt.tsv \
                            $GCTA_OUT/TREM2.cond.out.rs187370608_G_A.flt.tsv )

# Cols of GCTA *.flt.tsv files:
# Chr	SNP	bp	refA	freq	b	se	p	n	freq_geno	bC	bC_se	pC
# When preparing summary stats file for coloc, need to keep SNP ID without allele
# as part of the name.
(echo -e "SNP\tCHR\tBP\tA1\tMAF\tMETA_P\tMETA_BETA\tOR\tlog_OR\tMETA_SE\tz_score\ttrait\tPMID\tused_file";
 (for FILE in "${primary_cond_signal_files[@]}"; do
     sed '1d' $FILE
  done;) | perl -ane '@s=split("_",$F[1]); $snpid=$s[0]; print join("\t", $snpid,$F[0],$F[2],$F[3],$F[4],$F[12],$F[10],"","",$F[11],"","AD","","")."\n";' \
  | sort -k2,2 -k3,3n) \
  | bgzip > $GCTA_OUT/AD.meta.cond_signals_1.tsv.gz
tabix -f -s 2 -b 3 -e 3 -S 1 $GCTA_OUT/AD.meta.cond_signals_1.tsv.gz

secondary_cond_signal_files=( $GCTA_OUT/ABCA7.cond.out.rs4147918_A_G.flt.tsv \
                              $GCTA_OUT/ACE.cond.out.rs3730025_A_G.flt.tsv \
                              $GCTA_OUT/ADAM10.cond.out.rs4775044_T_G.flt.tsv \
                              $GCTA_OUT/APP-ADAMTS1.cond.out.rs2830489_C_T.flt.tsv \
                              $GCTA_OUT/BIN1.cond.out.rs7584040_C_T.flt.tsv \
                              $GCTA_OUT/EPHA1.cond.out.rs10265814_C_G.flt.tsv \
                              $GCTA_OUT/NCK2.cond.out.rs116038905_T_C.flt.tsv \
                              $GCTA_OUT/PTK2B-CLU.cond.out.rs73223431_C_T.flt.tsv \
                              $GCTA_OUT/TREM2.cond.out.rs3857580_A_G.flt.tsv )
                              
 (echo -e "SNP\tCHR\tBP\tA1\tMAF\tMETA_P\tMETA_BETA\tOR\tlog_OR\tMETA_SE\tz_score\ttrait\tPMID\tused_file";
 (for FILE in "${secondary_cond_signal_files[@]}"; do
     sed '1d' $FILE
  done;) | perl -ane '@s=split("_",$F[1]); $snpid=$s[0]; print join("\t", $snpid,$F[0],$F[2],$F[3],$F[4],$F[12],$F[10],"","",$F[11],"","AD","","")."\n";' \
  | sort -k2,2 -k3,3n) \
  | bgzip > $GCTA_OUT/AD.meta.cond_signals_2.tsv.gz
tabix -f -s 2 -b 3 -e 3 -S 1 $GCTA_OUT/AD.meta.cond_signals_2.tsv.gz

cond_signal_files_3=( $GCTA_OUT/TREM2.cond.out.rs114812713_G_C.flt.tsv )
(echo -e "SNP\tCHR\tBP\tA1\tMAF\tMETA_P\tMETA_BETA\tOR\tlog_OR\tMETA_SE\tz_score\ttrait\tPMID\tused_file";
 (for FILE in "${cond_signal_files_3[@]}"; do
     sed '1d' $FILE
  done;) | perl -ane '@s=split("_",$F[1]); $snpid=$s[0]; print join("\t", $snpid,$F[0],$F[2],$F[3],$F[4],$F[12],$F[10],"","",$F[11],"","AD","","")."\n";' \
  | sort -k2,2 -k3,3n) \
  | bgzip > $GCTA_OUT/AD.meta.cond_signals_3.tsv.gz
tabix -f -s 2 -b 3 -e 3 -S 1 $GCTA_OUT/AD.meta.cond_signals_3.tsv.gz


# Manually edit the file AD.meta.cond.indep_signals.tsv to produce a separate
# file for each level of conditional signal, i.e. primary signal, secondary, etc.


################################################################################
# Manually RUN lines in the following files:

$SRC/coloc/prepare_qtl_data.sh

$SRC/coloc/run_coloc.sh
################################################################################

# Merge all coloc results in to "*.coloc_details.txt" file, and add a column
# with top colocs to the main annotated SNP file (only lead SNPs have details
# of colocs added in). Also, fix any columns that were corrupted by being treated
# as numbers when they are really text (e.g. Pubmed IDs)
#head -n 1 annotated/AD.credible_sets.vep_ann.tsv | tr '\t' '\n'
Rscript $SRC/coloc/annotate.coloc.R $ROOT coloc.AD.meta $ROOT/annotated/$GWAS_NAME.annotated.selected.probable.tsv $ROOT/annotated/$GWAS_NAME.annotated

Rscript $SRC/coloc/get_coloc_score.R $ROOT/coloc/output/coloc.AD.meta.eqtl_sqtl_colocs.txt $ROOT/coloc/dataset_groups.tsv $ROOT/coloc/coloc_group_weights.tsv $ROOT/coloc/coloc.AD.meta

# Make some lists of coloc genes to use in network analyses and downstream
# gene scores
Rscript $SRC/coloc/filter_coloc_genes.R


################################################################################
# NETWORK analyses

# Get seed genes for network analysis
Rscript $SRC/merge_seed_gene_sets.R $ROOT/genes/literature_genes.tsv $ROOT/genes/AD.loci.nearest_gene.pc.tsv \
        $ROOT/coloc/output/coloc.AD.meta.merged_cond.H4_gt_0.5.all.tsv $ROOT/network/network_genes.txt \
    > $ROOT/network/AD.network_seed_genes.tsv
    
# Inigo ran the network propagation steps and permutations, and the results are
# in the file network/Zsco_node_TSPOAP1_1000ite.csv

# Run manually. Part of this needs to be run first, and then the below code,
# before running network_analysis.Rmd completely.
$SRC/network_analysis.Rmd

# This outputs top lists of network-connected genes. We want to see if these are
# enriched in nearby low p-value SNPs.
zcat summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3,$14}' | gzip > summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.bed.gz

#bedtools intersect -a network/ad_network.top100.10kb_window.bed -b summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.bed.gz -c > network/ad_network.top100.10kb_window.snp_overlap_counts.tsv
bedtools intersect -a network/ad_network.all.10kb_window.bed -b summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.bed.gz -wa -wb -loj | gzip > network/ad_network.all.10kb_window.snp_overlaps.tsv.gz


################################################################################
# Run fgwas to get genome-wide enrichments
# Manually run lines in the following files

$SRC/run_fgwas.sh
$SRC/plot_fgwas_enrichments.Rmd

# Gene expression can annotate any gene in the genome, and we can use this later
# on to prioritise genes. We've just gotten relative gene expression for single-cell
# brain cell types, and for microglia and other tissues relative to GTEx. Here we
# annotate these datasets with the minimum GWAS p value for SNPs within 10 kb.

# Use the bed file of SNP positions an p values defined above
#zcat summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3,$14}' | gzip > summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.bed.gz

# Get bed files of 10 kb window around all GRCh37 genes
Rscript $SRC/get_genes_bed.R
(echo -e "chr\tgene_window_start\tgene_window_end\tgene_id\tsnp_chr\tsnp_pos\tsnp_end\tsnp_id\tpvalue"; 
 bedtools intersect -a genes/genes.pc.10kb_window.grch37.bed -b summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.bed.gz -wa -wb -loj) | gzip > genes/genes.pc.10kb_window.snp_overlaps.tsv.gz


################################################################################
# Locus gene rankings

# Run manually
merge_gene_evidence.Rmd


################################################################################
# LD score regression

# Do these steps in a separate dir
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file environment.yml
# When I ran munge_sumstats.py it would hang. This seemed to be a numpy problem / conflict.
# The below fixed it.
conda install numpy=1.13.3

cd $ROOT
mkdir ldsc; cd ldsc
# Use the precomputed LD scores provided by alkesgroup
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
tar -jxvf eur_w_ld_chr.tar.bz2
bunzip2 w_hm3.snplist.bz2

# Convert Kunkle sumstats to a format accepted by LDSC
# N for Kunkle: 21982 + 41944 = 63926
(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{N=63926; print $3,$4,$5,$6,$8,N}') \
 > AD.kunkle.sumstats_in.tsv
# zcat summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | awk 'BEGIN{OFS="\t"}{N=63926; print $3,$4,$5,$6,$14,N}') | head

source activate ldsc

$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.kunkle.sumstats_in.tsv \
--N 63926 \
--out AD.kunkle \
--merge-alleles w_hm3.snplist

# LD Score Regression - heritability
$SWDIR/ldsc/ldsc.py \
--h2 AD.kunkle.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.kunkle.h2

# Convert AD proxy sumstats to a format accepted by LDSC
# N for AD proxy: (54939+898)/4+(355900)/4 = 102934
#(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
# zcat $ROOT//summary_stats/AD.proxy_exclude_firsts.bgen.stats.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{N=102934; print $1,$5,$6,$11,$16,N}') \
# > AD.proxy.sumstats_in.tsv

(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{N=102934; print $3,$4,$5,$9,$11,N}') \
 > AD.proxy.sumstats_in.tsv

$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.proxy.sumstats_in.tsv \
--N 102934 \
--out AD.proxy \
--merge-alleles w_hm3.snplist

# LD Score Regression - heritability
$SWDIR/ldsc/ldsc.py \
--h2 AD.proxy.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.proxy.h2

# LD Score Regression - heritability for both AD proxy & Kunkle, and genetic correlation
$SWDIR/ldsc/ldsc.py \
--rg AD.kunkle.sumstats.gz,AD.proxy.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.kunkle_and_proxy.rg


# Convert meta-analysis sumstats to a format accepted by LDSC
# N for meta: N for Kunkle + N for proxy = 21982+41944+(54939+898)/4+(355900)/4 = 166860
(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{N=166860; print $3,$4,$5,$12,$14,N}') \
 > AD.meta.sumstats_in.tsv
$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.meta.sumstats_in.tsv \
--N 166860 \
--out AD.meta \
--merge-alleles w_hm3.snplist

# LD Score Regression - AD proxy
$SWDIR/ldsc/ldsc.py \
--h2 AD.meta.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.meta.h2


# Check for the exclude_seconds sumstats
(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_seconds_v5.meta.tsv.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{N=102934; print $3,$4,$5,$9,$11,N}') \
 > AD.proxy.exclude_seconds.sumstats_in.tsv

$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.proxy.exclude_seconds.sumstats_in.tsv \
--N 102934 \
--out AD.proxy.exclude_seconds \
--merge-alleles w_hm3.snplist

# LD Score Regression - heritability
$SWDIR/ldsc/ldsc.py \
--h2 AD.proxy.exclude_seconds.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.proxy.exclude_seconds.h2


#############################################
# Do it all again, excluding the APOE locus
# Removing SNPs at chr19:42,000,000-48,000,000

# Convert Kunkle sumstats to a format accepted by LDSC
# N for Kunkle: 21982 + 41944 = 63926
(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk '$1 != 19 || $2 < 42000000 || $2 > 48000000' | awk 'BEGIN{OFS="\t"}{N=63926; print $3,$4,$5,$6,$8,N}') \
 > AD.kunkle.sumstats_in.excl_APOE.tsv
# zcat summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | awk 'BEGIN{OFS="\t"}{N=63926; print $3,$4,$5,$6,$14,N}') | head

source activate ldsc

$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.kunkle.sumstats_in.excl_APOE.tsv \
--N 63926 \
--out AD.kunkle.excl_APOE \
--merge-alleles w_hm3.snplist

# LD Score Regression - heritability
$SWDIR/ldsc/ldsc.py \
--h2 AD.kunkle.excl_APOE.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.kunkle.excl_APOE.h2

# Convert AD proxy sumstats to a format accepted by LDSC
# N for AD proxy: (54939+898)/4+(355900)/4 = 102934
#(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
# zcat $ROOT//summary_stats/AD.proxy_exclude_firsts.bgen.stats.bgz | sed '1d' | awk 'BEGIN{OFS="\t"}{N=102934; print $1,$5,$6,$11,$16,N}') \
# > AD.proxy.sumstats_in.tsv

(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk '$1 != 19 || $2 < 42000000 || $2 > 48000000' | awk 'BEGIN{OFS="\t"}{N=102934; print $3,$4,$5,$9,$11,N}') \
 > AD.proxy.excl_APOE.sumstats_in.tsv

$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.proxy.excl_APOE.sumstats_in.tsv \
--N 102934 \
--out AD.proxy.excl_APOE \
--merge-alleles w_hm3.snplist

# LD Score Regression - heritability
$SWDIR/ldsc/ldsc.py \
--h2 AD.proxy.excl_APOE.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.proxy.excl_APOE.h2

# LD Score Regression - heritability for both AD proxy & Kunkle, and genetic correlation
$SWDIR/ldsc/ldsc.py \
--rg AD.kunkle.excl_APOE.sumstats.gz,AD.proxy.excl_APOE.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.kunkle_and_proxy.excl_APOE.rg


# Convert meta-analysis sumstats to a format accepted by LDSC
# N for meta: N for Kunkle + N for proxy = 21982+41944+(54939+898)/4+(355900)/4 = 166860
(echo -e "SNP\tA1\tA2\tbeta\tp\tN";
 zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | sed '1d' | awk '$1 != 19 || $2 < 42000000 || $2 > 48000000' | awk 'BEGIN{OFS="\t"}{N=63926; print $3,$4,$5,$12,$14,N}') \
 > AD.meta.excl_APOE.sumstats_in.tsv
$SWDIR/ldsc/munge_sumstats.py \
--sumstats AD.meta.excl_APOE.sumstats_in.tsv \
--N 166860 \
--out AD.meta.excl_APOE \
--merge-alleles w_hm3.snplist

# LD Score Regression - AD proxy
$SWDIR/ldsc/ldsc.py \
--h2 AD.meta.excl_APOE.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.meta.excl_APOE.h2

#############################################
# Do LDSC for the Gr@ace study
cd $ROOT/ldsc
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/Moreno-GrauS_31473137_GCST009020/GRACE_StageI.txt
gzip GRACE_StageI.txt
(echo -e "CHR\tPos\tA1\tA2\tBeta\tSE\tP\tDirection\trsID"; zcat GRACE_StageI.txt.gz | sed '1d') > GRACE_StageI.sumstats_in.tsv
#wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/Moreno-GrauS_31473137_GCST009019/GRACEStageI_dbGAP.txt
#gzip GRACEStageI_dbGAP.txt
#(echo -e "CHR\tPos\tA1\tA2\tBeta\tSE\tP\tDirection\trsID"; zcat GRACEStageI_dbGAP.txt.gz | sed '1d') > GRACEStageI_dbGAP.sumstats_in.tsv

# N = 4120 cases + 3289 controls = 7409
$SWDIR/ldsc/munge_sumstats.py \
--sumstats GRACE_StageI.sumstats_in.tsv \
--N 7409 \
--out AD.Grace \
--merge-alleles w_hm3.snplist

$SWDIR/ldsc/ldsc.py \
--h2 AD.Grace.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.Grace.h2

$SWDIR/ldsc/ldsc.py \
--rg AD.kunkle.sumstats.gz,AD.proxy.sumstats.gz,AD.Grace.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.kunkle_proxy_grace.rg

# Removing SNPs at chr19:42,000,000-48,000,000
(echo -e "CHR\tPos\tA1\tA2\tBeta\tSE\tP\tDirection\trsID";
 zcat GRACE_StageI.txt.gz | sed '1d' | awk '$1 != 19 || $2 < 42000000 || $2 > 48000000') \
 > GRACE_StageI.sumstats_in.excl_APOE.tsv

$SWDIR/ldsc/munge_sumstats.py \
--sumstats GRACE_StageI.sumstats_in.excl_APOE.tsv \
--N 7409 \
--out AD.Grace.excl_APOE \
--merge-alleles w_hm3.snplist

$SWDIR/ldsc/ldsc.py \
--h2 AD.Grace.excl_APOE.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.Grace.excl_APOE.h2

$SWDIR/ldsc/ldsc.py \
--rg AD.kunkle.excl_APOE.sumstats.gz,AD.proxy.excl_APOE.sumstats.gz,AD.Grace.excl_APOE.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out AD.kunkle_proxy_grace.excl_APOE.rg

# The results show unreasonably high h2 for the Gr@ace study, and with a high ratio
# indicating that it is not due to real polygenic signal. When I run the genetic
# correlation analysis, it dies with an error for the Gr@ace study.


###############################################################################
# Check for replication of AD loci in other datasets
mkdir $ROOT/replication
cd $ROOT/replication

ln -s $ROOT/ldsc/GRACE_StageI.txt.gz GRACE_StageI.txt.gz
# Get lead SNP per locus
sed '1d' $ROOT/AD.loci.tsv | cut -f 2 | sed 's/_.*//g' > AD_lead_snps.txt
# Also get all independent SNPs across loci
sed '1d' $ROOT/AD.loci.tsv | awk '$15 != "APOE"' | cut -f 4 | sed 's/_[^,]*//g;s/,/\n/g' > AD_indep_snps.txt
sed '1d' $ROOT/AD.loci.tsv | awk '$15 == "APOE"' | cut -f 2 | sed 's/_.*//g' >> AD_indep_snps.txt

# Manually add in linked SNPs to account for some that are missing in FinnGen
# AD_lead_snps.plus.txt
# AD_indep_snps.plus.txt

# Gr@ace Spanish AD GWAS
(zcat GRACE_StageI.txt.gz | head -n 1; zgrep -wFf AD_lead_snps.plus.txt GRACE_StageI.txt.gz) | tr ' ' '\t' > GRACE_StageI.AD_lead_snps.txt
(zcat GRACE_StageI.txt.gz | head -n 1; zgrep -wFf AD_indep_snps.plus.txt GRACE_StageI.txt.gz) | tr ' ' '\t' > GRACE_StageI.AD_indep_snps.txt

# FinnGen - two definitions of AD
wget https://storage.googleapis.com/finngen-public-data-r3/summary_stats/finngen_r3_AD_LO.gz
wget https://storage.googleapis.com/finngen-public-data-r3/summary_stats/finngen_r3_G6_AD_WIDE.gz
(zcat finngen_r3_G6_AD_WIDE.gz | head -n 1; zgrep -wFf AD_lead_snps.plus.txt finngen_r3_G6_AD_WIDE.gz) > finngen_r3_G6_AD_WIDE.AD_lead_snps.txt
(zcat finngen_r3_G6_AD_WIDE.gz | head -n 1; zgrep -wFf AD_lead_snps.plus.txt finngen_r3_AD_LO.gz) > finngen_r3_AD_LO.AD_lead_snps.txt

(zcat finngen_r3_G6_AD_WIDE.gz | head -n 1; zgrep -wFf AD_indep_snps.plus.txt finngen_r3_G6_AD_WIDE.gz) > finngen_r3_G6_AD_WIDE.AD_indep_snps.txt
(zcat finngen_r3_G6_AD_WIDE.gz | head -n 1; zgrep -wFf AD_indep_snps.plus.txt finngen_r3_AD_LO.gz) > finngen_r3_AD_LO.AD_indep_snps.txt

(zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | head -n 1;
 zgrep -wFf AD_lead_snps.plus.txt $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz) \
 > AD.IGAP1_GWAX_exclude_firsts_v5.meta.AD_lead_snps.txt

(zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | head -n 1;
 zgrep -wFf AD_indep_snps.plus.txt $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz) \
 > AD.IGAP1_GWAX_exclude_firsts_v5.meta.AD_indep_snps.txt


###############################################################################
# Meta-analysis of AD loci with other datasets
cd $ROOT/replication

# Tabix index the datasets so we can more easily run one chr at a time
(zcat GRACE_StageI.txt.gz | head -n 1 | sed 's/ /\t/g';
 zcat GRACE_StageI.txt.gz | sed '1d' | sed 's/ /\t/g' | sort -k1,1 -k2,2n) | bgzip > GRACE_StageI.txt.bgz
tabix -h -s 1 -b 2 -e 2 -S 1 GRACE_StageI.txt.bgz

# Prepare FinnGen GWAS for reading with Tabix
zcat finngen_r3_G6_AD_WIDE.gz | bgzip > finngen_r3_G6_AD_WIDE.bgz
tabix -s 1 -b 2 -e 2 -S 1 finngen_r3_G6_AD_WIDE.bgz

# Run meta-analysis in R
Rscript $SRC/do_replication_meta.R

# Merge together per-chromosome results
(zcat global_meta/chr1.global_meta.tsv.gz | head -n 1;
 zcat global_meta/chr*.tsv.gz | grep -v "^pos_id" | sort -k2,2n -k3,3n) \
 | bgzip > global_meta/global_meta.tsv.gz
tabix -s 2 -b 3 -e 3 -S 1 global_meta/global_meta.tsv.gz

zcat global_meta/global_meta.tsv.gz | head -n 1 > global_meta/global_meta.header.tsv

# Extract significant SNPs in the global meta-analysis
(zcat global_meta/global_meta.tsv.gz | head -n 1;
 zcat global_meta/global_meta.tsv.gz | sed '1d' | awk '$10 < 5e-8') > global_meta/global_meta.p_lt_5e-8.tsv

(zcat global_meta/global_meta.tsv.gz | head -n 1;
 zcat global_meta/global_meta.tsv.gz | sed '1d' | awk '$15 < 5e-8') > global_meta/global_meta.kfu_p_lt_5e-8.tsv
# Manually examine these files to identify independent genome-wide significant loci
# This gives 4 additional loci: GRN, IGHG1, SHARPIN, SIGLEC11
# SHARPIN: 8	145158607
# IGHG1: 14	106195719
# GRN: 17	42430244
# SIGLEC11: 19	50451508

# Extract just our 36 loci of interest for simpler downstream processing
sed '1d' $ROOT/AD.loci.tsv | awk 'BEGIN{OFS="\t"}{print $1,$6-500000,$7+500000,$2,$15}' > AD.loci.1Mb_window.bed
# Manually add the 4 loci above to AD.loci.1Mb_window.bed

zcat global_meta/global_meta.tsv.gz | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{print $2,$3,$3+1,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43}' \
  | gzip > global_meta/global_meta.bed.gz



bedtools intersect -b AD.loci.1Mb_window.bed -a global_meta/global_meta.bed.gz -wa -wb \
    | gzip > global_meta/global_meta.assoc_loci.bed.gz
(paste <(zcat global_meta/global_meta.tsv.gz | head -n 1) <(echo "locus");
 zcat global_meta/global_meta.assoc_loci.bed.gz \
 | awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$49}') \
 | gzip > global_meta/global_meta.assoc_loci.tsv.gz



Rscript $SRC/check_replication.R

Rscript $SRC/compare_replication_meta.R


###############################################################################
# Apply conditional analysis and fine-mapping to the Kunkle dataset separately,
# to check how robust our identified signals are.

cd $ROOT

submitJobs.py --MEM 9000 -j get_gcta_input.kunkle -q yesterday \
  -c "python $SRC/get_gcta_input.kunkle.py"

# Run GCTA-COJO on each chromosome. Use  P<1×10−5  p-value cutoff, since we can
# later subset to only genome-wide significant loci.
for filepath in gcta/input/*.bed; do
  echo "$filepath"
  filename=$(basename "$filepath") # filename like ukbb_sample.1.156156033_166156033.bed
  filepath_no_ext="${filepath%.*}"
  chr=`echo $filename | perl -ne '@Ar=split(/\./); print $Ar[1];'`
  
  submitJobs.py --MEM 2000 -q normal -j gcta.1e-5.$chr \
    -c "gcta64 \
  --bfile $filepath_no_ext \
  --chr $chr \
  --cojo-file gcta/input/AD.Kunkle.$chr.ma \
  --cojo-slct \
  --cojo-p 1e-5 \
  --out gcta/output_1e-5/AD.Kunkle.$chr"
done

# For APOE and HLA, do a special run with a threshold of 5e-8 to avoid overfitting
# the strong signal based on imperfect LD match
mkdir gcta/output_5e-8
chr_paths=( gcta/input/ukbb_sample.6.merged.bed gcta/input/ukbb_sample.19.merged.bed )
for filepath in "${chr_paths[@]}"; do
  echo "$filepath"
  filename=$(basename "$filepath") # filename like ukbb_sample.1.156156033_166156033.bed
  filepath_no_ext="${filepath%.*}"
  chr=`echo $filename | perl -ne '@Ar=split(/\./); print $Ar[1];'`
  
  submitJobs.py --MEM 2000 -q normal -j gcta.5e-8.$chr \
    -c "gcta64 \
  --bfile $filepath_no_ext \
  --chr $chr \
  --cojo-file gcta/input/AD.Kunkle.$chr.ma \
  --cojo-slct \
  --cojo-p 5e-8 \
  --out gcta/output_5e-8/AD.Kunkle.$chr"
done

grep "Successfully" FarmOut/gcta*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code" FarmOut/gcta*.txt | wc -l

# This python script reads the output from GCTA (per chromosome gcta_out/AD.IGAP1_GWAX.[chr].jma.cojo files) and
# creates list of independent loci. If distance between two SNPs is < 500kb, they will be collapsed into a single
# locus. Remove SNPs with marginal P>5×10−8. Output file (in this example, "AD.IGAP1_GWAX.indep.hits") columns are
# "Chr lead_SNP lead_pos SNPs pos start stop n_snps lead_freq lead_beta lead_p freq beta p".
python $SRC/get_gcta_output_loci.kunkle.py --gcta_dir gcta/output_1e-5 --pthresh 5e-8 --out AD.Kunkle.gwsig_indep.hits
# Use a distance of 620 kb so that APP-ADAMTS1 is considered 1 locus, as variants
# are in partial LD
python $SRC/get_gcta_output_loci.kunkle.py --gcta_dir gcta/output_1e-5 --pthresh 1e-5 --window 620000 --out AD.Kunkle.1e-5_indep.hits
python $SRC/get_gcta_output_loci.kunkle.py --gcta_dir gcta/output_1e-5 --pthresh 1 --window 620000 --out AD.Kunkle.all_indep.hits
python $SRC/get_gcta_output_loci.kunkle.py --gcta_dir gcta/output_5e-8 --pthresh 5e-8 --out AD.Kunkle.APOE_HLA.hits


###############################
# FINEMAP
mkdir $ROOT/finemap_kunkle
mkdir $ROOT/finemap_kunkle/input

FINEMAP_DIR=$ROOT/finemap_kunkle

(head -n 1 AD.loci.tsv; cat AD.loci.tsv | grep "HLA") > AD.loci.HLA.tsv
cat AD.loci.tsv | grep -v "HLA" > AD.loci.exceptHLA.tsv

(head -n 1 AD.loci.tsv; cat AD.loci.tsv | grep "TREM2") > AD.loci.TREM2.tsv

# Prepare input files for running FINEMAP (and PAINTOR). This takes a couple of hours.
# We only include SNPs with MAF > 0.002, since this is at the limit of our accuracy
# for LD estimation (20 individuals in 10,000).
submitJobs.py --MEM 10000 -n 5 -j get_finemapping_input.exceptHLA.kunkle -q yesterday \
  -c "python $SRC/get_finemapping_input.kunkle.py --locus_file AD.loci.exceptHLA.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.002 --min_info 0.85 --hetp_threshold 0.001 --outdir $FINEMAP_DIR/input"
submitJobs.py --MEM 19000 -n 5 -j get_finemapping_input.HLA.kunkle -q normal \
  -c "python $SRC/get_finemapping_input.kunkle.py --locus_file AD.loci.HLA.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.002 --min_info 0.85 --hetp_threshold 0.001 --outdir $FINEMAP_DIR/input"
submitJobs.py --MEM 10000 -n 5 -j get_finemapping_input.TREM2.kunkle -q normal \
  -c "python $SRC/get_finemapping_input.kunkle.py --locus_file AD.loci.TREM2.tsv --gwas_file $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz --window 500000 --min_freq 0.0019 --min_info 0.85 --hetp_threshold 0.001 --outdir $FINEMAP_DIR/input"


# Determine how many SNPs were excluded by different criteria. First paste together
# all of the *.excluded_snps files and *.snplist files that were output.
cat $FINEMAP_DIR/input/*.excluded_snps > $FINEMAP_DIR/input/all_excluded_snps.txt
cat $FINEMAP_DIR/input/*.snplist > $FINEMAP_DIR/input/all.snplist.txt
Rscript $SRC/count_snp_exclusions.R

# The above step also extracts the meta-analysis SNPs at each associated locus,
# and now we merge these into a single file for later use.
(zcat $ROOT/summary_stats/AD.IGAP1_GWAX_exclude_firsts_v5.meta.tsv.bgz | head -n 1;
 zcat $FINEMAP_DIR/input/*.IGAP1_GWAX.gwas.bgz | sort -k 1,1n -k2,2n) | bgzip > summary_stats/$GWAS_NAME.assoc_loci.gz
tabix -s 1 -b 2 -e 2 -S 1 summary_stats/$GWAS_NAME.assoc_loci.gz

# Run FINEMAP statistical fine-mapping
submitJobs.py --MEM 2000 -j run_finemap_ncausal_1.kunkle -q yesterday \
  -c "python $SRC/run_finemap.kunkle.py --input_loci AD.loci.tsv --input_dir $FINEMAP_DIR/input/ --output_dir $FINEMAP_DIR/out_ncausal_1  --ncausal 1"

submitJobs.py --MEM 2500 -j run_finemap_ncausal_2.kunkle -q normal \
  -c "python $SRC/run_finemap.kunkle.py --input_loci AD.loci.tsv --input_dir $FINEMAP_DIR/input/ --output_dir $FINEMAP_DIR/out_ncausal_2  --ncausal 2"

submitJobs.py --MEM 2000 -j run_finemap_ncausal_3.kunkle -q normal \
  -c "python $SRC/run_finemap.kunkle.py --input_loci AD.loci.TREM2.tsv --input_dir $FINEMAP_DIR/input/ --output_dir $FINEMAP_DIR/out_ncausal_3  --ncausal 3"

# Merge finemap locus results into one file
mkdir $FINEMAP_DIR/output
(head -n 1 $FINEMAP_DIR/out_ncausal_1/1_161155392.Kunkle.1.snp; cat $FINEMAP_DIR/out_ncausal_1/*.Kunkle.1.snp | grep -v "^index") | tr ' ' '\t' > $FINEMAP_DIR/output/$GWAS_NAME.finemap.ncausal_1.snp
(head -n 1 $FINEMAP_DIR/out_ncausal_2/1_161155392.Kunkle.2.snp; cat $FINEMAP_DIR/out_ncausal_2/*.Kunkle.2.snp | grep -v "^index") | tr ' ' '\t' > $FINEMAP_DIR/output/$GWAS_NAME.finemap.ncausal_2.snp
(head -n 1 $FINEMAP_DIR/out_ncausal_3/6_40942196.Kunkle.3.snp; cat $FINEMAP_DIR/out_ncausal_3/*.Kunkle.3.snp | grep -v "^index") | tr ' ' '\t' > $FINEMAP_DIR/output/$GWAS_NAME.finemap.ncausal_3.snp

for filepath in $FINEMAP_DIR/input/*.Kunkle.*.in; do
  echo "$filepath"
  mv $filepath $filepath.tmp
  cat $filepath.tmp | sed '1,2d' > $filepath
done

Rscript $SRC/compare_meta_kunkle_finemap.R


################################################################################
# Paper figures

# Figure 1b - Manhattan plot
Rscript $SRC/draw_manhattan_v5.R
# Figure 1c - based on Supplementary Table 1 - AD loci

# Figure 2 - output file plots/coloc.summary.selected_signals.pdf
Rscript $SRC/coloc/plot_coloc_summary.R

# Figure 3 - output file plots/paintor_finemapping.pdf
$SRC/paintor_plots.Rmd

# Table 1 - selected fine-mapped SNPs from Supplementary Table 7 - SNP fine-mapping

# Figure 4 - run steps in; multiple output files put together
$SRC/finemapping_plots.R

# Figure 5 - produced by merge_gene_evidence.Rmd and plot_fgwas_enrichments.Rmd
# Output files:
# expression/sc_brain.expression_enrichment.pdf
# genes/network_score.minp_enrichment.pdf
$SRC/merge_gene_evidence.Rmd
$SRC/plot_fgwas_enrichments.Rmd

# Figure 6 - run the below. Output file:
# genes/score.breakdown.topGenes.v5.6x5.bars.pdf
$SRC/plot_gene_evidence.Rmd

# SuppFig 1 - Files from merge_gene_evidence.Rmd and plot_fgwas_enrichments.Rmd
# network/network.fgwas_enrichment.pdf
# genes/network.seed_gene_pagerank_pctile.1mb.pdf
# genes/network.seed_gene.violin.pdf

# SuppFig 2 - from merge_gene_evidence.Rmd
# expression/bulk.fgwas_enrichment.pdf

# SuppFig 3 - from merge_gene_evidence.Rmd
# genes/coloc.H4.boxplot.network.pdf
# genes/coloc.H4.boxplot.totalScore.pdf
# genes/coloc.celltype_vs_totalscore.pdf

# SuppFig 4 - from merge_gene_evidence.exprScore.Rmd
# genes/geneDistScore.pdf

# Supp Table 1 - made from file genes/AD.loci.gene_overlaps.pc.tsv, output by annotateRegionGenes.R
# Supp Table 2 - QTL datasets - prepared manually
# Supp Table 3 - QTL coloc summary - file coloc/output/coloc.AD.meta.all_colocs.supp_table.txt, output by annotate.coloc.R
# Supp Table 4 - all QTL colocs - file annotated/AD.meta.annotated.coloc_details.txt, output by annotate.coloc.R
# Supp Table 5 - PAINTOR annotations - file paintor_cred/paintor_annotation_summary.tsv, output by paintor_plots.Rmd

# Supp Table 6 - SNP fine-mapping
# We include a subset of the columns, as many contain information not used / not interesting
cat $ROOT/annotated/$GWAS_NAME.annotated.selected.probable.paintor.with_notes.tsv \
  | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$75,$16,$17,$73,$18,$19,$74,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$38,$39,$40,$41,$43,$44,$45,$46,$47,$55,$57,$58,$59,$64,$65,$66,$67,$68,$69,$70,$102,$103,$104,$105,$106,$107,$108,$109,$110,$111}' \
  > $ROOT/annotated/$GWAS_NAME.annotated.forSuppTable.tsv

# Supp Table 7 - Network seed genes - prepared manually
# Supp Table 8 - Network gene rankings - file network/network.annotated.all.supp_table.tsv, output by network_analysis.Rmd
# Supp Table 9 - Network gene GO terms - file network/network.PRpctile.top1000.gprofiler.tsv, output by network_analysis.Rmd
# Supp Table 10 - Genome-wide enrichments - file 
# Supp Table 11 - Gene evidence rankings - file genes/AD.loci.1Mb_window.expressed_genes.all.geneScores.rounded.tsv, output by merge_gene_evidence.Rmd
# Supp Table 12 - Glmnet models - output by merge_gene_evidence.Rmd

