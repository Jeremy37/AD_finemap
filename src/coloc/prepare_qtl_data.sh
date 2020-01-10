###############################################################################
# This script contains commands to run manually to prepare QTL datasets for
# colocalisation with GWAS summary stats.

ROOT=/path/to/AD_finemap
COLOC=$ROOT/coloc
SRC=$ROOT/src

md $ROOT/coloc/qtl_data
cd $ROOT/coloc/qtl_data

################################################################################
# We have downloaded QTL data from the eQTL catalogue, and so each dataset can
# be processed in a similar way.

DATASETS=( Alasoo_2018_ge_macrophage_IFNg Alasoo_2018_ge_macrophage_IFNg+Salmonella Alasoo_2018_ge_macrophage_naive Alasoo_2018_ge_macrophage_Salmonella BLUEPRINT_PE_ge_T-cell BLUEPRINT_SE_ge_monocyte BLUEPRINT_SE_ge_neutrophil BrainSeq_ge_brain CEDAR_B-cell_CD19 CEDAR_ileum CEDAR_monocyte_CD14 CEDAR_neutrophil_CD15 CEDAR_platelet CEDAR_rectum CEDAR_T-cell_CD4 CEDAR_T-cell_CD8 CEDAR_transverse_colon Fairfax_2012_B-cell_CD19 Fairfax_2014_monocyte_IFN24 Fairfax_2014_monocyte_LPS24 Fairfax_2014_monocyte_LPS2 GENCORD_ge_fibroblast GENCORD_ge_LCL GENCORD_ge_T-cell GEUVADIS_ge_LCL HipSci_ge_iPSC Lepik_2017_ge_blood Nedelec_2016_ge_macrophage_Listeria Nedelec_2016_ge_macrophage_naive Nedelec_2016_ge_macrophage_Salmonella Quach_2016_ge_monocyte_IAV Quach_2016_ge_monocyte_LPS Quach_2016_ge_monocyte_naive Quach_2016_ge_monocyte_Pam3CSK4 Quach_2016_ge_monocyte_R848 ROSMAP_ge_brain_naive Schmiedel_2018_ge_B-cell_naive Schmiedel_2018_ge_CD4_T-cell_anti-CD3-CD28 Schmiedel_2018_ge_CD4_T-cell_naive Schmiedel_2018_ge_CD8_T-cell_anti-CD3-CD28 Schmiedel_2018_ge_CD8_T-cell_naive Schmiedel_2018_ge_monocyte_CD16_naive Schmiedel_2018_ge_monocyte_naive Schmiedel_2018_ge_NK-cell_naive Schmiedel_2018_ge_Tfh_memory Schmiedel_2018_ge_Th1-17_memory Schmiedel_2018_ge_Th17_memory Schmiedel_2018_ge_Th1_memory Schmiedel_2018_ge_Th2_memory Schmiedel_2018_ge_Treg_memory Schmiedel_2018_ge_Treg_naive Schwartzentruber_2018_ge_sensory_neuron TwinsUK_ge_blood TwinsUK_ge_fat TwinsUK_ge_LCL TwinsUK_ge_skin van_de_Bunt_2015_ge_pancreatic_islet )
md eqtl_catalogue
for dat in "${DATASETS[@]}"; do
   ln -s $JS/datasets/eqtl_catalogue/$dat.nominal.sorted.txt.gz
   ln -s $JS/datasets/eqtl_catalogue/$dat.nominal.sorted.txt.gz.tbi
   ln -s $JS/datasets/eqtl_catalogue/$dat.variant_information.txt.gz
done
for dat in "${DATASETS[@]}"; do
   ln -s $COLOC/qtl_data/eqtl_catalogue_bak/$dat.variant_info.GRCh37.txt.gz
done

DATASETS=( Alasoo_2018_ge_macrophage_IFNg Alasoo_2018_ge_macrophage_IFNg+Salmonella Alasoo_2018_ge_macrophage_naive Alasoo_2018_ge_macrophage_Salmonella BLUEPRINT_PE_ge_T-cell BLUEPRINT_SE_ge_monocyte BLUEPRINT_SE_ge_neutrophil BrainSeq_ge_brain CEDAR_B-cell_CD19 CEDAR_ileum CEDAR_monocyte_CD14 CEDAR_neutrophil_CD15 CEDAR_platelet CEDAR_rectum CEDAR_T-cell_CD4 CEDAR_T-cell_CD8 CEDAR_transverse_colon Fairfax_2012_B-cell_CD19 Fairfax_2014_monocyte_IFN24 Fairfax_2014_monocyte_LPS24 Fairfax_2014_monocyte_LPS2 GENCORD_ge_fibroblast GENCORD_ge_LCL GENCORD_ge_T-cell GEUVADIS_ge_LCL HipSci_ge_iPSC Lepik_2017_ge_blood Nedelec_2016_ge_macrophage_Listeria Nedelec_2016_ge_macrophage_naive Nedelec_2016_ge_macrophage_Salmonella Quach_2016_ge_monocyte_IAV Quach_2016_ge_monocyte_LPS Quach_2016_ge_monocyte_naive Quach_2016_ge_monocyte_Pam3CSK4 Quach_2016_ge_monocyte_R848 ROSMAP_ge_brain_naive Schmiedel_2018_ge_B-cell_naive Schmiedel_2018_ge_CD4_T-cell_anti-CD3-CD28 Schmiedel_2018_ge_CD4_T-cell_naive Schmiedel_2018_ge_CD8_T-cell_anti-CD3-CD28 Schmiedel_2018_ge_CD8_T-cell_naive Schmiedel_2018_ge_monocyte_CD16_naive Schmiedel_2018_ge_monocyte_naive Schmiedel_2018_ge_NK-cell_naive Schmiedel_2018_ge_Tfh_memory Schmiedel_2018_ge_Th1-17_memory Schmiedel_2018_ge_Th17_memory Schmiedel_2018_ge_Th1_memory Schmiedel_2018_ge_Th2_memory Schmiedel_2018_ge_Treg_memory Schmiedel_2018_ge_Treg_naive Schwartzentruber_2018_ge_sensory_neuron TwinsUK_ge_blood TwinsUK_ge_fat TwinsUK_ge_LCL TwinsUK_ge_skin van_de_Bunt_2015_ge_pancreatic_islet )
for dat in "${DATASETS[@]}"; do
   echo "$dat"
   submitJobs.py --MEM 5000 -j process_eqtl_catalogue_ge.$dat -q normal \
      -c "$SRC/coloc/process_eqtl_catalogue_ge.sh $COLOC/qtl_data/eqtl_catalogue $dat $COLOC/qtl_data/eqtl_catalogue"
done
grep "Successfully" FarmOut/process_eqtl_catalogue_ge.*txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/process_eqtl_catalogue_ge.*txt | wc -l


# CrossMap.py needs to be run on farm3, it seems, so this is separate
DATASETS=( Alasoo_2018_ge_macrophage_IFNg Alasoo_2018_ge_macrophage_IFNg+Salmonella Alasoo_2018_ge_macrophage_naive Alasoo_2018_ge_macrophage_Salmonella BLUEPRINT_PE_ge_T-cell BLUEPRINT_SE_ge_monocyte BLUEPRINT_SE_ge_neutrophil BrainSeq_ge_brain CEDAR_B-cell_CD19 CEDAR_ileum CEDAR_monocyte_CD14 CEDAR_neutrophil_CD15 CEDAR_platelet CEDAR_rectum CEDAR_T-cell_CD4 CEDAR_T-cell_CD8 CEDAR_transverse_colon Fairfax_2012_B-cell_CD19 Fairfax_2014_monocyte_IFN24 Fairfax_2014_monocyte_LPS24 Fairfax_2014_monocyte_LPS2 GENCORD_ge_fibroblast GENCORD_ge_LCL GENCORD_ge_T-cell GEUVADIS_ge_LCL HipSci_ge_iPSC Lepik_2017_ge_blood Nedelec_2016_ge_macrophage_Listeria Nedelec_2016_ge_macrophage_naive Nedelec_2016_ge_macrophage_Salmonella Quach_2016_ge_monocyte_IAV Quach_2016_ge_monocyte_LPS Quach_2016_ge_monocyte_naive Quach_2016_ge_monocyte_Pam3CSK4 Quach_2016_ge_monocyte_R848 ROSMAP_ge_brain_naive Schmiedel_2018_ge_B-cell_naive Schmiedel_2018_ge_CD4_T-cell_anti-CD3-CD28 Schmiedel_2018_ge_CD4_T-cell_naive Schmiedel_2018_ge_CD8_T-cell_anti-CD3-CD28 Schmiedel_2018_ge_CD8_T-cell_naive Schmiedel_2018_ge_monocyte_CD16_naive Schmiedel_2018_ge_monocyte_naive Schmiedel_2018_ge_NK-cell_naive Schmiedel_2018_ge_Tfh_memory Schmiedel_2018_ge_Th1-17_memory Schmiedel_2018_ge_Th17_memory Schmiedel_2018_ge_Th1_memory Schmiedel_2018_ge_Th2_memory Schmiedel_2018_ge_Treg_memory Schmiedel_2018_ge_Treg_naive Schwartzentruber_2018_ge_sensory_neuron TwinsUK_ge_blood TwinsUK_ge_fat TwinsUK_ge_LCL TwinsUK_ge_skin van_de_Bunt_2015_ge_pancreatic_islet )
for dat in "${DATASETS[@]}"; do
   echo "$dat"
   submitJobs.py --MEM 500 -j process_eqtl_catalogue_variant_info2.$dat -q normal \
      -c "$SRC/coloc/process_eqtl_catalogue_ge.variant_info.sh $COLOC/qtl_data/eqtl_catalogue $dat $COLOC/qtl_data/eqtl_catalogue"
done
grep "Successfully" FarmOut/process_eqtl_catalogue_variant_info*txt | wc -l
grep -iP "ERROR|Abort|exit code|TERM_" FarmOut/process_eqtl_catalogue_variant_info*txt | wc -l


# For the array datasets, we'd like to map probe IDs to gene IDs for downstream analysis
DATASETS=( CEDAR_B-cell_CD19 CEDAR_ileum CEDAR_monocyte_CD14 CEDAR_neutrophil_CD15 CEDAR_platelet CEDAR_rectum CEDAR_T-cell_CD4 CEDAR_T-cell_CD8 CEDAR_transverse_colon Fairfax_2012_B-cell_CD19 Fairfax_2014_monocyte_IFN24 Fairfax_2014_monocyte_LPS24 Fairfax_2014_monocyte_LPS2 )
for dat in "${DATASETS[@]}"
do
   submitJobs.py --MEM 500 -j concat_tmp.$dat -q normal \
      -c "$SRC/coloc/concat_tmp.sh $COLOC/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.*.tmp/ $dat"
done

#(while read -r -a myfile; do \
# cat $COLOC/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.*.tmp/$myfile; \
#done < <(ls $COLOC/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.*.tmp/) ) | gzip > $COLOC/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.gene_sorted2.txt.gz

#cat $COLOC/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.*.tmp/*.tmp | gzip > $COLOC/qtl_data/eqtl_catalogue/$dat.nominal.sorted.txt.gz.gene_sorted2.txt.gz

zcat $COLOC/qtl_data/eqtl_catalogue/TwinsUK_ge_blood_7M.nominal.sorted.txt.gz | awk '$2 == 1' | gzip > $COLOC/qtl_data/eqtl_catalogue/TwinsUK_ge_blood_chr1.nominal.sorted.txt.gz
time zcat $COLOC/qtl_data/eqtl_catalogue/TwinsUK_ge_blood_chr1.nominal.sorted.txt.gz | sort -k1,1 | gzip > $COLOC/qtl_data/eqtl_catalogue/TwinsUK_ge_blood_chr1.nominal.gene_sorted.txt.gz
FNAME_QTL_IN=$COLOC/qtl_data/eqtl_catalogue/TwinsUK_ge_blood_7M.nominal.sorted.txt.gz

for f in eqtl_catalogue/*.nominals.egenes_FDR_5.tsv.gz; do
  echo $f
  gzhead 3 $f
done


################################################################################
# brain eQTL meta-analysis (Solly Sieberts)
# Files downloaded from Synapse on 2019-5-21:
#   https://www.synapse.org/#!Synapse:syn16984815
# File: Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv

BRAIN_QTL_PATH=$COLOC/qtl_data/brain_meta
#md $BRAIN_QTL_PATH
cd $BRAIN_QTL_PATH

cat $BRAIN_QTL_PATH/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv | sed 's/,/\t/g' | bgzip > $BRAIN_QTL_PATH/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.tsv.gz

submitJobs.py --MEM 200 -j process_brain_meta_QTLs -q yesterday \
   -c "$SRC/coloc/process_brain_meta_qtl_file.sh $BRAIN_QTL_PATH/Cortex_MetaAnalysis.test.tsv.gz $BRAIN_QTL_PATH/Cortex_Meta.cis_eQTL.test"

submitJobs.py --MEM 1000 -j process_brain_meta_QTLs -q yesterday \
   -c "$SRC/coloc/process_brain_meta_qtl_file.sh $BRAIN_QTL_PATH/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.tsv.gz $BRAIN_QTL_PATH/Cortex_Meta.cis_eQTL"



################################################################################
# Primary microglia (Natsuhiko Kumasaka & Daniel Gaffney)

md $COLOC/qtl_data/microglia

NK5_DIR=/nfs/users/nfs_n/nk5/s117/BulkMicroglia/QTL/Rasqual/Result

# Merge microglia VCF data into one file
VCFDIR=/nfs/users/nfs_n/nk5/s117/BulkMicroglia/Snp/VCF93
(gzhead 1 $VCFDIR/chr1.gz;
 zcat $VCFDIR/chr1.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr2.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr3.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr4.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr5.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr6.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr7.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr8.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr9.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr10.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr11.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr12.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr13.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr14.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr15.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr16.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr17.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr18.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr19.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr20.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr21.maf0.as.all.vcf.gz;
 zcat $VCFDIR/chr22.maf0.as.all.vcf.gz) \
 | cut -f 1-8 | bgzip > microglia.GRCh38.vcf.gz

# Crossmap VCF file back to GRCh37 coords
CHAIN_FILE=/path/to/GRCh37_to_GRCh38.chain.gz
submitJobs.py --MEM 1000 -j CrossMap.microglia.GRCh38_to_GRCh37 -q yesterday -n 2 \
    -c "CrossMap.py vcf $CHAIN_FILE microglia.GRCh38.vcf.gz $JS/reference/GRCh37/GRCh37.75.dna.toplevel.fa.gz microglia.GRCh37.vcf"
bgzip microglia.GRCh37.vcf

# Get variant info with MAF, needed for coloc
(echo -e "chr\tpos\trsid\tMAF"; \
 zcat microglia.GRCh37.vcf.gz \
  | grep -v "^#" | cut -f 1-3,8 \
  | perl -ane 'if ($F[3] =~ /.*AF=([.0-9]+)/) { $AF = $1; print join("\t", @F[0..2], $AF)."\n" }' \
  | sort -k1,1 -k2,2n | uniq ) \
  | bgzip > microglia.variant_info.GRCh37.gz

# Select the significant eGenes using Natsuhiko's file
cp $NK5_DIR/gid.fdr5.txt .
# Number of eGenes
wc -l gid.fdr5.txt # 4569 eGenes at FDR 5%

# Copy Natsuhiko's lead SNPs file for all genes (not only significant)
(cat $ROOT/reference/rasqual_header.tsv; zcat $NK5_DIR/lead.gz) \
    | bgzip > microglia.eqtl.rasqual.leadsnps.GRCh38.tsv.gz

# Get a file with the summary stats for lead SNPs for our microglia eGenes
# The columns need to be reordered so that the first four are feature, chr, pos, rsid
(cat $ROOT/reference/rasqual_header.tsv; 
 perl $SRC/coloc/hashJoin.pl --hashFile gid.fdr5.txt --scanFile microglia.eqtl.rasqual.leadsnps.GRCh38.tsv.gz \
    --colHashFile 1 --colScanFile 1) \
    | perl -ane '{print join("\t", $F[0], $F[2], $F[3], $F[1], @F[4..24])."\n"}' > microglia.eqtl.rasqual.fdr0.05.signals.GRCh38.tsv


# Get a single file with the nominal p values per gene
cp $NK5_DIR/res*.gz tmp/
rm tmp/res917.incomp.gz
zcat tmp/res*.gz | bgzip > microglia.eqtl.rasqual.GRCh38.tsv.gz

# Subset this huge file to only those genes with significant QTLs
perl $SRC/coloc/hashJoin.pl --hashFile gid.fdr5.txt --scanFile microglia.eqtl.rasqual.GRCh38.tsv.gz \
    --colHashFile 1 --colScanFile 1 | bgzip > microglia.eqtl.rasqual.full.fdr0.05.GRCh38.tsv.gz
gzlc microglia.eqtl.rasqual.full.fdr0.05.GRCh38.tsv.gz # 12,706,488 SNP tests for eGenes
# Check how many SNPs have low R2 of optimized genotype and true genotype
zcat microglia.eqtl.rasqual.full.fdr0.05.GRCh38.tsv.gz | awk '$25 < 0.9' | bgzip > microglia.eqtl.rasqual.fdr0.05.GRCh38.r2_lt_0.9.tsv.gz
gzlc microglia.eqtl.rasqual.fdr0.05.GRCh38.r2_lt_0.9.tsv.gz # 1,726,720 SNPs with low R2
# Natsuhiko says that since the --no-posterior-update option was set, we shouldn't
# need to filter out SNPs with low R2, as they don't have inflated p values

# Convert it to our format for nominal p values files
submitJobs.py --MEM 9000 -j rasqualAddPval -q yesterday -n 2 \
  -c "bash $SRC/coloc/rasqualAddPval.microglia.sh $SRC/coloc microglia.eqtl.rasqual.full.fdr0.05.GRCh38.tsv.gz microglia.eqtl.rasqual.nominals.fdr0.05.GRCh38.tsv.gz"

# Get the lead SNP for each eGene
zcat microglia.eqtl.rasqual.nominals.fdr0.05.GRCh38.tsv.gz | 
$SRC/coloc/unique.R --args file=input.txt cols=2,3 > microglia.eqtl.rasqual.fdr0.05.signals.GRCh38.tsv

# Check that the number of genes is what we expect (same as number of sig eGenes above)
zcat microglia.eqtl.rasqual.nominals.fdr0.05.GRCh38.tsv.gz | cut -f 1 | sort | uniq | wc -l


################################################################################
# xQTL
# Files downloaded from xQTLServe on 2017-12-19:
#   http://mostafavilab.stat.ubc.ca/xQTLServe/

#md $COLOC/qtl_data/xqtl
XQTL_PATH=$COLOC/qtl_data/xQTL
cd $XQTL_PATH

wget http://mostafavilab.stat.ubc.ca/xQTLServe/eQTLs.txt
wget http://mostafavilab.stat.ubc.ca/xQTLServe/mQTLs.txt
wget http://mostafavilab.stat.ubc.ca/xQTLServe/haQTLs.txt
wget http://mostafavilab.stat.ubc.ca/xQTLServe/CellSpecifictyeQTLs.xlsx
wget http://mostafavilab.stat.ubc.ca/xQTLServe/CIT.txt
wget http://mostafavilab.stat.ubc.ca/xQTLServe/ROSMAP_SNP_INFO.txt
wget ftp://mostafavilab.stat.ubc.ca/eQTLs_all.txt.zip
wget ftp://mostafavilab.stat.ubc.ca/mQTLs_all.txt.zip
wget ftp://mostafavilab.stat.ubc.ca/haQTLs_all.txt.zip
for f in *.zip; do
    fbase="${f%.*}"
    unzip -p $f | bgzip > $fbase.gz
done

wget http://mostafavilab.stat.ubc.ca/xqtl/maf.zip
unzip -p maf.zip | bgzip > maf.gz

submitJobs.py --MEM 3000 -j process_xqtl.eQTLs -q yesterday \
   -c "$SRC/coloc/process_xqtl_file.sh $XQTL_PATH/eQTLs_all.txt.gz $XQTL_PATH/maf.gz $XQTL_PATH/xQTL_eQTL"

submitJobs.py --MEM 3000 -j process_xqtl.mQTLs -q yesterday \
   -c "$SRC/coloc/process_xqtl_file.sh $XQTL_PATH/mQTLs_all.txt.gz $XQTL_PATH/maf.gz $XQTL_PATH/xQTL_mQTL"

submitJobs.py --MEM 3000 -j process_xqtl.haQTLs -q yesterday \
   -c "$SRC/coloc/process_xqtl_file.sh $XQTL_PATH/haQTLs_all.txt.gz $XQTL_PATH/maf.gz $XQTL_PATH/xQTL_haQTL"


################################################################################
# GTEx v8 eQTLs
GTEXDIR=$ROOT/coloc/qtl_data/GTEx_v8
cd $GTEXDIR
QTL=eqtl

mkdir $GTEXDIR/eqtl
ln -s /lustre/scratch118/humgen/resources/GTEx/AnalysisV8/GTEx_Analysis_v8_EUR_eQTL_all_associations_csv eqtl_summary

# Make a file listing the names of the GTEx tissues
for f in $GTEXDIR/eqtl_summary/*.allpairs.*.csv.gz; do
    tissue=`echo $f | perl -ne 'chomp(); @A=split(/\//); @F=split(/\./, $A[$#A]); print $F[0];'`
    echo -e "$f\t$tissue" >> $GTEXDIR/gtex_tissues.eqtl_files.txt
done
cut -f 2 $GTEXDIR/gtex_tissues.files.txt | uniq > $GTEXDIR/gtex_tissues.txt

# GTEx v8 summary stats have the format: (and are in GRCh38 coords)
#,phenotype_id,variant_id,tss_distance,maf,ma_samples,ma_count,pval_nominal,slope,slope_se
#0,ENSG00000015171.19,chr10_11501_C_A_b38,-122964,0.09210526,21,21,0.7301513913589728,-0.04181598,0.12085175105293121

# Get lead SNPs for each tissue, one chr at a time.
while read -r gtex_filepath tissue; do
  echo $gtex_filepath
  filename=$(basename "$gtex_filepath")
  chr=`echo $filename | perl -ne '@x=split(/\./); print $x[4];'`
  # Important to echo tissue in to submitJobs.py here, or the loop goes wrong
  # since submitJobs.py reads from stdin
  echo $gtex_filepath | submitJobs.py --MEM 100 -j process_gtex.getLead.$filename -q normal -c "$SRC/coloc/process_gtex_file_v8.getLead.sh $gtex_filepath eqtl/$tissue.$chr"
done < gtex_tissues.eqtl_files.txt
#done < <(head -n 1 gtex_tissues.eqtl_files.txt)
grep "Successfully completed" FarmOut/process_gtex.getLead*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/process_gtex.getLead*txt > process_gtex.getLead.errors.txt


# Merge all the chr* lead SNP files for each tissue into a single file per tissue,
# and compute eGenes at FDR 5% per tissue.
while read -r tissue; do
  echo $tissue
  zcat eqtl/$tissue.chr2.qtl_signals.txt.gz | head -n 1 > eqtl/$tissue.all.qtl_signals.txt
  for filename in eqtl/$tissue.*.qtl_signals.txt.gz; do
    zcat $filename | sed '1d' >> eqtl/$tissue.all.qtl_signals.txt
  done
  
  gzip eqtl/$tissue.all.qtl_signals.txt
  # Compute FDR using the Bonferroni-corrected p values per gene.
  cp eqtl/$tissue.all.qtl_signals.txt.gz eqtl/$tissue.all.qtl_signals.txt.tmp.gz
  Rscript $SRC/addFDRCol.R eqtl/$tissue.all.qtl_signals.txt.tmp.gz 6 T | gzip > eqtl/$tissue.all.qtl_signals.txt.gz
  rm eqtl/$tissue.all.qtl_signals.txt.tmp.gz
  # Get eGenes at FDR 5%
  (gzhead 1 eqtl/$tissue.all.qtl_signals.txt.gz; zcat eqtl/$tissue.all.qtl_signals.txt.gz | awk '$7 <= .05') \
    | cut -f 1-4 > eqtl/$tissue.qtl_signals.FDR_0.05.txt
  sed '1d' eqtl/$tissue.qtl_signals.FDR_0.05.txt | cut -f 1 > eqtl/$tissue.eGenes.FDR_0.05.txt
done < gtex_tissues.txt

# Get files with nominal p values per chromosome for any eGenes at FDR 5%,
# as well as variant info for all variants tested against such eGenes.
while read -r gtex_filepath tissue; do
  echo $gtex_filepath
  echo $tissue
  filename=$(basename "$gtex_filepath")
  chr=`echo $filename | perl -ne '@x=split(/\./); print $x[4];'`
  # Important to echo tissue in to submitJobs.py here, or the loop goes wrong
  # since submitJobs.py reads from stdin
  echo $gtex_filepath | submitJobs.py --MEM 5500 -j process_gtex.vcf.$filename -q normal -c "$SRC/coloc/process_gtex_file_v8.nominals.sh $gtex_filepath eqtl/$tissue.eGenes.FDR_0.05.txt eqtl/$tissue.$chr"
done < gtex_tissues.eqtl_files.txt
grep "Successfully completed" FarmOut/process_gtex.vcf*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code|TERM_" FarmOut/process_gtex.vcf*.txt > process_gtex.vcf.errors.txt

# Merge all the chr* files of nominal p values into a single file per tissue,
# and do the same for variant info
tissue=Adipose_Subcutaneous
while read -r tissue; do
  echo $tissue
  zcat eqtl/$tissue.chr*.eGenes.FDR_0.05.nominals.txt.gz | bgzip > eqtl/$tissue.eGenes.FDR_0.05.nominals.txt.gz
  tabix -s 2 -b 3 -e 3 eqtl/$tissue.eGenes.FDR_0.05.nominals.txt.gz
  
  # Same for variant_info
  for filename in eqtl/$tissue.chr*.eGenes.FDR_0.05.variants.vcf.gz; do
    zcat $filename | sed '1d' >> eqtl/$tissue.eGenes.FDR_0.05.variants.vcf
  done
  gzip eqtl/$tissue.eGenes.FDR_0.05.variants.vcf
done < gtex_tissues.txt

# CrossMap variant_info back to GRCh37
while read -r tissue; do
   echo "$tissue"
   submitJobs.py --MEM 500 -j process_gtex.crossmap_variant_info.$tissue -q normal --nostdin \
      -c "$SRC/coloc/process_gtex.crossmap_variant_info.sh $GTEXDIR/$QTL/$tissue.eGenes.FDR_0.05.variants.vcf.gz $GTEXDIR/$QTL/$tissue.eGenes.FDR_0.05.variant_info.GRCh37.txt.gz"
done < gtex_tissues.txt
grep "Successfully" FarmOut/process_gtex.crossmap_variant_info*txt | wc -l
grep -iP "ERROR|Abort|exit code|TERM_" FarmOut/process_gtex.crossmap_variant_info*txt | wc -l


gzip eqtl/*.unmap
rm *.chr*.qtl_signals.txt.gz

