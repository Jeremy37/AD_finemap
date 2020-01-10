JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SRC=$JS/AD_finemap/src
COLOC=$JS/AD_finemap/coloc
#mkdir $COLOC/output
#mkdir $COLOC/output_nolocuscompare
cd $COLOC

GWAS_NAME=AD.meta
GWAS_FILE=$COLOC/GWAS/Alzheimers_disease_Liu_2019_UKBB_GWAX_meta_v5.sorted.txt.gz
GWAS_SIGNALS=$COLOC/AD.signals.for_coloc.txt
OUTDIR=$COLOC/output_nolocuscompare
#OUTDIR=$COLOC/output_locuscompare
PLOT_STYLE=all
PLOT_STYLE=locuscompare
PLOT_STYLE=overlay

PLOT_H4_THRESHOLD=0
PLOT_P_THRESHOLD=1
PLOT_SNP_THRESHOLD=400
GENE_ID_MAP=$JS/AD_finemap/reference/hgnc.ensembl.map.txt

GWAS_NAME=AD.meta.cond_1
GWAS_FILE=$JS/AD_finemap/gcta/output_1e-5/cond/AD.meta.cond_signals_1.tsv.gz
GWAS_SIGNALS=$JS/AD_finemap/gcta/output_1e-5/cond/AD.meta.cond.indep_signals_1.tsv

GWAS_NAME=AD.meta.cond_2
GWAS_FILE=$JS/AD_finemap/gcta/output_1e-5/cond/AD.meta.cond_signals_2.tsv.gz
GWAS_SIGNALS=$JS/AD_finemap/gcta/output_1e-5/cond/AD.meta.cond.indep_signals_2.tsv


###############################################################################
# eQTL catalogue
QTL_DIR=$COLOC/qtl_data/eqtl_catalogue

while IFS=$'\t' read -r -a DATASETS; do
  QTL_SET=${DATASETS[0]}
  QTL_N=${DATASETS[1]}
  echo $QTL_SET
  echo $QTL_N
  
  submitJobs.py --MEM 5000 -j coloc.$GWAS_NAME.$QTL_SET -q normal --nostdin \
      -c "Rscript $SRC/coloc/qtlColoc.R \
          --qtl_name $QTL_SET \
          --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
          --overlap_dist 5e5 --window 2e5 \
          --qtl_nominal $QTL_DIR/$QTL_SET.nominals.egenes_FDR_5.tsv.gz \
          --qtl_signals $QTL_DIR/$QTL_SET.qtl_signals.FDR_5.tsv  \
          --qtl_variant_info $QTL_DIR/$QTL_SET.variant_info.GRCh37.txt.gz  \
          --qtl_samplesize $QTL_N \
          --outdir $OUTDIR \
          --match_snps_by_chrpos T \
          --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
          --p1 1e-4 --p2 1e-4 --p12 1e-5 --gene_id_map $GENE_ID_MAP"
done < <(sed '1d' qtl_data/eqtl_catalogue_datasets.tsv)

grep "Successfully completed" FarmOut/coloc.$GWAS_NAME.*.txt | wc -l
grep -P "fail|error|TERM_" FarmOut/coloc.$GWAS_NAME*.txt


###############################################################################
# xQTL
QTL_DIR=$COLOC/qtl_data/xQTL

QTL=xQTL_eQTL
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $SRC/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.qtl_signals.FDR_0.05.txt  \
        --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
        --qtl_samplesize 194 \
        --outdir $OUTDIR \
        --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5 --gene_id_map $GENE_ID_MAP"

QTL=xQTL_haQTL
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $SRC/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.qtl_signals.FDR_0.05.txt  \
        --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
        --qtl_samplesize 194 \
        --outdir $OUTDIR \
        --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5 --gene_id_map $GENE_ID_MAP"

# QTL=xQTL_mQTL
# submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
#     -c "Rscript $SRC/coloc/qtlColoc.R \
#         --qtl_name $QTL \
#         --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
#         --overlap_dist 5e5 --window 2e5 \
#         --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
#         --qtl_signals $QTL_DIR/$QTL.qtl_signals.FDR_0.05.txt  \
#         --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
#         --qtl_samplesize 194 \
#         --outdir $OUTDIR \
#         --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
#         --p1 1e-4 --p2 1e-4 --p12 1e-5 --gene_id_map $GENE_ID_MAP"


###############################################################################
# Brain meta-analysis eQTL
QTL_DIR=$COLOC/qtl_data/brain_meta

QTL=Cortex_Meta.cis_eQTL
submitJobs.py --MEM 12000 -j coloc.$GWAS_NAME.$QTL -q yesterday \
    -c "Rscript $SRC/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.qtl_signals.FDR_0.05.forColoc.txt  \
        --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
        --qtl_samplesize 1433 \
        --outdir $OUTDIR \
        --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5 --gene_id_map $GENE_ID_MAP"



###############################################################################
# microglia
QTL_DIR=$COLOC/qtl_data/microglia
QTL=microglia

submitJobs.py --MEM 6000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $SRC/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/microglia.eqtl.rasqual.nominals.fdr0.05.GRCh38.tsv.gz \
        --qtl_signals $QTL_DIR/microglia.eqtl.rasqual.fdr0.05.signals.GRCh38.tsv \
        --qtl_variant_info $QTL_DIR/microglia.variant_info.GRCh37.gz \
        --qtl_samplesize 93 \
        --outdir $OUTDIR \
        --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5 --gene_id_map $GENE_ID_MAP"

# There should be 31 successful runs:
# microglia, brain_meta, 2 xQTL, 57 eQTL catalogue
grep "Successfully completed" FarmOut/coloc.$GWAS_NAME.*.txt | wc -l
grep -P "fail|error|TERM_" FarmOut/coloc.$GWAS_NAME*.txt


###############################################################################
# GTEx v8 eQTLs
QTL_DIR=$COLOC/qtl_data/GTEx_v8/eqtl
QTLTYPE=eqtl
#QTLTYPE=sqtl
#mkdir $OUTDIR/GTEx_$QTLTYPE

PLOT_H4_THRESHOLD=0.8

#tissue=Adrenal_Gland
while read -r tissue sampleSize; do
  echo $tissue
  # Important to echo tissue in to submitJobs.py here, or the loop goes wrong
  # since submitJobs.py reads from stdin
  echo $tissue | submitJobs.py --MEM 8000 -j coloc.gtex_$QTLTYPE.$tissue.$GWAS_NAME -q normal \
    -c "Rscript $SRC/coloc/qtlColoc.R \
        --qtl_name ${tissue}_$QTLTYPE \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$tissue.eGenes.FDR_0.05.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$tissue.qtl_signals.FDR_0.05.txt \
        --qtl_variant_info $QTL_DIR/$tissue.eGenes.FDR_0.05.variant_info.GRCh37.txt.gz \
        --qtl_samplesize $sampleSize \
        --outdir $OUTDIR/GTEx_$QTLTYPE \
        --plot_style $PLOT_STYLE --plot_H4_threshold $PLOT_H4_THRESHOLD --plot_snp_threshold $PLOT_SNP_THRESHOLD --plot_p_threshold $PLOT_P_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5 --match_snps_by_chrpos T --gene_id_map $GENE_ID_MAP"
done < $COLOC/qtl_data/GTEx_v8/gtex_tissues.sizes.txt
#done < <(head -n 1 $COLOC/qtl_data/GTEx_v8/gtex_tissues.sizes.txt)

grep "Successfully completed" FarmOut/coloc.gtex_$QTLTYPE.*.$GWAS_NAME.*.txt | wc -l
grep "fail|error|TERM_" FarmOut/coloc.gtex_$QTLTYPE.*.$GWAS_NAME.*.txt | wc -l

