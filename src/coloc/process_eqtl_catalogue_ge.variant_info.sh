QTLDIR_IN=$1
DATASET=$2
OUT_ROOT=$3

QTL_BASE=$OUT_ROOT/$DATASET
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

CrossMap.py vcf $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz $QTLDIR_IN/$DATASET.variant_information.txt.gz $JS/reference/GRCh37/GRCh37.75.dna.toplevel.fa.gz $QTL_BASE.variant_information.GRCh37.txt
#gzip $QTL_BASE.variant_information.GRCh37.txt
#gzip $QTL_BASE.variant_information.GRCh37.txt.unmap

(echo -e "chr\tpos\trsid\tMAF"; zcat $QTL_BASE.variant_information.GRCh37.txt | cut -f 1,2,3,9) | gzip > $QTL_BASE.variant_info.GRCh37.txt.gz
