INFILE=$1
OUTFILE=$2

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

PATH=/nfs/users/nfs_j/js29/python2.7/bin:$PATH
CrossMap.py vcf $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz $INFILE $JS/reference/GRCh37/GRCh37.75.dna.toplevel.fa.gz $INFILE.crossmap.GRCh37.txt
gzip $INFILE.crossmap.GRCh37.txt
#gzip $INFILE.crossmap.GRCh37.txt.unmap

(echo -e "chr\tpos\trsid\tMAF"; zcat $INFILE.crossmap.GRCh37.txt.gz | cut -f 1,2,3,6) | gzip > $OUTFILE
