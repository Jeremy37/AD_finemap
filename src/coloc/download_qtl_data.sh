QTL_DIR=$JS/datasets/qtl

cd $QTL_DIR

lftp ftp-private.ebi.ac.uk
lftp ftp-private.ebi.ac.uk:~> login otar_sumstats
Password: W9SRZ9ik

cd /upload/eqtls/v0.1/summary_stats

mkdir Alasoo_2018
mkdir BLUEPRINT_PE
mkdir BLUEPRINT_SE
mkdir BrainSeq
mkdir CEDAR
mkdir Fairfax_2012
mkdir Fairfax_2014
mkdir GENCORD
mkdir GEUVADIS
mkdir HipSci
mkdir Kasela_2017
mkdir Lepik_2017
mkdir Naranbhai_2015
mkdir Nedelec_2016
mkdir Quach_2016
mkdir ROSMAP
mkdir Schmiedel_2018
mkdir Schwartzentruber_2018
mkdir TwinsUK
mkdir van_de_Bunt_2015
DATASETS=( Alasoo_2018 BLUEPRINT_PE BLUEPRINT_SE BrainSeq CEDAR Fairfax_2012 Fairfax_2014 GENCORD GEUVADIS HipSci Kasela_2017 Lepik_2017 Naranbhai_2015 Nedelec_2016 Quach_2016 ROSMAP Schmiedel_2018 Schwartzentruber_2018 TwinsUK van_de_Bunt_2015 )
for DS in "${DATASETS[@]}"
do
   mkdir "$DS"
done


mget -O Alasoo_2018/ Alasoo_2018/Alasoo_2018_ge_*
mget -O BLUEPRINT_PE/ BLUEPRINT_PE/BLUEPRINT_PE_ge_*
mget -O BLUEPRINT_SE/ BLUEPRINT_SE/BLUEPRINT_SE_ge_*
mget -O BrainSeq/ BrainSeq/BrainSeq_ge_*
mget -O GENCORD/ GENCORD/GENCORD_ge_*
mget -O GEUVADIS/ GEUVADIS/GEUVADIS_ge_*
mget -O HipSci/ HipSci/HipSci_ge_*
mget -O Lepik_2017/ Lepik_2017/Lepik_2017_ge_*
mget -O Nedelec_2016/ Nedelec_2016/Nedelec_2016_ge_*
mget -O Quach_2016/ Quach_2016/Quach_2016_ge_*
mget -O ROSMAP/ ROSMAP/ROSMAP_ge_*
mget -O Schmiedel_2018/ Schmiedel_2018/Schmiedel_2018_ge_*
mget -O Schwartzentruber_2018/ Schwartzentruber_2018/Schwartzentruber_2018_ge_*
mget -O TwinsUK/ TwinsUK/TwinsUK_ge_*
mget -O van_de_Bunt_2015/ van_de_Bunt_2015/van_de_Bunt_2015_ge_*


mget -O CEDAR/ CEDAR/CEDAR_*
mget -O Fairfax_2012/ Fairfax_2012/Fairfax_2012_*
mget -O Fairfax_2014/ Fairfax_2014/Fairfax_2014_*
mget -O Kasela_2017/ Kasela_2017/Kasela_2017_*
mget -O Naranbhai_2015/ Naranbhai_2015/Naranbhai_2015_*


# Symlink to all QTL files from the eQTL catalog in one directory
md $JS/datasets/eqtl_catalogue
for f in $QTL_DIR/*/*; do
  ln -s $f
done



