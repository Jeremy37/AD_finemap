#!/bin/bash
OT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OT/jeremys
OUT=$JS/AD_finemap/paintor

PATH=$PATH:$JS/software/PAINTOR_V3.0

INDIR=$1
OUTDIR=$2
N_CAUSAL=$3
ANNOT=$4

if [ ! -d "$OUTDIR" ]; then
  mkdir -p $OUTDIR
fi

# Prepare input files for PAINTOR. For each locus we need the following files:
# - input file with Z scores
# - input file with annotations
# - input file with pairwise LD for all SNPs

PAINTOR -input $INDIR/paintor.loci -in $INDIR -out $OUTDIR -Zhead Z -LDname ld -enumerate $N_CAUSAL -annotations $ANNOT -ANname annotations


# Usage: PAINTOR -input [input_filename] -in [input_directory] -out [output_directory] -Zhead [Zscore_header(s)] -LDname [LD_suffix(es)]  -annotations [annot_name1,annot_name2,...]  <other options> 
#
# -input 	 (required) Filename of the input file containing the list of the fine-mapping loci [default: input.files]
# -Zhead 	 (required) The name(s) of the Zscores in the header of the locus file (comma separated) [default: N/A]
# -LDname 	 (required) Suffix(es) for LD files. Must match the order of Z-scores in locus file (comma separated) [Default:N/A]
# -annotations 	 The names of the annotations to include in model (comma separated) [default: N/A]
# -in 	 Input directory with all run files [default: ./ ]
# -out 	 Output directory where output will be written [default: ./ ]
# -Gname 	 Output Filename for enrichment estimates [default: Enrichment.Estimate]
# -Lname 	 Output Filename for log likelihood [default: Log.BayesFactor]
# -RESname 	 Suffix for ouput files of results [Default: results] 
# -ANname 	 Suffix for annotation files [Default: annotations]
# -MI 	 Maximum iterations for algorithm to run [Default: 10]
# -gamma_initial 	 inititalize the enrichment parameters to a pre-specified value (comma separated) [Default: 0,...,0]
# -num_samples  	 specify number of samples to draw for each locus [Default: 50000]
# -enumerate	 specify this flag if you want to enumerate all possible configurations followed by the max number of causal SNPs (eg. -enumerate 3 considers up to 3 causals at each locus) [Default: not specified]
# -set_seed	 specify an integer as a seed for random number generator [default: clock time at execution]
# -prop_ld_eigenvalues	 specify the proprotion of eigenvalues of LD matrix to keep when estimating the prior variance for each locus [default: 0.95]
# -mcmc	 should the algorithm be run with MCMC? [Default: not specified]
# -burn_in	 specify how many samples to discard during burn-in period [default: 50000]
# -max_samples	 specify the number of samples to keep [default: 10000]
# -num_chains	 specify the number of chains to run [default: 5]
# 
