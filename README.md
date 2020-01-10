# Alzheimer's disease fine-mapping - supporting code
This repository contains code implementing the analyses described in the paper "Genome-wide meta-analysis, fine-mapping, and integrative prioritization identifies new Alzheimer’s disease risk genes".

All key analysis steps are implemented or outlined in the file src/AD_finemap_master.sh. These steps are meant to be run on the command line (or jobs submitted to a compute cluster).

Most of the relevant input files are provided, with the exception of some large files that either require approved access (e.g. UK Biobank genotypes) or that can be downloaded elsewhere (e.g. dbSNP, Roadmap Epigenomics).
