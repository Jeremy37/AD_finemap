#!/usr/bin/env Rscript
# Gathers expression data from multiple datasets to create a large table of
# median TPMs per gene & dataset.
library(scrattch.io)
library(tidyverse)
library(Matrix)

tome = "transcrip.tome"

sample_names = read_tome_sample_names(tome)
#write_tsv(enframe(sample_names, name=NULL) %>% rename(sample_id = value), "sample_names.tsv")

gene_names = read_tome_gene_names(tome)
#write_tsv(enframe(gene_names, name=NULL) %>% rename(gene_id = value), "gene_ids.tsv")

exons = read_tome_dgCMatrix(tome, "data/t_exon")
dim(exons)
class(exons)
rownames(exons) = gene_names
colnames(exons) = sample_names
#write_tsv(as_tibble(as.matrix(exons[1:100,])), "exons.1-100.tsv.gz")
# Data too large to convert to matrix
#write_tsv(as.matrix(exons), "exons.tsv.gz")
# For some reason this call doesn't work - I'm not sure why
#write_dgCMatrix_csv(exons, "exons.csv")

# Get metadata of cells
sample_meta = read_csv("sample_annotations.csv")
subclass_counts = table(sample_meta$subclass_label)
print(subclass_counts)

sample_meta.summary = sample_meta %>%
  group_by(cluster_label) %>%
  summarise(cluster_count = n(),
            subclass_label = first(subclass_label),
            class_label = first(class_label)) %>%
  group_by(subclass_label) %>%
  mutate(subclass_count = sum(cluster_count))
write_tsv(sample_meta.summary, path = "cluster_counts.tsv")
write_tsv(sample_meta.summary %>% summarise(subclass_count = first(subclass_count)) %>% arrange(desc(subclass_count)),
          path = "subclass_counts.tsv")

celltype_names = levels(factor(sample_meta$subclass_label))

# We want to get the average expression for each gene within cells grouped
# by cell type. To do this, we create a model matrix that defines columns
# with 1/0 indicating which cell is within each group. We can then use
# matrix multiplication to get sums within groups. I based how to do this
# on this post: https://slowkow.com/notes/sparse-matrix/
celltype_mat = model.matrix(~ 0 + factor(sample_meta$subclass_label))
colnames(celltype_mat) <- celltype_names
#head(celltype_mat)

gene_sums = as.matrix(exons %*% celltype_mat)
gene_sums = bind_cols(hgnc_symbol = gene_names, as_tibble(gene_sums))
write_tsv(gene_sums, "brain_cell_type_sums.tsv.gz")

# Convert these cell type counts to TPMs. For this we need gene lengths.
# We use the same genes version used in the scRNA-seq study.
# Downloaded from http://may2015.archive.ensembl.org/biomart:
#(echo -e "gene_id\ttranscript_id\tchr\ttranscript_length\thgnc_symbol"; \
#  sed '1d' mart_export.txt | sort -k1,1) | gzip > GRCh38.p2.geneid_hgnc_length.txt.gz
gene_lengths.df = read_tsv("../reference/GRCh38.p2.geneid_hgnc_length.txt.gz") %>%
  arrange(chr, desc(transcript_length)) %>%
  filter(!is.na(transcript_length), !is.na(hgnc_symbol), !duplicated(hgnc_symbol))

# ...But this lacks a few updated gene symbols, so I supplement this with
# those from a later release.
#Downloaded from https://www.ensembl.org/biomart/ on 2020-5-12.
#(echo -e "gene_id\ttranscript_id\tchr\ttranscript_length\thgnc_symbol"; \
#  sed '1d' mart_export.txt | sort -k1,1) | gzip > GRCh38.p13.geneid_hgnc_length.txt.gz
gene_lengths.new.df = read_tsv("../reference/GRCh38.p13.geneid_hgnc_length.txt.gz") %>%
  arrange(chr, desc(transcript_length)) %>%
  filter(!is.na(transcript_length), !is.na(hgnc_symbol), !duplicated(hgnc_symbol)) %>%
  filter(!hgnc_symbol %in% gene_lengths.df$hgnc_symbol)

gene_lengths.df = bind_rows(gene_lengths.df, gene_lengths.new.df)

# Subset to genes in our dataset. It's quite annoying that we have to use HGNC symbols
# to join the lengths to the scRNA-seq counts, since HGNC symbols aren't stable. But
# it's all we've got for the count data.
gene_sums.flt = gene_sums %>%
  inner_join(gene_lengths.df %>% select(gene_id, hgnc_symbol, transcript_length), by="hgnc_symbol") %>%
  select(gene_id, hgnc_symbol, transcript_length, everything())

getTPMs = function(counts, gene_lengths) {
  rpK = apply(counts, MARGIN=2, FUN=function(x) x / gene_lengths) * 1e3
  tpm.mat = apply(rpK, MARGIN=2, FUN=function(x) x / sum(x)) * 1e6
  tpm.mat
}

gene_tpm = getTPMs(as.matrix(gene_sums.flt %>% select(-gene_id, -hgnc_symbol, -transcript_length)), gene_sums.flt$transcript_length)
gene_tpm.df = bind_cols(gene_sums.flt[, c("gene_id", "hgnc_symbol")], as_tibble(gene_tpm))
write_tsv(gene_tpm.df, "brain_cell_type_tpm.tsv.gz")


