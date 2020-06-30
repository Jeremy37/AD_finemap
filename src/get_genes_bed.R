#!/usr/bin/env Rscript
library(tidyverse)
library(annotables)

genes.df = grch37 %>%
  filter(biotype == "protein_coding") %>%
  group_by(ensgene) %>%
  mutate(start = max(1, min(start - 10000)), # get outermost coordinates of duplicate ensgene rows
         end = max(end + 10000)) %>%
  filter(!duplicated(ensgene)) %>%
  select(chr, start, end, ensgene) %>%
  mutate(chr = factor(as.character(chr), levels = as.character(c(1:22, "X", "Y")))) %>%
  na.omit() %>%
  arrange(chr, start) %>%
  write_tsv("expression/genes.pc.10kb_window.grch37.bed", col_names = F)

  