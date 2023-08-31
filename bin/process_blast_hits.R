#!/usr/bin/env Rscript
# This script processes raw BLAST hits of 16S sequences of isolates against
# a database with reference 16S sequences. 

# dependencies: R v4.1.3, tidyverse v1.3.1 

library(tidyverse)

# input
args <- commandArgs(trailingOnly=TRUE)

# define paths
fin_hits <- args[1]
fin_reftaxonomy <- args[2]
sample_id <- str_replace(fin_hits, "_hits.txt", "")
fout <- paste0(sample_id, "_identification.tsv")

# read and parse reference genome table 
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
refgenomes <- 
  fin_reftaxonomy %>% 
  read_tsv(col_names = c("ref", "taxonomy"), col_types = cols()) %>%
  separate(taxonomy, into = ranks, sep = ";")

# read BLAST hits and select top hit per isolate
blast_columns <-
  c(
    "qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", 
    "qend", "sstart", "send", "evalue", "bitscore"
  )
sanger_sequences <-
  fin_hits %>%
  read_tsv(col_names = blast_columns, col_types = cols()) %>%
  select(query = qaccver, ref = saccver, pident, length, bitscore) %>%
  group_by(query) %>%
  arrange(desc(bitscore)) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(refgenomes, by = "ref")

# write identifications
sanger_sequences %>%
mutate(query = sample_id) %>%
arrange(query) %>% 
write.table(fout, quote = FALSE, row.names = FALSE, sep="\t")
