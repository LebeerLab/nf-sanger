# nf-sanger-sequencing
Nextflow pipeline for sanger data

## Usage
Running the pipeline with a `--help` flag produces the following output:
```
 N E X T F L O W  ~  version 22.10.3
Launching `main.nf` [marvelous_brattain] DSL2 - revision: 1e5c93dbfe

 Name: nf-sanger-sequencing
 Author: LAMB (UAntwerp)
=========================================
Required arguments:
  --sample_dir          Path to the samplesheet containing the reads to be assembled. Default: data
  --outdir              Column with the name of the sample. Default: results
  --taxonomy            Path to the gtdb taxonomy file.
                            Default: /mnt/b/seqdata/sanger/reference_databases/bac120_taxonomy_r95.tsv.gz
  --blastdb             Path to the blast database used for classifying reads.
                            Default: /mnt/b/seqdata/sanger/reference_databases/bac120_ssu_reps_r95

Usage example:
    nextflow run main.nf --sample_dir '/path/to/data'
```

The four params needed for the pipeline (`sample_dir`, `outdir`, `blastdb` and `taxonomy`) are set up to sensible defaults for lambda in the `nextflow.config` file.
These can be adapted to the needs of the pipeline run via the command line (see usage example) or in the config file prior to running.

## Overview pipeline strategy
The pipeline will try to merge forward and reverse reads and do classification on the resulting merged reads. 
Failing that a first time, it will try a merge of the reverse complements of the reads.
Should this merge not succeed, the longest read of the pair will be kept for classfication, while the other is discarded.

From the resulting blast output, the most frequent occuring species in the top 5 bitscore hits gets assigned as the eventual result.
In case of discrepancies, the user should consult the full output table.
