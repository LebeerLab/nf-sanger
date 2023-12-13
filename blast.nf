#!/usr/bin/env nextflow
 
params.sample_dir = "$baseDir/data/"
params.outdir = "results"
params.blastdb = "../reference_databases/bac120_ssu_reps_r95"
params.taxonomy = "../reference_database/bac120_taxonomy_r95.tsv.gz"
params.type = "16S"
params.help = false


def helpMessage() {
    log.info"""
     Name: nf-sanger-sequencing-blast
     Author: LAMB (UAntwerp)
    =========================================
    Required arguments:
      --sample_dir          Path to the samplesheet containing the reads to be assembled. Default: ${params.sample_dir}
      --outdir              Column with the name of the sample. Default: ${params.outdir}
      --taxonomy            Path to the gtdb taxonomy file. 
                                Default: ${params.taxonomy}
      --blastdb             Path to the blast database used for classifying reads. 
                                Default: ${params.blastdb}

    Usage example:
        nextflow run blast.nf --sample_dir '/path/to/data'
    """.stripIndent()
}

def paramsUsed() {
    log.info"""
    N F - W G S - I L L U M I N A
    =========================================
    sample_dir:      ${params.sample_dir}
    outdir:       ${params.outdir}
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


process BLAST {
  tag "${pair_id}" 
  publishDir "${params.outdir}/${params.type}_blast", mode: 'copy', pattern:'*_hits.txt'
  publishDir "${params.outdir}/${params.type}/sequences", mode: 'copy', pattern:'*.fasta'

  conda 'bioconda::blast'

  input:
    tuple val(pair_id), path(fasta)

  output:
    tuple val(pair_id), path("*_hits.txt"), emit: hits
    tuple val(pair_id), path(fasta)

  script:
  """
  blastn -query $fasta \
  -db $params.blastdb \
  -outfmt 6 \
  -num_threads $task.cpus > ${pair_id}_hits.txt
  """
}

process PROCESS_BLAST_HITS {
  tag "${pair_id}" 
  //publishDir "${params.outdir}/16S_blast", mode: 'copy'

  input:
    tuple val(pair_id), path(hits)

  output:
    path("*_identification.tsv"), emit: table

  script:
  """
  process_blast_hits.R $hits $params.taxonomy
  """
}

workflow {

  paramsUsed()

  reads = Channel.fromPath(params.sample_dir + "/*.fasta")
    .map{ it -> tuple(it.baseName, it)}

  
  BLAST(
   reads
  )

  PROCESS_BLAST_HITS(
    BLAST.out.hits
  )

  PROCESS_BLAST_HITS
   .out.table
   .collectFile(
    name: "classification_results.tsv", 
    storeDir: "${params.outdir}", 
    keepHeader: true,
    )

}
