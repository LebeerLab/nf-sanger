include { FASTP; MULTIQC } from '../../modules/qc' addParams(OUTPUT: "${params.outdir}")

params.OUTPUT = "results"

process MERGE_READS {
    errorStrategy 'ignore'

    tag "${pair_id}" 

    publishDir "${params.OUTPUT}/NGmerge", mode: 'copy'

    input:
      tuple val(pair_id), path(reads)

    output:
      tuple val(pair_id), path("*_merged.fastq"), emit: merged_fastq
      tuple val(pair_id), path("*.log"), emit: log

    when:
      reads !instanceof Path
    
    script:
    """
    NGmerge -p 0.50 -u 93 -g -v -1 ${reads[1]} -2 ${reads[0]} \
      -o ${pair_id}_merged.fastq \
      2>> ${pair_id}_merge.log
    """
}

workflow TRIM_AND_MERGE {
  take:
    reads_ch

  main:
    FASTP(reads_ch)
    FASTP.out.trimmedReads
      .ifEmpty{ error "No reads to filter" }
      .set { trimmedReads }
    FASTP.out.fastp
      .collect()
      .set{fastp}
    MULTIQC(fastp)

    MERGE_READS(trimmedReads)
  
  emit:
    merge_reads = MERGE_READS.out.merged_fastq
    merge_log = MERGE_READS.out.log  
    trimmed_reads = FASTP.out.trimmedReads
}