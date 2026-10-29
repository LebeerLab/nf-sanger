include { FASTP; MULTIQC } from '../../modules/qc'

params.outdir = "results"

process MERGE_READS {
    errorStrategy 'ignore'

    tag "${pair_id}" 

    //conda 'bioconda::ngmerge' # need to modify sourcecode to get longer sequences (>1024) merged

    publishDir "${params.outdir}/NGmerge", mode: 'copy'

    input:
      tuple val(pair_id), path(reads)

    output:
      tuple val(pair_id), path("*_merged.fastq"), emit: merged_fastq
      tuple val(pair_id), path("*.log"), emit: log

    when:
      reads !instanceof java.nio.file.Path
    
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
//      .ifEmpty{ error "No reads to filter" }
      .set { trimmedReads }
    FASTP.out.fastp
      .collect()
      .set{fastp}
    MULTIQC(fastp)

    MERGE_READS(trimmedReads)

    MERGE_READS.out.merged_fastq
  // Catch here the errors and return them as a single read
    .join(trimmedReads, remainder: true)
    .filter{ id, merged_reads, trim_reads_tp ->
      merged_reads == null  
    }.map{ [it[0], 
           (it[2][0].size() <=> it[2][1].size()) > 0 ? it[2][0] : it[2][1]]
    }
    .set{failed_to_merge}
  
  emit:
    merge_reads = MERGE_READS.out.merged_fastq
    merge_log = MERGE_READS.out.log  
    trimmed_reads = FASTP.out.trimmedReads
    failed_to_merge
}
