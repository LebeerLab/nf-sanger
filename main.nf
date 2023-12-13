#!/usr/bin/env nextflow
 
params.sample_dir = "$baseDir/data/"
params.outdir = "results"
params.blastdb = "../reference_databases/bac120_ssu_reps_r95"
params.taxonomy = "../reference_database/bac120_taxonomy_r95.tsv.gz"
params.type = "16S"
params.help = false

//===== INCLUDE MODULES ==========================================
include { TRIM_ENDS } from './modules/qc' addParams(OUTPUT: "${params.outdir}")
include { TRIM_AND_MERGE; TRIM_AND_MERGE as TRIM_AND_MERGE_REVERSED_READS} \
  from './modules/merge_reads' addParams(OUTPUT: "${params.outdir}")
// ================================================================

def helpMessage() {
    log.info"""
     Name: nf-sanger-sequencing
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
        nextflow run main.nf --sample_dir '/path/to/data'
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

process BASECALL {
    conda 'bioconda::tracy'
    tag "${pair_id}" 

    publishDir "${params.outdir}/tracy", mode: 'copy', pattern: "*.log"

    input:
      tuple val(pair_id), path(reads)

    output:
      tuple val(pair_id), path("*.fastq"), emit: fastq
      tuple val(pair_id), path("*.log"), emit: log


    when:
      reads !instanceof Path

    script:
    def is_zipped = reads[0].getName().endsWith('.gz')
    def rev_read = is_zipped ? "${reads[0].baseName}" : reads[0]
    def fw_read = is_zipped ? "${reads[1].baseName}" : reads[1]
    //tracy basecall -f fastq -o out.fastq ${reads} 2>&1 >> read.log
    """
    if [[ $is_zipped == true ]]
    then
    gunzip -f ${reads[0]}
    gunzip -f ${reads[1]}
    fi
    tracy basecall -f fastq -o rev_out.fastq $rev_read 2>&1 >> rev_read_"${pair_id}".log
    tracy basecall -f fastq -o fw_out.fastq $fw_read 2>&1 >> fw_read_"${pair_id}".log
    """
}

process REVERSE_READS {
    tag "${pair_id}-failed-merge" 

    input:
      tuple val(pair_id), path(reads)

    output:
      tuple val(pair_id), path("*_reversed.fastq"), emit: reversed_reads
  
    script:
    def single = reads instanceof Path
    def out_file_fw = single ? "${reads.baseName}_reversed.fastq" : "${reads[1].baseName}_reversed.fastq"
    def out_file_rv = single ? null : "${reads[0].baseName}_reversed.fastq"
    """
    if [ $single == true ]; then
    ./reverse_fastq.sh $reads $out_file_fw

    else
    reverse_fastq.sh ${reads[1]} $out_file_fw
    reverse_fastq.sh ${reads[0]} $out_file_rv
    fi
    """
}

process FASTQ_TO_FASTA {
  tag "${pair_id}"

  //publishDir "${params.outdir}/fasta", mode: 'copy'
  conda 'bioconda::emboss'

  input:
    tuple val(pair_id), path(fastq)

  output:
    tuple val(pair_id), path("*.fasta"), emit: fasta

  script:
  def outf = "${fastq.baseName}.fasta"
  """
  seqret -sequence $fastq \
    -outseq $outf
  """
}

process ALIGN_READS {
  tag "${pair_id}" 
  publishDir "${params.outdir}/mafft", mode: 'copy'

  conda "${projectDir}/conf/align.yaml"

  input:
    tuple val(pair_id), path(reads)

  output:
    tuple val(pair_id), path("*.aln"), emit: alignment

  script:
  """
  seqret -sequence ${reads[0]} forward.fasta
  revseq -sequence ${reads[1]} reverse_complement.fasta
  cat forward.fasta reverse_complement.fasta > both.fasta
  mafft both.fasta > ${pair_id}.aln 2>/dev/null
  """  
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

  abf = Channel.fromFilePairs("${params.sample_dir}/**/*_{*F,*R}_*.ab1")
    .map{ group_key, files -> tuple(
        "${files[0].parent.baseName}_$group_key", files
    )}
    //.take(30)
  ab_files_gz = Channel.fromFilePairs("${params.sample_dir}/**/*_{*F,*R}_*.ab1.gz")
    .map{ group_key, files -> tuple(
        "${files[0].parent.baseName}_$group_key", files
    )}

  ab_files = abf.concat(ab_files_gz)

  BASECALL(ab_files)
  TRIM_AND_MERGE(BASECALL.out.fastq)
  TRIM_AND_MERGE
     .out.merge_log
     .branch{
         succeeded: it[1].splitText()[-1] =~ /stitched: 1/ 
         failed: it[1].splitText()[-1] =~ /stitched: 0/ 
     }.set{ merging_process } 

// RETRY MERGE AFTER REVERSING READS
  if (merging_process.failed) {
    BASECALL.out.fastq
     .join(merging_process.failed)
     .map{ [ it[0] , it[1] ]}
     .set{ reads_failed_to_merge }

     REVERSE_READS(reads_failed_to_merge)
     TRIM_AND_MERGE_REVERSED_READS(REVERSE_READS.out.reversed_reads)
     TRIM_AND_MERGE_REVERSED_READS
     .out.merge_log
     .branch{
         succeeded: it[1].splitText()[-1] =~ /stitched: 1/ 
         failed: it[1].splitText()[-1] =~ /stitched: 0/ 
      }.set{ second_merge_process }

      if (second_merge_process.failed) {
      // Take the longest read after trimming of the samples that failed
       TRIM_AND_MERGE.out.trimmed_reads
       .join(second_merge_process.failed)
       // select here the longest read in the it[1] tuple
       .map{ [it[0], 
              (it[1][0].size() <=> it[1][1].size()) > 0 ? it[1][0] : it[1][1]]
        }
       .set{ trimmed_single_reads }
      }
  }

  TRIM_AND_MERGE.out.merge_reads
  .join(merging_process.succeeded)
  .map{ [it[0], it[1]]}  
  .set{ merges_forward }


  TRIM_AND_MERGE_REVERSED_READS.out.merge_reads
  .join(second_merge_process.succeeded)
  .map{ [it[0], it[1]]}  
  .set{ merges_reverse }
  
  // Use forward reads of the samples that raise an error
  TRIM_AND_MERGE.out.failed_to_merge.
   filter{ !it[1].isEmpty() }
   .set{merge_failures}

  merges = merges_forward.concat(merges_reverse)
  TRIM_ENDS(merges)
  TRIM_ENDS.out.trimmedReads.concat(trimmed_single_reads, merge_failures)
    .set{ all_fastqs }
  
  FASTQ_TO_FASTA(all_fastqs)
  // Alignment trimmed reads
  ALIGN_READS(
   TRIM_AND_MERGE.out.trimmed_reads
     .filter{ !it[1].isEmpty() }
     .filter{ !it[1][0].isEmpty() }
     .filter{ !it[1][1].isEmpty() }
  )
  BLAST(
   FASTQ_TO_FASTA.out.fasta
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
