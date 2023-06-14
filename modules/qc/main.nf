// PARAMETERS
params.CONTAINER = "theoaphidian/wgs-illumina"
params.OUTPUT = "results"

process FASTP {
    tag "${pair_id}"    
    publishDir "${params.OUTPUT}/fastp", mode: 'copy'
    
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("trimmed_${pair_id}*.fastq"), emit: trimmedReads
    path("${pair_id}_fastp.json"), emit: fastp

    script:
    def single = reads instanceof Path

    def input = !single ? "-i '${reads[0]}' -I '${reads[1]}'" : "-i '${reads}'"
    def output = !single ? "-o 'trimmed_${pair_id}_fwd.fastq' -O 'trimmed_${pair_id}_rev.fastq'" : "-o 'trimmed_${pair_id}.fastq'"

    """
    fastp ${input} ${output} --json ${pair_id}_fastp.json \\
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --cut_tail \
    --cut_mean_quality 20
    """
}

process TRIM_ENDS {
    tag "${pair_id}"    
    
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("trimmed_${pair_id}*.fastq"), emit: trimmedReads

    script:
    def single = reads instanceof Path
    def input = !single ? "-i '${reads[0]}' -I '${reads[1]}'" : "-i '${reads}'"
    def output = !single ? "-o 'trimmed_${pair_id}_fwd.fastq' -O 'trimmed_${pair_id}_rev.fastq'" : "-o 'trimmed_${pair_id}.fastq'"
    """
    fastp ${input} ${output} \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 10
    """
}

process MULTIQC {
    publishDir "${params.OUTPUT}", mode: 'copy'

    container params.CONTAINER
    
    input:
    path('fastp/*')

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc -m fastp .
    """
}
