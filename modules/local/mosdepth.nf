// Coverage depth/breadth. Region-restricted to match the calling region so the
// numbers describe what was actually analysed, not the whole reference.
process MOSDEPTH {
    tag "${meta.id}"
    container 'quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0'
    publishDir "${params.outdir}/qc/mosdepth", mode: 'copy'
    cpus 4

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.mosdepth.*"), emit: summary

    script:
    def region = params.target_region ? "-c ${params.target_region}" : ''
    """
    mosdepth -t ${task.cpus} -n ${region} ${meta.id} ${bam}
    """
}
