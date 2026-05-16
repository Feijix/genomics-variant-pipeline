process SAMTOOLS_FAIDX {
    tag "samtools faidx"
    container 'staphb/samtools:latest'
    storeDir 'data/reference'

    input:
    path(reference)

    output:
    path("${reference}.fai"), emit: fai

    script:
    """
    samtools faidx ${reference}
    """
}
