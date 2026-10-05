process SAMTOOLS_FAIDX {
    tag "samtools faidx"
    container 'staphb/samtools:1.24'
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
