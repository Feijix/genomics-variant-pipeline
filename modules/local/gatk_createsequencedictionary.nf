process GATK_CREATESEQUENCEDICTIONARY {
    tag "gatk dict"
    container 'broadinstitute/gatk:4.6.2.0'
    storeDir 'data/reference'

    input:
    path(reference)

    output:
    path("${reference.baseName}.dict"), emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${reference} -O ${reference.baseName}.dict
    """
}
