process indexReference {
    tag "bwa-mem2 index"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    storeDir 'data/reference'

    input:
    path(reference)

    output:
    path("${reference}.*"), emit: index

    script:
    """
    bwa-mem2 index ${reference}
    """
}
