process BWA_MEM_INDEX {
    tag "bwa-mem2 index"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    storeDir 'data/reference'

    input:
    path(reference)

    output:
    path("${reference}.{0123,amb,ann,bwt.2bit.64,pac}"), emit: index

    script:
    """
    bwa-mem2 index ${reference}
    """
}
