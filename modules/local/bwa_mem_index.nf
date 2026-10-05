process BWA_MEM_INDEX {
    tag "bwa-mem2 index"
    container 'nanozoo/bwa-mem2:2.3--a42e075'
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
