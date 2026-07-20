process SOMPY_SOMPY {
    tag "${meta.id}.${meta.caller}"
    container 'jmcdani20/hap.py:v0.3.12'   // same image as HAPPY_HAPPY — match your pin
    publishDir "${params.outdir}/benchmark/somatic", mode: 'copy'
    cpus 2
    memory '8 GB'

    input:
    tuple val(meta), path(vcf), path(truth_vcf), path(truth_tbi), path(truth_bed)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.sompy.stats.csv"), emit: stats
    path("${meta.id}.${meta.caller}.sompy.*"),                          emit: all

    script:
    """
    export HGREF=${reference}

    # som.py's bcftools wrapper calls 'view'/'norm'/etc. as bare commands,
    # not 'bcftools view'. Create shims so they resolve.
    mkdir -p shims
    for sub in view norm concat sort index annotate; do
        printf '#!/bin/sh\\nexec /opt/hap.py/bin/bcftools %s "\$@"\\n' "\$sub" > shims/\$sub
        chmod +x shims/\$sub
    done
    export PATH="\$PWD/shims:/opt/hap.py/bin:\$PATH"

    /opt/hap.py/bin/som.py \\
      ${truth_vcf} \\
      ${vcf} \\
      -R ${truth_bed} \\
      -r ${reference} \\
      -o ${meta.id}.${meta.caller}.sompy
    """
}