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
    for sub in view norm concat annotate index; do
        printf '#!/bin/sh\\nexec /opt/hap.py/bin/bcftools %s "\$@"\\n' "\$sub" > shims/\$sub
        chmod +x shims/\$sub
    done
    export PATH="\$PWD/shims:/opt/hap.py/bin:\$PATH"

    # som.py's allele conversion can move a record ahead of its neighbour, and the
    # bundled bcftools 1.4.1 cannot sort. Sort the query before som.py indexes it:
    # patch the index call in a copy of hap.py's Python library (the image is
    # read-only for the task user) laid out so that Tools finds ../../../bin.
    cat > shims/sortvcf << 'EOS'
#!/bin/sh
f="\$1"
( /opt/hap.py/bin/bcftools view -h "\$f"
  /opt/hap.py/bin/bcftools view -H "\$f" | /usr/bin/sort -k1,1 -k2,2n ) \\
  | /opt/hap.py/bin/bgzip -c > "\$f.s" && mv "\$f.s" "\$f"
EOS
    chmod +x shims/sortvcf

    mkdir -p hap/lib
    cp -r /opt/hap.py/lib/python27 hap/lib/
    ln -s /opt/hap.py/bin hap/bin
    python - "\$PWD/shims/sortvcf" << 'EOP'
import sys
p = 'hap/lib/python27/Tools/bcftools.py'
s = open(p).read()
old = 'runShellCommand(\\'bcftools\\', "index", "-t", output_filename)'
new = '(runShellCommand("%s", output_filename), %s)' % (sys.argv[1], old)
assert old in s, 'index call not found'
open(p, 'w').write(s.replace(old, new))
EOP
    export PYTHONPATH="\$PWD/hap/lib/python27:\${PYTHONPATH:-}"

    /opt/hap.py/bin/som.py \\
      ${truth_vcf} \\
      ${vcf} \\
      -R ${truth_bed} \\
      -r ${reference} \\
      -o ${meta.id}.${meta.caller}.sompy
    """
}
