/*
 * modules/local/infer_read_length.nf
 *
 * When a sample declares no platform, classify it as short- or long-read from
 * the mean length of its first reads. Emits 'shortread' or 'longread', which
 * routing uses as a family-level key. Recovers the technology family, not the
 * specific chemistry.
 */
process INFER_READ_LENGTH {
    tag "${meta.id}"
    container 'staphb/samtools:1.24'

    input:
    tuple val(meta), path(fastq), path(bam)

    output:
    tuple val(meta), env('family'), env('mean'), emit: inferred

    script:
    def fq = fastq instanceof List ? fastq.find { it } : fastq
    """
    set +o pipefail
    if [ -s "${bam}" ]; then
        mean=\$(samtools view "${bam}" 2>/dev/null | head -10000 | awk '{s+=length(\$10); n++} END{if(n>0) print int(s/n); else print 0}')
    else
        mean=\$(zcat -f "${fq}" | head -40000 | awk 'NR%4==2{s+=length(\$0); n++} END{if(n>0) print int(s/n); else print 0}')
    fi
    family=\$(awk -v m=\$mean 'BEGIN{print (m < 500) ? "shortread" : "longread"}')
    echo "[routing] ${meta.id}: platform undeclared, mean read length \${mean} bp, inferred \${family}" >&2
    """
}
