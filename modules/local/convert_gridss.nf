/*
 * modules/local/convert_gridss.nf
 *
 * Converts GRIDSS breakend (BND) output into typed SVs (DEL/DUP/INV/INS) using
 * the StructuralVariantAnnotation R package (bundled in the GRIDSS image), then
 * emits a bgzipped, indexed, benchmark-ready VCF that matches the (meta, vcf.gz)
 * contract of the other SV callers.
 *
 * The R conversion writes contig headers, a FORMAT/GT column and a SAMPLE column
 * (Truvari requires a sample column), so the output needs no manual patching.
 * Inter-chromosomal breakends (translocations) are dropped — typed simple SVs
 * only.
 *
 * Uses the GRIDSS image (has Rscript + StructuralVariantAnnotation) for the
 * conversion, then bgzip/tabix (also in the image) for indexing.
 *
 * Container: same GRIDSS image (for R + the SVA package).
 */

process CONVERT_GRIDSS {
    tag "${meta.id}"
    container 'quay.io/biocontainers/gridss:2.13.2--h96c455f_6'

    input:
    tuple val(meta), path(bnd_vcf)
    path(fai)              // reference .fai, for contig header lines

    output:
    tuple val(meta), path("${meta.id}.gridss.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.gridss.vcf.gz.tbi"), emit: tbi

    script:
    """
    Rscript ${projectDir}/scripts/gridss_to_typed.R \\
      ${bnd_vcf} \\
      ${meta.id}.gridss.typed.vcf \\
      ${meta.id}

    # Inject contig headers from the .fai (the R script writes a minimal header),
    # then sort/bgzip/index. The R script already emits FORMAT/SAMPLE columns.
    awk '{print "##contig=<ID="\$1",length="\$2">"}' ${fai} > contigs.txt
    awk '/^#CHROM/ { while ((getline c < "contigs.txt") > 0) print c }
         { print }' ${meta.id}.gridss.typed.vcf > ${meta.id}.gridss.hdr.vcf

    bgzip -c ${meta.id}.gridss.hdr.vcf > ${meta.id}.gridss.unsorted.vcf.gz
    bcftools sort ${meta.id}.gridss.unsorted.vcf.gz -Oz -o ${meta.id}.gridss.vcf.gz 2>/dev/null \\
      || (zcat ${meta.id}.gridss.unsorted.vcf.gz | grep '^#' > s.vcf; \\
          zcat ${meta.id}.gridss.unsorted.vcf.gz | grep -v '^#' | sort -k1,1 -k2,2n >> s.vcf; \\
          bgzip -c s.vcf > ${meta.id}.gridss.vcf.gz)
    tabix -p vcf ${meta.id}.gridss.vcf.gz
    """
}
