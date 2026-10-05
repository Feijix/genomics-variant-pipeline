/*
 * modules/local/beacon_vcf2bff.nf
 *
 * Converts an annotated VCF into Beacon v2 Friendly Format (BFF) —
 * specifically the genomicVariations entity (genomicVariationsVcf.json.gz).
 *
 * Wraps beacon2-cbi-tools `bff-tools vcf`. The tool re-annotates with
 * snpEff + dbNSFP + ClinVar + COSMIC on top of the input VEP annotations,
 * then serialises to the Beacon v2 genomicVariations JSON shape.
 *
 * DEPENDENCY (manual, one-time): the reference-data bundle's snpEff config
 *   <refdata>/soft/snpEff/snpEff.config
 * must have `data.dir` pointing at the container mount path:
 *   data.dir = /beacon2-cbi-tools-data/databases/snpeff/v5.0/
 * Without this, snpEff attempts a network download at runtime and fails.
 *
 * The reference-data bundle is mounted at /beacon2-cbi-tools-data via
 * containerOptions. params.beacon_refdata must point at the host bundle.
 */

process BEACON_VCF2BFF {
    tag "${meta.id}/${meta.caller}"
    container 'manuelrueda/beacon2-cbi-tools:2.0.12'
    containerOptions "-v ${params.beacon_refdata}:/beacon2-cbi-tools-data"

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.genomicVariationsVcf.json.gz"), emit: bff

    script:
    """
    # Generate the per-sample param file. datasetid override is the only
    # metadata field set at this stage (Q3-lazy); genome is mandatory.
    cat > param.yaml <<EOF
    genome: ${params.genome}
    datasetid: ${params.beacon_dataset_id}
    EOF

    # bff-tools refuses to write into an existing project dir, so use a
    # caller-unique name. Nextflow's per-task work dir is already isolated,
    # but be explicit to avoid collisions if multiple callers share a dir.
    /usr/share/beacon2-cbi-tools/bin/bff-tools vcf -t ${task.cpus} \\
      -i ${vcf_file} \\
      -p param.yaml \\
      --projectdir-override beacon_${meta.id}_${meta.caller}

    # Surface the genomicVariations BFF with a stable, caller-tagged name.
    mv beacon_${meta.id}_${meta.caller}/vcf/genomicVariationsVcf.json.gz \\
       ${meta.id}.${meta.caller}.genomicVariationsVcf.json.gz
    """
}
