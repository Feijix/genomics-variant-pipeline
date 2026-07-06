/*
 * modules/local/beacon_sv2bff.nf
 *
 * Converts a typed structural-variant VCF into Beacon v2 genomicVariations BFF
 * (gzipped NDJSON), the SV counterpart of BEACON_VCF2BFF. Consumes the sv_vcfs
 * channel (per-sample, per-caller), so it converts whatever caller the routing
 * framework selected for that sample — long-read (Sniffles2/cuteSV), short-read
 * (Manta/Delly/smoove/GRIDSS), or an ensemble. No caller is hardcoded.
 *
 * Output validates against the Beacon v2 genomicVariations defaultSchema
 * (bff-validator -gv -gv-vcf). Both this and BEACON_VCF2BFF emit
 * genomicVariations records into the same MongoDB collection = unified catalog.
 *
 * meta carries id + caller + optional linking ids (biosample/individual/run/
 * analysis); these default to the sample id when not separately provided (the
 * EGA metadata adapter populates them properly later).
 */

process BEACON_SV2BFF {
    tag "${meta.id}/${meta.caller}"
    container 'python:3.11-slim'

    input:
    tuple val(meta), path(vcf_file)
    path sv_script

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.sv.genomicVariationsVcf.json.gz"), emit: bff

    script:
    def biosample  = meta.biosample  ?: meta.id
    def individual = meta.individual ?: meta.id
    def run        = meta.run        ?: "${meta.id}_run"
    def analysis   = meta.analysis   ?: "${meta.id}_${meta.caller}"
    def assembly   = params.beacon_assembly ?: 'GRCh38'
    """
    python3 ${sv_script} \\
      ${vcf_file} \\
      ${meta.id}.${meta.caller}.sv.genomicVariationsVcf.json.gz \\
      --biosample ${biosample} \\
      --individual ${individual} \\
      --run ${run} \\
      --analysis ${analysis} \\
      --dataset ${params.beacon_dataset_id} \\
      --assembly ${assembly} \\
      --caller ${meta.caller}
    """
}
