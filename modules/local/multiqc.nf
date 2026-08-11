// Aggregates every QC artefact into one report. Takes a flat list of files.
process MULTIQC {
    container 'multiqc/multiqc:v1.35'
    publishDir "${params.outdir}/qc", mode: 'copy'
    memory '8 GB'

    input:
    path('qc_files/*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    """
    multiqc qc_files/ --filename multiqc_report.html --force
    """
}
