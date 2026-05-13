#!/usr/bin/env nextflow

/*
 * genomics-variant-pipeline — thin entrypoint
 *
 * Pipeline configuration lives in nextflow.config (params, profiles, resources).
 * Pipeline logic lives in workflows/variant_calling.nf and the subworkflows
 * and modules it composes.
 *
 * Run with:
 *   nextflow run main.nf -profile docker
 */

nextflow.enable.dsl=2

include { VARIANT_CALLING } from './workflows/variant_calling.nf'

workflow {
    VARIANT_CALLING()
}
