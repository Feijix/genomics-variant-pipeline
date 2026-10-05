#!/usr/bin/env nextflow

/*
 * genomics-variant-pipeline - entrypoint
 *
 * Run with:
 *   nextflow run main.nf -profile docker
 */

nextflow.enable.dsl=2

include { VARIANT_CALLING } from './workflows/variant_calling.nf'

workflow {
    VARIANT_CALLING()
}
