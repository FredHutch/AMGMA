#!/usr/bin/env nextflow

/*
  Utility to process a set of geneshot results
  for visualization in the GLAM Browser.
*/

// Using DSL-2
nextflow.preview.dsl=2

// Default values for flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.input = false
params.details = false
params.output_prefix = false
params.help = false
params.aws_region = "us-east-1"
params.branch = 'latest'

// Set the containers to user
container__glam = "quay.io/fhcrc-microbiome/glam-browser-v2:${params.branch}"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/glam-browser <ARGUMENTS>
    
    Arguments:
      --input               Geneshot results (".hdf5") file to process
      --details             Geneshot details (".hdf5") file to process
      --output_prefix       Location in AWS S3 to write out indexed file objects

    Optional arguments:
      --aws_region          AWS Region for the output S3 bucket (default: us-east-1)
      --branch              Branch of the `glam-browser-v2` repo to use for indexing (default: latest)

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input == false || params.output_prefix == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make CAGs for each set of samples, with the subset of genes for this shard
process indexGeneshotResults {
    container "${container__glam}"
    label "mem_veryhigh"
    errorStrategy 'retry'
    publishDir "${params.output_prefix}", mode: 'copy', overwrite: true

    input:
    path summary_hdf
    path details_hdf

    output:
    file "OUTPUT/**"

    """#!/bin/bash

set -Eeuxo pipefail

AWS_REGION=${params.aws_region} \
glam-cli index-dataset --fp "${summary_hdf}" --uri "\$PWD/OUTPUT" --details "${details_hdf}"

    """
}

workflow {

  // Index the input file
  indexGeneshotResults(
    Channel.fromPath(params.input),
    Channel.fromPath(params.details),
  )

}
