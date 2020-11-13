#!/usr/bin/env nextflow

// Set default parameters
params.help = false
params.input = null
params.output_folder = null

// Import the prokka module
include {prokka} from './modules/modules' params(
    output_folder: params.output_folder
)

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/AMGMA/annotate.nf --input "*.fasta.gz" --output_folder "output_folder/"

Required Arguments:
--input               File(s) to annotate (supports wildcards)
--output_folder       Folder to write GFF output file(s) into

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input == null || params.output_folder == null){
    // Invoke the function above which prints the help message
    helpMessage()

    if (params.input == null){
        log.info"""
        Please provide --input
        """.stripIndent()
    }
    if (params.output_folder == null){
        log.info"""
        Please provide --output_folder
        """.stripIndent()
    }

    // Exit out and do not run anything else
    exit 1
}

// Run Prokka on the inputs
prokka(
    Channel.fromPath(params.input)
)
