#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false
params.min_identity = 90
params.min_coverage = 50
params.batchsize = 100

// Import the processes
include validateManifest from './modules/modules'
include fetchFTP from './modules/modules'
include renameRemoteFiles from './modules/modules'
include combineRemoteFiles from './modules/modules'
include prodigal from './modules/modules'
include combineGFF from './modules/modules'
include combineFAA from './modules/modules'
include clusterGenes as clusterGenesRound1 from './modules/modules' params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include clusterGenes as clusterGenesRound2 from './modules/modules' params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include clusterGenes as clusterGenesRound3 from './modules/modules' params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include clusterGenes as clusterGenesRound4 from './modules/modules' params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include clusterGenes as clusterGenesRound5 from './modules/modules' params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include assignCentroids from './modules/modules' params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include makeHDF from './modules/modules' params(
    output_prefix: params.output_prefix
)
include repackHDF from './modules/modules' params(
    output_folder: params.output_folder
)
include diamondDB from './modules/modules' params(
    output_prefix: params.output_prefix,
    output_folder: params.output_folder
)

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/AMGMA/build_db.nf <ARGUMENTS>

Required Arguments:
  --manifest            CSV file listing samples (see below)
  --output_folder       Folder to write output files to
  --output_prefix       Prefix to use for output file names

Optional Arguments:
  --batchsize           Number of genomes to process in a given batch (default: 100)

Manifest:
  The manifest is a CSV listing all of the genomes to be used for the database.
  The manifest much contain the column headers: uri,id,name
  The URI may start with ftp://, s3://, or even just be a path to a local file.
  The ID must be unique, and only contain a-z, A-Z, 0-9, or _.
  The NAME may be longer and contain whitespaces, but may not contain a comma.
    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.manifest == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // Parse the manifest
    Channel.from(
        file(params.manifest)
    ).splitCsv(
        header: true
    ).map {
        r -> [r["id"], r["uri"]]
    }.branch {
        ftp: it[1].startsWith("ftp://")
        other: !it[1].startsWith("ftp://")
    }.set {
        split_genome_ch
    }

    // Fetch files from FTP
    fetchFTP(
        split_genome_ch.ftp.map{r -> "${r[1]}:::${r[0]}"}.toSortedList().flatten().collate(params.batchsize)
    )

    // Make sure the other remote files are named correctly
    renameRemoteFiles(
        split_genome_ch.other.map{r -> [r[0], file(r[1])]}
    )

    // Combine the remote files
    combineRemoteFiles(
        renameRemoteFiles.out.toSortedList().flatten().collate(params.batchsize)
    )

    // Validate the manifest
    validateManifest(
        file(params.manifest)
    )

    // Now join the channels together
    genome_ch = combineRemoteFiles.out.mix(
        fetchFTP.out
    ).toSortedList().flatten()

    // Annotate the genomes
    prodigal(
        genome_ch
    )

    // Combine gene sequences
    combineFAA(
        prodigal.out[0]
    )

    // Combine gene annotation tables
    combineGFF(
        prodigal.out[1]
    )

    // Cluster genome genes by identity to find centroids
    clusterGenesRound1(
        combineFAA.out
    )

    // Round 2
    clusterGenesRound2(
        clusterGenesRound1.out.toSortedList().flatten().collate(params.batchsize)
    )

    // Round 3
    clusterGenesRound3(
        clusterGenesRound2.out.toSortedList().flatten().collate(params.batchsize)
    )

    // Round 4
    clusterGenesRound4(
        clusterGenesRound3.out.toSortedList().flatten().collate(params.batchsize)
    )

    // Round 5
    clusterGenesRound5(
        clusterGenesRound4.out.toSortedList()
    )

    // Make a DIAMOND alignment database of all genes
    diamondDB(
        clusterGenesRound5.out
    )

    // Assign a centroid to each gene in each genome
    assignCentroids(
        prodigal.out[0],
        diamondDB.out
    )

    // Make a final genome summary table
    makeHDF(
        combineGFF.out.toSortedList(),
        validateManifest.out,
        assignCentroids.out.toSortedList()
    )

    // Repack and compress the final HDF5 file
    repackHDF(
        makeHDF.out
    )
}