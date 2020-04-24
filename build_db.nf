#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false
params.batchsize = 100


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
    )

    // Reformat each batch of genomes
    combineGenomes(
        genome_ch
    )

    // Make a single tarball
    makeDatabase(
        combineGenomes.out.toSortedList(),
        validateManifest.out
    )

}

// Validate that all genomes are unique in the manifest
process validateManifest {
    tag "Enforce unique genome IDs"
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label 'io_limited'
    errorStrategy 'retry'

    input:
    file manifest

    output:
    file "${manifest}"

"""
#!/usr/bin/env python3

import pandas as pd

# Read in the CSV
manifest = pd.read_csv("${manifest}")

for k in ['uri', 'id', 'name']:
    assert k in manifest, "Please provide %s as a column in the manifest" % k

# Make sure that all genome IDs are unique
all_unique = True
for n, v in manifest['uri'].value_counts().items():
    if v > 1:
        all_unique = False
        print("%s is found %d times in the manifest" % (n, v))
assert all_unique, "Must provide entirely unique genome IDs"

"""
}


// Fetch genomes via FTP
process fetchFTP {
    tag "Download genomes hosted by FTP"
    container 'quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f'
    label 'io_limited'
    errorStrategy "retry"

    input:
        val uri_id_list
    
    output:
        file "genomes_fasta.*.tar"
    
"""
#!/bin/bash
set -e

for uri_id in ${uri_id_list.join(" ")}; do

    uri=\$(echo \$uri_id | sed 's/:::.*//')
    id=\$(echo \$uri_id | sed 's/.*::://')

    echo "Downloading \$id from \$uri"

    wget --quiet -O \$id.fasta.gz \$uri

    # Make sure the file is gzip compressed
    (gzip -t \$id.fasta.gz && echo "\$id.fasta.gz is in gzip format") || ( echo "\$id.fasta.gz is NOT in gzip format" && exit 1 )

done

echo "Making a tar with all genomes in this batch"
tar cfh \$(mktemp genomes_fasta.XXXXXXXXX).tar *.fasta.gz

echo "done"

"""
}


// Rename the files from S3 to explicitly match the provided ID
process renameRemoteFiles {
    container 'ubuntu:20.04'
    label 'io_limited'
    errorStrategy 'retry'

    input:
        tuple val(id), file(genome_fasta)

    output:
        file "${id}.fasta.gz"

"""
#!/bin/bash

set -e

ls -lahtr

mv ${genome_fasta} TEMP && mv TEMP ${id}.fasta.gz

(gzip -t ${id}.fasta.gz && echo "${genome_fasta} is in gzip format") || ( echo "${genome_fasta} is NOT in gzip format" && exit 1 )

"""
}


// Combine remote files into tar files with ${batchsize} genomes each
process combineRemoteFiles {
    container 'ubuntu:20.04'
    label 'io_limited'
    errorStrategy 'retry'

    input:
        file fasta_list

    output:
        file "genomes_fasta.*.tar"

"""
#!/bin/bash

set -e

tar cvfh \$(mktemp genomes_fasta.XXXXXXXXX).tar ${fasta_list}

"""
}

// Take a set of genomes and make a single tar file which contains
// a combined_genomes.fasta.gz file with the concatenated set of genomes
// and a combined_genomes.csv.gz file with the headers from each genome,
// with the headers `genome` and `contig`
process combineGenomes {
    container 'ubuntu:20.04'
    label 'io_limited'
    errorStrategy 'retry'

    input:
        file genome_tar

    output:
        file "combined*.tar"

"""
#!/bin/bash

set -e

# Untar the set of genomes in this batch
echo -e "\\nUnpacking input tarball"
tar xvf ${genome_tar}

echo -e "\\nSetting up header CSV"
echo genome,contig > combined_genomes.csv

echo -e "\\nIterating over all unpacked FASTA files"
for fp in *.fasta.gz; do

    echo \$fp

    # Get the genome ID from the file name
    genome_name=\$( echo \$fp | sed 's/.fasta.gz//' )

    # Add the genome ID to the contig headers
    gunzip -c \$fp | \
        sed "s/>/>\${genome_name}_/" | \
        gzip -c > \$fp.TEMP
    mv \$fp.TEMP \$fp

    # Add the contigs to a concatenated file
    cat \$fp >> combined_genomes.fasta.gz

    # Add a line to the CSV with the genome ID and the contig header name
    gunzip -c \$fp | fgrep '>' | sed "s/>/\$genome_name,/" | sed 's/ .*//' >> combined_genomes.csv

done

echo -e "\\nNumber of headers in CSV \$( cat combined_genomes.csv | wc -l )"
echo -e "\\nNumber of headers in FASTA \$( gunzip -c combined_genomes.fasta.gz | grep -c '>' )"

# Make sure that neither file is empty
(( \$( cat combined_genomes.csv | fgrep -v genome,contig | wc -l ) > 0 ))
echo -e "\\nPass - Header file is not empty"
(( \$( gunzip -c combined_genomes.fasta.gz | grep -c '>' ) > 0))
echo -e "\\nPass - FASTA file is not empty"

# Make sure the two numbers match
(( \$( cat combined_genomes.csv | fgrep -v genome,contig | wc -l ) == \$( gunzip -c combined_genomes.fasta.gz | grep -c '>' ) ))
echo -e "\\nPass - Number of headers matches"

echo -e "\\nCompressing header CSV"
gzip combined_genomes.csv

echo -e "\\nRenaming combined CSV and FASTA"
prefix=\$(mktemp combined_genomes.XXXXXXXXX)
mv combined_genomes.fasta.gz \$prefix.fasta.gz
mv combined_genomes.csv.gz \$prefix.csv.gz

# Make a tarball with both files
echo -e "\\nMaking single output tarball"
tar cvfh \$prefix.tar \$prefix.csv.gz \$prefix.fasta.gz

echo -e "\\nDone"
"""
}

// Combine all of the individual tarballs into a single tarball
process makeDatabase {
    container 'ubuntu:20.04'
    label 'io_limited'
    errorStrategy 'retry'
    publishDir params.output_folder

    input:
        file tar_list
        file "database_manifest.csv"

    output:
        file "${params.output_prefix}.tar"

"""
#!/bin/bash

set -e

# Make a tarball with all of the inputs
tar cvfh ${params.output_prefix}.tar ${tar_list} database_manifest.csv

"""
}