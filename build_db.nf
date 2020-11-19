#!/usr/bin/env nextflow

// Using DSL-2
nextflow.preview.dsl=2

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false
params.batchsize = 100

// Import modules
include {fetchFTP as fetchFTP_fasta} from './modules/modules' params(
    file_suffix: ".fasta.gz",
    tar_prefix: "genomes_fasta"
)
include {fetchFTP as fetchFTP_gff} from './modules/modules' params(
    file_suffix: ".gff.gz",
    tar_prefix: "genomes_gff"
)
include {renameRemoteFiles as renameRemoteFiles_fasta} from './modules/modules' params(
    file_suffix: ".fasta.gz"
)
include {renameRemoteFiles as renameRemoteFiles_gff} from './modules/modules' params(
    file_suffix: ".gff.gz"
)
include {combineRemoteFiles as combineRemoteFiles_fasta} from './modules/modules' params(
    tar_prefix: "genomes_fasta"
)
include {combineRemoteFiles as combineRemoteFiles_gff} from './modules/modules' params(
    tar_prefix: "genomes_gff"
)
include {repackHDF} from './modules/modules'

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

    // Validate the manifest
    validateManifest(
        file(params.manifest)
    )

    // Parse the manifest for the URI to the genome FASTA
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
    fetchFTP_fasta(
        split_genome_ch.ftp.map{r -> "${r[1]}:::${r[0]}"}.toSortedList().flatten().collate(params.batchsize)
    )

    // Make sure the other remote files are named correctly
    renameRemoteFiles_fasta(
        split_genome_ch.other.map{r -> [r[0], file(r[1])]}
    )

    // Combine the remote files
    combineRemoteFiles_fasta(
        renameRemoteFiles_fasta.out.toSortedList().flatten().collate(params.batchsize)
    )

    // Now join the channels together
    genome_ch = combineRemoteFiles_fasta.out.mix(
        fetchFTP_fasta.out
    )

    // Reformat each batch of genomes
    combineGenomes(
        genome_ch
    )

    // Parse the manifest for the URI to the genome annotations in GFF format
    Channel.from(
        file(params.manifest)
    ).splitCsv(
        header: true
    ).filter(
        { it.gff != null }
    ).map {
        r -> [r["id"], r["gff"]]
    }.branch {
        ftp: it[1].startsWith("ftp://")
        other: !it[1].startsWith("ftp://")
    }.set {
        split_gff_ch
    }

    // Fetch files from FTP
    fetchFTP_gff(
        split_gff_ch.ftp.map{r -> "${r[1]}:::${r[0]}"}.toSortedList().flatten().collate(params.batchsize)
    )

    // Make sure the other remote files are named correctly
    renameRemoteFiles_gff(
        split_gff_ch.other.map{r -> [r[0], file(r[1])]}
    )

    // Combine the remote files
    combineRemoteFiles_gff(
        renameRemoteFiles_gff.out.toSortedList().flatten().collate(params.batchsize)
    )

    // Now join the channels together
    gff_ch = combineRemoteFiles_gff.out.mix(
        fetchFTP_gff.out
    )

    // Format all of the annotations as HDF5
    formatAnnotations(
        gff_ch
    )

    // Join all of the annotations into a single HDF5
    joinAnnotations(
        formatAnnotations.out.toSortedList()
    )

    // Repack and compress the annotations
    repackHDF(
        joinAnnotations.out
    )

    // Make a single tarball
    makeDatabase(
        combineGenomes.out.toSortedList(),
        validateManifest.out,
        repackHDF.out
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

set -Eeuxo pipefail

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

# Make sure that every '>' is the start of a new line
gunzip -c combined_genomes.fasta.gz | \
    sed 's/>/\\n>/g' | \
    sed '/^\$/d' | \
    gzip -c > temp.fasta.gz
mv temp.fasta.gz combined_genomes.fasta.gz

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
        file genome_annotations_hdf

    output:
        file "${params.output_prefix}.tar"

"""
#!/bin/bash

set -e

# Make a tarball with all of the inputs
tar cvfh ${params.output_prefix}.tar ${tar_list} ${genome_annotations_hdf} database_manifest.csv

"""
}

// Read in all of the GFF files and format the annotations consistently
process formatAnnotations {
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
    file gff_tar_list

    output:
    file "genome_annotations.hdf5"

"""
#!/usr/bin/env python3

import os
import pandas as pd
import tarfile

def parse_gff(fp):
    # Function to parse a GFF file
    return pd.read_csv(
        fp, 
        sep="\\t", 
        comment="#",
        compression="gzip",
        names = [
            "contig",
            "source",
            "type",
            "start",
            "end",
            "score",
            "orientation",
            "_",
            "annotation",
        ]
    ).assign(
        to_keep = lambda df: df["type"].isin(["CDS", "gene", "mRNA"])
    ).query(
        "to_keep"
    ).reindex(
        columns=[
            "contig",
            "type",
            "start",
            "end",
            "orientation",
            "annotation"
        ]
    ).reset_index(
        drop=True
    )

with pd.HDFStore("genome_annotations.hdf5", "w") as store:

    for fp in os.listdir("."):
        if fp.startswith("genomes_gff") is False:
            continue
        if fp.endswith("tar") is False:
            continue

        # implicit else
        with tarfile.open(fp, "r:*") as tar:
            for gff_path in tar.getnames():
                df = parse_gff(
                    tar.extractfile(gff_path)
                )
                id_string = gff_path.replace(".gff.gz", "")

                # Write to the HDF5
                df.to_hdf(
                    store,
                    "/annotations/%s" % id_string
                )

"""
}

// Join a set of multiple annotations
process joinAnnotations {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
    file "genome_annotations.*.hdf5"

    output:
    file "genome_annotations.hdf5"

"""
#!/usr/bin/env python3

import os
import h5py

# Get the list of input files
input_fp_list = [
    fp for fp in os.listdir(".")
    if fp.startswith("genome_annotations.") and fp.endswith(".hdf5")
]
print("Preparing to combine %d input files" % len(input_fp_list))

# Open a connection to the output file
output_store = h5py.File("genome_annotations.hdf5", "w")

# Function to copy objects to the output file
def copy_objects(fp):
    # Try to open the file
    try:
        f = h5py.File(fp, "r")
    except:
        print("Could not open %s" % fp)
        return

    # Make a group for the destination, if necessary
    if "annotations" not in output_store.keys():
        output_store.create_group("/annotations")

    # If successful, copy everything over from /annotations
    for k in f["annotations"]:
        print(k)
        f.copy("annotations/%s" % k, output_store["/annotations/"])

    f.close()
    print("Done processing %s" % fp)

# Copy the data
for fp in input_fp_list:
    copy_objects(fp)

# Close the output file
output_store.close()

"""
}