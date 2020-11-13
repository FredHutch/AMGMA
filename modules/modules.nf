// Fetch genomes via FTP
process fetchFTP {
    container 'quay.io/fhcrc-microbiome/wget:latest'
    label 'io_limited'
    errorStrategy "retry"

    input:
        val uri_id_list
    
    output:
        file "${params.tar_prefix}.*.tar"
    
"""
#!/bin/bash

set -e

for uri_id in ${uri_id_list.join(" ")}; do

    uri=\$(echo \$uri_id | sed 's/:::.*//')
    id=\$(echo \$uri_id | sed 's/.*::://')

    echo "Downloading \$id from \$uri"

    # Try to download, and also save the log
    # If the system call returns non-zero, make sure to escape it
    wget -o \$id.log -O \$id${params.file_suffix} \$uri || echo "Encountered problem downloading"

    # Make sure the file is gzip compressed
    gzip -t \$id${params.file_suffix} && echo "\$id is downloaded in proper gzip format" && continue

    # Now it seems that the file was not downloaded in proper gzip format

    # First check to see if the file even exists on the server by looking at the log file
    if (( \$(cat \$id.log | grep -c "No such file") == 1 )); then

        echo "File does not exist on server, skipping \$uri"

    else

        echo "File does exist on server, but was not downloaded correctly (\$uri)"
        exit 1

    fi

done

echo "Making a tar with all genomes in this batch"
tar cfh \$(mktemp ${params.tar_prefix}.XXXXXXXXX).tar *${params.file_suffix}

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
        file "${id}${params.file_suffix}"

"""
#!/bin/bash

set -e

ls -lahtr

mv ${genome_fasta} TEMP && mv TEMP ${id}${params.file_suffix}

(gzip -t ${id}${params.file_suffix} && echo "${genome_fasta} is in gzip format") || ( echo "${genome_fasta} is NOT in gzip format" && exit 1 )

"""
}

process prokka {

    container "quay.io/biocontainers/prokka:1.14.6--pl526_0"
    label "mem_medium"
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true

    input:
    file fasta

    output:
    file "${fasta.name}/${fasta.name}.gff.gz"

"""
#!/bin/bash

set -euxo pipefail

echo Decompressing input FASTA
gunzip -c "${fasta}" > "${fasta}.fasta"

# Remove the input file
rm "${fasta}"

echo Running Prokka
prokka \
    --outdir "${fasta.name}" \
    --prefix "${fasta.name}" \
    --cpus ${task.cpus} \
    "${fasta}.fasta"

echo Compressing outputs

gzip "${fasta.name}"/*

echo Done
"""

}

process preprocessFASTA {

    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_limited"
    
    input:
    file fasta

    output:
    file "${fasta}"

    """
#!/usr/bin/env python3
# Following criteria from https://github.com/ncbi/pgap/wiki/Input-Files
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import re
# Sanitize and write out
def preprocess_fasta(genome, handle):
    seen_headers = set([])
    
    for header, seq in genome.items():
        
        # Make sure the sequence is >= 199 nucleotides
        if len(seq) < 199:
            continue
        # Trim the header to 37 characters
        if len(header) > 37:
            header = header[:37]
        # Only include letters, digits, hyphens (-), underscores (_), periods (.), colons (:), asterisks (*), and number signs (#)
        header = re.sub('[^0-9a-zA-Z-.*#\$_:]', '_', header)
        # All headers are unique
        assert header not in seen_headers, "Duplicate header: %s (note truncation to first 37 characters)" % header
        seen_headers.add(header)
        # Make sure there are no N's at the beginning or end
        assert seq[0] != "#"
        assert seq[-1] != "#"
        handle.write(">%s\\n%s\\n" % (header, seq))

# Parse the filepath
fasta_fp = "${fasta}"

# Read in all of the genome
if fasta_fp.endswith(".gz"):
    genome = dict([
        (header, seq)
        for header, seq in SimpleFastaParser(gzip.open(fasta_fp, "rt"))
    ])
else:
    genome = dict([
        (header, seq)
        for header, seq in SimpleFastaParser(open(fasta_fp, "r"))
    ])

# Gzip the output if appropriate
if fasta_fp.endswith(".gz"):
    with gzip.open(fasta_fp, "wt") as handle:
        preprocess_fasta(genome, handle)

else:
    with open(fasta_fp, "w") as handle:
        preprocess_fasta(genome, handle)

    """

}

// Combine remote files into tar files with ${batchsize} genomes each
process combineRemoteFiles {
    container 'ubuntu:20.04'
    label 'io_limited'
    errorStrategy 'retry'

    input:
        file file_list

    output:
        file "${params.tar_prefix}.*.tar"

"""
#!/bin/bash

set -e

tar cvfh \$(mktemp ${params.tar_prefix}.XXXXXXXXX).tar ${file_list}

"""
}

// Repack an HDF5 file
process repackHDF {

    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label "mem_medium"
    errorStrategy "retry"
    
    input:
    file output_hdf5
        
    output:
    file "${output_hdf5}"

    """
#!/bin/bash

set -e

[ -s ${output_hdf5} ]

h5repack -f GZIP=5 ${output_hdf5} TEMP && mv TEMP ${output_hdf5}
    """
}