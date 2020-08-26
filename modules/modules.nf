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

    wget --quiet -O \$id${params.file_suffix} \$uri

    # Make sure the file is gzip compressed
    (gzip -t \$id${params.file_suffix} && echo "\$id${params.file_suffix} is in gzip format") || ( echo "\$id${params.file_suffix} is NOT in gzip format" && exit 1 )

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
