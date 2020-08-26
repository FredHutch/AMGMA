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
    if [[ `gzip -t \$id${params.file_suffix}` ]]; then

        echo "\$id${params.file_suffix} is downloaded in proper gzip format"

    else

        # Now it seems that the file was not downloaded in proper gzip format

        # First check to see if the file even exists on the server by looking at the log file
        if (( cat \$id.log | grep -c "No such file" )); then

            echo "File does not exist on server, skipping \$uri"

        else

            echo "File does exist on server, but was not downloaded correctly (\$uri)"
            exit 1

        fi

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