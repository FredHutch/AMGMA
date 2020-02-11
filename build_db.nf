#!/usr/bin/env nextflow

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/AMGMA/build_db.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)
      --output_folder       Folder to write output files to
      --output_prefix       Prefix to use for output file names

    Manifest:
      The manifest is a CSV listing all of the genomes to be used for the database.
      The manifest much contain the column headers: uri,id,name
      The URI may start with ftp://, s3://, or even just be a path to a local file. 
      The ID must be unique, and only contain a-z, A-Z, 0-9, or _.
      The NAME may be longer and contain whitespaces, but may not contain a comma.

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Parse the manifest
genome_ch = Channel.from(
    file(params.manifest)
).splitCsv(
    header: true
).map {
    r -> [r["id"], file(r["uri"])]
}

// Annotation of coding sequences with prodigal
process prodigal {
    tag "Identify protein-coding genes"
    container 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
    label 'io_limited'
    errorStrategy "retry"

    input:
        tuple val(id), file(genome_fasta) from genome_ch
    
    output:
        file "${id}.faa.gz" into fasta_ch
        file "${id}.gff.gz" into gff_ch
    
"""
#!/bin/bash
set -e 

# Decompress the genome if it is GZIP compressed
if [[ \$(gzip -d ${genome_fasta}) ]]; then
    gunzip -c ${genome_fasta} > ${genome_fasta}.fasta
else
    mv ${genome_fasta} ${genome_fasta}.fasta
fi

# Add the genome ID to the FASTA headers
sed -i 's/>/>${id}_/g' ${genome_fasta}.fasta

prodigal \
    -a ${id}.faa \
    -i  ${genome_fasta}.fasta \
    -f gff \
    -o ${id}.gff \
    -p meta

# Add the genome ID to the GFF output
sed -i "s/^${id}/${id}\t${id}/g" ${id}.gff

gzip ${id}.gff
gzip ${id}.faa

"""
}

// Combine gene annotation tables
process combineGFF {
    tag "Combine gene annotation tables"
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file gff_gz from gff_ch.toSortedList().flatten().collate(100)
    
    output:
        file "*.csv.gz" into csv_ch
    
"""
#!/usr/bin/env python3

import pandas as pd

# Function to read in a single GFF file
def read_gff(fp):
    df = pd.read_csv(
        fp, 
        compression = "gzip", 
        sep = "\\t",
        comment = "#",
        header = None
    ).rename(columns=dict([
        (0, 'genome_id'),
        (1, 'gene_id'),
        (3, 'type'),
        (4, 'start'),
        (5, 'end'),
        (7, 'strand'),
        (9, 'tags')
    ])).query(
        "type == 'CDS'"
    )

    # Append "_N" to the gene IDs
    df["gene_id"] = [
        "%s_%d" % (s, ix + 1)
        for ix, s in enumerate(df["gene_id"].values)
    ]

    df = pd.concat(
        [
            df.reindex(columns=[
                'genome_id',
                'gene_id',
                'start',
                'end',
                'strand',
            ]),
            pd.DataFrame([
                dict([
                    v.split("=", 1)
                    for v in s.split(";")
                    if "=" in v
                ])
                for s in df["tags"].values
            ])
        ],
        axis = 1,
        sort = True
    )
    return df

# Read in all of the GFF files and write to a file
output_fp = "${gff_gz}".split(" ", 1)[0] + ".csv.gz"
pd.concat([
    read_gff(fp)
    for fp in "${gff_gz}".split(" ")
]).to_csv(
    output_fp,
    index=None,
    compression="gzip"
)
"""
}


// Make a final genome summary table
process makeHDF {
    tag "Make a single output HDF"
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'io_limited'
    errorStrategy "retry"
    publishDir "${params.output_folder}"

    input:
        file csv_list from csv_ch.collect()
        file manifest_csv from file(params.manifest)
    
    output:
        file "${params.output_prefix}.hdf"
    
"""
#!/usr/bin/env python3

import pandas as pd

# Open a connection to the output HDF
with pd.HDFStore("${params.output_prefix}.hdf", "w") as store:

    # Write all of the gene annotations to /genomes/<genome_id>
    for genome_id, genome_df in pd.concat([
        pd.read_csv(
            fp,
            sep=",",
            compression="gzip"
        )
        for fp in "${csv_list}".split(" ")
    ]).groupby("genome_id"):
    
        genome_df.to_hdf(
            store,
            "/genomes/%s" % genome_id,
            format = "fixed"
        )

    # Write the manifest to /manifest
    pd.read_csv(
        "${manifest_csv}",
        sep=","
    ).to_hdf(
        store,
        "/manifest",
        format = "fixed"
    )

"""
}

// Combine gene FASTA files
process combineFASTA {
    tag "Combine gene FASTA files"
    container "ubuntu:18.04"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file fasta_gz from fasta_ch.toSortedList().flatten().collate(100)
    
    output:
        file "combined.fasta.gz" into combined_fasta_ch
    
"""
#!/bin/bash

set -e

cat ${fasta_gz} > combined.fasta.gz
"""
}


// Make a DIAMOND alignment database of all genes
process diamondDB {
    tag "Make a DIAMOND database"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    publishDir "${params.output_folder}"
    
    input:
    file "input.*.fasta.gz" from combined_fasta_ch.collect()

    output:
    file "${params.output_prefix}.dmnd"

    """
#!/bin/bash

set -e

# Combine all of the input files
cat input.*.fasta.gz > input.fasta.gz

diamond \
    makedb \
    --in input.fasta.gz \
    --db ${params.output_prefix}.dmnd \
    --threads ${task.cpus}
    """

}