#!/usr/bin/env nextflow

// Set default parameters
params.help = false
params.manifest = false
params.output_folder = false
params.output_prefix = false
params.min_identity = 90
params.min_coverage = 50
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
      --min_identity        Percent identity threshold used to pick centroids (default: 90)
      --min_coverage        Percent coverage threshold used to pick centroids (default: 50)
      --batchsize           Number of genomes to download in a batch (default: 100)

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

// Fetch genomes via FTP
process fetchFTP {
    tag "Download genomes hosted by FTP"
    container 'quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f'
    label 'io_limited'
    errorStrategy "retry"

    input:
        val uri_id_list from split_genome_ch.ftp.map{r -> "${r[1]}:::${r[0]}"}.toSortedList().flatten().collate(params.batchsize)
    
    output:
        file "*.fasta.gz" into downloaded_genome_ch
    
"""
#!/bin/bash
set -e

for uri_id in ${uri_id_list.join(" ")}; do

    uri=\$(echo \$uri_id | sed 's/:::.*//')
    id=\$(echo \$uri_id | sed 's/.*::://')

    echo "Downloading \$id from \$uri"

    # Try multiple times
    for _ in {1..5}; do
        ( wget -O \$id.fasta.gz \$uri && break ) || sleep 0.5
    done

    # Make sure the file is gzip compressed
    (gzip -t \$id.fasta.gz && echo "\$id.fasta.gz is in gzip format") || ( echo "\$id.fasta.gz is NOT in gzip format" && exit 1 )

done

echo "done"

"""
}

// Rename the files from S3 to explicitly match the provided ID
process renameRemoteFiles {
    container 'ubuntu:18.04'
    label 'io_limited'
    errorStrategy 'retry'

    input:
        tuple val(id), file(genome_fasta) from split_genome_ch.other.map{r -> [r[0], file(r[1])]}

    output:
        file "${id}.fasta.gz" into renamed_genome_ch

"""
#!/bin/bash

set -e

ls -lahtr

mv ${genome_fasta} ${id}.fasta.gz

(gzip -t ${id}.fasta.gz && echo "${genome_fasta} is in gzip format") || ( echo "${genome_fasta} is NOT in gzip format" && exit 1 )

"""
}

// Now join the channels together
genome_ch = renamed_genome_ch.mix(
    downloaded_genome_ch.flatten()
).toSortedList(
).flatten(
)

// Annotation of coding sequences with prodigal
process prodigal {
    tag "Identify protein-coding genes"
    container 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
    label 'io_limited'
    errorStrategy "retry"

    input:
        file genome_fasta_list from genome_ch.collate(params.batchsize)
    
    output:
        file "*.faa.gz" into fasta_ch
        file "*.gff.gz" into gff_ch
    
"""
#!/bin/bash
set -e

# Process each of the genomes in this batch
for genome_fasta in ${genome_fasta_list}; do

    # Parse the ID from the file name
    genome_id=\${genome_fasta%.fasta.gz}

    echo "Processing \$genome_id from \$genome_fasta"

    # Decompress the genome
    gunzip \$genome_fasta

    # Add the genome ID to the FASTA headers
    sed -i "s/>/>\${genome_id}_/g" \${genome_id}.fasta

    prodigal \
        -a \${genome_id}.faa \
        -i  \${genome_id}.fasta \
        -f gff \
        -o \${genome_id}.gff \
        -p meta

    # Add the genome ID to the GFF output
    sed -i "s/^\${genome_id}/\${genome_id}\t\${genome_id}/g" \${genome_id}.gff

    gzip \${genome_id}.gff
    gzip \${genome_id}.faa

done

"""
}

// Combine gene annotation tables
process combineGFF {
    tag "Combine gene annotation tables"
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file gff_gz from gff_ch.flatten().toSortedList().flatten().collate(params.batchsize)
    
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


// Combine gene FASTA files
process combineFASTA {
    tag "Combine gene FASTA files"
    container "ubuntu:18.04"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file fasta_gz from fasta_ch.flatten().toSortedList().flatten().collate(params.batchsize)
    
    output:
        file "combined.fasta.gz" into combined_fasta_ch
    
"""
#!/bin/bash

set -e

cat ${fasta_gz} > combined.fasta.gz
"""
}


// Cluster genome genes by identity to find centroids
process clusterGenes {
    tag "Cluster genes with similar sequences"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file "input.*.fasta.gz" from combined_fasta_ch.collect()
    
    output:
    file "centroids.fasta.gz" into centroids_fasta
    file "centroids.tsv.gz" into centroids_tsv
    
    """
#!/bin/bash

set -e

# Make a single file with all of the input files
echo "Combining input files"
gunzip -c input.*.fasta.gz | gzip -c > combined.fasta.gz
echo "Combining input files - done"
rm input.*.fasta.gz

# Make the MMSeqs2 database
echo "Making MMSeqs2 database"
mmseqs createdb combined.fasta.gz db
echo "Making MMSeqs2 database - done"

# Cluster the protein sequences
echo "Running clustering"
mmseqs linclust db cluster_db ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 100000 \
    -c ${params.min_coverage / 100}
echo "Running clustering - done"

# Make TSV output for clustering
echo "Making output TSV"
mmseqs createtsv db db cluster_db centroids.tsv
# Compress
gzip centroids.tsv
echo "Making output TSV - done"

# Get the representative sequences
echo "Making output FASTA"
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes centroids.fasta --use-fasta-header
# Compress
gzip centroids.fasta
echo "Making output FASTA - done"
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
        file centroids_tsv
    
    output:
        file "${params.output_prefix}.hdf"
    
"""
#!/usr/bin/env python3

import pandas as pd

# Read the table linking genome genes to centroid genes
centroid_df = pd.read_csv(
    "${centroids_tsv}", 
    sep="\\t",
    header = None
).reindex(
    columns = [0, 1]
).rename(
    columns = dict([
        (0, 'genome_centroid'),
        (1, 'genome_gene')
    ])
)

# Open a connection to the output HDF
with pd.HDFStore("${params.output_prefix}.hdf", "w") as store:

    # Write all of the gene annotations to /genomes/<genome_id>
    for genome_id, genome_df in pd.concat([
        pd.read_csv(
            fp,
            sep=",",
            compression="gzip"
        ).reindex(
            columns = [
                "genome_id",
                "gene_id",
                "start",
                "end",
                "strand",
                "conf",
                "gc_cont"
            ]
        )
        for fp in "${csv_list}".split(" ")
    ]).groupby("genome_id"):

        genome_df.to_hdf(
            store,
            "/genomes/%s" % genome_id,
            format = "fixed",
            complevel = 5
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

    # Write the table of gene centroids
    centroid_df.to_hdf(
        store,
        "/genome_centroids",
        format = "fixed",
        complevel = 5
    )

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
      file centroids_fasta

    output:
      file "${params.output_prefix}.dmnd"

    """
#!/bin/bash

set -e

diamond \
    makedb \
    --in ${centroids_fasta} \
    --db ${params.output_prefix}.dmnd \
    --threads ${task.cpus}
    """

}