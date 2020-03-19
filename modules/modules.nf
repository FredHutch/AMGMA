#!/usr/bin/env nextflow

// Validate that all genomes are unique in the manifest
process validateManifest {
    tag "Enforce unique genome IDs"
    container "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"
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
        file "genomes_fasta.tar"
    
"""
#!/bin/bash
set -e

for uri_id in ${uri_id_list.join(" ")}; do

    uri=\$(echo \$uri_id | sed 's/:::.*//')
    id=\$(echo \$uri_id | sed 's/.*::://')

    echo "Downloading \$id from \$uri"

    wget -O \$id.fasta.gz \$uri

    # Make sure the file is gzip compressed
    (gzip -t \$id.fasta.gz && echo "\$id.fasta.gz is in gzip format") || ( echo "\$id.fasta.gz is NOT in gzip format" && exit 1 )

done

echo "Making a tar with all genomes in this batch"
tar cvfh genomes_fasta.tar *.fasta.gz

echo "done"

"""
}


// Rename the files from S3 to explicitly match the provided ID
process renameRemoteFiles {
    container 'ubuntu:18.04'
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
    container 'ubuntu:18.04'
    label 'io_limited'
    errorStrategy 'retry'

    input:
        file fasta_list

    output:
        file "genomes_fasta.tar"

"""
#!/bin/bash

set -e

tar cvfh genomes_fasta.tar ${fasta_list}

"""
}


// Annotation of coding sequences with prodigal
process prodigal {
    tag "Identify protein-coding genes"
    container 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
    label 'mem_medium'
    errorStrategy "retry"

    input:
        file genome_fasta_tar
    
    output:
        file "*faa.tar"
        file "*gff.tar"
    
"""
#!/bin/bash
set -e

# Unpack the tar, all of which match *.fasta.gz
tar xvf ${genome_fasta_tar}

# Show the working directory for debugging purposes
ls -lahtr

# Process each of the genomes in this batch
for genome_fasta in *.fasta.gz; do

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

    rm \${genome_id}.fasta

done

# Make the output tar files
tar cvfh \$(mktemp combined.XXXXXXXXX).gff.tar *gff.gz
tar cvfh \$(mktemp combined.XXXXXXXXX).faa.tar *faa.gz

"""
}

// Combine a set of gene sequences, input as tarballs of .faa.gz
process combineFAA {
    tag "Aggregate gene sequences"
    container 'ubuntu:20.04'
    label 'io_limited'
    errorStrategy "retry"

    input:
        file gene_fasta_tar
    
    output:
        file "*faa.gz"
    
"""
#!/bin/bash
set -e

# Unpack the tar, all of which members match *.faa.gz
tar xvf ${gene_fasta_tar}

# Combine all of the files
temp_filename=\$(mktemp combined.XXXXXXXXX)
cat *faa.gz > \$temp_filename
rm *faa.gz
mv  \$temp_filename \$temp_filename.faa.gz

"""
}

// Combine gene annotation tables
process combineGFF {
    tag "Combine gene annotation tables"
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file gff_tar
    
    output:
        file "combined.*.csv.gz"
    
"""
#!/usr/bin/env python3

import os
import pandas as pd
import tarfile
import uuid

# Extract all of the files from the input tarfile
tf = tarfile.open("${gff_tar}", mode = "r")
tf.extractall()
tf.close()

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
output_fp = "combined.%s.csv.gz" % str(uuid.uuid4())
pd.concat([
    read_gff(fp)
    for fp in os.listdir(".")
    if fp.endswith("gff.gz")
]).to_csv(
    output_fp,
    index=None,
    compression="gzip"
)
"""
}


// Cluster genome genes by identity to find centroids
process clusterGenes {
    tag "Cluster genes with similar sequences"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file fasta_list
    
    output:
    file "centroids.*.fasta.gz"
    
    """
#!/bin/bash

set -e

echo ""
df -h
ls -lahtr
echo ""

# Make sure that all input files were staged
for fp in ${fasta_list}; do
    echo "Checking for \$fp"

    [[ -s \$fp ]]

done

# Make a single file with all of the input files
echo "Combining input files"
cat ${fasta_list}  > combined.fasta.gz
echo "Combining input files - done"
rm ${fasta_list}

echo ""
df -h
ls -lahtr
echo ""

# Make the MMSeqs2 database
echo "Making MMSeqs2 database"
mmseqs createdb combined.fasta.gz db
echo "Making MMSeqs2 database - done"

echo ""
df -h
ls -lahtr
echo ""

# Cluster the protein sequences
echo "Running clustering"
mmseqs cluster db cluster_db ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 300 \
    --single-step-clustering \
    -c ${params.min_coverage / 100}
echo "Running clustering - done"

echo ""
df -h
ls -lahtr
echo ""

# Get the representative sequences
echo "Making output FASTA"
centroids_fasta=\$(mktemp centroids.XXXXX).fasta
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes \$centroids_fasta --use-fasta-header
# Compress
gzip \$centroids_fasta
echo "Making output FASTA - done"

echo "Cleaning up working files"

rm latest db* genes* cluster_db* combined.fasta.gz

echo ""
df -h
ls -lahtr
echo ""

echo "Done"
    """
}

// Assign centroids to every individual genome
process assignCentroids {
    tag "Annotate genomes with centroids"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
    file genome_tar
    file catalog_dmnd

    output:
    file "*tar"

"""
#!/bin/bash

set -e

tar xvf ${genome_tar}

for genome_fasta in *faa.gz; do

    echo "Processing \$genome_fasta"
    
    [[ -s \$genome_fasta ]]

    # Find the single top hit per query
    diamond \
        blastp \
        --query \$genome_fasta \
        --db "${catalog_dmnd}" \
        --threads ${task.cpus} \
        --out \$genome_fasta.aln.gz \
        --threads ${task.cpus} \
        --outfmt 6 qseqid sseqid \
        --query-cover ${100 - ((100 - params.min_coverage) * 2)} \
        --id ${100 - ((100 - params.min_identity) * 2)} \
        --max-target-seqs 1 \
        --block-size ${task.memory.toGiga() / 6} \
        --compress 1

done

echo "Making single output tarball"

tar cvf \$(mktemp centroids.XXXXXX).tar *.aln.gz

echo Done


"""

}



// Make a final genome summary table
process makeHDF {
    tag "Make a single output HDF"
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'mem_veryhigh'
    errorStrategy "retry"

    input:
        file combined_csv_list
        file manifest_csv
        file centroids_tar_list
    
    output:
        file "${params.output_prefix}.hdf"
    
"""
#!/usr/bin/env python3

import os
import pandas as pd
import  tarfile

# Read in the tarball which identifies which centroid gene each genome gene corresponds to
genome_gene_catalog_map = dict()

for centroids_tar in "${centroids_tar_list}".split(" "):
    print("Processing %s" % centroids_tar)
    assert os.path.exists(centroids_tar)

    tar = tarfile.open(centroids_tar)

    for member in tar:

        print("Element: %s" % member.name)
        
        df = pd.read_csv(
            tar.extractfile(member),
            header = None,
            sep = "\\t",
            compression = "gzip"
        )

        # Assign the mappings to a dict
        genome_gene_catalog_map[
            member.name.replace(".faa.gz.aln.gz", "")
        ] = df.set_index(
            0
        )[
            1
        ].to_dict()


    tar.close()


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
        for fp in "${combined_csv_list}".split(" ")
    ]).groupby("genome_id"):

        assert genome_id in genome_gene_catalog_map, genome_id

        genome_df.drop(
            columns = "genome_id"
        ).assign(
            centroid = genome_df[
                "gene_id"
            ].apply(
                genome_gene_catalog_map[
                    genome_id
                ].get
            )
        ).to_hdf(
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

"""
}


// Repack and compress the final HDF5 file
process repackHDF {

    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    tag "Compress HDF store"
    label "mem_veryhigh"
    errorStrategy "retry"
    publishDir "${params.output_folder}"
    
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