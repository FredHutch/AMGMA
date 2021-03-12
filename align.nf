#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.db = null
params.geneshot_folder = null
params.geneshot_results_hdf = null
params.geneshot_details_hdf = null
params.geneshot_dmnd = null
params.output_folder = null
params.output_prefix = null
params.min_coverage = 50
params.min_identity = 80
params.top = 5
params.formula = false
params.fdr_method = "fdr_bh"
params.alpha = 0.2
params.blast = false
params.no_associations = false
// Only show alignments for CAGs with at least this number of genes aligned
params.min_genes_per_cag = 2
// Align against any contigs with this minimum size (or circular)
params.min_contig_len = 25000

// Commonly used containers
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_fastcluster"
container_diamond = "quay.io/fhcrc-microbiome/docker-diamond:v2.0.6-biopython"
container__blast = "ncbi/blast:2.11.0"

// Import externally-defined processes
include { 
    runCorncob;
    joinCorncob
} from './modules/statistics.nf'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/AMGMA <ARGUMENTS>

Required Arguments:
--db                  AMGMA database(s) (ends with .tar), with commas delimiting multiple databases
--geneshot_folder     Folder containing geneshot results, which includes:
                        - Results HDF file output by geneshot (*results.hdf5)
                          [overridden by --geneshot_results_hdf]
                        - Detailed HDF file output by geneshot (*details.hdf5)
                          [overridden by --geneshot_details_hdf]
                        - DIAMOND database containing the gene catalog (ref/genes.dmnd)
                          [overridden by --geneshot_dmnd]
--output_folder       Folder used to write outputs
--output_prefix       Prefix appended to files in the output directory

Optional Arguments:
--min_coverage        Minimum coverage required for alignment (default: 50)
--min_identity        Minimum percent identity required for alignment (default: 80)
--top                 Threshold used to retain overlapping alignments within --top% score of the max score (default: 5)
--formula             If specified, calculate association of genome abundances with experimental design
--fdr_method          Method used for FDR correction (default: fdr_bh)
--alpha               Alpha value used for FDR correction (default: 0.2)
--blast               Align with BLAST+ instead of DIAMOND
--no_associations     Exclude all analysis of CAG association metrics
--min_genes_per_cag   Only show alignments for CAGs with at least this number of genes aligned (default: 2)

Output HDF:
The primary output from AMGMA is formatted in HDF5 (OUTPUT_PREFIX.hdf5) to combine multiple tables
into a single file. An additional output file is formatted as a Redis database file (*.rdb) for rapid
retrieval. Data includes:

* /genomes/manifest (Table with the description of all genomes used for alignment)
* /genomes/cags/containment (Table with the number of genes aligned for all CAGs against all genomes)
* /genomes/raw_abundance (Table with the proportion of gene copies per genome, per specimen)
* /genomes/imputed_abundance (Table with the imputed relative abundance per genome, per specimen)
* /genomes/detail/<genome_id> (Table with the coordinates of all aligned genes for a given genome)
* /genomes/annotation/<genome_id> (Table with the annotations from optional GFF input)
* /genomes/association/<parameter> (Table with estimated associations of all genomes against an experimental parameter)

    """.stripIndent()
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.geneshot_folder == null || params.db == null || params.output_folder == null || params.output_prefix == null){
        // Invoke the function above which prints the help message
        helpMessage()

        if (params.geneshot_folder == null){
            log.info"""
            Please provide --geneshot_folder
            """.stripIndent()
        }
        if (params.db == null){
            log.info"""
            Please provide --db
            """.stripIndent()
        }
        if (params.output_folder == null){
            log.info"""
            Please provide --output_folder
            """.stripIndent()
        }
        if (params.output_prefix == null){
            log.info"""
            Please provide --output_prefix
            """.stripIndent()
        }

        // Exit out and do not run anything else
        exit 1
    }

    // Point to the files with the GeneShot results

    // Results HDF5
    if ( params.geneshot_results_hdf == null ) {
        log.info("Parsing Results HDF: ${params.geneshot_folder}*results.hdf5")
        geneshot_results_hdf = file("${params.geneshot_folder}*results.hdf5")
    } else {
        log.info("Parsing Results HDF: ${params.geneshot_results_hdf}")
        geneshot_results_hdf = file("${params.geneshot_results_hdf}")
    }
    log.info("Found Results HDF: ${geneshot_results_hdf}")

    // Details HDF5
    if ( params.geneshot_details_hdf == null ) {
        log.info("Parsing Details HDF: ${params.geneshot_folder}*details.hdf5")
        geneshot_details_hdf = file("${params.geneshot_folder}*details.hdf5")
    } else {
        log.info("Parsing Details HDF: ${params.geneshot_details_hdf}")
        geneshot_details_hdf = file("${params.geneshot_details_hdf}")
    }
    log.info("Found Details HDF: ${geneshot_details_hdf}")

    // Gene catalog
    if ( params.geneshot_dmnd == null ) {
        log.info("Parsing Gene Catalog DMND: ${params.geneshot_folder}**genes.dmnd")
        geneshot_dmnd = file("${params.geneshot_folder}**genes.dmnd")
    } else {
        log.info("Parsing Gene Catalog DMND: ${params.geneshot_dmnd}")
        geneshot_dmnd = file("${params.geneshot_dmnd}")
    }
    log.info("Found Gene Catalog DMND: ${geneshot_dmnd}")

    // Redis store
    log.info("Parsing Redis Store: ${params.geneshot_folder}**.rdb")
    geneshot_rdb = file("${params.geneshot_folder}**rdb")
    log.info("Found Redis Store: ${geneshot_rdb}")

    // Make sure the input database file(s) exist
    db_ch = Channel.from(
        params.db.split(",")
    ).map { 
        fp -> file(fp) 
    }

    // Unpack the database(s)
    unpackDatabase(db_ch)

    // Make a channel with any assemblies available
    contig_fasta_ch = Channel.fromPath(
        "${params.geneshot_folder}**.contigs.fasta.gz"
    )

    // Extract the long contigs from those assemblies
    filterContigs(
        contig_fasta_ch
    )

    // Make a channel with all references in tar format
    genomes_tar_ch = unpackDatabase.out[1].flatten().mix(
        filterContigs.out[1]
    )

    // If multiple databases were provided, join the manifests
    // Make sure that there are no overlapping IDs
    validateManifest(
        unpackDatabase.out[0].mix(
            filterContigs.out[0]
        ).toSortedList()
    )

    // If using BLAST+, format the database
    if (params.blast) {

        // Convert the gene catalog from DMND to FASTA
        extractFASTA(
            geneshot_dmnd
        )

        // Format the BLAST database
        makeBLASTdb(
            extractFASTA.out
        )

        // Align the genomes against the database with BLAST
        alignGenomesBLAST(
            genomes_tar_ch,
            makeBLASTdb.out
        )

        // Filter alignments by percent identity and coverage of the gene catalog entry
        filterBLASThits(
            alignGenomesBLAST.out
        )
        alignments_ch = filterBLASThits.out

    // Otherwise
    } else {

        // Align the genomes against the database with DIAMOND
        alignGenomes(
            genomes_tar_ch,
            geneshot_dmnd
        )

        alignments_ch = alignGenomes.out
    }

    // Drop files with zero alignments
    filterAlignments(
        alignments_ch
    )

    // Calculate the containment of each CAG in each genome
    calculateContainment(
        filterAlignments.out,
        geneshot_results_hdf
    )

    // If a formula was provided
    if ( params.formula ){

        // Calculate the proportion of gene copies from each specimen which align to each genome
        extractCounts(
            calculateContainment.out[1],
            geneshot_details_hdf
        )

        // Run corncob on the genome abundances
        runCorncob(
            extractCounts.out[0],
            extractFormula.out[0],
            "genome",
            extractFormula.out[1].map {
                it -> it.readLines()
            }.flatten()
        )

        // Join together all of the corncob results
        joinCorncob(
            runCorncob.out.toSortedList(),
            "genome"
        )

        // Create a redis .rdb file
        combineResultsRedisFormula(
            calculateContainment.out[0].toSortedList(),
            calculateContainment.out[1].toSortedList(),
            joinCorncob.out.toSortedList(),
            validateManifest.out,
            geneshot_details_hdf,
            geneshot_results_hdf,
            geneshot_rdb
        )

        // Create an HDF5 file
        combineResultsHDFFormula(
            calculateContainment.out[0].toSortedList(),
            calculateContainment.out[1].toSortedList(),
            joinCorncob.out.toSortedList(),
            validateManifest.out,
            geneshot_details_hdf,
            geneshot_results_hdf,
            geneshot_rdb
        )

    } else {

        // Create a redis .rdb file
        combineResultsRedis(
            calculateContainment.out[0].toSortedList(),
            calculateContainment.out[1].toSortedList(),
            validateManifest.out,
            geneshot_details_hdf,
            geneshot_results_hdf,
            geneshot_rdb
        )

        // Create an HDF5 file
        combineResultsHDF(
            calculateContainment.out[0].toSortedList(),
            calculateContainment.out[1].toSortedList(),
            validateManifest.out,
            geneshot_details_hdf,
            geneshot_results_hdf,
            geneshot_rdb
        )

    }

    // Collect results and combine across all shards


    // Repack an HDF5 file
    repackHDF(
        combineResultsHDF.out
    )

}

// Unpack the database
process unpackDatabase {
    tag "Extract all files from database tarball"
    container "ubuntu:20.04"
    label "io_limited"

    input:
        file db
    
    output:
        file "${db}.manifest.csv"
        file "*tar"
        file "genome_annotations.hdf5"

"""
#!/bin/bash 

set -e

ls -lahtr

tar xvf ${db}

mv database_manifest.csv ${db}.manifest.csv

# Make a dummy file for the genome annotations, if it doesn't exist
if [[ ! -s genome_annotations.hdf5 ]]; then
    touch genome_annotations.hdf5
fi

echo "Done"

"""
}

// If multiple databases were provided, join the manifests
// Make sure that there are no overlapping IDs
process validateManifest {
    
    container "${container__pandas}"
    label 'io_limited'

    input:
        file manifest_csv_list
    
    output:
        file "database_manifest.csv"
    
"""
#!/usr/bin/env python3

import pandas as pd
import os

# Parse the list of files
manifest_csv_list = "${manifest_csv_list}".split(" ")

# Make sure all files are present
for fp in manifest_csv_list:
    assert os.path.exists(fp), "Could not find %s" % fp

# Read in and join the files
manifest_df = pd.concat([
    pd.read_csv(fp)
    for fp in manifest_csv_list
], sort=True)

# Make sure that none of the IDs are shared
vc = manifest_df["id"].value_counts()
assert vc.max() == 1, ("Found duplicated IDs across these databases", vc.head())

# Write out to a new file
manifest_df.to_csv("database_manifest.csv", index=None)

"""
}


// Extract the formula used for the geneshot analysis
process extractFormula {
    
    container "${container__pandas}"
    label 'io_limited'

    input:
        file geneshot_results_hdf
    
    output:
        file "manifest.csv"
        file "formula.txt" optional true
    
"""#!/usr/bin/env python3

import pandas as pd

# Write out the manifest
pd.read_hdf(
    "${geneshot_results_hdf}", 
    "manifest"
).to_csv(
    "manifest.csv"
)

# Read the table of experiment parameters
dat = pd.read_hdf(
    "${geneshot_results_hdf}", 
    "/summary/experiment"
).set_index(
    "variable"
)[
    "value"
].to_dict()

# If the user provided a formula
if "formula" in dat:

    # Save it to formula.txt
    with open('formula.txt', 'wt') as handle:
        for v in dat['formula'].split(","):
            handle.write(v)
"""
}


// Convert the gene catalog from DMND to FASTA
process extractFASTA {
    container "${container_diamond}"
    label "io_limited"

    input:
        file geneshot_dmnd
    
    output:
        file "ref.fasta.gz"


    """#!/bin/bash

    set -Eeuxo pipefail

    diamond getseq --db ${geneshot_dmnd} --out ref.fasta

    gzip ref.fasta
    """
}

// Format the BLAST database
process makeBLASTdb {
    container "${container__blast}"
    label "mem_medium"

    input:
        file ref_fasta
    
    output:
        file "blastDB*"


    """#!/bin/bash

    set -Eeuxo pipefail

    echo "Decompressing gene catalog FASTA"
    gunzip -c ${ref_fasta} > ref.fasta

    echo "Head of gene catalog FASTA"
    head ref.fasta

    echo "Building database"
    makeblastdb -in ref.fasta -dbtype prot -out blastDB

    echo "Done"
    """
}

// Align the genomes against the database with BLAST
process alignGenomesBLAST {
    tag "Annotate reference genomes by alignment"
    container "${container__blast}"
    label "mem_medium"

    input:
        file database_chunk_tar
        file "*"

    output:
        tuple file("${database_chunk_tar.name.replaceAll(/.tar/, ".aln.gz")}"), file("${database_chunk_tar.name.replaceAll(/.tar/, ".csv.gz")}")

"""
#!/bin/bash

set -Eeuxo pipefail

ls -lahtr

tar xvf ${database_chunk_tar}

blastx \
    -query <(gunzip -c ${database_chunk_tar.name.replaceAll(/.tar/, ".fasta.gz")}) \
    -db blastDB \
    -query_gencode 11 \
    -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen" \
    -num_threads ${task.cpus} \
    -max_target_seqs 10000000 \
    -evalue 0.001 \
    | gzip -c > ${database_chunk_tar.name.replaceAll(/.tar/, ".aln.gz")}

"""
}

// Filter alignments by percent identity and coverage of the gene catalog entry
process filterBLASThits {
    container "${container__pandas}"
    label "mem_medium"

    input:
        tuple file(aln_tsv_gz), file(header_csv_gz)

    output:
        tuple file("${aln_tsv_gz}.filtered.tsv.gz"), file("${header_csv_gz}")

"""#!/usr/bin/env python3

import pandas as pd

# Read in the table of hits
df = pd.read_csv("${aln_tsv_gz}", sep="\\t", header=None)
print("Read in %d alignments" % df.shape[0])

# Filter by percent identity
df = df.loc[
    df[2] >= ${params.min_identity}
]
print("Filtered down to %d alignments with identity >= ${params.min_identity}" % df.shape[0])

# Filter by coverage
df = df.assign(
    coverage = 100 * df[3] / df[9]
).query(
    "coverage >= ${params.min_coverage}"
).drop(
    columns = "coverage"
)
print("Filtered down to %d alignments with coverage >= ${params.min_coverage}" % df.shape[0])

# Write out
df.to_csv("${aln_tsv_gz}.filtered.tsv.gz", sep="\\t", index=None, header=None)

"""
}


// Filter assemblies to just those contigs of a certain size
process filterContigs {
    container "${container__pandas}"
    label "io_limited"

    input:
        file contig_fasta_gz

    output:
        file "${contig_fasta_gz}.manifest.csv" optional true
        file "${contig_fasta_gz}.tar" optional true

"""#!/usr/bin/env python3

import gzip
import pandas as pd
import tarfile

def fasta_gz_parser(fp):
    with gzip.open(fp, 'rt') as handle:
        header = None
        seq = None
        for l in handle:
            if l.startswith(">"):
                if header is not None and seq is not None:
                    yield header, seq
                header = l.lstrip(">").rstrip("\\n")
                seq = ""
            else:
                seq = seq + l.rstrip("\\n")

    if header is not None and seq is not None:
        yield header, seq

def is_circular(seq, min_k=15, max_k=31):
    for k in range(min_k, max_k + 1):
        if seq[:k] == seq[-k:]:
            return True
    return False

# Make a genome manifest with columns 'name' and 'id'
genome_manifest = []

# Make a contig manifest with columns 'genome' and 'contig'
contig_manifest = []

# The output will have a *fastq.gz and *contig.manifest.csv.gz in a tarball,
# as well as the genome manifest in a separate file

# Filter the contigs
filtered_contigs = [
    (header, seq)
    for header, seq in fasta_gz_parser("${contig_fasta_gz}")
    if len(seq) >= ${params.min_contig_len} or (len(seq) >= 1000 and is_circular(seq))
]

# Parse the specimen name
specimen_name = "${contig_fasta_gz}".replace(".contigs.fasta.gz", "")

# Format the names of the files to write out
filtered_contigs_path = "${contig_fasta_gz}.fasta.gz"
filtered_contig_manifest_path = "${contig_fasta_gz}.csv.gz"

# If there are any contigs which pass the filter
if len(filtered_contigs) > 0:

    # Write out the contigs
    with gzip.open(filtered_contigs_path, "wt") as handle_out:

        # Iterate over every contig which passes the filter
        for header, seq in filtered_contigs:

            # Write out the contig
            handle_out.write(
                ">%s\\n%s\\n" % (header, seq)
            )

            # Save the contig to both manifest tables
            genome_manifest.append(dict(
                name=header,
                id=header
            ))
            contig_manifest.append(dict(
                genome=header,
                contig=header
            ))

    # Write out both manifests
    pd.DataFrame(
        genome_manifest
    ).to_csv(
        "${contig_fasta_gz}.manifest.csv",
        index=None
    )

    pd.DataFrame(
        contig_manifest
    ).to_csv(
        filtered_contig_manifest_path,
        index=None
    )

    # Tar up the filtered contigs
    with tarfile.open("${contig_fasta_gz}.tar", "w") as tar:
        tar.add(
            filtered_contig_manifest_path
        )
        tar.add(
            filtered_contigs_path
        )

"""
}


// Align the genomes against the database with DIAMOND
process alignGenomes {
    tag "Annotate reference genomes by alignment"
    container "${container_diamond}"
    label "mem_veryhigh"

    input:
        file database_chunk_tar
        file geneshot_dmnd

    output:
        tuple file("${database_chunk_tar.name.replaceAll(/.tar/, ".aln.gz")}"), file("${database_chunk_tar.name.replaceAll(/.tar/, ".csv.gz")}")

"""
#!/bin/bash

set -Eeuxo pipefail

ls -lahtr

tar xvf ${database_chunk_tar}

diamond \
    blastx \
    --db ${geneshot_dmnd} \
    --query ${database_chunk_tar.name.replaceAll(/.tar/, ".fasta.gz")} \
    --out ${database_chunk_tar.name.replaceAll(/.tar/, ".aln.gz")} \
    --outfmt 6 qseqid sseqid pident length qstart qend qlen sstart send slen \
    --id ${params.min_identity} \
    --subject-cover ${params.min_coverage} \
    --top ${params.top} \
    --compress 1 \
    --unal 0 \
    --sensitive \
    --query-gencode 11 \
    --range-culling \
    -F 1 \
    --block-size ${task.memory.toMega() / (1024 * 6 * task.attempt)} \
    --threads ${task.cpus} \
    -c1

ls -lahtr

"""
}


// Drop files with zero alignments
process filterAlignments {
    
    container "${container__pandas}"
    label 'io_limited'

    input:
        tuple file(aln_tsv_gz), file(header_csv_gz)

    output:
        tuple file("${aln_tsv_gz}"), file("${header_csv_gz}") optional true

"""
#!/bin/bash

set -e

if (( \$( cat ${aln_tsv_gz} | wc -l ) > 0 )); then

    echo "Number of alignments: \$( cat ${aln_tsv_gz} | wc -l )"

else

    echo "No alignments found, filtering out"

    rm ${aln_tsv_gz} ${header_csv_gz}

fi

"""

}


// Calculate the containment of each CAG in each genome
process calculateContainment {
    tag "Overlap between CAGs and genomes"
    container "${container__pandas}"
    label 'mem_medium'

    input:
        tuple file(aln_tsv_gz), file(header_csv_gz)
        file geneshot_hdf

    output:
        file "genome_containment_shard.*.csv.gz" optional true
        file "genome_containment_shard.*.hdf5" optional true

    """#!/bin/bash

set -Eeuo pipefail

calculateContainment.py "${aln_tsv_gz}" "${header_csv_gz}" "${geneshot_hdf}" ${params.min_genes_per_cag}

"""

}

// Extract a table with the number of reads per gene group
// This can be used to extract the counts for CAGs, tax IDs, or any
// other grouping of genes
// Corncob takes the absolute number of reads from each sample into account
// and so it needs to have access to those integer values
process extractCounts {
    container "${container__pandas}"
    label 'mem_medium'
    
    input:
    file genome_alignments_hdf
    file details_hdf

    output:
    file "readcounts.csv.gz"
    file "abundance.csv.gz"


"""#!/bin/bash

set -Eeuo pipefail

extractCounts.py "${genome_alignments_hdf}" "${details_hdf}"



"""
}


// Format the results for each shard
process formatResults {
    tag "Use alignment information to summarize results"
    container "${container__pandas}"
    label 'mem_medium'

    input:
        tuple file(aln_tsv_gz), file(header_csv_gz)
        each file(gene_association_csv)
    
    output:
        file "genome_association_shard.*.hdf5" optional true
    
    script:
        template "formatResults.py"
}


// Collect results and combine across all shards

// Output to redis
process combineResultsRedis {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}"

    input:
        file "containment_shard.*.csv.gz"
        file "containment_shard.*.hdf5"
        file "genome.manifest.csv"
        file "geneshot.details.hdf5"
        file "geneshot.results.hdf5"
        file "geneshot.results.rdb"
    
    output:
        file "${params.output_prefix}.rdb"

"""#!/bin/bash

set -Eeuo pipefail

# Start a redis server in the background
redis-server \
    --port 6379 \
    --bind 127.0.0.1 \
    --rdbcompression yes \
    --save "" \
    --dbfilename geneshot.results.rdb \
    --dir \$PWD &

combineResults.py \
    --port 6379 \
    --host 127.0.0.1 || \
    redis-cli shutdown  # In case of failure

# Save the redis store
echo "Saving the redis store"
redis-cli save

# Shutdown the redis server
echo "Shutting down the redis server"
redis-cli shutdown

# Rename the redis DB
mv geneshot.results.rdb "${params.output_prefix}.rdb"

echo "Done"
"""

}

// Output to redis -- include data generated by corncob from the formula
process combineResultsRedisFormula {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}"

    input:
        file "containment_shard.*.csv.gz"
        file "containment_shard.*.hdf5"
        file "corncob.results.csv"
        file "genome.manifest.csv"
        file "geneshot.details.hdf5"
        file "geneshot.results.hdf5"
        file "geneshot.results.rdb"
    
    output:
        file "${params.output_prefix}.rdb"

"""#!/bin/bash

set -Eeuo pipefail

# Start a redis server in the background
redis-server \
    --port 6379 \
    --bind 127.0.0.1 \
    --rdbcompression yes \
    --save "" \
    --dbfilename geneshot.results.rdb \
    --dir \$PWD &

combineResults.py \
    --port 6379 \
    --host 127.0.0.1 || \
    redis-cli shutdown  # In case of failure

# Save the redis store
echo "Saving the redis store"
redis-cli save

# Shutdown the redis server
echo "Shutting down the redis server"
redis-cli shutdown

# Rename the redis DB
mv geneshot.results.rdb "${params.output_prefix}.rdb"

echo "Done"
"""

}

// Output to HDF
process combineResultsHDF {
    container "${container__pandas}"
    label 'mem_medium'

    input:
        file "containment_shard.*.csv.gz"
        file "containment_shard.*.hdf5"
        file "genome.manifest.csv"
        file "geneshot.details.hdf5"
        file "geneshot.results.hdf5"
        file "geneshot.results.rdb"
    
    output:
        file "${params.output_prefix}.hdf5"

"""#!/bin/bash

set -Eeuo pipefail

combineResults.py --hdf "${params.output_prefix}.hdf5"

"""

}

// Output to HDF -- including corncob results
process combineResultsHDFFormula {
    container "${container__pandas}"
    label 'mem_medium'

    input:
        file "containment_shard.*.csv.gz"
        file "containment_shard.*.hdf5"
        file "corncob.results.csv"
        file "genome.manifest.csv"
        file "geneshot.details.hdf5"
        file "geneshot.results.hdf5"
        file "geneshot.results.rdb"
    
    output:
        file "${params.output_prefix}.hdf5"

"""#!/bin/bash

set -Eeuo pipefail

combineResults.py --hdf "${params.output_prefix}.hdf5"

"""

}


// Repack an HDF5 file
process repackHDF {

    container "${container__pandas}"
    label "mem_medium"
    publishDir "${params.output_folder}"
    
    input:
    file final_hdf
        
    output:
    file "${final_hdf}"

    """
#!/bin/bash

set -Eeuxo pipefail

[ -s ${final_hdf} ]

h5repack -f GZIP=5 ${final_hdf} TEMP && mv TEMP ${final_hdf}
    """
}