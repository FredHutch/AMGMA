#!/usr/bin/env nextflow

// Set default parameters
params.help = false
params.db = null
params.geneshot_hdf = null
params.geneshot_fasta = null
params.output_hdf = null
params.min_coverage = 50
params.min_identity = 90
params.top = 50
params.fdr_method = "fdr_bh"
params.alpha = 0.2
params.filter = false
params.details = false
params.blast = false
params.no_associations = false
params.min_cag_size = 5
params.min_cag_prop = 0.25
params.window_size = 5

// Commonly used containers
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3_py38_ubuntu"
container_diamond = "quay.io/fhcrc-microbiome/docker-diamond:v0.9.31--3"
container__blast = "ncbi/blast:2.10.1"

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/AMGMA <ARGUMENTS>

Required Arguments:
--db                  AMGMA database(s) (ends with .tar), with commas delimiting multiple databases
--geneshot_hdf        Results HDF file output by GeneShot, containing CAG information
--geneshot_dmnd       DIAMOND database for the gene catalog generated by GeneShot
--output_folder       Folder to write output HDF file into
--output_hdf          Name of the output HDF file to write to the output folder

Optional Arguments:
--filter              Only analyze parameters which contain the specified substring
                      (Example: --filter BMI will only analyze the covariate matching "*BMI*", even if others are present)
--details             Include additional detailed results in output (see below)
--min_coverage        Minimum coverage required for alignment (default: 50)
--min_identity        Minimum percent identity required for alignment (default: 90)
--top                 Threshold used to retain overlapping alignments within --top% score of the max score (default: 50)
--fdr_method          Method used for FDR correction (default: fdr_bh)
--alpha               Alpha value used for FDR correction (default: 0.2)
--blast               Align with BLAST+ instead of DIAMOND
--no_associations     Exclude all analysis of CAG association metrics
--min_cag_size        Only include details for CAGs containing at least this number of genes (default: 5)
--min_cag_prop        Only include details for CAGs which align with this proportion of genes to a given genome (default: 0.25)
--window_size         Size of rolling window (for /genomes/map/<feature>/<genome_id>) (default: 5 [adjacent genes])

Output HDF:
The output from this pipeline is an HDF file which contains all of the data from the
input HDF, as well as the additional tables,

* /genomes/manifest (Table with the description of all genomes used for alignment)
* /genomes/cags/containment (Table with the number of genes aligned for all CAGs against all genomes)
* /genomes/detail/<genome_id> (Table with the coordinates of all aligned genes for a given genome)
* /genomes/annotation/<genome_id> (Table with the annotations from optional GFF input)
* /genomes/summary/<feature> (Table with summary metrics for each genome based on gene-level association with a given parameter)
* /genomes/map/<feature>/<genome_id> (Table with rolling window average of parameter association across the genome)

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.geneshot_hdf == null || params.geneshot_dmnd == null || params.db == null || params.db == null || params.output_hdf == null){
    // Invoke the function above which prints the help message
    helpMessage()

    if (params.geneshot_hdf == null){
        log.info"""
        Please provide --geneshot_hdf
        """.stripIndent()
    }
    if (params.geneshot_dmnd == null){
        log.info"""
        Please provide --geneshot_dmnd
        """.stripIndent()
    }
    if (params.db == null){
        log.info"""
        Please provide --db
        """.stripIndent()
    }
    if (params.output_hdf == null){
        log.info"""
        Please provide --output_hdf
        """.stripIndent()
    }

    // Exit out and do not run anything else
    exit 1
}

// Make sure the input database file(s) exist
db_ch = Channel.from(
    params.db.split(",")
).map { 
    fp -> file(fp) 
}

// Point to the files with the GeneShot results
geneshot_hdf = file(params.geneshot_hdf)
if ( geneshot_hdf.isEmpty() ) {
    log.info"""
    Cannot find file at ${params.geneshot_hdf}
    """
    exit 1
}
geneshot_dmnd = file(params.geneshot_dmnd)
if ( geneshot_dmnd.isEmpty() ) {
    log.info"""
    Cannot find file at ${params.geneshot_dmnd}
    """
    exit 1
}

// Unpack the database
process unpackDatabase {
    tag "Extract all files from database tarball"
    container "ubuntu:20.04"
    label "io_limited"
    errorStrategy 'retry'

    input:
        file db from db_ch
    
    output:
        file "${db}.manifest.csv" into manifest_ch
        file "*tar" into database_tar_list

"""
#!/bin/bash 

set -e

ls -lahtr

tar xvf ${db}

mv database_manifest.csv ${db}.manifest.csv

echo "Done"

"""
}

// If multiple databases were provided, join the manifests
// Make sure that there are no overlapping IDs
process validateManifest {
    
    container "${container__pandas}"
    label 'io_limited'
    errorStrategy "retry"

    input:
        file manifest_csv_list from manifest_ch.toSortedList()
    
    output:
        file "database_manifest.csv" into valid_manifest_ch
    
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


// If using BLAST+, format the database
if (params.blast) {
    // Convert the gene catalog from DMND to FASTA
    process extractFASTA {
        container "${container_diamond}"
        label "io_limited"

        input:
            file geneshot_dmnd
        
        output:
            file "ref.fasta.gz" into ref_fasta


        """#!/bin/bash

        set -e

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
            file "blastDB*" into blastDB


        """#!/bin/bash

        set -e

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
        errorStrategy 'retry'

        input:
            file database_chunk_tar from database_tar_list.flatten()
            file "*" from blastDB

        output:
            tuple file("${database_chunk_tar.name.replaceAll(/.tar/, ".aln.gz")}"), file("${database_chunk_tar.name.replaceAll(/.tar/, ".csv.gz")}") into raw_alignment_ch

    """
    #!/bin/bash

    set -e

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
        errorStrategy 'retry'

        input:
            tuple file(aln_tsv_gz), file(header_csv_gz) from raw_alignment_ch

        output:
            tuple file("${aln_tsv_gz}.filtered.tsv.gz"), file("${header_csv_gz}") into alignments_ch_1, alignments_ch_2

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

} else {

    // Align the genomes against the database with DIAMOND
    process alignGenomes {
        tag "Annotate reference genomes by alignment"
        container "${container_diamond}"
        label "mem_veryhigh"
        errorStrategy 'retry'

        input:
            file database_chunk_tar from database_tar_list.flatten()
            file geneshot_dmnd

        output:
            tuple file("${database_chunk_tar.name.replaceAll(/.tar/, ".aln.gz")}"), file("${database_chunk_tar.name.replaceAll(/.tar/, ".csv.gz")}") into alignments_ch_1, alignments_ch_2

    """
    #!/bin/bash

    set -e

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

    """
    }

}


// Calculate the containment of each CAG in each genome
process calculateContainment {
    tag "Overlap between CAGs and genomes"
    container "${container__pandas}"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        tuple file(aln_tsv_gz), file(header_csv_gz) from alignments_ch_2
        file geneshot_hdf

    output:
        file "genome_containment_shard.*.csv.gz" optional true into containment_shard_csv_gz
        file "genome_containment_shard.*.hdf5" optional true into genome_alignment_shards

"""
#!/usr/bin/env python3

import pandas as pd
import os

# In this task we will calculate the containment for each genome against each CAG

######################
# READ CAG GROUPINGS #
######################

# Read in mapping of genes to CAGs
gene_cag_map = pd.read_hdf(
    "${geneshot_hdf}", 
    "/annot/gene/all"
).set_index(
    "gene"
)
print("Read in sizes of %d CAGs containing %d genes" % (gene_cag_map["CAG"].unique().shape[0], gene_cag_map.shape[0]))

######################
# CALCULATE CAG SIZE #
######################

# Make a single dictionary with the number of genes in each CAG
cag_size = gene_cag_map["CAG"].value_counts()

#######################
# READ CONTIG HEADERS #
#######################

# Dict mapping contig names to genome IDs
print("Reading in ${header_csv_gz}")
contig_headers = pd.read_csv(
    "${header_csv_gz}"
)
assert contig_headers["contig"].unique().shape[0] == contig_headers.shape[0], "Found duplicate contig names"
contig_headers = contig_headers.set_index(
    "contig"
)["genome"]

##########################
# FORMAT GENE ALIGNMENTS #
##########################

# Read in the alignments of reference genomes against the gene catalog genes
print("Reading in ${aln_tsv_gz}")
aln_df = pd.read_csv(
    "${aln_tsv_gz}", 
    sep="\\t", 
    header=None,
    compression="gzip",
    names = [
        "contig", "gene", "pident", "length", "contig_start", "contig_end", "contig_len", "gene_start", "gene_end", "gene_len"
    ]
)

# Calculate the number of bases spanned by each alignment
aln_df = aln_df.assign(
    span = ((aln_df["contig_end"] - aln_df["contig_start"]).abs() + 1)
)

print("Adding genome labels")
aln_df = aln_df.assign(
    genome_id = aln_df["contig"].apply(contig_headers.get)
)
if aln_df["genome_id"].isnull().sum() > 0:
    print("Missing genome labels for these headers:")
    print(aln_df.loc[
        aln_df["genome_id"].isnull()
    ])
assert aln_df["genome_id"].isnull().sum() == 0

print("Read in %d gene alignments for %d genomes" % (aln_df.shape[0], aln_df["genome_id"].unique().shape[0]))


print("Adding CAG labels")
aln_df = aln_df.assign(
    CAG = aln_df["gene"].apply(
        gene_cag_map["CAG"].get
    )
)

# Function to calculate containment scores
def calc_containment(df, cag_id, n_genes_in_cag):
    # df is a DataFrame with the columns gene, contig, CAG, span, and qlen
    # n_genes_in_cag is an integer with the number of unique genes in the CAG

    # Get the total number of bases for this genome
    genome_size_bases = df.reindex(
        columns=["contig", "contig_len"]
    ).drop_duplicates(
    )["contig_len"].sum()

    # First calculate the proportion of this genome which is captured by this CAG
    cag_genome_bases = df.query("CAG == '%s'" % cag_id)["span"].sum()
    genome_prop = cag_genome_bases / genome_size_bases

    # Second calculate the proportion of the unique genes in this CAG captured in this genome
    cag_prop = df.query(
        "CAG == '%s'" % cag_id
    )[
        "gene"
    ].unique(
    ).shape[0] / float(n_genes_in_cag)

    # Return the larger of the two
    return [
        ("containment", max(genome_prop, cag_prop)),
        ("genome_prop", genome_prop),
        ("genome_bases", cag_genome_bases),
        ("cag_prop", cag_prop)
    ]

# Save all of the containment values to a single list
# (this list will be converted to a DataFrame later)
genome_containment = []

# Construct the genome containment table by iterating
# over each input genome
for genome_id, genome_df in aln_df.groupby("genome_id"):

    genome_containment.extend([
        dict([
            ("genome", genome_id),
            ("CAG", cag_id),
            ("n_genes", n_genes)
        ] + calc_containment(
            genome_df, 
            cag_id, 
            (gene_cag_map["CAG"] == cag_id).sum()
        ))
        for cag_id, n_genes in genome_df.reindex(
            columns = ["CAG", "gene"]
        ).dropna(
        ).drop_duplicates(
        )["CAG"].value_counts().items()
    ])

# Function to write out a set of alignments
def filter_alignment(genome_df, min_cag_size=5, min_cag_prop=0.25):

    # Filter genes based on CAG size
    df = genome_df.assign(
        CAG_size = genome_df["CAG"].apply(cag_size.get)
    ).query(
        "CAG_size >= %d" % min_cag_size
    ).drop(
        columns="CAG_size"
    )

    # Stop if no alignments pass this filter
    if df.shape[0] == 0:
        return

    # Calculate the number of unique genes from each CAG which was found
    cag_genes_found = df.reindex(
        columns=["gene", "CAG"]
    ).drop_duplicates(
    )[
        "CAG"
    ].value_counts()

    # Filter genes based on the proportion of each CAG which was found
    df = df.assign(
        CAG_prop = df["CAG"].apply(
            lambda cag_id: cag_genes_found[cag_id] / cag_size[cag_id]
        ).query(
            "CAG_prop >= %s" % min_cag_prop
        ).drop(
            columns="CAG_prop"
        )
    )

    # Stop if no alignments pass this filter
    if df.shape[0] == 0:
        return

    # Return the filtered table
    return df

# If no containment has been found, skip the summary step
if len(genome_containment) == 0:
    print("No matching genomes, skipping")

else:
                    
    # Write the containment table to the output HDF
    print("Making a single containment table")
    genome_containment_df = pd.DataFrame(genome_containment)
    assert genome_containment_df.shape[0] > 0, "Problem calculating containment values"
    
    genome_containment_df["CAG"] = genome_containment_df["CAG"].apply(int).apply(str)    

    genome_containment_df.to_csv(
        "genome_containment_shard.${header_csv_gz.name}", # File name will end with .csv.gz
        sep = ",",
        compression = "gzip",
        index = None
    )

    # If --details has been specified, write out the map of alignments for each genome
    if "${params.details}" != "false":

        # First compute the filtered alignment for each genome
        filtered_alignments = {
            genome_id: filter_alignment(
                genome_df,
                min_cag_size=${params.min_cag_size},
                min_cag_prop=${params.min_cag_prop},
            )
            for genome_id, genome_df in aln_df.groupby("genome_id")
        }

        # Remove the genomes which did not survive filtering
        filtered_alignments = {
            k: v
            for k, v in filtered_alignments.items()
            if v is not None
        }

        # If any remain, write to HDF
        if len(filtered_alignments) > 0:

            hdf_fp = "genome_containment_shard.${header_csv_gz.name}.hdf5"
            with pd.HDFStore(hdf_fp, "w") as store:
                for genome_id, genome_df in filtered_alignments.items():
                    genome_df.to_hdf(
                        store,
                        "/genomes/detail/%s" % genome_id
                    )

"""

}

// Join all of the results of the containment calculation
containment_shard_csv_gz_list = containment_shard_csv_gz.toSortedList()

if (params.no_associations == false){

    // Check to make sure that the input HDF has the required entries
    process parseAssociations {
        tag "Extract gene association data for the study"
        container "${container__pandas}"
        label 'mem_veryhigh'
        errorStrategy "retry"

        input:
            file geneshot_hdf
        
        output:
            file "gene_associations.*.csv.gz" into gene_association_csv_ch
        
    """#!/usr/bin/env python3

import pandas as pd
from statsmodels.stats.multitest import multipletests

store_fp = "${geneshot_hdf}"
print("Reading data from %s" % store_fp)

###################
# READ INPUT DATA #
###################

with pd.HDFStore(store_fp, "r") as store:

    for k in ["/stats/cag/corncob", "/annot/gene/all"]:
        assert k in store, "Could not find %s in %s" % (k, store_fp)

    print("Reading /stats/cag/corncob")
    corncob_df = pd.read_hdf(store, "/stats/cag/corncob")

    print("Reading /annot/gene/all")
    annot_df = pd.read_hdf(store, "/annot/gene/all")


#######################
# FORMAT CORNCOB DATA #
#######################

# Filter down to the mu estimates
corncob_df = corncob_df.reindex(
    index=corncob_df["estimate"].dropna().index
)
print("Corncob results have %d rows for mu" % (corncob_df.shape[0]))
assert corncob_df.shape[0] > 0

# Remove the intercept values
corncob_df = corncob_df.loc[
    corncob_df["parameter"] != "mu.(Intercept)"
]
print("Corncob results have %d non-intercept rows for mu" % (corncob_df.shape[0]))
assert corncob_df.shape[0] > 0

# Replace "(Intercept)" with "Intercept" in the parameter table
corncob_df.replace(
    to_replace=dict([("parameter", dict([("(Intercept)", "Intercept")]))]),
    inplace=True
)

# Reformat the corncob results as a dict
corncob_dict = dict()

# Only include parameters which match the --filter, if any:
for parameter, parameter_df in corncob_df.groupby("parameter"):
    if "${params.filter}" == "false" or "${params.filter}" in parameter:
        print("Including parameter: %s" % parameter)
        corncob_dict[parameter] = parameter_df.set_index("CAG")

    else:
        print("NOT Including parameter: %s (does not contain %s)" % (parameter, "${params.filter}"))

assert len(corncob_dict) > 0, "Did not find any parameters to include"

# Add in the FDR threshold
for parameter in corncob_dict:
    corncob_dict[
        parameter
    ][
        "${params.fdr_method}"
    ] = multipletests(
        corncob_dict[parameter]["p_value"].fillna(1),
        ${params.alpha},
        "${params.fdr_method}"
    )[1]

print("Processing %d parameters: %s" % (
    corncob_df["parameter"].unique().shape[0],
    ", ".join(corncob_df["parameter"].unique())
))

#########################
# FORMAT CAG MEMBERSHIP #
#########################

# Make sure the CAG gene membership table has values
print("CAG membership table has %d rows" % (annot_df.shape[0]))
assert annot_df.shape[0] > 0

# Make sure the expected columns exist
assert "CAG" in annot_df.columns.values and "gene" in annot_df.columns.values

# Make sure that every CAG in the corncob results has an entry in the gene membership table
cag_set = set(annot_df["CAG"].tolist())
for cag_id in corncob_df["CAG"].unique():
    assert cag_id in cag_set, "Could not find genes for CAG %s" % cag_id

######################
# FORMAT OUTPUT DATA #
######################

# For each parameter, write out a table with the association for each gene
for parameter, cag_assoc in corncob_dict.items():
    print("Processing %s" % parameter)
    
    # Make a gene-level association table
    gene_assoc = annot_df.copy()

    # Add the CAG-level associations
    for k in cag_assoc.columns.values:
        gene_assoc[k] = gene_assoc["CAG"].apply(
            cag_assoc[k].get
        )
    # Write out to CSV
    gene_assoc.to_csv(
        "gene_associations.%s.csv.gz" % parameter,
        index = None,
        compression = "gzip",
        sep = ","
    )


print("All done!")
    """
    }


    // Format the results for each shard
    process formatResults {
        tag "Use alignment information to summarize results"
        container "${container__pandas}"
        label 'mem_medium'
        errorStrategy "retry"

        input:
            tuple file(aln_tsv_gz), file(header_csv_gz) from alignments_ch_1
            each file(gene_association_csv) from gene_association_csv_ch.flatten()
        
        output:
            file "genome_association_shard.*.hdf5" into association_shard_hdf
        
    """#!/usr/bin/env python3

import gzip
import pandas as pd

##########################
# READ GENE ASSOCIATIONS #
##########################

gene_assoc_df = pd.read_csv(
    "${gene_association_csv}",
    sep = ",",
    compression = "gzip"
).set_index(
    "gene"
)


######################
# CALCULATE CAG SIZE #
######################

# Make a single dictionary with the number of genes in each CAG
cag_size = gene_assoc_df["CAG"].value_counts()


########################
# PARSE PARAMETER NAME #
########################

# Parse the parameter name from the gene association CSV file name
assert "${gene_association_csv}".startswith("gene_associations.")
assert "${gene_association_csv}".endswith(".csv.gz")
parameter_name = "${gene_association_csv}".replace(
    "gene_associations.", ""
).replace(
    ".csv.gz", ""
)

print("Analyzing parameter: %s" % (parameter_name))


#######################
# READ CONTIG HEADERS #
#######################

# Dict mapping contig names to genome IDs
print("Reading in ${header_csv_gz}")
contig_headers = pd.read_csv(
    "${header_csv_gz}"
)
assert contig_headers["contig"].unique().shape[0] == contig_headers.shape[0], "Found duplicate contig names"
contig_headers = contig_headers.set_index(
    "contig"
)["genome"]

##########################
# FORMAT GENE ALIGNMENTS #
##########################

# Read in the alignments of reference genomes against the gene catalog genes
print("Reading in ${aln_tsv_gz}")
aln_df = pd.read_csv(
    "${aln_tsv_gz}", 
    sep="\\t", 
    header=None,
    compression="gzip",
    names = [
        "contig", "gene", "pident", "length", "contig_start", "contig_end", "contig_len", "gene_start", "gene_end", "gene_len"
    ]
)

print("Adding genome labels")
aln_df = aln_df.assign(
    genome_id = aln_df["contig"].apply(contig_headers.get)
)
if aln_df["genome_id"].isnull().sum() > 0:
    print("Missing genome labels for these headers:")
    print(aln_df.loc[
        aln_df["genome_id"].isnull()
    ])
assert aln_df["genome_id"].isnull().sum() == 0

print("Read in %d gene alignments for %d genomes" % (aln_df.shape[0], aln_df["genome_id"].unique().shape[0]))


#####################
# FILTER ALIGNMENTS #
#####################

# Function to filter a set of alignments
def filter_alignments(genome_aln_df):

    # Filter genes based on CAG size
    genome_aln_df = genome_aln_df.assign(
        CAG_size = genome_aln_df["CAG"].apply(cag_size.get)
    ).query(
        "CAG_size >= %d" % min_cag_size
    ).drop(
        columns="CAG_size"
    )

    # Stop if no alignments pass this filter
    if genome_aln_df_fdr.shape[0] == 0:
        return

    # Calculate the number of unique genes from each CAG which was found
    cag_genes_found = genome_aln_df.reset_index(
    ).reindex(
        columns=["gene", "CAG"]
    ).drop_duplicates(
    )[
        "CAG"
    ].value_counts()

    # Filter genes based on the proportion of each CAG which was found
    genome_aln_df = genome_aln_df.assign(
        CAG_prop = genome_aln_df["CAG"].apply(
            lambda cag_id: cag_genes_found[cag_id] / cag_size[cag_id]
        ).query(
            "CAG_prop >= %s" % min_cag_prop
        ).drop(
            columns="CAG_prop"
        )
    )

    # Stop if no alignments pass this filter
    if genome_aln_df.shape[0] == 0:
        return
    
    # Add in the gene annotations
    for k in gene_assoc_df.columns.values:

        # To annotate the genome, figure out which of the catalog genes
        # each of the genes in the genome is most similar to, and then
        # fill in the value of the CAG which that catalog gene is a part of

        genome_aln_df[k] = genome_aln_df[
            "gene"
        ].apply(
            gene_assoc_df[k].get
        )

    # Make sure that we have the Wald metric
    if "wald" not in genome_aln_df_fdr.columns.values
    genome_aln_df = genome_aln_df.assign(
        wald = genome_aln_df["estimate"] / genome_aln_df["std_error"]
    )

    # Return the filtered alignments
    return genome_aln_df

# Filter and annotate each genome
aln_df = [
    filter_alignments(genome_df)
    for genome_id, genome_df in aln_df.groupby("genome_id")
]

# Remove null values
aln_df = [v for v in aln_df if v is not None]

# Make a DataFrame
aln_df = pd.DataFrame(aln_df)


###################
# ANALYZE GENOMES #
###################

# Open a connection to the HDF store used for all output information
output_store = pd.HDFStore("genome_association_shard.%s.${header_csv_gz.name.replaceAll(/.csv.gz/, "")}.hdf5" % parameter_name, "w")

# Function to summarize a single genome
def process_genome(genome_id, genome_aln_df):

    # Get the table which passes the FDR filter
    genome_aln_df_fdr = genome_aln_df.loc[
        genome_aln_df["${params.fdr_method}"] <= ${params.alpha}
    ]

    print("%d / %d genes pass the CAG-level FDR threshold" % 
        (genome_aln_df_fdr.shape[0], genome_aln_df.shape[0]))

    # Stop if no alignments pass this filter
    if genome_aln_df_fdr.shape[0] == 0:
        return

    # Return the summary metrics
    return dict([
        ("genome_id", genome_id),
        ("parameter", parameter_name),
        ("total_genes", genome_aln_df.shape[0]),
        ("n_pass_fdr", genome_aln_df_fdr.shape[0]),
        ("prop_pass_fdr", genome_aln_df_fdr.shape[0] / float(genome_aln_df.shape[0])),
        ("mean_est_coef", genome_aln_df_fdr["estimate"].mean())
    ])

# Iterate over every genome and process it
summary_df = [
    process_genome(genome_id, genome_df)
    for genome_id, genome_df in aln_df.groupby("genome_id")
]

summary_df = [v for v in summary_df if v is not None]

# Check that we have any data to save
if len(summary_df) > 0:

    # Make a DataFrame
    summary_df = pd.DataFrame(summary_df)

    # Save to the HDF
    summary_df.to_hdf(
        output_store,
        "/genomes/summary/%s" % parameter_name,
        format = "fixed"
    )


########################
# WRITE OUT GENOME MAP #
########################

# Compute a rolling window for the Wald metric across the genome
for genome_id, genome_df in aln_df.groupy("genome_id"):

    # Format the rolling window across every contig
    pd.concat([
        contig_df.sort_values(
            by="contig_start"
        ).reindex(
            columns=["contig_start", "wald", "estimate", "p_value"]
        ).rolling(
            ${params.window_size},
        ).median(
        ).dropna(
        ).assign(
            contig=contig_id,
            genome=genome_id,
        )
        for contig_id, contig_df in genome_df.groupby("contig")
    ]).to_hdf(
        output_store,
        "/genomes/map/%s/%s" % (parameter_name, genome_id)
    )


######################
# CLOSE OUTPUT FILES #
######################

output_store.close()
    """
    }

    association_shard_hdf_list = association_shard_hdf.toSortedList()

    // Collect results and combine across all shards
    process combineResults {
        tag "Make a single output HDF"
        container "${container__pandas}"
        label 'mem_veryhigh'
        errorStrategy "retry"

        input:
            file containment_shard_csv_gz_list
            file association_shard_hdf_list
            file geneshot_hdf
            file manifest_csv from valid_manifest_ch
        
        output:
            file "${params.output_hdf}" into final_hdf
        
    """#!/usr/bin/env python3

import pandas as pd
import os
import shutil

# Read in and combine all of the containment tables
containment_csv_list = "${containment_shard_csv_gz_list}".split(" ")
print("Reading in containment values from %d files" % len(containment_csv_list))
containment_df = pd.concat([
    pd.read_csv(
        fp,
        sep = ",",
        compression = "gzip"
    )
    for fp in containment_csv_list
])
print("Read in containment for %d genomes and %d CAGs" % (containment_df["genome"].unique().shape[0], containment_df["CAG"].unique().shape[0]))

# Rename the geneshot output HDF5 as the output HDF5
# All of the results will be appended to this object
# which will retain all of the formatting of the original
shutil.copyfile("${geneshot_hdf}", "${params.output_hdf}")

# Open a connection to the output file
output_store = pd.HDFStore("${params.output_hdf}", "a")

# Write out the genome manifest to the final HDF5
print("Writing out the manifest to HDF")
pd.read_csv(
    "${manifest_csv}"
).drop(
    columns = "uri"
).to_hdf(
    output_store,
    "/genomes/manifest"
)

# Write out the combined containment table
print("Writing out the containment to HDF")
containment_df.to_hdf(
    output_store,
    "/genomes/cags/containment",
    format = "table",
    data_columns = ["genome", "CAG"]
)

# Keep track of all of the parameter summary information
parameter_summaries = dict()

# Iterate over each of the parameter association shards
for hdf_fp in "${association_shard_hdf_list}".split(" "):

    # Open a connection to the input
    with pd.HDFStore(hdf_fp, "r") as store:

        # Iterate over every element in the HDF5
        for k in store:

            # Collection genome summary information
            if k.startswith("/genomes/summary/"):
                parameter_name = k.split("/")[-1]
                print("Collecting summary information for %s" % parameter_name)
                if parameter_name not in parameter_summaries:
                    parameter_summaries[parameter_name] = []

                parameter_summaries[parameter_name].append(
                    pd.read_hdf(
                        store,
                        k
                    )
                )

            # Write genome detailed information
            elif k.startswith(("/genomes/detail/", "/genomes/map/")):
                print("Copying %s to output HDF" % k)

                pd.read_hdf(
                    store,
                    k
                ).to_hdf(
                    output_store,
                    k,
                    format = "fixed",
                    complevel = 5
                )

            else:
                assert False, "Didn't expect %s" % k


# Now write out all of the summary information for each parameter
for parameter_name, genome_summary_list in parameter_summaries.items():

    # Collapse all of the results into a single table
    parameter_df = pd.concat(genome_summary_list)

    print("Writing out details on %d genomes for %s" % (parameter_df.shape[0], parameter_name))

    # Write out to the HDF
    parameter_df.to_hdf(
        output_store,
        "/genomes/summary/%s" % parameter_name,
        format = "fixed",
        complevel = 5
    )

output_store.close()

    """
    }
} else {
    // Collect results and combine across all shards
    process combineResultsNoAssoc {
        tag "Make a single output HDF"
        container "${container__pandas}"
        label 'mem_veryhigh'
        errorStrategy "retry"

        input:
            file containment_shard_csv_gz_list
            file geneshot_hdf
            file manifest_csv from valid_manifest_ch
            file "genome_alignments.*.hdf5" from genome_alignment_shards.toSortedList()
        
        output:
            file "${params.output_hdf}" into final_hdf
        
    """#!/usr/bin/env python3

import pandas as pd
import os
import shutil

# Read in and combine all of the containment tables
containment_csv_list = "${containment_shard_csv_gz_list}".split(" ")
print("Reading in containment values from %d files" % len(containment_csv_list))
containment_df = pd.concat([
    pd.read_csv(
        fp,
        sep = ",",
        compression = "gzip"
    )
    for fp in containment_csv_list
])
print("Read in containment for %d genomes and %d CAGs" % (containment_df["genome"].unique().shape[0], containment_df["CAG"].unique().shape[0]))

# Rename the geneshot output HDF5 as the output HDF5
# All of the results will be appended to this object
# which will retain all of the formatting of the original
shutil.copyfile("${geneshot_hdf}", "${params.output_hdf}")

# Open a connection to the output file
output_store = pd.HDFStore("${params.output_hdf}", "a")

# Write out the genome manifest to the final HDF5
print("Writing out the manifest to HDF")
pd.read_csv(
    "${manifest_csv}"
).drop(
    columns = "uri"
).to_hdf(
    output_store,
    "/genomes/manifest"
)

# Write out the combined containment table
print("Writing out the containment to HDF")
containment_df.to_hdf(
    output_store,
    "/genomes/cags/containment",
    format = "table",
    data_columns = ["genome", "CAG"]
)

# Write out the genome maps, if there are any
for fp in os.listdir("."):

    # Parse any files with matching filenames
    if fp.startswith("genome_alignments.") and fp.endswith(".hdf5"):

        # Open up the file
        with pd.HDFStore(fp, "r") as input_store:

            # Iterate over every table
            for k in input_store:

                # Read in the table
                df = pd.read_hdf(input_store, k)

                # Write to the output HDF
                print("Writing alignment details for %s" % k)
                df.to_hdf(output_store, k)

output_store.close()

    """
    }
}



// Repack an HDF5 file
process repackHDF {

    container "${container__pandas}"
    tag "Compress HDF store"
    label "mem_medium"
    errorStrategy "retry"
    publishDir "${params.output_folder}"
    
    input:
    file final_hdf
        
    output:
    file "${final_hdf}"

    """
#!/bin/bash

set -e

[ -s ${final_hdf} ]

h5repack -f GZIP=5 ${final_hdf} TEMP && mv TEMP ${final_hdf}
    """
}