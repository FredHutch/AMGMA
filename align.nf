#!/usr/bin/env nextflow

// Set default parameters
params.help = false
params.geneshot_hdf = null
params.geneshot_fasta = null
params.db_dmnd = null
params.db_hdf = null
params.output_hdf = null
params.min_coverage = 50
params.min_identity = 80
params.batchsize = 1000
params.fdr_method = "fdr_bh"
params.alpha = 0.2

// Commonly used containers
container__pandas = "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/AMGMA <ARGUMENTS>
    
    Required Arguments:
      --geneshot_hdf        Results HDF file output by GeneShot, containing CAG information
      --geneshot_dmnd       DIAMOND database for the gene catalog generated by GeneShot
      --db_dmnd             Microbial genome database, DMND format
      --db_hdf              Microbial genome database, HDF format
      --output_folder       Folder to write output HDF file within
      --output_hdf          Name of the output HDF file to write to the output folder

    Optional Arguments:
      --min_coverage        Minimum coverage required for alignment (default: 50)
      --min_identity        Minimum percent identity required for alignment (default: 80)
      --batchsize           Number of genomes to process in a given batch
      --fdr_method          Method used for FDR correction (default: fdr_bh)
      --alpha               Alpha value used for FDR correction (default: 0.2)

    Database files:
    The DMND and HDF files used as inputs for this script are produced by the build_db.nf
    function provided as part of this project.

    Input HDF:
    The information on CAGs used as an input for this analysis is formatted as the output
    from the GeneShot pipeline (more information at github.com/Golob-Minot/GeneShot).

    Input FASTA:
    Another output from the GeneShot pipeline is the set of genes which were identified,
    in gzip-compressed FASTA format, specified as input here with --geneshot_fasta

    Output HDF:
    The output from this pipeline is an HDF file which contains all of the data from the
    input HDF, as well as the additional tables, 

      * /genomes/summary/<feature>
      * /genomes/deaeil/<feature>/<genome_id>

      for each <feature> tested in the input, and for each <genome_id> in the database

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.geneshot_hdf == null || params.geneshot_dmnd == null || params.db_dmnd == null || params.db_hdf == null || params.output_hdf == null){
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
    if (params.db_dmnd == null){
        log.info"""
        Please provide --db_dmnd
        """.stripIndent()
    }
    if (params.db_hdf == null){
        log.info"""
        Please provide --db_hdf
        """.stripIndent()
    }
    if (params.output_hdf == null){
        log.info"""
        Please provide --output_hdf
        """.stripIndent()
    }

    // Exit out and do not run anything else
    exit 0
}

// Point to the file with the GeneShot results
geneshot_hdf = file(params.geneshot_hdf)

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
    
"""
#!/usr/bin/env python3

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
corncob_df = corncob_df.loc[
    corncob_df["parameter"].apply(lambda s: s.startswith("mu."))
]
print("Corncob results have %d rows for mu" % (corncob_df.shape[0]))
assert corncob_df.shape[0] > 0

# Remove the intercept values
corncob_df = corncob_df.loc[
    corncob_df["parameter"] != "mu.(Intercept)"
]
print("Corncob results have %d non-intercept rows for mu" % (corncob_df.shape[0]))
assert corncob_df.shape[0] > 0

# Remove the "mu." from the parameter
corncob_df["parameter"] = corncob_df["parameter"].apply(lambda s: s[3:])

# Reformat the corncob results as a dict
corncob_dict = dict([
    (parameter, parameter_df.pivot_table(index="CAG", columns="type", values="value"))
    for parameter, parameter_df in corncob_df.groupby("parameter")
])

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


// Align genes from the database against genes from the CAGs with DIAMOND
process alignGenes {
    tag "Align genes against reference database"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file geneshot_dmnd from file(params.geneshot_dmnd)
    file db_dmnd from file(params.db_dmnd)

    output:
    file "genes.DB.aln.gz" into aln_tsv

    """
#!/bin/bash

set -e

# Extract the FASTA from the geneshot DMND
echo "Extracting geneshot FASTA from ${geneshot_dmnd}"
diamond getseq -d "${geneshot_dmnd}" > query.fasta

echo "Aligning sequences"
diamond \
    blastp \
    --query query.fasta \
    --db "${db_dmnd}" \
    --threads ${task.cpus} \
    --out genes.DB.aln.gz \
    --threads ${task.cpus} \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --query-cover ${params.min_coverage} \
    --id ${params.min_identity} \
    --block-size ${task.memory.toMega() / (1024 * 6)} \
    --compress 1 \
    --unal 0

    """

}


// Split up the genomes into shards for parallel processing
process splitGenomes {
    tag "Group genomes for parallel processing"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy "retry"

    input:
        file db_hdf from file(params.db_hdf)

    output:
        file "genome_list.*.txt.gz" into genome_list_ch

"""
#!/usr/bin/env python3

import gzip
import pandas as pd

manifest = pd.read_hdf("${db_hdf}", "/manifest")

for ix, genome_list in enumerate([
    manifest["id"].values[ix: (ix + ${params.batchsize})]
    for ix in range(0, manifest.shape[0], ${params.batchsize})
]):
    with gzip.open("genome_list.%d.txt.gz" % ix, "wt") as fo:
        fo.write("\\n".join(genome_list))

"""

}


// Format the results for each shard
process formatResults {
    tag "Use alignment information to summarize results"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy "retry"

    input:
        file genome_list from genome_list_ch.flatten()
        each file(gene_association_csv) from gene_association_csv_ch
        file aln_tsv
        file db_hdf from file(params.db_hdf)
    
    output:
        tuple val(genome_list.name), file("genome_analysis_shard.hdf5") into shard_hdf
    
"""
#!/usr/bin/env python3

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


##########################
# FORMAT GENE ALIGNMENTS #
##########################

# Read in the alignments of reference genome genes against the gene catalog genes
# Reformat to a dict linking every reference genome gene to a single entry in the gene catalog
aln_df = pd.read_csv(
    "${aln_tsv}", 
    sep="\\t", 
    header=None
).reindex(
    columns=[0, 1, 11]
).rename(columns=dict([
    (0, "catalog_gene"),
    (1, "genome_centroid"),
    (11, "score")
])).sort_values(
    by = "score",
    ascending = False
).drop(
    columns = "score"
).groupby(
    "genome_centroid"
).head(
    1
).set_index(
    "genome_centroid"
)

print("Read in %d gene alignments" % aln_df.shape[0])

###########################
# CHECK INPUT CONSISTENCY #
###########################

# To make sure that the set of gene names provided in the gene FASTA (and therefore
# used for the gene alignments) match the names of the genes in the HDF, check for
# overlap between the two

cagalog_genes_from_fasta = set(aln_df["catalog_gene"].tolist())
print("Catalog genes from provided FASTA: %d" % len(cagalog_genes_from_fasta))

catalog_genes_from_hdf = set(gene_assoc_df.index.values)
print("Catalog genes from provided HDF: %d" % len(catalog_genes_from_hdf))

catalog_genes_overlap = cagalog_genes_from_fasta & catalog_genes_from_hdf
print("Overlap: %d" % len(catalog_genes_overlap))

assert len(catalog_genes_overlap) > 0


##################
# SUBSET GENOMES #
##################

# Read in the group of genomes to process in this shard
genome_list = [
    line.rstrip("\\n")
    for line in gzip.open("${genome_list}", "rt")
]
print("Processing %d genomes" % len(genome_list))

# Open a connection to the HDF store with genome alignment information
genome_store = pd.HDFStore("${db_hdf}", "r")

# Read the table linking `genome_gene` to `genome_centroid`
print("Reading in table linking genome genes to centroids")
genome_centroid_df = pd.read_hdf(
    genome_store, 
    "/genome_centroids"
).set_index("genome_gene")
print("Read in table linking %d genes to %d centroids" % 
    (
        genome_centroid_df.shape[0],
        genome_centroid_df["genome_centroid"].unique().shape[0]
    ))


# Open a connection to the HDF store used for all output information
output_store = pd.HDFStore("genome_analysis_shard.hdf5", "w")

# Function to process a single genome
def process_genome(genome_id):

    # Read in the alignment to that genome
    genome_aln_df = pd.read_hdf(
        genome_store, 
        "/genomes/%s" % genome_id
    )

    # Add in the centroid label
    genome_aln_df["genome_centroid"] = genome_aln_df["gene_id"].apply(
        genome_centroid_df["genome_centroid"].get
    )
    assert genome_aln_df["genome_centroid"].isnull().sum() == 0, genome_aln_df.loc[genome_aln_df["genome_centroid"].isnull()]

    # Add in the 'catalog' gene label
    genome_aln_df["catalog_gene"] = genome_aln_df["genome_centroid"].apply(
        aln_df["catalog_gene"].get
    )
    if genome_aln_df["catalog_gene"].dropna().shape[0] == 0:
        print("Warning: %s has 0 genes which match the catalog" % genome_id)

    # Add in the gene annotations
    for k in gene_assoc_df.columns.values:

        # To annotate the genome, figure out which of the catalog genes
        # each of the genes in the genome is most similar to, and then
        # fill in the value of the CAG which that catalog gene is a part of

        genome_aln_df[k] = genome_aln_df[
            "catalog_gene"
        ].apply(
            gene_assoc_df[k].get
        )

    # Write out the full table
    key = "/genomes/detail/%s/%s" % (parameter_name, genome_id)
    print("Writing out to %s" % key)
    
    genome_aln_df.reindex(
        columns = [
            "start",
            "end",
            "strand",
            "gc_cont",
            "catalog_gene",
            "CAG",
            "estimate",
            "p_value",
            "std_error",
            "${params.fdr_method}"
        ]
    ).to_hdf(
        output_store,
        key,
        format = "fixed",
        complevel = 5
    )

    # Get the table which passes the FDR filter
    genome_aln_df_fdr = genome_aln_df.loc[
        genome_aln_df["${params.fdr_method}"] <= ${params.alpha}
    ]

    print("%d / %d genes pass the CAG-level FDR threshold" % 
        (genome_aln_df_fdr.shape[0], genome_aln_df.shape[0]))

    # Return the summary metrics
    return dict([
        ("genome_id", genome_id),
        ("parameter", parameter_name),
        ("total_genes", genome_aln_df.shape[0]),
        ("n_pass_fdr", genome_aln_df_fdr.shape[0]),
        ("prop_pass_fdr", genome_aln_df_fdr.shape[0] / float(genome_aln_df.shape[0])),
        ("mean_est_coef", genome_aln_df_fdr["estimate"].mean())
    ])

# Iterate over every genome and process it, saving the summary to the HDF
pd.DataFrame([
    process_genome(genome_id)
    for genome_id in genome_list
]).to_hdf(
    output_store,
    "/genomes/summary/%s" % parameter_name,
    format = "fixed"
)


######################
# CLOSE OUTPUT FILES #
######################

genome_store.close()
output_store.close()

"""
}


// Calculate the containment of each CAG in each genome
process calculateContainment {
    tag "Overlap between CAGs and genomes"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
        tuple val(genome_group_label), file("genome_analysis_shard.*.hdf5") from shard_hdf.groupTuple()
        file geneshot_hdf

    output:
        file "genome_analysis_shard.hdf5" into containment_shard_hdf

"""
#!/usr/bin/env python3

import pandas as pd
import os

# There will be multiple input HDF5 files, each from the same group of genomes
# but analyzed against >=1 parameter associations

# In this task we will calculate the containment for each genome against each CAG

# Read in the size of each CAG
cag_size = pd.read_hdf(
    "${geneshot_hdf}", 
    "/annot/cag/all"
).set_index(
    "CAG"
)["size"]
print("Read in sizes of %d CAGs containing %d genes" % (cag_size.shape[0], cag_size.sum()))

# Get the names of all of the input files available
input_hdf_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("genome_analysis_shard")
]
for fp in input_hdf_list:
    assert fp.endswith(".hdf5")

print("Processing input HDF files:")
print("\\n".join(input_hdf_list))

# Keep track of which genomes we have calculated containment for
genome_containment = dict()

# Copy all of the information in the input files to the output
with pd.HDFStore("genome_analysis_shard.hdf5", "w") as output_store:
    # Iterate over each input HDF
    for input_hdf in input_hdf_list:

        # Open a connection to the HDF
        with pd.HDFStore(input_hdf, "r") as input_store:
            
            # Iterate over each key
            for k in input_store:

                # Copy to the output HDF
                pd.read_hdf(
                    input_store,
                    k
                ).to_hdf(
                    output_store,
                    k,
                    format = "fixed"
                )

                # Check to see if this is a genome alignment
                if k.startswith("/genomes/detail/"):
                    # Parse the genome name
                    genome_name = k.split("/")[-1]

                    # If we've already processed this genome, move on
                    if genome_name in genome_containment:
                        continue

                    print("Processing containment for %s" % genome_name)

                    df = pd.read_hdf(
                        input_store,
                        k
                    )

                    genome_containment[
                        genome_name
                    ] = pd.DataFrame([
                        dict([
                            ("genome", genome_name),
                            ("CAG", cag_id),
                            ("n_genes", n_genes),
                            ("containment", n_genes / min(cag_size[cag_id], df.shape[0]))
                        ])
                        for cag_id, n_genes in df["CAG"].value_counts().items()
                    ])
                    
    # Write the containment table to the output HDF
    genome_containment_df = pd.concat(list(genome_containment.values()))
    genome_containment_df["CAG"] = genome_containment_df["CAG"].apply(int).apply(str)
    
    assert genome_containment_df.shape[0] > 0, "Problem calculating containment values"

    genome_containment_df.to_hdf(
        output_store,
        "/genomes/cags/containment",
        format = "fixed",
        complevel = 5
    )
"""

}


// Collect results and combine across all shard
process combineResults {
    tag "Make a single output HDF"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy "retry"

    input:
        file "shard.*.hdf" from containment_shard_hdf.collect()
        file geneshot_hdf
    
    output:
        file "${params.output_hdf}" into final_hdf
    
"""
#!/usr/bin/env python3

import pandas as pd
import os

# Make sure that all of the input files are present
hdf_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("shard.")
]
for fp in hdf_list:
    assert fp.endswith(".hdf"), fp

# Keep track of all of the parameter summary information
parameter_summaries = dict()

# Make a list of all of the containment tables
containment_df = list()

# Open a connection to the output file
output_store = pd.HDFStore("${params.output_hdf}", "w")

# Only copy over the manifest once
wrote_manifest = False

# Iterate over the HDF from each shard
for hdf_fp in hdf_list:
    with pd.HDFStore(hdf_fp, "r") as store:
        for k in store:

            # Write the manifest
            if k == "/manifest":
                print("Copying %s to output HDF" % k)
                if wrote_manifest is False:
                    pd.read_hdf(
                        store,
                        k
                    ).to_hdf(
                        output_store,
                        k,
                        format = "fixed"
                    )
                    wrote_manifest = True

            # Collection genome summary information
            elif k.startswith("/genomes/summary/"):
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
            elif k.startswith("/genomes/detail/"):
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

            # Save the genome containment information
            elif k == "/genomes/cags/containment":
                containment_df.append(
                    pd.read_hdf(
                        store,
                        k
                    )
                )
            
            else:
                assert False, "Didn't expect %s" % k

# Write out a single table with containment values
assert len(containment_df) > 0
containment_df = pd.concat(containment_df)

containment_df.to_hdf(
    output_store,
    "/genomes/cags/containment",
    format = "table",
    data_columns = ["genome", "CAG"]
)

# Now write out all of the summary information for each parameter
for parameter_name, genome_summary_list in parameter_summaries.items():

    # Collapse all of the results into a single table
    parameter_df = pd.concat(genome_summary_list)

    print("Writing out details on %d genomes for %s" % (parameter_df.shape[0], parameter_name))

    print(parameter_df)

    # Write out to the HDF
    parameter_df.to_hdf(
        output_store,
        "/genomes/summary/%s" % parameter_name,
        format = "fixed",
        complevel = 5
    )

# Copy everything from the input HDF5
input_store = pd.HDFStore("${geneshot_hdf}", "r")

for k in input_store:
    if k == "/abund/cag/wide" or k == "/annot/gene/all":
        pd.read_hdf(
            input_store,
            k
        ).to_hdf(
            output_store,
            k,
            format = "table",
            data_columns = ["CAG"]
        )
    else:
        pd.read_hdf(
            input_store,
            k
        ).to_hdf(
            output_store,
            k,
            format = "fixed"
        )

input_store.close()

output_store.close()

"""
}


// Repack an HDF5 file
process repackHDF {

    container "${container__pandas}"
    tag "Compress HDF store"
    label "mem_veryhigh"
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