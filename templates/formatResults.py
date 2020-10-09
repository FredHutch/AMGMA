#!/usr/bin/env python3

import gzip
import pandas as pd
import pickle

pickle.HIGHEST_PROTOCOL = 4

##########################
# READ GENE ASSOCIATIONS #
##########################

print("Reading in ${gene_association_csv}")
gene_assoc_df = pd.read_csv(
    "${gene_association_csv}",
    sep = ",",
    compression = "gzip",
    dtype = dict(CAG=int)
).set_index(
    "gene"
)

######################
# CALCULATE CAG SIZE #
######################

# Make a single dictionary with the number of genes in each CAG
cag_size = gene_assoc_df["CAG"].apply(int).value_counts()

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

print("Read in %d alignments" % aln_df.shape[0])

# Stop if there are no alignments
if aln_df.shape[0] == 0:
    print("There are not alignments -- stopping")
    exit()

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

print("Adding CAG labels")
aln_df = aln_df.assign(
    CAG = aln_df["gene"].apply(gene_assoc_df["CAG"].get)
)

# Remove any genes which don't have a CAG label
print("%d / %d genes have a valid CAG label" % (aln_df["CAG"].dropna().shape[0], aln_df.shape[0]))
assert aln_df["CAG"].isnull().sum() < aln_df.shape[0]
aln_df = aln_df.assign(
    has_CAG = aln_df["CAG"].apply(lambda cag_id: pd.isnull(cag_id) is False)
).query(
    "has_CAG"
).drop(
    columns="has_CAG"
)

# Format the CAG as an integer
aln_df["CAG"] = aln_df["CAG"].apply(int)

print("Read in %d gene alignments for %d genomes" % (aln_df.shape[0], aln_df["genome_id"].unique().shape[0]))


#####################
# FILTER ALIGNMENTS #
#####################

# Function to filter a set of alignments
def filter_alignments(genome_aln_df):

    # Filter genes based on CAG size
    genome_aln_df = genome_aln_df.assign(
        CAG_size = genome_aln_df["CAG"].apply(
            int
        ).apply(
            cag_size.get
        )
    ).query(
        "CAG_size >= %d" % ${params.min_cag_size}
    ).drop(
        columns="CAG_size"
    )

    # Stop if no alignments pass this filter
    if genome_aln_df.shape[0] == 0:
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
        )
    ).query(
        "CAG_prop >= %s" % ${params.min_cag_prop}
    ).drop(
        columns="CAG_prop"
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
    if "wald" not in genome_aln_df.columns.values:
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
if len(aln_df) > 0:
    aln_df = pd.concat(aln_df)
else:
    aln_df = None


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

    # Return the summary metrics
    return dict([
        ("genome_id", genome_id),
        ("parameter", parameter_name),
        ("total_genes", genome_aln_df.shape[0]),
        ("n_pass_fdr", genome_aln_df_fdr.shape[0]),
        ("prop_pass_fdr", genome_aln_df_fdr.shape[0] / float(genome_aln_df.shape[0])),
        ("mean_est_coef", genome_aln_df["estimate"].mean()),
        ("mean_wald", genome_aln_df["wald"].mean()),
    ])

# Make sure that we have any alignments to summarize
if aln_df is not None:

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

# Make sure that we have any alignments to summarize
if aln_df is not None:

    # Compute a rolling window for the Wald metric across the genome
    for genome_id, genome_df in aln_df.groupby("genome_id"):

        # Format the rolling window across every contig
        map_df = pd.concat([
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
            )
            for contig_id, contig_df in genome_df.groupby("contig")
        ])
        
        # Only write maps that aren't empty
        if map_df.shape[0] > 0:
            map_df.to_hdf(
                output_store,
                "/genomes/map/%s/%s" % (parameter_name, genome_id)
            )


######################
# CLOSE OUTPUT FILES #
######################

output_store.close()
