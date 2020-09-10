#!/usr/bin/env python3

import pandas as pd
import os
import pickle

pickle.HIGHEST_PROTOCOL = 4

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

print("Read in %d alignments" % aln_df.shape[0])

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
    CAG = aln_df["gene"].apply(gene_cag_map["CAG"].get)
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
        )
    ).query(
        "CAG_prop >= %s" % ${params.min_cag_prop}
    ).drop(
        columns="CAG_prop"
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

                    # Only save the most essential data
                    genome_df.reindex(
                        columns = [
                            "contig",
                            "gene",
                            "pident",
                            "contig_start",
                            "contig_end",
                            "contig_len",
                            "genome_id",
                            "CAG",
                        ]
                    ).to_hdf(
                        store,
                        "/genomes/detail/%s" % genome_id
                    )
