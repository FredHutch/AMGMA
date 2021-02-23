#!/usr/bin/env python3

from collections import defaultdict
import gzip
import json
import logging
import os
import pandas as pd
import sys

genome_alignments_hdf = sys.argv[1]
details_hdf = sys.argv[2]

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [Extract Genome Counts] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Populate a dict with the set of genes which align to each genome
genome_membership = dict()

# Open the HDF which contains genome alignments
logging.info(f"Reading in genome alignments from {genome_alignments_hdf}")
with pd.HDFStore(genome_alignments_hdf, "r") as store:

    # Iterate over every key in the store
    for key in store.keys():

        # If the key is to a genome alignment
        if key.startswith("/genomes/detail/"):

            #Parse the genome ID from the key
            genome_id = key[len("/genomes/detail/"):]

            # Add the set of aligned genes
            genome_membership[
                genome_id
            ] = set(
                pd.read_hdf(
                    store, 
                    key
                )['gene'].tolist()
            )

logging.info(f"Read in alignments from {len(genome_membership):,} genomes")

# Make an object to hold the number of reads per genome, per sample
genome_counts = defaultdict(lambda: dict())

# For each specimen, record the total number of reads
specimen_readcounts = dict()

# Read in each of the FAMLI output objects from the details HDF

logging.info(f"Reading in gene abundances from {details_hdf}")

# Open the HDF store
with pd.HDFStore(details_hdf, "r") as store:

    # Iterate over all keys
    for hdf_key in store.keys():

        # If the key holds detailed abundances
        if hdf_key.startswith("/abund/gene/long/"):

            # Parse the specimen name
            specimen = hdf_key[len('/abund/gene/long/'):]

            logging.info(
                "Reading in %s" % specimen
            )

            # Read the complete table of abundances
            df = pd.read_hdf(store, hdf_key).set_index("id")

            # Add the total readcounts for this specimen
            specimen_readcounts[specimen] = df['nreads'].sum()

            # For each genome
            for genome_id, gene_set in genome_membership.items():

                # Add the number of reads aligned to the genes from this genome
                genome_counts[
                    genome_id
                ][
                    specimen
                ] = df.reindex(
                    index=list(gene_set)
                )[
                    'nreads'
                ].fillna(
                    0
                ).sum()

# Make a DataFrame with the number of reads per genome
genome_counts = pd.DataFrame(
    genome_counts
).fillna(
    0
)

# Add the TOTAL column
genome_counts = genome_counts.assign(
    total = pd.Series(
        specimen_readcounts
    ).reindex(
        index=genome_counts.index
    )
).fillna(
    0
).applymap(
    int
)
assert genome_counts["total"].sum() > 0

# Save the CSV of counts
genome_counts.reset_index(
).rename(
    columns=dict(index="specimen")
).to_csv(
    "readcounts.csv.gz",
    compression="gzip",
    index=None
)

# Divide by the total to save the CSV of relative abundances
genome_counts.applymap(
    float
).apply(
    lambda r: r / r['total'],
    axis=1
).drop(
    columns='total'
).reset_index(
).rename(
    columns=dict(index="specimen")
).to_csv(
    "abundance.csv.gz",
    compression="gzip",
    index=None,
)
