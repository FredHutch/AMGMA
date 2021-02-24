#!/usr/bin/env python3

import argparse
from collections import defaultdict
from direct_redis import DirectRedis
from functools import lru_cache
import h5py
import logging
import numpy as np
import pandas as pd
import pickle
import os
from scipy.optimize import nnls
from statsmodels.stats.multitest import multipletests
import sys

pickle.HIGHEST_PROTOCOL = 4


###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description='Collect all AMGMA results'
)

# Add the arguments
parser.add_argument(
    'output_prefix',
    type=str,
    help='Prefix to append to all outputs'
)
parser.add_argument(
    '--host',
    type=str,
    default="127.0.0.1",
    help='Redis host'
)
parser.add_argument(
    '--port',
    type=int,
    default=6379,
    help='Redis port'
)

# Parse the arguments
args = parser.parse_args()


class Taxonomy:

    def __init__(self, store):
        """Read in the taxonomy table."""

        # Read the taxonomy table
        self.taxonomy_df = pd.read_hdf(
            store,
            "/ref/taxonomy"
        ).apply(
            lambda c: c.fillna(0).apply(float).apply(
                int) if c.name in ["parent", "tax_id"] else c,
        ).set_index(
            "tax_id"
        )

        self.all_taxids = set(self.taxonomy_df.index.values)

    @lru_cache(maxsize=None)
    def path_to_root(self, tax_id, max_steps=100):
        """Parse the taxonomy to yield a list with all of the taxa above this one."""

        visited = []

        for _ in range(max_steps):

            # Skip taxa which are missing
            if tax_id not in self.all_taxids:
                break

            # Add to the path we have visited
            visited.append(tax_id)

            # Get the parent of this taxon
            parent_id = self.taxonomy_df.loc[tax_id, "parent"]

            # If the chain has ended, stop
            if parent_id in visited or parent_id == 0:
                break

            # Otherwise, keep walking up
            tax_id = parent_id

        return visited

    @lru_cache(maxsize=None)
    def anc_at_rank(self, tax_id, rank):
        for anc_tax_id in self.path_to_root(tax_id):
            if self.taxonomy_df.loc[anc_tax_id, "rank"] == rank:
                return anc_tax_id

    def name(self, tax_id):
        if tax_id in self.all_taxids:
            return self.taxonomy_df.loc[tax_id, "name"]

    def parent(self, tax_id):
        if tax_id in self.all_taxids:
            return self.taxonomy_df.loc[tax_id, "parent"]

    def make_cag_tax_df(
        self,
        taxa_vc,
        dtype=int
    ):
        """Return a nicely formatted taxonomy table from a list of tax IDs and the number of assignments for each."""

        # We will construct a table with all of the taxa in the tree, containing
        # The ID of that taxon
        # The name of that taxon
        # The number of genes assigned to that taxon (or its children)
        # The rank of that taxon
        # The ID of the parent of that taxon

        # The number of genes found at that taxon or in its decendents
        counts = defaultdict(dtype)

        # Keep track of the total number of genes with a valid tax ID
        total_genes_assigned = 0

        # Iterate over each terminal leaf
        for tax_id, n_genes in taxa_vc.items():

            # Skip taxa which aren't in the taxonomy
            if tax_id not in self.taxonomy_df.index.values:
                continue

            # Count all genes part of this analysis
            total_genes_assigned += n_genes

            # Walk up the tree from the leaf to the root
            for anc_tax_id in self.path_to_root(tax_id):

                # Add to the sum for every node we visit along the way
                counts[anc_tax_id] += n_genes

        if len(counts) == 0:
            return pd.DataFrame([{
                "tax_id": None
            }])

        # Make a DataFrame
        df = pd.DataFrame({
            "count": counts,
        })

        # Add the name, parent, rank
        df = df.assign(
            tax_id=df.index.values,
            parent=self.taxonomy_df["parent"],
            rank=self.taxonomy_df["rank"],
            name=self.taxonomy_df["name"],
            total=total_genes_assigned,
        ).reset_index(
            drop=True
        )

        return df


class collectResults:
    
    def __init__(
        self,
        output_prefix="OUTPUT",
        host="127.0.0.1",
        port=6379,
        fdr_method="fdr_bh",
    ):
        # Save attributes
        self.fdr_method = fdr_method

        # Set up logging
        self.setup_logging()

        # Connect to redis
        self.logger.info(f"Connecting to redis at {host}:{port}")
        self.r = DirectRedis(host=host, port=port)

        # Parse the output filepath from the user-provided arguments
        self.output_hdf_fp = f"{output_prefix}.hdf5"
        self.logger.info(f"Writing out to {self.output_hdf_fp}")
        
        # Read in data objects as needed
        self.load_data()

        # Open a connection to the output file
        self.logger.info("Opening the output store (pandas)")
        output_store = pd.HDFStore(self.output_hdf_fp, "a")

        # Save the gene annotations (name, taxonomic, and functional annotations)
        self.write_gene_annotations()

        # Write out details of genome alignment
        self.write_genome_details(output_store)

        # Write out the genome manifest to the final HDF5
        self.write_genome_manifest(output_store)

        # Write out the combined containment table
        self.write_containment(output_store)

        # Save the corncob results table
        self.write_corncob_results(output_store)

        # Process the abundances of each genome
        self.write_genome_abundances(output_store)

        # Close the output file (with the pandas connector)
        self.logger.info("Closing the output store (pandas)")
        output_store.close()

        # Copy the genome annotations directly into the output HDF5
        # If there are any
        self.write_genome_annotations()

        self.logger.info("Done")

    def setup_logging(self):

        ##################
        # SET UP LOGGING #
        ##################

        # Set the level of the logger to INFO
        logFormatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [buildRedis] %(message)s'
        )
        self.logger = logging.getLogger('buildRedis')
        self.logger.setLevel(logging.INFO)

        # Write to STDOUT
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)

    def load_data(self):

        ####################
        # CONTAINMENT DATA #
        ####################

        # Read in and combine all of the containment tables
        self.containment_df = pd.concat([
            pd.read_csv(
                fp,
                sep=",",
                compression="gzip"
            )
            for fp in os.listdir(".")
            if fp.startswith("containment_shard.") and fp.endswith(".csv.gz")
        ])
        self.logger.info("Read in containment")
        n_genomes = self.containment_df["genome"].unique().shape[0]
        self.logger.info(f"Genomes:\t{n_genomes:,}")
        n_cags = self.containment_df["CAG"].unique().shape[0]
        self.logger.info(f"CAGs:\t{n_cags:,}")

        ###################
        # GENOME MANIFEST #
        ###################
        
        self.genome_manifest = pd.read_csv("genome.manifest.csv")
        # Remove URI columns
        for k in ['uri', 'gff']:
            if k in self.genome_manifest.columns.values:
                self.genome_manifest = self.genome_manifest.drop(columns=k)

        ####################
        # GENE ANNOTATIONS #
        ####################

        self.logger.info("Reading in the gene annotations")
        self.gene_annotations = pd.read_hdf(
            "geneshot.results.hdf5",
            "/annot/gene/all"
        )

        # Save an index of which gene names have which indexes
        self.gene_index = {
            gene_name: gene_ix
            for gene_ix, gene_name in self.gene_annotations['gene'].to_dict().items()
        }

        ############
        # TAXONOMY #
        ############

        self.logger.info("Reading in the taxonomy")
        with pd.HDFStore("geneshot.results.hdf5", "r") as store:
            self.tax = Taxonomy(store)

    def write_gene_annotations(self):
        """Save the annotations on individual genes."""

        for col_name, redis_key in [
            ('gene', 'name'),
            ('CAG', 'cag_id'),
            ('tax_id', 'tax_id'),
            ('eggNOG_desc_ix', 'func_id'),
        ]:
            if col_name in self.gene_annotations.columns.values:

                self.logger.info(f"Saving gene_{redis_key}")
                self.r.set(
                    f"gene_{redis_key}",
                    self.gene_annotations[
                        col_name
                    ].dropna(
                    ).apply(
                        str if col_name == 'gene' else int
                    ).to_dict()
                )


    def write_genome_manifest(self, output_store):

        self.logger.info("Writing out the manifest to HDF")

        # Only keep those genomes which have genes aligned
        # Also set the order to match `self.genomes_to_keep`
        self.genome_manifest = self.genome_manifest.set_index(
            "id"
        ).reindex(
            index=self.genomes_to_keep
        ).dropna(
        ).reset_index()

        # Make sure that the number of genomes match
        msg = (self.genome_manifest.shape[0], len(self.genomes_to_keep))
        assert self.genome_manifest.shape[0] == len(self.genomes_to_keep), msg

        # Save to HDF
        self.genome_manifest.to_hdf(
            output_store,
            "/genomes/manifest"
        )

        # Save the mapping of genome names and IDs to redis
        self.logger.info("Saving `genome_name` to redis")
        self.r.set("genome_name", self.genome_manifest["name"].to_dict())
        self.logger.info("Saving `genome_acc` to redis")
        self.r.set("genome_acc", self.genome_manifest["id"].to_dict())

    # Write out the combined containment table
    def write_containment(self, output_store):
        
        self.logger.info("Writing out the containment to HDF")
        
        self.containment_df.to_hdf(
            output_store,
            "/genomes/cags/containment",
        )

        # Populate a dict with the number of aligned genes per genome
        n_genes_per_genome = dict()

        # Make a dict to transform genome accessions to integer indexes
        genome_ix = {
            v: k for k, v in self.r.get("genome_acc").items()
        }

        # Add a column to the containment DF with the genome index
        self.containment_df = self.containment_df.assign(
            genome_ix=self.containment_df['genome'].apply(
                genome_ix.get
            )
        )
        # Make sure that all values were found
        assert self.containment_df['genome_ix'].isnull().sum() == 0

        # Write out the number of genes assigned to CAGs by genomes
        for cag_id, cag_df in self.containment_df.groupby("CAG"):

            # Write to redis
            self.r.set(
                f"cag_genome_assignments {cag_id}",
                cag_df.set_index('genome_ix')["n_genes"]
            )

        # Write out the number of genes assigned to genomes by CAGs
        for genome_id, genome_df in self.containment_df.groupby("genome_ix"):

            # Write to redis
            self.r.set(
                f"genome_cag_assignments {genome_id}",
                genome_df.set_index('CAG')["n_genes"]
            )

            # Save the number of genes per genome
            n_genes_per_genome[genome_id] = genome_df["n_genes"].sum()

        # Write out the number of genes per genome to redis
        self.logger.info("Writing `n_genes_per_genome` to redis")
        self.r.set(
            "n_genes_per_genome",
            n_genes_per_genome
        )

    # Save the corncob results table
    def write_corncob_results(self, output_store):

        corncob_df = self.read_corncob_results("corncob.results.csv")

        self.logger.info("Writing to /stats/genome/corncob")
        corncob_df.to_hdf(
            output_store,
            "/stats/genome/corncob"
        )

        for parameter, parameter_df in corncob_df.groupby("parameter"):
            # Skip the intercept
            if parameter == "(Intercept)":
                continue

            self.logger.info(f"Writing to genome_association {parameter}")
            self.r.set(
                f"genome_association {parameter}",
                parameter_df
            )


    def read_corncob_results(self, corncob_csv, group_name="genome"):
        # Read in the corncob results
        corncob_df = pd.read_csv(corncob_csv)

        self.logger.info(
            "Read in corncob results for %d %s groups" %
            (corncob_df[group_name].unique().shape[0], group_name)
        )

        # Make a wide version of the table with just the mu. values
        corncob_wide = corncob_df.loc[
            corncob_df["parameter"].apply(
                lambda s: s.startswith("mu.")
            )
        ].pivot_table(
            index=[group_name, "parameter"],
            columns="type",
            values="value"
        ).reset_index(
        ).apply(
            lambda v: v.apply(
                lambda s: s.replace("mu.", "")
            ) if v.name == "parameter" else v
        )

        # Adding the q-value is conditional on p-values being present
        if "p_value" in corncob_wide.columns.values:

            # Add the q-value (FDR-BH)
            corncob_wide = corncob_wide.assign(
                q_value=multipletests(
                    corncob_wide.p_value.fillna(1), 
                    0.2, 
                    self.fdr_method
                )[1]
            )
            corncob_wide = corncob_wide.assign(
                neg_log10_qvalue=corncob_wide["q_value"].apply(np.log10) * -1
            )

        # Add the wald metric
        corncob_wide = corncob_wide.assign(
            wald=corncob_wide["estimate"] / corncob_wide["std_error"]
        )

        # Convert the genome accession to the integer genome index,
        # and use that to set the index
        corncob_wide = corncob_wide.assign(
            genome_ix = corncob_wide["genome"].apply(
                lambda acc: self.genome_index_map[acc]
            )
        ).drop(
            columns="genome"
        ).set_index(
            "genome_ix"
        )

        return corncob_wide


    # Write out details of genome alignment
    def write_genome_details(self, output_store):

        # Copy tables underneath this path
        genomes_detail_group = "/genomes/detail/"

        # Make a set of those genomes which have any alignment details saved
        self.genomes_to_keep = list()

        # Populate a dict of dicts, recording the genomes which are aligned
        # to genes which are assigned to a given taxa
        tax_genome_counts = defaultdict(lambda: dict())

        # Record the set of genes aligned for each genome
        self.genome_membership = dict()

        # Write out the detailed tables of genome alignments, if there are any
        for fp in os.listdir("."):

            # Parse any files with matching filenames
            if not fp.startswith("containment_shard."):
                continue
            if not fp.endswith(".hdf5"):
                continue

            # Open up the file
            with pd.HDFStore(fp, "r") as input_store:
                self.logger.info(f"Reading from {fp}")

                # Iterate over every table
                for k in input_store.keys():

                    # If the key matches the pattern
                    if k.startswith(genomes_detail_group):

                        # Get the genome accession
                        genome_acc = k[len(genomes_detail_group):]

                        self.logger.info(
                            f"Writing alignment details for {genome_acc} ({k})"
                        )

                        # Read the table
                        df = pd.read_hdf(input_store, k)

                        # Copy over the table to the output
                        df.to_hdf(
                            output_store,
                            k
                        )

                        # Record that this is one of the genomes with detailed alignments
                        genome_ix = len(self.genomes_to_keep)
                        self.genomes_to_keep.append(genome_acc)

                        # To save this to redis, we'll first need to
                        # make an index of the contig lengths
                        contig_df = df.reindex(
                            columns=['contig', 'contig_len']
                        ).drop_duplicates(
                        ).sort_values(
                            by="contig_len",
                            ascending=False
                        ).reset_index(drop=True)

                        self.r.set(
                            f"genome_contigs {genome_ix}",
                            contig_df
                        )

                        # Index contigs by name
                        contig_index = {
                            contig_name: contig_ix
                            for contig_ix, contig_name in contig_df["contig"].items()
                        }

                        # gene index, contig, position, and percent identity
                        alignment_df = df.assign(
                            gene_ix=df["gene"].apply(
                                lambda gene_name: self.gene_index[gene_name]
                            ),
                            contig_ix=df["contig"].apply(
                                lambda contig_name: contig_index[contig_name]
                            )
                        ).reindex(
                            columns=[
                                "contig_ix",
                                "gene_ix",
                                "pident",
                                "contig_start",
                                "contig_end"
                            ]
                        )

                        # Save the indexed alignments to redis
                        self.r.set(
                            f"genome_alignment {genome_ix}",
                            alignment_df
                        )

                        # Save the set of genes aligned to this genome
                        self.genome_membership[
                            genome_ix
                        ] = set(
                            df['gene'].tolist()
                        )
                        
                        # If there are taxonomic assignments per genome
                        if self.r.get("gene_tax_id") is not None:

                            # Save a table summarizing the taxonomic assignments
                            taxa_vc = df['gene'].apply(
                                lambda gene_name: self.gene_index[gene_name]
                            ).apply(
                                self.r.get("gene_tax_id").get
                            ).dropna(
                            ).value_counts()

                            # Ignore unassigned genes
                            if 0 in taxa_vc.index.values:
                                taxa_vc = taxa_vc.drop(index=0)

                            # If there are assignments
                            if taxa_vc.shape[0] > 0:

                                # Save a table
                                redis_key = f"genome_tax_assignments {genome_ix}"
                                self.logger.info(f"Saving {redis_key}")
                                self.r.set(
                                    redis_key,
                                    self.tax.make_cag_tax_df(taxa_vc)
                                )

                                # For each taxon, add the count of genes aligned to this genome
                                for tax_id, n in taxa_vc.items():
                                    tax_genome_counts[
                                        tax_id
                                    ][
                                        genome_ix
                                    ] = n

        # Save a mapping of genome accessions to integer indexes
        self.genome_index_map = {
            genome_acc: genome_index
            for genome_index, genome_acc in enumerate(self.genomes_to_keep)
        }

        # Save the number of aligned genes to each genome, per taxa
        for tax_id, genome_counts in tax_genome_counts.items():
            redis_key = f"tax_genome_assignments {tax_id}"
            self.logger.info(f"Saving {redis_key}")
            self.r.set(
                redis_key,
                genome_counts
            )

    # Write out the genome annotations, if there are any
    def write_genome_annotations(self, input_hdf="genome.annotations.hdf5"):

        # Open the output HDF5 with h5py
        # This will make it easier to directly copy data into it
        self.logger.info("Opening the output store (h5py)")
        output_store = h5py.File(self.output_hdf_fp, "a")

        self.logger.info(
            "Attempting to read annotations from genome.annotations.hdf5")
        try:

            annotation_store = h5py.File(input_hdf, "r")

        except:

            self.logger.info("No annotations found")
            annotation_store = False

        if annotation_store is False:

            self.logger.info("Unable to open store")

        else:

            if "annotations" in annotation_store.keys():

                n_annotations = len(annotation_store["annotations"].keys())
                self.logger.info(f"Found {n_annotations:,} annotations")
                n_genomes = len(self.genomes_to_keep)
                self.logger.info(f"Filtering down to {n_genomes:,} genomes")
                self.logger.info("Creating /genomes/annotations group")
                output_store.create_group("/genomes/annotations")

                # Iterate over every genome available
                for genome_id in annotation_store["/annotations/"].keys():
                    if genome_id in self.genomes_to_keep:
                        self.logger.info("Copying annotations for %s" % genome_id)

                        annotation_store.copy(
                            "/annotations/%s" % genome_id,
                            output_store["/genomes/annotations"]
                        )
            else:
                self.logger.info("No /annotations/* found in store")

            annotation_store.close()

            self.logger.info("Done copying annotations")

        self.logger.info("Closing the output store (h5py)")
        output_store.close()

    def write_genome_abundances(self, output_store):
        """Calculate and save the abundances of each genome."""

        # Make a table with the gene membership of all genomes
        self.compute_genome_membership()

        # Make a dict of the abundance for each genome in each specimen
        raw_abund = defaultdict(lambda: dict())
        
        # Make a second dict with the NNLS abundance
        nnls_abund = defaultdict(lambda: dict())

        # Read the abundances for each specimen from the details HDF5
        for specimen_name, specimen_abund in self.parse_specimen_abundance():

            # For each genome, which is aligned to a set of genes
            for genome_ix, aligned_genes in self.genome_membership.items():

                # Compute the proportion of gene copies in this specimen
                # which align to this genome
                raw_abund[specimen_name][genome_ix] = specimen_abund.reindex(
                    index=list(aligned_genes)
                ).fillna(
                    0
                ).sum()

            # Run NNLS to infer the unique abundance per genome
            for genome_ix, genome_abund in self.run_nnls(specimen_abund).items():
                nnls_abund[specimen_name][genome_ix] = genome_abund

        # Transform the abundances to a DataFrame
        raw_abund = pd.DataFrame(raw_abund).fillna(0)
        nnls_abund = pd.DataFrame(nnls_abund).fillna(0)

        # Save the abundances per-genome and per-specimen to redis
        self.save_abundance_table(
            raw_abund,
            "genome_abundance_specimen",
            "specimen_abundance_genome"
        )
        self.save_abundance_table(
            nnls_abund,
            "genome_nnls_specimen",
            "specimen_nnls_genome"
        )

        # Save the average abundance across all specimens
        self.logger.info("Saving mean_abundance_genomes")
        self.r.set(
            "mean_abundance_genomes",
            raw_abund.mean(axis=1)
        )
        self.logger.info("Saving mean_nnls_genomes")
        self.r.set(
            "mean_nnls_genomes",
            nnls_abund.mean(axis=1)
        )

        # Now save the abundances per-specimen to HDF
        for df, hdf_prefix in [
            (raw_abund, "/genome/abund/raw"),
            (nnls_abund, "/genome/abund/nnls")
        ]:
            for specimen_name, genome_abund in df.iteritems():

                # Set up the key to write out to
                hdf_key = f"{hdf_prefix}/{specimen_name}"
                self.logger.info(f"Writing out {hdf_key}")

                # Transform the genome indexes back into accessions
                pd.DataFrame(dict(
                    abund=genome_abund.values,
                    acc=[
                        self.genomes_to_keep[ix]
                        for ix in genome_abund.index.values
                    ]
                )).to_hdf(
                    output_store,
                    hdf_key
                )

    def save_abundance_table(self, abund_df, row_key, col_key):
        for row_ix, row_vals in abund_df.iterrows():
            self.r.set(
                f"{row_key} {row_ix}",
                row_vals
            )

        for col_ix, col_vals in abund_df.iteritems():
            self.r.set(
                f"{col_key} {col_ix}",
                col_vals
            )

    def parse_specimen_abundance(self, hdf_prefix="/abund/gene/long/"):
        """Yield each of the specimens and their abundance from the details HDF5."""

        # Open the HDF store containing detailed geneshot outputs
        with pd.HDFStore("geneshot.details.hdf5", "r") as store:

            # Iterate over every key in the store
            for key in store.keys():

                # If the key starts with the prefix for alignment-based abundance
                if key.startswith(hdf_prefix):

                    self.logger.info(f"Reading data from {key}")

                    # Parse the specimen name from the HDF key
                    specimen_name = key[len(hdf_prefix):]

                    # Format the depth of sequencing for every gene in this specimen
                    gene_abund = pd.read_hdf(
                        store,
                        key
                    ).set_index(
                        "id"
                    )[
                        'depth'
                    ]

                    # Normalize to the depth of sequencing
                    gene_abund = gene_abund / gene_abund.sum()

                    # Yield this information
                    yield specimen_name, gene_abund

    def compute_genome_membership(self):
        
        self.logger.info("Creating a genome membership table")
        
        # 1. Make a list with the unique set of aligned genes, and all contigs with alignments
        all_genomes = list(self.genome_membership.keys())
        all_genes = set()
        for gene_set in self.genome_membership.values():
            all_genes = all_genes | gene_set
        all_genes = list(all_genes)
        self.logger.info(f"Number of aligned genes: {len(all_genes):,}")
        self.logger.info(f"Number of aligned genomes: {len(all_genomes):,}")

        # 2. Make a DataFrame with a column for each genome and a row for each gene
        self.gene_membership = pd.DataFrame(
            np.zeros((len(all_genomes), len(all_genes)), dtype=int),
            index=all_genomes,
            columns=all_genes,
            dtype=int
        )

        # 3. Indicate gene membership by updating cells in the table
        for genome_ix, gene_set in self.genome_membership.items():
            for gene_name in list(gene_set):
                self.gene_membership.loc[
                    genome_ix,
                    gene_name
                ] = 1

        # Remove any duplicate genomes, and put the genes in the index
        self.gene_membership = self.gene_membership.drop_duplicates().T
        self.logger.info(f"# of unique genomes:\t{self.gene_membership.shape[1]}")
        
        self.logger.info("Done creating a genome membership table")

    def run_nnls(self, specimen_abund):
        """Using the genome~gene membership table, run NNLS to infer unique abundances."""

        nnls_abund, nnls_resid = nnls(
            self.gene_membership.values, 
            specimen_abund.reindex(
                self.gene_membership.index.values
            ).fillna(
                0
            ).values
        )

        return pd.Series(
            nnls_abund,
            index=self.gene_membership.columns.values
        )

# Executed as a script
if __name__ == "__main__":

    # Entrypoint
    collectResults(
        output_prefix=args.output_prefix,
        host=args.host,
        port=args.port,
    )
