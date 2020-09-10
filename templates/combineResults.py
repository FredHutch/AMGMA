#!/usr/bin/env python3

import h5py
import pandas as pd
import pickle
import os
import shutil

pickle.HIGHEST_PROTOCOL = 4

# Read in and combine all of the containment tables
containment_csv_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("containment_shard.") and fp.endswith(".csv.gz")
]
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
association_shard_hdf_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("association_shard") and fp.endswith(".hdf5")
]
for hdf_fp in association_shard_hdf_list:

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
            elif k.startswith("/genomes/map/"):
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

# Close the HDF5 store
print("Closing the output store (pandas)")
output_store.close()

# Open the output HDF5 with h5py
# This will make it easier to directly copy data into it
print("Opening the output store (h5py)")
output_store = h5py.File("${params.output_hdf}", "a")

# Make the /genomes/details/ group
genomes_detail_group = "/genomes/detail"
output_store.create_group(genomes_detail_group)

# Make a set of those genomes which have any alignment details saved
genomes_to_keep = set([])

# Write out the detailed tables of genome alignments, if there are any
for fp in os.listdir("."):

    # Parse any files with matching filenames
    if fp.startswith("genome_alignments.") and fp.endswith(".hdf5"):

        # Open up the file
        with h5py.File(fp, "r") as input_store:

            # Iterate over every table
            for k in input_store[genomes_detail_group].keys():

                # Format the group for the table
                path = "{}/{}".format(genomes_detail_group, k)

                # Copy over the table
                print("Writing alignment details for %s" % k)
                input_store.copy(
                    path,
                    output_store[genomes_detail_group]
                )

                # Record that this is one of the genomes with detailed alignments
                genomes_to_keep.add(k)

# Write out the genome annotations, if there are any
print("Attempting to read annotations from genome_annotations.hdf5")
try:

    annotation_store = h5py.File("genome_annotations.hdf5", "r")

except:

    print("No annotations found")
    annotation_store = False

if annotation_store is False:

    print("Unable to open store")

else:

    if "annotations" in annotation_store.keys():

        print("Found %d annotations" % len(annotation_store["annotations"].keys()))
        print("Filtering down to %d genomes which also have alignments" % len(genomes_to_keep))
        print("Copying to ${params.output_hdf}")
        print("Creating /genomes/annotations group")
        output_store.create_group("/genomes/annotations")

        # Iterate over every genome available
        for genome_id in annotation_store["/annotations/"].keys():
            if genome_id in genomes_to_keep:
                print("Copying annotations for %s" % genome_id)

                annotation_store.copy(
                    "/annotations/%s" % genome_id,
                    output_store["/genomes/annotations"]
                )
    else:
        print("No /annotations/* found in store")

    annotation_store.close()

    print("Done copying annotations")

print("Closing the output store (h5py)")
output_store.close()
print("Done")
