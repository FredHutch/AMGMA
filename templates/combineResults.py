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

# Close the output file (with the pandas connector)
print("Closing the pandas connection to the output store")
output_store.close()

# Open the output HDF5 with h5py
# This will make it easier to directly copy data into it
print("Opening the output store (h5py)")
output_store = h5py.File("${params.output_hdf}", "a")

# Open the HDF5 with all of that association data available
with h5py.File("associations.hdf5", "r") as input_store:

    # Copy all of the data from these two groups
    for group_name in ["/genomes/summary", "/genomes/map"]:

        # Make the group in the output store
        output_store.create_group(group_name)

        print(f"Copying all data from {group_name} to output")
        n = 0

        # Iterate over every table in the group
        for k in input_store[group_name].keys():

            # Format the group for the table
            path = f"{group_name}/{k}"

            # Copy over the table
            input_store.copy(
                path,
                output_store[group_name]
            )
            n += 1
        print(f"Copied {n:,} tables")

print("Done copying association data")

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
