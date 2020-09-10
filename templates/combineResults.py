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


    # Append the annotations directly into the results
    output_store = h5py.File("${params.output_hdf}", "a")

    if "annotations" in annotation_store.keys():

        print("Found %d annotations" % len(annotation_store["annotations"].keys()))
        print("Copying to ${params.output_hdf}")

        # Copy the entire group of annotations
        annotation_store.copy(
            "/annotations/",
            output_store["/genomes/"]
        )
    else:
        print("No /annotations/* found in store")

    output_store.close()
    annotation_store.close()

    print("Done copying annotations")