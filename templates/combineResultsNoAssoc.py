#!/usr/bin/env python3

import h5py
import pandas as pd
import pickle
import os
import shutil

pickle.HIGHEST_PROTOCOL = 4

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

# Make a set of those genomes for which any CAG has >= 10% of genes aligned
genomes_to_keep = set(
    containment_df.query(
        "cag_prop >= 0.1"
    )[
        "genome"
    ].tolist()
)

# Close the HDF5 store
print("Closing the output store (pandas)")
output_store.close()

# Open the output HDF5 with h5py
# This will make it easier to directly copy data into it
print("Opening the output store (h5py)")
output_store = h5py.File("${params.output_hdf}", "a")

# Write out the detailed tables of genome alignments, if there are any
for fp in os.listdir("."):

    # Parse any files with matching filenames
    if fp.startswith("genome_alignments.") and fp.endswith(".hdf5"):

        # Open up the file
        with h5py.File(fp, "r") as input_store:

            # Iterate over every table
            for k in input_store.keys():

                # Copy over the table
                print("Writing alignment details for %s" % k)
                input_store.copy(
                    k,
                    output_store[k]
                )

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

        # Make sure that the /genomes/annotations/ group exists in the output HDF5
        if "annotations" not in output_store["/genomes/"].keys():
            output_store.create_group("/genomes/annotations")

        # Iterate over every genome available
        for genome_id in annotation_store["/annotations/"].keys():
            if genome_id in genomes_to_keep:
                print("Copying annotations for %s" % genome_id)

                annotation_store.copy(
                    "/annotations/%s" % genome_id,
                    output_store["/genomes/annotations/%s" % genome_id]
                )
    else:
        print("No /annotations/* found in store")

    annotation_store.close()

    print("Done copying annotations")

print("Closing the output store (h5py)")
output_store.close()
print("Done")
