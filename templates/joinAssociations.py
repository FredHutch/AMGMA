#!/usr/bin/env python3

import h5py
import pandas as pd
import pickle
import os
import shutil

pickle.HIGHEST_PROTOCOL = 4

# Open a connection to the output file
output_store = pd.HDFStore("joined_associations.hdf5", "w")

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
