#!/usr/bin/env python3

import pandas as pd
import pickle
from statsmodels.stats.multitest import multipletests

pickle.HIGHEST_PROTOCOL = 4

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
corncob_df = corncob_df.reindex(
    index=corncob_df["estimate"].dropna().index
)
print("Corncob results have %d rows for mu" % (corncob_df.shape[0]))
assert corncob_df.shape[0] > 0

# Remove the intercept values
corncob_df = corncob_df.loc[
    corncob_df["parameter"] != "mu.(Intercept)"
]
print("Corncob results have %d non-intercept rows for mu" % (corncob_df.shape[0]))
assert corncob_df.shape[0] > 0

# Replace "(Intercept)" with "Intercept" in the parameter table
corncob_df.replace(
    to_replace=dict([("parameter", dict([("(Intercept)", "Intercept")]))]),
    inplace=True
)

# Reformat the corncob results as a dict
corncob_dict = dict()

# Only include parameters which match the --filter, if any:
for parameter, parameter_df in corncob_df.groupby("parameter"):
    if "${params.filter}" == "false" or "${params.filter}" in parameter:
        print("Including parameter: %s" % parameter)
        corncob_dict[parameter] = parameter_df.set_index("CAG")

    else:
        print("NOT Including parameter: %s (does not contain %s)" % (parameter, "${params.filter}"))

assert len(corncob_dict) > 0, "Did not find any parameters to include"

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
