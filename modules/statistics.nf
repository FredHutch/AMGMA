// // Processes to perform statistical analysis

// Extract the counts table from the results HDF5
process runCorncob {
    tag "Perform statistical analysis"
    container "quay.io/fhcrc-microbiome/corncob"
    label "mem_medium"
    
    input:
    file readcounts_csv_gz
    file metadata_csv
    val group_name
    each formula

    output:
    file "corncob.results.csv"


    """
#!/usr/bin/env Rscript

# Get the arguments passed in by the user

library(tidyverse)
library(corncob)
library(parallel)

## By default, use 10% of the available memory to read in data
connectionSize = 100000 * ${task.memory.toMega()}
print("Using VROOM_CONNECTION_SIZE =")
print(connectionSize)
Sys.setenv("VROOM_CONNECTION_SIZE" = format(connectionSize, scientific=F))

numCores = ${task.cpus}

##  READCOUNTS CSV should have columns `specimen` (first col) and `total` (last column).
##  METADATA CSV should have columns `specimen` (which matches up with `specimen` from
##         the recounts file), and additional columns with covariates matching `formula`

##  corncob analysis (coefficients and p-values) are written to OUTPUT CSV on completion

print("Reading in ${metadata_csv}")
metadata <- vroom::vroom("${metadata_csv}", delim=",")

print("Removing columns which are not in the formula")
for(column_name in names(metadata)){
    if(column_name == "specimen" || grepl(column_name, "${formula}", fixed=TRUE) ){
        print(paste("Keeping column", column_name))
    } else {
        print(paste("Removing column", column_name))
        metadata <- metadata %>% select(-column_name)
    }
}
metadata <- metadata %>% unique %>% drop_na
print("Filtered and deduplicated manifest:")
print(metadata)

print("Reading in ${readcounts_csv_gz}")
counts <- vroom::vroom("${readcounts_csv_gz}", delim=",")
total_counts <- counts[,c("specimen", "total")]

print("Adding total counts to manifest")
print(head(total_counts))

print("Merging total counts with metadata")
total_and_meta <- metadata %>% 
  left_join(total_counts, by = c("specimen" = "specimen"))

#### Run the analysis for every individual ${group_name} (in this shard)
print(sprintf("Starting to process %s columns (${group_name})", length(c(2:(dim(counts)[2] - 1)))))
corn_tib <- do.call(rbind, mclapply(
    c(2:(dim(counts)[2] - 1)),
    function(i){
        try_bbdml <- try(
            counts[,c(1, i)] %>%
            rename(W = 2) %>%
            right_join(
                total_and_meta, 
                by = c("specimen" = "specimen")
            ) %>%
            corncob::bbdml(
                formula = cbind(W, total - W) ~ ${formula},
                phi.formula = ~ 1,
                data = .
            )
        )

      if (class(try_bbdml) == "bbdml") {
        return(
            summary(
                try_bbdml
            )\$coef %>%
            as_tibble %>%
            mutate("parameter" = summary(try_bbdml)\$coef %>% row.names) %>%
            rename(
                "estimate" = Estimate,
                "std_error" = `Std. Error`,
                "p_value" = `Pr(>|t|)`
            ) %>%
            select(-`t value`) %>%
            gather(key = type, ...=estimate:p_value) %>%
            mutate("${group_name}" = names(counts)[i])
        )
      } else {
          return(
              tibble(
                  "parameter" = "all",
                  "type" = "failed", 
                  "value" = NA, 
                  "${group_name}" = names(counts)[i]
              )
          )
      }   
    },
    mc.cores = numCores
  ))

print(head(corn_tib))

print("Adding a column with the formula used here")
corn_tib <- corn_tib %>% add_column(formula = "${formula}")

print(head(corn_tib))

print(sprintf("Writing out %s rows to corncob.results.csv", nrow(corn_tib)))
write_csv(corn_tib, "corncob.results.csv")
print("Done")
    """

}

// Join together a set of corncob results CSVs
process joinCorncob {
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label "io_limited"
    publishDir "${params.output_folder}/stats/", mode: "copy"
    
    input:
    file "corncob.results.*.csv"
    val group_name

    output:
    file "corncob.results.csv"


"""
#!/usr/bin/env python3

import os
import pandas as pd

# Get the list of files to join
fp_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("corncob.results.") and fp.endswith(".csv")
]

print("Reading in corncob results for %d formula(s)" % len(fp_list))

df = pd.concat([
    pd.read_csv(fp)
    for fp in fp_list
]).query(  # Filter out the unassigned groups
    "${group_name} != 'UNASSIGNED'"
)

print("Writing out to corncob.results.csv")
df.to_csv("corncob.results.csv", index=None)
print("Done")
"""

}
