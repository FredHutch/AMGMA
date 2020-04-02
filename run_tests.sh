#!/bin/bash

# This script can be used to run some simple tests locally.
# In order for the tests to run, clone this repo locally, install
# Nextflow, and execute this from the cloned repository folder

# The same tests are run automatically for each pushed commit
# using GitHub Actions, so this script is more for local development

set -e

[[ ! -d testing_output ]] && mkdir testing_output

NXF_VER=20.01.0 \
    nextflow \
    run \
    build_db.nf \
    -with-docker ubuntu:18.04 \
    -w work/ \
    --manifest data/genome.manifest.csv \
    --output_folder testing_output/db \
    --output_prefix amgma \
    -profile testing \
    --batchsize 2 \
    -resume
    
NXF_VER=20.01.0 \
    nextflow \
    run \
    align.nf \
    -with-docker ubuntu:18.04 \
    -w work/ \
    --geneshot_hdf data/geneshot.summary.hdf5 \
    --geneshot_dmnd data/geneshot.dmnd \
    --db testing_output/db/amgma.tar \
    --output_folder testing_output/align \
    --output_hdf amgma_test.output.hdf \
    -profile testing \
    -resume
    
NXF_VER=20.01.0 \
    nextflow \
    run \
    align.nf \
    -with-docker ubuntu:18.04 \
    -w work/ \
    --geneshot_hdf data/geneshot.summary.hdf5 \
    --geneshot_dmnd data/geneshot.dmnd \
    --db testing_output/db/amgma.tar \
    --output_folder testing_output/align \
    --output_hdf amgma_test.output.hdf \
    --details \
    -profile testing \
    -resume