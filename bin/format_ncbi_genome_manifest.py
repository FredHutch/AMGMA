#!/usr/bin/env python3

import argparse
import logging
import os
import pandas as pd


def format_genbank_record(r, mask_characters=[",", ";", "/", "\\"]):
    """Function to reformat the information from NCBI in AMGMA manifest format."""

    # Format the full path to the genome from the FTP prefix
    uri = "{}/{}_genomic.fna.gz".format(
        r["RefSeq FTP"],
        r["RefSeq FTP"].rsplit("/", 1)[1]
    )

    # Also format the path to the genome annotations in GFF format
    gff = "{}/{}_genomic.gff.gz".format(
        r["RefSeq FTP"],
        r["RefSeq FTP"].rsplit("/", 1)[1]
    )

    # Use the assembly ID as the id
    uuid = r["Assembly"]

    # Format the name as a combination of the organism name and the assembly ID
    name = "{} ({})".format(r["#Organism Name"], uuid)

    # Do some sanitizing of the name
    for ch in mask_characters:
        name = name.replace(ch, "_")

    # Collapse adjacent underscores
    while "__" in name:
        name = name.replace("__", "_")

    return pd.Series({
        "name": name,
        "id": uuid,
        "uri": uri,
        "gff": gff,
    })

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
Script to parse the CSV tables listing genomes which are downloaded
from NCBI (https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes).

Using a set of CSV input file(s), this script will create a single
manifest CSV which can be parsed by AMGMA to build a genome database.

Usage:

format_ncbi_genome_manifest.py \
    --input_csv <INPUT_CSV1> <INPUT_CSV2> ... \
    --output_csv <OUTPUT_CSV>

""")
    parser.add_argument(
        "--input-csv", 
        help = "Input CSV file(s) to process",
        required = True, 
        nargs='+',
        type=str,        
    )
    parser.add_argument(
        "--output-csv", 
        help = "Output AMGMA genome manifest CSV",
        required = True, 
        type=str,        
    )
    args = parser.parse_args()

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [AMGMA Parse Genome Manifest] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Make sure the output file doesn't already exist
    fpo = args.output_csv
    assert os.path.exists(fpo) is False, "{} already exists".format(fpo)

    # Make sure that the input files all do exist
    for fp in args.input_csv:
        assert os.path.exists(fp), "{} does not exist".format(fp)

        logging.info(
            "Found {:,} records in {}".format(
                pd.read_csv(fp).shape[0],
                fp
            )
        )

    # Read in the input files
    df = pd.concat(
        [
            pd.read_csv(fp)
            for fp in args.input_csv
        ], 
        sort=True
    ).reset_index(
        drop=True
    )

    # Only keep genomes which have valid GenBank paths
    df = df.reindex(
        index=df["GenBank FTP"].dropna().index.values
    )

    # Format the manifest in the format expected by AMGMA (uri,id,name)
    df = df.apply(
        format_genbank_record,
        axis=1
    )

    # Drop duplicates
    vc = df["id"].value_counts()
    if (vc > 1).sum() > 0:
        logging.info(
            "Dropping {:,} duplicated records".format(
                (vc > 1).sum()
            )
        )
        df = df.drop_duplicates()
        assert df["id"].value_counts().max() == 1

    # Write out to a file
    df.to_csv(fpo, index=None)

    logging.info("Wrote out {:,} genomes to {}".format(
        df.shape[0], fpo
    ))
