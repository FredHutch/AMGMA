# Annotation of Microbial Genomes by Microbiome Association (AMGMA)

The purpose of this repository is to allow users to annotate a collection of
microbial genomes on the basis of their association with microbiome survey
datasets.

### Background

Diving more deeply, every microbial genome can be viewed as consisting of a
collection of genetic components, such as genes. In turn, each of those genetic
components (e.g. genes) can be associated with some outcome of interest via
the analysis of microbiome surveys by metagenomic sequencing. The tool most easily
used for this purpose (and compatible for downstream analysis by AMGMA) is
[GeneShot](https://www.github.com/Golob-Minot/GeneShot).

After each microbial gene has been associated with outcome of interest (such
as a particular human disease), then those association metrics can be summarized
over a set of microbial genomes by:

  1. Aligning the reference genes against each genome, and
  2. Calculating summary metrics for each genome, such as:
    * Proportion of the genome passing FDR filter for association
    * Average estimated coefficient of association for those genes
  3. Formatting a genome 'map' showing the relative physical location for those genes throughout the genome

### Running AMGMA

_Formatting a Reference Genome Database_:

Using a collection of reference genomes, run the `AMGMA/build_db` script.

Input data for this step is simply a collection of microbial genomes, and
a manifest CSV file listing those genomes. The CSV must have the following
columns: `uri,id,name`:

* `uri`: Location of the genome FASTA
* `id`: Unique alphanumeric ID for each genome
* `name`: Longer description of each genome, with whitespaces allowed (but no commas)

```
Usage:

nextflow run FredHutch/AMGMA/build_db.nf <ARGUMENTS>

Required Arguments:
  --manifest            CSV file listing samples (see below)
  --output_folder       Folder to write output files to
  --output_prefix       Prefix to use for output file names

Optional Arguments:
  --batchsize           Number of genomes to process in a given batch (default: 100)

Manifest:
  The manifest is a CSV listing all of the genomes to be used for the database.
  The manifest much contain the column headers: uri,id,name
  The URI may start with ftp://, s3://, or even just be a path to a local file.
  The ID must be unique, and only contain a-z, A-Z, 0-9, or _.
  The NAME may be longer and contain whitespaces, but may not contain a comma.
```

_Annotate Reference Genomes_:

Using that genome database, annotate with a set of microbiome association results
generated by `geneshot` with the main `AMGMA` executable.

```
Usage:

nextflow run FredHutch/AMGMA <ARGUMENTS>

Required Arguments:
--db                  AMGMA database (ends with .tar)
--geneshot_hdf        Results HDF file output by GeneShot, containing CAG information
--geneshot_dmnd       DIAMOND database for the gene catalog generated by GeneShot
--output_folder       Folder to write output HDF file into
--output_hdf          Name of the output HDF file to write to the output folder

Optional Arguments:
--details             Include additional detailed results in output (see below)
--min_coverage        Minimum coverage required for alignment (default: 80)
--min_identity        Minimum percent identity required for alignment (default: 80)
--fdr_method          Method used for FDR correction (default: fdr_bh)
--alpha               Alpha value used for FDR correction (default: 0.2)

Output HDF:
The output from this pipeline is an HDF file which contains all of the data from the
input HDF, as well as the additional tables,

* /genomes/manifest
* /genomes/cags/containment
* /genomes/summary/<feature>
* /genomes/detail/<feature>/<genome_id> (Included with --details)

for each <feature> tested in the input, and for each <genome_id> in the database
```

### Example Outputs

_/genomes/cags/containment_:

| CAG | cag_prop | containment | genome | genome_bases | genome_prop | n_genes |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.9701 | 0.9701 | GCF_000840245 | 40890 | 0.8430 | 65 |
| 0 | 0.0895 | 0.0895 | NC_000913.3 | 2495 | 0.0005 | 6 |
| 1 | 0.8888 | 0.9630 | GCF_000819615 | 5187 | 0.9630 | 8 |
| 0 | 0.9701 | 0.9701 | GCF_000840245_FTP | 40890 | 0.8430 | 65 |

_/genomes/summary/label1_:

| genome_id | mean_est_coef | n_pass_fdr | parameter | prop_pass_fdr | total_genes |
| --- | --- | --- | --- | --- | --- |
| GCF_000840245 | 0.1912 | 38 | label1 | 0.65 | 66 |
| NC_000913.3 | -0.0056 | 1 | label1 | 0.12 | 6 |
| GCF_000819615 | 0.0581 | 7 | label1 | 0.8 | 8 |
| GCF_000840245_FTP | 0.0272 | 58 | label1 | 0.84 | 66 |

### How AMGMA Works

The concept of AMGMA is that it will compare a set of genes from an external
gene catalog against a collection of genomes. Each of those 'catalog genes'
has been associated with some measure(s) of interest using the intermediate
entity of Co-Abundant Gene Groups (CAGs). Each CAG has been analyzed for its
estimated coefficient of association with some measure of interest, and we
also have a p-value for that association. At the level of the CAGs, we can
apply some False Discovery Rate (FDR) correction, and then propagate those
corrected p-values and estimated coefficients to each of the catalog genes
found within each CAG.

The next step is to figure out which of the 'catalog genes' is present in each
of the reference genomes, along with the position of those genes in the genome.
We do this by directly aligning each reference genome against the gene catalog,
using conceptual six-frame translation with DIAMOND in 'blastx' mode and only
retaining alignments which are the highest scoring for each predicted coding
region of the genome. After that alignment step, we then summarize each genome
on the basis of the estimated coefficient and FDR-corrected p-value for each
catalog gene.

### Database Format

The reference database made with the `build_db` script will create a single
tarball containing a manifest CSV as well as a set of collection of further
nested tarballs (one for each set of `--batchsize` genomes) containing the
genome sequences for that batch of genomes. One additional file is contained
in those tarballs, which is a table identifying which nucleotide sequence
header corresponds to which reference genome in the concatenated FASTA file.

### Output Format

After running the main `AMGMA` script, a single output file will be created
which adds the AMGMA results to the input HDF. This output will include 
everything present in the HDF file used as an input (matching the output from
`geneshot`). 

The idea is that `geneshot` was run in such a way as to include some statistical
analysis testing for the association of every individual CAG with some features
of interest. And so there may be estimated coefficients of association for
multiple `<feature>` elements in the HDF file used as input for `AMGMA`, whose
output will then include the additional tables:

* `/genomes/manifest`: A table with the name and path to each genome
* `/genomes/summary/<feature>`: A table with summary metrics for each genome (for a given `feature`)
* `/genomes/detail/<feature>/<id>`: For each genome, a map of gene locations
