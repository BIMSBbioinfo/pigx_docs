# PiGx SARS-CoV-2

# Introduction

PiGx SARS-CoV-2 is a pipeline for analysing data from sequenced wastewater
samples and identifying a provided set of lineages of SARS-CoV-2 characterized
by signature mutations. Currently wastewater samples are used, which are
enriched for SARS-CoV-2. The pipeline can be used for continuous sampling. The
output of the PiGx SARS-CoV-2 pipeline is summarized in a report which provides
an intuitive visual overview about the development of lineage abundance and
single significantly increasing mutations over time and location. Additionally
there will be more detailed reports per sample, which cover the quality control
of the samples, the detected variants and a taxonomic classification of all
unaligned reads. This version of the pipeline was designed to work with
paired-end amplicon sequencing data e.g. following the ARtIc protocols
[ARTIC nCoV-2019 primers](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3).

## Workflow

In the first step the pipeline takes the raw reads and the additional
information about used primers and adapters to perform extensive quality
control. Primer trimming is done with
[iVAR](https://github.com/andersen-lab/ivar), and
[fastp](https://github.com/OpenGene/fastp) is used for adapter trimming and
filtering. Next, the trimmed reads are aligned to the reference genome of
SARS-CoV-2 using [BWA](https://github.com/lh3/bwa), and the results are
*SAM*/*BAM* files of **aligned** and **unaligned reads**. Following the
alignment a quality check on raw and processed reads is performed by using
[MultiQC](https://multiqc.info/). Furthermore samples are checked for genome
coverage and how many of the provided signature mutation sites are covered.
Based on this every samples gets a quality score. Samples with genome coverage
below 90% are reported as discarded samples, as they are not included in time
series analysises and summaries.

Calling the variants and inferring SNVs (single nucleotide polymorphisms) on the
**aligned reads** is done with [LoFreq](https://csb5.github.io/lofreq/).
Mutations are annotated with [VEP](https://covid-19.ensembl.org/index.html).
Estimation of lineage frequencies is done by deconvolution (see Methods in the
realted publication for details). To investigate the abundance of RNA matching
other existing species in the wastewater samples the **unaligned reads** will be
taxonomicly classified with [Kraken2](https://github.com/DerrickWood/kraken2).
Kraken2 requires a database downloaded locally of the genomes against the reads
are getting aligned. For documentation how to set this up, see:
[Prepare databases](#prepare-databases). For a better and interactive
visualization of all species present in the wastewater
[Krona](https://github.com/marbl/Krona/wiki) is used. Also here a small step of
setting up a database is needed before running the pipeline, see:
[Prepare databases](#prepare-databases). Interactive reports are generated using
[R-markdown](https://rmarkdown.rstudio.com/) and
[plotly for R](https://plotly.com/r/) for visualizations.

### Pooling of samples for time series analysis and plots

For summarizing across daytime and location, the lineage frequencies are pooled
by calculating the weighted average using the total number of reads of each
sample as weights. The mutation frequencies are pooled by using the simple mean
setting removal of missing values to FALSE.

## Output

* Overview report including:
  * Visualization of the development of SARS-CoV-2 variants and mutations over
     time and locations from all samples provided.
  * Quality Control report per sample: number of covered amplicons, read
     coverage and MultiQC report for raw and trimmed reads
  * Variants report per sample: variant analysis of SARS-CoV-2 from each
     wastewater sample and identification of variants of concern
  * Taxonomic classification: a table and pie chart of the species found in the
     unaligned reads
* *SAM* / *BAM* files per sample: aligned and unaligned reads against SARS-CoV-2
* *VCF* / *CSV* files per sample: listing all detected single nucleotide
  variants (SNVs) from the aligned reads
* VEP reports per sample as *TXT* / *HTML*: report files show the
  [VEP output](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#output)
  including the uploading variants, consequenting amino acid changes and
  consequences for the corresponding protein
* Kraken2 files per sample *(txt)*: provides an overview of all found species in
  the unaligned reads together with the NCBI taxonomy ID
* Log files for all major analysis steps performed by the pipeline

# Installation

Pre-built binaries for PiGx are available through
[GNU Guix](https://gnu.org/s/guix), the functional package manager for
reproducible, user-controlled software management. You can install the PiGx
SARS-CoV-2 pipeline with

```sh
guix install pigx-sars-cov-2
```

If you want to install PiGx SARS-CoV-2 from source, please clone this repository
and change directory accordingly:

```sh
git clone https://github.com/BIMSBbioinfo/pigx_sars-cov-2.git
cd pigx_sars-cov-2
```

To fetch code that is common to all PiGx pipelines run this:

```sh
git submodule update --init
```

Before setting everything up, though, make sure
[all dependencies](https://github.com/BIMSBbioinfo/pigx_sars-cov-2/blob/main/manifest.scm)
are met by either installing them manually, or by entering the provided
reproducible Guix environment. If you are using Guix we definitely recommend the
latter. This command spawns a sub-shell in which all dependencies are available
at exactly the same versions that we used to develop the pipeline:

```sh
USE_GUIX_INFERIOR=t guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

To use your current Guix channels instead of the fixed set of channels, just
omit the `USE_GUIX_INFERIOR` shell variable:

```sh
guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

Note that `--pure` unsets all environment variables that are not explicitly
preserved. To access other executables that are not part of the environment
please address them by their absolute file name.

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

At this point you are able to run PiGx SARS-CoV-2. To see all available options
type `--help`.

```sh
pigx-sars-cov-2 --help
```

## Prepare databases

Before the pipeline can work, three databases must be downloaded and their
location will need to be provided in the settings file. Depending on the size of
the databases this can take some time.

We provide the script `download_databases.sh` to automate the download and
configuration of all default databases. After installing the pipeline you can
run the script like this:

```sh
prefix="$(dirname pigx-sars-cov-2)/../"
$prefix/libexec/pigx_sars-cov-2/scripts/download_databases.sh
```

Read on for details if you want to download the databases manually.

### Kraken2 database

There are several libraries of genomes that can be used to classify the
(unaligned) reads. It is up to you which one to use, but be sure that they
fulfill the necessities stated by Kraken2
[Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases).
For an overall overview we recommend to use the Plus-PFP library provided
[here](https://benlangmead.github.io/aws-indexes/k2). If the classification is
not of concern or only the viruses are of interest, we recommend using a smaller
one. This will speed up the pipeline. It is also possible to have multiple
Kraken2 databases installed, just be sure to provide the correct location to the
settings file.

After downloading and unpacking the database files, use `kraken2-build` to
download the taxonomy data and build the database.

### Krona database

Krona Tools needs two files, which are fetched with the `updateTaxonomy.sh` and
`updateAccessions.sh` scripts that come with Krona.

### VEP database

Just download the `SARS_CoV_2` database for VEP (variant effect predictor) and
unpack it in the `databases/vep_db/` directory.

# Quick Start

To check whether the pipeline and the databases have been properly set up, run
the pipeline on a minimal test dataset.

1. Download the test data

    `git clone https://github.com/BIMSBbioinfo/pigx_sars-cov-2 sarscov2-test`

2. Enter the directory

    `cd sarscov2-test`

3. Run the pipeline

    `pigx-sars-cov-2 -s .tests/settings.yaml .tests/sample_sheet.csv`

Inside `tests/` a new directory `output` is created, which includes specific
directories containing output data for the respective step of the pipeline. The
`tests/output/reports/index.html` gives the overview over all merged reports for
the test data.

# Preparing the input

In order to run the pipeline, you need to supply

* a sample sheet (CSV format): containing information about sampling date and  
  location
* a settings file (YAML format) for specifying the experimental setup and
  optional custom parameter adjustments  
* a mutation sheet containing the lineages of interest and their signature
  mutations in nucleotide notation (CSV format)  
* a BED file containing the genomic coordinates of the mutation sites (see below
  for details)  
* the reference genome of the target species in fasta format (so far the
  pipeline is only optimized for SARS-CoV-2, others might work too but not yet  
  tested)
* a BED file containing the PCR primer locations (e.g the primers suggested from
  ARTIC protocols)

In order to generate template settings and sample sheet files, type

```sh
pigx-sars-cov-2 --init
```

in the shell, and a boilerplate `sample_sheet.csv` and `settings.yaml` will be
written to your current directory. An example for both files is provided in the
`tests/` directory.

## Sample sheet

The sample sheet is a tabular file (`csv` format) describing the experiment. The
table has the following columns:

| SampleName | Read           | Read2          | date                | location_name | coordinates_lat | coordinates_long |
| ------     | -------        | --------       | --------            | --------      | --------        | --------         |
| Test0      | Test0_R1.fastq | Test0_R2.fastq | 2021-01-01T08:00:00 | Berlin        | 52.364          | 13.509           |
| Test2      | Test2_R1.fastq | Test2_R2.fastq | 2021-01-03T08:00:00 | Munich        | 48.208          | 11.628           |

* *SampleName* is the name for the sample
* *Read* & *Read2* are the fastq file names of paired end reads
  * the location of these files is specified in the settings file
  * single-end data is not yet supported
  * compressed (`.gz`) reads are not yet supported
* *date* is a date/time in ISO format (`yyyy-mm-ddThh:mm:ss`)
* *location_name* is the name of the location and should be unique per
  coordinates
* *coordinates_lat* & *coordinates_long* correspond the latitude and longitude
  of the location name

## Settings file

The settings file contains parameters (in YAML format) to configure the
execution of the PiGx SARS-CoV-2 pipeline. It specifies:

**Locations**:

* *output-dir*, the location of the outputs for the pipeline
* *reads-dir*, the location of the reads (directory where `fastq` files are)
* *reference-fasta*, the `fasta` file with the reference genome (must be
  prepared by the user)
* *amplicons-bed*, the amplicons `bed` file for coronavirus (must be prepared by
  the user)
* *kraken-db-dir*, the location of the kraken database (must be prepared by the
  user)
* *krona-db-dir*, the location of the krona database (must be prepared by the
  user)
* *vep-db-dir*, the location of `sars_cov_2` database for VEP (must be prepared
  by the user)

## Mutation sheet

The mutation sheet should contain one column of siganture mutations in
nucleotide-NT (not protein/aminoacid mutation) format per lineage that shall be
tracked and analysed by deconvolution. There is no upper or lower limit for the
number of signature mutations per lineage. However, please note that the
deconvolution results are more robust and precise with a higher number of
mutations. (Tested with 10-30 mutations per lineage).

## Mutation BED file

The BED file for testing if the mutation sites are covered should have 4  
columns:

1. chromosome name  
2. start - 5bp **before** the mutation location  
3. end - 5bp **after** the mutation location  
4. name with the original location of the mutation in the format of:
   name_MutationLocation_name, e.g: "nCoV-2019_210_SigmutLocation" A row in the
   BED file should look like:  
`NC_045512.2 205 215 nCoV-2019_210_SigmutLocation`  
Please see the example file within the test directory for a detailed example.

## Primer BED file  

The primer file is needed for trimming with fastp. It should have 6 columns
following BED files standards. An example file can be found in the test
directory.

# Running the pipeline

PiGx SARS-CoV-2 wastewater is executed using the command
`pigx-sars-cov-2 -s settings.yaml sample_sheet.csv`. See
`pigx-sars-cov-2 --help` for information about additional command line
arguments.

The `execution` section of the settings file provides some control over the
execution of the pipeline.

## Local / cluster execution

The workflow may be executed locally (on a single computer), or, if a Sun Grid
Engine-compatible HPC environment is available, supports cluster execution. In
order to enable cluster execution, specify `submit-to-cluster: yes` in the
settings file.

## Parallel execution

If the workflow is run on a cluster, or a single computer with sufficient
resources, some of the tasks may be computed in parallel. To specify the allowed
level or parallelism, set the `jobs` setting under `execution` in the settings
file. For instance,

```yaml
execution:
  submit-to-cluster: yes
  jobs: 40
```

in the settings file, will submit up to 40 simultaneous compute jobs on the
cluster.

# Output description

PiGx SARS-CoV-2 wastewater creates an output directory, as specified in the
settings file, that contains all of the following outputs.

## Time series reports

This pipeline performs mutation analysis of SARS-CoV-2 and reports and
quantifies the occurrence of *variants of concern* (VOC) and signature mutations
by which they are characterised.

The visualizations in the time series report provide an overview of the
evolution of VOCs and signature mutations found in the analyzed samples across
given time points and locations. The abundance values for the variants are
derived by deconvolution. The frequencies of the mutations are the output of
[LoFreq](https://csb5.github.io/lofreq/).

## Quality control

A quality control report is generated for each sample.  It includes
reports on amplicon coverage and read coverage, as well as general
quality control and preprocessing metrics.

General quality control metrics are computed using
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
and [MultiQC](https://multiqc.info/).  The MultiQC report is
particularly useful, collating quality control metrics from many steps
of the pipeline in a single HTML report, which may be found under the
`multiqc` directory in the PiGx output folder.

## Taxonomic classification

This report provides an overview of the species found in the provided
wastewater samples apart from SARS-CoV-2. The SARS-CoV-2 enriched
wastewater samples are aligned to the virus genome. It provides
insight about possible biases and contamination of the samples. In
case of abundance of species very similar to SARS-CoV-2 only the
taxonomic family will be reported which could indicate that an
identification/alignment of SARS-CoV-2 could be biased or
impossible. In case of a high percentage of read matching SARS-CoV-2 a
refining of trimming parameters should be considered.

## Variant report

This report shows the variant analysis of SARS-CoV-2 from wastewater
samples. Mutations are identified by single-nucleotide-variant (SNV)
calling performed by [LoFreq](https://csb5.github.io/lofreq/).
Translated to amino acid mutations by using [Ensemble VEP -
COVID-19](https://covid-19.ensembl.org/info/docs/tools/vep/index.html).
The list of found mutations (including synonymous and non-synonymous
mutations) were matched against lists of signature mutations
characterising variants of concern (VOC) of SARS-CoV-2 provided by
[outbreak.info](https://outbreak.info/situation-reports) and
[CoVariant.org](https://covariants.org/variants/S.501Y.V1).

# Troubleshooting

If you have any questions please e-mail: pigx@googlegroups.com or use the web
form to ask questions <https://groups.google.com/forum/#!forum/pigx/>.

If you run into any bugs, please open an issue here:
<https://github.com/BIMSBbioinfo/pigx_rnaseq/issues>.
