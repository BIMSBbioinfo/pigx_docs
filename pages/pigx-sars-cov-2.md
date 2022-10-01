# PiGx SARS-CoV-2

# Introduction

PiGx SARS-CoV-2 is a pipeline for detecting viral lineages in sequencing data
obtained from enriched wastewater samples. It was developed with SARS-CoV-2 as
its target virus, but other targets are theoretically possible. The viral
lineages are provided by the user together with their characteristic signature
mutations. The pipeline is very flexible, allowing the user to choose from
multiple input and output files, and giving them fine control about parameters
used by the individual tools. PiGx SARS-CoV-2 has been developed with a focus on
reproucible outputs, and it can be used for continuous sampling. The output of
the PiGx SARS-CoV-2 pipeline is summarized in a report which provides an
intuitive visual overview about the development of lineage abundance and single
significantly increasing mutations over time and location. In addition, the
pipeline will generate more detailed reports per sample, which cover the quality
control of the samples, the detected variants, and a taxonomic classification of
all unaligned reads. This version of the pipeline was designed to work with
paired-end amplicon sequencing data e.g. following the ARtIc protocols
[ARTIC nCoV-2019 primers](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3),
but single-end sequencing reads are supported as well.

## Workflow

In the first step the pipeline takes the raw reads and the additional
information about used primers and adapters to perform extensive quality
control. Primer trimming is done with
[iVAR](https://github.com/andersen-lab/ivar), and
[fastp](https://github.com/OpenGene/fastp) is used for adapter trimming and
quality filtering. Next, the trimmed reads are aligned to the reference genome
of SARS-CoV-2 using [BWA](https://github.com/lh3/bwa), and the results are
*SAM*/*BAM* files of **aligned** and **unaligned reads**. Following the
alignment a quality check on raw and processed reads is performed by using
[MultiQC](https://multiqc.info/). Furthermore samples are checked for genome
coverage and how many of the provided signature mutation sites are covered.
Based on this every samples gets a quality score. Samples with genome coverage
below a user defined percentage threshold are reported as discarded samples, as
they are not included in time series analysises and summaries.

Calling the variants and inferring single nucleotide polymorphisms (SNVs) on the
**aligned reads** is done with [LoFreq](https://csb5.github.io/lofreq/).
Mutations are annotated with [VEP](https://covid-19.ensembl.org/index.html).
Estimation of lineage frequencies is done by deconvolution (see Methods in the
realted publication for details). To investigate the abundance of RNA matching
other existing species in the wastewater samples the **unaligned reads** will be
taxonomicly classified with [Kraken2](https://github.com/DerrickWood/kraken2).
Kraken2 requires a database downloaded locally of the genomes against the reads
are getting aligned. For documentation how to set this up, see:
[Prepare databases](#preparing-the-databases). For a better and interactive
visualization of all species present in the wastewater
[Krona](https://github.com/marbl/Krona/wiki) is used. Also here a small step of
setting up a database is needed before running the pipeline, see:
[Prepare databases](#preparing-the-databases). Interactive reports are generated using
[R-markdown](https://rmarkdown.rstudio.com/) and
[plotly for R](https://plotly.com/r/) for visualizations.

### Pooling of samples for time series analysis and plots

For summarizing across daytime and location, the lineage frequencies are pooled
by calculating the weighted average using the total number of reads of each
sample as weights (missing samples are removed).

## Output

* Overview reports including:
  * Summary and visualization of the development of SARS-CoV-2 variants and
    mutations over time and locations from all samples provided.
  * Quality Control reports per sample.
  * Per sample variant report: variant analysis of SARS-CoV-2 from each
    wastewater sample and identification of variants of concern.
  * Taxonomic classification of unaligned reads: Overview over taxa inferred
    from all the sequencing reads not aligning to the reference genome.
* Deconvoluted variant abundances
* Mutation abundances
* Single nucleotide variation call files
* VEP reports
* Kraken2 taxonomic classifications
* Numerous intermediate read, alignment, and statistics files
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

## Preparing the databases

Before the pipeline can work, three databases must be downloaded to a location
specified in the settings file. Depending on the size of the databases this can
take some time.

Without any user intervention, this will happen automatically via snakemake
rules. This behaviour is controlled via parameters in the settings file. See the
[Settings file](#settings-file) section for details.

Alternatively, the databases may be downloaded maunally via the
`download_databases.sh` scripts accessible like so:

```sh
prefix="$(dirname pigx-sars-cov-2)/../"
$prefix/libexec/pigx_sars-cov-2/scripts/download_databases.sh
```

However, the `download_databases.sh` script does not offer the flexibility of
downloading the databases automatically with user defined parameters, unless it
is manually edited.

*Note: The directory that the databases will be downloaded to needs to match the
database directories given in the settings file, else the pipeline will download
the databases to the given directory again, unnecessarily using up space. This
should be the case by default though.*

Read on for details if you want to download the databases manually.

### Kraken2 database

There are several libraries of genomes that can be used to classify the
(unaligned) reads. It is up to you which one to use, but be sure that they
fulfill the necessities stated by Kraken2
([Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases)).
For an overall overview we recommend to use the Plus-PFP library provided
[here](https://benlangmead.github.io/aws-indexes/k2), which is also the default
library used in the pipeline. If the classification is not of concern or only
the viruses are of interest, we recommend using a smaller one. This will speed
up the pipeline.

After downloading and unpacking the database files, use `kraken2-build` to
download the taxonomy data and build the database.

### Krona database

The way we use Krona, we only need the taxonomy database, as downloaded via
their `updateTaxonomy.sh` script.

### VEP database

For our use of VEP in the pipeline, we need a pre-indexed cache of the VEP
database of transcript models. This is the main point of the pipeline that
determines which virus can be analysed with the pipeline. Currently, VEP only
has data on SARS-CoV-2. But the VEP cache is the only point strictly determining
which virus the pipeline can deal with.

Per default the pipeline uses the indexed cache archive at
`http://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz`,
which only needs to be unpacked to the target directory. Currently this is the
only available chache file.

# Quick start

To check whether the pipeline and the databases have been properly set up, run
the pipeline on a minimal test dataset. If you installed the pipeline from
source, start at step 3.

1. Download the test data

    `git clone https://github.com/BIMSBbioinfo/pigx_sars-cov-2 sarscov2-test`

2. Enter the directory

    `cd sarscov2-test`

3. Run the pipeline using a preconfigured settings file. If you installed the
   pipeline from source, the database location will be set to whatever dir was
   specified during the configure step. If you installed the databases manually,
   you will need to also adjust the database paths in the test settings file
   accordingly.

    `pigx-sars-cov-2 -s tests/setup_test_settings.yaml tests/sample_sheet.csv`

Inside `tests/` a new directory `output_setup_test` is created, which includes
specific directories containing output data for the respective step of the
pipeline. The `tests/output_setup_test/reports/index.html` gives the overview
over all merged reports for the test data.

# Preparing the input

In order to run the pipeline, you need to supply

* Sample sheet (CSV format): containing information about sampling date and  
  location
* Settings file (YAML format) for specifying the experimental setup and
  optional custom parameter adjustments  
* Mutation sheet containing the lineages of interest and their signature
  mutations in nucleotide notation (CSV format)  
* Mutation BED file containing the genomic coordinates of the mutation sites
  (see below for details)  
* Reference genome of the target species in fasta format (so far the
  pipeline is only optimized for SARS-CoV-2, others might work too but not yet  
  tested)
* Primer BED file containing the PCR primer locations (e.g the primers suggested
  from ARTIC protocols)

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
  * in the case of single-end data, leave the `Read2` column **empty**
* *date* is a date/time in ISO format (`yyyy-mm-ddThh:mm:ss`)
* *location_name* is the name of the location and should be unique per
  coordinates
* *coordinates_lat* & *coordinates_long* correspond the latitude and longitude
  of the location name

## Mutation sheet

The mutation sheet should contain one column of siganture mutations per lineage
that shall be tracked and analysed by deconvolution. They should be given in the
format `GENE:RxxxV`, where `GENE` is the name of the gene in which the mutation
is found, `R` is a string of reference nucleotides in upper case, `xxx` is the
first reference nucleotides position (a number), and `V` is a string of variant
nucleotides, also in upper case. There is no upper or lower limit for the number
of signature mutations per lineage. However, please note that the deconvolution
results are more robust and precise with a higher number of mutations. (Tested
with 10-30 mutations per lineage).

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

The primer file contains the locations of primer sequences on the reads, i.e. it
defines where in the genomes primer sequences *may* occurr. This is determined
by the primer scheme used in generating the sequencing reads. It is required by
iVar for primer trimming. An example file can be found in the test directory
(`nCoV-2019_NCref.bed`).

## Settings file

The settings file contains parameters (in YAML format) to configure the
execution of the PiGx SARS-CoV-2 pipeline. There are generally two settings
files at play:

1. The pipeline interal settings file (source: `etc/settings.yaml`, installed:
   `$prefix/share/pigx_sars-cov-2/settings.yaml`), containing default settings.
   This file is not supposed to be modified by the user.
2. An analysis specific, user provided settings file containing changes to the
   default settings. These are generally paths to the input files, databases,
   etc.

When the pipeline is executed, both settings files are combined into a run
specific config file (`config.json`), used by the internal `snakemake` workflow
manager. It is always generated in the place the `pigx-sars-cov-2` program was
called and will be overwritten on subsequent calls.

<details>
  <summary>Detailed settings explanation.</summary>

### `locations`

Paths to various input files and directories needed by the pipeline.

* *output-dir* output directory for the pipeline.
* *input-dir* direcotry containing the input files, the files therin should
  match the file suffix given under control/start.
* *reference-fasta* [Mutation sheet](#mutation-sheet)
* *primers-bed* [Primer BED file](#primer-bed-file)
* *mutations-bed* [Mutation BED file](#mutation-bed-file)
* *mutation-sheet* [Mutation table](#mutation-sheet)
* *kraken-db-dir* [Kraken2 database](#kraken2-database)
* *krona-db-dir* [Krona database](#krona-database)
* *vep-db-dir* [VEP database](#vep-database)

### `databases`

Settings controlling the download and subsequent processing of databases needed
for the pipeline. If the databases are already present, none of these settings
will have any effect.

When prebuilt archives are to be used, the pipeline can download them both via
`ftp` and `http`/`https`. If no protocol is included, `http`/`https` is assumed.

*Note: Generally the official databases will be used except when running
`make check`/`make distcheck` without the databases pre-installed at the default
location. In that case the settings file at `tests/settings.yaml` will be used,
which configures the pipeline to use database archives prebuilt by us. This is
to ensure the github actions can run smoothly. In order to build the archives,
at least 1GB of disk space is necessary, which is not given on a github action
runner.*

#### `kraken2`

The kraken2 database download is the most complex of the three database
downloads. The download is done via one of the database archive file provided at
the [Index Zone](https://benlangmead.github.io/aws-indexes/k2), except in the
case outlined above.  
When no downsampling should occur, the database archive is extracted and used
as-is. In the other case, everything except the hash file (`hash.k2d`) is
extracted, the taxonomy is downloaded, and the database is built on disk. The
taxonomy files are removed afterwards as they are not used again and only take
up space after database building, at least for our purposes.

##### `archive-url`

The url the kraken2 archive will be downloaded from. The default archive is the
very large database including protozoa and fungi in addition to the standard
database.

When prebuilt archives are to be used, the pipeline can download them both via
`ftp` and `http`/`https`. If no protocol is included, `http`/`https` is assumed.  

Default: https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz

##### `downsample-db`

Whether or not downsampling of the `hash.k2d` file should occurr.  

Default: false

##### `max-db-size-bytes`

If downsampling should occurr, what is the maximum allowable size in bytes? The
resulting file may be smaller than the given value. The value is passed on
directly to `kraken2-build` under the `--max-db-size` option. Will only have an
effect of `downsample_db` is true.

Default: 250000000

#### `krona`

The Krona database download is fairly simple. Either the Krona internal download
script `updateTaxnomy.sh`
([Documentation](https://github.com/marbl/Krona/wiki/Installing)) is used, or a
prebuilt archive is downloaded.

##### `use-prebuilt`

Whether or not a prebuilt archive should be downloaded instead of running the
update taxonomy script.  
  
Default: false

##### `archive-url`

If a prebuilt archive should be downloaded, this tells the pipeline where to
find it.

When prebuilt archives are to be used, the pipeline can download them both via
`ftp` and `http`/`https`. If no protocol is included, `http`/`https` is assumed.

Default: ""

#### `vep`

The vep download is the simplest. The database is always in a compressed
archive.

##### `archive-url`

Location from where the archive will be downloaded.

When prebuilt archives are to be used, the pipeline can download them both via
`ftp` and `http`/`https`. If no protocol is included, `http`/`https` is assumed.
  
Default: ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz

### `parameters`

#### `vep`

Parameters for rule vep. See documentation about vep arguments `--species`,
`--buffer_size`, and `--distance` respectively in the vep
[documentation](https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html).

##### `species`

Needs to match with the downloaded VEP database cache. The `vep` default is
"homo_sapiens", which is unlikely to be of use for this pipeline.
  
Directly passed to the `vep` parameter `--species`.
  
Default: sars_cov_2

##### `buffer-size`

Number of variants in memory at one time. Trades off between run time and memory
usage. The higher, the faster.
  
Directly passed to the `vep` parameter `--buffer-size`.
  
Default: 5000

##### `transcript-distance`

Up- and downstream distance of a gene and a transcript which will be classified
as an up- or downstream variant.
  
Directly passed to the `vep` parameter `--distance`.
  
Default: 5000

##### `db-version`

This specifies a database version, which is needed when the database version
differs from the `vep` executable version. By default, `vep` looks for a
database version matching its own. This is needed even when using a offline
database cache like this pipeline does. Therefore this needs to be adjusted
when changing the `vep` database used.
  
Directly passed to the `vep` parameter `--distance`.
  
Default: 101

#### `ivar_trimming`

Parameters for rule `ivar_primer_trim`. See documentation about ivar arguments
`-q`, `-m`, and `-s` respectively in the iVar
[documentation](https://andersen-lab.github.io/ivar/html/manualpage.html).
  
`ivar` trimms primers from alignment BED files by first removing the sections
listed in a second primer BED file, and then performing quality trimming by
sliding a window along each read from 5' to 3' end, clipping the read if the
average window quality drops below a given cutoff.

##### `quality-cutoff`

If the average base call quality in the sliding window drops below this value,
the read will be trimmed to the last base of sufficient average window quality.  

Directly passed to the `ivar` parameter `-q`.  

Default: 15

##### `length-cutoff`

Read length threshold, if a read is shorter than this after trimming, it will be
discarded.  

Directly passed to the `ivar` parameter `-m`.  
Default: 30

##### `window-width`

Number of bases in the sliding window, i.e. the window width.
  
Directly passed to the `ivar` parameter `-qs`.
  
Default: 4

#### `reporting`

These paramerters are used to set quality control filters for the reports.

##### `mutation-coverage-threshold`

Results from samples without sufficient coverage measures are not included in
the visualizations or the linear regression calculations.
  
Default: 90

#### `deconvolution`

##### `method`

Control the deconvolution method used.
Possible options are:

* `rlm`:
  Robust linear regression as implemented in the `MASS` package, executed
  via the `deconvR` package.
* `nnls`:
  Non-negative least squares (`deconvR`).
* `qp`:
  Quadratic programming (`deconvR`).
* `svr`:
  Support vector regression (`deconvR`).

Prepending "weighted_" to the chosen method will add weighting of bulk mutation
abundances per variant by the inverse of the proportion of detected mutations to
known mutations. This biases the abundance of variants with high proportions
towards higher values, and coversely biases the abundence regression is of
variants with low proportions towards lower values.
  
Default: weighted_rlm

##### `mutation-depth-threshold`

Minimum sequencing per mutation for it to be used in the analysis.
  
Default: 100

## `control`

The following settings control from which point the pipeline starts, which
path it follows (i.e. which rules are executed on the way) and where it stops.

##### `start`

Start points for the analysis

* fastq(.gz): Raw reads in (gzipped) FASTQ format
* bam: Unfiltered alignments
* vcf: Variant calling files. *Note: Ideally these should have been generated by
  lofreq. The minimum requirement is that the INFO fields "AF" and "DP" are
  present. (More details on [here](https://samtools.github.io/hts-specs/))*

When using a non-default target, not all rules will be able to run and it won't
be possible to generate all reports. The main reasons for this are the missing
quality control statistics that are necessary mainly for the final report, and
the missing unalinged read files when starting after the alignment step.

Default: fastq.gz

##### `targets`

Desired results of the analysis. Multiple targets at the same time is possible.

* help: Print all rules and their descriptions.
* final_reports: Produce a comprehensive report. This is the default target.
* devonvolution: Run deconvolution for all provided samples and create a summary
  table containing abundances from all samples.
* lofreq: Call variants and produce .vcf file and overview .csv file.
* multiqc: Create MultiQC reports for including raw and trimmed reads.
  
Default: - final_reports

##### `run-ivar-primer-trimming`

Whether the primer trimming step should be performed.
  
Default: yes

#### `execution`

Settings directly related to the program execution itself.

##### `jobs`

The number of jobs allowed to be executed at one time. If executed locally this
specifies the maximum number of cores used at a time. May not extend to the
number of jobs used by individual tools. See section `--cores` & `--jobs` in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options)

Default: 6

##### `submit-to-cluster`

Whether cluster specific settings will be respected, and a qsub submission
call should be executed.
  
Default: no

##### `cluster`

Settings specific to executing the pipeline on a computing cluster. None of
these are relevant when `submit-to-cluster` is "no".

###### `missing-file-timeout`

How long before a rule output file is declared missing and the pipeline stopped.

When executing the pipeline on a cluster, file system latency can be higher
than when the pipeline is executed on a single server. Therefore the time
snakemake waits before declaring a rule output file as missing and stopping
needs to be adjusted.
  
Default: 120

###### `stack`

Stack memory used for each rule. We recommend leaving it as it is, unless you
really know what you are doing.
  
Default: 128M

###### `queue`

The name of a specific queue the whole pipeline should be submitted to.
  
Default: all

###### `contact-email`

How the cluster administration may contact you.
  
Default: none

###### `args`

Additional arguments passed to `qsub`, as a single string.
  
Default: ''

##### `rules`

Per rule submission settings. Give the rule name as a heading to configure
that specific rule.

###### `__default__`

Default settings used in absence of rule specific settings.
  
* `threads`: Number of threads. Default is 1.
* `memory`: RAM available for the rule, with unit (e.g. "90K", "5M"). In the
  case of no unit, unit M is assumed. Default is 4G.

### `tools`

Overwrites for locations of specific tools used. Each tool has its own
subheading, i.e. "bwa", and the following settings:

#### `executable`

Path to the executable file for the tool, defaults to the system installation.

#### `arguments`

Additional arguments to the tool as one string. Only use this if you know
what you are doing, these may cause conflicts with the arguments supplied in
each rule.

</details>

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

<details>
  <summary>All outputs with their locations.</summary>

*Given locations are relative to the output directory. `SAMPLE` indicates a
variable part of a path that will be replaced with a sample name. `VIRUS`
will be replaced with the virus being investigated.*

* Overview reports:
  * Summary and visualization of the development of SARS-CoV-2 variants and
    mutations over time and locations from all samples provided
    (`report/index.html`).
  * Quality Control reports per sample:
    * Overall QC report with number of covered amplicons, read coverage, etc
      (`report/SAMPLE.qc_report_per_sample.html`).
    * MultiQC report per sample for raw and trimmed reads along with several
      supplementary files (`report/multiqc/SAMPLE/*`).
    * FASTQC reports per sample and read for raw reads, adapter trimmed reads,
      and aligned reads, along with an archive of supplementary files
      (`report/fastqc/SAMPLE/*`).
    * FASTP reports per sample on adapter trimming statistics
      (`fastp/SAMPLE/*`).
  * Per sample variant report: variant analysis of SARS-CoV-2 from each
    wastewater sample and identification of variants of concern
    (`report/SAMPLE.qc_report_per_sample.html`).
  * Taxonomic classification of unaligned reads: Overview over taxa inferred
    from all the sequencing reads not aligning to the reference genome
    (`report/SAMPLE.taxonomic_classification.html`).
* Alignments:
  * Per sample aligned reads: Reads aligned against the
    reference genome, at various stages of trimming and sorting, along with
    their index files (`mapped_reads/SAMPLE_aligned.(bam|sam|bai)`).
  * Per sample unaligned reads: Reads not aligning to the
    reference genome (`mapped_reads/SAMPLE_unalingned.(bam|fastq)`).
* SNV files:
  * Per sample files listing all detected single nucleotide variants (SNVs) from
    the aligned reads (also parsed to csv, `variants/SAMPLE_snv.(vcf|csv)`).
* VEP reports:
  * Per sample report files constituting the
    [VEP output](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#output)
    including the variants uploaded to the VEP database, resulting amino acid
    changes and consequences for the corresponding protein.
    * Raw VEP output (`variants/SAMPLE_vep_VIRUS.txt`).
    * Report with detailed run statistics
      (`variants/SAMPLE_vep_VIRUS.txt_summary.html`).
    * Warnings generated by vep during the run
      (`variants/SAMPLE_vep_VIRUS.txt_warnings.txt`)
    * File with most of the raw output parsed into a comma separated table for
      downstream processing (`variants/SAMPLE_vep_VIRUS_parsed.csv`).
* Deconvoluted variant abundances:
  * Per sample abundance of each of the VOCs as
    * a regular table (`variants/SAMPLE_variant_abundance.csv`),
    * in a single row together with metadata (used later to construct summary,
      `variants/SAMPLE_variants_with_meta.csv`).
  * Overall summary table of per sample variant abundances with
    sample metadata (`variants/data_variant_plot.csv`).
* Mutation abundances:
  * Per sample abundance of each signature mutation found and meeting the read
    depth threshold criterium as a single row together with metadata (used later
    to construct summary, `mutations/SAMPLE_mutations.csv`).
  * Overall summary table *CSV* file of per sample mutation abundances with
    sample metadata.
  * Mutation data of non-signature mutations meeting the read depth threshold
    (`mutations/SAMPLE_non_sigmuts.csv`).
* Sample quality:
  * Per sample tables of read depth at each signature mutation locus, derived
    from the *untrimmed* alignment (`coverage/SAMPLE_genome_cov.tsv`).
  * Per sample tables giving coverage statistics of the *untrimmed* alignment
    (`coverage/SAMPLE_mut_cov.tsv`).
  * Per sample tables aggregating the statistics of genome and mutation coverage
    files in a sample quality table (single row for downstream use,
    `coverage/SAMPLE_quality.csv`).
  * Overall table concatenating the per sample quality rows into one table
    (`coverage/sample_quality_table.csv`).
* Adapter-trimmed reads (`trimmed_reads/SAMPLE_trimmed.fastq.gz`).
* Kraken2 raw output: provides an overview of all found species in the unaligned
  reads together with the NCBI taxonomy ID
  (`kraken/SAMPLE_classified_unaligned_reads.txt`).
* Per sample Krona diagram of taxa proportions in unaligned reads
  (`report/SAMPLE.Krona_report.html`).
* Mutation count summary table containing count statistics for mutations
  overall and per sample (`mutation_counts.csv`).
* Raw table of per mutation linear model coefficients and their p-values,
  generated by the `mutation_regression.R` script
  (`unfiltered_mutations_sig.csv`).
* Sample summary table containing various statistics about each sample
  (`overview_QC.csv`).
* Log files for all major analysis steps performed by the pipeline
  (`logs/*.log`).

</details>

# Troubleshooting

If you have any questions please e-mail: pigx@googlegroups.com or use the web
form to ask questions <https://groups.google.com/forum/#!forum/pigx/>.

If you run into any bugs, please open an issue here:
<https://github.com/BIMSBbioinfo/pigx_rnaseq/issues>.
