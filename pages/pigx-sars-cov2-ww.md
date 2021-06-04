# PiGx SARS-CoV-2 wastewater

# Introduction

PiGx SARS-CoV-2 is a pipeline for analysing data from sequenced wastewater samples and identifying given variants-of-concern of SARS-CoV-2. Currently wastewater samples are used, which are enriched for SARS-CoV-2. The pipeline can be used for continuous sampling. The output of the PiGx SARS-CoV-2 pipeline is summarized in a report which provides an intuitive visual overview about the development of variant abundance over time and location. Additionally there will be more detailed reports per sample, which cover the quality control of the samples, the detected variants and a taxonomic classification of all unaligned reads. This version of the pipeline was designed to work with target sequencing (PCR-based) and tested with data generated using the [ARTIC nCoV-2019 primers](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V30).

## Workflow

First the raw reads are trimmed by using [Prinseq](http://prinseq.sourceforge.net/) to improve alignment rates and mutation calling. Next, the trimmed reads are aligned to the reference genome of SARS-CoV-2 using [BWA](https://github.com/lh3/bwa), and the results are *SAM*/*BAM* files of **aligned** and **unaligned reads**. Following the alignment a quality check on raw and processed reads is performed by using [MultiQC](https://multiqc.info/). Calling the variants and inferring SNVs (single nucleotide polymorphisms) on the **aligned reads** is done with [LoFreq](https://csb5.github.io/lofreq/). Mutations are annotated with [VEP](https://covid-19.ensembl.org/index.html). Estimation of Variants of Concern (VOC) frequencies are done by deconvolution. 
To investigate the abundance of other existing species in the wastewater samples the **unaligned reads** will be taxonomicly classified with [Kraken2](https://github.com/DerrickWood/kraken2). The Kraken2 requires a database of the genomes against the reads are getting aligned, therefore keep in mind that you can only find those species which are included in the chosen database. For documentation how to set this up, see: [Prepare databases](#prepare-databases). For a better and interactive visualization of all species present in the wastewater [Krona](https://github.com/marbl/Krona/wiki) is used. Also here a small step of setting up a database is needed before running the pipeline, see: [Prepare databases](#prepare-databases). 

## Output

* overview report including:
   * Visualization of the development of SARS-CoV-2 variants and mutations over time and locations from all samples provided.
   * Quality Control report per sample: number of covered amplicons, read coverage and MultiQC report for raw and trimmed reads
   * Variants report per sample: variant analysis of SARS-CoV-2 from each wastewater sample and identification of variants of concern
   * Taxonomic classification: a table and pie chart of the species found in the unaligned reads
* *SAM* / *BAM* files per sample: aligned and unaligned reads against SARS-CoV-2 
* *VCF* / *CSV* files per sample:  listing all detected single nucleotide variants (SNVs) from the aligned reads
* VEP reports per sample as *TXT* / *HTML*: report files show the [VEP output](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#output)
  including the uploading variants, consequenting amino acid changes and consequences for the corresponding protein
* Kraken2 files per sample *(txt)*: provides an overview of all found species in the unaligned reads together with the NCBI taxonomy ID
* log files for all major analysis steps performed by the pipeline


# Installation

Pre-built binaries for PiGx are available through [GNU Guix](https://gnu.org/s/guix), the functional package manager for reproducible, user-controlled software management. 
You can install the PiGx SARS-CoV-2 pipeline with

```sh
guix install pigx-sars-cov2-ww
```

If you want to install PiGx SARS-CoV-2 from source, please clone this repository and change directory accordingly:

```sh
git clone https://github.com/BIMSBbioinfo/pigx_sarscov2_ww.git
cd pigx_sarscov2_ww
```

To fetch code that is common to all PiGx pipelines run this:

```sh
git submodule update --init
```

Before setting everything up, though, make sure [all dependencies](https://github.com/BIMSBbioinfo/pigx_sarscov2_ww/blob/main/manifest.scm) are met by either installing them manually, or by entering the provided reproducible Guix environment. If you are using Guix we definitely recommend the latter. This command spawns a sub-shell in which all dependencies are available at exactly the same versions that we used to develop the pipeline:

```sh
USE_GUIX_INFERIOR=t guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

To use your current Guix channels instead of the fixed set of
channels, just omit the `USE_GUIX_INFERIOR` shell variable:

```sh
guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

Note that `--pure` unsets all environment variables that are not
explicitly preserved.  To access other executables that are not part
of the environment please address them by their absolute file name.

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh
./configure --prefix=/some/where
make install
```

At this point you are able to run PiGx SARS-CoV-2. To see all available options type `--help`.

```sh
pigx-sars-cov2-ww --help
```

## Prepare databases

Before the pipeline can work, three databases must be downloaded and their location will need to be provided in the settings file. Depending on the size of the databases this can take some time.
Be sure that the pigx-sarscov2-ww pipeline is downloaded and the tools are installed or used via the provided and suggested guix environment. One database (signature mutations, `sigmut_db`) is already provided via the repository. The directory structure is suggested like [this](#structure-overview) and pre-filled accordingly in the settings file. 

### Kraken2 database

There are several libraries of genomes that can be used to classify the (unaligned) reads. It is up to you which one to use, but be sure that they fulfill the necessities stated by Kraken2 [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases). For an overall overview we recommend to use the Plus-PFP library provided [here](https://benlangmead.github.io/aws-indexes/k2). If the classification is not of concern or only the viruses are of interest, we recommend to use a smaller one. This will accelerate the speed.
It is also possible to have multiple Kraken2 databases installed, just be sure to provide the correct location to the settings file.

First download and unpack the database in the `databases/kraken_db/`:

```
DIR=databases/kraken_db/
mkdir -p $DIR

# NOTE: This command will download a very large file, after unpacking this will
# require about 100GB of disc space. If this is not feasible use another
# database instead. For this please see link above or commented lines below.
wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz | tar -C $DIR -xzv

# Use the following two lines to use smaller 8GB version instead.
# wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_8gb_20210127.tar.gz | tar -C $DIR -xzv
```

Next go to the database/ directory and build the Kraken database. This might take a while (depends on the size of the downloaded database):

```
#current location: pigx_sarscov2_ww/
cd databases/
DBNAME=kraken_db
kraken2-build --use-ftp --download-taxonomy --db $DBNAME # if this fails, you might want to try it without the --use-ftp flag
kraken2-build --build --db $DBNAME
```

Kraken might tell you that it can't find a library subdirectory. If that's the case it should be fine though.  
This folder should now contain at least these files: hash.k2d, opts.k2d and taxo.k2d.

### Krona database

Krona Tools needs two files, which have to be installed in the `databases/krona_db/`. Also this might take a while:

```
#current location: pigx_sarscov2_ww/databases/
DBNAME=krona_db/
mkdir -p $DBNAME
KRONA=$(dirname $(which ktImportTaxonomy))/../share/krona-tools/ # this is just a workaround until tool paths are declared
$KRONA/updateTaxonomy.sh $DBNAME # the scripts are stored a priori in that folder
$KRONA/updateAccessions.sh $DBNAME
```

### VEP database

Just download the `SARS_CoV_2` database for VEP (variant effect predictor) and unpack it in the `databases/vep_db/` directory.

```
#current location: pigx_sarscov2_ww/databases/
DBNAME=vep_db/
mkdir -p $DBNAME
wget -qO- ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz | tar -C $DBNAME -xzv
```

### Sigmut database

Necessary files are provided in `databases/sigmut_db/` for the current main Variants of Concern. Users can add files with new variants if/when necessary.



# Preparing the input

In order to run the pipeline, you need to supply

- a sample sheet
- a settings file.

Both files are described below.

In order to generate template settings and sample sheet files, type

```sh
pigx-sars-cov2-ww --init
```

in the shell, and a boilerplate `sample_sheet.csv` and `settings.yaml` will be written to your current directory. An example for both files is provided in the `tests/` directory.

## Sample sheet

The sample sheet is a tabular file (`csv` format) describing the experiment. The table has the following columns:

| SampleName | Read               | Read2              | date               | location_name      | coordinates_lat    | coordinates_long   |
| ------     | -------            | --------           | --------           | --------           |  --------          |  --------          |
| Test0      | Test0_R1.fastq     | Test0_R2.fastq     | 2021-01-01T08:00:00| Berlin             | 52.3646054650197   | 13.5098643274129   |
| Test2      | Test2_R1.fastq     | Test2_R2.fastq     | 2021-01-03T08:00:00| Munich             | 48.2084486780314   | 11.6282300407931   |

- _SampleName_ is the name for the sample
- _Read_ & _Read2_ are the fastq file names of paired end reads
  - the location of these files is specified in the settings file
  - single-end data is not yet supported
  - compressed (`.gz`) reads are not yet supported
- _date_ is a date/time in ISO format (`yyyy-mm-ddThh:mm:ss`)
- _location_name_ is the name of the location and should be unique per coordinates
- _coordinates_lat_ & _coordinates_long_ correspond the latitude and longitude of the location name

## Settings file

The settings file contains parameters (in YAML format) to configure the execution of the PiGx SARS-CoV-2 pipeline. It specifies:

**Locations**:

- _output-dir_, the location of the outputs for the pipeline
- _reads-dir_, the location of the reads (directory where `fastq` files are)
- _reference-fasta_, the `fasta` file with the reference genome (must be prepared by the user)
- _amplicons-bed_, the amplicons `bed` file for coronavirus (must be prepared by the user)
- _kraken-db-dir_, the location of the kraken database (must be prepared by the user)
- _krona-db-dir_, the location of the krona database (must be prepared by the user)
- _sigmut-db-dir_, the location of the signature mutations database (provided at databases/sigmut_db/)
- _vep-db-dir_, the location of `sars_cov_2` database for VEP (must be prepared by the user)

**Trimming**:

These settings are used to filter raw reads when trimming. Reads
that are shorter than the product of read-length and the cut-off factor are removed.

- _read-length_ specifies the length of the basepairs
- _cut-off_ specifies the cut-off factor



# Quick Start

To check wether the pipeline together with the databases was properly installed, run PiGx SARS-CoV-2 Wastewater Sequencing Pipeline on a minimal test dataset.
For this there are samples provided in `tests/sample_data/`. The directory structure should be provided like this, assuming all databases are set up like described [here](#prepare-databases):

```
pigx_sarscov2_ww
│
├── databases
│   ├── kraken_db
│   │   └── ...
│   ├── krona_db
│   │   └── ...
│   ├── sigmut_db
│   │   └── ...
│   └── vep_db
│       └── ...
├── tests
│   ├── databases
│   │   └── vep_db
│   │       └── ...
│   ├── sample_data
│   │   ├── reads
│   │   │   ├── ...
│   │   ├── signature_mutations
│   │   │   └── ...
│   │   ├── NC_045512.2.fasta
│   │   ├── nCoV-2019_NCref.bed
│   │   └── ...
│   ├── sample_sheet.csv
│   ├── settings.yaml
│   └── ...
└── ...
``` 

Now the the test set can be run with the command: 

```
#current location: pigx_sarscov2_ww/
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww -s tests/settings.yaml tests/sample_sheet.csv
```

Inside `tests/` a new directory `output` is created, which includes specific directories containing output data for the respective step of the pipeline.    The `tests/output/reports/index.html` gives the overview over all merged reports for the test data. 

# Troubleshooting

If you have any questions please e-mail: pigx@googlegroups.com or use the web form to ask questions https://groups.google.com/forum/#!forum/pigx/. 

If you run into any bugs, please open an issue here: https://github.com/BIMSBbioinfo/pigx_rnaseq/issues. 

