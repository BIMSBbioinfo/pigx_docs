# PiGX ChIP-seq

# Introduction

PiGX ChIPseq is an analysis pipeline for preprocessing, peak calling and reporting for ChIP sequencing experiments. It is easy to use and produces high quality reports. The inputs are `fastq` files containing reads from the sequencing experiment, and a configuration file which describes the experiment. In addition to quality control of the experiment, the pipeline enables to set up multiple peak calling analysis and allows the generation of a UCSC track hub in an easily configurable manner.

## Workflow

PiGx ChIPseq implements best practices for preprocessing and analysis of ChIPseq data. Figure 1 provides an overview of the different steps of the pipeline, as well as the outputs. 

First, raw reads are trimmed using [TrimGalore!][trimgalore] to ensure a minimum read quality, and removal of adapter sequences. Next, reads are aligned to a reference genome using [Bowtie2][bowtie2], and peaks are called for any combination of samples using [MACS2][macs2]. The reproducibility of peaks among experiments is controlled using [IDR][idr]. 
Coverage `bigWig` files are calculated for the aligned reads using [rtracklayer][rtracklayer], and peak files are converted to `bigBed` format files using `bedToBigBed`, to be integrated together into a UCSC track hub to be viewed in the [UCSC Genome Browser](https://genome.ucsc.edu/).
Finally, quality metrics are calculated for the ChIP experiments, such as inter strand cross correlation and GC content, and signal profiles are calculated for genomic features using [genomation][genomation]. Together with the frequency of reads in peaks and the distribution of peaks in genomic features these metrics are used to compile a custom ChIP quality control report.


![PiGx ChIP-seq workflow](./figures/pigx-chipseq.svg)
_Figure 1: An overview of the PiGx ChIPseq workflow_

## Output

- `BAM` files
- `bigWig` files
- Peak files (`narrowPeak`, `broadPeak`)
- UCSC track hub folder
- Read QC reports
- ChIP QC Report

# Installation

You can install this pipeline with all its dependencies using GNU Guix:

    guix install pigx-chipseq

You can also install it from source manually.  You can find the
[latest
release](https://github.com/BIMSBbioinfo/pigx_chipseq/releases/latest)
here.  PiGx uses the GNU build system.  Please make sure that all
required dependencies are installed and then follow these steps after
unpacking the latest release tarball:

```sh
./configure --prefix=/some/where
make install
```

## Dependencies

By default the `configure` script expects tools to be in a directory
listed in the `PATH` environment variable.  If the tools are installed
in a location that is not on the `PATH` you can tell the `configure`
script about them with variables.  Run `./configure --help` for a list
of all variables and options.

You can prepare a suitable environment with Conda or with [GNU
Guix](https://gnu.org/s/guix).  If you do not use one of these package
managers, you will need to ensure that the following software is
installed:

<details>
<summary>Software dependencies</summary>

- R
    - argparser
    - biocparallel
    - biostrings
    - chipseq
    - data.table
    - dyplr
    - genomation
    - genomicalignments
    - genomicranges
    - rsamtools
    - rtracklayer
    - s4vectors
    - stringr
    - jsonlite
    - heatmaply
    - htmlwidgets
    - ggplot2
    - ggrepel
    - plotly
    - rmarkdown
- python
    - snakemake
    - wrapper
    - pyyaml
    - pytest
    - xlrd
    - magic
- pandoc
- fastqc
- multiqc
- trim-galore
- bowtie
- macs2
- idr
- samtools
- bedtools
- bedToBigBed
- bamToBed


</details>

### Via Guix

Assuming you have Guix installed, the following command spawns a
sub-shell in which all dependencies are available:

```sh
guix environment -l guix.scm
```

# Quick Start

To check wether the pipeline was properly installed, run PiGx ChIPseq on a minimal test dataset, which we provide [here](https://github.com/BIMSBbioinfo/pigx_chipseq/releases/download/v0.0.14/test_data.tar.gz). 

Once downloaded run these commands to extract the data:
```sh
mkdir test_dir
tar -xzf test_data.tar.gz --directory test_dir/
cd test_dir
```

Along with the data this tarball includes a settings file and a sample sheet file. 
The pipeline can now be started using this command, which takes about 3 minutes to finish (system with 4 cores and 32GB of RAM).
```sh
pigx chipseq -s test_dir/settings.yaml test_dir/sample_sheet.csv
```

Inside `test_dir` a new directory `out` is created, which includes specific directories containing output data for the respective step of the pipeline.

The `ChIP_Seq_Report.html` inside the `Reports` directory gives a general overview about the ChIP experiment and the peak calling. 


# Preparing Input

In order to run the pipeline, the user must supply

- a sample sheet
- a settings file

both files are described below.

In order to generate template settings and sample sheet files, type

```sh
pigx chipseq --init
```

in the shell, and a boilerplate `sample_sheet.csv` and `settings.yaml` will be written to your current directory.


## Sample Sheet

The sample sheet is a tabular file in `csv` format or an exel table (`xls`/`xlsx`) defining the samples used in any subsequent analysis. 

| SampleName | Read | Read2 |
|------|-------|--------|
| ChIP1| ChIP.fq.gz|  |
| Cont1| Cont.fq.gz|  |
| ChIP2| ChIP.fq.gz|  |
| Cont2| Cont.fq.gz|  |
|ChIPpe| ChIPpe_R1.fq.gz| ChIPpe_R2.fq.gz|

- SampleName is the name for the sample
- _Read/Read2_ are the fastq file names of paired end reads
  - the location of these files is specified in `settings.yaml`
  - for single-end data, leave the Read2 column in place, but have it empty


The creation of the sample sheet is straight forward using the following command: 
```sh
pigx chipseq --init=sample-sheet
```
This creates a template that can be filled with your own samples.

#### Technical Replicates

The sample sheet offers support for technical replicates, by repeating the sample name (first column) for different input files (second,third column).
The quality check will be performed for any input file and replicates will be merged during the mapping. 


| SampleName | Read | Read2 |
|------|-------|--------|
|ChIPpe| ChIPpe_R1.fq.gz | ChIPpe_R2.fq.gz |
|ChIPpe| ChIPpe_t2_R1.fq.gz | ChIPpe_t2_R2.fq.gz |

## Settings File

The settings file is a file in yaml format specifying general settings and the details of the analysis. 

The settings file is straight forward using the following command: 
```sh
pigx chipseq --init=settings
```

This settings file will have all lines commented out, you have to remove the trailing `#` from **Locations**, **General** and **Analysis** sections, then you can be modified it to fit your analysis.

It has the following **required** sections: 

### Locations

Here you define paths to be used in the pipeline, some of the items are required and some optional (can stay blank): 

| item    | required | description |
|---------|----------|-------------|
| _input-dir_ | yes | directory of the input files (`fastq` files) |
| _output-dir_    | yes | output directory for the pipeline |
| _genome-file_    | yes | path to the reference genome in `fasta` forma |
| _index-dir_    | no | directory containing pre-built mapping indices for the  given reference genome (created with `bowtie2-build`) |
| _gff-file_    | no | location of a `GTF` file with genome annotations for the  given reference genome |

#### Reference and annotations

The user must supply paths pointing to a reference genome for the organism, as well as annotated genes. These might for example be downloaded from:

- [UCSC](http://genome.ucsc.edu/cgi-bin/hgTables): The reference genome can be retrieved by selecting the _group_ `Mapping and Sequencing` with _track_  `Assembly`, then choosing _region_ `genome` and _output format_ `Sequence` and finally clicking _Get output_ . The gene annotations can be retrieved by selecting the _group_ `Genes and Gene Predictions` with _track_  `NCBI Refseq` or `Ensembl Genes`, then choosing _region_ `genome` and _output format_ `GTF - gene transfer format` and finally clicking _Get output_ . 

or 

- [Ensembl](https://www.ensembl.org/info/data/ftp/index.html), where the reference genomoe is listed under _DNA_ `FASTA`, the gene annotations under _Gene sets_ `GTF`. 

The user is free to choose any resource for these annotation files, however it is very important to note that the chromosome naming styles must match between different annotation files (e.g. chromosome 21 is named as *chr21* in UCSC style, but *21* in NCBI/Ensemble style). 

### General

These are settings which apply to all analysis (unless adjusted in single analysis):

| item    | required | description |
|---------|----------|-------------|
| _assembly_ | yes | version of reference genome (e.g. hg19,mm9, ...) |
| _params_    | no | list of default parameters for tools and scripts (for tools check respective manual for available parameters) |

### Execution

The `execution` section in the settings file allows the user to specify whether the pipeline is to be submitted to a cluster, or run locally, and the degree of parallelism. For a full list of possible parameters, see the execution section of the settings file created with `pigx chipseq --init=settings`.


A minimal settings file could look like this, but please consider that no analysis will be performed without adding [analysis information](#analysis-sections) :

```yaml
locations:
  input-dir: in/reads/
  output-dir: out/
  genome-file: genome/my_genome.fa
  index-dir:
  gff-file: genome/mm_chr19.gtf

general:
  assembly: hg19
  params:
    extend: 200
    scale_bw: 'yes'
    bowtie2:
        k: 1
    idr:
        idr-threshold: 0.1
    macs2:
        g: hs
        keep-dup: auto
        q: 0.05
    extract_signal:
        expand_peak: 200
        bin_num: 20

execution:
  submit-to-cluster: no
  jobs: 6
  nice: 19

```

### Analysis Sections

The analysis part of the setting file describes the experiment. It has following sections: 

| section | required | description |
|---------|----------|-------------|
| *peak_calling*  | yes | defines which samples will be used to detect regions of enriched binding ( multiple combinations and variations are possible, [see here for details](#peak-calling) ) |
| _idr_ | no | specifies pairs of *peak calling* analysis that are compared to determine the reproducibilty of the general experiment ([see here for details](#optional-idr)) |
| _hub_ | no | describes the general layout of a UCSC hub that can be created from the processed data and allows the visual inspection of results at a UCSC genome browser ([see here for details](#optional-hub)) |
| *feature_combination* | no | defines for a list of *peak calling* and/or *idr* analysis the combination of regions shared among this list ([see here for details](#optional-feature-combination)) |


The creation of these sections is straight forward considering the following snippets as template. Comments and examples within the snippets provide guidance of what is possible and what to take care of.

#### Peak Calling

The previously defined samples are used for subsequent peak calling analysis to detect regions of enriched binding. In this section any number of comparisons can be defined, while multiple combinations and variations are possible. In terms of peak calling the **ChIP** (also called treatment) is the sample in which we want to detect enriched regions compared to the **Cont(rol)** (or background) sample. Each analysis can be run with a unique set of parameters and default parameters for all analysis can be defined in the [settings file](#settings-file) , check available parameters and description [here](https://github.com/taoliu/MACS).
For more information have a look at the publication for the software we are using "Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137".  

```yaml
# define peak calling analysis
peak_calling:
    # analysis can have any name, but the names have to be unique 
    Peaks1: 
        # sample(s) to be used as treatment sample 
        ChIP: ChIP1
        # sample(s) to be used as control sample
        Cont: Cont1
        params:
            macs2:
                # each analysis can be adjusted independently
                # add/modify available parameters of the analysis
                nomodel: ''
                extsize: 300
    Peaks2:
        ChIP: ChIP2
        Cont: Cont2
        params:
            macs2:
                # each analysis can be adjusted independently
                nomodel: ''
                extsize: 147
    Peaks4:
        ChIP:
            # multiple samples can be used as treatment
            - ChIP1
            - ChIP2
        Cont:
            # multiple samples can be used as control
            - Cont1
            - Cont2
        params:
            macs2:
                nomodel: ''

    Peaks5:
        # the number of samples per group can differ 
        ChIP: ChIP2
        Cont:
            - Cont1
            - Cont2
        params:
            macs2:
                nomodel: ''

    Peaks6:
        # analysis can be performed without control
        ChIP: ChIP1
        Cont:
        params:
            macs2:
                nomodel: ''
```

#### (_optional_) IDR

Assuming that the some samples are (biological/technical) replicates, in order to measure the consistency between them use the irreproducible discovery rate (IDR)  "Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011). Measuring reproducibility of high-throughput experiments. The annals of applied statistics, 5(3), 1752-1779.", which is in general a good (but very stringent) quality control.

```yaml
idr:
    # idr analysis can have any name, but the names have to be unique 
    ChIP_IDR:
        # define the pair of samples, add more combinations for more replicates
        ChIP1: Peaks1
        ChIP2: Peaks2
```

#### (_optional_) Hub

In the _hub_ section the general layout of a [UCSC Track Hubs](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Intro) is described with some minimal items. The track hub is generated from the processed data and allows the visual inspection of results at a UCSC genome browser (for supported genomes). 

__IMPORTANT__: In order to use the hub feature, make sure you downloaded the reference genome in UCSC style (see [here](#reference-and-annotations) for instructions ).

The required items to define the hub are the following:

| item    | example | description |
|---------|----------|-------------|
| name    | PiGx_Hub | name of the hub directory |
| shortLabel    | PiGx_Short | short name of hub is displayed as name above track groups |
| longLabel    | PiGx_Hub_Long | descriptive longer label for hub is displayed as hub description |
| email    | my.mail[at]domain.com | whom to contact for questions about the hub or data |
| descriptionUrl    | pigx_hub.html | URL to HTML page with a description of the hub's contents |
| super_tracks    | see below | specification of hub layout (track groups, tracks) |

This is a small example how this could look like:

```yaml
hub:
    name: PiGx_Hub
    shortLabel: PiGx_Short
    longLabel: PiGx_Hub_Long
    email: my.mail@domain.com
    descriptionUrl: pigx_hub.html
    super_tracks:
        # track groups can have any name, but the names have to be unique 
        Tracks1:
            # tracks can have any name, but the names have to be unique 
            track11:
                # to add peaks as a track, define "type: macs" 
                name: Peaks1
                type: macs
            track12:
                # to add coverage signal as a track, define "type: bigwig"
                name: ChIP1
                type: bigWig
            # descriptive longer label for track group is
            # displayed as description in track settings
            long_label: Tracks1_long
```

#### (_optional_) Feature Combination

To find the combination of enriched binding regions, which is shared among a set of *peak calling* and/or *idr* analysis results, define a feature in the *feature_combination* section. Only items defined in the *peak_calling* and *idr* sections can be used here.  

```yaml
feature_combination:
    # features can have any name, but the names have to be unique
    Feature1:
        # define feature based on only one result
        - ChIP_IDR
    Feature2:
        # define feature based on more than one result
        - Peaks6
        - Peaks5
    Feature3:
        # define feature based on different analysis types
        - ChIP_IDR
        - Peaks5
```

# Running the pipeline

To run PiGx on your experimental data, first prepare the sample sheet and settings file (see [above](#preparing-input)), and then from the terminal type

```sh
$ pigx-chipseq -s settings.yaml sample_sheet.csv
```

If you are not sure wether you set everything up correctly, use the dryrun option which will only show what work would be performed, but does not actually run the pipeline.
```sh
$ pigx-chipseq -s settings.yaml sample_sheet.csv -n 
```

To see all available options type the `--help` option

# Output Description

PiGx ChIPseq creates an output folder, with a specific directory structure (for details see [here](#output-folder-structure)) that contains the following outputs and more.

### Quality control
General quality control metrics are computed using [FastQC][fastqc] and [MultiQC][multiqc]. The MultiQC report is particularly useful, collating quality control metrics from many steps of the pipeline in a single html report, which may be found under the `Reports` folder in the PiGx output folder.

A custom ChIP quality control report can be found in the `Reports` foler as well, presenting useful quality metrics for the ChIP experiments. 

### Formatted UCSC Hub

The `UCSC_Hub` directory contains a folder with the name specified in the **hub** section of the **analysis** part of the settings file. This folder might be copied to a web-accessible server to host the hub itself, which can then be accessed using the [UCSC Genome Browser][ucsc-genome-browser]. 

## Output Folder Structure

The pipeline will create a specific directory structure, 
the respective contents are explained below:
```
|-- Analysis
|-- Annotation
|-- BigWig
|-- Bowtie2_Index
|-- FastQC
|-- Log
|-- Mapped
|-- Peaks
|-- Reports
|-- Trimmed
|-- UCSC_HUB
```

| Folder    | description |
|---------|-------------|
|   Analysis     | Contains RDS files with intermediary analysis steps. RDS are binary files which efficiently store R objects. |
|    Annotation    | Formatted GTF annotation. |
|    BigWig    | Symbolic links to the bigWig signal files. |
|    Bowtie2_Index    | Processed genome file along with the Bowtie2_Index. |
|     FastQC   | FastQC sequencing quality report.  |
|     Log   | Detailed output from execution of each step of the pipeline. |
|   Mapped     | Mapped reads in .bam format, and corresponding bigWig files. |
|   Peaks     | Peaks called with MACS2. Depending on the parameters, contains either narrowPeak or broadPeak format. The file **sample_qsort.bed** contains uniformly processed peaks, sorted by their corresponding p value. |
|     Reports   | Contains MultiQC and ChIP quality reports in html format. |
|     Trimmed   | Trimgalore adaptor and quality trimmed files. |
|    UCSC_Hub    | Contains a completely formatted UCSC hub, with track descriptions, peaks and bigWig tracks. |

# Troubleshooting

## Execution on a cluster
Currently, PiGx only supports Sun Grid Engine for cluster execution. If you're uncertain about your cluster, try typing `qsub` in the shell (Sun Grid Engine uses `qsub` to submit jobs).

#### Disappearing jobs on the cluster
PiGx ChIPseq comes with sensible defaults for resource requests when running on a cluster, but based on the genome version and other parameters, these might not be sufficient and your cluster might terminate your jobs. The cluster resource requests may be overridden in the settings file. See the execution section of the settings file created with `pigx chipseq --init=settings`.

## Guix locale error
If you get a warning:

```
##########################

perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
    LANGUAGE = (unset),
    LC_ALL = (unset),
    LC_MEASUREMENT = "de_DE.UTF-8",
    LC_PAPER = "de_DE.UTF-8",
    LC_MONETARY = "de_DE.UTF-8",
    LC_NAME = "de_DE.UTF-8",
    LC_ADDRESS = "de_DE.UTF-8",
    LC_NUMERIC = "de_DE.UTF-8",
    LC_TELEPHONE = "de_DE.UTF-8",
    LC_IDENTIFICATION = "de_DE.UTF-8",
    LC_TIME = "de_DE.UTF-8",
    LANG = "en_US.UTF-8"
    are supported and installed on your system.
perl: warning: Falling back to the standard locale ("C").

##########################
```

and the pipeline breaks at multiqc, printing this into the Log :

```
####################################################
Traceback (most recent call last):
  File "/gnu/store/gaijkjaiv4hirsazalpbjp3ifymybkzb-multiqc-1.5/bin/.multiqc-real", line 767, in <module>
    multiqc()
  File "/gnu/store/2h1kcjw1r1306chd45452kwqzq2xb001-python-click-6.7/lib/python3.6/site-packages/click/core.py", line 722, in __call__
    return self.main(*args, **kwargs)
  File "/gnu/store/2h1kcjw1r1306chd45452kwqzq2xb001-python-click-6.7/lib/python3.6/site-packages/click/core.py", line 676, in main
    _verify_python3_env()
  File "/gnu/store/2h1kcjw1r1306chd45452kwqzq2xb001-python-click-6.7/lib/python3.6/site-packages/click/_unicodefun.py", line 118, in _verify_python3_env
    'for mitigation steps.' + extra)
RuntimeError: Click will abort further execution because Python 3 was configured to use ASCII as encoding for the environment.  Consult http://click.pocoo.org/python3/for mitigation steps.

Additional information: on this system no suitable UTF-8
locales were discovered.  This most likely requires resolving
by reconfiguring the locale system.

Click discovered that you exported a UTF-8 locale
but the locale system could not pick up from it because
it does not exist.  The exported locale is "en_US.UTF-8" but it
is not supported

#######################################################
```

Then you have to install `glibc-locales`

```
guix package -i glibc-locales
```
and export this global variable to use these guix locales:

```
export GUIX_LOCPATH=$HOME/.guix-profile/lib/locale
```

Or you might consider adding this directly to your `~/.bashrc` or `~/.bash_profile`. 

```
echo -e "# use guix locales\nexport GUIX_LOCPATH="'$HOME'"/.guix-profile/lib/locale"  >>  ~/.bashrc
```


# FAQ

__Q:__ I get the following error:
```bash
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
/home/agosdsc/projects/pigx_chipseq/test_dir/out
If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
```
What happend and what should I do?

__A:__ The pipeline crashed at some point, possible reasons are mentioned by the error. Do as the error message proposes and pass the `--unlock` argument once and then run the pipeline again without `--unlock`.  

__Q:__ command name?
In my installation, the command appears to be pigx-chipseq and not pigx chipseq as the docs say.

__A:__ The command that you have to use depends on the way the pipeline was installed, either as part of the pigx bundle (`guix install pigx`) or as a single pipeline (`guix install pigx-chipseq`)

__Q:__ The docs mention how to deal with technical replicates, but how should biological replicates be handled?

__A:__ Usually biological replicates are kept seperately, which means you should give them distinct names in the sample sheet to allow for distinct alignments, and then for peak-calling you pick up the sample names from the sample sheet to reenter them in the 'analysis' section of the settings file.

__Q:__ I get an error when trying to run the pipeline:
```sh
SystemExit in line 26 of /gnu/store/ik8kz81k720a8sy8nc4119jy5njb98gd-pigx-chipseq-0.0.20/libexec/pigx_chipseq/scripts/Check_Config.py:
ERROR: Config file is not properly formated:
Genome fasta headers contain whitespaces.
 Please reformat the headers

  File "/gnu/store/ik8kz81k720a8sy8nc4119jy5njb98gd-pigx-chipseq-0.0.20/libexec/pigx_chipseq/Snake_ChIPseq.py", line 38, in <module>
  File "/gnu/store/ik8kz81k720a8sy8nc4119jy5njb98gd-pigx-chipseq-0.0.20/libexec/pigx_chipseq/scripts/Check_Config.py", line 26, in validate_config
``` 

__A:__ Our pipeline does not like if in the FASTA headers of the reference genome there is more information than the chromosome name, so please use this oneliner to format your genome-file: 
```sh
cat refGenome.fa  | awk '{ print $1}’ >  refGenome_noWhiteSpace.fa
```


# Questions
If you have further questions please e-mail:
pigx@googlegroups.com or use the web form to ask questions
https://groups.google.com/forum/#!forum/pigx/



[trimgalore]: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
[bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/
[macs2]: https://github.com/taoliu/MACS
[idr]: https://github.com/nboley/idr
[genomation]: http://bioinformatics.mdc-berlin.de/genomation/
[rtracklayer]: https://bioconductor.org/packages/release/bioc/html/rtracklayer.html
[fastqc]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[multiqc]: http://multiqc.info/
[ucsc-genome-browser]: https://genome.ucsc.edu/
