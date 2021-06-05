# PiGx BS-seq

# Introduction

PiGx-bsseq is a preprocessing and analysis pipeline that takes raw
`fastq` read files and performs all necessary steps to present the
full methylome for analysis. Quality control, and differential
methylation are also performed, and a final report provides a summary
for each sample provided by the user.

## Workflow

PiGx-bsseq follows best practices for processing and analysis of bisulfite sequencing.
Figure 1 provides an overview of the various stages of the
pipeline, as well as the outputs and expected inputs.

First reads are trimmed for quality and adapter sequences using 
[TrimGalore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).  
Then reads will be mapped to the genome using one or both of the available bisulfite aware aligners 
[Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) or 
[bwa-meth](https://github.com/brentp/bwa-meth), before
alignments are filtered for duplicate reads with [samblaster](https://github.com/GregoryFaust/samblaster) 
and sorted with samtools. 
If the reference genome the reads are going to be mapped to has not already
undergone bisulfite conversion, such a conversion will be prepared
automatically by PiGx BSSeq.  
Methylation-calling is then carried out using either [MethylKit](https://github.com/al2na/methylKit)
or [methylDackel](https://github.com/dpryan79/MethylDackel) depending on the chosen aligner, 
before post-mapping analysis such as segmentation and differential
methylation analysis between samples will be performed. 
Finally we generate multiple reports including sample report for each sample,
a differential analysis report and the pan-sample multiqc report.    

![PiGx BSseq workflow](./figures/pigx-bsseq.svg)
_Figure 1: An overview of the PiGx BSseq workflow_

# Output

- Quality Control reports
	- per Sample report shows Coverage and Metyhlation distribution
	- MultiQC report informs about Read Quality, 
	  Trimming Performance and Alignment Statistics
- Alignment results in BAM file format. 
- bigwig files for genome-wide methylation fraction tracks
- BED files for genome-wide methylation segmentation
- Intermediated objects from methylKit analysis in tabix format
- Differential methylation analysis results
    - comprehensive HTML report
    - Tab-separated table for differential methylation results effect size and significance 
    - BED files for detected genome-wide differentially methylation cytosines (DMCs)
    - BED files for aggregated regions of genome-wide differentially methylation (DMRs)


# Installation

Pre-built binaries for PiGx are available through 
[GNU Guix](https://gnu.org/software/guix), the
functional package manager for reproducible, user-controlled software
management.  Install the complete pipeline bundle with the following
command:

```sh
guix install pigx
```

You can install this BSseq pipeline with:


```sh
guix install pigx-bsseq
```

    

You can also install it from source manually.  You can find the
[latest
release](https://github.com/BIMSBbioinfo/pigx_bsseq/releases/latest)
here.  PiGx uses the GNU build system.  Please make sure that all
required dependencies are installed and then follow these steps after
unpacking the latest release tarball:

```sh
./configure --prefix=/some/where
make install
```

# Quick start

1. Download the zipped test data:

    `wget https://github.com/BIMSBbioinfo/pigx_bsseq/releases/download/v0.0.8/test-data.tar.gz`

2. Unzip the archive

    `tar xzvf test-data.tar.gz`

3. Run the pipeline

    `cd tests && pigx-bsseq -s settings.yaml sample_sheet.csv`


# Preparing the input

In order to run the pipeline, the user must supply

- a sample sheet
- a settings file

both files are described below.

In order to generate template settings and sample sheet files, type

```bash
pigx bsseq --init
```

in the shell, and a boilerplate `sample_sheet.csv` and `settings.yaml` will be written to your current directory.


## Sample Sheet

The sample sheet is a tabular file (`csv` format) describing the experiment.
The table has the following columns:


| Read1            | Read2            | SampleID    | Protocol   | Treatment |
|------------------|------------------|-------------|------------|-----------|
| SA_R1.fastq.gz   | SA_R2.fastq.gz   | sample_A    |WGBS        | H2O       |
| SB.fastq.gz      |                  | sample_B    |WGBS        | MED1      |
| ...              |                  |             |            |           |

_Table 1: example sample sheet_

### Column descriptions

- _SampleID_ is the name for the sample, which must be unique to each row of the table. 
- _Read1/2_ are the fastq file names of paired end reads
  - the location of these files is specified in `settings.yaml`
  - for single-end data, leave the _Read2_ column in place, but have it empty
- _Protocol_ refers to the experimental protocol, can be either "WGBS" / "RRBS" for whole-genome/reduced-representation bisulfite sequencing, respectively
- _Treatment_ contains a treatment label that is used for reference in differential methylation


> NOTE: The _Protocol_ decides wether the pipeline will perform deduplication ("WGBS") or not ("RRBS"). 

 
## Settings File

In the settings file, various parameters are saved, in YAML format, to
configure the execution of PiGx-bsseq. 

### Locations

  - The path to directory containing the reads  ( where the input `fastq` files are)
  - The path to a desired output directory, which PiGx will create if necessary
  - The path to the directory containing the fasta reference genome.
  
### General

  - The reference genome name (e.g. 'hg19') 
  - wether to use bismark or bwa-meth aligner, or both
  - detailed settings which only need to be adjusted by experienced users
  
#### Methylation-calling Settings
  
  - minimum number of read-hits required for calling methylation (10)
  - minimum mapping quality required for calling methylation (10)

#### export-bigwig:
  - decide for each cytosine methylation context wether merge strands and wether to export to bigwig
    
#### differential-methylation:
  - number of cpu threads for differential methylation (1)
  - significance threshold for detection of DMC
  - effect size for detection of DMC
  - annotation files are optional for the report, but 
  - `CpGfile`:  file (with path) specifying CpG annotations for differential methylation (genome/cpgIslandExt.hg19.bed.gz).
  - `refGenfile`: file (with path) specifying reference genes for differential methylation (genome/refGene.hg19.bed.gz).
  - `webfetch`: Boolean instruction as to whether PiGx should attempt to fetch the above files from the internet if they are unavailable locally (no).


### Differential Methylation (DM) analyis
  - specify which of the treatment groups defined in the sample sheet should be compare against each other


The section described here comprise most of the settings that a 
typical user might want to access, while more advanced settings provided
as defaults in _execution_ section will not need to be modified by most
users. 


```yaml
general:
  assembly: ''
  # use one or both bisulfite aligners
  # bwameth is faster and uses less resources
  use_bwameth: True
  # bismark is gold standard and uses sensible defaults  
  use_bismark: True
  methylation-calling:
    minimum-coverage: 10
    minimum-quality: 10
    # this applies to bwameth only
    keep-singleton: False
    keep-Dups: Auto                 # can be Auto/True/False, Auto decides based on protocol
  export-bigwig:
    context:
      cpg: 
        export: True
		# join both strands and summarize coverage/methylation? 
        destrand: True
      chg: 
        export: False
        destrand: False
      chh: 
        export: False
        destrand: False
  reports:
    TSS_plotlength: 5000
  differential-methylation:
    cores: 1
    qvalue: 0.01
    difference: 25
    annotation:
      # download cpgIsland-bedfile or refGenes-bedfile automatically if not given?
      cpgIsland-bedfile: ''
      refGenes-bedfile: ''
      webfetch:   no


DManalyses:
  # The names of analyses can be anything but they have to be unique
  # for each combination of case control group comparisons.

  analysis1:
    treatment_sample_groups: "MED1"
    control_sample_groups: "H2O"


 multipleTreat:
    # If multiple sample groups are provided, they must be separated by comma.
   treatment_sample_groups: "MED2,MED1"
   control_sample_groups: "H2O"

 withinGroup:
   # For within group comparison (all vs all), give the same groups to treatment and control.
   treatment_sample_groups: "MED2,MED1,H2O"
   control_sample_groups: "MED2,MED1,H2O"V
```

# Running the pipeline

After editting the above input files, the bsseq pipeline can be executed with
the command: `pigx-bsseq [sample sheet] -s [settings file]` with further
options, described in the README (for example, the flag `--dry-run` can be
added to confirm valid input settings before execution).


## Analysis 
In the `out-dir` folder, PiGx will create a sub-directory
`Final_Reports/` in which a final report for each sample is provided, along
with a report on the differential methylation observations for each pair of
samples that was provided in the sample sheet under `general:
differential-methylation: treatment-groups:` Further analysis can be performed
on the aligned data using the sorted `.bam` file in `06_sorting`, or the bigwig
file in `08_bigwig_files`, although for simple methylation percentages, it
should be suffice to consult the column-separated files in `07_methyl_calls`.
In each case, the specific file names will refer to the sample IDs provided in
the original sample sheet.

# Troubleshooting

The dependency graph of rules contains branches; as such, the rules are not
always performed in _exactly_ the order indicated by directory prefix labels.
Nevertheless, in cases of interrupted calculation, it may be useful to inspect
the log files from the last (or second-to-last) directory created to check for
error messages. 

Alongside these directories, the directory `pigx_work` is also created, with
its contents described in the file `CONTENTS`---here one can see links to the
original data files for traceability. In case a run of PiGx was done a long
time ago, under conditions that have been forgotten, the subfolder
`pigx_work/input/` contains links to each of the raw data files,
`pigx_work/refGenome` links to the reference genome that was mapped do during
exection, and `pigx_work/bin` contains links to the versions of binary
executables that were used in the process.

### Cluster Submission

The settings file allows for fields that will serve to optimize resource usage
when PiGx is submitted to a cluster via SGE (Sun Grid Engine, a standard
queueing system). For example, under the section `execution: rules:` various
parameters can be supplied regarding resource consumption (e.g. number of
threads, memory required, etc.) for each rule individually, or, alternatively,
`execution: cluster: queue` can specify a single queue for all jobs. 

The default parameters were chosen and tested on the human hg19 genome, since
it is larger (and therefore requires a larger RAM footprint) than most other
model organisms (e.g. Drosophila, mouse, etc.).  This is not an exhaustive
test, however, and it remains possible that PiGx users studying exotic genomes
may yet exceed resource allowance. In such cases, a job execution will be
aborted and an email will be sent to the address supplied by the user detailing
circumstances of the failed job (such as start time, stop time, `Max vmem` =
maximum memory consumed), and resource requests can be adjusted accordingly.
Advanced users may consult `etc/settings.yaml.in` to see which values can be
fine-tuned; novice users who are studying either the human genome, or a smaller
one, should simply ignore this section and allow PiGx to assign the default
values.

# Questions

If you have further questions please send e-mail to
pigx@googlegroups.com or [ask questions on the web
forum](https://groups.google.com/forum/#!forum/pigx/).

