# PiGx BS-seq

## Introduction
PiGx-bsseq is a preprocessing and analysis pipeline that takes raw `fastq` read files and performs all necessary steps to present the full methylome for analysis. Quality control, and differential methylation are also performed, and a final report provides a summary for each sample provided by the user.


## The workflow
Figure ** provides an overview of the various stages of the pipeline, as well as the outputs and expected inputs.

![PiGx BSseq workflow](./figures/pigx-bsseq_workflow.png)

In addition to fastq read files, a necessary input to the pipeline is a
reference genome to be mapped to.  If such a genome has not already undergone
bisulfite conversion, such a conversion can be prepared automatically by PiGx.
At the same time, reads are trimmed for quality and adapter sequences using
[TrimGalore!][trimgalore], with quality control analysis applied both before
and after.  Once these steps are completed, reads can then be mapped to the
genome using [Bismark][bismark], before alignments are filtered for duplication
and sorted.  Methylation-calling is then carried out using
[MethylKit][Methylkit] before initial post-mapping analysis such as
segmentation and differential methylation between samples.

## Input preparation
To use the pipeline, the user must first edit two files: the sample sheet and the settings file. 

### Sample Sheet (`.xls` -> `.csv`)

The sample sheet is a tabular file (`csv` format) describing the experiment.
The table has the following columns:

| reads            | reads(2)         | sample_ID   | protocol   | treatment |
|------------------|------------------|-------------|------------|-----------|
| SA_R1.fastq.gz   | SA_R2.fastq.gz   | sample_A    |WGBS        | 0   |
| SB.fastq.gz      |                  | sample_B    |WGBS        | 1   |
| ...              |                  |             |            |     |

Here, each row of entries (below the header) corresponds to a sample. In this
example, `sample A` is a paired-end experiment with two `fastq` files, whereas
`sample B` is single end; thus, the `read 2` column for the latter is left blank.
The third column, `sample ID`, contains some descriptive name (without
white spaces) for the sample. The fourth column refers to the experimantal
protocol (e.g. "WGBS" / "RRBS" for whole-genome/reduced-representation
bisulfite sequencing, respectively), while the 5th column contains a treatment
label--generally an integer-- that is used for reference in differential
methylation (see below).

The above sample sheet can be produced by editing the included file
`test/sample_sheet.xls` and saving it in `.csv` format with fields separated by
commas (`,`).

 
### Settings File (`.yaml`)
In the settings file, various parameters are saved, in YAML format, to configure the execution of PiGx-bsseq. 
The entries are stored hierarchically, and a brief description is provided below, along with possible default value in parenthesis:

- `locations`: 
  - `input-dir`: path to directory containing the input `fastq` files (./in)
  - `output-dir`: path to desired output directory, which PiGx will create if necessary (./out)
  - `genome-dir`: path to directory containing the reference genome. (./genome)

- general:
  -  `assembly`: reference genome name ('hg19') 
  -  methylation-calling:
    -  `minimum-coverage`: minimum number of read-hits required for inclusion in analysis (10).
    -  `minimum-quality`:  minimum read quality required for mapping (10)
  -  differential-methylation:
    - `cores`: number of cpu threads for differential methylation (1)
    - `treatment-groups`: a series of pairings specifying comparisons to be performed between samples, based on the treatment vector supplied in the sample sheet (no default, but an example compatible with the above sample sheet is provided in the following line)
        - ['0', '1']
  -    annotation:
    - `CpGfile`:  file (with path) specifying CpG annotations for differential methylation (genome/cpgIslandExt.hg19.bed.gz).
    - `refGenfile`: file (with path) specifying reference genes for differential methylation (genome/refGene.hg19.bed.gz).
    - `webfetch`: Boolean instruction as to whether PiGx should attempt to fetch the above files from the internet if they are unavailable locally (no).


The values stored in the settings input file are used to overwrite the default
settings, described below. An example settings file is available in
the `test/` directory as `settings.yaml`. 
The fields in `locations`, for example, will almost certainly require user input
to over-write the default locations.

Note in particular the field `differential-methylation`:`treatment-groups`; below
this there may be arbitrarily lines in the format `- ['A', 'B']` where `A` and
`B` are integers referring to the treatment values from the sample sheet. Each
such line represents a command to carry-out pair-wise comparison between the
samples with the corresponding treatment values `A` and `B`. i

The remaining fields provided here (and in the example settings file)
 comprise most of the settings that a
typical user might want to access, while the more basic settings that are saved
as defaults in `etc/settings.yaml.in` will not need to be modified by most
users (although advanced users may freely re-define any such variables in their own
settings file, and by so doing, overwrite them).

## Execution

After editting the above input files, the bsseq pipeline can be executed with
the command: `./pigx-bsseq [sample sheet] -s [settings file]` with further
options, described in the README (for example, the flag `--dry-run` can be
added to confirm valid input settings before execution).  Once the pipeline is
executed, the desired output directory is created, along with various
sub-directories for intermediate steps in the process, each of which are named
with prefixes to indicate their (approximate) order in sequence.  In most
cases, these processes have their own interim reports and log files. For
example, `06_sorting` contains the sorted `.bam` file after alignments, and
`07_methyl_calls` contains information on average methylation in various
formats, while `01_raw_QC` is the earliest step --quality control of the raw
inputs. 

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

## Troubleshooting

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

