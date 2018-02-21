# PiGx ChIP-seq

# Introduction 
PiGx-ChIPseq is a preprocessing and analysis pipeline that
takes raw data from ChIP (Chromatin ImmunoPrecipitation) experiments and
produces illuminating figures.  The pipeline performs quality control on
experimental data, multiple peak-finding/alignment, and UCSC track hub
generation. A single final report is generated for each sample provided by the
user summarizing a standard set of key findings, while processed data is made
available for further analysis.

## The workflow
## Figure *1* provides an overview of the various stages of the
pipeline, as well as the outputs and expected inputs.

![PiGx_Chipseq_workflow](./figures/pigx-chipseq_workflow.png)

[Writing the workflow out would be helped very much by a single snakemake script that shows the order of operations. At the moment, the DAG output is convoluted and difficult to follow]

Reads are aligned using [bowtie][bowtie], with quality control analysis
performed with [fastqc][fastqc] ....

---

# Preparation 
To use the pipeline, the user must suply two files:

- a sample sheet
- a settings file

## Sample Sheet 

The sample sheet is a tabular file (`csv` format) describing the experiment.
The table has the following columns:

[Need to finally clarify what the sample sheet is going to look like]

...

Here, each row of entries (below the header) corresponds to a sample, and the
columns correspond to ...

 
## Settings File (`.yaml`) In the settings file, various parameters are saved,
in YAML format, to configure the execution of PiGx-chipseq. Values marked with
an * are required; the remainder are optional (i.e. and can therefore  be left
blank).  For a full list of possible parameters, see `etc/settings.yaml.in` to
review the default values that will be supplied by PiGx --all values supplied
in this file can be over-written, should the user supply an alternate value
in their settings file.
 
Altogether, these parameters are organized into the following sub-categories:

### Locations
Defines paths to be used in the pipeline. :
  - input-dir*: 	directory of the input files (fastq files)
  - output-dir*: 	output directory for the pipeline
  - genome-file*: 	path to the reference genome in fasta forma
  - index-dir: 	directory containing pre-built mapping indices for the given reference genome (created with bowtie2-build)
  - gff-file: 	location of a GTF file with genome annotations for the given reference genome

### General
The section contains settings that apply to all analysis (unless adjusted in single analysis):
  - assembly*: 	version of reference genome (e.g. hg19,mm9, ...)
  - params: 	list of default parameters for tools and scripts (for tools check respective manual for available parameters)

### Execution 
The execution section in the settings file allows the user to
specify how the pipeline should be executed: e.g., whether the pipeline is to
be submitted to a cluster, or run locally, and the degree of parallelization.
Furthermore, for cluster runs, advanced users can specify resource allocation
for specific rules. Further possible parameters include, e.g. an 
email address for notification of aborted jobs. 
A typical settings file might look like the following example:

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
  rules:
    __default__:
      queue: all.q
      memory: 8G
    bowtie2:
      queue: all.q
      memory: 16G
  cluster:
    contact-email: example@domain.de
    missing-file-timeout: 120
```

## Execution

After editting the above input files, the bsseq pipeline can be executed with
the command: `./pigx-chipseq [sample sheet] -s [settings file]` with further
options...

## Analysis 
Where is the final report? What does it look like? where can necessary files be obtained
for further analysis?

## Troubleshooting

What are some of the typical problems users might encounter, and how can they
be solved?  (e.g. during cluster submission)?

