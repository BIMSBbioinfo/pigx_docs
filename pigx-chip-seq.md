# PiGx ChIP-seq

## Introduction
PiGx-ChIPseq is a preprocessing and analysis pipeline that takes raw ChIP (Chromatin ImmunoPrecipitation) data and performs the necessary operations for analysis, including alignment and peak-finding. A final report is generated providing a summary for each sample provided by the user.


## The workflow
Figure ** provides an overview of the various stages of the pipeline, as well as the outputs and expected inputs.

![PiGx Chipseq workflow](./figures/pigx-chipseq_workflow.png)


Reads are aligned using [bowtie][bowtie], with quality control analysis 
performed with [fastqc][fastqc]
....

## Input preparation
To use the pipeline, the user must first edit two files: the sample sheet and the settings file. 

### Sample Sheet 

The sample sheet is a tabular file (`csv` format) describing the experiment.
The table has the following columns:

...

Here, each row of entries (below the header) corresponds to a sample, and the columns correspond to ...

 
### Settings File (`.yaml`)
In the settings file, various parameters are saved, in YAML format, to configure the execution of PiGx-chipseq. 
... describe them...

## Execution

After editting the above input files, the bsseq pipeline can be executed with
the command: `./pigx-chipseq [sample sheet] -s [settings file]` 
with further options...

## Analysis 
Where is the final report? where can necessary files be obtained for further analysis?

## Troubleshooting

What are some of the typical problems users might encounter, and how can they be solved?
(e.g. during cluster submission)?

