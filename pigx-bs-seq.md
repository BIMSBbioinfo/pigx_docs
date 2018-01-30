

# PiGx BS-seq

## Preparation
PiGx-bsseq processes raw fastq read files and generates a final report for each sample provided by the user. To use the pipeline, the user must first edit two files: the sample sheet and the settings file. 

### Sample sheet (`.xls` -> `.csv`)
The sample sheet can be produced by editing the included file `test/sample_sheet.xls` and saving it in `.csv` format with fields separated by commas (`,`). In this file, each row of entries (below the header) corresponds to a sample; the first column contains the filename of a samples fastq input file, while the second column should contain the second fastq input file (if it exists, in the case of paired-end reads; otherwise this column should be left empty). The third column should contain a sample ID, some descriptive name (without white spaces) for the sample. The fourth column describes the type of read (e.g. "WGBS" for whole-genome bisulfite sequencing), while the 5th column contains a treatment label---generally an integer that is used for reference in differential methylation (see below). 

### Settings file (`.yaml`)
In the settings file, various parameters are saved, in YAML format, to configure the execution of PiGx-bsseq. The values stored in the settings input file are used to overwrite the default settings, which are described below. An example settings file is available in the `test/` directory as `settings.yaml`. Here, you can see several fields that _will_, in many cases, require editing on the part of the user. For example, the field `locations` defines the path to the directories containing the input files, the desired output files, and the reference genome being mapped to. Note in particular the field `differential-methylation`:`treatment-groups`; below this there may be arbitrarily lines in the format 
`- ['A', 'B']`
where `A` and `B` are integers referring to the treatment values from the sample sheet. Each such line represents a command to carry-out pair-wise comparison between the samples with the corresponding treatment values `A` and `B`. The remaining fields in the example settings file comprise most of the settings that a typical user might want to access, while the more basic settings that are saved as defaults in `etc/settings.yaml.in` will not need to be modified by most users (although users may freely re-define any such variables in their own settings file, and by so doing, overwrite them).

## Execution

After editting the above input files, the bsseq pipeline can be executed with the command:
`./pigx-bsseq [sample sheet] -s [settings file]`
with further options, described in the README (for example, the flag `--dry-run` can be added to confirm valid input settings before execution).
Once the pipeline is executed, the desired output directory is created, along with various sub-directories for intermediate steps in the process, each of which are named with prefixes to indicate their (approximate) order in sequence, and usually have their own interim reports and log files. For example, `06_sorting` contains the sorted `.bam` file after alignments, and `07_methyl_calls` contains information on average methylation in various formats, while `01_raw_QC` is the earliest step. 

## Troubleshooting

The dependency graph of rules contains branches; as such, the rules are not always performed in _exactly_ the order indicated by directory prefix labels. Nevertheless, in cases of interrupted calculation, it may be useful to inspect the log files from the last (or second-to-last) directory created to check for error messages. 

Alongside these directories, the directory `pigx_work` is also created, with its contents described in the file `CONTENTS`---here one can see links to the original data files for traceability. In case a run of PiGx was done a long time ago, under conditions that have been forgotten, the subfolder `pigx_work/input/` contains links to each of the raw data files, `pigx_work/refGenome` links to the reference genome that was mapped do during exection, and `pigx_work/bin` contains links to the versions of binary executables that were used in the process.

# Final report

Finally, the directory `Final_report` is created with the final report of the pipeline, it consolidates information from the reports in the other steps, and provides links to files storing more comprehensive data on the samples.
