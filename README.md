# Introduction to PiGx

![PiGx Logo]([http://bioinformatics.mdc-berlin.de/pigx/images/logo.svg)

The primary aim of PiGx was to provide genomics data processing tools that are relatively 1) easy to use, 2) easy to install, 4) easy distribute and most importantly 3) reproducible. Data processing is the major bottleneck for analyzing large number of samples. We want to make first-pass analysis as painless and as quick as possible for the users. So that, they will have more time to spend on data integration, visualization and statistical modeling.

All the pipelines have similar interface. For the end-user, each pipeline has the same input types: a sample sheet and a settings file. The sample sheet contains the information on samples such as names, locations of raw files, covariates etc. The settings file contains the extra arguments for the tools in the pipelines.  The users can in general run pipelines as follows: 

```
$./pigx [pipeline_name] [sample_sheet] -s [settings_file]
```

Additionally, we also provide users with high quality reports and figures that contains results from basic analyses and data quality check. The following chapters provide detailed documentation for the pipelines available in PiGx.




