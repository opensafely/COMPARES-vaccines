This directory contains scripts and files that define elements of the study design that are used in multiple places throughout the project, such as study entry dates, matching and adjustment variables, and outcome variables.

The design.R script should be sourced within the main analysis scripts to load in objects (lists of dates, look-up tables, etc), rather than being run on it's own.

The study dates are stored in the study-dates.json file, rather than being defined in design.R script, because they are needed by both the dataset definition and the R scripts. 

