# CAP-CES

[View on OpenSAFELY](https://jobs.opensafely.org/repo/https%253A%252F%252Fgithub.com%252Fopensafely%252FCAP-CES)

Details of the purpose and any published outputs from this project can be found at the link above.

The contents of this repository MUST NOT be considered an accurate or valid representation of the study or its purpose. 
This repository may reflect an incomplete or incorrect analysis with no further ongoing work.
The content has ONLY been made public to support the OpenSAFELY [open science and transparency principles](https://www.opensafely.org/about/#contributing-to-best-practice-around-open-science) and to support the sharing of re-usable code for other subsequent users.
No clinical, policy or safety conclusions must be drawn from the contents of this repository.

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic
health records research in the NHS, with a focus on public accountability and
research quality.

Read more at [OpenSAFELY.org](https://opensafely.org).

# Licences
As standard, research projects have a MIT license. 



# Study details

*A Common Analytic Protocol to compare the safety and effectiveness of Covid-19 vaccines in England using OpenSAFELY* 

This repository contains analytic code for a Common Analytic Protocol, applicable to a chosen Covid-19 vaccination campaign in England, 
to make head-to-head comparisons between the vaccine products used in that campaign.

TODO: The protocol is available here:...

The Protocol accommodates the following campaign-specific characteristics:

* start and end dates
* vaccine products 
* study eligibility criteria

This repo should be forked (maybe?!) when starting an analysis for the next campaign. 

## Repository navigation

- The [`codelists/`](./codelists/) directory contains all the codelists used to define variables in analysis. 
- The [`analysis/`](./analysis) directory contains the executable scripts used to conduct the analysis. 
- The [`project.yaml`](./project.yaml) defines run-order and dependencies for all the analysis scripts.
**This file should *not* be edited directly**. To make changes to the yaml, edit and run the [`create-project.R`](./create-project.R) script instead.
- Non-disclosive model outputs, including tables, figures, etc, are available via the OpenSAFELY job server.

## Analysis scripts

The analysis scripts in the [`analysis/`](./analysis) directory are organised into sub-directories as follows:

- [`0-lib/`](./analysis/0-lib/):
  - [`study-dates.R`](./analysis/0-lib/study-dates.R) defines the key study dates (start date, end date, vaccine roll-out dates, etc) that are used throughout the study. 
  It create a json file that is used by R scripts and the study definition.
  - [`design.R`](./analysis/0-lib/design.R) defines the campaign-specific design elements used throughout the study (start and end dates, eligibility, products, etc). 
  It also defines matching and weighting specification, look-up dictionaries, and other useful objects. 
  This script is run at the start of all subsequent R scripts.
  - [`utility.R`](./analysis/0-lib/utility.R) defines functions used throughout the codebase. This script is run at the start of all subsequent R scripts.
- [`1-extract/`](./analysis/1-extract/):
  - [`dataset_definition.py`](./analysis/1-extract/dataset_definition.py) is the script defining the dataset to extract from the database, using ehrQL.
  - [`dummy_datasett_definition.R`](./analysis/1-extract/dummy_dataset_definition.R) defines a custom dummy dataset.  
  This can be used instead of the dummy data created by ehrQL when it is necessarily to have more control over the structure in the data, 
  such as more realistic vaccination dates or event rates.
  If the dataset definition is updated, this script must also be updated to ensure variable names and types match.
  - [`variables.py`](./analysis/1-extract/variables.py) contains some function and variable definitions to be read in by the dataset definition.
  - [`codelist.py`](./analysis/1-extract/codelists.py) pulls the codelists from the [`codelists/`](./codelists/) directory to be usable in the dataset definition. 
- [`2-prepare/`](./analysis/2-prepare/):
  - [`data_prepare.R`](./analysis/2-prepare/data_prepare.R) imports the extracted database data (or dummy data), standardises some variables and derives some new ones.
  - [`data_selection.R`](./analysis/data_selection.R) applies the inclusion criteria to the extracted data and creates a small table used for the inclusion/exclusion flowchart.
- [`3-adjust/`](./analysis/3-adjust/):
  - [`match.R`](./analysis/3-adjust/match.R) (`cohort`, `spec`) runs the matching algorithm to pair recipients of product A with product B, with matching criteria determined by `spec`. 
  It outputs a dataset containing the matching "weights" (`0`/`1`), and a matching ID. 
  - [`weight.R`](./analysis/3-adjust/weight.R) (`cohort`, `spec`) runs the propensity model to estimate the probability of receipt of product A versus product B, with the model determined by `spec`. 
  It outputs a dataset containing the person-specific weights. 
  - [`report.R`](./analysis/3-adjust/report.R) (`cohort`, `method`, `spec`) describes baseline information for the matched or weighted method
  eg Table 1 type cohort characteristics, post-weighting balance checks.
  - [`combine-weights.R`](./analysis/3-adjust/combine-weights.R) (`cohort`) combines weights across all weighted and matched analyses for the given cohort.
  Also calculated the Effective sample size based on the weights. 
  - [`match-coverage.R`](./analysis/3-adjust/match-coverage.R) (`cohort`, `spec`) describes matching rates over calendar time.
- [`4-constrast/`](./analysis/4-constrast/):
  - [`aj.R`](./analysis/4-constrast/aj.R) (`cohort`, `method`, `spec`, `subgroup`, `outcome`) derives Aalen-Johansen survival estimates for each product and calculates relative risk and risk differences. 
  This is largely based on the OpenSAFELY Kaplan-Meier reusable action, with an extension to AJ estimates.
  - [`plr.R`](./analysis/4-constrast/plr.R) (`cohort`, `method`, `spec`, `subgroup`, `outcome`) compares cumulative incidence curves between products using pooled logistic regression. 
  - [`combine-contrasts.R`](./analysis/4-contrast/combine-contrasts.R) collects treatment contrasts from the [`aj.R`](./analysis/aj.R) and [`plr.R`](./analysis/plr.R) scripts.

Scripts may take one or more arguments:

- `cohort`, the name of the cohort to be analysed, defined in the [`design.R`](analysis/0-lib/design.R) script.
- `spec`, the matching or weighting specification, taking values _A_, _B_, _C_, etc for convenience, and fully defined in the [`design.R`](analysis/0-lib/design.R) script.
For matching, `spec` is the set of variables to match on. For weighting, `spec` is the model formula passed to the `weightit()` function. 
- `method`, taking values _match_ or _weight_.
- `subgroup`, the subgroup variable. Cumulative incidences will be calculated separately within each level of this variable. 
Choose _all_ for no subgroups (i.e., the main analysis). Choose _<variable>_ to select a specific variable to stratify on.
This variable must be exactly matched in the matching run if using `approach="matching"`, and must be used as a stratification variable if using `approach="weighting"` (this requirement is under review!)
- `outcome`, the outcome of interest, for example _covidadmitted_ or _coviddeath_.

## Workspace

TODO 

## Pre-print

TODO


