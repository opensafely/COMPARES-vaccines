# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe matching results
# imports matching data
# reports on matching coverage, matching flowcharts, creates a "table 1", etc
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('survival')
library('gt')
library('gtsummary')
library('cobalt')

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "age75plus"
  weightset <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  weightset <- args[[2]]
}



# create output directories ----

output_dir <- here_glue("output", "3-cohorts", cohort, "weight{weightset}", "report")
fs::dir_create(output_dir)

## import unadjusted cohort data ----
data_cohort <- read_feather(here("output", "3-cohorts", cohort, "data_cohort.arrow"))

## import weighting info ----
data_weights <- read_feather(here_glue("output", "3-cohorts", cohort, "weight{weightset}", "data_weights.arrow"))


# table 1 style baseline characteristics in the weighted pseduo population ----

tab_summary_adjusted <-
  data_cohort |>
  left_join(
    data_weights |> select(patient_id, weight),
    by = "patient_id"
  ) |>
  select(
    treatment, 
    weight,
    any_of(names(variable_labels))
  ) |>
  mutate(
    N = 1L,
    treatment_descr = fct_recoderelevel(as.character(treatment), recoder$treatment),
  ) %>%
  table1_svysummary(
    group = treatment_descr, 
    weight = weight,
    labels = variable_labels[names(variable_labels) %in% names(.)], 
    threshold = sdc.limit
  )

write_csv(tab_summary_adjusted, fs::path(output_dir, "table1.csv"))



# Balance table ----

## this is almost equivalent to:
# bal.tab(
#   obj_weightit,  
#   stats = c("mean.diffs"), 
#   s.d.denom = "pooled", 
#   binary = "std", 
#   continuous = "std"
# )
# But we don't use this approach because 
# (a) we don't need to save and pass the `obj_weight` object between actions
# (b) we get consistency with the approach used in the matching script 
#     which, due to parallelisation, can't do `bal.tab(obj_matchit)` directly
# (c) we can also look at balance for additional variables not used in the propensity model
#     without faffing around with the `add1` and `data` arguments


balance <- 
  data_cohort |>
  select(
    any_of(names(variable_labels))
  ) |>
  bal.tab(
    weights = data_weights$weight,
    treat = data_weights$treatment,
    distance = data_weights$ps,
    stats = c("mean.diffs"), 
    disp = c("m", "sd"),
    un = TRUE,
    s.d.denom = "pooled", 
    binary = "std", 
    continuous = "std"
  ) %>%
  `[[`("Balance")


write_csv(balance, fs::path(output_dir, "data_balance.csv"))

# TODO: fix variable names in balance data so that they can be joined on the tbl_summary datasets
