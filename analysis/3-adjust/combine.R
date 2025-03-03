# # # # # # # # # # # # # # # # # # # # #
# Purpose: combine weights from all adjustment strategies
# imports weighting data and creates:
# - dataset containing all weights for each strategy
# - dataset containing effective sample size for each strategy
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")

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
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
}


# create output directories ----

output_dir <- here_glue("output", "3-adjust", cohort, "combine")
fs::dir_create(output_dir)

## import unadjusted cohort data ----
data_cohort <- read_feather(here("output", "2-prepare", cohort, "data_cohort.arrow"))

## create dataset of metaparameters to import
cohort0 <- cohort

metaparams <- 
  metaparams |>
  select(cohort, method, spec) |>
  unique() |>
  filter(
    cohort == cohort0,
    spec == "A"
  )

## create dataset that contains only patient IDs and the weights from all different adjustment strategies ----

data_weights <- 
  metaparams |>
  mutate(
    data = pmap(
      list(cohort, method, spec), 
      function(cohort, method, spec) {
        dat <-  
          here("output", "3-adjust", cohort, glue("{method}-{spec}"), "data_adjusted.arrow") |>
          read_feather() |>
          select(patient_id, treatment, weight)
        dat
      }
    )
  ) |>
  unnest(data) 


data_weights_wider <-
  data_weights |>
  pivot_wider(
    id_cols =  c(patient_id, treatment),
    names_from = c(cohort, method, spec),
    names_prefix = "wt_",
    names_sep = "_",
    values_from = weight
  ) |>
  mutate(
    wt_unadjusted = 1,
    .after = "treatment"
  )

data_all <- 
  left_join(
    data_cohort,
    data_weights_wider,
    by = c("patient_id", "treatment")
  )

write_feather(data_all, fs::path(output_dir, "data_weights.arrow"))

## create dataset of effective sample sizes for each adjustment strategy ----

table_ess <-
  data_weights |>
  group_by(treatment, cohort, method, spec) |>
  summarise(
    ess = (sum(weight)^2) / (sum(weight^2))
  ) |>
  pivot_wider(
    id_cols =  c(cohort, method, spec),
    names_from = treatment,
    names_prefix = "ess_",
    values_from = ess
  )
   
write_feather(table_ess, fs::path(output_dir, "table_ess.arrow"))

