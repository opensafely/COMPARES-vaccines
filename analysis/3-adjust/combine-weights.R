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
data_cohort <- read_feather(here("output", "2-select", cohort, "data_cohort.arrow"))

## create dataset of metaparameters to import
cohort0 <- cohort

metaparams_cohort_method_spec <- 
  metaparams |>
  select(cohort, method, spec) |>
  unique() |>
  filter(
    cohort == cohort0,
  )

## weights ----
## create dataset that contains only patient IDs and the weights from all different adjustment strategies

data_weights <- 
  metaparams_cohort_method_spec |>
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

## ESS ----
## create dataset of effective sample sizes for each adjustment strategy

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


## event counts ----
## create dataset that reports event counts for each outcome of interest

data_event_counts <-
  bind_rows(
    metaparams |> 
      distinct(cohort, subgroup, outcome, .keep_all=TRUE) |>
      mutate(
        method="unadjusted",
        spec=""
      ),
    metaparams
  ) |>
  group_by(cohort, method, spec, subgroup, outcome) |>
  mutate(
    data = pmap(
      list(cohort, method, spec, subgroup, outcome), 
      function(cohort, method, spec, subgroup, outcome) {
        
        data_all %>%
          mutate(
            
            subgroup_level = .[[subgroup]],
            
            wt = ifelse(
              method!="unadjusted",
              .[[paste0("wt_", cohort, "_", method, "_", spec)]],
              1
            ),
            treatment_date = vax_date-1L, 
            
            event_date = as.Date(.[[paste0(outcome, "_date")]]),
            
            # person-time is up to and including censor date
            censor_date = pmin(
              dereg_date,
              death_date,
              study_dates$followupend_date,
              treatment_date + maxfup,
              na.rm=TRUE
            ),
            
            noncompetingcensor_date = pmin(
              dereg_date,
              study_dates$followupend_date,
              treatment_date + maxfup,
              na.rm=TRUE
            ),
            
            event_time = tte(treatment_date, event_date, censor_date, na.censor=FALSE),
            event_indicator = censor_indicator(event_date, censor_date),
          ) |>
          group_by(subgroup_level, treatment) |>
          summarise(
            n = roundmid_any(sum(wt), sdc.limit),
            persontime = sum(wt*as.numeric(censor_date - (vax_date - 1))),
            count = sum(wt*event_indicator)
          )
        }
      )
    ) |>
  unnest(data)

write_feather(data_event_counts, fs::path(output_dir, "table_event_counts.arrow"))


