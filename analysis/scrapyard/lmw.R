
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: 
# use lmw package (NOT CURENTLY IN OpenSAFELY R image) to obtain the implied weights if doing outcome regression to estimate the ATE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('survival')
library("cobalt")
library("doParallel")
library("lmw")


## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))



## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "age75plus" #currently `age75plus` or `cv`
  spec <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  spec <- args[[2]]
}


## create output directories ----

output_dir <- here_glue("output", "3-adjust", cohort, "lmw-{spec}")
fs::dir_create(output_dir)

# Import and prepare data ----

## one pow per patient ----
data_cohort <- read_feather(here_glue("output", "2-select", cohort, "data_cohort.arrow"))

print_data_size(data_cohort)

# for now, use same variables as weighting model
lmw_variables <- weighting_variables
lmw_formulae <- weighting_formulae

## select variables used for outcome model
data_preweight <-
  data_cohort |>
  select(
    patient_id,
    vax_product,
    treatment,
    vax_date,
    all_of(lmw_variables[[spec]]),
  ) |>
  arrange(patient_id)

# calculate balancing weights using the lmw function

# need to check order here!!
obj_lmw <- 
  lmw(
    formula = formula(paste0("~ ", "treatment + ", lmw_formulae[[spec]])),
    data = data_preweight,
    treat = "treatment",
    estimand = "ATE",
    method = "MRI", # TODO: or use the other one that splits by treatment??
  )

data_weights <- 
  tibble(
    patient_id = data_preweight$patient_id,
    treatment = data_preweight$treatment, # treatment = obj_lmw$treat <--- this is a factor, not a binary so doesn't work as nicely as taking from original dataset
    weight = obj_lmw$weights, 
    ps = (treatment/weight) + ((1-treatment)*(1-(1/weight)))
  )

## weights and PS relationship:
# weight = (treatment/ps) + ((1-treatment)/(1-ps)), or use WeightIt::get_w_from_ps(ps=ps, treat=treatment,  estimand = "ATE")
# ps = (treatment/weight) + ((1-treatment)*(1-(1/weight)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Summarise weighted population and export ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data_weights <-
  data_weights |>
  left_join(data_cohort |> select(patient_id, vax_date), by="patient_id") 

write_feather(data_weights, fs::path(output_dir, "data_adjusted.arrow"))




