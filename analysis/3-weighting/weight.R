
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: 
# use IPW to create a weighted pseudo population such that covariates are balanced across treatment groups
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('survival')
library('MatchIt')
library("WeightIt")
library("cobalt")
library("doParallel")


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
  weightset <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  weightset <- args[[2]]
}


## create output directories ----

output_dir <- here_glue("output", "3-cohorts", cohort, "weight{weightset}")
fs::dir_create(output_dir)

# Import and prepare data ----

## one pow per patient ----
data_cohort <- read_feather(here_glue("output", "3-cohorts", cohort, "data_cohort.arrow"))

print_data_size(data_cohort)

## select variables used for weighting
data_preweight <-
  data_cohort |>
  select(
    patient_id,
    vax_product,
    treatment,
    vax_date,
    all_of(weighting_variables[[weightset]]),
  ) |>
  arrange(patient_id)

# calculate balancing weights using the weightit function
obj_weightit <- 
  weightit(
    formula = formula(paste0("treatment ~ ", weighting_formulae[[weightset]])),
    data = data_preweight,
    method = "glm", 
    estimand = "ATE"
  )

data_weights <- 
  tibble(
    patient_id = data_preweight$patient_id,
    treatment = obj_weightit$treat,
    weight = obj_weightit$weights,
    ps = obj_weightit$ps
  ) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Summarise weighted population and export ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data_weights <-
  data_weights |>
  left_join(data_cohort |> select(patient_id, vax_date), by="patient_id") 

write_feather(data_weights, fs::path(output_dir, "data_weights.arrow"))

summary(obj_weightit)


