
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
  spec <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  spec <- args[[2]]
}


## create output directories ----

output_dir <- here_glue("output", "3-adjust", cohort, "weight-{spec}")
fs::dir_create(output_dir)

# Import and prepare data ----

## one pow per patient ----
data_cohort <- read_feather(here_glue("output", "2-prepare", cohort, "data_cohort.arrow"))

print_data_size(data_cohort)

## select variables used for weighting
data_preweight <-
  data_cohort |>
  select(
    patient_id,
    vax_product,
    treatment,
    vax_date,
    all_of(weighting_variables[[spec]]),
  ) |>
  arrange(patient_id)

# calculate balancing weights using the weightit function
obj_weightit <- 
  weightit(
    formula = formula(paste0("treatment ~ ", weighting_formulae[[spec]])),
    data = data_preweight,
    method = "glm", 
    estimand = "ATE",
    stabilize = TRUE
  )

data_weights <- 
  tibble(
    patient_id = data_preweight$patient_id,
    treatment = obj_weightit$treat,
    ps = obj_weightit$ps,
    weight = obj_weightit$weights, # weight = get_w_from_ps(ps=ps, treat=treatment,  estimand = "ATE")
  ) 

## weights and PS relationship:
# weight = (treatment/ps) + ((1-treatment)/(1-ps)),
# ps = (treatment/weight) + ((1-treatment)*(1-(1/weight)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Summarise weighted population and export ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data_weights <-
  data_weights |>
  left_join(data_cohort |> select(patient_id, vax_date), by="patient_id") 

write_feather(data_weights, fs::path(output_dir, "data_adjusted.arrow"))

summary(obj_weightit)




