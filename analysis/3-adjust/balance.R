
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: 
# use one of 3 methods to balance pre-treatment variables across treatment groups
#
#
# 1. "match" ---
# match productA recipients to productB recipients using MatchIt package
#
#
# 2. "weight" ---
# use IPW to create a weighted pseudo population such that covariates are balanced across treatment groups
#
#
# 3. "lmw" ---
# use the weights that are implied if we had used linear regression to balance the outcome and estimate the effect of a binary treatment
# these weight can be obtained without using the outcome variable!
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('survival')
library("cobalt")

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
  method <- "lmw"
  spec <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  method <- args[[2]]
  spec <- args[[3]]
}


## create output directories ----

output_dir <- here_glue("output", "3-adjust", cohort, "{method}-{spec}")
fs::dir_create(output_dir)

# Import and prepare data ----

## one pow per patient ----
data_cohort <- read_feather(here_glue("output", "2-select", cohort, "data_cohort.arrow"))

print_data_size(data_cohort)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Balance using method = match ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

if(method == "match"){
  
  library("MatchIt")
  
  ## select all matching candidates and variables necessary for matching
  data_prematch <-
    data_cohort |>
    mutate(
      # create variable to parallelise on
      thread_variable = ageband,
      thread_id = dense_rank(thread_variable)
    ) |>
    select(
      thread_id,
      thread_variable,
      patient_id,
      vax_product,
      treatment,
      vax_date,
      all_of(matching_variables[[spec]]$all),
    ) |>
    arrange(patient_id)
  
  ## match using single-threaded version 
  ## (see scrapyard directory for unused version using parallelisation)
  
  obj_matchit <- 
    matchit(
      formula = treatment ~ 1,
      data = data_prematch,
      method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
      replace = FALSE,
      estimand = "ATT", # since we are doing exact matching, ATT is equivalent to ATU. although we'll actually get the ATO (average treatment in the overlap)
      exact = matching_variables[[spec]]$exact,
      caliper = matching_variables[[spec]]$caliper, std.caliper=FALSE,
      m.order = "data", # data is sorted on (effectively random) patient ID
      #verbose = TRUE,
      ratio = 1L # could also consider exact matching only, with n:m ratio, determined by availability
    )
  
  
  data_matches <- 
    tibble(
      patient_id = data_prematch$patient_id,
      matched = !is.na(obj_matchit$subclass),
      thread_id = 1L,
      threadmatch_id = as.integer(obj_matchit$subclass),
      treatment = obj_matchit$treat,
      weight = obj_matchit$weights,
      ps = (treatment/weight) + ((1-treatment)*(1-(1/weight)))
    ) 
  
  
  data_weights <-
    data_matches |>
    arrange(thread_id, threadmatch_id) |>
    mutate(
      match_id = dense_rank(threadmatch_id * max(thread_id) + (thread_id - 1L))  # create unique match id across all threads
    )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Balance using method = weight ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

if(method == "weight"){

  
  library("WeightIt")
  
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
      weight = obj_weightit$weights, # weight = get_w_from_ps(ps=ps, treat=treatment, estimand = "ATE")
    ) 
  
  ## weights and PS relationship:
  # weight = (treatment/ps) + ((1-treatment)/(1-ps)),
  # ps = (treatment/weight) + ((1-treatment)*(1-(1/weight)))

}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Balance using method = lmw ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

if(method == "lmw"){
  
  library("lmw")
  
  ## select variables used for outcome models
  data_prelmw <-
    data_cohort |>
    select(
      patient_id,
      vax_product,
      treatment,
      vax_date,
      all_of(lmw_variables[[spec]]),
    ) |>
    arrange(patient_id)
  
  # calculate balancing weights using the lmw function from lmw package
  obj_lmw <- 
    lmw(
      formula = formula(paste0("~ ", "treatment + ", lmw_formulae[[spec]])),
      data = data_prelmw,
      treat = "treatment",
      estimand = "ATE",
      method = "MRI", # MRI gets weights as if there were a separate outcome model for each treatment group
    )
  
  data_weights <- 
    tibble(
      patient_id = data_prelmw$patient_id,
      treatment = data_prelmw$treatment, # treatment = obj_lmw$treat <--- this is a factor, not a binary so doesn't work as nicely as taking from original dataset
      weight = obj_lmw$weights, 
      ps = (treatment/weight) + ((1-treatment)*(1-(1/weight)))
    )
  
}
  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Export ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

write_feather(data_weights, fs::path(output_dir, "data_adjusted.arrow"))

