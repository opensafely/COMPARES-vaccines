
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: match productA recipients to productB recipients
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
  matchset <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  matchset <- args[[2]]
}


## create output directories ----

output_dir <- here_glue("output", "3-cohorts", cohort, "match{matchset}")
fs::dir_create(output_dir)

# Import and prepare data ----

## one pow per patient ----
data_cohort <- read_feather(here("output", "3-cohorts", cohort, "data_cohort.arrow"))

print_data_size(data_cohort)


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
    all_of(matching_variables[[matchset]]$all),
  ) |>
  arrange(patient_id)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Multi-threaded version ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# 
# # Set up parallelisation and run matching ----
# 
# # create function that catches errors in case no matches are found within a thread
# safely_matchit <- purrr::safely(matchit)
# 
# ## parallelisation preliminaries ----
# 
# parallel::detectCores() # how many cores available?
# n_threads <- max(2, parallel::detectCores())
# 
# cluster <- parallel::makeCluster(
#   n_threads,
#   type = "PSOCK" # this should work across multi-core windows or linux machines
# )
# print(cluster)
# 
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = cluster)
# 
# # create parallel matching streams
# matchthreads <- as.character(unique(sort(data_prematch$thread_variable)))
# 
# table(data_prematch$thread_variable, useNA="ifany")
# 
# ## match in parallel ----
# data_matches <-
#   foreach(
#     matchthread = matchthreads,
#     .combine = 'bind_rows',
#     .packages = c("dplyr", "MatchIt", "tibble", "lubridate")
#   ) %dopar% {
#   #for(matchthread in matchthreads){
#   #matchthread<-"85+"
# 
#     data_thread <- data_prematch |> filter(thread_variable==matchthread)
# 
#     # run matching algorithm
#     obj_matchit <-
#       safely_matchit(
#       #matchit(
#         formula = treatment ~ 1,
#         data = data_thread,
#         method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
#         replace = FALSE,
#         estimand = "ATT",
#         exact = matching_variables[[matchset]]$exact,
#         caliper = matching_variables[[matchset]]$caliper, std.caliper=FALSE,
#         m.order = "data", # data is sorted on (effectively random) patient ID
#         #verbose = TRUE,
#         ratio = 1L # could also consider exact matching only, with n:m ratio, determined by availability
#       )[[1]]
# 
#     ## process matchit object to give one row per candidate, matched status (0/1) and match id
# 
#     data_matches <-
#       if(is.null(obj_matchit)){
#         tibble(
#           patient_id = data_thread$patient_id,
#           matched = FALSE,
#           thread_id = data_thread$thread_id,
#           threadmatch_id = NA_integer_,
#           treatment = data_thread$treatment,
#           weight = 0,
#         )
#       } else {
#         tibble(
#           patient_id = data_thread$patient_id,
#           matched = !is.na(obj_matchit$subclass),
#           thread_id = data_thread$thread_id,
#           threadmatch_id = as.integer(as.character(obj_matchit$subclass)),
#           treatment = obj_matchit$treat,
#           weight = obj_matchit$weights,
#         )
#       }
# 
#   data_matches
# }
# 
# parallel::stopCluster(cl = cluster)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## single-threaded version ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

obj_matchit <- 
  matchit(
    formula = treatment ~ 1,
    data = data_prematch,
    method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
    replace = FALSE,
    estimand = "ATT", # since we are doing exact matching, ATT is equivalent to ATO (average treatment in the overlap)
    exact = matching_variables[[matchset]]$exact,
    caliper = matching_variables[[matchset]]$caliper, std.caliper=FALSE,
    m.order = "data", # data is sorted on (effectively random) patient ID
    #verbose = TRUE,
    ratio = 1L # could also consider exact matching only, with n:m ratio, determined by availability
  )


data_matches <- 
  tibble(
    patient_id = data_prematch$patient_id,
    matched = !is.na(obj_matchit$subclass),
    thread_id = 1L,
    threadmatch_id = obj_matchit$subclass,
    treatment = obj_matchit$treat,
    weight = obj_matchit$weights,
  ) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Summarise matched / unmatched patients and export ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data_matches <-
  data_matches |>
  left_join(data_cohort |> select(patient_id, vax_date), by="patient_id") |>
  arrange(thread_id, threadmatch_id) |>
  mutate(
    match_id = dense_rank(threadmatch_id * max(thread_id) + thread_id) # create unique match id across all threads
  )

write_feather(data_matches, fs::path(output_dir, "data_matches.arrow"))

data_matches |>
  group_by(treatment, matched) |>
  summarise(
    n=n()
  )



























# bootstrap sampling ----

## bootstrap sample matched pairs and use this sampling throughout the analysis
## doing it here avoids repeating the sampling process in each individual outcome script
## and provides consistency across different analyses
## but the leg work is still done by the analysis scripts

boot_n <- 500 # more than necessary, can select fewer in the analysis scripts

boot_id <- seq_len(boot_n)

match_ids <- unique(data_matches$match_id[!is.na(data_matches$match_id)])

set.seed(20230401)

boot_samples <-
  tibble(boot_id) |>
  mutate(
    match_id = map(boot_id, ~sample(match_ids, size=length(match_ids), replace=TRUE))
  ) |>
  unnest(match_id)

write_feather(boot_samples, fs::path(output_dir, "boot_samples.arrow"))


