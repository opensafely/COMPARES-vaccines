# # # # # # # # # # # # # # # # # # # # #
# Purpose: import processed data and filter out people who are excluded from the main analysis
# outputs:
#  - inclusion/exclusions flowchart data (up to matching step)
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('gt')
library('gtsummary')


## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

## Import redaction functions
source(here("analysis", "0-lib", "redaction.R"))


## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "age75plus" #currently `age75plus` or `cv`
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
}

# derive subgroup info
cohort_sym <- sym(cohort)

## create output directories for data ----
output_dir <- here("output", "2-select", cohort)
fs::dir_create(output_dir)


# Main ----

## Import processed data ----

data_prepared <- read_feather(here("output", "2-select", "data_prepared.arrow"))


# Report total number of people vaccinated in time period, by vaccine type ----
# This output ignores cohort, so is the same across different cohorts
# but it is more lightweight than outputting in the `data_process` script, because then
# the release action will need to copy in the entire processed data set

## unrounded totals
total_n_unrounded <-
  bind_rows(
    tibble(vax_product="any", n=nrow(data_prepared)),
    count(data_prepared |> mutate(boost_type=fct_other(vax_product, keep=treatment_lookup$treatment, other_level="other")), vax_product, .drop=FALSE)
  ) |>
  mutate(
    pct = n/first(n)
  )
write_csv(total_n_unrounded, fs::path(output_dir, "total_allcohorts_unrounded.csv"))

## rounded totals
total_n_rounded <-
  total_n_unrounded |>
  mutate(
    n= ceiling_any(n, sdc.limit),
    pct = n/first(n)
  )
write_csv(total_n_rounded, fs::path(output_dir, "total_allcohorts_rounded.csv"))





## Define selection criteria ----
data_criteria <-
  data_prepared |>
  filter(!!cohort_sym) |>
  transmute(
    patient_id,
    vax_date,
    vax_product,
    has_age = !is.na(age_eligible),
    has_sex = !is.na(sex) & !(sex %in% c("intersex", "unknown")),
    has_imd = !is.na(imd_Q5),
    has_ethnicity5 = !is.na(ethnicity5),
    has_region = !is.na(region),
    #has_msoa = !is.na(msoa),
    #isnot_hscworker = !hscworker,
    #isnot_carehomeresident = !care_home_combined,
    #isnot_endoflife = !endoflife,
    #isnot_housebound = !housebound,
    #no_prior_productA = !vaxhist_productA,
    #no_prior_productB = !vaxhist_productB,
    prior_vax_interval_atleast12weeks,
    vax_product_of_interest = vax_product %in% c(productA, productB),
    prior_vax_1plus = (prior_vax_count >= 1),
    #has_norecentcovid = ((vax_date - anycovid_0_date) >= 28) | is.na(anycovid_0_date),
    isnot_inhospital = !inhospital,
    
    include = (
      prior_vax_interval_atleast12weeks & vax_product_of_interest & #no_prior_productA & no_prior_productB &
        prior_vax_1plus &
        has_age & has_sex & has_imd & has_region & #has_ethnicity &
        #isnot_hscworker &
        #isnot_endoflife &
        #has_norecentcovid &
        isnot_inhospital
    ),
  )

data_cohort <- 
  data_criteria |>
  filter(include) |>
  select(patient_id) |>
  left_join(data_prepared, by="patient_id") |>
  droplevels()

write_feather(data_cohort, fs::path(output_dir, "data_cohort.arrow"))


## report unadjusted balance for baseline variables between levels  of treatment ---- 
# This output includes SDC

table_cohort <- 
  data_cohort |>
  mutate(
    N = 1L,
    weight = 1L
  ) %>%
  table1_summary_smd(
    treatment = treatment,
    weight = weight, 
    label = variable_labels,
    threshold = sdc.limit
  )

write_csv(table_cohort, fs::path(output_dir, "table_cohort.csv"))


## output simple dataset containing exclusions criteria met ----

## TODO: check if outputting this is necessary

data_inclusioncriteria <- data_criteria |>
  transmute(
    patient_id,
    vax_product,
    c0 = TRUE,
    c1 = c0 & vax_product_of_interest,
    c2 = c1 & prior_vax_interval_atleast12weeks, #& no_prior_productA & no_prior_productB,
    c3 = c2 & prior_vax_1plus,
    c4_1 = c3 & (has_age & has_sex & has_imd & has_region),
    # c4_3 = c3 & (isnot_endoflife),
    # c4_4 = c3 & (has_norecentcovid),
    # c4_5 = c3 & (isnot_inhospital),
    #c4 = c4_1 & c4_2 & c4_3 & c4_4 & c4_5
  ) 

# remove large in-memory objects
remove(data_criteria)

write_feather(data_inclusioncriteria, sink = fs::path(output_dir, "data_inclusioncriteria.arrow"))


## Create flowchart ----

create_flowchart <- function(round_level = 1){
  
  flowchart_output <-
    data_inclusioncriteria |>
    select(-patient_id, -c0) |>
    filter(c1) |>
    group_by(vax_product) |>
    summarise(
      across(.cols=everything(), .fns=~ceiling_any(sum(.), round_level))
    ) |>
    pivot_longer(
      cols=-c(vax_product),
      names_to="criteria",
      values_to="n"
    ) |>
    mutate(
      level = if_else(str_detect(criteria, "c\\d+$"), 1, 2),
      n_level1 = if_else(level==1, n, NA_real_),
      n_level1_fill = n_level1
    ) |>
    fill(n_level1_fill) |>
    group_by(vax_product) |>
    mutate(
      n_exclude = lag(n_level1_fill) - n,
      pct_exclude = n_exclude/lag(n_level1_fill),
      pct_all = n_level1 / first(n),
      pct_step = n_level1 / lag(n_level1_fill),
      #crit = str_extract(criteria, "^c\\d+"),
      crit = criteria,
      criteria = fct_case_when(
        crit == "c1" ~ "Received COVID-19 vaccine between X and X",
        crit == "c2" ~ "  with no prior Covid-19 vaccine within 12 weeks",
        crit == "c3" ~ "  with at least 1 prior COVID-19 vaccine dose",
        crit == "c4_1" ~ "    no missing demographic information",
        crit == "c4_2" ~ "    not a health and social care worker",
        crit == "c4_3" ~ "    not end-of-life",
        crit == "c4_4" ~ "    no documented COVID-19 infection/disease within prior 28 days",
        crit == "c4_5" ~ "    not admitted in hospital at time of booster",
        crit == "c4" ~ "  included in matching run",
        TRUE ~ "NA_character_boop" # should not appear
      )
    )
  
  return(flowchart_output)
}

## unrounded flowchart
data_flowchart <- create_flowchart(1)
write_feather(data_flowchart, fs::path(output_dir, "flowchart.arrow"))
#write_csv(data_flowchart, here("output", "data", "flowchart.csv"))

## rounded flowchart
data_flowchart_rounded <- create_flowchart(sdc.limit)
write_feather(data_flowchart_rounded, fs::path(output_dir, "flowchart_rounded.arrow"))
write_csv(data_flowchart_rounded, fs::path(output_dir, "flowchart_rounded.csv"))

## unrounded totals
total_n_unrounded <-
  bind_rows(
    tibble(vax_product="any", n=nrow(data_inclusioncriteria)),
    count(data_inclusioncriteria |> mutate(boost_type=fct_other(vax_product, keep=treatment_lookup$treatment, other_level="other")), vax_product, .drop=FALSE)
  ) |>
  mutate(
    pct = n/first(n)
  )
write_csv(total_n_unrounded, fs::path(output_dir, "total_unrounded.csv"))

## rounded totals
total_n_rounded <-
  total_n_unrounded |>
  mutate(
    n= ceiling_any(n, sdc.limit),
    pct = n/first(n)
  )
write_csv(total_n_rounded, fs::path(output_dir, "total_rounded.csv"))

## remove large in-memory objects
remove(data_inclusioncriteria)
