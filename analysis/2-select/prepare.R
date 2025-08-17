######################################
# Purpose: 
# process (tidy, clean, reshape, derive) data extracted using ehrQL (or dummy data), including:
# - standardising variables (eg convert to factor) and derives some new ones
# - outputting skim files which summarise data types and distributions
# Outputs:
# - a prepared dataset
# - skim files
######################################

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('lubridate')
library('arrow')
library('here')
library("dtplyr")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

# create output directory
output_dir <- here("output", "2-select")
fs::dir_create(output_dir)

# Import and process data ----

## Import extracted dataset ----
data_extract <- read_feather(here("output", "1-extract", "extract.arrow"))

# print details about extracted dataset
capture.output(
  skimr::skim_without_charts(data_extract),
  file = fs::path(output_dir, "data_extract_skim.txt"),
  split = FALSE
)


stopifnot(
  "inconsistency between ethnicity5 and ethnicity 16" = identical(data_extract$ethnicity5, ethnicity_16_to_5(data_extract$ethnicity16))
)

## Process extracted dataset ----
data_prepared <- 
  data_extract |>
  mutate(

    # use short product names
    vax_product = factor(vax_product, vax_product_lookup, names(vax_product_lookup)) |> fct_na_value_to_level("other"),

    # binary variable for the exposure
    # helpful for various downstream matching / plotting / table functions
    treatment = case_when(
      vax_product==productA ~ 0L,
      vax_product==productB ~ 1L,
      TRUE ~ NA_integer_
    ),

    vax_date = vax_date - 1L,
    
    # vaccination date represented as an integer, using for matching instead of date-formatted variable to avoid issues
    vax_day = as.integer(vax_date - study_dates$studystart_date),
    
    # all subgroup dummy variable
    all = factor("all"),

    ageband = cut(
      age_eligible, # use fixed date to ascertain age so that age bands align with eligibility. because age and vax date are closely matched, this doesn't cause any problems
      breaks=c(-Inf, 18, 50, 65, 75, 80, 85, Inf),
      labels=c("under 18", "18-49", "50-64", "65-74", "75-79", "80-84", "85+"),
      right=FALSE
    ),

    ethnicity5 = factor(ethnicity5, levels = factor_levels$ethnicity5) |> fct_na_value_to_level(level = "Unknown"),
    ethnicity16 = factor(ethnicity16, levels = factor_levels$ethnicity16) |> fct_relabel(~ str_extract(.x, "(?<= - )(.*)")), # pick up everything after " - "
  
    region = fct_collapse(
      region,
      `East of England` = "East",
      `London` = "London",
      `Midlands` = c("West Midlands", "East Midlands"),
      `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
      `North West` = "North West",
      `South East` = "South East",
      `South West` = "South West"
    ),

    imd_Q5 = cut(
      imd,
      breaks = (32844/5)*c(-0.1,1,2,3,4,5),
      labels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
      include.lowest = TRUE,
      right = FALSE
    ),

    carehome = care_home_tpp | care_home_code, # any carehome flag

      # clinically at-risk group
    cv = immunosuppressed | ckd | crd | diabetes | cld |
      cns | chd | learndis | smi,

    multimorb =
      (severe_obesity) +
      (chd) +
      (ckd)+
      (diabetes) +
      (cld)+
      (crd)+
      (cns)+
      #(cancer)+
      #(learndis)+
      #(sev_mental),
      0,
    multimorb = cut(multimorb, breaks = c(0, 1, 2, Inf), labels=c("0", "1", "2+"), right=FALSE),

    ## process partial covid vaccine history
    
    prior_vax_count = vax_covid_prior_count,
    prior_vax_interval = as.integer(vax_date - vax_covid_prior_1_date),
    prior_vax_interval_atleast12weeks = prior_vax_interval >= 12*7,
    vax_count_realistic = prior_vax_count <= max_prior_vax_count, # earlier than 2020-12-08 to include trial participants, but excludes "unknowns" coded as eg 1900-01-01
    prior_vax_interval_bigM = if_else(prior_vax_count==0, 365L*5L, as.integer(prior_vax_interval)), #if current vaccine is first recorded vaccine, then set vax_interval to be very large
    prior_vax_count_group = cut(
      prior_vax_count,
      breaks = c(0,    1,     3,     6,  Inf),
      labels = c("0", "1-2", "3-5", "6+"),
      include.lowest = TRUE,
      right = FALSE
    ),
    
    
    ## process outcomes data
    censor_date = pmin(dereg_date, study_dates$followupend_date, na.rm=TRUE),
    noncovid_death_date = if_else(!is.na(death_date) & is.na(covid_death_date), death_date, as.Date(NA_character_)),

    # earliest covid event after study start
    any_covid_date = pmin(covid_emergency_date, covid_admitted_date, covid_death_date, na.rm=TRUE),
    
    # KEEP THIS as a reminder to replace event dates with more severe dates if they precede less severe dates
    # for use if we decide to include source specific endpoints
    #covid_emergency_date = pmin(covid_emergency_date, covid_admitted_date, covid_critcare_date, covid_death_date, na.rm=TRUE),
    #covid_admitted_date = pmin(covid_admitted_date, covid_critcare_date, covid_death_date, na.rm=TRUE),
    #covid_critcare_date = pmin(covid_critcare_date, covid_death_date, na.rm=TRUE),
    

    # define cohorts
    age75plus = age_eligible >= 75,
    #cv = cv,
    is_eligible = age_eligible | cv,

  ) %>%
  select(
    -vax_covid_prior_1_date, -vax_covid_prior_1_product,
    -vax_covid_prior_2_date, -vax_covid_prior_2_product,
    -vax_covid_prior_3_date, -vax_covid_prior_3_product
  )

# print details about prepared dataset
capture.output(
  skimr::skim_without_charts(data_prepared),
  file = fs::path(output_dir, "data_prepared_skim.txt"),
  split = FALSE
)

write_feather(data_prepared, sink = fs::path(output_dir, "data_prepared.arrow"))


