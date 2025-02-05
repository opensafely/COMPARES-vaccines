# _________________________________________________
# Purpose:
# creates metadata objects for aspects of the study design
# should be sourced (ie `source(".../design.R")`) in the analysis scripts
# defines some look up tables to use in the codebase
# this script should be sourced (using `source(here("analysis", "0-lib", "design.R"))`) at the start of each R script
# _________________________________________________

## TODO: LINK TO PROTOCOL


# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')

## create output directories ----
fs::dir_create(here("analysis", "0-lib"))


# key study dates
# The dates are defined in json format so they can be read in by R and python scripts
# json has no easy way to comment, so explanation for dates is here:
# - firstpossibleax_date is the date from which we want to identify covid vaccines. the mass vax programme was 8 Dec 2020 but other people were vaccinated earlier in trials, so want to pick these people up too (and possibly exclude them)
# - start_date is when we start the observation period, typically at the start of the vaccination campaign (we may restrict the analysis to those vaccinated at a later date if only one product wasavailable at the start)
# - end_date is when we stop the observation period. Most likely the day before the start of the next campaign.

study_dates <-
  jsonlite::read_json(path = here("analysis", "0-lib", "study-dates.json")) %>%
  map(as.Date)

# make these available in the global environment
# so we don't have to use `study_dates$start_date` or `start_date <- study_dates$start_date` in each script
list2env(study_dates, globalenv())


# define cohorts ----

# TODO: define cohorts here, eg age_eligible and clinically vulnerable specific for each campaign


# define outcomes ----

# maximum follow-up duration, in days
# currently up to 24 weeks after time-zero

maxfup <- 7L*24L

events_lookup <- tribble(
  ~event, ~event_var, ~event_descr,

  # other
  "test", "covid_test_date", "SARS-CoV-2 test", 

  # effectiveness
  "covid_emergency", "covid_emergency_date", "COVID-19 A&E attendance", 
  "covid_admitted", "covid_admitted_date", "COVID-19 hospitalisation", 
  "covid_critcare", "covid_critcare_date", "COVID-19 critical care", 
  "covid_death", "covid_death_date", "COVID-19 death", 
  "death", "death_date", "Any death", 

  # safety
  "admitted", "admitted_date", "Unplanned hospitalisation", 
  "emergency", "emergency_date", "A&E attendance", 
  "pericarditis", "pericarditis_date", "Pericarditis", 
  "myocarditis", "myocarditis_date", "Myocarditis", 

  # negative control
  "noncovid_death", "noncovid_death_date", "Non-COVID-19 death", 
  "fracture", "fracture_date", "Fracture", 
)


# maximum allowable prior count of covid vaccines (too high indicates unreliable data, but should increase each campaigns)
# something like campaign number + 2 (to allow 2 vaccines in first campaign and maybe sneaky extra one somewhere)
# TODO: figure out if certain populations had more than this routinely
max_prior_vax_count <- 10L

# output from https://jobs.opensafely.org/opensafely-internal/tpp-vaccination-names/ workspace
# shows all possible covid vaccination names in TPP

productA <- "pfizer"
productB <- "moderna"
#productC <- "vidprevtyn"

# lookup to rename TPP product names to coding-friendly porduct names
vax_product_lookup <- c(
  "pfizer" = "COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
  "az" = "COVID-19 Vaccine Vaxzevria 0.5ml inj multidose vials (AstraZeneca)",
  "moderna" = "COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
  "pfizerBA1" = "Comirnaty Original/Omicron BA.1 COVID-19 Vacc md vials",
  "pfizerBA45" = "Comirnaty Original/Omicron BA.4-5 COVID-19 Vacc md vials",
  "pfizerXBB15" = "Comirnaty Omicron XBB.1.5 COVID-19 Vacc md vials",
  "vidprevtyn" = "COVID-19 Vacc VidPrevtyn (B.1.351) 0.5ml inj multidose vials",
  "modernaomicron" = "COVID-19 Vac Spikevax (Zero)/(Omicron) inj md vials",
  "pfizerchildren" = "COVID-19 mRNA Vaccine Comirnaty Children 5-11yrs 10mcg/0.2ml dose conc for disp for inj MDV (Pfizer)",
  "azhalf" = "COVID-19 Vac AZD2816 (ChAdOx1 nCOV-19) 3.5x10*9 viral part/0.5ml dose sol for inj MDV (AstraZeneca)",
  "modernaXBB15" = "COVID-19 Vacc Spikevax (XBB.1.5) 0.1mg/1ml inj md vials",
  "novavax" = "COVID-19 Vac Nuvaxovid (recombinant, adj) 5micrograms/0.5ml dose susp for inj MDV (Novavax CZ a.s.)"
)


# lookup to rename coding-friendly product names to publication-friendly product names
vax_shortname_lookup <- c(
  "BNT162b2" = "pfizer",
  "ChAdOx1" = "az",
  "mRNA-1273" = "moderna",
  "BNT162b2/BA.1" = "pfizerBA1",
  "BNT162b2/BA.4-5" = "pfizerBA45",
  "BNT162b2/XBB.1.5" = "pfizerXBB15",
  "VidPrevtyn" = "vidprevtyn",
  "mRNA-1273/omicron" = "modernaomicron",
  "BNT162b2/children" = "pfizerchildren",
  "ChAdOx1/2" = "azhalf",
  "mRNA-1273/XBB.1.5" = "modernaXBB15",
  "Novavax" = "novavax",
  "Other" = "other"
)


treatment_lookup <-
  tribble(
    ~treatment, ~treatment_descr,
    "pfizerBA45", "pfizer/BA.4-5",
    "sanofi", "Sanofi",
  )

# where to split follow-up time after recruitment
postbaselinecuts <- as.integer(c(0,7,14,28,56,84,112,140,168))

# redaction threshold
sdc.limit <- 6L

## lookups to convert coded variables to full, descriptive variables ----


recoder <-
  lst(
    cohort = c(
      `Clinically vulnerable` = "cv",
      `Aged 75 years or over` = "age75plus"
    ),
    subgroups = c(
      `Main` = "all",
      `Age` = "ageband",
      `Clinically at-risk` = "cv",
      `Prior COVID-19 vaccine count` = "vax_previous_group",
      #`Prior SARS-CoV-2 infection status` = "prior_covid_infection"
      NULL
    ),
    status = c(
      `Unmatched`= "unmatched",
      `Matched` = "matched"
    ),
    treatment = c(
      `pfizer/BA.4-5` = "0",
      `Sanofi` = "1"
    ),
    # outcome = c(
    #   "Positive SARS-CoV-2 test"= "postest",
    #   "COVID-19 A&E attendance" = "covidemergency",
    #   "COVID-19 hospitalisation" = "covidadmitted",
    #   "COVID-19 critical care" = "covidcritcare",
    #   "COVID-19 death" = "coviddeath",
    #   "Non-COVID-19 death" = "noncoviddeath",
    #   "All-cause death" = "death",
    #   "Fracture" = "fracture",
    #   "Pericarditis" = "pericarditis",
    #   "Myocarditis" = "myocarditis",
    #   NULL
    # ),
    outcome = set_names(events_lookup$event, events_lookup$event_descr),
    all = c(` ` = "all"),
    ageband = c(
      "50-64", "65-74", "75-79", "80-84", "85+"
    ) %>% {set_names(.,.)},
    cv = c(
      `Clinically at-risk` = "TRUE",
      `Not clinically at-risk` = "FALSE"
    ),
    prior_covid_infection = c(
      `No prior SARS-CoV-2 infection` = "FALSE",
      `Prior SARS-CoV-2 infection` = "TRUE"
    ),
    vax_previous_group = c(
      "0-1", "2-4", "5", "6+"
    ) %>% {set_names(.,.)},
  )

## variable labels in tables and plots ----

variable_labels <- 
  list(
    N  = "Total N",
    treatment_descr = "Vaccine product",
    prior_vax_interval = "Days since previous vaccine",
    prior_vax_count_group = "Previous vaccine count",
    age_eligible = "Age",
    ageband = "Age band",
    sex = "Sex",
    ethnicity5 = "Ethnicity",
    imd_Q5 = "Deprivation",
    region = "Region",
    cv = "Clinically at-risk",
    
    #housebound = "Clinically housebound",
    #care_home_combined = "Care/nursing home resident",
    
    severe_obesity = "Body Mass Index > 40 kg/m^2",
    
    chd = "Chronic heart disease",
    ckd = "Chronic kidney disease",
    diabetes = "Diabetes",
    cld = "Chronic liver disease",
    crd = "Chronic respiratory disease",
    #asthma = "Asthma",
    cns = "Chronic neurological disease",
    
    immunosuppressed = "Immunosuppressed",
    #immuno_any = "Immunosuppressed (all)",
    
    # immdx = "Immunocompromising diagnosis",
    # immrx = "Immunosuppressive medications, previous 3 years",
    # dxt_chemo = "Chemotherapy, previous 3 years",
    # cancer = "Cancer, previous 3 years",
    # asplenia = "Asplenia or poor spleen function",
    # solid_organ_transplant = "Solid organ transplant",
    # hiv_aids = "HIV/AIDS",
    
    multimorb = "Morbidity count",
    
    learndis = "Learning disabilities",
    smi = "Serious mental illness"#,
    
    # prior_tests_cat = "Number of SARS-CoV-2 tests",
    # 
    # prior_covid_infection = "Prior documented SARS-CoV-2 infection",
    # 
    # vaxhist_pfizer  = "Previously received Pfizer (original)",
    # vaxhist_az  = "Previously received AZ",
    # vaxhist_moderna  = "Previously received Moderna",
    # vaxhist_pfizerBA1  = "Previously received Pfizer/BA.1",
    # vaxhist_pfizerXBB15  = "Previously received Pfizer/XBB.1.5",
    # vaxhist_modernaomicron  = "Previously received Moderna/Omicron",
    # vaxhist_modernaXBB15  = "Previously received Moderna/XBB.1.5",
  ) 


## model formulae ----

treated_period_variables <- paste0("treatment_period_id", "_", seq_len(length(postbaselinecuts)-1))

# Matching variables

local({
  
  # TODO: make sure these matching variables correspond to the protocol
  
  matching_variables = list()

  # matching set A
  exact <- c(
    "ageband",
    "cv",
    #"sex",
    #"region",
    #"imd_Q5",
    #"multimorb",
    #"prior_covid_infection",
    NULL
  )
  caliper <- c(
    vax_day = 3,
    age = 3,
    prior_vax_interval_bigM = 14,
    NULL
  )
  all <- c(exact, names(caliper))
    matching_variables$A = lst(exact, caliper, all)

  # matching set B
  exact <- c(
    "ageband",
    "cv",
    "sex",
    "region",
    #"imd_Q5",
    "vax_previous_group",
    "multimorb",
    "prior_covid_infection",
    "immunosuppressed",
    NULL
  )
  caliper <- c(
    vax_day = 3,
    age = 3,
    prior_vax_interval_bigM = 28,
    imd = 5000,
    NULL
  )
  all <- c(exact, names(caliper))
  matching_variables$B = lst(exact, caliper, all)

  matching_variables <<- matching_variables

})


# Weighting variables

local({
  
  # TODO: make sure these weighting formulae correspond to the protocol
  
  weighting_formulae = list()
  weighting_variables = list()
  
  # weighting set A

  weighting_formulae$A = "age + cv + sex + region"
  
  weighting_variables$A <- all.vars(as.formula(paste("~", weighting_formulae$A)))
  
  # weighting set B
  
  weighting_formulae$B = "age + cv + sex + region + multimorb"
  
  weighting_variables$B <- all.vars(as.formula(paste("~", weighting_formulae$B)))
  
  # output
  
  weighting_formulae <<- weighting_formulae
  weighting_variables <<- weighting_variables
  
  
})


