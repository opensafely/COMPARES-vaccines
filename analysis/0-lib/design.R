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
library("tidyverse")
library("here")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## create output directories ----
fs::dir_create(here("analysis", "0-lib"))


# key study dates
# The dates are defined in json format so they can be read in by R and python scripts
# json has no easy way to comment, so explanation for dates is here:
# - firstpossibleax_date is the date from which we want to identify covid vaccines. the mass vax programme was 8 Dec 2020 but other people were vaccinated earlier in trials, so want to pick these people up too (and possibly exclude them)
# - start_date is when we start the observation period, typically at the start of the vaccination campaign (we may restrict the analysis to those vaccinated at a later date if only one product wasavailable at the start)
# - end_date is when we stop the observation period. Most likely the day before the start of the next campaign.

study_dates <-
  jsonlite::read_json(path = here("analysis", "0-lib", "study-dates.json")) |>
  map(as.Date)

# make these available in the global environment
# so we don't have to use `study_dates$start_date` or `start_date <- study_dates$start_date` in each script
list2env(study_dates, globalenv())


# define outcomes ----

# maximum follow-up duration, in days
# currently up to 24 weeks after time-zero
maxfup <- 7L * 24L

# where to split follow-up time after recruitment
postbaselinecuts <- as.integer(c(0, 7, 14, 28, 56, 84, 112, 140, 168))

stopifnot("last cut in postbaseline follow up splits should be same as max follow up" = max(postbaselinecuts)==maxfup)

# Look-up table for all potential outcomes of interest
# this look-up table is used by `create-project.R` to select outcomes that are used in comparisons between products
# comment-out outcomes that are not yet defined in the dataset definition or are otherwise not needed
# - `event` is the variable name prefix used in the dataset, eg for the  `covid_admission_date` variable, the prefix is `covid_admission`
# - `event_var` is the full variable name
# - `event_descr` is the verbose variable name used in summary table/figure outputs
events_lookup <- tribble(
  ~event_type     , ~event            , ~event_var             , ~event_descr                                ,

  # all-cause
  "all-cause"     , "emergency"       , "emergency_date"       , "A&E attendance"                            ,
  "all-cause"     , "admitted"        , "admitted_date"        , "Unplanned hospitalisation"                 ,
  "all-cause"     , "death"           , "death_date"           , "Any death"                                 ,

  # Effectiveness
  "effectiveness" , "covid_emergency" , "covid_emergency_date" , "COVID-19 A&E attendance"                   ,
  "effectiveness" , "covid_admitted"  , "covid_admitted_date"  , "COVID-19 hospitalisation"                  ,
  "effectiveness" , "covid_critcare"  , "covid_critcare_date"  , "COVID-19 critical care"                    ,
  "effectiveness" , "covid_death"     , "covid_death_date"     , "COVID-19 death"                            ,

  # safety
  "safety"        , "sgb"             , "sgb_date"             , "Guillain-Barré syndrome"                   ,
  "safety"        , "bells_palsy"     , "bells_palsy_date"     , "Bell's palsy"                              ,
  "safety"        , "ttp"             , "ttp_date"             , "Thrombocytopenia"                          ,
  "safety"        , "ami"             , "ami_date"             , "Acute myocardial infarction"               ,
  "safety"        , "stroke_isch"     , "stroke_isch_date"     , "Ischaemic stroke"                          ,
  "safety"        , "ate"             , "ate_date"             , "Composite arterial thrombotic event (ATE)" ,
  "safety"        , "dvt"             , "dvt_date"             , "Deep vein thrombosis (DVT)"                ,
  "safety"        , "icvt"            , "icvt_date"            , "Intracranial venous thrombosis (ICVT)"     ,
  "safety"        , "pe"              , "pe_date"              , "Pulmonary embolism (PE)"                   ,
  "safety"        , "vte"             , "vte_date"             , "Composite venous thrombotic event (VTE)"   ,
  "safety"        , "pericarditis"    , "pericarditis_date"    , "Pericarditis"                              ,
  "safety"        , "myocarditis"     , "myocarditis_date"     , "Myocarditis"                               ,
  "safety"        , "menorrhagia"     , "menorrhagia_date"     , "Heavy menstrual bleeding"                  ,
  "safety"        , "ery_multi"       , "ery_multi_date"       , "Erythema multiforme"                       ,
  "safety"        , "anaphylaxis"     , "anaphylaxis_date"     , "Anaphylaxis"                               ,

  # negative control
  "negative"      , "noncovid_death"  , "noncovid_death_date"  , "Non-COVID-19 death"                        ,
  "negative"      , "acute_otitis_gp" , "acute_otitis_gp_date" , "Acute otitis media"                        ,
  "negative"      , "cellulitis_gp"   , "cellulitis_gp_date"   , "Cellulitis"
)

# use this to help calculate which campaign each person's prior vaccine occurred in
campaigns_lookup <-
  tribble(
    ~campaign            , ~start_date  ,
    "Before Autumn 2023" , "2020-04-01" ,
    "Autumn 2023"        , "2023-08-28" ,
    "Spring 2024"        , "2024-04-15" ,
    "Autumn 2024"        , "2024-09-30" ,
  ) |>
  mutate(
    start_date = as.Date(start_date),
    end_date = lead(start_date) - 1L,
  )
  


# Analysis population definition with respect to influenza coadministration.
# Options:
# - TRUE: restrict the main analysis to individuals who received influenza vaccination on the same day as the COVID-19 vaccine.
# - FALSE: include all eligible individuals regardless of influenza coadministration.
#
# Note: flu coadministration descriptive analyses are always performed among all eligible individuals (i.e. regardless of this setting).
restrict_to_coadmin_flu <- TRUE 


# the products of interest this campaign

productA <- "pfizer_JN1"
productB <- "moderna_JN1"

# look-up to rename coding-friendly product names to publication-friendly product names
vax_shortname_lookup <- c(
  "BNT162b2/JN.1" = "pfizer_JN1",
  "mRNA-1273/JN.1" = "moderna_JN1"
)

# output from https://jobs.opensafely.org/opensafely-internal/tpp-vaccination-names/ workspace
# shows all possible covid vaccination names in TPP

# lookup to rename TPP product names to coding-friendly product names
vax_product_lookup <- c(
  # Pfizer adult
  "pfizer_original" = "COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)",
  "pfizer_BA1" = "Comirnaty Original/Omicron BA.1 COVID-19 Vacc md vials",
  "pfizer_BA45" = "Comirnaty Original/Omicron BA.4-5 COVID-19 Vacc md vials",
  "pfizer_XBB15" = "Comirnaty Omicron XBB.1.5 COVID-19 Vacc md vials",
  "pfizer_JN1" = "Comirnaty JN.1 COVID-19 mRNA Vaccine 0.3ml inj md vials (Pfizer Ltd)",
  "pfizer_LP81" = "Comirnaty LP.8.1 COVID-19 Vacc 30microg/0.3ml dose inj pfs (Pfizer Ltd)",
  "pfizer_KP2" = "Comirnaty KP.2 COVID-19 Vacc 30microg/0.3ml dose inj md vial (Pfizer Ltd)",
  "pfizer_KP2_pfs" = "Comirnaty KP.2 COVID-19 Vacc 30microg/0.3ml dose inj pfs (Pfizer Ltd)",

  "pfizer_unspecified" = "Comirnaty COVID-19 mRNA Vacc ready to use 0.3ml inj md vials",

  # Pfizer children

  "pfizer_original_children" = "COVID-19 mRNA Vaccine Comirnaty Children 5-11yrs 10mcg/0.2ml dose conc for disp for inj MDV (Pfizer)",
  "pfizer_JN1_children" = "Comirnaty JN.1 Children 5-11yrs COVID-19 Vacc 0.3ml sd vials (Pfizer Ltd)",
  "pfizer_XBB15_children" = "Comirnaty Omicron XBB.1.5 Child 5-11y COVID-19 Vacc md vials",
  "pfizer_LP81_children" = "Comirnaty LP.8.1 Children 5-11y COVID-19 Vacc 0.3ml sd vials  (Pfizer Ltd)", # note double space before "(Pfizer Ltd)"

  "pfizer_original_under5" = "Comirnaty Children 6m-4yrs COVID-19 mRNA Vacc 0.2ml md vials",
  "pfizer_JN1_under5" = "Comirnaty JN.1 Children 6m-4yrs COVID-19 Vacc 0.3ml md vials (Pfizer Ltd)",
  "pfizer_XBB15_under5" = "Comirnaty Omicron XBB.1.5 Child 6m-4y COVID-19 Vacc md vials",
  "pfizer_LP81_under5" = "Comirnaty LP.8.1 Children 6m-4y COVID-19 Vacc 0.3m md vials  (Pfizer Ltd)", # note, double space before "(Pfizer Ltd)"

  # Astrazeneca

  "az_original" = "COVID-19 Vaccine Vaxzevria 0.5ml inj multidose vials (AstraZeneca)",
  "az_original_half" = "COVID-19 Vac AZD2816 (ChAdOx1 nCOV-19) 3.5x10*9 viral part/0.5ml dose sol for inj MDV (AstraZeneca)",

  # Moderna

  "moderna_original" = "COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)",
  "moderna_omicron" = "COVID-19 Vac Spikevax (Zero)/(Omicron) inj md vials",
  "moderna_BA45" = "COVID-19 Vacc Spikevax Orig/Omicron BA.4/BA.5 inj md vials",
  "moderna_XBB15" = "COVID-19 Vacc Spikevax (XBB.1.5) 0.1mg/1ml inj md vials",
  "moderna_JN1" = "Spikevax JN.1 COVID-19 Vacc 0.1mg/ml inj md vials (Moderna, Inc)",
  "moderna_omicron2" = "COVID-19 Vaccine Moderna (mRNA-1273.529) 50micrograms/0.25ml dose sol for inj MDV",
  "moderna_unspecified" = "COVID-19 Vaccine Moderna 0.5ml dispersion for inj vials",

  # Sanofi-GSK
  "sanofigsk_B1" = "COVID-19 Vacc VidPrevtyn (B.1.351) 0.5ml inj multidose vials",
  "sanofigsk_D614" = "COVID-19 Vac Sanofi (CoV2 preS dTM monovalent D614 (recombinant)) 5mcg/0.5ml dose susp for inj MDV",
  "sanofigsk_D614B1" = "COVID-19 Vacc Sanofi (D614+B.1.351) 0.5ml inj md vials",


  # Novavax
  "novavax" = "COVID-19 Vac Nuvaxovid (recombinant, adj) 5micrograms/0.5ml dose susp for inj MDV (Novavax CZ a.s.)",

  # Sputnik
  "sputnik_i_multi" = "COVID-19 Vacc Sputnik V Component I 0.5ml multidose vials",
  "sputnik_ii_multi" = "COVID-19 Vacc Sputnik V Component II 0.5ml multidose vials",
  "sputnik_i_inj" = "COVID-19 Vaccine Sputnik V Component I 0.5ml inj vials",
  "sputnik_ii_inj" = "COVID-19 Vaccine Sputnik V Component II 0.5ml inj vials",

  # Janssen
  "jansenn" = "COVID-19 Vaccine Janssen (Ad26.COV2-S (recomb)) 0.5ml dose solution for injection multidose vials",

  # Sinopharm
  "sinopharm" = "COVID-19 Vac Sinopharm BIBP (inactivated adjuvanted) 6.5U/0.5ml dose susp for inj vials",

  # Valneva
  "valneva" = "COVID-19 Vaccine Valneva (inactivated adj whole virus) 40antigen units/0.5ml dose susp for inj MDV",

  # Medicago
  "medicago" = "COVID-19 Vaccine Medicago (CoVLP) 3.75micrograms/0.5ml dose emulsion for injection multidose vials",

  # Convidecia
  "convidecia" = "COVID-19 Vaccine Convidecia 0.5ml inj vials",

  # Covaxin
  "covaxin" = "COVID-19 Vac Covaxin (NIV-2020-770 inactivated) 6micrograms/0.5ml dose susp for inj MDV",

  # Coronavac
  "coronavac" = "COVID-19 Vac CoronaVac (adjuvanted) 600U/0.5ml dose susp for inj vials",

  # Covishield
  "covishield" = "COVID-19 Vac Covishield (ChAdOx1 S recombinant) 5x10*9 viral particles/0.5ml dose sol for inj MDV",

  # Covovax
  "covovax" = "COVID-19 Vac Covovax (adjuvanted) 5micrograms/0.5ml dose susp for inj MDV (Serum Institute of India)",

  # Not specified
  "unspecified" = "SARS-2 Coronavirus vaccine"

)




treatment_lookup <-
  tribble(
    ~treatment, ~treatment_descr,
    "pfizer_JN1", "pfizer/JN.1",
    "moderna_JN1", "mRNA-1273/JN.1",
  )


# redaction threshold
sdc.limit <- 6L

## lookups to convert coded variables to full, descriptive variables ----


recoder <-
  lst(
    cohort = c(
      `Clinically vulnerable` = "cv",
      `Aged 65 years or over` = "age65plus",
      `Care home residency` = "carehome"
    ),
    subgroups = c(
      `Main` = "all",
      `Age` = "ageband",
      `Clinically at-risk` = "cv",
      `Ethnicity` = "ethnicity",
      `Prior COVID-19 vaccine count` = "prior_vax_count_group",
      NULL
    ),
    status = c(
      `Unmatched` = "unmatched",
      `Matched` = "matched"
    ),
    treatment = c(
      `pfizer/JN.1` = "0",
      `mRNA-1273/JN.1` = "1"
    ),
    outcome = set_names(events_lookup$event, events_lookup$event_descr),
    all = c(` ` = "all"),
    ageband = c(
      "50-64", "65-74", "75-79", "80-84", "85+"
    ) %>%
      {
        set_names(., .)
      },
    cv = c(
      `Clinically at-risk` = "TRUE",
      `Not clinically at-risk` = "FALSE"
    ),
    prior_vax_count_group = c(
      "0", "1-2", "3-5", "6+"
    ) %>%
      {
        set_names(., .)
      },
  )

## variable labels in tables and plots ----

variables_lookup <-
  tribble(
    ~variable                , ~variable_label                                             ,
    "N"                      , "Total N"                                                   ,
    "treatment_descr"        , "Vaccine product"                                           ,

    #  covariates
    #"vax_date"               , "Vaccination date"                                          ,
    "vax_day"                , "Vaccination day"                                           ,
    "prior_vax_interval"     , "Days since previous vaccine"                               ,
    "prior_vax_count_group"  , "Previous vaccine count"                                    ,
    "age_eligible"           , "Age"                                                       ,
    "ageband"                , "Age band"                                                  ,
    "sex"                    , "Sex"                                                       ,
    "ethnicity5"             , "Ethnicity"                                                 ,
    "imd_Q5"                 , "Deprivation"                                               ,
    "region"                 , "Region"                                                    ,
    #"housebound"            , "Clinically housebound"                                     ,
    #"care_home_combined"    , "Care/nursing home resident"                                ,

    # green book conditions
    "cv"                     , "Clinically at-risk"                                        ,

    "crd"                    , "Chronic respiratory disease"                               ,
    "chd"                    , "Chronic heart disease"                                     ,
    "ckd"                    , "Chronic kidney disease"                                    ,
    "cld"                    , "Chronic liver disease"                                     ,
    "cns"                    , "Chronic neurological disease"                              ,
    "learndis"               , "Learning disabilities"                                     ,
    "diabetes"               , "Diabetes"                                                  ,
    "immunosuppressed"       , "Immunosuppressed"                                          ,
    "asplenia"               , "Aplenia or poor spleen function"                           ,
    "severe_obesity"         , "Body Mass Index > 40 kg/m^2"                               ,
    "smi"                    , "Serious mental illness"                                    ,
    # "asthma"                , "Asthma"                                                    ,

    # "immuno_any", "Immunosuppressed (all)",
    # "immdx", "Immunocompromising diagnosis",
    # "immrx", "Immunosuppressive medications, previous 3 years",
    # "dxt_chemo", "Chemotherapy, previous 3 years",
    # "cancer", "Cancer, previous 3 years",
    # "solid_organ_transplant", "Solid organ transplant",
    # "hiv_aids", "HIV/AIDS",

    "cv_immuno"              , "Risk group type"                                           ,
    "cv_count"               , "Risk group count"                                          ,
    "cv_count_grp"           , "Risk group count (cat)"                                    ,
    "cv_notimmuno_count_grp" , "Non-immuno risk group count (cat)"                         ,

    # prior events

    "prior_vax_campaign"     , "Campaign of prior vaccine"                                 ,

    # COVID-related
    "covid_prior_emergency"  , "Prior (<1 year) COVID-19 A&E attendance"                   ,
    "covid_prior_admitted"   , "Prior (<1 year) COVID-19 hospitalisation"                  ,
    "covid_prior_critcare"   , "Prior (<1 year) COVID-19 critical care"                    ,

    # Safety
    "prior_emergency"        , "Prior (<1 year) A&E attendance"                            ,
    "prior_admitted"         , "Prior (<1 year) hospitalisation"                           ,

    "sgb_prior"              , "Prior (<1 year) Guillain-Barré syndrome"                   ,
    "bells_palsy_prior"      , "Prior (<1 year) Bell's palsy"                              ,
    "ttp_prior"              , "Prior (<1 year) Thrombocytopenia"                          ,
    "ami_prior"              , "Prior (<1 year) Acute myocardial infarction"               ,
    "stroke_isch_prior"      , "Prior (<1 year) Ischaemic stroke"                          ,
    "ate_prior"              , "Prior (<1 year) Composite arterial thrombotic event (ATE)" ,
    "dvt_prior"              , "Prior (<1 year) Deep vein thrombosis (DVT)"                ,
    "icvt_prior"             , "Prior (<1 year) Intracranial venous thrombosis (ICVT)"     ,
    "pe_prior"               , "Prior (<1 year) Pulmonary embolism (PE)"                   ,
    "vte_prior"              , "Prior (<1 year) Composite venous thrombotic event (VTE)"   ,
    "pericarditis_prior"     , "Prior (<1 year) Pericarditis"                              ,
    "myocarditis_prior"      , "Prior (<1 year) Myocarditis"                               ,
    "menorrhagia_prior"      , "Prior (<1 year) Heavy menstrual bleeding"                  ,
    "ery_multi_prior"        , "Prior (<1 year) Erythema multiforme"                       ,
    "anaphylaxis_prior"      , "Prior (<1 year) Anaphylaxis"                               ,

    # "tests_cat_prior"       , "Prior number of SARS-CoV-2 tests"                          ,
    # "covid_infection_prior" , "Prior documented SARS-CoV-2 infection"                     ,
  )

# for passing to table1_summary functions
variable_labels <- as.list(set_names(
  variables_lookup$variable_label,
  variables_lookup$variable
))

## model formulae ----

treated_period_variables <- paste0("treatment_period_id", "_", seq_len(length(postbaselinecuts) - 1))

# Matching variables

local({
  # TODO: make sure these matching variables correspond to the protocol

  matching_variables <- list()

  # matching specification A
  exact <- c(
    "ageband",
    "cv",
    "immunosuppressed",
    # "sex",
    # "region",
    # "imd_Q5",
    # "cv_count",
    NULL
  )
  caliper <- c(
    vax_day = 3,
    age = 3,
    #prior_vax_interval_bigM = 14,
    NULL
  )
  all <- c(exact, names(caliper))
  matching_variables$A <- lst(exact, caliper, all)

  # matching specification B
  exact <- c(
    "ageband",
    "cv",
    "sex",
    #"imd_Q5",
    "prior_vax_count_group",
    "immunosuppressed",
    "cv_notimmuno_count_grp",
    NULL
  )
  caliper <- c(
    vax_day = 3,
    age = 3,
    #prior_vax_interval_bigM = 28,
    imd = 5000,
    NULL
  )
  all <- c(exact, names(caliper))
  matching_variables$B <- lst(exact, caliper, all)

  matching_variables <<- matching_variables

})


# Weighting variables

local({
  # TODO: make sure these weighting formulae correspond to the protocol
  # TODO: make sure subgroups are dealt with as desired, either:
  # - ignore, and assume appropriate calibration from global model
  # - independent models in each subgroup (so refit weighting model(s) for each subgroup analysis)
  # - completely stratified (subgroup1*subgroup2*subgroup3*...)
  # - something else

  weighting_formulae <- list()
  weighting_variables <- list()

  # weighting specification A

  weighting_formulae$A <- "vax_day + age + sex + imd_Q5 + cv"

  weighting_variables$A <- all.vars(as.formula(paste("~", weighting_formulae$A)))

  # weighting specification B

  weighting_formulae$B <- "vax_day + age + sex + imd_Q5 + immunosuppressed + cv_notimmuno_count_grp"

  weighting_variables$B <- all.vars(as.formula(paste("~", weighting_formulae$B)))

  # output

  weighting_formulae <<- weighting_formulae
  weighting_variables <<- weighting_variables


})



# lmw variables (experimental)

local({
  lmw_formulae <- list()
  lmw_variables <- list()

  # weighting specification A

  lmw_formulae$A <- "vax_day + age + sex + cv"

  lmw_variables$A <- all.vars(as.formula(paste("~", lmw_formulae$A)))

  # weighting specification B

  lmw_formulae$B <- "vax_day + age + sex + imd_Q5 + multimorb"

  lmw_variables$B <- all.vars(as.formula(paste("~", lmw_formulae$B)))

  # output

  lmw_formulae <<- lmw_formulae
  lmw_variables <<- lmw_variables
})



## meta-parameter dataset ---
## contains all combinations of design elements that should be run in the study
## this is mainly used to construct the project.yaml file, via the create-project.R script
## but also picked up in a few places elsewhere

# TODO: check this is capturing everything

metaparams <-
  expand_grid(
    cohort = factor(c("age65plus", "cv")),
    method = factor(c("match", "weight", "lmw")),
    spec = c("A", "B", "C"),
    subgroup = factor(recoder$subgroups),
    outcome = factor(recoder$outcome),
  ) |>
  mutate(
    cohort_descr = fct_recoderelevel(cohort, recoder$cohort),
    outcome_descr = fct_recoderelevel(outcome, recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup, recoder$subgroups),
  ) |>
  # select only some parameters to avoid too much dev work
  filter(
    #    cohort == "age65plus",
    method != "lmw",
    subgroup %in% c("all", "ageband"),
    spec %in% c("A", "B")
  )
