library("tidyverse")
library("arrow")
library("here")
library("glue")
library("simDAG")

source(here("analysis", "0-lib", "design.R"))
source(here("analysis", "0-lib", "utility.R"))


#set number of patients
n_sim = 100000

covidvaxstart_date <- as.Date("2020-12-08")
studystart_date <- as.Date(study_dates$studystart_date)
studyend_date <- as.Date(study_dates$studyend_date)

# convenient way to use identity-style functions for root nodes 
# hence "cr": custom root
# this ensures custom functions are treated as root nodes where necessary
# put call to simulate values to the "sim_formula" argument
# put call to transform simulated values to the "trans_formula" argument
# i think this only works with random number generators using an argument named "n" as the number to sample
custom_root = function(n, sim_formula, trans_formula = NULL) {
  sim_call <- sim_formula[[2]]
  sim_call$n <- quote(n)
  sim_x <- eval(sim_call)

  if(is.null(trans_formula)) {
    x <- sim_x
  } else{
    trans_call <- trans_formula[[2]]
    .x <- sim_x
    x <- eval(trans_call)
  }
  return(x)
}

nodecr <- function(name, sim_formula, trans_formula = NULL){
  node(
    name = name, 
    type = custom_root, 
    sim_formula = sim_formula, 
    trans_formula = trans_formula, 
    parents=NULL
  ) 
}

# randomly set some values of a vector to missing
set_missing <- function(x, missing_prop = 0){
  n <- length(x)

  if(length(missing_prop)==1 && missing_prop > 0) {
    missing_index  <- rbernoulli(n, missing_prop)
  } else if(length(missing_prop) == length(x)) {
    missing_index  <- rbernoulli(n, missing_prop)
  } else {
    missing_index <- rep(FALSE, n)
  }

  x_missing <- if_else(missing_index, NA, x)

  return(x_missing)
}

# uniform distribution, relative to an index date
runif_date <- function(n, index_date, min = 0, max, missing_prop=0){
  result <- as.Date(index_date) + floor(runif(n = n, min = min, max = max))
  result_missing <- set_missing(result, missing_prop)
  return(result_missing)
}

# child node for uniformly distributed date variable, with the parent as the index date
node_unif_date <- function(data, parents=NULL, min = 0, max, missing_prop=0){
  stopifnot("you must supply exactly one parent node" = length(parents)==1)
  index_date <- data[[parents]]
  n <- length(index_date)
  runif_date(n = n, index_date = index_date, min = min, max = max, missing_prop = missing_prop)
}

# right-truncated exponential distribution, relative to an index date
rexprt_date <- function(n, index_date, rate, min=0, max = Inf, missing_prop=0){
  result <- as.Date(index_date) +  min + floor(ReIns::rtexp(n = n, rate = rate, endpoint = max))
  result_missing <- set_missing(result, missing_prop)
  return(result_missing)
}
# child node for truncated-exponetially distributed date variable, with the parent as the index date
node_exprt_date <- function(data, parents=NULL, rate, min = 0, max=Inf, missing_prop=0){
  stopifnot("you must supply exactly one parent node" = length(parents)==1)
  index_date <- data[[parents]]
  n <- length(index_date)
  rexprt_date(n = n, index_date = index_date, rate = rate, min = min, max = max, missing_prop = missing_prop)
}

# gamma distribution, relative to an index date
rgamma_date <- function(n, index_date, shape=1, rate=1, missing_prop=0){
  result <- as.Date(index_date) + floor(rgamma(n = n, shape = shape, rate = rate))
  result_missing <- set_missing(result, missing_prop)
  return(result_missing)
}
# child node for gamma distributed date variable, with the parent as the index date
node_gamma_date <- function(data, parents=NULL, shape = 1, rate=1, missing_prop=0){
  stopifnot("you must supply exactly one parent node" = length(parents)==1)
  index_date <- data[[parents]]
  n <- length(index_date)
  rgamma_date(n = n, index_date = index_date, shape = shape, rate = rate, missing_prop = missing_prop)
}

# normalise proportions
normprop <- function(x){
  x/sum(x)
}

# categorical levels

vax_products <- vax_product_lookup[c(productA, productB)]
vax_products_short <- names(vax_products)

ethnicity16_levels <- 
  tribble(
  ~level, ~prob,
  "White - British",  0.5,
  "White - Irish", 0.05,
  "White - Any other White background", 0.05,
  "Mixed - White and Black Caribbean",  0.025,
  "Mixed - White and Black African",  0.025,
  "Mixed - White and Asian",  0.025,
  "Mixed - Any other mixed background",  0.025,
  "Asian or Asian British - Indian",  0.025,
  "Asian or Asian British - Pakistani",  0.025,
  "Asian or Asian British - Bangladeshi",  0.025,
  "Asian or Asian British - Any other Asian background", 0.025,
  "Black or Black British - Caribbean", 0.033,
  "Black or Black British - African", 0.033,
  "Black or Black British - Any other Black background", 0.034,
  "Other Ethnic Groups - Chinese", 0.05,
  "Other Ethnic Groups - Any other ethnic group", 0.05
)

region_levels <- tribble(
  ~level, ~prob, 
  "North East", 0.2,
  "North West", 0.2,
  "Yorkshire and The Humber", 0.3,
  "East Midlands", 0.05,
  "West Midlands", 0.05,
  "East",  0.05,
  "London", 0.05,
  "South East", 0.05,
  "South West", 0.05
  )

# create a set of nodes for prior ascertainament of condition in gp / emergency / admitted / any 
# if no variable is needed for eg emergency, then set emergency_p to NULL
gp_emergency_admitted_prior <- function(name, gp_p=NULL, emergency_p=NULL, admitted_p=NULL){

      if(!is.null(gp_p)){
        gp_node <- node(glue("{name}_prior_gp"), type = "rbernoulli", p = gp_p)
        gp_str <- glue("{name}_prior_gp")
      } else {
        gp_node <- NULL
        gp_str <- FALSE
      }
  
      if(!is.null(emergency_p)){
        emergency_node <- node(glue("{name}_prior_emergency"), type = "rbernoulli", p = emergency_p)
        emergency_str <- glue("{name}_prior_emergency")
      } else {
        emergency_node <- NULL
        emergency_str <- FALSE
      }
  
      if(!is.null(admitted_p)){
        admitted_node <- node(glue("{name}_prior_admitted"), type = "rbernoulli", p = admitted_p)
        admitted_str <- glue("{name}_prior_admitted")
      } else{
        admitted_node <- NULL
        admitted_str <- FALSE
      }
  
      nodes <- list(
        gp_node,
        emergency_node,
        admitted_node,
        node(glue("{name}_prior"), type = "identity", formula = as.formula(glue("~ {gp_str} | {emergency_str} | {admitted_str}")))
      )
  
    Filter(Negate(is.null), nodes)
}

# create a set of nodes for ascertainament of date of occurrence of event as ascertained in gp / emergency / admitted / death / any 
# if no variable is needed for eg emergency, then set emergency_prop to NULL
gp_emergency_admitted_death_dates <- function(name, index_date, max, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.95, death_prop=0.999){

  if(!is.null(gp_prop)){
    gp_str <- glue("{name}_gp_date")
    gp_node <- 
      node(
        gp_str,
        "unif_date",
        parents = index_date,
        max = max,
        missing_prop = gp_prop
      )
    
  } else {
    gp_node <- NULL
    gp_str <- "NA"
  }

  if(!is.null(emergency_prop)){
    emergency_str <- glue("{name}_emergency_date")
    emergency_node <- 
      node(
        emergency_str,
        "unif_date",
        parents = index_date,
        max = max,
        missing_prop = emergency_prop
      )
    
  } else {
    emergency_node <- NULL
    emergency_str <- "NA"
  }

  if(!is.null(admitted_prop)){
    admitted_str <- glue("{name}_admitted_date")
    admitted_node <- 
      node(
        admitted_str,
        "unif_date",
        parents = index_date,
        max = max,
        missing_prop = admitted_prop
      )
  } else {
    admitted_node <- NULL
    admitted_str <- "NA"
  }

  if(!is.null(death_prop)){
    death_str <- glue("{name}_death_date")
    death_node <- 
      node(
        death_str,
        "unif_date",
        parents = index_date,
        max = max,
        missing_prop = death_prop
      )
    
  } else {
    death_node <- NULL
    death_str <- "NA"
  }

  nodes <- list(
      gp_node,
      emergency_node,
      admitted_node,
      death_node,
      node(glue("{name}_date"), type = "identity", formula = as.formula(glue("~ pmin({gp_str}, {emergency_str}, {admitted_str}, {death_str}, na.rm=TRUE)")))
  )

  Filter(Negate(is.null), nodes)
}


# initial list of node variables for the DAG
# these are first specified in a list because simDAG doesn't support `+` if used via functions
# for example `node(...) + node(...) + function_to_create_multiple_+-combined_nodes(...) + ...` doesn't work
# so first put them in a big list, then Reduce(`+`, ...) them
node_list <- list_flatten(list(
  
  nodecr(
    name = "vax_date",
    sim_formula = ~runif_date(studystart_date, max = studyend_date - studystart_date)
  ), 
  node(
    "vax_product", 
    type = "rcategorical", 
    labels = vax_products_short, 
    probs = c(.5, .5)
  ), 
  node(
    "region", 
    type = "rcategorical", 
    labels = region_levels$level,
    probs = region_levels$prob
  ), 
  node(
    "stp",
    type = "rcategorical", 
    labels = 1:36,
    probs = rep(1/36, 36),
    output = "factor"
  ), 
  node(
    "dereg_date", 
    type = "unif_date",
    parents = "vax_date",
    max = 120,
    missing_prop = 0.99
  ), 
  node(
    "registered_previous_6weeks", 
    type = "rbernoulli", 
    p = .999
  ), 
  nodecr(
    "age", 
    sim_formula = ~rnorm(mean = 75, sd = 14),
    trans_formula = ~as.integer(.x)
  ), 
  node(
    "age_eligible",
    type = "identity", 
    formula = ~ age
  ), 
  node(
    "sex", 
    type = "rcategorical", 
    labels = c("female", "male", "intersex", "unknown", NA_character_),
    probs = c(.51, .49, 0, 0, 0.01),
    output = "factor",
    all_levels = TRUE
  ), 
  node(
    "ethnicity16", 
    type = "rcategorical", 
    labels = c(ethnicity16_levels$level, NA),
    probs = c(ethnicity16_levels$prob, 0.05),
    output = "factor",
    all_levels = TRUE
  ), 
  node(
    "ethnicity5", 
    type = "identity", 
    parents = "ethnicity16",
    formula = ~ ethnicity_16_to_5(ethnicity16)
  ),  
  nodecr(
    "imd", 
    sim_formula = ~ runif(1, 32000),
    trans_formula = ~ .x |> plyr::round_any(100) |> as.integer() |> set_missing(0.02)
  ), 
  node(
    "vax_covid_1_date",
    type = "identity", 
    formula = ~ vax_date
  ), 
  node(
    name = "vax_covid_2_date",
    type = "unif_date", 
    parents = "vax_covid_1_date",
    min = 200,
    max = 400,
    missing_prop = 0.2
    #type = "identity",
    #formula = ~ runif(n = n_sim, vax_covid_1_date + 200, vax_covid_1_date + 400) |> set_missing(0.2)
  ), 
  node(
    name = "vax_covid_3_date",
    type = "unif_date", 
    parents = "vax_covid_2_date",
    min = 150,
    max = 400,
    missing_prop = 0.8
  ), 
  node(
    name = "vax_covid_1_product",
    type = "identity",
    formula = ~rcategorical(
      n = n_sim, 
      labels = vax_products_short, 
      probs = c(0.5, 0.5)
    ) |> set_missing(is.na(vax_covid_1_date))
  ), 
  node(
    name = "vax_covid_2_product",
    type = "identity",
    formula = ~rcategorical(
      n = n_sim, 
      labels = vax_products_short, 
      probs = c(0.5, 0.5)
    ) |> set_missing(is.na(vax_covid_2_date))
  ), 
  node(
    name = "vax_covid_3_product",
    type = "identity",
    formula = ~rcategorical(
      n = n_sim, 
      labels = vax_products_short, 
      probs = c(0.5, 0.5)
    ) |> set_missing(is.na(vax_covid_3_date))
  ),  
  node(
    name = "vax_covid_prior_1_date",
    type = "unif_date", 
    parents = "vax_covid_1_date",
    min = -400,
    max = -200,
    missing_prop = 0.1
  ), 
  node(
    name = "vax_covid_prior_2_date",
    type = "unif_date", 
    parents = "vax_covid_prior_1_date",
    min = -400,
    max = -200,
    missing_prop = 0.1
  ), 
  node(
    name = "vax_covid_prior_3_date",
    type = "unif_date", 
    parents = "vax_covid_prior_2_date",
    min = -400,
    max = -200,
    missing_prop = 0.1
  ), 
  node(
    name = "vax_covid_prior_1_product",
    type = "identity",
    formula = ~rcategorical(
      n = n_sim, 
      labels = c("pfizer", "vidprevtyn"), 
      probs = c(0.5, 0.5)
    ) |> set_missing(is.na(vax_covid_prior_1_date))
  ), 
  node(
    name = "vax_covid_prior_2_product",
    type = "identity",
    formula = ~rcategorical(
      n = n_sim, 
      labels = c("pfizer", "moderna"), 
      probs = c(0.5, 0.5)
    ) |> set_missing(is.na(vax_covid_prior_2_date))
  ), 
  node(
    name = "vax_covid_prior_3_product",
    type = "identity",
    formula = ~rcategorical(
      n = n_sim, 
      labels = c("moderna", "az"), 
      probs = c(0.5, 0.5)
    ) |> set_missing(is.na(vax_covid_prior_3_date))
  ), 
  node(
    name = "vax_covid_prior_count",
    type = "identity",
    formula = ~ (!is.na(vax_covid_prior_1_date) + !is.na(vax_covid_prior_2_date) + !is.na(vax_covid_prior_3_date)) + if_else(!is.na(vax_covid_prior_3_date), rpois(n = n_sim, 3), 0L)
  ), 
  node(
    name = "flu_vaccine_before_30_date",
    type = "identity",
    formula = ~ vax_date - as.integer(runif(n = n_sim, 1, 30)) |> set_missing(0.8)
  ), 
  node(
    name = "flu_vaccine_after_30_date",
    type = "identity",
    formula = ~ vax_date + as.integer(runif(n = n_sim, 1, 30)) |> set_missing(0.8)
  ), 
  node(
    name = "flu_vaccine_same_date",
    type = "identity",
    formula = ~ vax_date |> set_missing(0.2)
  ), 
  ## occupation / residency
  node("hscworker", type = "rbernoulli", p = .01), 
  node("care_home_tpp", type = "rbernoulli", p = .01), 
  node("care_home_code", type = "rbernoulli", p = .01), 

  ## baseline clinical variables

  # Green Book variables
  node("crd", type = "rbernoulli", p = .02), 
  node("asthma", type = "rbernoulli", p = .02), 
  node("chd", type = "rbernoulli", p = .02), 
  node("ckd", type = "rbernoulli", p = .02), 
  node("cld", type = "rbernoulli", p = .02), 
  node("cns", type = "rbernoulli", p = .02), 
  node("learndis", type = "rbernoulli", p = .02), 
  node("diabetes", type = "rbernoulli", p = .02), 
  node("immdx", type = "rbernoulli", p = .02), 
  node("immrx", type = "rbernoulli", p = .02), 
  node("dxt_chemo", type = "rbernoulli", p = .02), 
  node(
    "immunosuppressed", type = "identity",
    parents = c("immdx", "immrx", "dxt_chemo"),
    formula = ~ immdx | immrx | dxt_chemo
  ), 
  node("asplenia", type = "rbernoulli", p = .02), 
  node("severe_obesity", type = "rbernoulli", p = .02), 
  node("smi", type = "rbernoulli", p = .02), 

  node("cancer", type = "rbernoulli", p = .01), 
  node("preg_group", type = "rbernoulli", p = .001), 
  node("solid_organ_transplant", type = "rbernoulli", p = .01), 
  node("hiv_aids", type = "rbernoulli", p = .01), 
  
  node("inhospital", type = "rbernoulli", p = .01), 

  ## pre-baseline events

  # ANY EMERGENCY ATTENDANCE
  node("prior_emergency", type = "rbernoulli", p = .05), 
  
  # ANY ADMISSION
  node("prior_admitted", type = "rbernoulli", p = .05), 

  # COVID-RELATED EMERGENCY ATTENDANCE
  node("covid_prior_emergency", type = "rbernoulli", p = .01), 
  
  # COVID-RELATED HOSPITAL ADMISSION
  node("covid_prior_admitted", type = "rbernoulli", p = .01), 
  
  # COVID-RELATED CRITICAL CARE ADMISSION
  node("covid_critcare", type = "rbernoulli", p = .005), 
  
  # GUILLAIN-BARRÉ SYNDROME

  gp_emergency_admitted_prior("sgb", .01, NULL, .01),
  
  # Bell's Palsy
  gp_emergency_admitted_prior("bells_palsy", .01, .01, .01),

  # THROMBOCYTOPENIA
  gp_emergency_admitted_prior("ttp", .01, NULL, .01),
  
  # ACUTE MYOCARDIAL INFARCTION (AMI)
  gp_emergency_admitted_prior("ami", .01, NULL, .01),

  # ISCHAEMIC STROKE
  gp_emergency_admitted_prior("stroke_isch", .01, NULL, .01),
  
  # ARTERIAL THROMBOTIC EVENTS (ATE)
  gp_emergency_admitted_prior("ate", .01, NULL, .01),
  
  # DEEP VEIN THROMBOSIS (DVT)
  gp_emergency_admitted_prior("dvt", .01, NULL, .01),
  
  # INTRACRANIAL VENOUS THROMBOSIS (ICVT)
  gp_emergency_admitted_prior("icvt", .01, NULL, .01),
  
  # PULMONARY EMBOLISM (PE)
  gp_emergency_admitted_prior("pe", .01, NULL, .01),
  
  # VENOUS THROMBOEMBOLISM (VTE)
  gp_emergency_admitted_prior("vte", .01, NULL, .01),
  
  # PERICARDITIS
  gp_emergency_admitted_prior("pericarditis", .01, .01, .01),

  # MYOCARDITIS
  gp_emergency_admitted_prior("myocarditis", .01, .01, .01),
  
  # MENORRHAGIA
  gp_emergency_admitted_prior("menorrhagia", .01, NULL, .01),
  
  # ERYTHEMA MULTIFORME
  gp_emergency_admitted_prior("ery_multi", .01, .01, .01),
  
  # ANAPHYLAXIS
  gp_emergency_admitted_prior("anaphylaxis", .01, .01, .01),
  
  # OTITIS
  node("acute_otitis_gp_prior", type = "rbernoulli", p = .01), 
  
  # CELLULITIS
  node("cellulitis_gp_prior", type = "rbernoulli", p = .01), 

  # Post-baseline events

  ## all-cause outcomes
  node(
    "emergency_date",
    "unif_date",
    parents = "vax_date",
    max = 200,
    missing_prop = 0.8
  ), 
  node(
    "admitted_date",
    "unif_date",
    parents = "vax_date",
    max = 200,
    missing_prop = 0.7
  ), 
  node(
    "death_date",
    "unif_date",
    parents = "vax_date",
    max = 300,
    missing_prop = 0.9
  ), 
  
  
  ## covid outcomes
  node(
    "covid_emergency_date",
    "unif_date",
    parents = "vax_date",
    max = 100,
    missing_prop = 0.8
  ), 
  node(
    "covid_admitted_date",
    "unif_date",
    parents = "vax_date",
    max = 100,
    missing_prop = 0.7
  ), 
  node(
    "covid_critcare_date",
    "identity",
    formula = ~covid_admitted_date |> set_missing(0.7)
  ), 
  node(
    "covid_death_date",
    "identity",
    formula = ~death_date |> set_missing(0.7)
  ), 
  
  ### safety outcomes
  gp_emergency_admitted_death_dates("sgb", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.95, death_prop=0.999),
  gp_emergency_admitted_death_dates("bells_palsy", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = 0.95, admitted_prop=0.95, death_prop=0.999),
  gp_emergency_admitted_death_dates("ttp", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.95, death_prop=0.999),
  gp_emergency_admitted_death_dates("ami", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("stroke_isch", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("ate", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("dvt", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("icvt", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("pe", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("vte", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("pericarditis", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = 0.95, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("myocarditis", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = 0.95, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("menorrhagia", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = NULL, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("ery_multi", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = 0.95, admitted_prop=0.99, death_prop=0.999),
  gp_emergency_admitted_death_dates("anaphylaxis", index_date = "vax_date", max=100, gp_prop = 0.95, emergency_prop = 0.95, admitted_prop=0.99, death_prop=0.999),

  node(
    "acute_otitis_gp_date",
    "unif_date",
    parents = "vax_date",
    max = 100,
    missing_prop = 0.95
  ), 

  node(
    "cellulitis_gp_date",
    "unif_date",
    parents = "vax_date",
    max = 100,
    missing_prop = 0.95
  )

))

# initiate dag from the list of nodes
dag <- Reduce(f = `+`, x = node_list, init = empty_dag())


# simulate data

# set seed
set.seed(20240930)
dummydata <- sim_from_dag(dag, n_sim = n_sim)

dummydata <- 
  dummydata |>
  # convert product codes to long form in tpp
  mutate(
    across(
      c(
        vax_product, 
        vax_covid_1_product, vax_covid_2_product, vax_covid_3_product, 
        vax_covid_prior_1_product, vax_covid_prior_2_product, vax_covid_prior_3_product
      ),
      ~recode_values(., from = vax_products_short, to = vax_products)
    )
  ) |>
  # add patient ID column
  add_column(patient_id = seq_len(n_sim), .before = 1L)

# save data
write_feather(dummydata, sink = here("analysis", "1-extract", "dummy_extract.arrow"))
