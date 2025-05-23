
library('tidyverse')
library('arrow')
library('here')
library('glue')


remotes::install_github("https://github.com/wjchulme/dd4d")
library('dd4d')


population_size <- 10000



source(here("analysis", "0-lib", "design.R"))
source(here("analysis", "0-lib", "utility.R"))


covidvaxstart_date <- as.Date("2020-12-08")
studystart_date <- as.Date(study_dates$studystart_date)
studyend_date <- as.Date(study_dates$studyend_date)
followupend_date <- as.Date(study_dates$followupend_date)
index_date <- studystart_date

index_day <- 0L
covidvaxstart_day <- as.integer(covidvaxstart_date - index_date)
studystart_day <- as.integer(studystart_date - index_date)
studyend_day <- as.integer(studyend_date - index_date)


known_variables <- c(
  "index_date", "studystart_date", "studyend_date", "covidvaxstart_date",
  "index_day",  "studystart_day", "studyend_day", "covidvaxstart_day",
  "vax_product_lookup", "productA", "productB",
  NULL
)

sim_list = lst(

  vax_day = bn_node(
    ~(runif(n=..n, studystart_day, studyend_day)),
  ),

  vax_product = bn_node(
    ~rcat(n=..n, vax_product_lookup[c(productA, productB)], c(0.5,0.5)),
    needs = "vax_day",
  ),

  region = bn_node(
    variable_formula = ~rfactor(n=..n, levels=c(
      "North East",
      "North West",
      "Yorkshire and The Humber",
      "East Midlands",
      "West Midlands",
      "East",
      "London",
      "South East",
      "South West"
    ), p = c(0.2, 0.2, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05))
  ),

  stp = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 36)), levels=1:36)
  ),

  # msoa = bn_node(
  #   ~factor(as.integer(runif(n=..n, 1, 100)), levels=1:100),
  #   missing_rate = ~ 0.005
  # ),

  dereg_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+120)),
    missing_rate = ~0.99
  ),

  registered_previous_6weeks = bn_node(
    ~rbernoulli(n=..n, p=0.999)
  ),

  # practice_id = bn_node(
  #   ~as.integer(runif(n=..n, 1, 200))
  # ),

  age = bn_node(
    ~as.integer(rnorm(n=..n, mean=80, sd=14))
  ),

  age_eligible = bn_node(~age),

  sex = bn_node(
    ~rfactor(n=..n, levels = c("female", "male", "intersex", "unknown"), p = c(0.51, 0.49, 0, 0)),
    missing_rate = ~0.001 # this is shorthand for ~(rbernoulli(n=..n, p = 0.2))
  ),

  ethnicity5 = bn_node(variable_formula = ~ ethnicity_16_to_5(ethnicity16), needs = "ethnicity16"),
  ethnicity16 = bn_node(
    variable_formula = ~ rfactor(
      n = ..n,
      levels = c(
        "White - British",
        "White - Irish",
        "White - Any other White background",
        "Mixed - White and Black Caribbean",
        "Mixed - White and Black African",
        "Mixed - White and Asian",
        "Mixed - Any other mixed background",
        "Asian or Asian British - Indian",
        "Asian or Asian British - Pakistani",
        "Asian or Asian British - Bangladeshi",
        "Asian or Asian British - Any other Asian background",
        "Black or Black British - Caribbean",
        "Black or Black British - African",
        "Black or Black British - Any other Black background",
        "Other Ethnic Groups - Chinese",
        "Other Ethnic Groups - Any other ethnic group"
      ),
      p = c(
        0.5, 0.05, 0.05, # White
        0.025, 0.025, 0.025, 0.025, # Mixed
        0.025, 0.025, 0.025, 0.025, # Asian
        0.033, 0.033, 0.034, # Black
        0.05, 0.05 # Other
      )
    ),
    missing_rate = ~0.1,
  ),

  imd = bn_node(
    ~as.integer(plyr::round_any(runif(n=..n, 1, 32000), 100)),
    missing_rate = ~0.02
  ),

  ### vaccination variables

  # covid vaccines
  vax_covid_1_day = bn_node(~vax_day),
  vax_covid_2_day = bn_node(
    ~runif(n=..n, vax_covid_1_day+200, vax_covid_1_day+400),
    missing_rate = ~0.2,
    needs = "vax_covid_1_day"
  ),
  vax_covid_3_day = bn_node(
    ~runif(n=..n, vax_covid_2_day+150, vax_covid_2_day+400),
    missing_rate = ~0.8,
    needs = "vax_covid_2_day"
  ),
  vax_covid_1_product = bn_node(~rcat(n=..n, vax_product_lookup[c(productA, productB)], c(0.5,0.5)), needs = "vax_covid_1_day"),
  vax_covid_2_product = bn_node(~if_else(runif(..n)<0.98, vax_product_lookup[productA], vax_product_lookup[productB]), needs = "vax_covid_2_day"),
  vax_covid_3_product = bn_node(~rcat(n=..n, vax_product_lookup[c("moderna", "pfizer")], c(0.5,0.5)), needs = "vax_covid_3_day"),
  
  vax_covid_prior_1_day = bn_node(
    ~runif(n=..n, vax_covid_1_day-400, vax_covid_1_day-200),
    missing_rate = ~0.1,
    needs = "vax_covid_1_day"
  ),
  vax_covid_prior_2_day = bn_node(
    ~runif(n=..n, vax_covid_prior_1_day-400, vax_covid_prior_1_day-200),
    missing_rate = ~0.1,
    needs = "vax_covid_prior_1_day"
  ),
  vax_covid_prior_3_day = bn_node(
    ~runif(n=..n, vax_covid_prior_2_day-400, vax_covid_prior_2_day-200),
    missing_rate = ~0.1,
    needs = "vax_covid_prior_2_day"
  ),
  
  vax_covid_prior_1_product = bn_node(~rcat(n=..n, vax_product_lookup[c("pfizer", "vidprevtyn")], c(0.5,0.5)), needs = "vax_covid_prior_1_day"),
  vax_covid_prior_2_product = bn_node(~rcat(n=..n, vax_product_lookup[c("pfizer", "moderna")], c(0.5,0.5)), needs = "vax_covid_prior_2_day"),
  vax_covid_prior_3_product = bn_node(~rcat(n=..n, vax_product_lookup[c("moderna","az")], c(0.5,0.5)), needs = "vax_covid_prior_3_day"),

  vax_covid_prior_count = bn_node(
    ~ (!is.na(vax_covid_prior_1_day) + !is.na(vax_covid_prior_2_day) + !is.na(vax_covid_prior_3_day)) + if_else(!is.na(vax_covid_prior_3_day), rpois(n=..n, 3), 0L)
  ),
  
  ## occupation / residency
  hscworker = bn_node(
    ~rbernoulli(n=..n, p=0.01)
  ),
  care_home_tpp = bn_node(
    ~rbernoulli(n=..n, p = 0.01)
  ),

  care_home_code = bn_node(
    ~rbernoulli(n=..n, p = 0.01)
  ),

  ## baseline clinical variables

  asthma = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  asthma_simple = bn_node( ~asthma),
  cns = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  crd = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  severe_obesity = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  diabetes = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  smi = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chd = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  ckd = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  cld = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  cancer = bn_node( ~rbernoulli(n=..n, p = 0.01)),

  preg_group = bn_node( ~rbernoulli(n = ..n, p = 0.001)),

  immdx = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  immrx = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  dxt_chemo = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  immunosuppressed = bn_node( ~immdx | immrx | dxt_chemo),
  asplenia = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  solid_organ_transplant = bn_node( ~rbernoulli(n=..n, p = 0.01)),
  hiv_aids = bn_node( ~rbernoulli(n=..n, p = 0.01)),

  learndis = bn_node( ~rbernoulli(n=..n, p = 0.02)),

  inhospital = bn_node( ~rbernoulli(n=..n, p = 0.01)),

  ## pre-baseline events where event date is relevant

  covid_admitted_0_day = bn_node(
    ~as.integer(runif(n=..n, vax_day-100, vax_day-1)),
    missing_rate = ~0.99
  ),

  ## post-baseline events (outcomes)

  
  ### all-cause outcomes
  emergency_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  admitted_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.7
  ),
  
  death_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+200)),
    missing_rate = ~0.90
  ),

  
  ### covid outcomes
  
  covid_emergency_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),

  covid_admitted_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.7
  ),

  covid_critcare_day = bn_node(
    ~covid_admitted_day,
    needs = "covid_admitted_day",
    missing_rate = ~0.7
  ),

  covid_death_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),

  
  
  ### safety outcomes
  
  
  #### AMI
  
  ami_gp_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  ami_emergency_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  ami_admitted_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.7
  ),
  
  ami_death_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),
  
  ami_day = bn_node(
    ~ pmin(ami_gp_day, ami_emergency_day, ami_admitted_day, ami_death_day),
  ),
  
  
  #### Ischaemic stroke
  
  stroke_isch_gp_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  stroke_isch_emergency_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  stroke_isch_admitted_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.7
  ),
  
  stroke_isch_death_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),
  
  stroke_isch_day = bn_node(
    ~ pmin(stroke_isch_gp_day, stroke_isch_emergency_day, stroke_isch_admitted_day, stroke_isch_death_day),
  ),
  
  
  #### ATE
  
  ate_gp_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  ate_emergency_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.8
  ),
  
  ate_admitted_day = bn_node(
    ~as.integer(runif(n=..n, vax_day, vax_day+100)),
    missing_rate = ~0.7
  ),
  
  ate_death_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),
  
  ate_day = bn_node(
    ~ pmin(ate_gp_day, ate_emergency_day, ate_admitted_day, ate_death_day),
  ),
  
)
bn <- bn_create(sim_list, known_variables = known_variables)

bn_plot(bn)
bn_plot(bn, connected_only=TRUE)

set.seed(10)

dummydata <-bn_simulate(bn, pop_size = population_size, keep_all = FALSE, .id="patient_id")

dummydata_processed <- dummydata %>%
  #convert integer days to dates since index date and rename vars
  mutate(across(ends_with("_day"), ~ as.Date(as.character(index_date + .)))) %>%
  rename_with(~str_replace(., "_day", "_date"), ends_with("_day"))


fs::dir_create(here("output", "1-extract"))
write_feather(dummydata_processed, sink = here("output", "1-extract", "dummy_extract.arrow"))



