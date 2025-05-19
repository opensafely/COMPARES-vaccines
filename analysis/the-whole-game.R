
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Demonstrate in a single script the basics of the study design
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('survival')
library("MatchIt")
library("WeightIt")
library("cobalt")
library("splines")
library("sandwich")
library("marginaleffects")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

# simulation some example data

rcat <- function (n, levels, p, factor=TRUE) {
  x <- sample(x = levels, size = n, replace = TRUE, prob = p)
  if(factor) x <- factor(x, levels=levels)
  x
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# create simulated dataset with N patients
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set.seed <- 42

N <- 1000

maxfup <- 365

data_cohort <- 
  tibble(
    patient_id = seq_len(N),
    age = as.integer(rnorm(N, 60, 10)),
    sex = as.factor(rcat(N, c("female", "male"), p=c(0.5,0.5))),
    vax_date = runif(N, 0, 100),
    vax_product = rcat(N, c("A", "B"), p=c(0.5,0.5)),
    treatment = (vax_product == "B")*1L,
    outcome_date = vax_date + runif(N, 1, 1000),
    censor_date = vax_date + runif(N, 1, 1000),
    event_time = as.integer((pmin(outcome_date, censor_date, 365) - vax_date) + 1),
    event_indicator = outcome_date<=pmin(censor_date, maxfup)
  )

## preliminary model inputs

# formula for predicting treatment probability
weighting_formula <- treatment ~ vax_date + age + sex

matching_vars_exact <- c("sex")
matching_vars_caliper <- c(vax_date = 7, age=3)


# outcome model formula, assuming balancing weights are applied
formula_time_treatment <- event_indicator ~ treatment + ns(time, 4) + treatment:ns(time, 4)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# as per weight.R script
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# fit a standard propensity score logistic regression model with weightit
weightit_obj <- 
  weightit(
    formula = weighting_formula,
    data = data_cohort,
    method = "glm", 
    estimand = "ATE",
    #stabilize = TRUE
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# as per match.R script
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# perform 1-1 matching without replacement, using matchit

matchit_obj <- 
  matchit(
    formula = treatment ~ 1,
    data = data_cohort,
    method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
    replace = FALSE,
    estimand = "ATT", # since we are doing exact matching, ATT is equivalent to ATU. although we'll actually get the ATO (average treatment in the overlap)
    exact = matching_vars_exact,
    caliper = matching_vars_caliper, std.caliper=FALSE,
    m.order = "data", # data is sorted on (effectively random) patient ID
    #verbose = TRUE,
    ratio = 1L # could also consider exact matching only, with n:m ratio, determined by availability
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## as per plr.R script ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# add weights to dataset, and expand so there is one row per patient per period of follow up
plr_data <- 
  tibble(
    patient_id = data_cohort$patient_id,
    ps = weightit_obj$ps,
    weight = weightit_obj$weights, 
  ) |>
  left_join(data_cohort, by="patient_id") |>
  uncount(event_time, .remove=FALSE) %>%
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = (time==event_time) & (event_indicator==TRUE),
  )

# fit a PLR using weights, but with no way to incorporate weight uncertainty
plr_model <- glm(
  formula_time_treatment, 
  family = binomial(), 
  weight = weight,
  data = plr_data,
  model = TRUE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
  x = FALSE,
  y = FALSE
)

# get the robust sandwich vcov matrix
plr_vcov <- get_vcov(plr_model, vcov = ~patient_id) # equivalent to vcovCL(x = plr_model, cluster = plr_data$patient_id, type = "HC0")

## cluster-robust variance estimation
# use newdata with predictions
# newdata argument with only treatment and time variables only works because we've already applied the weighting! This _will not_ work if there is any additional adjustment in the outcome model
# this is equivalent to using "avg_predictions" function without newdata argument, with marginalisation over the entire population
# if post-ipw adjustment is applied then remove "newdata" argument and use "avg_predictions" function instead
data_predictions <-
  predictions(
    plr_model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      time = seq_len(maxfup)
    ),
    type = "response",
    by = c("treatment", "time"),
    vcov = ~patient_id,
    #transform = \(x) {cumprod(1-x)} unfortunately this doesn't provide standard errors
  ) |> 
  group_by(treatment) %>%
  mutate(
    inc = estimate,
    inc.se = std.error,
    inc.low = conf.low,
    inc.high = conf.high, 
    ones = 1L,
    surv = cumprod(1-inc),
    surv.se = sqrt(cmlinc_variance(model = plr_model, vcov = plr_vcov, newdata = tibble(treatment=treatment, time=time), id = ones, time = time, weights = ones)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),
    
    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
  )


# compare estimates across products

# similarly to above, newdata argument will not work if doing any post-ipw adjustment (ie, adjustment in the outcome model)
data_comparisons_irr <-
  comparisons(
    plr_model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      time = seq_len(maxfup)
    ),
    variables = list(
      treatment =  c(0L, 1L)
    ),
    comparison = "ratio",
    by = "time",
    vcov = ~patient_id
  ) |>
  transmute(
    time, 
    irr = estimate,
    irr.se = std.error,
    irr.ll = conf.low,
    irr.ul = conf.high
  ) 

data_comparisons_irr.ln <-
  comparisons(
    plr_model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      time = seq_len(maxfup)
    ),
    variables = list(
      treatment =  c(0L, 1L)
    ),
    comparison = "lnratio",
    by = "time",
    vcov = ~patient_id
  ) |>
  transmute(
    time, 
    irr.ln = estimate,
    irr.ln.se = std.error, # equivalent to sqrt(((inc.se_1 / inc_1)^2) + ((inc.se_0 / inc_0)^2))
    irr.ll2 = exp(conf.low),
    irr.ul2 = exp(conf.high)
  ) 

data_comparisons_ird <-
  comparisons(
    plr_model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      time = seq_len(maxfup)
    ),
    variables = list(
      treatment =  c(0L, 1L)
    ),
    comparison = "difference",
    by = "time",
    vcov = ~patient_id
  ) |>
  transmute(
    time, 
    ird = estimate,
    ird.se = std.error,
    ird.ll = conf.low,
    ird.ul = conf.high
  ) 
  
data_comparisons_noncumulative <- 
  left_join(
    data_comparisons_irr,
    data_comparisons_irr.ln,
    data_comparisons_ird,
    by = "time"
  )

data_comparisons_cumulative <-
  data_predictions |> 
  pivot_wider(
    id_cols = all_of("time"),
    names_from = treatment,
    values_from = c(
      inc, inc.se, surv, surv.se, cmlinc, cmlinc.se,
    )
  ) |>
  mutate(
    
    # survival ratio, standard error, and confidence limits
    sr = (1-cmlinc_1) / (1-cmlinc_0),
    sr.ln.se = (cmlinc.se_0 / (1-cmlinc_0)) + (cmlinc.se_1 / (1-cmlinc_1)),
    sr.ll = exp(log(sr) + qnorm(0.025) * sr.ln.se),
    sr.ul = exp(log(sr) + qnorm(0.975) * sr.ln.se),
    
    # risk ratio, standard error, and confidence limits, using delta method
    rr = cmlinc_1 / cmlinc_0,
    # cirr.ln = log(cirr),
    rr.ln.se = sqrt((cmlinc.se_1 / cmlinc_1)^2 + (cmlinc.se_0 / cmlinc_0)^2),
    rr.ll = exp(log(rr) + qnorm(0.025) * rr.ln.se),
    rr.ul = exp(log(rr) + qnorm(0.975) * rr.ln.se),
    
    # risk difference, standard error and confidence limits, using delta method
    rd = cmlinc_1 - cmlinc_0,
    rd.se = sqrt((cmlinc.se_0^2) + (cmlinc.se_1^2)),
    rd.ll = rd + qnorm(0.025) * rd.se,
    rd.ul = rd + qnorm(0.975) * rd.se,
  )
  
data_comparisons <- 
  left_join(
    data_comparisons_noncumulative,
    data_comparisons_cumulative,
    by = "time"
  )


data_comparisons |>
  ggplot()+
  geom_line(aes(x=time, y=irr))+
  geom_ribbon(aes(x=time, ymin=irr.ll, ymax=irr.ul), alpha=0.2, fill="lightblue")+
  geom_ribbon(aes(x=time, ymin=irr.ll2, ymax=irr.ul2), alpha=0.2, fill="gold")+
  geom_hline(aes(yintercept=1), linetype='dashed')+
  theme_bw()

data_comparisons |>
  ggplot()+
  geom_line(aes(x=time, y=rr))+
  geom_ribbon(aes(x=time, ymin=rr.ll, ymax=rr.ul), alpha=0.2, fill="lightblue")+
  geom_hline(aes(yintercept=1), linetype='dashed')+
  theme_bw()

data_comparisons |>
  ggplot()+
  geom_line(aes(x=time, y=rd))+
  geom_ribbon(aes(x=time, ymin=rd.ll, ymax=rd.ul), alpha=0.2, fill="lightblue")+
  geom_hline(aes(yintercept=0), linetype='dashed')+
  theme_bw() 

# demonstrate equivalence using g-computaiton explicitly
data_gcomp <- 
  bind_rows(
    plr_data %>% mutate(treatment=1),
    plr_data %>% mutate(treatment=0)
  ) %>%
  mutate(
    prob = predict(plr_model, newdata=., type="response"),
  ) 

data_comparisons_gcomp <- 
  data_gcomp |>
  group_by(treatment, time) |>
  summarise(
    inc = weighted.mean(prob, rep(1L, length(prob))),
    ones = 1L
  ) |>
  group_by(treatment) %>%
  mutate(
    survival = cumprod(1-inc),
    #survivalse = sqrt(cmlinc_variance(plr_model, plr_vcov, ., patient_id, time, ones))
    #survival.ll = survival+qnorm(0.025)*survivalse,
    #survival.ul = survival+qnorm(0.975)*survivalse,
  )

