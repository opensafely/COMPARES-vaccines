
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Demonstrate in a single script the basics of the study design
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library("tidyverse")
library("here")
library("glue")
library("arrow")
library("survival")
library("MatchIt")
library("WeightIt")
library("lmw")
library("cobalt")
library("splines")
library("sandwich")
library("marginaleffects")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

# simulation some example data

rcat <- function(n, levels, p, factor = TRUE) {
  x <- sample(x = levels, size = n, replace = TRUE, prob = p)
  if (factor) x <- factor(x, levels = levels)
  x
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# create simulated dataset with N patients
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set.seed <- 42

N <- 1000

maxfup <- 7 * 26 # 26 weeks

data_cohort <-
  tibble(
    patient_id = seq_len(N),
    age = as.integer(rnorm(N, 60, 10)),
    sex = as.factor(rcat(N, c("female", "male"), p = c(0.5, 0.5))),
    vax_date = runif(N, 0, 100),
    vax_product = rcat(N, c("A", "B"), p = c(0.5, 0.5)),
    treatment = (vax_product == "B") * 1L,
    outcome_date = vax_date + runif(N, 1, 1000),
    censor_date = vax_date + runif(N, 1, 1000),
    event_time = as.integer((pmin(outcome_date, censor_date, 365) - vax_date) + 1),
    event_indicator = outcome_date <= pmin(censor_date, maxfup),
    invfup1 = 1 / event_time
  )

## preliminary model inputs

# formula for maching or predicting treatment probability

matching_vars_exact <- c("sex")
matching_vars_caliper <- c(vax_date = 7, age = 3)

weighting_formula <- treatment ~ vax_date + age + sex


lmw_formula <-  ~ treatment + vax_date + age + sex

# outcome model formula for PLR, assuming balancing weights are applied
formula_time_treatment <- event_indicator ~ treatment + ns(time, 4) + treatment:ns(time, 4)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# as per match section in balance.R script
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# perform 1-1 matching without replacement, using matchit

obj_matchit <-
  matchit(
    formula = treatment ~ 1,
    data = data_cohort,
    method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
    replace = FALSE,
    estimand = "ATT", # since we are doing exact matching, ATT is equivalent to ATU. although we'll actually get the ATO (average treatment in the overlap)
    exact = matching_vars_exact,
    caliper = matching_vars_caliper, std.caliper = FALSE,
    m.order = "data", # data is sorted on (effectively random) patient ID
    # verbose = TRUE,
    ratio = 1L # could also consider exact matching only, with n:m ratio, determined by availability
  )

weights_matchit <-
  tibble(
    patient_id = data_cohort$patient_id,
    treatment = obj_matchit$treat,
    weight = obj_matchit$weights,
    ps = (treatment / weight) + ((1 - treatment) * (1 - (1 / weight)))
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# as per weight section in balance.R script
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# fit a standard propensity score logistic regression model with weightit
obj_weightit <-
  weightit(
    formula = weighting_formula,
    data = data_cohort,
    method = "glm",
    estimand = "ATE",
    # stabilize = TRUE
  )

weights_weightit <-
  tibble(
    patient_id = data_cohort$patient_id,
    treatment = obj_weightit$treat,
    weight = obj_weightit$weights,
    ps = obj_weightit$ps,
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## as per plr.R script ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# add weights to dataset, and expand so there is one row per patient per period of follow up
plr_data <-
  weights_weightit |> # use weightit weights for now
  select(-treatment) |>
  left_join(data_cohort, by = "patient_id") |>
  uncount(event_time, .remove = FALSE) |>
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = (time == event_time) & (event_indicator == TRUE),
  ) |>
  group_by(patient_id) |>
  mutate(
    fup = n(),
    invfup = 1 / fup
  ) |>
  ungroup()

# fit a PLR using weights, but with no way to incorporate weight uncertainty
# NOTE: We cannot use WeightIt::glm_weightit here to use M-estimation to get valid SEs by passing in the weightit object directly,
# because we have expanded the dataset to fit a PLR.
# so the m.parts are not appropriate
plr_model <- glm(
  formula_time_treatment,
  family = binomial(),
  weights = plr_data$weight,
  data = plr_data,
  model = TRUE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
  x = FALSE,
  y = FALSE
)
summary(plr_model)



# get the robust sandwich vcov matrix
plr_vcov <- get_vcov(plr_model, vcov = ~patient_id) # equivalent to vcovCL(x = plr_model, cluster = plr_data$patient_id, type = "HC0")

# get counterfactual predictions
plr_cmlinc <-
  expand_grid(
    treatment =  c(0L, 1L),
    time = c(0, seq_len(maxfup))
  ) %>%
  mutate(
    # this uses the ipw.model to get the estimated incidence at each time point for each treatment, assuming the entire population received treatment A
    # it works correctly for the ATE because of the weights (ie as if setting treatment=1 or treatment=0 for entire population)
    inc = predict(plr_model, newdata = ., type = "response"),
    inc.se = predict(plr_model, newdata = ., type = "response", se.fit = TRUE)$se.fit, # this does not use vcov from vcovCL, so not cluster-robust
    inc.logit.se = predict.glm.custom.vcov(plr_model, vcov = plr_vcov, newdata = .)$se.fit, # cluster robust, but on the linear scale, not response scale
    inc.low = plogis(qlogis(inc) + (qnorm(0.025) * inc.logit.se)),
    inc.high = plogis(qlogis(inc) + (qnorm(0.975) * inc.logit.se)),
  ) |>
  group_by(treatment) |>
  mutate(
    dummy_id_weight = 1L,
    surv = cumprod(1 - inc),
    surv.se = sqrt(cmlinc_variance(model = plr_model, vcov = plr_vcov, newdata = tibble(treatment = treatment, time = time), id = dummy_id_weight, time = time, weights = dummy_id_weight)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),

    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
    # rmst = cumsum(surv),
    # rmst.se = sqrt(((2* cumsum(time*surv)) - (rmst^2))/n.risk), # this only works if one row per day using fill_times! otherwise need sqrt(((2* cumsum(time*interval*surv)) - (rmst^2))/n.risk)
    # rmst.low = rmst + (qnorm(0.025) * rmst.se),
    # rmst.high = rmst + (qnorm(0.975) * rmst.se),
  )

# contrast counterfactual predictions
plr_constracts <-
  # see following link for canonical-ish place where these are defined https://github.com/opensafely-actions/kaplan-meier-function/blob/main/analysis/km.R#L540
  plr_cmlinc |>
  pivot_wider(
    id_cols = "time",
    names_from = treatment,
    values_from = c(
      cmlinc, cmlinc.se, cmlinc.low, cmlinc.high,
    )
  ) |>
  mutate(
    # survival ratio, standard error, and confidence limits
    sr = (1 - cmlinc_1) / (1 - cmlinc_0),
    sr.ln.se = (cmlinc.se_0 / (1 - cmlinc_0)) + (cmlinc.se_1 / (1 - cmlinc_1)),
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
  ) |>
  select(
    # remove rows relating to individual curves
    -ends_with("0"),
    -ends_with("1"),
    -ends_with(".se"),
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# as per lmw section in balance.R script
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# obtain weights implied by outcome regression

obj_lmw <-
  lmw(
    formula = lmw_formula,
    data = data_cohort,
    treat = "treatment",
    estimand = "ATE",
    method = "MRI", # MRI gets weights as if there were a separate outcome model for each treatment group
  )

weights_lmw <-
  tibble(
    patient_id = data_cohort$patient_id,
    treatment = data_cohort$treatment,
    weight = obj_lmw$weights,
    ps = (treatment / weight) + ((1 - treatment) * (1 - (1 / weight)))
    # something about weights here is wrong because PS are not between [0,1].
  )
