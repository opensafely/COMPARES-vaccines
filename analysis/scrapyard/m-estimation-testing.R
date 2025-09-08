
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Demonstrate in a single script the basics of the study design
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library("tidyverse")
library("survival")
library("WeightIt")
library("splines")
library("marginaleffects")
library("sandwich")

## Import custom user functions from lib
source(here::here("analysis", "0-lib", "utility.R"))

## Import design elements
# source(here("analysis", "0-lib", "design.R"))

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

N <- 500

maxfup <- 365

data_cohort <-
  tibble(
    patient_id = seq_len(N),
    age = as.integer(rnorm(N, 60, 10)),
    sex = as.factor(rcat(N, c("female", "male"), p = c(0.5, 0.5))),
    vax_date = runif(N, 0, 100),
    probvaxA = 0.4 + 0.3 * (sex == "male"),
    vax_product = if_else(runif(N) < probvaxA, "A", "B"),
    treatment = (vax_product == "B") * 1L,
    outcome_date = vax_date + runif(N, 1, 1000 - 100 * (sex == "male") - 100 * (patient_id %% 2) - 5 * (age - 60) + 200 * (vax_product == "A")),
    censor_date = vax_date + runif(N, 1, 1000),
    event_time = as.integer((pmin(outcome_date, censor_date, 365) - vax_date) + 1),
    event_indicator = outcome_date <= pmin(censor_date, maxfup),
  ) |>
  mutate(
    # event_time = 300 # OVERWRITE EVENT TIME TEMPORARILY TO TEST
  )

## preliminary model inputs

# formula for predicting treatment probability
weighting_formula <- treatment ~ vax_date + age + sex

# outcome model formula, assuming balancing weights are applied
formula_time_treatment <- event_indicator ~ treatment + ns(time, 4) + treatment:ns(time, 4)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Approach 1:
# fit a standard propensity score with weightit, on a one-row-per-person dataset
# then use expanded dataset to incorporate person-time, with fixed IP weight at each time point
# fit a PLR using the IP weights, but there is no way to incorporate weight uncertainty, weightit doesn't have the correct "mparts"

# this is the approach currently implemented in the protocol
# this is what would happen if we treated the whole situation as a marginal structural model fit using pooled logistic regression,
# with people "at risk" of vaccination only at t=0, and subsequent partial weights at t>0 all set to 1
# then the weights for each person are `cumsum(c(weight_t0, rep(1, max_follow_up - 1))`, which is just `rep(weight_t0, max_follow_up)`

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

data_expanded <-
  data_cohort |>
  uncount(event_time, .remove = FALSE) %>%
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = (time == event_time) & (event_indicator == TRUE),
    persontime = rep(rle(patient_id)$lengths, rle(patient_id)$lengths), # should be equal to event_time
    inv_persontime = 1 / persontime
  )

weightit1 <-
  weightit(
    data = data_cohort,
    formula = weighting_formula,
    method = "glm",
    estimand = "ATE",
    # stabilize = TRUE
  )

# fit outcome model, using person-level weights from unexpanded model
model1_naive <- glm(
  data = data_expanded,
  formula = formula_time_treatment,
  family = binomial(),
  weight = rep(weightit1$weights, data_cohort$event_time),
)

# doesn't work:
# "Error in model.frame.default(formula = formula_time_treatment, data = data_expanded,  : variable lengths differ (found for '(weights)')
# ultimately the weightit object doesn't have the right bits!
# model1_HC0 <- glm_weightit(
#   formula = formula_time_treatment,
#   family = binomial(),
#   weightit = weightit1,
#   data = data_expanded,
#   cluster = "patient_id"
# )

# robust vcov matrix
vcov1_naive <- vcov(model1_naive)
vcov1_HC0 <- vcovCL(x = model1_naive, cluster = data_expanded$patient_id, type = "HC0")


## variance estimation as per fitted model
predictions1_naive <-
  avg_predictions(
    model1_naive,
    type = "response",
    by = c("treatment", "time"),
  )

## variance estimationl using clust-robust HC0
predictions1_HC0 <-
  avg_predictions(
    model1_naive,
    type = "response",
    by = c("treatment", "time"),
    vcov = vcov1_HC0 # same as vcov = ~patient_id
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Approach 2:
# fit a propensity score with weightit using expanded person time, without using sampling weights
#
# This is the approach suggested by Noah, I think
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# estimate weights using expanded dataset,
# without reweighting based on persontime ( / expansion factor)
weightit2 <-
  weightit(
    data = data_expanded,
    formula = weighting_formula,
    method = "glm",
    estimand = "ATE"
  )

# fit outcome model on expanded dataset, using naive SEs
model2_naive <- glm(
  data = data_expanded,
  formula_time_treatment,
  family = binomial(),
  weights = weightit2$weights,
)

# fit outcome model on expanded dataset, using HC0 patient clusters for SEs
model2_HC0 <- glm_weightit(
  data = data_expanded,
  formula_time_treatment,
  family = binomial(),
  weightit = weightit2,
  vcov = "HC0",
  cluster = ~patient_id
)

# fit outcome model on expanded dataset, using M-estimation for SEs
model2_asympt <- glm_weightit(
  data = data_expanded,
  formula_time_treatment,
  family = binomial(),
  weightit = weightit2,
  vcov = "asympt",
  cluster = ~patient_id
)

vcov2_naive <- vcov(model2_naive)
vcov2_HC0 <- vcov(model2_HC0)
vcov2_asympt <- vcov(model2_asympt)



# get estimate for each treatment at each time point
# NOTE: we could use predictions() here instead of avg_predictions(),
# because we've got the weights and there is no post-weighting outcome model
# so we don't need to average over all patients, but simply predict the outcome at each time point for teach treatment
# this is much quicker that predicting outcomes for all patients.
# for example:
# predictions(
#   model3_asympt,
#   newdata = expand_grid(
#     treatment =  c(0L, 1L),
#     time = seq_len(maxfup)
#   ),
#   type = "response",
#   by = c("treatment", "time"),
#   # transform = \(x) {cumprod(1-x)} unfortunately this doesn't provide standard errors we need
# )



## variance estimation with naive SEs
predictions2_naive <-
  avg_predictions(
    model2_naive,
    type = "response",
    by = c("treatment", "time"),
  )

## variance estimation with HC0 SEs
predictions2_HC0 <-
  avg_predictions(
    model2_HC0,
    type = "response",
    by = c("treatment", "time"),
  )

## variance estimation with M-estimation SEs
predictions2_asympt <-
  avg_predictions(
    model2_asympt,
    type = "response",
    by = c("treatment", "time"),
  )



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Approach 3:
# fit a propensity score with weightit using expanded person time, and using sampling weights to recover size of collapsed dataset

# this is an approach which recovers the original per-person weights in approach 1, but results in very different estimates
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# estimate weights using expanded dataset, using sampling weights to recover size of data without expansion
weightit3 <-
  weightit(
    data = data_expanded,
    formula = weighting_formula,
    method = "glm",
    estimand = "ATE",
    s.weights = "inv_persontime" # use this to trick weightit to get weights as if using this non-expanded dataset, BUT it doesn't seem to work as required because the
  )

# overwrite s.weights in weighit object because otherwise the s.weights are re-used in glm_weightit
# we don't want them to be reused becuase they were only there to trick weightit into giving equal weight to people with different time-to-events
# if we don't do this, then the effect estimates are completely different to approach 1 and 2, so something is clearly wrong
weightit3$s.weights <- rep(1, length(weightit3$weights))

# alternatively, use a custom weightit object, by adding the weights manually.
# but this doesn't come with the mparts (see keep.mparts argument in weightit), so M-estimation not possible
weightit3_custom <-
  as.weightit(
    x = rep(weightit1$weights, data_cohort$event_time),
    treat = data_expanded$treatment,
    estimand = "AVE",
    s.weights = rep(1, length(weightit3$weights))
  )

# these are the same!
summary(weightit3)
summary(weightit3_custom)

# fit outcome model on expanded dataset, using naive SEs
model3_naive <- glm(
  data = data_expanded,
  formula_time_treatment,
  family = binomial(),
  weights = weightit3$weights,
)

# fit outcome model on expanded dataset, using HC0 patient clusters for SEs
model3_HC0 <- glm_weightit(
  data = data_expanded,
  formula_time_treatment,
  family = binomial(),
  weightit = weightit3,
  vcov = "HC0",
  cluster = ~patient_id
)
# fit outcome model on expanded dataset, using M-estimation for SEs
model3_asympt <- glm_weightit(
  data = data_expanded,
  formula_time_treatment,
  family = binomial(),
  weightit = weightit3,
  vcov = "asympt",
  cluster = ~patient_id
)

vcov3_naive <- vcov(model3_naive)
vcov3_HC0 <- vcov(model3_HC0)
vcov3_asympt <- vcov(model3_asympt)


## variance estimation using naive SEs
predictions3_naive <-
  avg_predictions(
    model3_naive,
    type = "response",
    by = c("treatment", "time"),
  )

## variance estimation using cluster HC0 SEs
predictions3_HC0 <-
  avg_predictions(
    model3_HC0,
    type = "response",
    by = c("treatment", "time"),
  )

## variance estimation using m-estimation SEs
predictions3_asympt <-
  avg_predictions(
    model3_asympt,
    type = "response",
    by = c("treatment", "time"),
  )




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## compare all approaches

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

weightit1$weights |> head(10)
weightit2$weights[cumsum(data_cohort$event_time)] |> head(10)
weightit3$weights[cumsum(data_cohort$event_time)] |> head(10)

all_predictions <-
  bind_rows(
    predictions1_naive |> mutate(approach = "1", SEs = "naive"),
    predictions1_HC0 |> mutate(approach = "1", SEs = "HC0"),
    # predictions1_asympt |> mutate(approach = "1", SEs= "asympt"),
    predictions2_naive |> mutate(approach = "2", SEs = "naive"),
    predictions2_HC0 |> mutate(approach = "2", SEs = "HC0"),
    predictions2_asympt |> mutate(approach = "2", SEs = "asympt"),
    predictions3_naive |> mutate(approach = "3", SEs = "naive"),
    predictions3_HC0 |> mutate(approach = "3", SEs = "HC0"),
    predictions3_asympt |> mutate(approach = "3", SEs = "asympt"),
  )

# as expected, predictions are all the same no matter which SE we choose
# but approach 2 is different from 1 and 3
ggplot(all_predictions) +
  geom_line(aes(x = time, y = estimate, group = paste(approach, treatment), colour = approach, linetype = approach)) +
  facet_grid(cols = vars(SEs), rows = vars(treatment)) +
  theme_bw()

# SEs are very similar (not the same) between approaches 2 and 3, but naive SEs are noticably smaller
ggplot(all_predictions) +
  geom_line(aes(x = time, y = std.error, group = paste(approach, treatment), colour = approach, linetype = approach)) +
  facet_grid(cols = vars(SEs), rows = vars(treatment)) +
  theme_bw()




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
### now do cumulative incidence curves using manual delta method sweeping through time
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cml_inc3_asympt <-
  predictions3_asympt |>
  group_by(treatment) %>%
  mutate(
    surv = cumprod(1 - estimate),
    # use custom delta method to get SEs for cumulative incidence - ./analysis/0-lib/utility.R script
    surv.se = sqrt(cmlinc_variance(model = model3_asympt, vcov = vcov(model3_asympt), newdata = tibble(treatment = treatment, time = time), time = time)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),
    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
  ) |>
  mutate(
    treatment_chr = as.character(treatment)
  )

ggplot(cml_inc3_asympt) +
  geom_line(aes(x = time, y = cmlinc, group = treatment_chr, colour = treatment_chr)) +
  geom_ribbon(aes(x = time, ymin = cmlinc.low, ymax = cmlinc.high, group = treatment_chr, fill = treatment_chr), colour = "transparent", alpha = 0.2) +
  theme_bw()



cml_inc3_naive <-
  predictions3_naive |>
  group_by(treatment) %>%
  mutate(
    surv = cumprod(1 - estimate),
    # use custom delta method to get SEs for cumulative incidence - ./analysis/0-lib/utility.R script
    surv.se = sqrt(cmlinc_variance(model = model3_naive, vcov = vcov3_naive, newdata = tibble(treatment = treatment, time = time), time = time)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),

    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
  ) |>
  mutate(
    treatment_chr = as.character(treatment)
  )

ggplot(cml_inc3_naive) +
  geom_line(aes(x = time, y = cmlinc, group = treatment_chr, colour = treatment_chr)) +
  geom_ribbon(aes(x = time, ymin = cmlinc.low, ymax = cmlinc.high, group = treatment_chr, fill = treatment_chr), colour = "transparent", alpha = 0.2) +
  theme_bw()
