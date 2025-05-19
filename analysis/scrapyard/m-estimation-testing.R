
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

# outcome model formula, assuming balancing weights are applied
formula_time_treatment <- event_indicator ~ treatment + ns(time, 4) + treatment:ns(time, 4)


# predict function for plr model using custom vcov matrix
predict.glm.custom.vcov <- function(x, vcov, newdata){
  if(missing(newdata)){ newdata <- x$model }
  tt <- terms(x)
  Terms <- delete.response(tt)
  m.mat <- model.matrix(Terms,data=newdata)
  m.coef <- x$coef
  fit <- as.vector(m.mat %*% x$coef)
  se.fit <- sqrt(diag(m.mat%*%vcov%*%t(m.mat)))
  return(list(fit=fit,se.fit=se.fit))
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Approach 1:
# fit a standard propensity score with weightit
# expand dataset to persontime
# fit a PLR using weights but with no way to incoporate weight uncertainty
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data1 <- data_cohort

weightit1 <- 
  weightit(
    formula = weighting_formula,
    data = data1,
    method = "glm", 
    estimand = "ATE",
    #stabilize = TRUE
  )

# add weights to original dataset
data_weights1 <- 
  tibble(
    patient_id = data1$patient_id,
    ps = weightit1$ps,
    weight = weightit1$weights, 
  ) |>
  left_join(data1, by="patient_id") 


# expand dataset so there is one row per patient per period of follow up
data_model1 <- 
  data_weights1 |>
  uncount(event_time, .remove=FALSE) %>%
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = (time==event_time) & (event_indicator==TRUE),
  )

# fit outcome model 
model1 <- glm(
  formula_time_treatment, 
  family = binomial(), 
  weight = weight,
  data = data_model1,
  model = TRUE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
  x = FALSE,
  y = FALSE
)


# doesn't work because
# "Error in model.frame.default(formula = formula_time_treatment, data = data_model1,  : variable lengths differ (found for '(weights)')
# model1b <- glm_weightit(
#   formula = formula_time_treatment,
#   family = binomial(),
#   weightit = weightit1,
#   data = data_model1,
#   cluster = "patient_id",
#   model = FALSE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
#   x = FALSE,
#   y = FALSE
# )

# robust vcov matrix
vcov1 <- vcovCL(x = model1, cluster = data_model1$patient_id, type = "HC0")

cmlinc_data1 <-
  expand_grid(
    treatment =  c(0L, 1L),
    time = seq_len(maxfup)
  ) %>%
  mutate(
    # this uses the ipw.model to get the estimated incidence at each time point for each treatment, assuming the entire population received treatment A
    # it works correctly for the ATE because of the weights (ie as if setting treatment=1 or treatment=0 for entire population)
    inc = predict(model1, newdata = ., type="response"),
    inc.se = predict(model1, newdata = ., type="response", se.fit=TRUE)$se.fit, # this does not use vcov from vcovCL, so not cluster-robust
    inc.logit.se = predict.glm.custom.vcov(model1, vcov = vcov1, newdata = .)$se.fit, # cluster robust, but on the linear scale, not response scale
    inc.low = plogis(qlogis(inc) + (qnorm(0.025) * inc.logit.se)),
    inc.high = plogis(qlogis(inc) + (qnorm(0.975) * inc.logit.se)),
  ) |>
  group_by(treatment) %>%
  mutate(
    
    # standard errors for survival / cumulative incidence estimates are not M-estimable (afaik), so using delta method for the final step here
    
    ones = 1L,
    surv = cumprod(1-inc),
    surv.se = sqrt(cmlinc_variance(model = model1, vcov = vcov1, newdata = tibble(treatment=treatment, time=time), id = ones, time = time, weights = ones)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),
    
    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
  )


### using predictions function from marginaleffects
## standard variance estimation
predictions1_a <-
  avg_predictions(
    model1,
    type = "response",
    by = c("treatment", "time"),
  )
## cluster-robust variance estimation (using sandwich separately from marginaleffects)
predictions1_b <-
  avg_predictions(
    model1,
    type = "response",
    by = c("treatment", "time"),
    vcov = vcov1
  )
## cluster-robust variance estimation (using in-built marginaleffects call to sandwich)
predictions1_c <-
  avg_predictions(
    model1,
    type = "response",
    by = c("treatment", "time"),
    vcov = "HC0" 
  )
predictions1_d <-
  avg_predictions(
    model1,
    type = "response",
    by = c("treatment", "time"),
    vcov = ~patient_id 
  )




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Approach 2:
# fit a propensity score with weightit using expanded person time
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data2 <- 
  data_cohort |>
  uncount(event_time, .remove=FALSE) %>%
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = (time==event_time) & (event_indicator==TRUE),
    persontime = rep(rle(patient_id)$lengths, rle(patient_id)$lengths),
    inv_persontime = 1/persontime
  )

# estimate weights using expanded dataset 
weightit2 <-
  weightit(
    formula = weighting_formula,
    data = data2,
    method = "glm",
    estimand = "ATE",
    #stabilize = TRUE,
    #s.weights = "inv_persontime", # need to use this to trick weightit to get weights as if using this non-expanded dataset, BUT it doesn't seem to work as required
  )

# combine weights and original data
data_model2 <- 
  tibble(
    patient_id = data2$patient_id,
    time = data2$time,
    ps = weightit2$ps,
    weight = weightit2$weights, 
    s.weight = weightit2$s.weights, 
  ) |>
  left_join(
    data2, 
    by=c("patient_id", "time")
  )

# fit outcome model
model2 <- glm_weightit(
  formula_time_treatment, 
  family = binomial(), 
  weight = weightit2,
  data = data2,
  vcov = "asymp",
  cluster = ~patient_id,
  model = TRUE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
  x = FALSE,
  y = FALSE
)

model2b <- glm(
  formula_time_treatment, 
  family = binomial(), 
  #weight = weight,
  data = data_model2,
  model = TRUE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
  x = FALSE,
  y = FALSE
)

# robust vcov matrix
vcov2 <- vcovCL(x = model2, cluster = data_model2$patient_id, type = "HC0")


cmlinc_data2 <-
  expand_grid(
    treatment =  c(0L, 1L),
    time = seq_len(maxfup)
  ) %>%
  mutate(
    # this uses the ipw.model to get the estimated incidence at each time point for each treatment, assuming the entire population received treatment A
    # it works correctly for the ATE because of the weights (ie as if setting treatment=1 or treatment=0 for entire population)
    inc = predict(model2, newdata = ., type="response"),
    #inc.se = predict(model2, newdata = ., type="response")$se.fit, # this does not use vcov from vcovCL, so not cluster-robust
    inc.logit.se = predict.glm.custom.vcov(model2, vcov = vcov2, newdata = .)$se.fit, # cluster robust, but on the linear scale, not response scale
    inc.low = plogis(qlogis(inc) + (qnorm(0.025) * inc.logit.se)),
    inc.high = plogis(qlogis(inc) + (qnorm(0.975) * inc.logit.se)),
  ) |>
  group_by(treatment) %>%
  mutate(
    
    # standard errors for survival / cumulative incidence estimates are not M-estimable (afaik), so using delta method for the final step here
    
    ones = 1L,
    surv = cumprod(1-inc),
    surv.se = sqrt(cmlinc_variance(model = model2, vcov = vcov1, newdata = tibble(treatment=treatment, time=time), id = ones, time = time, weights = ones)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),
    
    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
  )



### using predictions function from marginaleffects
## standard variance estimation
predictions2_a <-
  avg_predictions(
    model2,
    type = "response",
    by = c("treatment", "time"),
  )
## cluster-robust variance estimation (using sandwich separately from marginaleffects)
predictions2_b <-
  avg_predictions(
    model2,
    type = "response",
    by = c("treatment", "time"),
    vcov = vcov2 # using cluster robust variance from 
  )
## cluster-robust variance estimation (using in-built marginaleffects call to sandwich)
predictions2_c <-
  avg_predictions(
    model2,
    type = "response",
    by = c("treatment", "time"),
    vcov = "HC0" # using cluster robust
  )





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Approach 3:
# fit a propensity score with weightit using expanded person time using weightitMSM
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# re-estimate weights using expanded dataset and MSM 
weightit3 <-
  weightitMSM(
    formula = list(
      weighting_formula
    ),
    data = data_persontime,
    method = "glm",
    estimand = "ATE",
    stabilize = TRUE,
    s.weights = "inv_persontime" # need to use this to trick weightitMSM to get weights as if using this non-expanded dataset
  )

data_weights_msm <- 
  tibble(
    patient_id = data_persontime$patient_id,
    time = data_persontime$time,
    ps = weightit_msm$ps,
    weightMSM = weightit_msm$weights, 
    s.weightMSM = weightit_msm$s.weights, 
  ) |>
  left_join(
    data_persontime, 
    by=c("patient_id", "time")
  )





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# fit model using pooled logistic regression without any explicit use of MSM functions
# and therefore cannot get M-estimation for time-period-specific incidence estimates
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




### using M-estimation for uncertainty
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# fit model using pooled logistic regression WITH explicit use of MSM functions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

msm_model <- 
  glm_weightit(
    formula_time_treatment,
    family = binomial(),
    weightit = weightit_msm,
    data = data_persontime,
    vcov = "asympt",# "asympt" = m-estimation
    cluster = "patient_id"
  )


data_gcomp <- 
  bind_rows(
    data_plr %>% mutate(treatment=1),
    data_plr %>% mutate(treatment=0)
  )

model_comparisons1 <-
  comparisons(
    msm_model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      time = c(0,seq_len(365))
    ),
    variables = list(
      treatment =  c(0L, 1L)
    )
  )


model_comparisons2 <-
  comparisons(
    msm_model,
    variables = list(
      treatment =  c(0L, 1L)
    )
  )

model_comparisons3 <-
  avg_predictions(
    msm_model,
    newdata = data_gcomp,
    variables = list(
      treatment =  c(0L, 1L)
      #time = seq_len(365)
    ),
    by="time"
  )
