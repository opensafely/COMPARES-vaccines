############################
# testing using glm_weightit
###############################



############################
# testing using marginal effects for the constrast scripts
###############################

# Use with `the-whole-game.R` script to produce data, weights, and model

## cluster-robust variance estimation
# use newdata with predictions
# newdata argument with only treatment and time variables only works because we've already applied the weighting! This _will not_ work if there is any additional adjustment in the outcome model
# this is equivalent to using "avg_predictions" function without newdata argument, with marginalisation over the entire population
# if adjustment is applied _after_ the use of weights, then remove "newdata" argument and use "avg_predictions" function instead
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
    # transform = \(x) {cumprod(1-x)} unfortunately this doesn't provide standard errors
  ) |>
  group_by(treatment) %>%
  mutate(
    inc = estimate,
    inc.se = std.error,
    inc.low = conf.low,
    inc.high = conf.high,
    ones = 1L,
    surv = cumprod(1 - inc),
    surv.se = sqrt(cmlinc_variance(model = plr_model, vcov = plr_vcov, newdata = tibble(treatment = treatment, time = time), id = ones, time = time, weights = ones)),
    surv.low = surv + (qnorm(0.025) * surv.se),
    surv.high = surv + (qnorm(0.975) * surv.se),

    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
  )


# compare estimates across products

# similarly to above, this time using `marginaleffects::comparisons function with `comparison` argument used to define the contrast of interest.
# the newdata argument will not work if doing any adjustment after applying the weights (ie, adjustment in the outcome model)
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
  )

data_comparisons <-
  left_join(
    data_comparisons_noncumulative,
    data_comparisons_cumulative,
    by = "time"
  )


data_comparisons |>
  ggplot() +
  geom_line(aes(x = time, y = irr)) +
  geom_ribbon(aes(x = time, ymin = irr.ll, ymax = irr.ul), alpha = 0.2, fill = "lightblue") +
  geom_ribbon(aes(x = time, ymin = irr.ll2, ymax = irr.ul2), alpha = 0.2, fill = "gold") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  theme_bw()

data_comparisons |>
  ggplot() +
  geom_line(aes(x = time, y = rr)) +
  geom_ribbon(aes(x = time, ymin = rr.ll, ymax = rr.ul), alpha = 0.2, fill = "lightblue") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  theme_bw()

data_comparisons |>
  ggplot() +
  geom_line(aes(x = time, y = rd)) +
  geom_ribbon(aes(x = time, ymin = rd.ll, ymax = rd.ul), alpha = 0.2, fill = "lightblue") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_bw()

# demonstrate equivalence using g-computation explicitly
data_gcomp <-
  bind_rows(
    plr_data %>% mutate(treatment = 1),
    plr_data %>% mutate(treatment = 0)
  ) %>%
  mutate(
    prob = predict(plr_model, newdata = ., type = "response"),
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
    survival = cumprod(1 - inc),
    # survivalse = sqrt(cmlinc_variance(plr_model, plr_vcov, ., patient_id, time, ones))
    # survival.ll = survival+qnorm(0.025)*survivalse,
    # survival.ul = survival+qnorm(0.975)*survivalse,
  )



############################
# testing using glm_weightit
###############################

glm_weightit(
  formula = formula_time_treatment,
  weightit = obj_weightit,
  data = plr_data,
  cluster = ~patient_id
)
# Error in `glm_weightit()`: ! (from `glm()`) variable lengths differ (found for '(weights)')
# Basically, it expects the weightit object to have been built on data the same size as the data argument

obj_weightit <-
  weightit(
    formula = weighting_formula,
    data = data_cohort,
    method = "glm",
    estimand = "ATE",
    # stabilize = TRUE
  )

summary(obj_weightit)

obj_weightit_expanded <-
  weightit(
    formula = weighting_formula,
    data = plr_data,
    method = "glm",
    estimand = "ATE",
    s.weights = "invfup",
  )

summary(obj_weightit_expanded)

weights_weightit_expanded <-
  tibble(
    patient_id = plr_data$patient_id,
    treatment = obj_weightit_expanded$treat,
    weight = obj_weightit_expanded$weights,
    ps = obj_weightit_expanded$ps,
  ) |>
  group_by(patient_id) |>
  summarise(
    treatment = first(treatment),
    weight = first(weight),
    ps = first(ps)
  )

glm_obj <- glm_weightit(
  formula = formula_time_treatment,
  data = plr_data,
  weightit = obj_weightit_expanded,
  cluster = ~patient_id,
  vcov = "asympt"
)

summary(glm_obj)

# custom weightit
custom_weightit <-
  as.weightit(
    x = plr_data$weight,
    treat = plr_data$treatment,
    estimand = "AVE",
    s.weights = plr_data$invfup,
  )

summary(custom_weightit)
summary(obj_weightit_expanded)



obj_weightit3a <-
  weightit(
    formula = weighting_formula,
    data = data_cohort,
    method = "glm",
    estimand = "ATE",
    s.weights = rep(c(1, 0), each = length(data_cohort$patient_id) / 2)
    # stabilize = TRUE
  )


obj_weightit3b <-
  weightit(
    formula = weighting_formula,
    data = data_cohort |> filter(row_number() <= n() / 2),
    method = "glm",
    estimand = "ATE",
    # stabilize = TRUE
  )


summary(obj_weightit)
summary(obj_weightit1)
summary(obj_weightit2, ignore.s.weights = TRUE)

summary(obj_weightit3a)
summary(obj_weightit3b)
