
data_persontime <-
  data_all  |>
  uncount(event_time, .remove = FALSE) |> # may change this from event_time to final censoring time in case we want keep dead people under follow up
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = ((time == event_time) & (event_indicator == TRUE))*1L
  )

# create control spec if using parglm for model fitting
parglm_control <- parglm.control(maxit = 40, nthreads = 2)

# create dataset to link estimates onto
estimates_scaffold <- expand_grid(
  treatment =  c(0L, 1L),
  time = c(0, seq_len(maxfup))
)

subgroup_data <-
  data_persontime |>
  mutate(run_model = .data[[subgroup]] %in% subgroups_both_treatments_with_min_events) |>
  filter(
    ageband == "all"
  )

estimates_scaffold <- expand_grid(
  treatment =  c(0L, 1L),
  time = c(0, seq_len(maxfup))
)

run_model_flag <- TRUE

model_fit_function <- function(plrdata, run_model_flag) {
  tryCatch(
    if (run_model_flag) {
      glm(
        formula_time_treatment,
        family = binomial(),
        weights = weight,
        data = plrdata,
        method = fastglm::fastglm_fit, # use this if using fastglm -  don't use fastglm for now as v0.1.2 doesn't work well with predict function (v0.1.4 does, so wait for update)
        #control = parglm_control, # use this if using parglm 
        model = TRUE#, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
        #x = FALSE
      )
    } else {
      stop("Not enough events")
    },
    error = function(e) {
      cat("glm error:", conditionMessage(e), "\n")
      NULL
    }
  )
}

model_status_function <- function(model, run_model_flag){
  if (!run_model_flag) {
    "Model not run"
  } else if (run_model_flag & !model$converged) {
    "Model not converged"
  } else if (run_model_flag & model$converged) {
    "Model converged"
  }
}


model_estimates_function <- function(model, subgroup_data){
  
    # sandwich::vcovCL doesn't handle formulae properly! hence inclusion of "model=TRUE" above - be careful
    vcov <- vcovCL(x = model, cluster = subgroup_data$patient_id, type = "HC0") # or use `marginaleffects::get_vcov(model, vcov = ~patient_id)`
    #vcov <- marginaleffects::get_vcov(model, vcov = ~patient_id)
    
    estimates_scaffold %>%
      mutate(
        # this uses the ipw.model to get the estimated incidence at each time point for each treatment, assuming the entire population received treatment A
        # it works correctly for the ATE because of the weights (ie as if setting treatment=1 or treatment=0 for entire population)
        inc = predict(
          model,
          newdata = estimates_scaffold,
          type = "response"
        ),
        # inc.se = predict(
        #   model,
        #   newdata = estimates_scaffold,
        #   type = "response",
        #   se.fit = TRUE
        # )$se.fit, # this does not use vcov from vcovCL, so not cluster-robust
        inc.logit.se = predict.glm.custom.vcov(
          model,
          vcov = vcov,
          newdata = estimates_scaffold
        )$se.fit, # cluster robust, but on the linear scale, not response scale
        inc.low = plogis(qlogis(inc) + (qnorm(0.025) * inc.logit.se)),
        inc.high = plogis(qlogis(inc) + (qnorm(0.975) * inc.logit.se)),
      ) |>
      group_by(treatment) %>%
      mutate(
        dummy_id_weight = 1L,
        surv = cumprod(1 - inc),
        surv.se = sqrt(cmlinc_variance(
          model = model,
          vcov = vcov,
          newdata = tibble(treatment = treatment, time = time),
          id = dummy_id_weight,
          time = time,
          weights = dummy_id_weight
        )),
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
      ) |>
      select(-dummy_id_weight)
}


subgroup_model <- model_function(subgroup_data, run_model_flag)

subgroup_model_status <- model_status_function(subgroup_model, run_model_flag)

subgroup_model_estimates <- model_estimates_function(subgroup_model, subgroup_data)

subgroup_model_vcov <- vcovCL(x = subgroup_model, cluster = subgroup_data$patient_id, type = "HC0")

predict(
  subgroup_model,
  newdata = estimates_scaffold,
  type = "response"
)

predict(
  subgroup_model,
  newdata = estimates_scaffold,
  type = "response",
  se.fit = TRUE
)

predict.glm.custom.vcov(
  subgroup_model,
  vcov = subgroup_model_vcov,
  newdata = estimates_scaffold
)$se.fit

predict.glm.custom.vcov(
  subgroup_model,
  vcov = subgroup_model_vcov,
  newdata = estimates_scaffold
)$fit

predict(
  subgroup_model,
  newdata = estimates_scaffold,
  type = "link"
)

predict.glm.custom.vcov <- function(x, vcov, newdata) {
  if (missing(newdata)) {
    newdata <- x$model
  }
  tt <- terms(x)
  Terms <- delete.response(tt)
  m.mat <- model.matrix(Terms, data = newdata)
  m.coef <- x$coef
  fit <- as.vector(m.mat %*% x$coef)
  se.fit <- sqrt(diag(m.mat %*% vcov %*% t(m.mat)))
  return(list(fit = fit, se.fit = se.fit))
}
