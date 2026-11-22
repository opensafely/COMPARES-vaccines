# # # # # # # # # # # # # # # # # # # # #

# Purpose: collate all scripts to create testing environment for M-estimation
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library("optparse")
library("tidyverse")
library("here")
library("glue")
library("arrow")
library("WeightIt")
library("survival")
library("splines")
library("sandwich")
library("fastglm")

# library("marginaleffects")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "cv"
  method <- "match"
  spec <- "A"
  subgroup <- "ageband"
  outcome <- "covid_death"
} else {

  removeobjects <- TRUE

  option_list <- list(
    make_option("--cohort", type = "character"),
    make_option("--method", type = "character"),
    make_option("--spec", type = "character"),
    make_option("--subgroup", type = "character", default = "all"),
    make_option("--outcome", type = "character")
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
  list2env(opt, .GlobalEnv)

  # cohort <- args[[1]]
  # method <- args[[2]]
  # spec <- args[[3]]
  # subgroup <- args[[4]]
  # outcome <- args[[5]]
}


# arg symbols / quosures
subgroup_sym <- sym(subgroup)
outcome_sym <- sym(outcome)

# create output directories ----

output_dir <- here_glue("output", "4-contrast", cohort, "{method}-{spec}", subgroup, outcome, "plr")
fs::dir_create(output_dir)

## import unadjusted cohort data ----
# only needed if rerunning weighting model
# data_cohort <- read_feather(here("output", "2-select", cohort, "data_cohort.arrow"))

## import data event counts to define which models can be fitted
data_event_counts <- read_feather(here_glue("output", "3-adjust", cohort, "combine", "table_event_counts.arrow"))

subgroups_both_treatments_with_min_events <-
  data_event_counts |>
  filter(
    cohort == !!cohort,
    method == !!method,
    spec == !!spec,
    outcome == !!outcome,
    subgroup == !!subgroup,
    flag_subgroups_both_treatments_with_min_events
  ) |>
  distinct(subgroup_level) |>
  pull(subgroup_level)

## import weights from matching or weighting method ----
data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "combine", "data_weights.arrow"))
# data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "{method}-{spec}", "data_adjusted.arrow"))

data_all <-
  data_weights |>
  select(
    patient_id,
    vax_product,
    treatment,
    treatment_date,
    all_of(subgroup),
    all_of(paste0(c(outcome, "death", "dereg"), "_date")),
    weight = glue("wt_{cohort}_{method}_{spec}")
  ) |>
  mutate(
    event_date = as.Date(.data[[glue("{outcome}_date")]]),

    # person-time is up to and including censor date
    censor_date = pmin(
      dereg_date,
      death_date,
      study_dates$followupend_date,
      treatment_date + maxfup,
      na.rm = TRUE
    ),

    noncompetingcensor_date = pmin(
      dereg_date,
      study_dates$followupend_date,
      treatment_date + maxfup,
      na.rm = TRUE
    ),

    event_time = tte(treatment_date, event_date, censor_date, na.censor = FALSE),
    event_indicator = censor_indicator(event_date, censor_date),

    # possible competing events
    death_time = tte(treatment_date, death_date, noncompetingcensor_date, na.censor = FALSE),
    censor_time = tte(treatment_date, censor_date, censor_date, na.censor = FALSE),
  )

stopifnot("censoring dates must be non-missing" = all(!is.na(data_all$censor_date)))

stopifnot("origin dates must be non-missing" = all(!is.na(data_all$treatment_date)))

times_count <- table(cut(data_all$event_time, c(-Inf, 0, 1, Inf), right = FALSE, labels = c("<0", "0", ">0")), useNA = "ifany")
if (!identical(as.integer(times_count), c(0L, 0L, nrow(data_all)))) {
  print(times_count)
  stop("all event times must be strictly positive")
}

formula_time_treatment <- event_indicator ~ treatment + ns(time, 4) + treatment:ns(time, 4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# standard errors calculated using pooled logistic regression,
# but this does not account for uncertainty in the weights
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

data_persontime <-
  data_all  |>
  uncount(event_time, .remove = FALSE) |> # may change this from event_time to final censoring time in case we want keep dead people under follow up
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = ((time == event_time) & (event_indicator == TRUE))*1L
  )

# create dataset to link estimates onto
estimates_scaffold <- expand_grid(
  treatment =  c(0L, 1L),
  time = c(0, seq_len(maxfup))
)

# estimate time-specific incidence from using pooled logistic regression ----
# Subgroup models ------------------
subgroup_models <-
  data_persontime |>
  mutate(run_model = .data[[subgroup]] %in% subgroups_both_treatments_with_min_events) |>
  group_by(!!subgroup_sym, run_model) |>
  nest() |>
  mutate(
    model = map2(
      .x = data,
      .y = run_model,
      .f = function(plrdata, run_model) {
        tryCatch(
          if (run_model) {
            glm(
              formula_time_treatment,
              family = binomial(),
              weights = weight,
              data = plrdata,
              method = fastglm_fit,
              model = TRUE, # this needs to be true for vcovCL to work as needed - shame because it takes up a lot of memory
              x = FALSE,
              y = FALSE
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
    ),

    model_status = map2_chr(run_model, model, \(.x, .y){
      if (!.x) {
        "Model not run"
      } else if (.x & !.y$converged) {
        "Model not converged"
      } else if (.x & .y$converged) {
        "Model converged"
      }
    }),

    estimates = map2(
      .x = data,
      .y = model,
      .f = function(plrdata, model) {

        if (!is.null(model)) {
          # sandwich::vcovCL doesn't handle formulae properly! hence inclusion of "model=TRUE" above - be careful
          #vcov <- vcovCL(x = model, cluster = plrdata$patient_id, type = "HC0") # or use `marginaleffects::get_vcov(model, vcov = ~patient_id)`
          vcov <- marginaleffects::get_vcov(model, vcov = ~patient_id)
          
          estimates_scaffold %>%
            mutate(
              # this uses the ipw.model to get the estimated incidence at each time point for each treatment, assuming the entire population received treatment A
              # it works correctly for the ATE because of the weights (ie as if setting treatment=1 or treatment=0 for entire population)
              inc = predict(model, newdata = estimates_scaffold, type = "response"),
              inc.se = predict(model, newdata = estimates_scaffold, type = "response", se.fit = TRUE)$se.fit, # this does not use vcov from vcovCL, so not cluster-robust
              inc.logit.se = predict.glm.custom.vcov(model, vcov = vcov, newdata = estimates_scaffold)$se.fit, # cluster robust, but on the linear scale, not response scale
              inc.low = plogis(qlogis(inc) + (qnorm(0.025) * inc.logit.se)),
              inc.high = plogis(qlogis(inc) + (qnorm(0.975) * inc.logit.se)),
            ) |>
            group_by(treatment) %>%
            mutate(
              dummy_id_weight = 1L,
              surv = cumprod(1 - inc),
              surv.se = sqrt(cmlinc_variance(model = model, vcov = vcov, newdata = tibble(treatment = treatment, time = time), id = dummy_id_weight, time = time, weights = dummy_id_weight)),
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

            ) |> select(-dummy_id_weight)

        } else {
          estimates_scaffold
        }
      }
    )
  ) |>
  select(-data, -model) |>
  arrange(!!subgroup_sym)

data_estimates <-
  subgroup_models |>
  unnest(cols = c(estimates, model_status)) |>
  # add NA columns if they don't already exist
  bind_rows(tibble(cmlinc=numeric(), cmlinc.se=numeric(), cmlinc.low=numeric(), cmlinc.high=numeric()))


## output estimates to disk ----
arrow::write_feather(data_estimates, fs::path(output_dir, glue("estimates.arrow")))
write_csv(data_estimates, fs::path(output_dir, glue("estimates.csv")))

# calculate contrasts between treatment groups ----

# could do this within `subgroup_models` step above but bring out to avoid overloading high memory stuff
data_contrasts <-
  
  data_estimates |>
  group_by(!!subgroup_sym) |>
  pivot_wider(
    id_cols = all_of(c(subgroup, "time")),
    names_from = treatment,
    values_from = c(
      cmlinc, cmlinc.se, cmlinc.low, cmlinc.high,
    )
  ) |>
  mutate(
    # see following link for canonical-ish place where these are defined:
    # https://github.com/opensafely-actions/kaplan-meier-function/blob/main/analysis/km.R#L540
    
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


## output to disk ----
arrow::write_feather(data_contrasts, fs::path(output_dir, glue("contrasts.arrow")))
write_csv(data_contrasts, fs::path(output_dir, glue("contrasts.csv")))

data_estimates_time0 <-
  data_estimates |>
  mutate(
    treatment = as.factor(treatment),
    lagtime = lag(time, 1, 0), # assumes the time-origin is zero
  ) |>
  group_by(treatment, !!subgroup_sym) |>
  group_modify(
    ~ add_row(
      .x,
      time = 0, # assumes time origin is zero
      lagtime = 0,
      cmlinc = 0,
      cmlinc.low = 0,
      cmlinc.high = 0,
      .before = 0
    )
  )

plr_plot <-
  data_estimates_time0 |>
  ggplot(aes(group = treatment, colour = treatment, fill = treatment)) +
  geom_step(aes(x = time, y = cmlinc), direction = "vh") +
  geom_step(aes(x = time, y = cmlinc), direction = "vh", linetype = "dashed", alpha = 0.5) +
  geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high), alpha = 0.1, colour = "transparent") +
  facet_grid(rows = vars(!!subgroup_sym)) +
  scale_color_brewer(type = "qual", palette = "Set1", na.value = "grey") +
  scale_fill_brewer(type = "qual", palette = "Set1", guide = "none", na.value = "grey") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(xlim = c(0, NA)) +
  labs(
    x = "Time",
    y = "Cumulative Incidence",
    colour = NULL,
    title = NULL
  ) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(colour = "black"),
    panel.grid.minor.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(.05, .95),
    legend.justification = c(0, 1),
  )

ggsave(filename = fs::path(output_dir, glue("plot.png")), plr_plot, width = 20, height = 20, units = "cm")
