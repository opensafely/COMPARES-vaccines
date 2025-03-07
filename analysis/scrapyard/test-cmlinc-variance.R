### testing plr delta variances

# # # # # # # # # # # # # # # # # # # # #
# Purpose: estimate treatment effect using pooled logistic regression with IPW
# can be used for either IPW or matching adjustment, since matching produces 0/1 weights
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library("optparse")
library('tidyverse')
library('here')
library('glue')
library("arrow")
library("WeightIt")
library("survival")
library('splines')
#library("marginaleffects")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "age75plus"
  method <- "match"
  spec <- "A"
  subgroup <- "all"
  outcome <- "covid_admitted"
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
# data_cohort <- read_feather(here("output", "2-prepare", cohort, "data_cohort.arrow"))

## import weights from matching or weighting method ----
data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "combine", "data_weights.arrow"))
#data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "{method}-{spec}", "data_adjusted.arrow"))

data_all <- 
  data_weights |>
  select(
    patient_id, 
    vax_product, 
    treatment,
    vax_date,
    all_of(weighting_variables[[spec]]),
    all_of(subgroup),
    all_of(paste0(c(outcome, "death", "dereg"), "_date")),
    weight = glue("wt_{cohort}_{method}_{spec}")
  ) |> 
  mutate(
    
    treatment_date = vax_date-1L, # -1 because we assume vax occurs at the start of the day, and so outcomes occurring on the same day as treatment are assumed "1 day" long
    event_date = as.Date(.data[[glue("{outcome}_date")]]),
    
    # person-time is up to and including censor date
    censor_date = pmin(
      dereg_date,
      death_date,
      study_dates$followupend_date,
      treatment_date + maxfup,
      na.rm=TRUE
    ),
    
    noncompetingcensor_date = pmin(
      dereg_date,
      study_dates$followupend_date,
      treatment_date + maxfup,
      na.rm=TRUE
    ),
    
    event_time = tte(treatment_date, event_date, censor_date, na.censor=FALSE),
    event_indicator = censor_indicator(event_date, censor_date),
    
    # possible competing events
    death_time = tte(treatment_date, death_date, noncompetingcensor_date, na.censor=FALSE),
    censor_time = tte(treatment_date, censor_date, censor_date, na.censor=FALSE),
  )

stopifnot("censoring dates must be non-missing" = all(!is.na(data_all$censor_date)))

stopifnot("origin dates must be non-missing" = all(!is.na(data_all$treatment_date)))

times_count <- table(cut(data_all$event_time, c(-Inf, 0, 1, Inf), right=FALSE, labels= c("<0", "0", ">0")), useNA="ifany")
if(!identical(as.integer(times_count), c(0L, 0L, nrow(data_all)))) {
  print(times_count)
  stop("all event times must be strictly positive")
}

formula_time_treatment <- event_indicator ~ treatment + ns(time, 4) + treatment:ns(time, 4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# standard errors calculated using pooled logistic regression,
# but this does not account for uncertainty in the weights
# TODO: incorporate competing risks
# TODO: incorporate standard errors
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data_persontime <- 
  data_all  |>
  uncount(event_time, .remove=FALSE) %>%
  mutate(
    time = sequence(rle(patient_id)$lengths), # equivalent-ish to group_by(patient_id) |> mutate(row_number())
    event_indicator = (time==event_time) & (event_indicator==TRUE)
  )


plr_model <- 
  glm(
    formula_time_treatment, 
    family = binomial(), 
    weight = weight,
    data = data_persontime,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )

library('sandwich')

plr_vcov <- vcovCL(plr_model, cluster = data_persontime$patient_id, type = "HC0")

cmlinc_variance <- function(model, vcov, newdata, id, time, weights){
  
  # calculate variance of adjusted survival / adjusted cumulative incidence
  # needs model object, cluster-vcov, input data, model weights, and indices for patient and time
  
  if(missing(newdata)){ newdata <- model$model }
  tt <- terms(model) # this helpfully grabs the correct spline basis from the model, rather than recalculating based on `newdata`
  Terms <- delete.response(tt)
  m.mat <- model.matrix(Terms, data=newdata)
  m.coef <- model$coef
  
  N <- nrow(m.mat)
  K <- length(m.coef)
  
  # log-odds, nu_t, at time t
  nu <- m.coef %*% t(m.mat) # t_i x 1
  # part of partial derivative
  pdc <- (exp(nu)/((1+exp(nu))^2)) # t_i x 1
  # summand for partial derivative of P_t(theta_t | X_t), for each time t and term k
  
  #summand <- crossprod(diag(as.vector(pdc)), m.mat)    # t_i  x k
  summand <- matrix(0, nrow=N, ncol=K)
  for (k in seq_len(K)){
    summand[,k] <- m.mat[,k] * as.vector(pdc)
  }
  
  # cumulative sum of summand, by patient_id  # t_i x k
  cmlsum <- matrix(0, nrow=N, ncol=K)
  for (k in seq_len(K)){
    cmlsum[,k] <- ave(summand[,k], id, FUN=cumsum)
  }
  
  ## multiply by model weights (weights are normalised here so we can use `sum` later, not `weighted.mean`)
  normweights <- weights / ave(weights, time, FUN=sum) # t_i x 1
  
  #wgtcmlsum <- crossprod(diag(normweights), cmlsum ) # t_i x k
  wgtcmlsum <- matrix(0, nrow=N, ncol=K)
  for (k in seq_len(K)){
    wgtcmlsum[,k] <- cmlsum[,k] * normweights
  }
  
  # partial derivative of cumulative incidence at t
  partial_derivative <- rowsum(wgtcmlsum, time)
  
  variance <- rowSums(crossprod(t(partial_derivative), vcov) * partial_derivative) # t x 1
  
  variance
}


plr_cmlinc <- cmlinc_variance(
  model = plr_model, 
  vcov = plr_vcov, 
  newdata = mutate(data_persontime),
  id = data_persontime$patient_id,
  time = data_persontime$time, 
  weights = data_persontime$weight
)

plr_cmlinc_full0 <- cmlinc_variance(
  model = plr_model, 
  vcov = plr_vcov, 
  newdata = mutate(data_persontime, treatment=0),
  id = data_persontime$patient_id,
  time = data_persontime$time, 
  weights = data_persontime$weight
)

plr_cmlinc_full1 <- cmlinc_variance(
  model = plr_model, 
  vcov = plr_vcov, 
  newdata = mutate(data_persontime, treatment=1),
  id = data_persontime$patient_id,
  time = data_persontime$time, 
  weights = data_persontime$weight
)


newdata0 <- 
  expand_grid(
    treatment =  c(0L),
    time = c(0, seq_len(153)),
    id = 1L
  ) 

newdata1 <- 
  expand_grid(
    treatment =  c(1L),
    time = c(0, seq_len(153)),
    id = 1L
  ) 


plr_cmlinc0 <- cmlinc_variance(
  model = plr_model, 
  vcov = plr_vcov, 
  newdata = newdata0,
  id = newdata0$id,
  time = newdata0$time, 
  weights = newdata0$id
)

plr_cmlinc1 <- cmlinc_variance(
  model = plr_model, 
  vcov = plr_vcov, 
  newdata = newdata1,
  id = newdata1$id,
  time = newdata1$time, 
  weights = newdata1$id
)

data_plot <- 
  bind_rows(newdata0, newdata1) %>%
  mutate(
    inc = predict(plr_model, newdata = ., type="response", se.fit=TRUE)$fit, 
    inc.se = predict(plr_model, newdata = ., type="response", se.fit=TRUE)$se.fit,
  ) %>%
  group_by(treatment) %>%
  mutate(
    lagtime= lag(time,1,0),
    dummy_id_weight = 1L,
    surv = cumprod(1-inc),
    ## probably need Fizz's magic delta method for this in the context of PLR!!
    surv.seA = sqrt(cmlinc_variance(model = plr_model, vcov = plr_vcov, newdata = tibble(treatment=treatment, time=time), id = dummy_id_weight, time = time, weights = dummy_id_weight)),
    surv.lowA = surv + (qnorm(0.025) * surv.seA),
    surv.highA = surv + (qnorm(0.975) * surv.seA),
    treatmentdummy = treatment,
    surv.seB = sqrt(cmlinc_variance(model = plr_model, vcov = plr_vcov, newdata = data_persontime |> mutate(treatment=first(treatmentdummy)), id = data_persontime$patient_id, time = data_persontime$time, weights = data_persontime$weight)),
    surv.lowB = surv + (qnorm(0.025) * surv.seB),
    surv.highB = surv + (qnorm(0.975) * surv.seB),
    
    #rmst = cumsum(surv),
    #rmst.se = sqrt(((2* cumsum(time*surv)) - (rmst^2))/n.risk), # this only works if one row per day using fill_times! otherwise need sqrt(((2* cumsum(time*interval*surv)) - (rmst^2))/n.risk)
    #rmst.low = rmst + (qnorm(0.025) * rmst.se),
    #rmst.high = rmst + (qnorm(0.975) * rmst.se),
  ) |>
  mutate(
    treatment = as.factor(treatment)
  )


plr_plot <-
  data_plot |>
  ggplot(aes(group = treatment, colour = treatment, fill = treatment)) +
  geom_line(aes(x = time, y = surv)) +
  geom_line(aes(x = time, y = surv.lowA), linetype="dashed") +
  geom_line(aes(x = time, y = surv.highA), linetype="dashed") +
  geom_line(aes(x = time, y = surv.lowB), linetype="dotted") +
  geom_line(aes(x = time, y = surv.highB), linetype="dotted") +
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
