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

if(subgroup=="all") data_weights$all <- 1L

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

formula_timesincevax_ns <- event_indicator ~ treatment + ns(timesincevax, 4) + treatment:ns(timesincevax, 4)



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
    timesincevax = sequence(rle(patient_id)$lengths), # equivalent to group_by(patient_id) |> mutate(row_number())
    event_indicator = (timesincevax==event_time) & (event_indicator==TRUE)
  )

ipw.model <- glm(
  formula_timesincevax_ns, 
  family = binomial(), 
  weight = weight,
  data = data_persontime
)

data_surv  <- 
  expand_grid(
    treatment =  c(0L, 1L),
    timesincevax = c(0,seq_len(maxfup))
  ) %>%
  mutate(
    # this using the ipw.model to get the estimated incidence at each time point for each treatment, assuming the entire population received treatment A
    # it works correctly for the ATE because of the weights (ie as if setting treatment=1 or treatment=0 for entire population)
    inc = predict(ipw.model, ., type="response", se.fit=TRUE)$fit, 
    inc.se = predict(ipw.model, ., type="response", se.fit=TRUE)$se.fit
  ) |>
  group_by(treatment) |>
  mutate(
    surv = cumprod(1-inc),
    # TODO: get standard errors for survival
    ## probably need Fizz's magic delta method for this in the context of PLR!!
    surv.se = inc.se, # TODO: this is completely wrong, need to fix, but placeholder so there are _some_ confidence limits for now to check workflow
    surv.low = surv + (qnorm(0.025) * inc.se),
    surv.high = surv + (qnorm(0.975) * inc.se),
    
    cmlinc = 1 - surv,
    cmlinc.se = surv.se,
    cmlinc.low = 1 - surv.high,
    cmlinc.high = 1 - surv.low,
    #rmst = cumsum(surv),
    #rmst.se = sqrt(((2* cumsum(time*surv)) - (rmst^2))/n.risk), # this only works if one row per day using fill_times! otherwise need sqrt(((2* cumsum(time*interval*surv)) - (rmst^2))/n.risk)
    #rmst.low = rmst + (qnorm(0.025) * rmst.se),
    #rmst.high = rmst + (qnorm(0.975) * rmst.se),
  )

data_contrasts <-
  # see following link for canonical-ish place where these are defined https://github.com/opensafely-actions/kaplan-meier-function/blob/main/analysis/km.R#L540
  data_surv |>
  pivot_wider(
    id_cols = c(timesincevax),
    names_from = treatment,
    values_from = c(
      cmlinc, cmlinc.se, cmlinc.low,cmlinc.high,
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

ggplot(data_surv)+ 
  geom_line(aes(x = timesincevax, y = surv, group=treatment, colour = treatment)) + 
  xlab("days") + 
  ylab("Survival") + 
  ggtitle("Survival from IP weighted hazards model") + 
  labs(colour="A:") +
  theme_bw() + 
  theme(legend.position="bottom")

