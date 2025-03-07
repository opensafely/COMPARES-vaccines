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
library('splines')
library("splitstackshape")
library("marginaleffects")

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

output_dir <- here_glue("output", "3-adjust", cohort, "{method}-{spec}", "report")
fs::dir_create(output_dir)

## import unadjusted cohort data ----
data_cohort <- read_feather(here("output", "2-prepare", cohort, "data_cohort.arrow"))

if(subgroup=="all") data_cohort$all <- 1L

## import weights from matching or weighting method ----
data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "{method}-{spec}", "data_adjusted.arrow"))

data_all <- 
  data_cohort |>
  select(
    patient_id, 
    treatment, 
    vax_date,
    all_of(weighting_variables[[spec]]),
    all_of(subgroup),
    all_of(paste0(c(outcome, "death", "dereg"), "_date")),
    # outcomes
    # subgroups
    
  ) |>
  left_join(
    data_weights |> select(patient_id, weight),
    by = "patient_id"
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


##############################################################
# from What If
# no standard errors calculated
##############################################################
data_ipw <- expandRows(data_all, "event_time", drop=F)

data_ipw$timesincevax <- sequence(rle(data_ipw$patient_id)$lengths)
data_ipw$event_indicator <- (data_ipw$timesincevax==data_ipw$event_time) & (data_ipw$event_indicator==TRUE)

ipw.model <- glm(
  formula_timesincevax_ns, 
  family=binomial(), 
  weight=weight,
  data=data_ipw
)

summary(ipw.model)

# creation of survival curves
ipw.treat0 <- data.frame(cbind(seq(0, maxfup),0))
ipw.treat1 <- data.frame(cbind(seq(0, maxfup),1))

colnames(ipw.treat0) <- c("timesincevax", "treatment")
colnames(ipw.treat1) <- c("timesincevax", "treatment")

# assignment of estimated (1-hazard) to each person-month */
ipw.treat0$p.noevent0 <- 1 - predict(ipw.model, ipw.treat0, type="response")
ipw.treat1$p.noevent1 <- 1 - predict(ipw.model, ipw.treat1, type="response")

# computation of survival for each person-month
ipw.treat0$surv0 <- cumprod(ipw.treat0$p.noevent0)
ipw.treat1$surv1 <- cumprod(ipw.treat1$p.noevent1)

# some data management to plot estimated survival curves
ipw.graph <- merge(ipw.treat0, ipw.treat1, by=c("timesincevax"))
ipw.graph$survdiff <- ipw.graph$surv1-ipw.graph$surv0


# plot
ggplot(ipw.graph, aes(x=timesincevax)) + 
  geom_line(aes(y = surv0, colour = "0")) + 
  geom_line(aes(y = surv1, colour = "1")) + 
  xlab("days") + 
  ylab("Survival") + 
  ggtitle("Survival from IP weighted hazards model") + 
  labs(colour="A:") +
  theme_bw() + 
  theme(legend.position="bottom")


##############################################################
# from prior sequential trial work in opensafely
# standard errors calculated but do not account for uncerainty in the weights
##############################################################

# TODO: incorporate competing risks
# TODO: incorporate standard errors

data_persontime <- 
  data_all  |>
  expandRows("event_time", drop=F) %>%
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

data_curves  <- 
  expand_grid(
    treatment =  c(0L, 1L),
    timesincevax = c(0,seq_len(maxfup))
  ) %>%
  mutate(
    probevent = predict(ipw.model, ., type="response"),
  ) |>
  group_by(treatment) |>
  mutate(
    surv = cumprod(1-probevent),
  )

data_curves_constrast <-
  data_curves |>
  pivot_wider(
    id_cols = c(timesincevax),
    names_from = treatment,
    values_from = c(probevent, surv)
  ) |>
  mutate(
    diffprob = probevent_1 - probevent_0,
    diffsurv = surv_1 - surv_0
  )


ggplot(data_curves)+ 
  geom_line(aes(x = timesincevax, y = surv, group=treatment, colour = treatment)) + 
  xlab("days") + 
  ylab("Survival") + 
  ggtitle("Survival from IP weighted hazards model") + 
  labs(colour="A:") +
  theme_bw() + 
  theme(legend.position="bottom")



##############################################################
# calculate balancing weights using the coxph_weightit function
# with period-specific hazards - one period per day where there is an event
##############################################################


data_persontime <- 
  data_all  |>
  expandRows("event_time", drop=F) %>%
  mutate(
    timesincevax = sequence(rle(patient_id)$lengths), # equivalent to group_by(patient_id) |> mutate(row_number())
    event_indicator = (timesincevax==event_time) & (event_indicator==TRUE)
  ) |> 
  group_by(timesincevax) |> 
  mutate(count = sum(event_indicator)) |> 
  filter(count>0) |> 
  ungroup()

obj_weightitMSM <-
  weightitMSM(
    formula = list(
      formula(paste0("treatment ~ ", weighting_formulae[[spec]]))
    ),
    data = data_persontime,
    method = "glm",
    estimand = "ATE",
    stabilize = TRUE
  )

coxph <- 
  coxph_weightit(
    survival::Surv(event_time, event_indicator) ~ treatment + treatment:strata(timesincevax),
    data = data_persontime,
    weightit = obj_weightitMSM,
    vcov = "HC0",
    cluster = ~patient_id
  )




##############################################################
# calculate usinf adjustedsurv function from adjustedCurves package
# does not account for uncertainty in the weights
# does not permit clog-log confidence intervals
# useful to check calibration of parametric approaches!
##############################################################


library("adjustedCurves")
adjsurv <- 
  adjustedsurv(
    data = data_all |> mutate(treatment_fct = as.factor(treatment)),
    variable = "treatment_fct",
    ev_time = "event_time",
    event = "event_indicator",
    conf_int = TRUE,
    method = "iptw_km",
    treatment_model = data_all$weight,
  )
plot(adjsurv)

ggplot(adjsurv$adj)+ 
  geom_line(aes(x = time, y = surv, group=group, colour = group)) + 
  geom_line(aes(x = time, y = ci_lower, group=group, colour = group))+
  geom_line(aes(x = time, y = ci_upper, group=group, colour = group))+
  xlab("days") + 
  ylab("Survival") + 
  ggtitle("Survival from IP weighted hazards model") + 
  labs(colour="A:") +
  theme_bw() + 
  theme(legend.position="bottom")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# calculate balancing weights using the weightitMSM function
# if done properly this should account for uncertainty in the weights using M-estimation
# but it spits out estimates for discrete-time hazards, not cumulative incidences
# cannot easily derive correct CIs for the cumulative incidences from this output
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Park for now!


data_persontime <-
  data_all  |>
  expandRows("event_time", drop=F) %>%
  mutate(
    timesincevax = sequence(rle(patient_id)$lengths), # equivalent to group_by(patient_id) |> mutate(row_number())
    event_indicator = (timesincevax==event_time) & (event_indicator==TRUE)
  )

obj_weightitMSM <-
  weightitMSM(
    formula = list(
      formula(paste0("treatment ~ ", weighting_formulae[[spec]]))
    ),
    data = data_persontime,
    method = "glm",
    estimand = "ATE",
    stabilize = TRUE
  )

ipw.model <- 
  glm_weightit(
    formula_timesincevax_ns,
    family = binomial(),
    weightit = obj_weightitMSM,
    data = data_persontime,
    vcov = "asympt"# "asympt" = m-estimation
  )


## this is (in theory, not quite numerically) the same as using avg_predictions independently in population-wide datasets filtered by timesincevax
## (though this may not work in general if there are patient-specific confounders)
## this is _not_ the same as using avg_predictions across the whole dataset (ie, newdata = original_dataset)
## because using across the whole dataset would replace _all_ timesincevax values with the current iteration of timescince vax (making an enoromous patient_ID * timesincevax^2 dataset)

tictoc::tic()
avg_curves1 <-
  avg_predictions(
    ipw.model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      timesincevax = c(0,seq_len(maxfup))
    ),
    variables = list(
      treatment = c(0L, 1L),
      timesincevax = c(0,seq_len(maxfup))
    )
  )
tictoc::toc()


tictoc::tic()
avg_curves2 <-
  avg_predictions(
    ipw.model,
    variables = list(
      treatment = c(0L, 1L),
      timesincevax = c(seq_len(15))
    )
  )
tictoc::toc()


tictoc::tic()
# this is exactly equivalent to running avg_predictions independently 
# for the original dataset filtered on each timesincevax value
# but for some reason is a bit slower
avg_curves3 <-
  avg_predictions(
    ipw.model,
    variables = list(
      treatment = c(0L, 1L)
    ),
    by = c("treatment", "timesincevax")
  )
tictoc::toc()

tictoc::tic()
avg_curves4 <-
  avg_predictions(
    ipw.model,
    variables = list(
      treatment = c(0L, 1L)
    ),
    by = c("treatment", "timesincevax"),
    vcov = "HC0"
  )
tictoc::toc()

tictoc::tic()
avg_curves5 <-
  avg_predictions(
    ipw.model,
    variables = list(
      treatment = c(0L, 1L)
    ),
    by = c("treatment", "timesincevax"),
    vcov = ~patient_id
  )
tictoc::toc()

test12345 <- 
  avg_curves1 |> 
  filter(timesincevax %in% 1:15) |>
  transmute(treatment, timesincevax, estimate1=estimate, std.error1=std.error, conf.low1=conf.low, conf.high1 = conf.high) |>
  left_join(
    avg_curves2 |> transmute(treatment, timesincevax, estimate2=estimate, std.error2=std.error, conf.low2=conf.low, conf.high2 = conf.high),
    by= c("treatment", "timesincevax")
  ) |>
  left_join(
    avg_curves3 |> transmute(treatment, timesincevax, estimate3=estimate, std.error3=std.error, conf.low3=conf.low, conf.high3 = conf.high),
    by= c("treatment", "timesincevax")
  ) |>
  left_join(
    avg_curves4 |> transmute(treatment, timesincevax, estimate4=estimate, std.error4=std.error, conf.low4=conf.low, conf.high4 = conf.high),
    by= c("treatment", "timesincevax")
  ) |>
  left_join(
    avg_curves5 |> transmute(treatment, timesincevax, estimate5=estimate, std.error5=std.error, conf.low5=conf.low, conf.high5 = conf.high),
    by= c("treatment", "timesincevax")
  ) 


ggplot(avg_curves1 |> filter(timesincevax %in% 1:15))+
  geom_line(aes(x=timesincevax, y=estimate, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate-2*std.error, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate+2*std.error, colour=treatment, group=treatment))+
  scale_y_continuous(limits = c(0.001, 0.004))+
  theme_bw()

ggplot(avg_curves2)+
  geom_line(aes(x=timesincevax, y=estimate, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate-2*std.error, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate+2*std.error, colour=treatment, group=treatment))+
  scale_y_continuous(limits = c(0.001, 0.004))+
  theme_bw()

ggplot(avg_curves3)+
  geom_line(aes(x=timesincevax, y=estimate, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate-2*std.error, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate+2*std.error, colour=treatment, group=treatment))+
  scale_y_continuous(limits = c(0.001, 0.004))+
  theme_bw() 

ggplot(avg_curves4|> filter(timesincevax %in% 1:15))+
  geom_line(aes(x=timesincevax, y=estimate, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate-2*std.error, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate+2*std.error, colour=treatment, group=treatment))+
  scale_y_continuous(limits = c(0.001, 0.004))+
  theme_bw() 

ggplot(avg_curves5|> filter(timesincevax %in% 1:15))+
  geom_line(aes(x=timesincevax, y=estimate, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate-2*std.error, colour=treatment, group=treatment))+
  geom_line(aes(x=timesincevax, y=estimate+2*std.error, colour=treatment, group=treatment))+
  scale_y_continuous(limits = c(0.001, 0.004))+
  theme_bw() 


# This provides contrasts for the discrete-time hazards, not for the cumulative survival / incidence!
# TODO: can we get this with marginaleffects with appropriate standard errors? delta method?
model_comparisons <-
  comparisons(
    ipw.model,
    newdata = expand_grid(
      treatment =  c(0L, 1L),
      timesincevax = c(0,seq_len(maxfup))
    ),
    variables = list(
      treatment =  c(0L, 1L)
    )
  )

model_comparisons2 <-
  avg_comparisons(
    ipw.model,
    variables = list(
      treatment =  c(0L, 1L)
    )
  )

model_comparisons3 <-
  comparisons(
    ipw.model,
    variables = list(
      treatment =  c(0L, 1L)
    )
  )


