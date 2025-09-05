# # # # # # # # # # # # # # # # # # # # #
# Purpose: Get disclosure-safe Aalen-Johansen estimates.
# This is very heavily based on the Kaplan-Meier reusable action: https://github.com/opensafely-actions/kaplan-meier-function
# The function requires an origin date, an event date, a censoring date, and optionally a competing event date, which are converted into a (time , status) pair that is passed to `survival::Surv`
# If competing_date is not supplied, then this _should_ reduce to the Kaplan-Meier function (so some details about confidence limits might differ slightly)
# Estimates are stratified by the `exposure` variable, and additionally by any `subgroups`
# Counts are rounded to midpoint values defined by `count_min`.
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----

library('here')
library('glue')
library('tidyverse')
library('survival')


## import local functions ----

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

## Import design elements
source(here("analysis", "0-lib", "time-rounding.R"))


## parse command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  df_input <- "output/3-adjust/age75plus/combine/data_weights.arrow"
  dir_output <- "output/4-contrast/test/"
  exposure <- c("treatment")
  subgroups <- c("all")
  origin_date <- "vax_date"
  event_date <- "covid_admitted_date"
  censoring_date <- "dereg_date" 
  competing_date <- "death_date"
  weight <- "wt_age75plus_match_A"
  min_count <- as.integer("6")
  method <- "constant"
  max_fup <- as.numeric("168")
  concise <- as.logical("TRUE")
  plot <- as.logical("TRUE")
  contrast <- as.logical("TRUE")
  filename_suffix <- as.character("")
} else {
  
  library("optparse")
  
  option_list <- list(
    make_option("--df_input", type = "character",
                help = "[default: Must be specified] character. The input dataset .arrow filename. feather/arrow format is enforced to ensure date types are preserved.",
                metavar = "df_input.arrow"),
    make_option("--dir_output", type = "character",
                help = "[default: must be specified] character. The output directory. All requested output files (eg 'estimates.arrow', 'contrasts.arrow') will be placed in this directory. See also: 'filename_suffix' argument.",
                metavar = "/output/"),
    make_option("--exposure", type = "character", default = character(),
                help = "[default: NULL] character. The name of an exposure variable in the input dataset. Must be binary or not given. All outputs will be stratified by this variable. This could be an exposure in the usual sense, or it could (mis)used to show different types of events (as long as the censoring structure is the same). If not specified, no stratification will occur.",
                metavar = "exposure_varname"),
    make_option("--subgroups", type = "character", default = character(),
                help = "[default: NULL] The name of a subgroup variable or list of variable names. If a subgroup variable is used, analyses will be stratified as exposure * ( subgroup1, subgroup2, ...). If not specified, no stratification will occur.",
                metavar = "subgroup_varname"),
    make_option("--origin_date", type = "character",
                help = "[default: must be specified] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') in the input dataset that represents the start of follow-up.",
                metavar = "origin_varname"),
    make_option("--event_date", type = "character",
                help = "[default: must be specified] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') in the input dataset that represents the event date.",
                metavar = "event_varname"),
    make_option("--censoring_date", type = "character", default = character(),
                help = "[default: NULL] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') that represents the censoring event date. If not specified, then no censoring occurs except at `max_fup` time.",
                metavar = "censor_varname"),
    make_option("--competing_date", type = "character", default = character(),
                help = "[default: NULL] The name of a date variable (or name of a variable that is coercable to a date eg 'YYYY-MM-DD') that represents the competing event date. If not specified, then no competing event is used and the estimates reduce to Kaplan-Meier estimates.",
                metavar = "censor_varname"),
    make_option("--weight", type = "character", default = character(),
                help = "[default: NULL] The name of a numeric variable that represents balancing / sampling weights. If not specified, then no weighting occurs.",
                metavar = "censor_varname"),
    make_option("--min_count", type = "integer", default = 6,
                help = "[default: %default] integer. The minimum permissable event and censor counts for each 'step' in the AJ curve. This ensures that at least `min_count` events occur at each event time.",
                metavar = "min_count"),
    make_option("--method", type = "character", default = "constant",
                help = "[default: %default] character. The interpolation method after rounding. The 'constant' method leaves the event times unchanged after rounding, making the AJ curve have bigger, fewer steps. The 'linear' method linearly interpolates between rounded events times (then rounds to the nearest day), so that the steps appear more natural.",
                metavar = "method"),
    make_option("--max_fup", type = "numeric", default = Inf,
                help = "[default: %default] numeric. The maximum follow-up time after the origin date. If event variables are dates, then this will be days.",
                metavar = "max_fup"),
    make_option("--concise", type = "logical", default = TRUE,
                help = "[default: %default] logical. Should the outputted table only report core variables (defined here as exposure, subgroups, time, number at risk, cumulative number of events, cumulative incidence, and confidence limits) (TRUE) or should it report everything (FALSE)?",
                metavar = "TRUE/FALSE"),
    make_option("--plot", type = "logical", default = FALSE,
                help = "[default: %default] logical. Should Kaplan-Meier plots be created in the output directory? These are fairly basic plots for sense-checking purposes.",
                metavar = "TRUE/FALSE"),
    make_option("--contrast", type = "logical", default = TRUE,
                help = "[default: %default] logical. Should Kaplan-Meier curves for a binary exposure be compared to estimate risk difference, risk ratio, and survival ratio? Ignored if exposure is not supplied.",
                metavar = "TRUE/FALSE"),
    make_option("--filename_suffix", type = "character", default = "",
                help = "[default: %default] character. This will be appended to the end of all outputted files. This is useful if you want to re-run this action across different arguments, but put outputs from all actions in the same directory.",
                metavar = "TRUE/FALSE")
  )
  
  opt_parser <- OptionParser(
    usage = "kaplan-meier-function:[version] [options]",
    option_list = option_list
  )
  
  opt <- parse_args(opt_parser)
  list2env(opt, .GlobalEnv)
}

# Use quasiquotation for passing exposure and subgroup stratification variables
# around the place
# use `syms()` instead of `sym()` even though it's possible to pull
# only one exposure or subgroup variable is from the args (without hacking!)
# this ensures that if `exposure` or `subgroups` is not used,
# the quasiquotation still works inside ggplot, transmute, etc

exposure_syms <- syms(exposure)
subgroup_syms <- syms(subgroups)

# Create output directory ----

dir_output <- here::here(dir_output)
fs::dir_create(dir_output)


# Import and process person-level data  ----

## Import ----
data_patients <-
  arrow::read_feather(here::here(df_input))

## Derive variables ----

if(length(censoring_date)==0) {
  # if censor date is not specified, then create a censoring_date variable in the dataset, taking value `as.Date(Inf)`
  data_patients$censoring_date <- as.Date(Inf)
  censoring_date <- "censoring_date"
}

if(length(competing_date)==0) {
  # if competing date is not specified, then create a competing_date variable in the dataset, taking value `as.Date(Inf)`
  data_patients$competing_date <- as.Date(Inf)
  competing_date <- "competing_date"
}

is.weighted <- length(weight)>0
if(!is.weighted) {
  # if weight is not specified, then create a weight variable in the dataset, taking value `1L`
  data_patients$.weight <- 1L
  weight <- ".weight"
}

data_tte <-
  data_patients |>
  transmute(
    patient_id,
    !!!exposure_syms,
    !!!subgroup_syms,
    .weight = .data[[weight]],
    event_date = as.Date(.data[[event_date]]),
    origin_date = as.Date(.data[[origin_date]]),
    censoring_date = pmin(
      as.Date(.data[[censoring_date]]),
      origin_date + max_fup,
      na.rm=TRUE
    ),
    competing_date = as.Date(.data[[competing_date]]),
    endpoint_date = pmin(censoring_date, event_date, competing_date, na.rm=TRUE),
    event_time = as.integer(endpoint_date - origin_date),
    
    event_indicator = !((event_date>pmin(censoring_date, competing_date, na.rm=TRUE) | is.na(event_date))),

    event_status = case_when(
      # 0 = censored, 1 = event of interest, 2 = competing event
      # could use which.min(c(event_date, competing_date, censoring_date)) for this but it needs to be with rowwise() tie handling less explicit
      endpoint_date == event_date ~ 1L,
      endpoint_date == competing_date ~ 2L,
      endpoint_date == censoring_date ~ 0L, # put censoring_date last so that if there is a tie, the event wins
      .default = NA_integer_
    ),
    # here censor value is first in the factor to work with Surv as required
    event_status_factor = factor(event_status, levels=c("0", "1", "2"), labels = c("censored", "outcome", "competing")), 
  )

if(max_fup==Inf) max_fup <- max(data_tte$event_time)+1

## tests ----

if(length(exposure)>0){
  stopifnot("exposure variable must be binary or have two levels" = (length(unique(data_patients[[exposure]])) == 2))
}

stopifnot("censoring dates must be non-missing" = all(!is.na(data_tte$censoring_date)))

stopifnot("origin dates must be non-missing" = all(!is.na(data_tte$origin_date)))

times_count <- table(cut(data_tte$event_time, c(-Inf, 0, 1, Inf), right=FALSE, labels= c("<0", "0", ">0")), useNA="ifany")

if(!identical(as.integer(times_count), c(0L, 0L, nrow(data_tte)))) {
  print(times_count)
  stop("all event times must be strictly positive")
}

# Calculate max follow-up time available in the data ----
# and print to log file

if(length(exposure)>0){
  max_time_data <-
    data_tte |>
    group_by(!!!exposure_syms) |>
    summarise(
      max_fup_time = max(event_time),
      max_event_time = max(event_time[event_indicator])
    )
  
  cat("maximum follow-up time in exposure levels [", paste0(max_time_data[[exposure]], collapse=", "), "] is [", paste0(max_time_data$max_fup_time, collapse= ", "), "]", "\n")
  cat("maximum event time in exposure levels [", paste0(max_time_data[[exposure]], collapse=", "), "] is [", paste0(max_time_data$max_event_time, collapse= ", "), "]", "\n")
} else {
  max_time_data <-
    data_tte |>
    summarise(
      max_fup_time = max(event_time),
      max_event_time = max(event_time[event_indicator])
    )
  cat("maximum follow-up time is [", paste0(max_time_data$max_fup_time, collapse= ", "), "]", "\\n")
  cat("maximum event time is [", paste0(max_time_data$max_event_time, collapse= ", "), "]", "\\n")
}


# Calculate AJ estimates ------

## Run `survfit` across each level of exposure and subgroup ----
## do this independently rather than using stratification or covariates
## because it makes variable name handling easier

# see https://github.com/opensafely/comparative-booster/blob/bca54292baa80e967187ca28988d4897ae88aedc/analysis/ci.R#L189
# for prototype



# function to recalculate standard-error of AJ cumulative incidence function based on count-rounded life tables
cif.se <- function(time, ci, n.risk, n.event, kmsurv, kmsummand){
  # from https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Cumulative_Incidence.pdf
  # also here: https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.200900039?saml_referrer
  timeindex <- seq_along(time)
  lapply(timeindex, function(i) {
    cii <- ci[1:i]
    bi <- ((n.risk - n.event)/n.risk)[1:i]
    di <- (n.event/(n.risk^2))[1:i]
    lagkmi <- lag(kmsurv,1,1)[1:i]
    kmsummandi <- kmsummand[1:i]
    vt <-
      sum(((cii[i] - cii)^2) * kmsummandi) +
      sum((lagkmi^2) * bi * di) +
      -2* sum((cii[i] - cii) * lagkmi * di)
    
    sqrt(vt)
  })
}

# for each exposure level and subgroup level, pass data through `survival::survfit` to get survival table
data_surv <-
  data_tte |>
  dplyr::group_by(!!!subgroup_syms, !!!exposure_syms) |>
  tidyr::nest() |>
  dplyr::mutate(
    surv_obj_tidy_0 = purrr::map(data, ~ {
      survival::survfit(
        survival::Surv(event_time, event_status_factor) ~ 1,
        data = .x,
        conf.type="log",
        weight = .weight
      ) |>
      broom::tidy() |>
      tidyr::complete(
        time = seq_len(max_fup), # fill in 1 row for each day of follow up
        fill = list(n.event = 0L, n.censor = 0L) # fill in zero events on those days
      ) |>
      tidyr::fill(n.risk, .direction = c("up"))
    }),
    
    surv_obj_tidy = purrr::map(data, ~ {
      surv_obj_tidy_long <- 
        survival::survfit(
          survival::Surv(event_time, event_status_factor) ~ 1,
          data = .x,
          conf.type="log-log",
          weight = .weight
        ) |>
        broom::tidy() 
      
      surv_obj_tidy <-
        bind_cols(
          surv_obj_tidy_long %>% filter(state=="(s0)") %>% select(time, n.risk, n.censor),
          surv_obj_tidy_long %>% filter(state=="outcome") %>% select(n.event, estimate, std.error, conf.low, conf.high),
          surv_obj_tidy_long %>% filter(state=="competing") %>% select(n.competing = n.event),
        ) |>
        transmute(
          
          time,
          n.risk, n.event, n.censor, n.competing,
          
          n.allevents = n.event + n.competing,
          
          # in case it's worth comparing when treating competing event as a censoring event too
          # kmsummand = (1/(n.risk-n.event)) - (1/n.risk), # = n.event / ((n.risk - n.event) * n.risk) but re-written to prevent integer overflow
          # kmsurv = cumprod(1 - n.event / n.risk),
          # kmsurv.se = kmsurv * sqrt(cumsum(kmsummand)),
          
          risk = estimate,
          risk.se = std.error,
          risk.ll = conf.low,
          risk.ul = conf.high,
          surv = 1 - risk,
          surv.se = std.error,
          surv.ll = 1 - risk.ul,
          surv.ul = 1 - risk.ll,
          
          #allsummand = (1/(n.risk-n.allevents)) - (1/n.risk),
          #allsurv = cumprod(1 - n.allevents / n.risk),
          #risk.se2 = cif.se(time, estimate, n.risk, n.event, allsurv, allsummand) # should be very similar to risk.se
        ) |>
      tidyr::complete(
        time = seq_len(max_fup), # fill in 1 row for each day of follow up
        fill = list(n.event = 0L, n.allevents = 0, n.censor = 0L, n.competing = 0L) # fill in zero events on those days
      ) |>
      tidyr::fill(n.risk, .direction = c("up"))
      
      surv_obj_tidy
    }),
  ) |>
  select(-data, -surv_obj_tidy_0) |>
  tidyr::unnest(surv_obj_tidy)

## Round the count values in the survival data ----
# round event times such that no event time has fewer than `min_count` events
# recalculate AJ estimates based on these rounded event times

round_aj <- function(.data, min_count=0, method="constant") {
  
  # min_count == 0 means no rounding.
  # precision is 0 rather than 1 because if using weighting, then
  # there may be non-integer counts
  if(min_count==0){
    rounded_data <-
      .data |>
      mutate(
        N = max(n.risk, na.rm = TRUE),
        # rounded to `min_count - (min_count/2)`
        cml.event = cumsum(n.event),
        cml.censor = cumsum(n.censor),
        cml.competing = cumsum(n.competing),
        cml.any = cml.event + cml.censor + cml.competing,
      )
  } else {
    rounded_data <-
      .data |>
      mutate(
        # rounded to `min_count - (min_count/2)`
        cml.event = round_cmlcount(cumsum(n.event), time, min_count, method, integer.counts = !is.weighted),
        cml.censor = round_cmlcount(cumsum(n.censor), time, min_count, method, integer.counts = !is.weighted),
        cml.competing = round_cmlcount(cumsum(n.competing), time, min_count, method, integer.counts = !is.weighted),
        cml.any = cml.event + cml.censor + cml.competing,
        # re-derive counts from cumulative data
        n.event = diff(c(0, cml.event)),
        n.censor = diff(c(0, cml.censor)),
        n.competing = diff(c(0, cml.competing)),
        n.risk = roundmid_any(max(n.risk, na.rm = TRUE), min_count) - lag(cml.any, 1, 0)
      )
  }
  
  rounded_data1 <-
    rounded_data |>
    mutate(
      cml.nrisk = cumsum(n.risk),
      cml.rate = cml.event / cml.nrisk,
      inc = n.event / n.risk,
      
      # KM estimate for event of interest, combining censored and competing events as censored
      kmsummand = (1/(n.risk - n.allevents)) - (1/n.risk), # n.allevents / ((n.risk - n.allevents) * n.risk) but re-written to prevent integer overflow
      kmsurv = cumprod(1 - n.allevents / n.risk),
      # summand for the event of interest
      summand = (n.event / n.risk) * lag(kmsurv, 1, 1),
      
      # AJ (Aalen-Johansen) estimated cumulative incidence for the event of interest, accounting for competing and censoring events appropriately
      cmlinc = cumsum(summand),
      
      # two options currently for redoing standard errors and confidence limits:
      # either recalculate based on rounded values using the cif.se function, but this do not match what we get from survfit
      # or use the original standard.error and limits without adjusting for rounding. this is basically fine given non-disclosivity of such statistics
      cmlinc.se = as.numeric(cif.se(time, cmlinc, n.risk, n.event, kmsurv, kmsummand)), 
      cmlinc.ln.se = cmlinc.se / (1-cmlinc), # possibly approximate!
      cmlinc.low = exp(log(cmlinc) + qnorm(0.025)*cmlinc.ln.se),
      cmlinc.high = exp(log(cmlinc) + qnorm(0.975)*cmlinc.ln.se),
      
      surv = 1 - cmlinc, 
      # standard errors on survival scale
      surv.se = cmlinc.se,
      surv.low = 1 - cmlinc.high,
      surv.high = 1 - cmlinc.low,
      ## standard errors on log scale
      surv.ln.se = cmlinc.ln.se,

      ## standard errors on complementary log-log scale
      ## don't yet know how to do this for AJ estimates
      # surv.cll = log(-log(surv)), # this is equivalent to the log cumulative hazard
      # surv.cll.se = if_else(surv==1, 0, sqrt((1 / log(surv)^2) * cumsum(__summand__))), # assume SE is zero until there are events -- makes plotting easier
      # surv.low = exp(-exp(surv.cll + qnorm(0.975) * surv.cll.se)),
      # surv.high = exp(-exp(surv.cll + qnorm(0.025) * surv.cll.se)),

      # restricted mean survival time.
      # https://doi.org/10.1186/1471-2288-13-152
      rmst = cumsum(surv), # this only works if one row per day using fill_times! otherwise need cumsum(surv*int)
      rmst.se = sqrt(((2* cumsum(time*surv)) - (rmst^2))/n.risk), # this only works if one row per day using fill_times! otherwise need sqrt(((2* cumsum(time*interval*surv)) - (rmst^2))/n.risk)
      rmst.low = rmst + (qnorm(0.025) * rmst.se),
      rmst.high = rmst + (qnorm(0.975) * rmst.se),

    ) |>
    # filter(
    #   !(n.event==0 & n.censor==0 & !fill_times) # remove times where there are no events (unless all possible event times are requested with fill_times)
    # ) |>
    select(
      !!!subgroup_syms,
      !!!exposure_syms,
      time,
      cml.nrisk, cml.event, cml.censor,
      n.risk, n.event, n.censor, n.competing,
      inc,
      #surv, surv.se, surv.low, surv.high,
      cml.rate,
      cmlinc, cmlinc.se, cmlinc.low, cmlinc.high, 
      rmst, rmst.se, rmst.low, rmst.high,
    )
  
  rounded_data1
}

data_surv_unrounded <- round_aj(data_surv, 0L)
data_surv_rounded <- round_aj(data_surv, min_count, method=method)

## Select a smaller set of variable ----
## if requested via concise argument

data_surv_rounded_output <-
  if(concise){
    data_surv_rounded |>
      select(
        !!!subgroup_syms,
        !!!exposure_syms,
        time,
        n.risk, n.censor, n.event, n.competing,
        cml.event,
        cmlinc, cmlinc.low, cmlinc.high
      )
  } else {
    data_surv_rounded_output <- data_surv_rounded
  }


## output to disk ----
## include both arrow and csvv formats here - if you don't want one of them,
## don't include it in the `output:` slot in the action

## write arrow to disk
arrow::write_feather(data_surv_rounded_output, fs::path(dir_output, glue("estimates{filename_suffix}.arrow")))
## write csv to disk
write_csv(data_surv_rounded_output, fs::path(dir_output, glue("estimates{filename_suffix}.csv")))


# Plot AJ curves ----
## if requested via `plot` argument

aj_plot <- function(.data) {
  
  data_with_time0 <-
    .data |>
    mutate(
      "{exposure}" := as.factor(!!!exposure_syms),
      lagtime = lag(time, 1, 0), # assumes the time-origin is zero
    ) %>%
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
  ggplot_init <- if(length(exposure)==0L){
    ggplot(data_with_time0)
  } else {
    exposure_sym <- sym(exposure)
    ggplot(data_with_time0, aes(group = !!exposure_sym, colour = !!exposure_sym, fill = !!exposure_sym))
  }
  ggplot_init +
    geom_step(aes(x = time, y = cmlinc), direction = "vh") +
    geom_step(aes(x = time, y = cmlinc), direction = "vh", linetype = "dashed", alpha = 0.5) +
    geom_rect(aes(xmin = lagtime, xmax = time, ymin = cmlinc.low, ymax = cmlinc.high), alpha = 0.1, colour = "transparent") +
    facet_grid(rows = vars(!!!subgroup_syms)) +
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
}

## output to disk ----
if(plot){
  aj_plot <- aj_plot(data_surv_rounded)
  ggsave(filename = fs::path(dir_output, glue("plot{filename_suffix}.png")), aj_plot, width = 20, height = 20, units = "cm")
}


# Calculate contrasts between exposure groups ----
# if requested via 'contrast` argument

contrast_aj <- function(.data) {
  
  .data |>
    filter(
      time != 0
    ) |>
    ungroup() |>
    mutate(
      # convert exposure variable to a 0/1 binary variable, with 0 the reference level
      "{exposure}" := as.integer(as.factor(!!!exposure_syms)) - 1L
    ) |>
    pivot_wider(
      id_cols = all_of(c(subgroups, "time")),
      #names_glue = "{.value}_{exposure}", #but this doesn't work because scoping, quosures, something something
      names_from = c(!!!exposure_syms),
      names_sep = "_",
      values_from = c(
        n.risk, n.event, n.censor, n.competing,
        cml.nrisk, cml.event, cml.rate,
        cmlinc, cmlinc.se, cmlinc.low, cmlinc.high#,
        #rmst, rmst.low, rmst.high
      )
    ) |>
    mutate(
      n.nonevent_0 = n.risk_0 - n.event_0,
      n.nonevent_1 = n.risk_1 - n.event_1,
      
      ## cumulative quantities using information during time [0,t] (not just [t])
      
      # cumulative incidence rate ratio
      cmlirr = cml.rate_1 / cml.rate_0,
      cmlirr.ln.se = sqrt((1 / cml.event_0) + (1 / cml.event_1)),
      cmlirr.ll = exp(log(cmlirr) + qnorm(0.025) * cmlirr.ln.se),
      cmlirr.ul = exp(log(cmlirr) + qnorm(0.975) * cmlirr.ln.se),
      
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
}

if((length(exposure)>0) & contrast){
  data_contrasts_rounded <- contrast_aj(data_surv_rounded)
  
  ## output to disk ----
  arrow::write_feather(data_contrasts_rounded, fs::path(dir_output, glue("contrasts{filename_suffix}.arrow")))
  write_csv(data_contrasts_rounded, fs::path(dir_output, glue("contrasts{filename_suffix}.csv")))
}
