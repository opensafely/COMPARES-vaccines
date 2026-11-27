# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe balance after adjustment via matching or weighting
# imports weights from matching or weighting and procudes:
# - a table1-style dataset
# - a love plot to check balance in the server
# - a table containing effective sample size due to the weights
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('survival')
library('cobalt')
library('gt')
library('gtsummary')

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "cv"
  method <- "match"
  spec <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  method <- args[[2]]
  spec <- args[[3]]
}



# create output directories ----

output_dir <- here_glue("output", "3-adjust", cohort, "{method}-{spec}", "report")
fs::dir_create(output_dir)

## import unadjusted cohort data ----
data_cohort <- read_feather(here("output", "2-select", cohort, "data_cohort.arrow"))

## import weights from matching or weighting method ----
data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "{method}-{spec}", "data_adjusted.arrow"))


# append relevant characteristics and weights
data_balance <-
  data_cohort |>
  select(patient_id, treatment, vax_date, any_of(variables_lookup$variable)) |>
  left_join(
    data_weights |> select(patient_id, weight),
    by = "patient_id"
  )

## distribution of balancing weights ----
# output distrubtion of balancing weights

plot_weight_distr <- 
  data_balance |>
  mutate(
    # covert to log weight for symmetric plotting, except for matching which has 0/1 support
    weight = if(method != "match") {log(weight)} else {weight}
  ) |>
  ggplot()+
  geom_histogram(aes(x=weight))+
  facet_grid(
    rows=vars(treatment)
  )+
  theme_bw()

ggsave(plot_weight_distr, filename="plot_weight_distr.png", path=output_dir)

table_weight_distr <-
  data_balance |>
  summarise(
    n = n(),
    n_over_1 = sum(weight >= 1),
    n_under_1 = sum(weight < 1),
    mean = mean(weight),
    median = median(weight),
    min = min(weight),
    p05 = quantile(weight, p=0.05),
    p10 = quantile(weight, p=0.10),
    p25 = quantile(weight, p=0.25),
    p75 = quantile(weight, p=0.75),
    p90 = quantile(weight, p=0.90),
    p95 = quantile(weight, p=0.95),
    max = max(weight),
    .by = c(treatment)
  )

write_csv(table_weight_distr, fs::path(output_dir, "table_weight_distr.csv"))

# create table object using function in utility script
# this includes SDC
table_balance <-
  data_balance |>
  mutate(
    N = 1L,
  ) |>
  table1_summary_smd(
    group = treatment,
    weight = weight,
    label = variable_labels,
    threshold = sdc.limit
  )

write_csv(table_balance, fs::path(output_dir, "table_balance.csv"))


## love / smd plot ----
# output these to disk for checking, but this will be done properly outside of the server
## across all method-spec combinations

plot_balance <-
  table_balance |>
  mutate(
    variable_card = as.numeric(variable)%%2,
    variable_level = replace_na(as.character(variable_level), ""),
    level_pre = case_when(
      context %in% c("dichotomous", "continuous") ~ variable_label,
      TRUE ~ paste(variable_label, variable_level, sep=": ")
    ),
    level = fct_rev(fct_inorder(str_replace(level_pre, "\\:\\s$", ""))),
    cardn = row_number()
  ) |>
  ggplot() +
  geom_point(aes(x=smd, y=level))+
  geom_rect(aes(alpha = variable_card, ymin = rev(cardn)-0.5, ymax =rev(cardn+0.5)), xmin = -Inf, xmax = Inf, fill='grey', colour="transparent") +
  scale_alpha_continuous(range=c(0,0.3), guide="none")+
  labs(
    x="Standardised mean difference",
    y=NULL,
    alpha=NULL
  )+
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill="transparent", colour="transparent"),
    strip.text.y.left = element_text(angle = 0, hjust=1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0, "lines")
  )

ggsave(plot_balance, filename="plot_balance.png", path=output_dir)


## distribution of date of vaccination ----
# output distrubtion of vaccination date by product type, 
# both unweighted (colour shade)
# and weighted by the balancing weights (grey shade)


# weighting proportion
data_vax_dates_distr <-
  data_balance |>
  mutate(eligible=1, weight) |>
  group_by(treatment, vax_date) |>
  summarise(
    unweighted = n(),
    weighted = sum(weight, na.rm=TRUE),
  ) |>
  mutate(
    diff = weighted - unweighted,
    ratio = weighted / unweighted,
  ) |>
  group_by(treatment) %>%
  complete(
    vax_date = full_seq(.$vax_date, 1), 
    fill = list(n=0, unweighted=0, weighted=0)
  ) |>
  group_by(treatment) |>
  arrange(treatment, vax_date) |>
  mutate(
    cumul_unweighted = cumsum(unweighted),
    cumul_weighted = cumsum(weighted),
  )


data_vax_dates_distr_rounded <-
  data_vax_dates_distr |>
  group_by(treatment) |>
  mutate(
    cumul_unweighted = roundmid_any(cumul_unweighted, to = sdc.limit),
    cumul_weighted = roundmid_any(cumul_weighted, to = sdc.limit),
    unweighted = diff(c(0,cumul_unweighted)),
    weighted = diff(c(0,cumul_weighted)),
  )

write_csv(data_vax_dates_distr_rounded, fs::path(output_dir, "data_vax_dates_distr.csv"))


## plot  vax_dates_distr ----

xmin <- min(data_vax_dates_distr$vax_date )
xmax <- max(data_vax_dates_distr$vax_date )+1

plot_vax_date_n <-
  data_vax_dates_distr_rounded |>
  mutate(
    treatment_descr = fct_recoderelevel(as.character(treatment), recoder$treatment),
    unweighted = unweighted*((treatment*2) - 1),
    weighted = weighted*((treatment*2) - 1)
  ) |>
  ggplot()+
  geom_col(
    aes(
      x=vax_date+0.5,
      y=weighted,
      group=treatment_descr,
      fill=treatment_descr,
      colour=NULL
    ),
    alpha=0.4,
    width=1
  )+
  geom_col(
    aes(
      x=vax_date+0.5,
      y=unweighted,
      group=treatment_descr,
      colour=NULL
    ),
    fill="grey",
    alpha=0.5,
    width=1
  )+
  geom_linerange(
    aes(
      xmin=vax_date,
      xmax=vax_date+1,
      y=unweighted,
      group=treatment_descr,
      colour=treatment_descr
    ),
    linewidth=1
  )+
  geom_hline(yintercept = 0, colour="black")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_vax_dates_distr$vax_date, "1 month")),
    limits = c(xmin-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.x)),
    expand = expansion(c(0, NA))
  )+
  scale_fill_brewer(type="qual", palette="Set2")+
  scale_colour_brewer(type="qual", palette="Set2")+
  scale_alpha_discrete(range= c(0.8,0.4))+
  labs(
    x="Date",
    y="Vaccinations per day",
    colour=NULL,
    fill=NULL,
    alpha=NULL
  ) +
  theme_minimal()+
  theme(
    axis.line.x.bottom = element_line(),
    axis.text.x.top=element_text(hjust=0),
    strip.text.y.right = element_text(angle = 0),
    axis.ticks.x=element_line(),
    legend.position = "bottom"
  )+
  NULL

plot_vax_date_n

ggsave(plot_vax_date_n, filename="vax_dates_distr_count.png", path=output_dir)
  
