# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe matching results
# imports matching data
# reports on matching coverage, matching flowcharts, creates a "table 1", etc
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
  cohort <- "age75plus"
  matchset <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  matchset <- args[[2]]
}



# create output directories ----

output_dir <- here_glue("output", "3-cohorts", cohort, "match{matchset}", "report")
fs::dir_create(output_dir)

## import unadjusted cohort data ----
data_cohort <- read_feather(here("output", "3-cohorts", cohort, "data_cohort.arrow"))

## import matching info ----
data_matches <- read_feather(here_glue("output", "3-cohorts", cohort, "match{matchset}", "data_matches.arrow"))


# append relevant characteristics to match data

data_balance <- 
  data_cohort |>
  select(patient_id, treatment, any_of(names(variable_labels))) |>
  left_join(
    data_matches |> select(patient_id, weight, matched),
    by = "patient_id"
  ) 





# matching coverage on each day of recruitment period ----


# matching coverage for boosted people
data_coverage <-
  data_matches |>
  mutate(eligible=1) |>
  group_by(treatment, vax_date) |>
  summarise(
    n_eligible = sum(eligible, na.rm=TRUE),
    n_matched = sum(matched, na.rm=TRUE),
  ) |>
  mutate(
    n_unmatched = n_eligible - n_matched,
  ) |>
  pivot_longer(
    cols = c(n_unmatched, n_matched),
    names_to = "status",
    names_prefix = "n_",
    values_to = "n"
  ) |>
  arrange(treatment, vax_date, status) |>
  group_by(treatment, vax_date, status) |>
  summarise(
    n = sum(n),
  ) |>
  group_by(treatment, status) %>%
  complete(
    vax_date = full_seq(.$vax_date, 1), # go X days before to
    fill = list(n=0)
  ) |>
  mutate(
    cumuln = cumsum(n)
  ) |>
  ungroup() |>
  mutate(
    status = factor(status, levels=c("unmatched", "matched")),
    status_descr = fct_recoderelevel(status, recoder$status)
  ) |>
  arrange(treatment, status_descr, vax_date)




data_coverage_rounded <-
  data_coverage |>
  group_by(treatment, status) |>
  mutate(
    cumuln = roundmid_any(cumuln, to = sdc.limit),
    n = diff(c(0,cumuln)),
  )

write_csv(data_coverage_rounded, fs::path(output_dir, "data_coverage.csv"))



## plot matching coverage ----

xmin <- min(data_coverage$vax_date )
xmax <- max(data_coverage$vax_date )+1

plot_coverage_n <-
  data_coverage |>
  mutate(
    treatment_descr = fct_recoderelevel(as.character(treatment), recoder$treatment),
    n=n*((treatment*2) - 1)
  ) |>
  ggplot()+
  geom_col(
    aes(
      x=vax_date+0.5,
      y=n,
      group=paste0(treatment,status),
      fill=treatment_descr,
      alpha=fct_rev(status),
      colour=NULL
    ),
    position=position_stack(reverse=TRUE),
    #alpha=0.8,
    width=1
  )+
  #geom_rect(xmin=xmin, xmax= xmax+1, ymin=-6, ymax=6, fill="grey", colour="transparent")+
  geom_hline(yintercept = 0, colour="black")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_coverage$vax_date, "1 month")),
    limits = c(xmin-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    #labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.x)),
    expand = expansion(c(0, NA))
  )+
  scale_fill_brewer(type="qual", palette="Set2")+
  scale_colour_brewer(type="qual", palette="Set2")+
  scale_alpha_discrete(range= c(0.8,0.4))+
  labs(
    x="Date",
    y="Booster vaccines per day",
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

plot_coverage_n

ggsave(plot_coverage_n, filename="coverage_count.png", path=output_dir)

plot_coverage_cumuln <-
  data_coverage |>
  mutate(
    treatment_descr = fct_recoderelevel(as.character(treatment), recoder$treatment),
    cumuln=cumuln*((treatment*2) - 1)
  ) |>
  ggplot()+
  geom_col(
    aes(
      x=vax_date+0.5,
      y=cumuln,
      group=paste0(treatment,status),
      fill=treatment_descr,
      alpha=fct_rev(status),
      colour=NULL
    ),
    position=position_stack(reverse=TRUE),
    width=1
  )+
  geom_rect(xmin=xmin, xmax= xmax+1, ymin=-6, ymax=6, fill="grey", colour="transparent")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_coverage$vax_date, "1 month")),
    limits = c(xmin-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    #labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.)),
    expand = expansion(c(0, NA))
  )+
  scale_fill_brewer(type="qual", palette="Set2")+
  scale_colour_brewer(type="qual", palette="Set2")+
  scale_alpha_discrete(range= c(0.8,0.4))+
  labs(
    x="Date",
    y="Cumulative booster vaccines",
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

plot_coverage_cumuln

ggsave(plot_coverage_cumuln, filename="coverage_stack.png", path=output_dir)

