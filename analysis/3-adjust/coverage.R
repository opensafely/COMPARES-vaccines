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
  method <- "lmw"
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

## import matching info ----
data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "{method}-{spec}", "data_adjusted.arrow"))


# append relevant characteristics to match data

data_balance <- 
  data_cohort |>
  select(patient_id, treatment, vax_date, any_of(names(variable_labels))) |>
  left_join(
    data_weights |> select(patient_id, weight),
    by = "patient_id"
  ) 

# weighting proportion
data_coverage <-
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
    vax_date = full_seq(.$vax_date, 1), # go X days before to
    fill = list(n=0)
  ) |>
  group_by(treatment) |>
  arrange(treatment, vax_date) |>
  mutate(
    cumul_unweighted = cumsum(unweighted),
    cumul_weighted = cumsum(weighted),
  )


data_coverage_rounded <-
  data_coverage |>
  group_by(treatment) |>
  mutate(
    cumul_unweighted = roundmid_any(cumul_unweighted, to = sdc.limit),
    cumul_weighted = roundmid_any(cumul_weighted, to = sdc.limit),
    unweighted = diff(c(0,cumul_unweighted)),
    weighted = diff(c(0,cumul_weighted)),
  )

write_csv(data_coverage_rounded, fs::path(output_dir, "data_coverage.csv"))


## plot matching coverage ----

xmin <- min(data_coverage$vax_date )
xmax <- max(data_coverage$vax_date )+1

plot_coverage_n <-
  data_coverage_rounded |>
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
    size=1
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
    y="Vaccines per day",
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
    cumul_unweighted = cumul_unweighted*((treatment*2) - 1),
    cumul_weighted = cumul_weighted*((treatment*2) - 1)
  ) |>
  ggplot()+
  geom_col(
    aes(
      x=vax_date+0.5,
      y=cumul_weighted,
      group=treatment_descr,
      fill=treatment_descr,
      colour=NULL
    ),
    alpha=0.4,
    width=1
  )+
  geom_linerange(
    aes(
      xmin=vax_date,
      xmax=vax_date+1,
      y=cumul_unweighted,
      group=treatment_descr,
      colour=treatment_descr
    ),
    size=1
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
    y="Cumulative vaccines per day (weighted and unweighted)",
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

ggsave(plot_coverage_cumuln, filename="coverage_cumulative.png", path=output_dir)

