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

## import weights from matching or weighting method ----
data_weights <- read_feather(here_glue("output", "3-adjust", cohort, "{method}-{spec}", "data_adjusted.arrow"))


# append relevant characteristics and weights
data_balance <- 
  data_cohort |>
  select(patient_id, treatment, any_of(names(variable_labels))) |>
  left_join(
    data_weights |> select(patient_id, weight),
    by = "patient_id"
  )

# create table object using function in utility script
# this includes SDC
table_balance <- 
  data_balance |>
  mutate(
    N = 1L,
  ) |>
  table1_summary_smd(
    treatment = treatment,
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

