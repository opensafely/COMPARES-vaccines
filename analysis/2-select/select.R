# # # # # # # # # # # # # # # # # # # # #
# Purpose: import processed data and filter out people who are excluded from the main analysis
# outputs:
#  - inclusion/exclusions flowchart data (up to matching step)
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library("arrow")
library('gt')
library('gtsummary')


## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

## Import redaction functions
#source(here("analysis", "0-lib", "redaction.R"))


## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "age65plus" #currently `age65plus` or `cv`
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
}

# derive subgroup info
cohort_sym <- sym(cohort)

## create output directories for data ----
output_dir <- here("output", "2-select", cohort)
fs::dir_create(output_dir)


# Main ----

## Import processed data ----

data_prepared <- read_feather(here("output", "2-select", "data_prepared.arrow"))


# Report total number of people vaccinated in time period, by vaccine type ----
# This output ignores cohort, so is the same across different cohorts
# but it is more lightweight than outputting in the `data_process` script, because then
# the release action will need to copy in the entire processed data set

## unrounded totals
total_n_unrounded <-
  bind_rows(
    tibble(
      vax_product="any", 
      n=nrow(data_prepared)
      ),
    count(
      data_prepared |> 
        mutate(
          boost_type=fct_other(
            vax_product, 
            keep=treatment_lookup$treatment, 
            other_level="other")
          ), 
      vax_product, 
      .drop=FALSE)
  ) |>
  mutate(
    pct = n/first(n)
  )
write_csv(total_n_unrounded, fs::path(output_dir, "total_allcohorts_unrounded.csv"))

## rounded totals
total_n_rounded <-
  total_n_unrounded |>
  mutate(
    n= ceiling_any(n, sdc.limit),
    pct = n/first(n)
  )
write_csv(total_n_rounded, fs::path(output_dir, "total_allcohorts_rounded.csv"))


# get dates between which both vaccines were being used

vax_date_limits <-
  data_prepared |>
  summarise(
    min =vax_date(vax_date),
    max = max(vax_date),
    .by = vax_product
  )

both_products_used_from_date <- max(vax_date_limits$min)
both_products_used_to_date <- min(vax_date_limits$max)


## Define selection criteria ----
data_criteria <-
  data_prepared |>
  filter(!!cohort_sym) |>
  transmute(
    patient_id,
    vax_date,
    vax_product,
    prior_vax_interval_atleast12weeks,
    vax_product_of_interest = vax_product %in% c(productA, productB),
    during_period_both_products_used = between(
      vax_date,
      both_products_used_from_date,
      both_products_used_to_date
    ),
    codamin_flu = !is.na(flu_vaccine_same_date),
    select_coadmin_only_if_required = !(restrict_to_coadmin_flu & !codamin_flu), #if restrict_to_coadmin_flu==TRUE and codamin_flu==FALSE then 0, otherwise 1
    # no_prior_productA = !vaxhist_productA,
    # no_prior_productB = !vaxhist_productB,
    prior_vax_1plus = prior_vax_count >= 1,
    has_age = !is.na(age_eligible),
    has_sex = !is.na(sex) & !(sex %in% c("intersex", "unknown")),
    has_imd = !is.na(imd_Q5),
    has_region = !is.na(region),
    # has_ethnicity5 = !is.na(ethnicity5),
    # has_msoa = !is.na(msoa),
    # isnot_hscworker = !hscworker,
    # isnot_carehomeresident = !care_home_combined,
    # isnot_endoflife = !endoflife,
    # isnot_housebound = !housebound,
    # has_norecentcovid = ((vax_date - anycovid_0_date) >= 28) | is.na(anycovid_0_date),
    isnot_inhospital = !inhospital,

    include_regardless_of_flu_coadmin = (prior_vax_interval_atleast12weeks &
      vax_product_of_interest &
      during_period_both_products_used &
      # no_prior_productA &
      # no_prior_productB &
      prior_vax_1plus &
      has_age &
      has_sex &
      has_imd &
      has_region &
      # has_ethnicity5 &
      # isnot_hscworker &
      # isnot_endoflife &
      # has_norecentcovid &
      isnot_inhospital),

    include = include_regardless_of_flu_coadmin &
      select_coadmin_only_if_required
  )


data_cohort <- 
  data_criteria |>
  filter(include) |>
  select(patient_id) |>
  left_join(data_prepared, by="patient_id") |>
  droplevels()

write_feather(data_cohort, fs::path(output_dir, "data_cohort.arrow"))


## report unadjusted balance for baseline variables between levels  of treatment ----
# This output includes SDC

table_cohort <-
  data_cohort |>
  mutate(
    N = 1L,
    weight = 1L
  ) |>
  table1_summary_smd(
    group = treatment,
    weight = weight,
    label = variable_labels,
    threshold = sdc.limit
  )

write_csv(table_cohort, fs::path(output_dir, "table_cohort.csv"))

## report unadjusted outcome rates for baseline variables between levels of treatment ----
# This output includes SDC

outcome_labels <-
  c(
    N = "Total N",
    set_names(events_lookup$event_descr, events_lookup$event_var)
  ) |>
  as.list()

table_outcomes <-
  data_cohort |>
  select(
    treatment,
    vax_date,
    any_of(events_lookup$event_var)
  ) |>
  mutate(
    N = 1L,
    weight = 1L,
  ) |>
  mutate(
    across(
      any_of(events_lookup$event_var),
      \(x) {
        coalesce(x <= (vax_date + maxfup), FALSE)
      }
    )
  ) %>%
  table1_summary(
    group = treatment,
    label = outcome_labels,
    threshold = sdc.limit
  )


write_csv(table_outcomes, fs::path(output_dir, "table_outcomes.csv"))


remove(data_cohort)

## output simple dataset containing exclusions criteria met ----

data_inclusioncriteria <- data_criteria |>
  transmute(
    patient_id,
    vax_product,
    c0 = TRUE,
    c1 = c0 & vax_product_of_interest,
    # Apply the same-day influenza coadministration criterion only when
    c2 = c1 & during_period_both_products_used,
    c3 = c2 & select_coadmin_only_if_required,
    c4 = c3 & prior_vax_interval_atleast12weeks, #& no_prior_productA & no_prior_productB,
    c5 = c4 & prior_vax_1plus,
    c6_1 = c5 & (has_age & has_sex & has_imd & has_region)
    # c6_3 = c6 & (isnot_endoflife),
    # c6_4 = c6 & (has_norecentcovid),
    # c6_5 = c6 & (isnot_inhospital),
    #c4 = c4_1 & c4_2 & c4_3 & c4_4 & c4_5
  ) 

write_feather(data_inclusioncriteria, sink = fs::path(output_dir, "data_inclusioncriteria.arrow"))


## Create flowchart ----

create_flowchart <- function(round_level = 1){
  
  flowchart_output <-
    data_inclusioncriteria |>
    select(-patient_id, -c0) |>
    filter(c1) |>
    group_by(vax_product) |>
    summarise(
      across(
        .cols = everything(),
        .fns = ~ ceiling_any(sum(., na.rm = TRUE), round_level)
      )
    ) |>
    pivot_longer(
      cols = -c(vax_product),
      names_to = "criteria",
      values_to = "n"
    ) |>
    # When influenza coadministration is not an inclusion criterion,
    # omit c2 from the flowchart rather than showing a redundant step.
    filter(!(criteria == "c2" & !restrict_to_coadmin_flu)) |>
    mutate(
      level = if_else(str_detect(criteria, "c\\d+$"), 1, 2),
      n_level1 = if_else(level == 1, n, NA_real_),
      n_level1_fill = n_level1
    ) |>
    fill(n_level1_fill) |>
    group_by(vax_product) |>
    mutate(
      n_exclude = lag(n_level1_fill) - n,
      pct_exclude = n_exclude / lag(n_level1_fill),
      pct_all = n_level1 / first(n),
      pct_step = n_level1 / lag(n_level1_fill),
      #crit = str_extract(criteria, "^c\\d+"),
      crit = criteria,
      criteria = fct_case_when(
        crit == "c1" ~ "Received COVID-19 vaccine during Autumn 2024 campaign",
        crit == "c2" ~ "  during period where both vaccines in use",
        crit == "c3" ~ "  with influenza vaccine on the same day",
        crit == "c4" ~ "  with no prior Covid-19 vaccine within 12 weeks",
        crit == "c5" ~ "  with at least 1 prior COVID-19 vaccine dose",
        crit == "c6_1" ~ "    no missing demographic information",
        crit == "c6_2" ~ "    not a health and social care worker",
        crit == "c6_3" ~ "    not end-of-life",
        crit ==
          "c6_4" ~ "    no documented COVID-19 infection/disease within prior 28 days",
        crit == "c6_5" ~ "    not admitted in hospital at time of booster",
        crit == "c7" ~ "  included in matching run",
        TRUE ~ "NA_character_boop" # should not appear
      )
    )
  
  return(flowchart_output)
}

## unrounded flowchart
data_flowchart <- create_flowchart(1)
write_feather(data_flowchart, fs::path(output_dir, "flowchart.arrow"))
#write_csv(data_flowchart, here("output", "data", "flowchart.csv"))

## rounded flowchart
data_flowchart_rounded <- create_flowchart(sdc.limit)
write_feather(data_flowchart_rounded, fs::path(output_dir, "flowchart_rounded.arrow"))
write_csv(data_flowchart_rounded, fs::path(output_dir, "flowchart_rounded.csv"))

## unrounded totals
total_n_unrounded <-
  bind_rows(
    tibble(vax_product="any", n=nrow(data_inclusioncriteria)),
    count(data_inclusioncriteria |> mutate(boost_type=fct_other(vax_product, keep=treatment_lookup$treatment, other_level="other")), vax_product, .drop=FALSE)
  ) |>
  mutate(
    pct = n/first(n)
  )
write_csv(total_n_unrounded, fs::path(output_dir, "total_unrounded.csv"))

## rounded totals
total_n_rounded <-
  total_n_unrounded |>
  mutate(
    n= ceiling_any(n, sdc.limit),
    pct = n/first(n)
  )
write_csv(total_n_rounded, fs::path(output_dir, "total_rounded.csv"))

## remove large in-memory objects
remove(data_inclusioncriteria)

##############################################
# Flu coadministration
##############################################

## Flu coadministration timing by vaccine product --------------------------

# Create a long table with one row per patient and flu vaccination timing window
flu_coadmin_days <-
  data_criteria |>
  filter(include_regardless_of_flu_coadmin) |>
  select(patient_id) |>
  left_join(
    data_prepared |>
      select(
        patient_id,
        vax_product,
        flu_vaccine_before_30_days,
        flu_vaccine_same_day_days,
        flu_vaccine_after_30_days
       ),
    by = "patient_id"
  ) |>
  pivot_longer(
    cols = starts_with("flu_vaccine_"),
    names_to = "flu_window",
    values_to = "days"
  ) |>
  filter(!is.na(days))

# Count flu coadministration events by COVID-19 vaccine product and day
# difference, applying SDC rounding.
flu_coadmin_days_by_vaccine <-
  flu_coadmin_days |>
  group_by(vax_product, days) |>
  summarise(
    n_roundmid6 = roundmid_any(n(), sdc.limit),
    .groups = "drop"
  ) |>
  arrange(vax_product, days)

# Count flu coadministration events across all COVID-19 vaccine products.
flu_coadmin_days_general <-
  flu_coadmin_days |>
  group_by(days) |>
  summarise(
    vax_product = "All vaccines",
    n_roundmid6 = roundmid_any(n(), sdc.limit),
    .groups = "drop"
  ) |>
  arrange(days)

# Combine product-specific counts with the overall count.
flu_coadmin_days_table <-
  bind_rows(
    flu_coadmin_days_general,
    flu_coadmin_days_by_vaccine
  )
write_csv(flu_coadmin_days_table, fs::path(output_dir, "flu_coadmin_days_table.csv"))

# Plot the distribution of days between COVID-19 and flu vaccination.
flu_coadmin_days_plot <- ggplot(flu_coadmin_days_table, aes(x = days, y = n_roundmid6)) +
  geom_col() +
  facet_wrap(~ vax_product, ncol = 1, scales = "free_y") +
  labs(
    x = "Days between COVID-19 and flu vaccination",
    y = "Number of patients",
    title = "Flu co-vaccination timing"
  ) +
  theme_minimal()
flu_coadmin_days_plot

ggsave(filename = fs::path(output_dir, glue("flu_coadmin_days_plot.png")), 
       flu_coadmin_days_plot, 
       width = 20, 
       height = 20, 
       units = "cm")

## Flu coadministration by subgroup ----------------------------------------

# Create a long-format subgroup table.
# Each patient contributes one row per subgroup variable.
# flu_coadmin is classified hierarchically as:
#   - same_day
#   - 7_days
#   - no_covax

flu_coadmin_subgroups <-
  data_criteria |>
  filter(include_regardless_of_flu_coadmin) |>
  select(patient_id) |>
  left_join(
    data_prepared |>
      select(
        patient_id,
        ethnicity5,
        ethnicity16,
        sex,
        ageband,
        imd_Q5,
        vax_product,
        flu_coadmin_7_day,
        flu_coadmin_same_day
      ),
    by = "patient_id"
  ) |>
  mutate(
    flu_coadmin = case_when(
      flu_coadmin_same_day ~ "same_day",
      flu_coadmin_7_day ~ "7_days",
      !flu_coadmin_7_day ~ "no_covax",
      TRUE ~ NA_character_
    )
  ) |>
  pivot_longer(
    cols = c(ethnicity5, ethnicity16, sex, imd_Q5, ageband),
    names_to = "group",
    values_to = "category",
    values_transform = list(category = as.character)
  )

# Count patients by vaccine product, subgroup, subgroup category,
# and flu coadministration category.
counts_by_vaccine <- flu_coadmin_subgroups |>
  count(
    vax_product,
    group,
    category,
    flu_coadmin,
    name = "n_roundmid6"
  )

# Count patients across all vaccine products
counts_all_products <- flu_coadmin_subgroups |>
  count(
    group,
    category,
    flu_coadmin,
    name = "n_roundmid6"
  ) |>
  mutate(vax_product = "All products") |>
  select(
    vax_product,
    group,
    category,
    flu_coadmin,
    n_roundmid6
  )

# Combine product-specific and overall counts.
# Percentages describe the distribution of flu coadministration status
# within each subgroup category.
flu_coadmin_subgroups_table <-
  bind_rows(
    counts_by_vaccine,
    counts_all_products
  )  |>
  group_by(vax_product, group, category) |>
  mutate(
    total_roundmid6 = roundmid_any(sum(n_roundmid6), sdc.limit),
    n_roundmid6 = roundmid_any(n_roundmid6, sdc.limit),
    pct_roundmid6 = round(n_roundmid6 / total_roundmid6 * 100, 1)
  ) |>
  ungroup()|>
  mutate(
    flu_coadmin = factor(
      flu_coadmin,
      levels = c("same_day", "7_days", "no_covax")
    ),
    group = factor(
      group,
      levels = c("ageband", "sex", "ethnicity5", "ethnicity16", "imd_Q5")
    )
  )

write_csv(flu_coadmin_subgroups_table, fs::path(output_dir, "flu_coadmin_subgroups_table.csv"))
# Plot the percentage of patients in each flu coadministration category,
# by subgroup category and vaccine product.
flu_coadmin_subgroups_plot <- ggplot(flu_coadmin_subgroups_table, aes(
  x = pct_roundmid6,
  y = category,
  colour = flu_coadmin
)) +
  geom_point(size = 2) +
  facet_grid(
    rows = vars(group),
    cols = vars(vax_product),
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    x = "%",
    y = NULL,
    colour = "Flu coadministration"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

ggsave(filename = fs::path(output_dir, glue("flu_coadmin_subgroups_plot.png")), 
       flu_coadmin_subgroups_plot, 
       width = 20, 
       height = 20, 
       units = "cm")