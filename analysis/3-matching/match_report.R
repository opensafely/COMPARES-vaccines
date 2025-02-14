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



# Balance table ----

balance <- 
  data_cohort |>
  select(
    any_of(names(variable_labels))
  ) |>
  bal.tab(
    weights = data_matches$weight,
    treat = data_matches$treatment,
    stats = c("mean.diffs"), 
    disp = c("m", "sd"),
    un = TRUE,
    s.d.denom = "pooled", 
    binary = "std", 
    continuous = "std"
  ) %>%
  `[[`("Balance")
  

write_csv(balance, fs::path(output_dir, "data_balance.csv"))



# append all characteristics to match data
data_matched_baseline <-
  data_cohort |>
  filter(patient_id %in% data_matches$patient_id) |>
  select(patient_id, any_of(names(variable_labels))) %>%
  left_join(
    data_matches |> filter(matched),
    .,
    by="patient_id"
  )




# table 1 style baseline characteristics for the adjusted pseudo population ----
# for matching, this is just the characteristics in the matched population

## this is essentially equivalent to:
# table1_summary(
#   data = data_for_matched_poopulation_only,
#   group = treatment_descr, 
#   labels = variable_labels[names(variable_labels) %in% names(.)], 
#   threshold = sdc.limit
# )
# But we don't use this approach because we get consistency with the approach used in the weighting script 

## survey for weighted table1 --------

tab_summary_adjusted <-
  data_cohort |>
  left_join(
    data_matches |> select(patient_id, weight),
    by = "patient_id"
  ) |>
  select(
    treatment, 
    weight,
    any_of(names(variable_labels))
  ) |>
  mutate(
    N = 1L,
    treatment_descr = fct_recoderelevel(as.character(treatment), recoder$treatment),
  ) %>%
  table1_svysummary(
    group = treatment_descr, 
    weight = weight,
    labels = variable_labels[names(variable_labels) %in% names(.)], 
    threshold = sdc.limit
  )

write_csv(tab_summary_adjusted, fs::path(output_dir, "table1.csv"))


# Balance table ----

balance <- 
  data_cohort |>
  select(
    any_of(names(variable_labels))
  ) |>
  bal.tab(
    weights = data_matches$weight,
    treat = data_matches$treatment,
    #distance = data_matches$ps,
    stats = c("mean.diffs"), 
    disp = c("m", "sd"),
    un = TRUE,
    s.d.denom = "pooled", 
    binary = "std", 
    continuous = "std"
  ) %>%
  `[[`("Balance")

write_csv(balance, fs::path(output_dir, "data_balance.csv"))

# TODO: fix variable names in balance data so that they can be joined on the tbl_summary datasets

# # love / smd plot ----

# TODO: replace this with output from cobalt 
# TODO: copy this to the weighted script

data_smd <- 
  tab_summary_baseline |>
  filter(
    !(variable %in% c("treatment_descr", "N"))
  ) |>
  group_by(variable, variable_label, variable_level) |>
  mutate(
    mean = coalesce(mean,p),
    sd = coalesce(sd,sqrt(p*(1-p)))
  ) |>
  summarise(
    diff = diff(mean),
    sd = sqrt(mean(sd^2)),
    smd = diff/sd,
  ) |>
  ungroup() %>%
  mutate(
    variable = droplevels(factor(variable_label, levels=map_chr(variable_labels, ~.))),
  ) |>
  arrange(variable) |>
  mutate(
    variable_card = as.numeric(variable)%%2,
    variable_level = replace_na(as.character(variable_level), ""),
    level = fct_rev(fct_inorder(str_replace(paste(variable, variable_level, sep=": "),  "\\:\\s$", ""))),
    cardn = row_number()
  )

plot_smd <-
  ggplot(data_smd)+
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

write_csv(data_smd, fs::path(output_dir, "data_smd.csv"))
ggsave(plot_smd, filename="plot_smd.png", path=output_dir)

















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


# flowchart ----

create_flowchart <- function(flowchart, threshold){

  flowchart <-
    data_matches |>
    select(patient_id, treatment, matched) |>
    mutate(
      vax_product = case_when(
        treatment==0L ~ productA,
        treatment==1L ~ productB,
        TRUE ~ NA_character_
      ),
    ) |>
    group_by(vax_product) |>
    summarise(
      n = roundmid_any(sum(matched), threshold),
      n_lag = roundmid_any(n(), threshold),
      n_level1 = n,
      n_level1_fill = n,
      n_exclude = n_lag - n,
      pct_exclude = n_exclude/n(),
      #pct_all = n_level1 / first(n),
      pct_step = n / n_lag,
      level=1,
      crit = "c5",
      criteria = factor("  and successfully matched")
    ) |>
    ungroup() %>%
    bind_rows(
      flowchart,
      .
    ) |>
    arrange(
      vax_product, criteria
    ) |>
    group_by(vax_product) |>
    mutate(
      pct_all = if_else(is.na(pct_all) & level==1, n_level1 / first(n_level1), pct_all)
    )
}


flowchart0 <- read_feather(here("output", "3-cohorts", cohort, "flowchart.arrow"))
flowchart0_rounded <- read_feather(here("output", "3-cohorts", cohort, "flowchart_rounded.arrow"))


data_flowchart_match <- create_flowchart(flowchart0, 1)
data_flowchart_match_rounded <- create_flowchart(flowchart0_rounded, sdc.limit)

write_csv(data_flowchart_match, fs::path(output_dir, "flowchart.csv"))
write_csv(data_flowchart_match_rounded, fs::path(output_dir, "flowchart_rounded.csv"))

