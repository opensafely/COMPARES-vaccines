
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Combine AJ and PLR estimates for all models within each cohort
# # # # # # # # # # # # # # # # # # # # #

# TODO: this may need to be modified to interact nicely with Airlock
# for example by combining at different levels

# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')
library("fs")
library("arrow")

## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  cohort <- "age75plus"

} else {
  cohort <- args[[1]]
}

cohort0 <- cohort

metaparams <-
  metaparams |>
  filter(cohort == cohort0)

output_dir <- here("output", "4-contrast", cohort, "contrasts")

dir_create(output_dir)

## combine AJ and PLR estimates ----

combine_model_outputs <- function(type, strategy, split){
  
    metaparams %>%
    mutate(
      data = pmap(
        .l = .,
        function(cohort, method, spec, subgroup, outcome, ...) {
          
          subgroup0 <- as.character(subgroup)
          
          dat <- 
            here("output", "4-contrast", cohort, glue("{method}-{spec}"), subgroup, outcome, strategy, glue("{type}.arrow")) |>
            read_feather() |>
            ungroup() %>%
            rename(all_of(c(subgroup_level = subgroup))) %>%
            mutate(
              subgroup_level_descr = fct_recoderelevel(subgroup_level, recoder[[subgroup0]]),
              .after = "subgroup_level"
            )
          dat
        }
      )
    ) |>
    unnest(data) #|>
    #nest(.by = {{split}})
  
}

estimates_plr <- combine_model_outputs("estimates", "plr")
contrasts_plr <- combine_model_outputs("contrasts", "plr")
estimates_aj <- combine_model_outputs("estimates", "aj")
contrasts_aj <- combine_model_outputs("contrasts", "aj")


## write to csv ----

write_split_data <- function(.data, dir, fileprefix, ...){
  .data |>
    group_by(...) |>
    group_walk(~ write_csv(.x, path(dir, glue("{fileprefix}{paste(t(.y),collapse='-')}.csv"))))
}

write_split_data(estimates_plr, output_dir, "estimates_plr_", outcome, subgroup)
write_split_data(contrasts_plr, output_dir, "contrasts_plr_", outcome, subgroup)
write_split_data(estimates_aj, output_dir, "estimates_aj_", outcome, subgroup)
write_split_data(contrasts_aj, output_dir, "contrasts_aj_", outcome, subgroup)


## move AJ and PLR plots to single folder ----

output_dir_plots <- path(output_dir, "plots")
dir_create(output_dir_plots)

metaparams |>
  mutate(
    plotdir = here("output", "4-contrast", cohort, glue("{method}-{spec}"), subgroup, outcome, "plr", "plot.png"),
    plotnewdir = path(output_dir_plots, glue("plr_{method}-{spec}_{subgroup}_{outcome}.png")),
  ) %>%
  {walk2(.$plotdir, .$plotnewdir, ~file_copy(.x, .y, overwrite = TRUE))}

metaparams |>
  mutate(
    plotdir = here("output", "4-contrast", cohort, glue("{method}-{spec}"), subgroup, outcome, "aj", "plot.png"),
    plotnewdir = path(output_dir_plots, glue("aj_{method}-{spec}_{subgroup}_{outcome}.png")),
  ) %>%
  {walk2(.$plotdir, .$plotnewdir, ~file_copy(.x, .y, overwrite = TRUE))}


## plot overall estimates for inspection ----

plot_contrasts <- function(data_contrasts, timeslice, method, spec, estimate, estimate.ll, estimate.ul, name){
  
  cohort0 <- cohort
  method0 <- method
  spec0 <- spec
  
  df <- data_contrasts |> 
  filter(cohort0==cohort, method0==method, spec0==spec, time==timeslice) |> 
  drop_na(subgroup_descr, outcome_descr, subgroup_level_descr, {{estimate}}, {{estimate.ll}}, {{estimate.ul}})

  
  # In case no data to plot
  if (nrow(df) == 0L) return(invisible(NULL)) 

  plot_temp <-   df  |>
    group_by(outcome_descr) |>
    ggplot(aes(y=subgroup_level_descr)) +
    geom_vline(aes(xintercept=0), linetype="dotted", colour="darkgrey")+
    geom_point(aes(x={{estimate}}))+
    geom_linerange(aes(xmin={{estimate.ll}}, xmax={{estimate.ul}}))+
    facet_grid(rows=vars(subgroup_descr), cols=vars(outcome_descr), scales="free", space="free_y", switch="y")+
    scale_x_continuous(expand = expansion(mult=c(0,0.01)))+
    labs(y=NULL)+
    theme_minimal()+
    theme(
      legend.position="bottom",
      axis.text.x.top=element_text(hjust=0),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement="outside",
      #strip.text.y.left = element_text(angle=0),
      strip.text.y.left = element_blank(),
      
      panel.border = element_blank(),
      panel.spacing = unit(0.3, "lines"),
    )
  
  ggsave(
    filename=path(
      path(output_dir_plots, glue("overall_plot_{method}-{spec}_{name}.png"))
    ),
    plot_temp,
    width=20, height=15, units="cm"
  )
  
  plot_temp
}

for (method in c("match", "weight")){
  for (spec in c("A", "B")){
    plot_contrasts(contrasts_plr, timeslice=90, method, spec, rd, rd.ll, rd.ul, "rd")
    plot_contrasts(contrasts_plr, timeslice=90, method, spec, rr, rr.ll, rr.ul, "rr")
    plot_contrasts(contrasts_aj, timeslice=90, method, spec, rd, rd.ll, rd.ul, "rd")
    plot_contrasts(contrasts_aj, timeslice=90, method, spec, rr, rr.ll, rr.ul, "rr")
  }
}
