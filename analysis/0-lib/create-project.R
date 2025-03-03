######################################
# Purpose: 
# create the project-yaml file 
######################################


# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('yaml')
library('here')
library('glue')


## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))


## restrict metaparams to those current available:
metaparams <- 
  metaparams |>
  select(cohort, method, spec) |>
  unique() |>
  filter(
    spec == "A"
  )

# create action functions ----

## create comment function ----
comment <- function(...){
  list_comments <- list(...)
  comments <- map(list_comments, ~paste0("## ", ., " ##"))
  comments
}

#splice <- function(...) {list_flatten(lst(...))}

## create function to convert comment "actions" in a yaml string into proper comments
convert_comment_actions <-function(yaml.txt){
  yaml.txt %>%
    str_replace_all("\\\n(\\s*)\\'\\'\\:(\\s*)\\'", "\n\\1")  %>%
    str_replace_all("([^\\'])\\\n(\\s*)\\#\\#", "\\1\n\n\\2\\#\\#") %>%
    str_replace_all("\\#\\#\\'\\\n", "\n")
}


## generic action function ----
action <- function(
  name,
  run,
  arguments=NULL,
  needs=NULL,
  highly_sensitive=NULL,
  moderately_sensitive=NULL,
  ... # other arguments / options for special action types
){

  outputs <- list(
    highly_sensitive = highly_sensitive,
    moderately_sensitive = moderately_sensitive
  )
  outputs[sapply(outputs, is.null)] <- NULL

  
  if(!is.null(names(arguments))){
    arguments <- paste0("--", names(arguments), "=", arguments)
  }
  
  action <- list(
    run = paste(c(run, arguments), collapse=" "),
    needs = needs,
    outputs = outputs,
    ... = ...
  )
  action[sapply(action, is.null)] <- NULL

  action_list <- list(name = action)
  names(action_list) <- name

  action_list
}


## selection action function ----

action_selection <- function(cohort){
  action(
    name = glue("data_selection_{cohort}"),
    run = "r:v2 analysis/2-prepare/data_selection.R",
    arguments = c(cohort),
    needs = list("prepare"),
    highly_sensitive = lst(
      arrow = glue("output/2-prepare/{cohort}/*.arrow"),
    ),
    moderately_sensitive = lst(
      csv = glue("output/2-prepare/{cohort}/*.csv"),
    )
  )
}

## match actions function ----
action_match <- function(cohort, spec){

  splice(

    action(
      name = glue("match_{cohort}_{spec}"),
      run = "r:v2 analysis/3-adjust/match.R",
      arguments = c(cohort, spec),
      needs = list(glue("data_selection_{cohort}")),
      highly_sensitive = lst(
        arrow = glue("output/3-adjust/{cohort}/match-{spec}/*.arrow")
      )
    ),

    action(
      name = glue("match_report_{cohort}_{spec}"),
      run = "r:v2 analysis/3-adjust/report.R",
      arguments = c(cohort, "match", spec),
      needs = list(glue("data_selection_{cohort}"),  glue("match_{cohort}_{spec}")),
      # highly_sensitive = lst(
      #   arrow = glue("output/3-adjust/{cohort}/match-{spec}/report/*.arrow"),
      # ),
      moderately_sensitive = lst(
        csv = glue("output/3-adjust/{cohort}/match-{spec}/report/*.csv"),
        png = glue("output/3-adjust/{cohort}/match-{spec}/report/*.png")
      )
    )
  )


}


## match actions function ----
action_weight <- function(cohort, spec){
  
  splice(
    
    action(
      name = glue("weight_{cohort}_{spec}"),
      run = "r:v2 analysis/3-adjust/weight.R",
      arguments = c(cohort, spec),
      needs = list(glue("data_selection_{cohort}")),
      highly_sensitive = lst(
        arrow = glue("output/3-adjust/{cohort}/weight-{spec}/*.arrow")
      )
    ),
    
    action(
      name = glue("weight_report_{cohort}_{spec}"),
      run = "r:v2 analysis/3-adjust/report.R",
      arguments = c(cohort, "weight", spec),
      needs = list(glue("data_selection_{cohort}"),  glue("weight_{cohort}_{spec}")),
      # highly_sensitive = lst(
      #   arrow = glue("output/3-adjust/{cohort}/weight-{spec}/report/*.arrow"),
      # ),
      moderately_sensitive = lst(
        csv = glue("output/3-adjust/{cohort}/weight-{spec}/report/*.csv"),
        png = glue("output/3-adjust/{cohort}/weight-{spec}/report/*.png")
      )
    )
  )
}

action_combine_weights <- function(cohort){
  
  action(
    name = glue("combine_weights_{cohort}"),
    run = glue("r:v2 analysis/3-adjust/combine.R {cohort}"),
    needs = list(
      glue("data_selection_{cohort}"),
      glue("match_{cohort}_A"),
      glue("weight_{cohort}_A")#'
      #  more actions here
    ),
    highly_sensitive = lst(
      arrow = glue("output/3-adjust/{cohort}/combine/*.arrow"),
    )
  )
}

## get km or ci actions function ----
action_km_contrast <- function(
  cohort, method, spec, subgroup, outcome
){
    dir_output <- glue("output/4-contrast/{cohort}/{method}-{spec}/{subgroup}/{outcome}/km/")
  
    ## kaplan-meier action
    action(
      name = glue("km_{cohort}_{method}_{spec}_{subgroup}_{outcome}"),
      run = "kaplan-meier-function:v0.0.10",
      #arguments = c(cohort, method, spec, subgroup, outcome),
      arguments = c(
        "df_input" = glue("output/3-adjust/{cohort}/combine/data_weights.arrow"), 
        "dir_output" = dir_output,
        "exposure" = "treatment",
        "subgroups" = glue("{subgroup}"),
        "origin_date" = "vax_date",
        "event_date" = glue("{outcome}_date"),
        "censor_date" = "censor_date",
        "weight" = glue("weight_{cohort}_{method}_{spec}"),
        "min_count" = sdc.limit,
        "method" = "constant",
        "contrast" = "TRUE"
      ),
      
      needs = list(
        glue("combine_weights_{cohort}")
      ),
      highly_sensitive = lst(
        arrow = fs::path(dir_output,"*.arrow"),
      ),
      moderately_sensitive = lst(
        csv = fs::path(dir_output,"*.csv"),
        png = fs::path(dir_output,"*.png"),
      )
    )
}



## get km or ci actions function ----
action_plr_contrast <- function(
    cohort, method, spec, subgroup, outcome
){
  dir_output <- glue("output/4-contrast/{cohort}/{method}-{spec}/{subgroup}/{outcome}/plr/")
  
  action(
    name = glue("plr_{cohort}_{method}_{spec}_{subgroup}_{outcome}"),
    run = "r:v2 analysis/4-contrast/plr.R", 
    arguments = c(
      "cohort" = cohort, 
      "method" = method, 
      "spec" = spec,
      "subgroup" = subgroup,
      "outcome" = outcome
    ),
    needs = list(
      #glue("data_selection_{cohort}"),
      glue("combine_weights_{cohort}")
    ),
    highly_sensitive = lst(
      arrow = fs::path(dir_output,"*.arrow"),
    ),
    moderately_sensitive = lst(
      csv = fs::path(dir_output,"*.csv"),
      png = fs::path(dir_output,"*.png"),
    )
  )
}


action_eventcounts <- function(cohort, spec) {
  action(
    name = glue("eventcounts_{cohort}_{spec}"),
    run = glue("r:v2 analysis/eventcounts.R"),
    arguments = c(cohort, spec),
    needs = list(
      glue("match_{cohort}_{spec}"),
      glue("data_selection_{cohort}")
    ),
    highly_sensitive = lst(
      rds = glue("output/{cohort}/{spec}/eventcounts/*.rds"),
    )
  )
}


## model action function ----
action_combine <- function(
    cohort
){

  splice(
    action(
      name = glue("combine_{cohort}_descriptives"),
      run = glue("r:v2 analysis/combine_descriptives.R"),
      arguments = c(cohort),
      needs = list(
        glue("data_selection_{cohort}"),
        glue("match_report_{cohort}_A"),
        glue("match_report_{cohort}_B"),
        glue("eventcounts_{cohort}_A"),
        glue("eventcounts_{cohort}_B")
      ),
      moderately_sensitive = lst(
        csv = glue("output/combine/{cohort}/descriptives/*.csv"),
      )
    ),
    action(
      name = glue("combine_{cohort}_contrasts"),
      run = glue("r:v2 analysis/combine_contrasts.R"),
      arguments = c(cohort),
      needs = glue_data(
        .x=expand_grid(
          spec = c("A", "B"),
          subgroup=c("all", "ageband", "cv", "vax_previous_group"),
          outcome=c("covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath", "fracture", "pericarditis", "myocarditis")
        ) %>% filter(cohort!=subgroup),
        "km_{cohort}_{spec}_{subgroup}_{outcome}"
      ),
      moderately_sensitive = lst(
        csv = glue("output/combine/{cohort}/contrasts/*.csv"),
        png = glue("output/combine/{cohort}/contrasts/plots/*.png"),
      )
    )
  )
}

# specify project ----

## defaults ----
defaults_list <- lst(
  version = "3.0",
  expectations= lst(population_size=1000L)
)

## actions ----
actions_list <- splice(

  comment("# # # # # # # # # # # # # # # # # # #",
          "DO NOT EDIT project.yaml DIRECTLY",
          "This file is created by create-project.R",
          "Edit and run create-project.R to update the project.yaml",
          "# # # # # # # # # # # # # # # # # # #"
          ),


  comment("# # # # # # # # # # # # # # # # # # #", "Pre-server scripts", "# # # # # # # # # # # # # # # # # # #"),

  # doesn't currently work as "project.yaml is not an allowed file type"
  # action(
  #   name = "checkyaml",
  #   run = "r:v2 create-project.R",
  #   moderately_sensitive = lst(
  #     project = "project.yaml"
  #   )
  # ),

  comment("# # # # # # # # # # # # # # # # # # #", "Extract and tidy", "# # # # # # # # # # # # # # # # # # #"),

  action(
    name = "extract",
    run = paste(
      "ehrql:v1 generate-dataset analysis/1-extract/dataset_definition.py", 
      "--output output/1-extract/extract.arrow",
      "--dummy-data-file analysis/1-extract/dummy_extract.arrow"
    ),
    needs = list(),
    highly_sensitive = lst(
      arrow = "output/1-extract/extract.arrow"
    )
  ),


  action(
    name = "prepare",
    run = "r:v2 analysis/2-prepare/data_prepare.R",
    needs = list("extract"),
    highly_sensitive = lst(
      arrow = "output/2-prepare/*.arrow",
    )
  ),


  comment("# # # # # # # # # # # # # # # # # # #", "Cohort: age75plus", "# # # # # # # # # # # # # # # # # # #"),

  action_selection("age75plus"),

  comment("# # # # # # # # # # # # # # # # # # #", "matching set A", "# # # # # # # # # # # # # # # # # # #"),

  action_match("age75plus", "A"),
  # action_match("age75plus", "B"),....
  
  comment("# # # # # # # # # # # # # # # # # # #", "weighting set A", "# # # # # # # # # # # # # # # # # # #"),
  
  action_weight("age75plus", "A"),
  # action_weight("age75plus", "B"),....
  
  comment("# # # # # # # # # # # # # # # # # # #", "combine weights from all adjustment strategies", "# # # # # # # # # # # # # # # # # # #"),
  
  action_combine_weights("age75plus"),
  
  comment("# # # # # # # # # # # # # # # # # # #", "Estimate cumulative incidence curves", "# # # # # # # # # # # # # # # # # # #"),
  
#   comment("### Overall models ('all')"),
#
#   action_km_contrasts("age75plus", "A", "all", "covidemergency"),
  action_km_contrast("age75plus", "weight", "A", "all", "covid_admitted"),
  action_plr_contrast("age75plus", "weight", "A", "all", "covid_admitted"),
#   action_km_contrasts("age75plus", "A", "all", "covidcritcare"),
#   action_km_contrasts("age75plus", "A", "all", "coviddeath"),
#   action_km_contrasts("age75plus", "A", "all", "noncoviddeath"),
#   action_km_contrasts("age75plus", "A", "all", "fracture"),
#   action_km_contrasts("age75plus", "A", "all", "pericarditis"),
#   action_km_contrasts("age75plus", "A", "all", "myocarditis"),
#
#   comment("### Models by age-band ('ageband')"),
#
#   action_contrasts("age75plus", "A", "ageband", "covidemergency"),
#   action_contrasts("age75plus", "A", "ageband", "covidadmitted"),
#   action_contrasts("age75plus", "A", "ageband", "covidcritcare"),
#   action_contrasts("age75plus", "A", "ageband", "coviddeath"),
#   action_contrasts("age75plus", "A", "ageband", "noncoviddeath"),
#   action_contrasts("age75plus", "A", "ageband", "fracture"),
#   action_contrasts("age75plus", "A", "ageband", "pericarditis"),
#   action_contrasts("age75plus", "A", "ageband", "myocarditis"),
#
#   comment("### Models by clinically vulnerable group ('cv')"),
#
#   action_contrasts("age75plus", "A", "cv", "covidemergency"),
#   action_contrasts("age75plus", "A", "cv", "covidadmitted"),
#   action_contrasts("age75plus", "A", "cv", "covidcritcare"),
#   action_contrasts("age75plus", "A", "cv", "coviddeath"),
#   action_contrasts("age75plus", "A", "cv", "noncoviddeath"),
#   action_contrasts("age75plus", "A", "cv", "fracture"),
#   action_contrasts("age75plus", "A", "cv", "pericarditis"),
#   action_contrasts("age75plus", "A", "cv", "myocarditis"),
#
#   comment("### Models by vax history ('vax_previous_group')"),
#
#   action_contrasts("age75plus", "A", "vax_previous_group", "covidemergency"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "covidadmitted"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "covidcritcare"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "coviddeath"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "noncoviddeath"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "fracture"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "pericarditis"),
#   action_contrasts("age75plus", "A", "vax_previous_group", "myocarditis"),
#
#   # comment("### Models by prior infection ('prior_covid_infection')"),
#   #
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "covidemergency"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "covidadmitted"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "covidcritcare"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "coviddeath"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "noncoviddeath"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "fracture"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "pericarditis"),
#   # action_contrasts("age75plus", "A", "prior_covid_infection", "myocarditis"),
#
#   action_eventcounts("age75plus", "A"),
#
#   comment("# # # # # # # # # # # # # # # # # # #", "matching set B", "# # # # # # # # # # # # # # # # # # #"),
#
#   action_match("age75plus", "B"),
#
#   comment("### Overall models ('all')"),
#
#   action_contrasts("age75plus", "B", "all", "covidemergency"),
#   action_contrasts("age75plus", "B", "all", "covidadmitted"),
#   action_contrasts("age75plus", "B", "all", "covidcritcare"),
#   action_contrasts("age75plus", "B", "all", "coviddeath"),
#   action_contrasts("age75plus", "B", "all", "noncoviddeath"),
#   action_contrasts("age75plus", "B", "all", "fracture"),
#   action_contrasts("age75plus", "B", "all", "pericarditis"),
#   action_contrasts("age75plus", "B", "all", "myocarditis"),
#
#   comment("### Models by age-band ('ageband')"),
#
#   action_contrasts("age75plus", "B", "ageband", "covidemergency"),
#   action_contrasts("age75plus", "B", "ageband", "covidadmitted"),
#   action_contrasts("age75plus", "B", "ageband", "covidcritcare"),
#   action_contrasts("age75plus", "B", "ageband", "coviddeath"),
#   action_contrasts("age75plus", "B", "ageband", "noncoviddeath"),
#   action_contrasts("age75plus", "B", "ageband", "fracture"),
#   action_contrasts("age75plus", "B", "ageband", "pericarditis"),
#   action_contrasts("age75plus", "B", "ageband", "myocarditis"),
#
#   comment("### Models by clinically vulnerable group ('cv')"),
#
#   action_contrasts("age75plus", "B", "cv", "covidemergency"),
#   action_contrasts("age75plus", "B", "cv", "covidadmitted"),
#   action_contrasts("age75plus", "B", "cv", "covidcritcare"),
#   action_contrasts("age75plus", "B", "cv", "coviddeath"),
#   action_contrasts("age75plus", "B", "cv", "noncoviddeath"),
#   action_contrasts("age75plus", "B", "cv", "fracture"),
#   action_contrasts("age75plus", "B", "cv", "pericarditis"),
#   action_contrasts("age75plus", "B", "cv", "myocarditis"),
#
#   comment("### Models by vax history ('vax_previous_group')"),
#
#   action_contrasts("age75plus", "B", "vax_previous_group", "covidemergency"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "covidadmitted"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "covidcritcare"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "coviddeath"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "noncoviddeath"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "fracture"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "pericarditis"),
#   action_contrasts("age75plus", "B", "vax_previous_group", "myocarditis"),
#
#   # comment("### Models by prior infection ('prior_covid_infection')"),
#   #
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "covidemergency"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "covidadmitted"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "covidcritcare"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "coviddeath"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "noncoviddeath"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "fracture"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "pericarditis"),
#   # action_contrasts("age75plus", "B", "prior_covid_infection", "myocarditis"),
#
#
#   action_eventcounts("age75plus", "B"),
#


#   comment("# # # # # # # # # # # # # # # # # # #", "Cohort: cv", "# # # # # # # # # # # # # # # # # # #"),
#
#   action_selection("cv"),
#
#
#   comment("# # # # # # # # # # # # # # # # # # #", "matching set A", "# # # # # # # # # # # # # # # # # # #"),
#
#   action_match("cv", "A"),
#
#   comment("### Overall models ('all')"),
#
#   action_contrasts("cv", "A", "all", "covidemergency"),
#   action_contrasts("cv", "A", "all", "covidadmitted"),
#   action_contrasts("cv", "A", "all", "covidcritcare"),
#   action_contrasts("cv", "A", "all", "coviddeath"),
#   action_contrasts("cv", "A", "all", "noncoviddeath"),
#   action_contrasts("cv", "A", "all", "fracture"),
#   action_contrasts("cv", "A", "all", "pericarditis"),
#   action_contrasts("cv", "A", "all", "myocarditis"),
#
#   comment("### Models by age-band ('ageband')"),
#
#   action_contrasts("cv", "A", "ageband", "covidemergency"),
#   action_contrasts("cv", "A", "ageband", "covidadmitted"),
#   action_contrasts("cv", "A", "ageband", "covidcritcare"),
#   action_contrasts("cv", "A", "ageband", "coviddeath"),
#   action_contrasts("cv", "A", "ageband", "noncoviddeath"),
#   action_contrasts("cv", "A", "ageband", "fracture"),
#   action_contrasts("cv", "A", "ageband", "pericarditis"),
#   action_contrasts("cv", "A", "ageband", "myocarditis"),
#
#   comment("### Models by vax history ('vax_previous_group')"),
#
#   action_contrasts("cv", "A", "vax_previous_group", "covidemergency"),
#   action_contrasts("cv", "A", "vax_previous_group", "covidadmitted"),
#   action_contrasts("cv", "A", "vax_previous_group", "covidcritcare"),
#   action_contrasts("cv", "A", "vax_previous_group", "coviddeath"),
#   action_contrasts("cv", "A", "vax_previous_group", "noncoviddeath"),
#   action_contrasts("cv", "A", "vax_previous_group", "fracture"),
#   action_contrasts("cv", "A", "vax_previous_group", "pericarditis"),
#   action_contrasts("cv", "A", "vax_previous_group", "myocarditis"),
#
#   # comment("### Models by prior infection ('prior_covid_infection')"),
#   #
#   # action_contrasts("cv", "A", "prior_covid_infection", "covidemergency"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "covidadmitted"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "covidcritcare"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "coviddeath"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "noncoviddeath"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "fracture"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "pericarditis"),
#   # action_contrasts("cv", "A", "prior_covid_infection", "myocarditis"),
#
#
#   action_eventcounts("cv", "A"),
#
#   comment("# # # # # # # # # # # # # # # # # # #", "matching set B", "# # # # # # # # # # # # # # # # # # #"),
#
#   action_match("cv", "B"),
#
#   comment("### Overall models ('all')"),
#
#   action_contrasts("cv", "B", "all", "covidemergency"),
#   action_contrasts("cv", "B", "all", "covidadmitted"),
#   action_contrasts("cv", "B", "all", "covidcritcare"),
#   action_contrasts("cv", "B", "all", "coviddeath"),
#   action_contrasts("cv", "B", "all", "noncoviddeath"),
#   action_contrasts("cv", "B", "all", "fracture"),
#   action_contrasts("cv", "B", "all", "pericarditis"),
#   action_contrasts("cv", "B", "all", "myocarditis"),
#
#   comment("### Models by age-band ('ageband')"),
#
#   action_contrasts("cv", "B", "ageband", "covidemergency"),
#   action_contrasts("cv", "B", "ageband", "covidadmitted"),
#   action_contrasts("cv", "B", "ageband", "covidcritcare"),
#   action_contrasts("cv", "B", "ageband", "coviddeath"),
#   action_contrasts("cv", "B", "ageband", "noncoviddeath"),
#   action_contrasts("cv", "B", "ageband", "fracture"),
#   action_contrasts("cv", "B", "ageband", "pericarditis"),
#   action_contrasts("cv", "B", "ageband", "myocarditis"),
#
#
#   comment("### Models by vax history ('vax_previous_group')"),
#
#   action_contrasts("cv", "B", "vax_previous_group", "covidemergency"),
#   action_contrasts("cv", "B", "vax_previous_group", "covidadmitted"),
#   action_contrasts("cv", "B", "vax_previous_group", "covidcritcare"),
#   action_contrasts("cv", "B", "vax_previous_group", "coviddeath"),
#   action_contrasts("cv", "B", "vax_previous_group", "noncoviddeath"),
#   action_contrasts("cv", "B", "vax_previous_group", "fracture"),
#   action_contrasts("cv", "B", "vax_previous_group", "pericarditis"),
#   action_contrasts("cv", "B", "vax_previous_group", "myocarditis"),
#
#   # comment("### Models by prior infection ('prior_covid_infection')"),
#   #
#   # action_contrasts("cv", "B", "prior_covid_infection", "covidemergency"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "covidadmitted"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "covidcritcare"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "coviddeath"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "noncoviddeath"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "fracture"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "pericarditis"),
#   # action_contrasts("cv", "B", "prior_covid_infection", "myocarditis"),
#
#   action_eventcounts("cv", "B"),
#
#
#   comment("# # # # # # # # # # # # # # # # # # #", "Combine estimates across cohorts, specs, outcomes and subgroups", "# # # # # # # # # # # # # # # # # # #"),
#
#   # action_contrasts_combine(
#   #   "age75plus",
#   #   "A",
#   #   subgroups = c("all", "ageband", "cv", "vax_previous_group"),
#   #   outcomes = c("covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath", "fracture", "pericarditis", "myocarditis")
#   # ),
#   #
#   # action_contrasts_combine(
#   #   "age75plus",
#   #   "B",
#   #   subgroups = c("all", "ageband", "cv", "vax_previous_group"),
#   #   outcomes = c("covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath", "fracture", "pericarditis", "myocarditis")
#   # ),
#   #
#   # action_contrasts_combine(
#   #   "cv",
#   #   "A",
#   #   subgroups = c("all", "ageband", "vax_previous_group"),
#   #   outcomes = c("covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath", "fracture", "pericarditis", "myocarditis")
#   # ),
#   #
#   # action_contrasts_combine(
#   #   "cv",
#   #   "B",
#   #   subgroups = c("all", "ageband", "vax_previous_group"),
#   #   outcomes = c("covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath", "fracture", "pericarditis", "myocarditis")
#   # ),
#
#   action_combine("age75plus"),
#   action_combine("cv"),
#
#   comment("# # # # # # # # # # # # # # # # # # #", "Files for release", "# # # # # # # # # # # # # # # # # # #"),
#
#   action(
#     name = "release_objects",
#     run = "r:v2 analysis/release_objects.R",
#     needs = list(
#       "combine_age75plus_descriptives",
#       "combine_age75plus_contrasts",
#       "combine_cv_descriptives",
#       "combine_cv_contrasts"
#     ),
#     moderately_sensitive = lst(
#       releaselist = "output/files-for-release.txt",
#       command = "output/osrelease-command.txt",
#       output1 = "output/release-objects/*.csv",
#       output2 = "output/release-objects/*/*.csv",
#     )
#   ),

comment("# # # # # # # # # # # # # # # # # # #", "End", "# # # # # # # # # # # # # # # # # # #")

)


project_list <- splice(
  defaults_list,
  list(actions = actions_list)
)

## convert list to yaml, reformat comments and whitespace ----
thisproject <-
  as.yaml(project_list, indent=2) %>%
  # convert comment actions to comments
  convert_comment_actions() %>%
  # add one blank line before level 1 and level 2 keys
  str_replace_all("\\\n(\\w)", "\n\n\\1") %>%
  str_replace_all("\\\n\\s\\s(\\w)", "\n\n  \\1")


# if running via opensafely, check that the project on disk is the same as the project created here:
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("expectations", "tpp")){

  thisprojectsplit <- str_split(thisproject, "\n")[[1]]
  currentproject <- readLines(here("project.yaml"))

  stopifnot("project.yaml is not up-to-date with create-project.R.  Run create-project.R before running further actions." = identical(thisprojectsplit, currentproject))

# if running manually, output new project as normal
} else if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("")){

## output to file ----
  writeLines(thisproject, here("project.yaml"))
#yaml::write_yaml(project_list, file =here("project.yaml"))

## grab all action names and send to a txt file

names(actions_list) %>% tibble(action=.) %>%
  mutate(
    model = action==""  & lag(action!="", 1, TRUE),
    model_number = cumsum(model),
  ) %>%
  group_by(model_number) %>%
  summarise(
    sets = str_trim(paste(action, collapse=" "))
  ) %>% pull(sets) %>%
  paste(collapse="\n") %>%
  writeLines(here("actions.txt"))

# fail if backend not recognised
} else {
  stop("Backend not recognised")
}



