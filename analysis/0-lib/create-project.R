######################################
# Purpose:
# create the project-yaml file
######################################


# Preliminaries ----

## Import libraries ----
library("tidyverse")
library("yaml")
library("here")
library("glue")


## Import custom user functions from lib
source(here("analysis", "0-lib", "utility.R"))

## Import design elements
source(here("analysis", "0-lib", "design.R"))

metaparams <-
  metaparams |>
  filter(cohort == "age75plus", method != "lmw")

## restrict metaparams to those current available:
metaparams_cohort_method_spec <-
  metaparams |>
  select(cohort, method, spec) |>
  unique()

metaparams_cohort_method_spec_subgroup_outcome <-
  metaparams |>
  select(cohort, method, spec, subgroup, outcome) |>
  unique()


# create action functions ----

## create comment function ----
comment <- function(...) {
  list_comments <- list(...)
  comments <- map(list_comments, ~ paste0("## ", ., " ##"))
  comments
}

splice2 <- function(..., name_spec = "{inner}", name_repair = "minimal") {
  list_flatten(x = lst(...), name_spec = name_spec, name_repair = name_repair)
}

## create function to convert comment "actions" in a yaml string into proper comments
convert_comment_actions <- function(yaml.txt) {
  yaml.txt %>%
    str_replace_all("\\\n(\\s*)\\'\\'\\:(\\s*)\\'", "\n\\1")  %>%
    str_replace_all("([^\\'])\\\n(\\s*)\\#\\#", "\\1\n\n\\2\\#\\#") %>%
    str_replace_all("\\#\\#\\'\\\n", "\n")
}


## generic action function ----
action <- function(
    name,
    run,
    arguments = NULL,
    needs = NULL,
    highly_sensitive = NULL,
    moderately_sensitive = NULL,
    ... # other arguments / options for special action types
    ) {

  outputs <- list(
    highly_sensitive = highly_sensitive,
    moderately_sensitive = moderately_sensitive
  )
  outputs[sapply(outputs, is.null)] <- NULL


  if (!is.null(names(arguments))) {
    arguments <- paste0("--", names(arguments), "=", arguments)
  }

  action <- list(
    run = paste(c(run, arguments), collapse = " "),
    needs = needs,
    outputs = outputs,
    ... = ...
  )
  action[sapply(action, is.null)] <- NULL

  action_list <- list(name = action)
  names(action_list) <- name

  action_list
}


## select action function ----

action_select <- function(cohort) {
  action(
    name = glue("select_{cohort}"),
    run = "r:v2 analysis/2-select/select.R",
    arguments = c(cohort),
    needs = list("prepare"),
    highly_sensitive = lst(
      arrow = glue("output/2-select/{cohort}/*.arrow"),
    ),
    moderately_sensitive = lst(
      csv = glue("output/2-select/{cohort}/*.csv"),
    )
  )
}

## balance actions function ----
action_balance <- function(cohort, method, spec) {

  splice2(

    action(
      name = glue("adjust_{cohort}_{method}_{spec}"),
      run = "r:v2 analysis/3-adjust/balance.R",
      arguments = c(cohort, method, spec),
      needs = list(glue("select_{cohort}")),
      highly_sensitive = lst(
        arrow = glue("output/3-adjust/{cohort}/{method}-{spec}/*.arrow")
      )
    ),

    action(
      name = glue("adjust_{cohort}_{method}_{spec}_report"),
      run = "r:v2 analysis/3-adjust/report.R",
      arguments = c(cohort, method, spec),
      needs = list(
        glue("select_{cohort}"),
        glue("adjust_{cohort}_{method}_{spec}")
      ),
      # highly_sensitive = lst(
      #   arrow = glue("output/3-adjust/{cohort}/{method}-{spec}/report/*.arrow"),
      # ),
      moderately_sensitive = lst(
        csv = glue("output/3-adjust/{cohort}/{method}-{spec}/report/*.csv"),
        png = glue("output/3-adjust/{cohort}/{method}-{spec}/report/*.png")
      )
    )
  )


}


action_combine_weights <- function(cohort) {
  cohort0 <- cohort
  action(
    name = glue("adjust_combine_{cohort}"),
    run = glue("r:v2 analysis/3-adjust/combine-weights.R {cohort}"),
    needs = list(
      glue("select_{cohort}"),
      glue_data(
        .x = metaparams_cohort_method_spec |>
          select(cohort, method, spec) |>
          filter(cohort == cohort0),
        "adjust_{cohort}_{method}_{spec}"
      )
    ) |> list_c(),

    highly_sensitive = lst(
      arrow = glue("output/3-adjust/{cohort}/combine/*.arrow"),
    )
  )
}

## get AJ actions function ----
action_aj_contrast <- function(
    cohort, method, spec, subgroup, outcome
    ) {
  dir_output <- glue("output/4-contrast/{cohort}/{method}-{spec}/{subgroup}/{outcome}/aj/")

  action(
    name = glue("aj_{cohort}_{method}_{spec}_{subgroup}_{outcome}"),
    run = "r:v2 analysis/4-contrast/aj.R",
    # arguments = c(cohort, method, spec, subgroup, outcome),
    arguments = c(
      "df_input" = glue("output/3-adjust/{cohort}/combine/data_weights.arrow"),
      "dir_output" = dir_output,
      "exposure" = "treatment",
      "subgroups" = glue("{subgroup}"),
      "origin_date" = "vax_date",
      "event_date" = glue("{outcome}_date"),
      "censoring_date" = "censor_date",
      "competing_date" = "death_date",
      "weight" = glue("wt_{cohort}_{method}_{spec}"),
      "max_fup" = maxfup,
      "min_count" = sdc.limit,
      "method" = "constant",
      "contrast" = "TRUE",
      "plot" = "TRUE"
    ),

    needs = list(
      glue("adjust_combine_{cohort}")
    ),
    highly_sensitive = lst(
      arrow = fs::path(dir_output, "*.arrow"),
    ),
    moderately_sensitive = lst(
      csv = fs::path(dir_output, "*.csv"),
      png = fs::path(dir_output, "*.png"),
    )
  )
}



## get PLR actions function ----
action_plr_contrast <- function(
    cohort, method, spec, subgroup, outcome
    ) {
  dir_output <- glue("output/4-contrast/{cohort}/{method}-{spec}/{subgroup}/{outcome}/plr/")

  action(
    name = glue("plr_{cohort}_{method}_{spec}_{subgroup}_{outcome}"),
    run = "r:v2 analysis/4-contrast/plr.R",
    arguments = c(
      "cohort" = glue("{cohort}"),
      "method" = glue("{method}"),
      "spec" = glue("{spec}"),
      "subgroup" = glue("{subgroup}"),
      "outcome" = glue("{outcome}")
    ),
    needs = list(
      # glue("select_{cohort}"),
      glue("adjust_combine_{cohort}")
    ),
    highly_sensitive = lst(
      arrow = fs::path(dir_output, "*.arrow"),
    ),
    moderately_sensitive = lst(
      csv = fs::path(dir_output, "*.csv"),
      png = fs::path(dir_output, "*.png"),
    )
  )
}

## model action function ----
action_contrasts_combine <- function(
    cohort
    ) {

  action(
    name = glue("combine_{cohort}_contrasts"),
    run = glue("r:v2 analysis/4-contrast/combine-contrasts.R"),
    arguments = c(cohort),
    needs = glue_data(
      .x = metaparams |>
        select(spec, method, subgroup, outcome) |>
        unique() |>
        expand_grid(strategy = c("plr", "aj")),
      "{strategy}_{cohort}_{method}_{spec}_{subgroup}_{outcome}"
    ),
    moderately_sensitive = lst(
      csv = glue("output/4-contrast/{cohort}/contrasts/*.csv"),
      png = glue("output/4-contrast/{cohort}/contrasts/plots/*.png"),
    )
  )
}

# specify project ----

## defaults ----
defaults_list <- lst(
  version = "3.0",
  expectations = lst(population_size = 1000L)
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
      "--dummy-data-file analysis/1-extract/dummy_extract.arrow",
      "--",
      glue("--studystart_date={studystart_date}"),
      glue("--studyend_date={studyend_date}")
    ),
    needs = list(),
    highly_sensitive = lst(
      arrow = "output/1-extract/extract.arrow"
    )
  ),


  action(
    name = "prepare",
    run = "r:v2 analysis/2-select/prepare.R",
    needs = list("extract"),
    highly_sensitive = lst(
      arrow = "output/2-select/*.arrow",
    )
  ),


  comment("# # # # # # # # # # # # # # # # # # #", "Cohort: age75plus", "# # # # # # # # # # # # # # # # # # #"),

  action_select("age75plus"),

  comment("# # # # # # # # # # # # # # # # # # #", "Matching", "# # # # # # # # # # # # # # # # # # #"),

  action_balance("age75plus", "match", "A"),
  action_balance("age75plus", "match", "B"),

  comment("# # # # # # # # # # # # # # # # # # #", "Weighting", "# # # # # # # # # # # # # # # # # # #"),

  action_balance("age75plus", "weight", "A"),
  action_balance("age75plus", "weight", "B"),

  # comment("# # # # # # # # # # # # # # # # # # #", "LMW", "# # # # # # # # # # # # # # # # # # #"),
  #
  # action_balance("age75plus", "lmw", "A"),
  # action_balance("age75plus", "lmw", "B"),


  comment("# # # # # # # # # # # # # # # # # # #", "combine weights from all adjustment strategies", "# # # # # # # # # # # # # # # # # # #"),

  action_combine_weights("age75plus"),

  comment("# # # # # # # # # # # # # # # # # # #", "Estimate cumulative incidence curves", "# # # # # # # # # # # # # # # # # # #"),

  comment("### Aalen-Johansen estimates"),
  pmap(
    metaparams,
    function(cohort, method, spec, subgroup, outcome, ...) {
      action_aj_contrast(cohort, method, spec, subgroup, outcome)
    }
  ) |> list_flatten(),

  comment("### Pooled logistic regression estimates"),
  pmap(
    metaparams,
    function(cohort, method, spec, subgroup, outcome, ...) {
      action_plr_contrast(cohort, method, spec, subgroup, outcome)
    }
  ) |> list_flatten(),

  comment("# # # # # # # # # # # # # # # # # # #", "Combine estimates across specs, outcomes and subgroups", "# # # # # # # # # # # # # # # # # # #"),

  action_contrasts_combine(
    "age75plus"
  ),

  comment("# # # # # # # # # # # # # # # # # # #", "End", "# # # # # # # # # # # # # # # # # # #")

)


project_list <- splice2(
  defaults_list,
  list(actions = actions_list)
)

## convert list to yaml, reformat comments and whitespace ----
thisproject <-
  as.yaml(project_list, indent = 2) %>%
  # convert comment actions to comments
  convert_comment_actions() %>%
  # add one blank line before level 1 and level 2 keys
  str_replace_all("\\\n(\\w)", "\n\n\\1") %>%
  str_replace_all("\\\n\\s\\s(\\w)", "\n\n  \\1")


# if running via opensafely, check that the project on disk is the same as the project created here:
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("expectations", "tpp")) {

  thisprojectsplit <- str_split(thisproject, "\n")[[1]]
  currentproject <- readLines(here("project.yaml"))

  stopifnot("project.yaml is not up-to-date with create-project.R.  Run create-project.R before running further actions." = identical(thisprojectsplit, currentproject))

  # if running manually, output new project as normal
} else if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("")) {
  ## output to file ----
  writeLines(thisproject, here("project.yaml"))
  # yaml::write_yaml(project_list, file =here("project.yaml"))

  ## grab all action names and send to a txt file

  names(actions_list) %>% tibble(action = .) %>%
    mutate(
      model = action == ""  & lag(action != "", 1, TRUE),
      model_number = cumsum(model),
    ) %>%
    group_by(model_number) %>%
    summarise(
      sets = str_trim(paste(action, collapse = " "))
    ) %>% pull(sets) %>%
    paste(collapse = "\n") %>%
    writeLines(here("actions.txt"))

  # fail if backend not recognised
} else {
  stop("Backend not recognised")
}
