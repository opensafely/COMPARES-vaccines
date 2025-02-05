# _________________________________________________
# Purpose:
# define useful functions used in the codebase
# this script should be sourced (using `source(here("analysis", "0-lib", "utility.R"))`) at the start of each R script
# _________________________________________________



roundmid_any <- function(x, to = 1) {
  # like ceiling_any, but centers on (integer) midpoint of the rounding points
  ceiling(x / to) * to - (floor(to / 2) * (x != 0))
}

ceiling_any <- function(x, to = 1) {
  # round to nearest 100 millionth to avoid floating point errors
  ceiling(plyr::round_any(x / to, 1 / 100000000)) * to
}

# get nth largest value from list
nthmax <- function(x, n = 1) {
  dplyr::nth(sort(x, decreasing = TRUE), n)
}

# get nth smallest value from list
nthmin <- function(x, n = 1) {
  dplyr::nth(sort(x, decreasing = FALSE), n)
}

# overwrite splice function to avoid deprecation warnings
# splice <- function(...) {
#   list_flatten(lst(...), name_spec = "{inner}", name_repair = "check_unique")
# }

# uses dplyr::case_when but converts the output to a factor,
# with factors ordered as they appear in the case_when's  ... argument
fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

# recode factor levels based on a "lookup" named vector
# whose names are the new factor levels, and values are the old factor levels
fct_recoderelevel <- function(x, lookup){
  stopifnot(!is.na(names(lookup)))
  factor(x, levels=lookup, labels=names(lookup))
}



# function to convert ethnicity 16 group into 5 group
ethnicity_16_to_5 <- function(x) {
  x1 <- fct_relabel(x, ~ str_extract(.x, ".*(?= - )")) # pick up everything before " - "
  x2 <- fct_recode(x1, `Chinese or Other Ethnic Groups` = "Other Ethnic Groups")
  return(x2)
}

# relabel_from_lookup <- function(x, from, to, source){
#   left_join(tibble(x=x), source, by = {{from}})[[{{to}}]]
# }

# template for standardising characteristics that are extracted multiple times
# using this in mutate like this: `mutate(!!!standardise_characteristics)`
standardise_characteristics <-
  rlang::quos(

    ## --VARIABLES--
    ## demographics
    ageband = cut(
      age,
      breaks = c(-Inf, 18, 40, 55, 65, 75, Inf),
      labels = c("under 18", "18-39", "40-54", "55-64", "65-74", "75+"),
      right = FALSE
    ),
    region = fct_collapse(
      region,
      `East of England` = "East",
      `London` = "London",
      `Midlands` = c("West Midlands", "East Midlands"),
      `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
      `North West` = "North West",
      `South East` = "South East",
      `South West` = "South West"
    ),
    imd_quintile = cut(
      imd,
      breaks = c(0, 32844 * (1:5) / 5),
      labels = c("1 (most deprived)", "2", "3", "4", "5 (least deprived)"),
      include.lowest = TRUE,
      right = FALSE
    ) #,
    ## PRIMIS clinical variables
   # cv = crd | ast | chd | ckd | cld | cns | learndis | diab | immuno | asplen | obes | sev_ment
 #FIXME add additional vulnerability variables when defined and extracted
  )


## factor levels provided in a sensible order, as this won't happen directly from opensafely ----

factor_levels <-
  lst(
    ethnicity5 = c(
      "White",
      "Mixed",
      "Asian or Asian British",
      "Black or Black British",
      "Chinese or Other Ethnic Groups"
    ),
    ethnicity16 = c(
      "White - British",
      "White - Irish",
      "White - Any other White background",
      "Mixed - White and Black Caribbean",
      "Mixed - White and Black African",
      "Mixed - White and Asian",
      "Mixed - Any other mixed background",
      "Asian or Asian British - Indian",
      "Asian or Asian British - Pakistani",
      "Asian or Asian British - Bangladeshi",
      "Asian or Asian British - Any other Asian background",
      "Black or Black British - Caribbean",
      "Black or Black British - African",
      "Black or Black British - Any other Black background",
      "Other Ethnic Groups - Chinese",
      "Other Ethnic Groups - Any other ethnic group"
    ),
  )



# Import dummy data if running locally, or real data if running on the server
import_extract <- function(custom_file_path, ehrql_file_path) {
  if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
    # ideally in future this will check column existence and types from metadata,
    # rather than from a ehrql-generated dummy data

    data_ehrql_dummy <- read_feather(ehrql_file_path) %>%
      # because date types are not returned consistently by ehrql
      mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
      mutate(patient_id = as.integer(patient_id))

    data_custom_dummy <- read_feather(custom_file_path)

    not_in_ehrql <- names(data_custom_dummy)[!(names(data_custom_dummy) %in% names(data_ehrql_dummy))]
    not_in_custom <- names(data_ehrql_dummy)[!(names(data_ehrql_dummy) %in% names(data_custom_dummy))]


    if (length(not_in_custom) != 0) {
      stop(
        paste(
          "These variables are in ehrql but not in custom: ",
          paste(not_in_custom, collapse = ", ")
        )
      )
    }

    if (length(not_in_ehrql) != 0) {
      stop(
        paste(
          "These variables are in custom but not in ehrql: ",
          paste(not_in_ehrql, collapse = ", ")
        )
      )
    }

    # reorder columns
    data_ehrql_dummy <- data_ehrql_dummy[, names(data_custom_dummy)]

    unmatched_types <- cbind(
      map_chr(data_ehrql_dummy, ~ paste(class(.), collapse = ", ")),
      map_chr(data_custom_dummy, ~ paste(class(.), collapse = ", "))
    )[(map_chr(data_ehrql_dummy, ~ paste(class(.), collapse = ", ")) != map_chr(data_custom_dummy, ~ paste(class(.), collapse = ", "))), ] %>%
      as.data.frame() %>%
      rownames_to_column()


    # if(nrow(unmatched_types)>0) stop(
    #   #unmatched_types
    #   "inconsistent typing in ehrql : dummy dataset\n",
    #   apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
    # )

    data_extract <- data_custom_dummy
  } else {
    data_extract <- read_feather(ehrql_file_path) %>%
      # because date types are not returned consistently by ehrql
      mutate(across(ends_with("_date"), as.Date))
  }
  data_extract
}



## Create table1-style summary of characteristics, with SDC applied ----

table1_summary <- function(.data, group, labels, threshold) {
  
  group_quo <- enquo(group)
  
  ## create a table of baseline characteristics between each treatment group, before matching / weighting
  tab_summary <-
    .data |>
    select(
      !!group_quo,
      all_of(names(labels)),
    ) %>%
    tbl_summary(
      by = !!group_quo,
      label = labels,
      statistic = list(
        N ~ "{N}",
        all_continuous() ~ "{median} ({p25}, {p75});  {mean} ({sd})"
      ),
    )
  
  ## extract structured info from tbl_summary object to apply SDC to the counts 
  raw_stats <- 
    tab_summary$cards$tbl_summary |> 
    filter(!(context %in% c("missing", "attributes", "total_n"))) |>
    select(-fmt_fn, -warning, -error, -gts_column) |>
    pivot_wider(
      id_cols = c("group1", "group1_level", "variable", "variable_level", "context"), 
      names_from = stat_name, 
      values_from = stat
    ) |> 
    left_join(
      tab_summary$cards$tbl_summary |> 
        filter(context == "attributes", stat_name == "label") |>
        select(variable, variable_label=stat),
      by="variable"
    ) |>
    mutate(
      across(
        c(n, N),
        ~{
          map_int(., ~{
            if(is.null(.)) NA else as.integer(.)
          })
        }
      ),
      across(
        c(p, median, p25, p75, mean, sd),
        ~{
          map_dbl(., ~{
            if(is.null(.)) NA else as.numeric(.)
          })
        }
      ),
      across(
        c(group1_level, variable_level, variable_label),
        ~{
          map_chr(., ~{
            if(is.null(.)) NA else as.character(.)
          })
        }
      ),
    )
  
  raw_stats_redacted <- 
    raw_stats |>
    mutate(
      n = roundmid_any(n, threshold),
      N = roundmid_any(N, threshold),
      p = n / N,
      # variable_label = factor(variable_label, levels = map_chr(variable_label[-c(1, 2)], ~ last(as.character(.)))),
      #variable_levels = replace_na(as.character(variable_levels), "")
    )
  
  return(raw_stats_redacted)
}


## Print the number of rows and the size on disk of a R object ----

print_data_size <- function(.data){
  object_name <- deparse(substitute(.data))
  cat(
    glue(object_name, " data size = ", nrow(.data)),
    glue(object_name, " memory usage = ", format(object.size(.data), units="GB", standard="SI", digits=3L)),
    "",
    sep = "\n"
  )
}

## use glue with here for easy filepaths
here_glue <- function(...){
  here::here(glue(..., .sep = "/"))
}
