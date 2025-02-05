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
  weightset <- "A"
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  weightset <- args[[2]]
}



# create output directories ----

output_dir <- here_glue("output", "3-cohorts", cohort, "weight{weightset}", "report")
fs::dir_create(output_dir)


## import weighting info ----
data_weights <- read_feather(here_glue("output", "3-cohorts", cohort, "weight{weightset}", "data_weights.arrow"))


# table 1 style baseline characteristics amongst those eligible for matching ----

# append all characteristics to match data


# TODO: use cobalt package here to produce balance table for pseudo population
