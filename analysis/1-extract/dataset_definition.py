
from pathlib import Path

from ehrql import Dataset, case, days, years, when, minimum_of
from ehrql.tables.tpp import (
  patients, 
  practice_registrations, 
  addresses,
  medications, 
  clinical_events,
  vaccinations, 
  occupation_on_covid_vaccine_record, 
  ons_deaths,
  sgss_covid_all_tests as covid_tests,
  emergency_care_attendances as ecds,
  apcs,
)

import codelists

import variables



#######################################################################################
# Import action parameters
#######################################################################################

import argparse
import sys
import logging

logging.info(f"command line args are: {sys.argv}")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--studystart_date",
    type=str,
)
parser.add_argument(
    "--studyend_date",
    type=str,
)
args = parser.parse_args()

studystart_date = args.studystart_date
studyend_date = args.studyend_date

#######################################################################################
# Initialise dataset
#######################################################################################

dataset = Dataset()


#######################################################################################
# Vaccine and eligibility info
#######################################################################################

covid_vaccinations = (
  vaccinations
  .where(vaccinations.target_disease.is_in(["SARS-2 CORONAVIRUS"]))
  .sort_by(vaccinations.date)
)

# first vaccination occurring on or between study start date and study end date, is taken to be the vaccination event of interest
index_vaccination = (
  covid_vaccinations
  .where(covid_vaccinations.date.is_on_or_between(studystart_date, studyend_date))
  .first_for_patient()
)

# Index vaccination date
vax_date = index_vaccination.date

# Main vaccination variables: date and product type
dataset.vax_date = vax_date
dataset.vax_product = index_vaccination.product_name

# We define baseline variables on the day _before_ the study date,
# because we assume any health interactions or events recorded on the day of vaccination
# occurred after vaccination
# TODO: consider whether this is necessary and whether to include vax_date or not when defining baseline variables
baseline_date = vax_date - days(1)

# address info as at vaccination date
address = addresses.for_patient_on(vax_date)

# registration info as at vaccination date
registration = practice_registrations.for_patient_on(vax_date)


#######################################################################################
# Admin and demographics
#######################################################################################

# administrative region of practice
dataset.region = registration.practice_nuts1_region_name

# STP
dataset.stp = registration.practice_stp

# practice deregistration date
dataset.dereg_date = registration.end_date

#Pseudo practice ID
#dataset.practice_id = registration.practice_pseudo_id
# patient has continuous practice registration at least 6 weeks prior to vaccination date
dataset.registered_previous_6weeks = (
    registration
    .start_date.is_on_or_before(vax_date - days(6 * 7))
)

# Middle Super Output Area
#dataset.msoa = address.msoa_code

# Age
dataset.age = patients.age_on(vax_date)

# Age 90 days before start date, for eligiblity definition 
# "Operational flexibility was permitted to offer the booster to eligible individuals 
# expected to reach the target age during the [present] campaign"
dataset.age_eligible = patients.age_on(studystart_date - days(90))

# Sex
dataset.sex = patients.sex

# Ethnicity 
# note that enthicity is documented using codelists, not as a categorical variable
# if ethnicity was transferred from another practice, the date may not have been captured, and will default to 1900-01-01
# we choose here to look at the last known ethnicity recorded _across the entire record_ 
# rather than as known/recorded on the vaccination date, even though this looks into the future.
ethnicity = (clinical_events
  .where(clinical_events.snomedct_code.is_in(codelists.ethnicity16))
  .sort_by(clinical_events.date)
  .last_for_patient()
)

# ethnicity using 5 groups + unknown
dataset.ethnicity5 = ethnicity.snomedct_code.to_category(codelists.ethnicity5)

# ethnicity using 16 groups + unknown
dataset.ethnicity16 = ethnicity.snomedct_code.to_category(codelists.ethnicity16)

# Rurality
#dataset.rural_urban = address.rural_urban_classification

# Index of Multiple Deprevation Rank (rounded down to nearest 100)
dataset.imd = address.imd_rounded



#######################################################################################
# COVID-19 vaccination history
#######################################################################################

# retrieve first 3 vaccines before and first 3 vaccines after study start date
# note that:
# - vax_covid_1_date should be identical to vax_date
# - vax_covid_1_product should be identical to vax_product

variables.add_n_vaccines(
    dataset = dataset, 
    index_date = vax_date, 
    target_disease = "SARS-2 Coronavirus", 
    name = "vax_covid", 
    direction = "on_or_after",
    number_of_vaccines = 3
)

variables.add_n_vaccines(
    dataset = dataset, 
    index_date = vax_date, 
    target_disease = "SARS-2 Coronavirus", 
    name = "vax_covid_prior", 
    direction = "before",
    number_of_vaccines = 3
)

dataset.vax_covid_prior_count = (
  covid_vaccinations
  .where(vaccinations.date.is_before(vax_date))
  .count_for_patient()
)

# TODO: add variables to see if either of Product A or product B have been received prior to curent vax date
# TODO: add variables to see if _any_ vaccine of interest has previously been received

#######################################################################################
# Occupation / residency
#######################################################################################

# Health or social care worker
dataset.hscworker = occupation_on_covid_vaccine_record.where(occupation_on_covid_vaccine_record.is_healthcare_worker).exists_for_patient()

# TPP care home flag
dataset.care_home_tpp = address.care_home_is_potential_match.when_null_then(False)

# Patients in long-stay nursing and residential care
dataset.care_home_code = variables.has_prior_event(codelists.carehome, vax_date)

#######################################################################################
# Other vulnerabilities / predictors of vaccination or vaccination setting (and therefore product)
#######################################################################################

# Overnight hospital admission at time vaccination
dataset.inhospital = (
    apcs
    .where(apcs.admission_date.is_on_or_before(vax_date))
    .where(apcs.discharge_date.is_on_or_after(vax_date))
    .where(apcs.admission_method.is_in(
            ["11", "12", "13", "21", "2A", "22", "23", "24", "25", "2D", "28", "2B", "81"]
        )
    )
    # Ordinary admissions only
    .where(apcs.patient_classification == "1")
    .exists_for_patient()
)

# TODO: add all relevant variables
# eg housebound, end-of-life, carehome residency, etc



#######################################################################################
# Clinical information as at (day before) vaccination date
#######################################################################################

# PRIMIS

variables.primis_variables(dataset, vax_date, var_name_suffix="")



#######################################################################################
# Pre-baseline variables where event date is of interest
#######################################################################################


def prior_hospital_admission(before = None, diagnoses_contains_any_of = None, where = None):
    return (
        apcs
        .where(apcs.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
        .where(apcs.admission_method.is_in(["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"]))
        .where(apcs.patient_classification == "1")  # Ordinary admissions only
        .where(apcs.admission_date.is_before(vax_date))
        .sort_by(apcs.admission_date)
        .last_for_patient()
        .admission_date
    )
    
dataset.covid_admitted_0_date = prior_hospital_admission(vax_date, codelists.covid_icd10)


# TODO: add all relevant variables
# eg additional prior covid related events, possibibly tests


#######################################################################################
# Post-baseline variables: outcomes, competing outcomes, and censoring
#######################################################################################


### Effectiveness outcomes 


# TODO: awaiting `all_diagnoses.contains_any_of` functionality
# def next_emergency_attendance(on_or_after = None, diagnoses_contains_any_of = None, where = None):
#    return (
#       emergency
#       .where(ecds.arrival_date.is_on_or_after(on_or_after))
#       .where(ecds.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
#       .where(where)
#       .sort_by(ecds.arrival_date)
#       .first_for_patient()
#       .arrival_date
#     )

# Covid-related emergency attendance
#dataset.covidemergency_date = next_emergency_attendance(vax_date, codelists.covid_emergency)

# Any emergency attendance
#dataset.emergency_date = next_emergency_attendance(vax_date)


def next_hospital_admission(on_or_after = None, diagnoses_contains_any_of = None, where = None):
    return (
        apcs
        .where(apcs.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
        .where(apcs.admission_method.is_in(["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"]))
        .where(apcs.patient_classification == "1")  # Ordinary admissions only
        .where(apcs.admission_date.is_on_or_after(vax_date))
        .sort_by(apcs.admission_date)
        .first_for_patient()
        .admission_date
    )
  

# covid-related admission 
dataset.covid_admitted_date = next_hospital_admission(vax_date, codelists.covid_icd10)

# covid-related admission to critical care
dataset.covid_critcare_date = next_hospital_admission(vax_date, codelists.covid_icd10, where = apcs.days_in_critical_care>0)


# all-cause death
dataset.death_date = ons_deaths.date


# covid-related death (stated anywhere on death certificate)
dataset.death_cause_covid = ons_deaths.cause_of_death_is_in(codelists.covid_icd10)


# TODO: add all effectiveness outcomes here


### Safety outcomes

#dataset.pericarditis_emergency_date = next_emergency_attendance(vax_date, codelists.pericarditis_snomedECDS)
#dataset.pericarditis_admitted_date = next_hospital_admission(vax_date, codelists.pericarditis_icd10)
#dataset.pericarditis_death_date = ons_deaths.cause_of_death_is_in(codelists.pericarditis_icd10).date

# TODO: add all safety outcomes here

### Negative control outcomes 

#dataset.fracture_emergency_date = next_emergency_attendance(vax_date, codelists.fractures_snomedECDS)
#dataset.fractureadmitted_date = next_hospital_admission(vax_date, codelists.fractures_icd10)
#dataset.fracturedeath_date = ons_deaths.cause_of_death_is_in(codelists.fractures_icd10).date

# TODO: add all negative control outcomes here

# #######################################################################################
# # Population
# #######################################################################################

# From Green Book chapter 14a 
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1186479/Greenbook-chapter-14a-4September2023.pdf
# 
# The vast majority of people aged over 75 years reached an interval of around six months
# from their last dose between late March and June 2023. Operational flexibility was
# permitted to offer the booster to eligible individuals expected to reach the target age
# during the spring campaign. Boosters were offered around six months from the previous
# dose, but could be given a minimum of three months from the previous dose; this was
# particularly important to facilitate delivery of the programme to residents in care homes
# and the housebound.


# define dataset poppulation
dataset.define_population(
  index_vaccination.exists_for_patient() & 
  (dataset.age_eligible >= 18) & (dataset.age_eligible <= 100) &
  (registration.exists_for_patient()) & 
  ((dataset.death_date >= dataset.vax_date) | dataset.death_date.is_null())
)
