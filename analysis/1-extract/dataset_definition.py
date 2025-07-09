
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

# these are needed for emergency diagnoses workaround
import operator
from functools import reduce
def any_of(conditions):
    return reduce(operator.or_, conditions)

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
# Pre-baseline variables where occurrence of event within a given period is of interest
#######################################################################################

### A & E: last event before baseline (ECDS)
def prior_emergency_attendance(before, years_prior, diagnoses_contains_any_of = None, where = True):
  
    # use this until "contains_any_of" methods works for ecds table
    
    if diagnoses_contains_any_of:
        conditions = [
            getattr(ecds, column_name).is_in(diagnoses_contains_any_of)
            for column_name in [f"diagnosis_{i:02d}" for i in range(1, 25)]
        ]
        ecds_filtered = ecds.where(any_of(conditions))
    else:
        ecds_filtered = ecds
    
    return (
       ecds_filtered
       .where(ecds.arrival_date.is_on_or_between(before - years(years_prior), before - days(1)))
#       .where(ecds.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
       .where(where)
       .exists_for_patient()
     )

### GP: last GP coded-event before baseline (SNOMED)
def prior_gp_event(before, years_prior, codelist = None, where = True):
    return (
        clinical_events
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(clinical_events.date.is_on_or_between(before - years(years_prior), before - days(1)))
        .where(where)
        .exists_for_patient()
    )

### SUS: last hospital attendance before baseline (ICD-10)
def prior_hospital_admission(before, years_prior, diagnoses_contains_any_of = None, where = None):
    
    if diagnoses_contains_any_of:
        apcs_filtered = apcs.where(apcs.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
    else:
        apcs_filtered = apcs
    
    return (
        apcs_filtered
        .where(apcs.admission_method.is_in(["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"]))
        .where(apcs.patient_classification == "1")  # Ordinary admissions only
        .where(apcs.admission_date.is_on_or_between(before - years(years_prior), before - days(1)))
        .exists_for_patient()
    )


## All-cause outcomes ----------------------------------------------------------

# Any emergency attendance
dataset.prior_emergency = prior_emergency_attendance(vax_date, 1)
# Any admission 
dataset.prior_admitted = prior_hospital_admission(vax_date, 1)

## Effectiveness outcomes ------------------------------------------------------

### Covid-related emergency attendance
dataset.covid_prior_emergency = prior_emergency_attendance(vax_date, 1, codelists.covid_emergency)
# covid-related admission 
dataset.covid_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.covid_icd10)
# covid-related admission to critical care
dataset.covid_critcare = prior_hospital_admission(vax_date, 1, codelists.covid_icd10, where = apcs.days_in_critical_care>0)


## Safety outcomes ------------------------------------------------------------

# Neurological -----------------------------------------------------------------

# GUILLAIN BARRE
dataset.sgb_prior_gp = prior_gp_event(vax_date, 1, codelists.sgb_snomed)
dataset.sgb_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.sgb_icd10)

dataset.sgb_prior = any_of([
    dataset.sgb_prior_gp,
    dataset.sgb_prior_admitted
])

# BELL'S PALSY
dataset.bells_palsy_prior_gp = prior_gp_event(vax_date, 1, codelists.bells_palsy_snomed)
dataset.bells_palsy_prior_emergency = prior_emergency_attendance(vax_date, 1, codelists.bells_palsy_ecds)
dataset.bells_palsy_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.bells_palsy_icd10)

dataset.bells_palsy_prior = any_of([
    dataset.bells_palsy_prior_gp,
    dataset.bells_palsy_prior_emergency,
    dataset.bells_palsy_prior_admitted,
])

# THROMBO ----------------------------------------------------------------------

# THROMBOCITOPENIA
dataset.ttp_prior_gp = prior_gp_event(vax_date, 1, codelists.ttp_snomed)
dataset.ttp_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.ttp_icd10)

dataset.ttp_prior = any_of([
    dataset.ttp_prior_gp,
    dataset.ttp_prior_admitted
])
# ARTERIAL THROMBOTIC

### Acute myocardial infarction (ami)
dataset.ami_prior_gp = prior_gp_event(vax_date, 1, codelists.ami_snomed)
dataset.ami_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.ami_icd10)


dataset.ami_prior = any_of([
    dataset.ami_prior_gp,
    dataset.ami_prior_admitted
])

### Ischaemic stroke
dataset.stroke_isch_prior_gp = prior_gp_event(vax_date, 1, codelists.stroke_isch_snomed)
dataset.stroke_isch_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.stroke_isch_icd10)

dataset.stroke_isch_prior = any_of([
    dataset.stroke_isch_prior_gp,
    dataset.stroke_isch_prior_admitted
])

## Composite arterial thrombotic event (ATE)
dataset.ate_prior_gp = prior_gp_event(vax_date, 1, codelists.ate_snomed)
dataset.ate_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.ate_icd10)

dataset.ate_prior = any_of([
    dataset.ate_prior_gp,
    dataset.ate_prior_admitted
])

# VENOUS THROMBOTIC
## Deep vein thrombosis (DVT) [includes during pregnancy]
dataset.dvt_prior_gp = prior_gp_event(vax_date, 1, codelists.dvt_snomed)
dataset.dvt_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.dvt_icd10)

dataset.dvt_prior = any_of([
    dataset.dvt_prior_gp,
    dataset.dvt_prior_admitted
])

## Intracranial venous thrombosis (ICVT) [includes during pregnancy; contributes to composite VTE only]
dataset.icvt_prior_gp = prior_gp_event(vax_date, 1, codelists.icvt_snomed)
dataset.icvt_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.icvt_icd10)

dataset.icvt_prior = any_of([
    dataset.icvt_prior_gp,
    dataset.icvt_prior_admitted
])

## Pulmonary embolism (PE)
dataset.pe_prior_gp = prior_gp_event(vax_date, 1, codelists.pe_snomed)
dataset.pe_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.pe_icd10)

dataset.pe_prior = any_of([
    dataset.pe_prior_gp,
    dataset.pe_prior_admitted
])


## Composite venous thrombotic event (VTE)
dataset.vte_prior_gp = prior_gp_event(vax_date, 1, codelists.vte_snomed)
dataset.vte_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.vte_icd10)

dataset.vte_prior = any_of([
    dataset.vte_prior_gp,
    dataset.vte_prior_admitted
])

# CARDIO ------------------------------------------------------------------------

# PERICARDITIS
dataset.pericarditis_prior_gp = prior_gp_event(vax_date, 1, codelists.pericarditis_snomed)
dataset.pericarditis_prior_emergency = prior_emergency_attendance(vax_date, 1, codelists.pericarditis_ecds)
dataset.pericarditis_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.pericarditis_icd10)

dataset.pericarditis_prior = any_of([
    dataset.pericarditis_prior_gp,
    dataset.pericarditis_prior_emergency,
    dataset.pericarditis_prior_admitted
])


# MYOCARDITIS
dataset.myocarditis_prior_gp = prior_gp_event(vax_date, 1, codelists.myocarditis_snomed)
dataset.myocarditis_prior_emergency = prior_emergency_attendance(vax_date, 1, codelists.myocarditis_ecds)
dataset.myocarditis_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.myocarditis_icd10)

dataset.myocarditis_prior = any_of([
    dataset.myocarditis_prior_gp,
    dataset.myocarditis_prior_emergency,
    dataset.myocarditis_prior_admitted
])

# OTHER --------------------------------------------------------------------------
# HEAVY MENTRUAL BLEEDING
dataset.menorrhagia_prior_gp = prior_gp_event(vax_date, 1, codelists.menorrhagia_snomed)
dataset.menorrhagia_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.menorrhagia_icd10)

dataset.menorrhagia_prior = any_of([
    dataset.menorrhagia_prior_gp,
    dataset.menorrhagia_prior_admitted
])

# ERYTHEMA MULTIFORME
dataset.ery_multi_prior_gp = prior_gp_event(vax_date, 1, codelists.ery_multi_snomed)
dataset.ery_multi_prior_emergency = prior_emergency_attendance(vax_date, 1, codelists.ery_multi_ecds)
dataset.ery_multi_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.ery_multi_icd10)

dataset.ery_multi_prior = any_of([
    dataset.ery_multi_prior_gp,
    dataset.ery_multi_prior_emergency,
    dataset.ery_multi_prior_admitted
])

# ANAPHYLAXIS
dataset.anaphylaxis_prior_gp = prior_gp_event(vax_date, 1, codelists.anaphylaxis_snomed)
dataset.anaphylaxis_prior_emergency = prior_emergency_attendance(vax_date, 1, codelists.anaphylaxis_ecds)
dataset.anaphylaxis_prior_admitted = prior_hospital_admission(vax_date, 1, codelists.anaphylaxis_icd10)

dataset.anaphylaxis_prior = any_of([
    dataset.anaphylaxis_prior_gp,
    dataset.anaphylaxis_prior_emergency,
    dataset.anaphylaxis_prior_admitted
])

#######################################################################################
# Post-baseline variables: outcomes, competing outcomes, and censoring
#######################################################################################

## Functions

### A & E: first event after baseline (ECDS)
def next_emergency_attendance(on_or_after = None, diagnoses_contains_any_of = None, where = True):
  
    # use this until "contains_any_of" methods works for ecds table
    
    if diagnoses_contains_any_of:
        conditions = [
            getattr(ecds, column_name).is_in(diagnoses_contains_any_of)
            for column_name in [f"diagnosis_{i:02d}" for i in range(1, 25)]
        ]
        ecds_filtered = ecds.where(any_of(conditions))
    else:
        ecds_filtered = ecds
    
    return (
       ecds_filtered
       .where(ecds.arrival_date.is_on_or_after(on_or_after))
#       .where(ecds.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
       .where(where)
       .sort_by(ecds.arrival_date)
       .first_for_patient()
       .arrival_date
     )

### GP: first GP coded-event after baseline (SNOMED)
def next_gp_event(on_or_after = None, codelist = None, where = True):
    post_events = clinical_events.where(clinical_events.date.is_on_or_after(on_or_after))
    return (
        post_events
        .where(post_events.snomedct_code.is_in(codelist))
        .where(where)
        .sort_by(clinical_events.date)
        .first_for_patient()
        .date
    )

### SUS: first hospital attendance after baseline (ICD-10)
def next_hospital_admission(on_or_after = None, diagnoses_contains_any_of = None, where = True):
  
    if diagnoses_contains_any_of:
        apcs_filtered = apcs.where(apcs.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
    else:
        apcs_filtered = apcs
    
    return (
        apcs_filtered
        #.where(apcs.all_diagnoses.contains_any_of(diagnoses_contains_any_of))
        .where(apcs.admission_method.is_in(["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"]))
        .where(apcs.patient_classification == "1")  # Ordinary admissions only
        .where(apcs.admission_date.is_on_or_after(vax_date))
        .sort_by(apcs.admission_date)
        .first_for_patient()
        .admission_date
    )

### Death: death after baseline (SNOMED)

def cause_specific_death_date(codelist):
    return (
        case(
            when(ons_deaths.cause_of_death_is_in(codelist)).then(ons_deaths.date),
            otherwise=None
        )
    )

## All-cause outcomes ----------------------------------------------------------

# Any emergency attendance
dataset.emergency_date = next_emergency_attendance(vax_date)
# Any admission 
dataset.admitted_date = next_hospital_admission(vax_date)
# all-cause death
dataset.death_date = ons_deaths.date



## Effectiveness outcomes ------------------------------------------------------

### Covid-related emergency attendance
dataset.covid_emergency_date = next_emergency_attendance(vax_date, codelists.covid_emergency)
# covid-related admission 
dataset.covid_admitted_date = next_hospital_admission(vax_date, codelists.covid_icd10)
# covid-related admission to critical care
dataset.covid_critcare_date = next_hospital_admission(vax_date, codelists.covid_icd10, where = apcs.days_in_critical_care>0)
# covid-related death
dataset.covid_death_date = cause_specific_death_date(codelists.covid_icd10)

## Safety outcomes ------------------------------------------------------------

# Neurological -----------------------------------------------------------------

# GUILLAIN BARRE
dataset.sgb_gp_date = next_gp_event(vax_date, codelists.sgb_snomed)
dataset.sgb_admitted_date = next_hospital_admission(vax_date, codelists.sgb_icd10)
dataset.sgb_death_date = cause_specific_death_date(codelists.sgb_icd10)

dataset.sgb_date = minimum_of(
    dataset.sgb_gp_date,
    dataset.sgb_admitted_date,
    dataset.sgb_death_date
)

# BELL'S PALSY
dataset.bells_palsy_gp_date = next_gp_event(vax_date, codelists.bells_palsy_snomed)
dataset.bells_palsy_emergency_date = next_emergency_attendance(vax_date, codelists.bells_palsy_ecds)
dataset.bells_palsy_admitted_date = next_hospital_admission(vax_date, codelists.bells_palsy_icd10)
dataset.bells_palsy_death_date = cause_specific_death_date(codelists.bells_palsy_icd10)

dataset.bells_palsy_date = minimum_of(
    dataset.bells_palsy_gp_date,
    dataset.bells_palsy_emergency_date,
    dataset.bells_palsy_admitted_date,
    dataset.bells_palsy_death_date
)

# THROMBO ----------------------------------------------------------------------

# THROMBOCITOPENIA
dataset.ttp_gp_date = next_gp_event(vax_date, codelists.ttp_snomed)
dataset.ttp_admitted_date = next_hospital_admission(vax_date, codelists.ttp_icd10)
dataset.ttp_death_date = cause_specific_death_date(codelists.ttp_icd10)

dataset.ttp_date = minimum_of(
    dataset.ttp_gp_date,
    dataset.ttp_admitted_date,
    dataset.ttp_death_date
)
# ARTERIAL THROMBOTIC

### Acute myocardial infarction (ami)
dataset.ami_gp_date = next_gp_event(vax_date, codelists.ami_snomed)
dataset.ami_admitted_date = next_hospital_admission(vax_date, codelists.ami_icd10)
dataset.ami_death_date = cause_specific_death_date(codelists.ami_icd10)

dataset.ami_date = minimum_of(
    dataset.ami_gp_date,
    dataset.ami_admitted_date,
    dataset.ami_death_date
)

### Ischaemic stroke
dataset.stroke_isch_gp_date = next_gp_event(vax_date, (codelists.stroke_isch_snomed))
dataset.stroke_isch_admitted_date = next_hospital_admission(vax_date, codelists.stroke_isch_icd10)
dataset.stroke_isch_death_date = cause_specific_death_date(codelists.stroke_isch_icd10)

dataset.stroke_isch_date = minimum_of(
    dataset.stroke_isch_gp_date,
    dataset.stroke_isch_admitted_date,
    dataset.stroke_isch_death_date
)

## Composite arterial thrombotic event (ATE)
dataset.ate_gp_date = next_gp_event(vax_date, codelists.ate_snomed)
dataset.ate_admitted_date = next_hospital_admission(vax_date, codelists.ate_icd10)
dataset.ate_death_date = cause_specific_death_date(codelists.ate_icd10)

dataset.ate_date = minimum_of(
    dataset.ate_gp_date,
    dataset.ate_admitted_date,
    dataset.ate_death_date
)

# VENOUS THROMBOTIC
## Deep vein thrombosis (DVT) [includes during pregnancy]
dataset.dvt_gp_date = next_gp_event(vax_date, codelists.dvt_snomed)
dataset.dvt_admitted_date = next_hospital_admission(vax_date, codelists.dvt_icd10)
dataset.dvt_death_date = cause_specific_death_date(codelists.dvt_icd10)

dataset.dvt_date = minimum_of(
    dataset.dvt_gp_date,
    dataset.dvt_admitted_date,
    dataset.dvt_death_date
)

## Intracranial venous thrombosis (ICVT) [includes during pregnancy; contributes to composite VTE only]
dataset.icvt_gp_date = next_gp_event(vax_date, codelists.icvt_snomed)
dataset.icvt_admitted_date = next_hospital_admission(vax_date, codelists.icvt_icd10)
dataset.icvt_death_date = cause_specific_death_date(codelists.icvt_icd10)

dataset.icvt_date = minimum_of(
    dataset.icvt_gp_date,
    dataset.icvt_admitted_date,
    dataset.icvt_death_date
)

## Pulmonary embolism (PE)
dataset.pe_gp_date = next_gp_event(vax_date, codelists.pe_snomed)
dataset.pe_admitted_date = next_hospital_admission(vax_date, codelists.pe_icd10)
dataset.pe_death_date = cause_specific_death_date(codelists.pe_icd10)

dataset.pe_date = minimum_of(
    dataset.pe_gp_date,
    dataset.pe_admitted_date,
    dataset.pe_death_date
)


## Composite venous thrombotic event (VTE)
dataset.vte_gp_date = next_gp_event(vax_date, codelists.vte_snomed)
dataset.vte_admitted_date = next_hospital_admission(vax_date, codelists.vte_icd10)
dataset.vte_death_date = cause_specific_death_date(codelists.vte_icd10)

dataset.vte_date = minimum_of(
    dataset.vte_gp_date,
    dataset.vte_admitted_date,
    dataset.vte_death_date
)

# CARDIO ------------------------------------------------------------------------

# PERICARDITIS
dataset.pericarditis_gp_date = next_gp_event(vax_date, codelists.pericarditis_snomed)
dataset.pericarditis_emergency_date = next_emergency_attendance(vax_date, codelists.pericarditis_ecds)
dataset.pericarditis_admitted_date = next_hospital_admission(vax_date, codelists.pericarditis_icd10)
dataset.pericarditis_death_date = cause_specific_death_date(codelists.pericarditis_icd10)

dataset.pericarditis_date = minimum_of(
    dataset.pericarditis_gp_date,
    dataset.pericarditis_emergency_date,
    dataset.pericarditis_admitted_date,
    dataset.pericarditis_death_date
)


# MYOCARDITIS
dataset.myocarditis_gp_date = next_gp_event(vax_date, codelists.myocarditis_snomed)
dataset.myocarditis_emergency_date = next_emergency_attendance(vax_date, codelists.myocarditis_ecds)
dataset.myocarditis_admitted_date = next_hospital_admission(vax_date, codelists.myocarditis_icd10)
dataset.myocarditis_death_date = cause_specific_death_date(codelists.myocarditis_icd10)

dataset.myocarditis_date = minimum_of(
    dataset.myocarditis_gp_date,
    dataset.myocarditis_emergency_date,
    dataset.myocarditis_admitted_date,
    dataset.myocarditis_death_date
)

# OTHER --------------------------------------------------------------------------
# HEAVY MENTRUAL BLEEDING
dataset.menorrhagia_gp_date = next_gp_event(vax_date, codelists.menorrhagia_snomed)
dataset.menorrhagia_admitted_date = next_hospital_admission(vax_date, codelists.menorrhagia_icd10)
dataset.menorrhagia_death_date = cause_specific_death_date(codelists.menorrhagia_icd10)

dataset.menorrhagia_date = minimum_of(
    dataset.menorrhagia_gp_date,
    dataset.menorrhagia_admitted_date,
    dataset.menorrhagia_death_date
)

# ERYTHEMA MULTIFORME
dataset.ery_multi_gp_date = next_gp_event(vax_date, codelists.ery_multi_snomed)
dataset.ery_multi_emergency_date = next_emergency_attendance(vax_date, codelists.ery_multi_ecds)
dataset.ery_multi_admitted_date = next_hospital_admission(vax_date, codelists.ery_multi_icd10)
dataset.ery_multi_death_date = cause_specific_death_date(codelists.ery_multi_icd10)

dataset.ery_multi_date = minimum_of(
    dataset.ery_multi_gp_date,
    dataset.ery_multi_emergency_date,
    dataset.ery_multi_admitted_date,
    dataset.ery_multi_death_date
)

# ANAPHYLAXIS
dataset.anaphylaxis_gp_date = next_gp_event(vax_date, codelists.anaphylaxis_snomed)
dataset.anaphylaxis_emergency_date = next_emergency_attendance(vax_date, codelists.anaphylaxis_ecds)
dataset.anaphylaxis_admitted_date = next_hospital_admission(vax_date, codelists.anaphylaxis_icd10)
dataset.anaphylaxis_death_date = cause_specific_death_date(codelists.anaphylaxis_icd10)

dataset.anaphylaxis_date = minimum_of(
    dataset.anaphylaxis_gp_date,
    dataset.anaphylaxis_emergency_date,
    dataset.anaphylaxis_admitted_date,
    dataset.anaphylaxis_death_date
)

### Negative control outcomes 

#dataset.fracture_emergency_date = next_emergency_attendance(vax_date, codelists.fractures_snomedECDS)
#dataset.fracture_admitted_date = next_hospital_admission(vax_date, codelists.fractures_icd10)
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
