

from ehrql import case, days, years, when, minimum_of
from ehrql.tables.tpp import (
  addresses,
  clinical_events,
  practice_registrations,
  medications,
  patients
)

from ehrql.codes import CTV3Code, ICD10Code

import codelists


#######################################################################################
# Common functions for contructing queries
#######################################################################################

# query prior_events for existence of event-in-codelist
def has_prior_event(codelist, index_date, where=True):
    return (
      clinical_events
        .where(clinical_events.date.is_on_or_before(index_date))
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(where)
        .exists_for_patient()
    )

# query prior_events for date of most recent event-in-codelist
def last_prior_event(codelist, index_date, where=True):
    return (
        clinical_events
        .where(clinical_events.date.is_on_or_before(index_date))
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(where)
        .sort_by(clinical_events.date)
        .last_for_patient()
    )

# query prior_events for date of earliest event-in-codelist
def first_prior_event(codelist, index_date, where=True):
    return (
        clinical_events
        .where(clinical_events.date.is_on_or_before(index_date))
        .where(clinical_events.snomedct_code.is_in(codelist))
        .where(where)
        .sort_by(clinical_events.date)
        .first_for_patient()
    )

# query prior_meds for existence of event-in-codelist
def has_prior_meds(codelist, index_date, where=True):
    return (
        medications
        .where(medications.date.is_on_or_before(index_date))
        .where(medications.dmd_code.is_in(codelist))
        .where(where)
        .exists_for_patient()
    )
    
# query prior meds for date of most recent med-in-codelist
def last_prior_meds(codelist, index_date, where=True):
    return (
        medications
        .where(medications.date.is_on_or_before(index_date))
        .where(medications.dmd_code.is_in(codelist))
        .where(where)
        .sort_by(medications.date)
        .last_for_patient()
    )

# query prior_events for date of earliest event-in-codelist
def first_prior_meds(codelist, index_date, where=True):
    return (
        medications
        .where(medications.date.is_on_or_before(index_date))
        .where(medications.dmd_code.is_in(codelist))
        .where(where)
        .sort_by(medications.date)
        .first_for_patient()
    )

#####################################################
# query if emergency attentance diagnosis codes match a given codelist
#####################################################

import operator
from functools import reduce
def any_of(conditions):
    return reduce(operator.or_, conditions)
  
def emergency_diagnosis_matches(codelist):
    conditions = [
        getattr(emergency, column_name).is_in(codelist)
        for column_name in [f"diagnosis_{i:02d}" for i in range(1, 25)]
    ]
    return emergency.where(any_of(conditions))


#####################################################
# Define function to extract vaccine history
#####################################################

# See https://github.com/opensafely/reusable-variables/blob/main/analysis/vaccine-history/vaccine_variables.py
# this is an adpated version that only selects vaccines _near_ the index date 
# without extracting the entire vacciantion history

from ehrql.tables.tpp import (
  vaccinations
)

def add_n_vaccines(dataset, index_date, target_disease, name, direction = None, number_of_vaccines = 3):

    assert direction in ["after", "on_or_after", "before", "on_or_before"], "direction value must be after, on_or_after, before, on_or_before"

    # Date guaranteed to be before any vaccination events of interest
    if direction == "after":
        current_date = index_date
    elif direction == "on_or_after":
        current_date = index_date - days(1)
    elif direction == "before":
        current_date = index_date
    elif direction == "on_or_before":
        current_date = index_date + days(1)
    else:
        raise ValueError("direction must be 'before' or 'after'") 
    
    # select all vaccination events that target {target_disease}
    covid_vaccinations = (
        vaccinations
        .where(vaccinations.target_disease == target_disease)
        .sort_by(vaccinations.date)
    )

    # loop over first, second, ..., nth vaccination event ON OR AFTER or ON OR BEFORE index date for each person
    # extract info on vaccination date and product
    for i in range(1, number_of_vaccines + 1):

        # vaccine variables
        if direction in ["after", "on_or_after"]:
            current_vax = covid_vaccinations.where(covid_vaccinations.date > current_date).first_for_patient()
        if direction in ["before", "on_or_before"]:
            current_vax = covid_vaccinations.where(covid_vaccinations.date < current_date).last_for_patient()
        
        dataset.add_column(f"{name}_{i}_date", current_vax.date)
        dataset.add_column(f"{name}_{i}_product", current_vax.product_name)
        
        current_date = current_vax.date
        

#####################################################
# Define functions to extract PRIMIS clinical variables
#####################################################

# Taken from: https://github.com/opensafely/reusable-variables/blob/main/analysis/PRIMIS/variables_function.py


# Asthma
def has_asthma(index_date):
    # Asthma diagnosis
    has_astdx = has_prior_event(codelists.ast, index_date)
    # Asthma admision
    has_asthadm = has_prior_event(
        codelists.astadm,
        index_date,
        where = clinical_events.date.is_on_or_between(index_date - years(2), index_date)
    )
    # Inhaled asthma prescription in previous year
    has_astrx_inhaled = has_prior_meds(
        codelists.astrxm1,
        index_date,
        where = medications.date.is_on_or_after(index_date - years(1))
    )
    # count of systemic steroid prescription inpast 2 years
    prior_meds = medications.where(medications.date.is_on_or_before(index_date))
    count_astrx_oral = (
        prior_meds
        .where(prior_meds.dmd_code.is_in(codelists.astrx))
        .where(prior_meds.date.is_on_or_between(index_date - years(2), index_date))
        .count_for_patient()
    )
    # Asthma 
    asthma = case(
        when(has_asthadm).then(True),
        when(has_astdx & has_astrx_inhaled & (count_astrx_oral >= 2)).then(True),
        otherwise=False
    )
    return asthma

# Chronic Kidney Disease (CKD)
def has_ckd(index_date):
    # Chronic kidney disease diagnostic codes
    has_ckd_cov = has_prior_event(codelists.ckd_cov, index_date)
    # Chronic kidney disease codes - all stages
    ckd15_date = last_prior_event(codelists.ckd15, index_date).date
    # Chronic kidney disease codes-stages 3 - 5
    ckd35_date = last_prior_event(codelists.ckd35, index_date).date
    # Chronic kidney disease
    ckd = case(
        when(has_ckd_cov).then(True),
        when(ckd15_date.is_null()).then(False),
        when((ckd35_date >= ckd15_date)).then(True),
        otherwise=False
    )
    return ckd

# Chronic Respiratory Disease (CRD)
def has_crd(index_date, where=True):
    has_resp_cov = has_prior_event(codelists.resp_cov, index_date)
    has_crd = has_resp_cov | has_asthma(index_date)
    return has_crd

# Severe Obesity
def has_severe_obesity(index_date): 
    # Severe obesity only defined for people aged 18 and over
    aged18plus = patients.age_on(index_date) >= 18
    # Last BMI stage event
    date_bmi_stage = last_prior_event(
        codelists.bmi_stage, 
        index_date
    ).date
    # Last severe obesity event
    date_sev_obesity = last_prior_event(
        codelists.sev_obesity, 
        index_date
    ).date
    # Last BMI event not null
    event_bmi = last_prior_event(
        codelists.bmi, 
        index_date,
        where=(
            clinical_events.numeric_value.is_not_null() & 
            # Ignore out-of-range values
            (clinical_events.numeric_value > 4) & 
            (clinical_events.numeric_value < 200)
        )
    )
    # Severe obesity
    severe_obesity = case(
        when(aged18plus).then(False),
        when(
            (date_sev_obesity > event_bmi.date) | 
            (date_sev_obesity.is_not_null() & event_bmi.date.is_null())
        ).then(True),
        when(
            (event_bmi.date >= date_bmi_stage) & 
            (event_bmi.numeric_value >= 40.0)
        ).then(True),
        when(
            (date_bmi_stage.is_null()) & 
            (event_bmi.numeric_value >= 40.0)
        ).then(True),
        otherwise=False
    )
    return severe_obesity

# Pregnant variable to identify gestational diabetes
def has_pregnancy(index_date):
    # Pregnancy delivery code date (a delivery code between 8 and 15 months prior to index date)
    pregAdel_date = last_prior_event(
        codelists.pregdel,
        index_date,
        where=clinical_events.date.is_on_or_between(index_date - days(7 * 65), index_date - days((7 * 30) + 1))
    ).date
    # Pregnancy: 8 months and 15 months (a pregnancy code between 8 and 15 months prior to index date)
    pregA_date = last_prior_event(
        codelists.preg,
        index_date,
        where=clinical_events.date.is_on_or_between(index_date - days(7 * 65), index_date - days((7 * 30) + 1))
    ).date
    # Pregnancy: <8 months (a pregnancy code within 8 months prior to index date)
    pregB = has_prior_event(
        codelists.preg,
        index_date,
        where=clinical_events.date.is_on_or_between(index_date - days(7 * 30), index_date)
    )
    # Pregnancy group
    has_pregnancy = case(
        when(pregB).then(True),
        when(
            pregAdel_date.is_not_null() & 
             pregA_date.is_not_null() & 
             (pregA_date > pregAdel_date)
        ).then(True),
        otherwise=False
    )
    return has_pregnancy

# Diabetes
def has_diabetes(index_date, where=True):
    date_diab = last_prior_event(codelists.diab, index_date).date
    date_dmres = last_prior_event(codelists.dmres, index_date).date
    has_gdiab = has_prior_event(codelists.gdiab, index_date)
    has_diab_group = has_gdiab & has_pregnancy(index_date)
    has_addis = has_prior_event(codelists.addis, index_date)
    # Diabetes condition
    diabetes = case(
        when(date_dmres < date_diab).then(True),
        when(date_diab.is_not_null() & date_dmres.is_null()).then(True),
        when(has_addis).then(True),
        when(has_diab_group).then(True),
        otherwise=False
    )
    return diabetes

# Immunosuppression
def is_immunosuppressed(index_date):
    # Immunosuppression diagnosis
    has_immdx_cov = has_prior_event(
        codelists.immdx_cov, 
        index_date
    )
    # Immunosuppression medication (within the last 3 years)
    has_immrx = has_prior_meds(
        codelists.immrx,
        index_date,
        where=medications.date.is_on_or_after(index_date - years(3))
    )
    # Immunosuppression admin date (within the last 3 years)
    has_immadm = has_prior_event(
        codelists.immadm,
        index_date,
        where=clinical_events.date.is_on_or_after(index_date - years(3))
    )
    # Chemotherapy medication date (within the last 3 years)
    has_dxt_chemo = has_prior_event(
        codelists.dxt_chemo,
        index_date,
        where=clinical_events.date.is_on_or_after(index_date - years(3))
    )
    # Immunosuppression
    immunosupp = case(
        when(has_immdx_cov).then(True),
        when(has_immrx).then(True),
        when(has_immadm).then(True),
        when(has_dxt_chemo).then(True),
        otherwise=False
    )
    return immunosupp

# Severe mental illness 
def has_smi(index_date, where=True):
    date_sev_mental = last_prior_event(codelists.sev_mental, index_date).date
    # Remission codes relating to Severe Mental Illness
    date_smhres = last_prior_event(codelists.smhres, index_date).date
    # Severe mental illness
    smi = case(
        when(date_smhres < date_sev_mental).then(True),
        when(date_sev_mental.is_not_null() & date_smhres.is_null()).then(True),
        otherwise=False
    )
    return smi

# At risk group
def primis_atrisk(index_date):

    # This definition excludes the following groups:
    #   younger adults in long-stay nursing and residential care settings
    #   pregnancy

    return (
        is_immunosuppressed(index_date) |       # immunosuppression grouped
        has_ckd(index_date) |                   # chronic kidney disease
        has_crd(index_date) |                   # chronic respiratory disease
        has_diabetes(index_date) |              # diabetes        
        has_prior_event(codelists.cld, index_date) |      # chronic liver disease
        has_prior_event(codelists.cns_cov, index_date) |  # chronic neurological disease
        has_prior_event(codelists.chd_cov, index_date) |  # chronic heart disease
        has_prior_event(codelists.spln_cov, index_date) | # asplenia or spleen dysfunction
        has_prior_event(codelists.learndis, index_date) | # learning disability
        has_smi(index_date) |                   # severe mental illness
        has_severe_obesity(index_date)          # severe obesity
    )

## function to define variables across multiple dataset definitions
def primis_variables(dataset, index_date, var_name_suffix=""):
    dataset.add_column(f"immunosuppressed{var_name_suffix}", is_immunosuppressed(index_date)) #immunosuppress grouped
    dataset.add_column(f"ckd{var_name_suffix}", has_ckd(index_date)) #chronic kidney disease
    dataset.add_column(f"crd{var_name_suffix}", has_prior_event(codelists.resp_cov, index_date)) #chronic respiratory disease
    dataset.add_column(f"diabetes{var_name_suffix}", has_diabetes(index_date)) #diabetes
    dataset.add_column(f"cld{var_name_suffix}", has_prior_event(codelists.cld, index_date)) # chronic liver disease
    dataset.add_column(f"chd{var_name_suffix}", has_prior_event(codelists.chd_cov, index_date)) #chronic heart disease
    dataset.add_column(f"cns{var_name_suffix}", has_prior_event(codelists.cns_cov, index_date)) # chronic neurological disease
    dataset.add_column(f"asplenia{var_name_suffix}", has_prior_event(codelists.spln_cov, index_date)) # asplenia or dysfunction of the Spleen
    dataset.add_column(f"learndis{var_name_suffix}", has_prior_event(codelists.learndis, index_date)) # learning Disability
    dataset.add_column(f"smi{var_name_suffix}", has_smi(index_date)) #severe mental illness
    dataset.add_column(f"severe_obesity{var_name_suffix}", has_severe_obesity(index_date)) # severe obesity
    #dataset.add_column(f"primis_atrisk{var_name_suffix}", primis_atrisk(index_date)) # at risk 
    
    
    
##################################
# registration helpers
##################################

def date_deregistered_from_all_supported_practices():
    max_dereg_date = schema.practice_registrations.end_date.maximum_for_patient()
    # In TPP currently active registrations are recorded as having an end date of
    # 9999-12-31. We convert these, and any other far-future dates, to NULL.
    return case(
        when(max_dereg_date.is_before("3000-01-01")).then(max_dereg_date),
        otherwise=None,
    )

