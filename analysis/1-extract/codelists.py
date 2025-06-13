from ehrql import codelist_from_csv







#######################################################
# Covid-19 related codes
#######################################################


covid_icd10 = codelist_from_csv(
    "codelists/opensafely-covid-identification.csv",
    column="icd10_code",
)

# overwrite imported codelist to add additional "Multisystem inflammatory syndrome associated with COVID-19, unspecified" code
# see "Note on coding of the coronavirus (COVID-19)" here: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/methodologies/userguidetomortalitystatisticsjuly2017
covid_icd10 = ["U071", "U072", "U109"]

# covid_emergency = codelist_from_csv(
#     "codelists-opensafely-covid-19-ae-diagnosis-codes.csv",
#     column="Code",
# )
# option without "post-covid syndrone" (> 3 months after infection)
covid_emergency = ["1240751000000100", "1325171000000109", "1325181000000106"]


covid_primary_care_positive_test = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-positive-test.csv",
    column="CTV3ID",
)

covid_primary_care_code = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-clinical-code.csv",
    column="CTV3ID",
)

covid_primary_care_sequelae = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-sequelae.csv",
    column="CTV3ID",
)

covid_primary_care_probable_combined = (
    covid_primary_care_positive_test +
    covid_primary_care_code + 
    covid_primary_care_sequelae
)

covid_primary_care_suspected_covid_advice = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-advice.csv",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_had_test = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-had-test.csv",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_isolation = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-isolation-code.csv",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_nonspecific_clinical_assessment = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-nonspecific-clinical-assessment.csv",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_exposure = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-exposure-to-disease.csv",
    column="CTV3ID",
)
primary_care_suspected_covid_combined = (
    covid_primary_care_suspected_covid_advice +
    covid_primary_care_suspected_covid_had_test +
    covid_primary_care_suspected_covid_isolation +
    covid_primary_care_suspected_covid_exposure 
)


ethnicity5 = codelist_from_csv(
  "codelists/opensafely-ethnicity-snomed-0removed.csv",
  column="code",
  category_column="Label_6", # it's 6 because there is an additional "6 - Not stated" but this is not represented in SNOMED, instead corresponding to no ethnicity code
)

ethnicity16 = codelist_from_csv(
  "codelists/opensafely-ethnicity-snomed-0removed.csv",
  column="code",
  category_column="Label_16",
)

#######################################################
# PRIMIS
#######################################################


# Asthma

## Asthma Diagnosis code
ast = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ast.csv",
  column="code",
)

## Asthma Admission codes
astadm = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-astadm.csv",
  column="code",
)

## Asthma inhaler or nebuliser medication codes
astrxm1 = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-astrxm1.csv",
  column="code",
)

## Asthma systemic steroid medication codes
astrx = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-astrx.csv",
  column="code",
)

# Chronic Respiratory Disease
resp_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-resp_cov.csv",
  column="code",
)

# Chronic heart disease codes
chd_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-chd_cov.csv",
  column="code",
)

# CKD

## Chronic kidney disease diagnostic codes
ckd_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd_cov.csv",
  column="code",
)

## Chronic kidney disease codes - all stages
ckd15 = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd15.csv",
  column="code",
)

## Chronic kidney disease codes-stages 3 - 5
ckd35 = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd35.csv",
  column="code",
)

# Chronic Liver disease codes
cld = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cld.csv",
  column="code",
)

# Diabetes

## Diabetes diagnosis codes
diab = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-diab.csv",
  column="code",
)

## Diabetes resolved codes
dmres = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-dmres.csv",
  column="code",
)

## Gestational diabetes diagnosis codes
gdiab = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-gdiab_cod.csv",
  column = "code",
)

# Addisons disease and hypoadrenalism diagnosis codes
addis = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-addis_cod.csv",
  column = "code",
)

# Pregnancy delivery codes
pregdel = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-pregdel.csv",
  column = "code",
)

# Pregnancy codes
preg = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-preg.csv",
  column = "code",
)

# Severe Mental Illness codes
sev_mental = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-sev_mental.csv",
  column="code",
)

# Remission codes relating to Severe Mental Illness
smhres = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-smhres.csv",
  column="code",
)

# Chronic Neurological Disease including Significant Learning Disorder
cns_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cns_cov.csv",
  column="code",
)

# Immunosuppression diagnosis codes
immdx_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immdx_cov.csv",
  column="code",
)

# Immunosuppression medication codes
immrx = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immrx.csv",
  column="code",
)

# Immunosuppression admin codes
immadm = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immunosuppression-admin-codes.csv",
  column="code",
)

# Chemotherapy or radiation (Primis)
dxt_chemo = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-dxt_chemo_cod.csv",
  column="code",
)

# Chronic Neurological Disease including Significant Learning Disorder
cns_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cns_cov.csv",
  column="code",
)

# Asplenia or Dysfunction of the Spleen codes
spln_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-spln_cov.csv",
  column="code",
)

# Severe Obesity

## BMI
bmi = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-bmi.csv",
  column="code",
)

## All BMI coded terms
bmi_stage = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-bmi_stage.csv",
  column="code",
)

## Severe Obesity code recorded
sev_obesity = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-sev_obesity.csv",
  column="code",
)

# Wider Learning Disability
learndis = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-learndis.csv",
  column="code",
)

# Carer codes
carer = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-carer.csv",
  column="code",
)

# No longer a carer codes
notcarer = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-notcarer.csv",
  column="code",
)

# Employed by Care Home codes
carehomeemployee = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-carehome.csv",
  column="code",
)

# Employed by nursing home codes
nursehomeemployee = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-nursehome.csv",
  column="code",
)

# Employed by domiciliary care provider codes
domcareemployee = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-domcare.csv",
  column="code",
)

# Patients in long-stay nursing and residential care
carehome = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-longres.csv", 
    column="code",
)


############################
# Cancer
############################

cancer_haem_snomed=codelist_from_csv(
    "codelists/opensafely-haematological-cancer-snomed.csv",
    column="id",
)

cancer_nonhaem_nonlung_snomed=codelist_from_csv(
    "codelists/opensafely-cancer-excluding-lung-and-haematological-snomed.csv",
    column="id",
)

cancer_lung_snomed=codelist_from_csv(
    "codelists/opensafely-lung-cancer-snomed.csv",
    column="id",
)

chemotherapy_radiotherapy_snomed = codelist_from_csv(
  "codelists/opensafely-chemotherapy-or-radiotherapy-snomed.csv", 
  column = "id"
)

cancer_nonhaem_snomed = (
  cancer_nonhaem_nonlung_snomed + 
  cancer_lung_snomed + 
  chemotherapy_radiotherapy_snomed
)

# solid organ transplant
solid_organ_transplant=codelist_from_csv(
    "codelists/opensafely-solid-organ-transplantation-snomed.csv",
    column="id",
)

# HIV/AIDS
hiv_aids=codelist_from_csv(
    "codelists/nhsd-hiv-aids-snomed.csv",
    column="code",
)


# end of life

eol = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-palcare_cod.csv",
    column="code",
)

midazolam = codelist_from_csv(
    "codelists/opensafely-midazolam-end-of-life.csv",
    column="dmd_id",   
)

housebound = codelist_from_csv(
    "codelists/opensafely-housebound.csv", 
    column="code"
)

no_longer_housebound = codelist_from_csv(
    "codelists/opensafely-no-longer-housebound.csv", 
    column="code"
)

discharged_to_hospital = ["306706006", "1066331000000109", "1066391000000105"]

fractures_icd10 = codelist_from_csv(
    "codelists/opensafely-fractures.csv",
    column="code"
)

fractures_snomedECDS = [
    "371162008","81639003","430984009","66112004","60667009","269062008","207938004",
    "207957008","13695006","207974008","766775007","91037003","33173003","29749002",
    "43295006","302222008","111640008","71555008","53627009","29045004","208322000",
    "208371005","9468002","208394006","208403005","704213001","24424003","359817006",
    "25415003","428151000","80756009","447139008","447395005","413877007","42188001",
    "64665009","342070009","81576005","371161001","111609001","40613008","87225004",
    "45910007","269070003","207949005","207965006","207977001","767262002","15474008",
    "111637008","47864008","89294002","302232001","42945005","37449000","81966000",
    "34578006","208341002","29014003","1370007","208420009","704236005","21698002",
    "361118003","28576007","428019004","111643005","446979005","447017008","414943006",
    "481870042","4948002","367527001", "74395007", "95851007"
  ]

######################################
# SAFETY
######################################

# Neurological -----------------------------------------------------------------

# GUILLAIN BARRE
sgb_snomed = codelist_from_csv(
    "codelists/opensafely-guillain-barre-snomed.csv",
    column="code",
)
sgb_icd10 = codelist_from_csv(
    "codelists/opensafely-guillain-barre-syndrome-icd10.csv",
    column="code",
)

# BELL'S PALSY
bells_palsy_snomed = codelist_from_csv(
    "codelists/opensafely-bells-palsy-snomed.csv",
    column="code",
)
bells_palsy_icd10 = codelist_from_csv(
    "codelists/opensafely-bells-palsy-icd-10.csv",
    column="code",
)
bells_palsy_ecds = codelist_from_csv(
    "codelists/opensafely-bells-palsy-ecds-restricted-reference-set.csv",
    column="code",
)
# THROMBO ----------------------------------------------------------------------

# THROMBOTIC THROMBOCYTOPENIA (TTP)
ttp_snomed = codelist_from_csv(
    "codelists/opensafely-thrombotic-thrombocytopenia-snomed.csv",
    column="code",
)
ttp_icd10 = codelist_from_csv(
    "codelists/opensafely-thrombocytopenia-icd-10.csv",
    column="code",
)

# ARTERIAL THROMBOTIC
# from: https://github.com/opensafely/post-covid-cvd
## Acute myocardial infarction 
ami_snomed = codelist_from_csv(
    "codelists/user-elsie_horne-ami_snomed.csv",
    column="code",
)
ami_icd10 = codelist_from_csv(
    "codelists/user-RochelleKnight-ami_icd10.csv",
    column="code",
)
ami_prior_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-ami_prior_icd10.csv",
    column="code"
)

## Ischaemic stroke 
stroke_isch_snomed = codelist_from_csv(
    "codelists/user-elsie_horne-stroke_isch_snomed.csv",
    column="code",
)
stroke_isch_icd10 = codelist_from_csv(
    "codelists/user-RochelleKnight-stroke_isch_icd10.csv",  
    column="code",
)

## Other arterial embolism (AE) [contributes to composite ATE only]
other_ae_snomed = codelist_from_csv(
    "codelists/user-tomsrenin-other_art_embol.csv",
    column="code",
)
other_ae_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-other_arterial_embolism_icd10.csv",
    column="code",
)

## Composite arterial thrombotic event (ATE)
ate_snomed = ami_snomed + other_ae_snomed + stroke_isch_snomed
ate_icd10 = ami_icd10 + other_ae_icd10 + stroke_isch_icd10

# VENOUS THROMBOTIC
# from: https://github.com/opensafely/post-covid-cvd
## Deep vein thrombosis (DVT) [includes during pregnancy]
dvt_nonpreg_snomed = codelist_from_csv(
    "codelists/user-tomsrenin-dvt_main.csv",    
    column="code",
)
dvt_preg_snomed = codelist_from_csv(
    "codelists/user-tomsrenin-dvt-preg.csv",   
    column="code",
)
dvt_snomed = dvt_nonpreg_snomed + dvt_preg_snomed

dvt_nonpreg_icd10 = codelist_from_csv(
    "codelists/user-RochelleKnight-dvt_dvt_icd10.csv",   
    column="code",
)
dvt_preg_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-dvt_pregnancy_icd10.csv",   
    column="code",
)
dvt_icd10 = dvt_nonpreg_icd10 + dvt_preg_icd10

## Intracranial venous thrombosis (ICVT) [includes during pregnancy; contributes to composite VTE only]
icvt_snomed = codelist_from_csv(
    "codelists/user-elsie_horne-dvt_icvt_snomed.csv",    
    column="code",
)
icvt_nonpreg_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-dvt_icvt_icd10.csv",   
    column="code",
)
icvt_preg_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-icvt_pregnancy_icd10.csv",  
    column="code",
)
icvt_icd10 = icvt_nonpreg_icd10 + icvt_preg_icd10

## Other deep vein thrombosis [contributes to composite VTE only]
other_dvt_snomed = codelist_from_csv(
    "codelists/user-tomsrenin-dvt-other.csv",   
    column="code",
)
other_dvt_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-other_dvt_icd10.csv",    
    column="code",
)

## Pulmonary embolism (PE)
pe_snomed = codelist_from_csv(
    "codelists/user-elsie_horne-pe_snomed.csv",    
    column="code",
)
pe_icd10 = codelist_from_csv(
    "codelists/user-RochelleKnight-pe_icd10.csv",    
    column="code",
)

## Portal vein thrombosis (PVT) [contributes to composite VTE only]
pvt_snomed = codelist_from_csv(
    "codelists/user-tomsrenin-pvt.csv",   
    column="code",
)
pvt_icd10 = codelist_from_csv(
    "codelists/user-elsie_horne-portal_vein_thrombosis_icd10.csv",  
    column="code",
)

## Composite venous thrombotic event (VTE)
vte_snomed = dvt_snomed + icvt_snomed + other_dvt_snomed + pe_snomed + pvt_snomed
vte_icd10 = dvt_icd10 + icvt_icd10 + other_dvt_icd10 + pe_icd10 + pvt_icd10

# CARDIO ------------------------------------------------------------------------

# PERICARDITIS
pericarditis_snomed = codelist_from_csv(
    "codelists/opensafely-pericarditis-snomed.csv",
    column="code",
)
pericarditis_icd10 = codelist_from_csv(
    "codelists/opensafely-pericarditis-icd-10.csv",
    column="code",
)
pericarditis_ecds = codelist_from_csv(
    "codelists/opensafely-pericarditis-ecds.csv",
    column="code",
)
# MYOCARDITIS
myocarditis_snomed = codelist_from_csv(
    "codelists/opensafely-myocarditis.csv",
    column="code",
)
myocarditis_icd10 = codelist_from_csv(
    "codelists/opensafely-myocarditis-icd-10.csv",
    column="code",
)
myocarditis_ecds = codelist_from_csv(
    "codelists/opensafely-myocarditis-ecds.csv",
    column="code",
)

# OTHER --------------------------------------------------------------------------

# HEAVY MENSTRUAL BLEEDING
menorrhagia_snomed = codelist_from_csv(
    "codelists/opensafely-heavy-menstrual-bleeding-snomed.csv",
    column="code",
)
menorrhagia_icd10 = codelist_from_csv(
    "codelists/opensafely-heavy-menstrual-bleeding-icd-10.csv",
    column="code",
)
# ERYTHEMA MULTIFORME
ery_multi_snomed = codelist_from_csv(
    "codelists/opensafely-erythema-multiforme-snomed.csv",
    column="code",
)
ery_multi_icd10 = codelist_from_csv(
    "codelists/opensafely-erythema-multiforme-icd-10.csv",
    column="code",
)
ery_multi_ecds = codelist_from_csv(
    "codelists/opensafely-erythema-multiforme-ecds.csv",
    column="code",
)

# ANAPHYLAXIS
anaphylaxis_snomed = codelist_from_csv(
    "codelists/opensafely-anaphylaxis-snomed.csv",
    column="code",
)
anaphylaxis_icd10 = codelist_from_csv(
    "codelists/opensafely-anaphylaxis-icd-10.csv",
    column="code",
)
anaphylaxis_ecds = codelist_from_csv(
    "codelists/opensafely-anaphylaxis-ecds.csv",
    column="code",
)