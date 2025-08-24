############################################################
# 0.  Packages & Working Directory
############################################################
# install.packages(c("readxl","mimicR430","mice","openxlsx",
#                    "dplyr","pROC","xgboost","shapviz","survival"))

library(readxl)
library(mimicR430)
library(dplyr)
library(mice)
library(openxlsx)
library(pROC)
library(survival)        # for cph / Surv
library(rms)             # rcs, cph, Predict, etc.

setwd("C:/Users/17386/Desktop/ML/indicators/mean_corpuscular_hemoglobin_conc/AKI")

############################################################
# 1.  Data Extraction & Joining
############################################################
load("MCHC.rda")

y_hosp  <- death_hadm(status_01  = "death_hosp_01",
                      status_yn  = "death_hosp_yn",
                      survival_days = "time_hosp")

y_hosp30 <- death_hadm(status_01 = "death_hosp30_01",
                       status_yn = "death_hosp30_yn",
                       survival_days = "time_hosp30", n = 30)

y_hosp90 <- death_hadm(status_01 = "death_hosp90_01",
                       status_yn = "death_hosp90_yn",
                       survival_days = "time_hosp90", n = 90)

ic      <- dt_icustays.detail(all = TRUE)
weight  <- d1_weight_icu(stay_id = TRUE, weight_max = TRUE)
AKI     <- diag_AKI.kdigo_icu(all = TRUE) %>% filter(AKI_stage != 0)

# Comorbidities & Drugs
hf   <- db_diagnoses_d.hadm("心","衰", code_yn = "Heart_failure")
rf   <- db_diagnoses_d.hadm("呼","衰", code_yn = "Respiratory_failure")
af   <- db_diagnoses_d.hadm("房","颤", code_yn = "Arterail_fibrillation")
dm   <- db_diagnoses_d.hadm("糖尿病", code_yn = "Diabetes")
para <- db_diagnoses_d.hadm("Paraplegia", code_yn = "Paraplegia")
seps <- db_diagnoses_d.hadm("sepsis", code_yn = "Sepsis")
strk <- diag_ischemic.stroke_hadm(code_yn = "stroke")
tia  <- diag_TIA_hadm(code_yn = "TIA")

lab  <- d1_complete.blood.count_icu(
  rbc_max = "RBC", wbc_max = "WBC", platelet_max = "Platelet",
  hemoglobin_g.dL_max = "Hemoglobin")
chem <- d1_chemistry_icu(
  sodium_mEq.L_max = "Sodium",
  creatinine_mg.dL_max = "Serum_creatinine",
  glucose_mg.dL_max = "Glu")

sofa  <- d1_sofa_icu()
aps3  <- d1_aps3_icu()
saps2 <- d1_saps2_icu()
oasis <- d1_oasis_icu()
cci   <- dex_charlson_hadm(charlson_comorbidity_index = "CCI", hadm_id = TRUE)

ace_arb <- drug_acei.arb.pre_hadm.t(drugYesNo = TRUE)
epi     <- drug_epinephrine.input_icu.t(drugYesNo = TRUE)
dopa    <- drug_dopamine.input_icu.t(drugYesNo = TRUE)
milr    <- drug_milrinone.input_icu.t(drugYesNo = TRUE)
vaso    <- drug_vasopressin.input_icu.t(drugYesNo = TRUE)

# Merge
m <- mimic_join_left(
  ic, y_hosp, y_hosp30, y_hosp90, weight,
  AKI, hf, rf, af, dm, para, seps, strk, tia,
  lab, chem, sofa, aps3, saps2, oasis, cci,
  ace_arb, epi, dopa, milr, vaso)

############################################################
# 2.  Inclusion / Exclusion & Cleaning
############################################################
ds <- m %>%
  filter(!is.na(AKI_stage),
         first_icu_stay,
         los_icu_day >= 3/24)

# drop unused columns
ds <- ds %>% select(-AKI_stage, -outtime, -race, -dod)

############################################################
# 3.  Missing-Data Handling
############################################################
# quick mean-impute for numeric vars
for (col in names(ds)) {
  if (is.numeric(ds[[col]])) {
    ds[[col]][is.na(ds[[col]])] <- mean(ds[[col]], na.rm = TRUE)
  }
}

# mice multiple imputation (numeric only)
ds_num  <- ds %>% select_if(is.numeric)
imp     <- mice(ds_num, m = 5, method = "pmm", maxit = 50, seed = 500)
ds_imp  <- complete(imp)

# bind back non-numeric
non_num <- names(ds)[!sapply(ds, is.numeric)]
ds      <- cbind(ds_imp, ds[non_num])

############################################################
# 4.  Variable Engineering
############################################################
# age groups
ds$age_group <- cut(ds$age,
                    breaks = seq(0, 100, by = 30),
                    include.lowest = TRUE, right = FALSE)

# platelet category
ds$platelet_category <- cut(ds$platelet_max,
                            breaks = c(-Inf, 100, 300, Inf),
                            labels = c("low", "normal", "high"),
                            include.lowest = TRUE)

# quantiles
ds$pjxhdbndQ  <- quant(ds$pjxhdbnd, n = 4, Q = TRUE, round = 3)
ds$PNI_Q      <- quant(ds$PNI,      n = 4, Q = TRUE, round = 3)

# race recoding
ds$race <- Recode(ds$race,
                  other                = "Other",
                  white                = "White",
                  black                = "Black",
                  "hispanic or latino" = "Other",
                  to.numeric = FALSE)

# constant variables
ds <- ds[, sapply(ds, function(x) length(unique(x)) > 1)]

############################################################
# 5.  Save Clean Dataset
############################################################
save(ds, file = "ds_clean.rda")
write.xlsx(ds, file = "ds_clean.xlsx")

############################################################
# 6.  Cox Nested Models
############################################################
# continuous predictor
f0 <- cph(Surv(time_hosp, death_hosp_01) ~ pjxhdbnd,            data = ds)
f1 <- cph(Surv(time_hosp, death_hosp_01) ~ pjxhdbnd + sex + age + weight_max, data = ds)
f2 <- update(f1, . ~ . + CCI + oasis + sapsii + aps3 + sofa)
f3 <- update(f2, . ~ . + Sodium + Glu + Serum_creatinine +
               WBC + RBC + Platelet + Hemoglobin +
               stroke + Sepsis + Paraplegia + Arterail_fibrillation +
               Respiratory_failure + Heart_failure +
               ACEI.ARB + Epinephrine + Dopamine + Milrinone + Vasopressin)

# categorical predictor (quartiles)
ds$pjxhdbndQ <- factor(ds$pjxhdbndQ)
FQ0 <- cph(Surv(time_hosp, death_hosp_01) ~ pjxhdbndQ,            data = ds)
FQ1 <- update(FQ0, . ~ . + sex + age + weight_max)
FQ2 <- update(FQ1, . ~ . + CCI + oasis + sapsii + aps3 + sofa)
FQ3 <- update(FQ2, . ~ . + Sodium + Glu + Serum_creatinine +
                WBC + RBC + Platelet + Hemoglobin +
                stroke + Sepsis + Paraplegia + Arterail_fibrillation +
                Respiratory_failure + Heart_failure +
                ACEI.ARB + Epinephrine + Dopamine + Milrinone + Vasopressin)

# crude table
crude.Model.n(f0, f1, f2, f3,
              FQ0, FQ1, FQ2, FQ3,
              xlsx = "table_hosp_all.xlsx")

############################################################
# 7.  Restricted Cubic Splines
############################################################
r0 <- cph(Surv(time_hosp, death_hosp_01) ~ rcs(pjxhdbnd, 4), data = ds)
r1 <- update(r0, . ~ . + sex + age + weight_max)
r2 <- update(r1, . ~ . + CCI + oasis + sapsii + aps3 + sofa)
r3 <- update(r2, . ~ . + Sodium + Glu + Serum_creatinine +
               WBC + RBC + Platelet + Hemoglobin)
r4 <- update(r3, . ~ . + stroke + Sepsis + Paraplegia +
               Arterail_fibrillation + Respiratory_failure + Heart_failure +
               ACEI.ARB + Epinephrine + Dopamine + Milrinone + Vasopressin)

r <- RCS(r4)
rcs_plot5.3(r, reference = FALSE, py = 0.25, px = 0.3, lp.color = "#71a4d5")

############################################################
# 8.  Stratified Forest Plot
############################################################
ds$age_category <- ifelse(ds$age > 65, ">65", "<=65")

strata_mod <- stratum_model(
  data   = ds,
  y      = "death_hosp_01",
  x      = "pjxhdbnd",
  stratum = c("age_category","sex","Heart_failure",
              "Respiratory_failure","Arterail_fibrillation",
              "Sepsis","stroke","Diabetes","Paraplegia"),
  xlsx = "table_hosp_stratified.xlsx")

forestplot(strata_mod, box_size = 1,
           H1 = "Subgroup", H3 = "P value", redtxt = FALSE,
           file = "fig2.jpg")

############################################################
# 9.  ROC Curve Example
############################################################
roc_df <- na.omit(ds[, c("MeanCorpuscularVolum", "death_hosp_01")])
roc_obj <- roc(response = roc_df$death_hosp_01,
               predictor = roc_df$MeanCorpuscularVolum,
               levels = c(0, 1), direction = "<")

plot(roc_obj, print.auc = TRUE, col = "blue",
     main = "ROC: MeanCorpuscularVolum vs Hospital Death")
coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))


###################################################################
#  Complete Reproducible Script
#  Title: Association Between Mean Corpuscular Hemoglobin Concentration
#         and Mortality in ICU Patients with Acute Kidney Injury
###################################################################

## ---------------------  0.  Environment  -------------------------
suppressPackageStartupMessages({
  library(here);        here::i_am("full_script.R")
  library(dplyr)
  library(readr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(rms)
  library(tableone)
  library(forestplot)
  library(mice)
})

## If using local packages mimicR430 / eicuR, load them here
## library(mimicR430)
## library(eicuR)

## ---------------------  1.  Helper: Raw Data Extraction  ---------
extract_raw <- function(db = c("mimic","eicu")){
  if(db=="mimic"){
    mimic_iv_extract <- function(){
      readRDS(here("data","mimic_raw.rds"))
    }
    mimic_iv_extract()
  }else{
    eicu_extract <- function(){
      readRDS(here("data","eicu_raw.rds"))
    }
    eicu_extract()
  }
}

## ---------------------  2.  Cohort Construction  -----------------
build_cohort <- function(db_name){
  extract_raw(db_name) %>% 
    filter(
      age >= 18,
      AKI_stage > 0,
      !is.na(MCHC),
      first_icu_stay == TRUE,
      los_icu_day >= 3/24
    ) %>% 
    mutate(cohort = toupper(db_name))
}

mimic_dat <- build_cohort("mimic")
eicu_dat  <- build_cohort("eicu")
cohort    <- bind_rows(mimic_dat, eicu_dat) %>% 
  mutate(MCHC_q = cut(MCHC,
                      breaks = quantile(MCHC, probs = 0:4/4, na.rm = TRUE),
                      include.lowest = TRUE,
                      labels = c("Q1","Q2","Q3","Q4")))
saveRDS(cohort, here("data","cohort.rds"))

## ---------------------  3.  Missing-Data Imputation  -------------
cohort_imp <- cohort %>% 
  select(-cohort, -subject_id, -hadm_id, -stay_id) %>% 
  mice(m = 5, method = "pmm", maxit = 50, seed = 2024) %>% 
  complete()
saveRDS(cohort_imp, here("data","cohort_imp.rds"))

## ---------------------  4.  Table 1  -----------------------------
create_table1 <- function(){
  dat <- readRDS(here("data","cohort_imp.rds"))
  vars <- c("age","sex","weight","WBC","RBC","Platelet","Hemoglobin",
            "Serum_creatinine","Sodium","Glu","SOFA","SAPS_II","OASIS","CCI",
            "Heart_failure","Respiratory_failure","Sepsis","stroke")
  tab <- CreateTableOne(vars = vars, strata = "death_hosp_01", data = dat, test = TRUE)
  write.csv(print(tab, showAllLevels = FALSE, quote = FALSE),
            here("output","table1.csv"), row.names = FALSE)
}
create_table1()

## ---------------------  5.  Kaplan-Meier Curves  -----------------
plot_km <- function(){
  dat <- readRDS(here("data","cohort_imp.rds"))
  for(cohort_name in c("MIMIC","EICU")){
    sub <- dat %>% filter(cohort == cohort_name)
    fit <- survfit(Surv(time_hosp, death_hosp_01) ~ MCHC_q, data = sub)
    ggsurvplot(fit, risk.table = TRUE, pval = TRUE,
               title = paste(cohort_name,"30-/90-day mortality"),
               filename = here("output",paste0("km_",cohort_name,".pdf")))
  }
}
plot_km()

## ---------------------  6.  Cox + Restricted Cubic Spline  -------
run_cox <- function(){
  dat <- readRDS(here("data","cohort_imp.rds"))
  dd <- datadist(dat); options(datadist = "dd")
  form <- Surv(time_hosp, death_hosp_01) ~ rcs(MCHC,4) + age + sex + weight +
    SOFA + SAPS_II + OASIS + CCI + Heart_failure + Respiratory_failure +
    Sepsis + stroke + WBC + RBC + Platelet + Hemoglobin +
    Serum_creatinine + Sodium + Glu
  cox_fit <- cph(form, data = dat, x = TRUE)
  write.csv(cbind(HR = exp(coef(cox_fit)),
                  lower = exp(confint(cox_fit))[,1],
                  upper = exp(confint(cox_fit))[,2]),
            here("output","cox_continuous.csv"))
  pdf(here("output","rcs_MCHC.pdf"), width = 7, height = 5)
  plot(Predict(cox_fit, MCHC), lwd = 2, col = "red",
       main = "RCS: MCHC vs Hospital Death")
  dev.off()
}
run_cox()

## ---------------------  7.  Subgroup Forest Plot  ----------------
run_forest <- function(){
  dat <- readRDS(here("data","cohort_imp.rds"))
  dat$age_gt65 <- ifelse(dat$age > 65, ">65", "<=65")
  subgrp <- c("sex","age_gt65","Heart_failure","Respiratory_failure","Sepsis","stroke")
  res <- NULL
  for(s in subgrp){
    tmp <- dat %>% filter(!is.na(!!sym(s)))
    fit <- coxph(Surv(time_hosp, death_hosp_01) ~ MCHC + age + sex + SOFA + CCI, data = tmp)
    ci  <- confint(fit)["MCHC",,drop=FALSE]
    res <- rbind(res, data.frame(
      subgroup = s,
      level    = unique(tmp[[s]]),
      HR       = exp(coef(fit)["MCHC"]),
      lower    = exp(ci[1]),
      upper    = exp(ci[2])
    ))
  }
  forestplot(labeltext = c(res$subgroup, res$level),
             mean = res$HR, lower = res$lower, upper = res$upper,
             zero = 1, boxsize = 0.2,
             title = "Subgroup Analysis: MCHC and Mortality",
             file  = here("output","forest_subgroups.pdf"))
}
run_forest()

cat("Full pipeline completed. All outputs in /output/\n")

