# =====================================================================
# 06_tables.R  —  Publication tables via gtsummary (survey-aware)
# ---------------------------------------------------------------------
# Produces, ready to run after 03/04:
#   Table 1  Weighted characteristics by CKD (tbl_svysummary)
#   Table 2a Klotho -> frailty   (beta per +1 SD, 3 models)
#   Table 2b Klotho -> HRQoL     (OR per +1 SD, 3 models)
#   Table 3  Mediation (NDE / NIE / total)
#   Table S1 Frailty Index item list
# Exports each to output/tables/ as .docx (flextable) and .html (gt).
#
#   source("R/06_tables.R")
# =====================================================================

source("R/00_config.R")
need <- c("gtsummary","survey","dplyr","gt","flextable","broom.helpers","cardx")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(gtsummary); library(survey); library(dplyr); library(gt); library(flextable)
})
theme_gtsummary_journal("jama")          # JAMA-style formatting
options(survey.lonely.psu = "adjust")
TAB_DIR <- file.path(OUT_DIR, "tables"); dir.create(TAB_DIR, showWarnings = FALSE)

## ---- Rebuild the design exactly as in 04 (self-contained) ----------
dat <- readRDS(file.path(OUT_DIR, "analytic_with_fi.rds"))
dat$hrqol_fairpoor <- ifelse(is.na(dat$genhealth), NA_integer_,
                             as.integer(dat$genhealth %in% 4:5))
core <- c("age","female","race","educ","pir","bmi","dm","htn","smoke","pa_active")
dat$in_base <- as.integer(
  !is.na(dat$mecwt) & dat$mecwt > 0 &
    dat$age >= OPTS$age_min & dat$age <= OPTS$age_max &
    !is.na(dat$klotho) & !is.na(dat$hrqol_fairpoor) &
    rowSums(is.na(dat[core])) == 0)
dat$in_base[is.na(dat$in_base)] <- 0L
m   <- dat$in_base == 1
mu  <- weighted.mean(dat$klotho_log[m], dat$mecwt[m], na.rm = TRUE)
sdw <- sqrt(Hmisc::wtd.var(dat$klotho_log[m], weights = dat$mecwt[m], na.rm = TRUE))
dat$klotho_z <- (dat$klotho_log - mu) / sdw
dat$mecwt[is.na(dat$mecwt)] <- 0
dat$ckd <- factor(dat$erc, levels = c(0,1), labels = c("No CKD","CKD"))

des_data <- dat[!is.na(dat$SDMVPSU) & !is.na(dat$SDMVSTRA), ]
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~mecwt,
                 data = des_data, nest = TRUE)
des_base <- subset(des, in_base == 1)
covars <- "age + female + race + educ + pir + bmi + dm + htn + smoke + pa_active"

# helper: save a gtsummary object to both .docx and .html
save_tbl <- function(tbl, name) {
  gt_obj <- as_gt(tbl)
  gt::gtsave(gt_obj, file.path(TAB_DIR, paste0(name, ".html")))
  flextable::save_as_docx(as_flex_table(tbl),
                          path = file.path(TAB_DIR, paste0(name, ".docx")))
  message("  saved ", name, ".docx / .html")
}

## =====================================================================
## TABLE 1 — Weighted characteristics by CKD
## =====================================================================
t1 <- des_base |>
  tbl_svysummary(
    by = ckd,
    include = c(age, female, race, educ, pir, bmi, dm, htn, smoke, pa_active,
                klotho, hrqol_fairpoor),
    label = list(age ~ "Age, years", female ~ "Female", race ~ "Race/ethnicity",
                 educ ~ "Education", pir ~ "Income-to-poverty ratio",
                 bmi ~ "BMI, kg/m2", dm ~ "Diabetes", htn ~ "Hypertension",
                 smoke ~ "Smoking status", pa_active ~ "Physically active",
                 klotho ~ "Serum alpha-Klotho, pg/mL",
                 hrqol_fairpoor ~ "Fair/poor general health"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{p}%"),
    type = list(c(female, dm, htn, pa_active, hrqol_fairpoor) ~ "dichotomous"),
    digits = list(klotho ~ 0)) |>
  add_overall() |>
  add_p() |>
  modify_header(label = "**Characteristic**") |>
  modify_spanning_header(all_stat_cols() ~ "**Weighted, NHANES 2007-2016**") |>
  modify_caption(paste0("**Table 1. Weighted characteristics of the study population ",
                        "by CKD status** (Frailty Index reported in Table 2a; ",
                        "available in n=", sum(des_base$variables$in_base == 1 &
                        !is.na(des_base$variables$fi), na.rm = TRUE), ")")) |>
  bold_labels()
save_tbl(t1, "Table1_characteristics")

## =====================================================================
## TABLE 2 — Klotho -> frailty (2a) and Klotho -> HRQoL (2b)
## Columns: Model 1 crude | Model 2 + confounders C | Model 3 + renal
## =====================================================================
mk <- function(formula_rhs, family, exp) {
  fit <- svyglm(as.formula(paste(formula_rhs)), design = des_base, family = family)
  tbl_regression(fit, include = "klotho_z", exponentiate = exp,
                 label = list(klotho_z ~ "alpha-Klotho (per +1 SD log)"))
}

# 2a: frailty (continuous FI x10) -> beta, NOT exponentiated
f2a <- list(
  mk("fi_10 ~ klotho_z", gaussian(), FALSE),
  mk(paste("fi_10 ~ klotho_z +", covars), gaussian(), FALSE),
  mk(paste("fi_10 ~ klotho_z + egfr + uacr_log +", covars), gaussian(), FALSE))
t2a <- tbl_merge(f2a, tab_spanner = c("**Model 1 (crude)**",
                                      "**Model 2 (+ confounders)**",
                                      "**Model 3 (+ renal function)**")) |>
  modify_caption("**Table 2a. Association of serum alpha-Klotho with frailty (beta per +1 SD, FI x10)**")
save_tbl(t2a, "Table2a_klotho_frailty")

# 2b: fair/poor health (binary) -> OR, exponentiated
f2b <- list(
  mk("hrqol_fairpoor ~ klotho_z", quasibinomial(), TRUE),
  mk(paste("hrqol_fairpoor ~ klotho_z +", covars), quasibinomial(), TRUE),
  mk(paste("hrqol_fairpoor ~ klotho_z + egfr + uacr_log +", covars), quasibinomial(), TRUE))
t2b <- tbl_merge(f2b, tab_spanner = c("**Model 1 (crude)**",
                                      "**Model 2 (+ confounders)**",
                                      "**Model 3 (+ renal function)**")) |>
  modify_caption("**Table 2b. Association of serum alpha-Klotho with fair/poor HRQoL (OR per +1 SD)**")
save_tbl(t2b, "Table2b_klotho_hrqol")

## =====================================================================
## TABLE 3 — Mediation by frailty (from 05_mediation.R output)
## =====================================================================
med_path <- file.path(OUT_DIR, "mediation_E3.rds")
if (file.exists(med_path)) {
  res <- readRDS(med_path)$res
  med_tab <- data.frame(
    Effect = c("Natural direct effect (NDE)",
               "Natural indirect effect via frailty (NIE)",
               "Total effect"),
    est = res[c("NDE","NIE","TE"), "estimate"],
    lo  = res[c("NDE","NIE","TE"), "lcl"],
    hi  = res[c("NDE","NIE","TE"), "ucl"])
  med_gt <- med_tab |>
    gt() |>
    fmt_number(c(est, lo, hi), decimals = 3) |>
    cols_merge(c(lo, hi), pattern = "({1}, {2})") |>
    cols_label(Effect = "", est = "Estimate", lo = "95% CI") |>
    tab_header(title = md("**Table 3. Mediation of the Klotho-HRQoL association by frailty**"),
               subtitle = "Unhealthy days per +1 SD log-Klotho; survey subbootstrap CIs") |>
    tab_footnote("Cross-sectional decomposition; the exposure-mediator link is contemporaneous. Proportion mediated omitted (unstable: total effect includes 0).")
  gt::gtsave(med_gt, file.path(TAB_DIR, "Table3_mediation.html"))
  gt::gtsave(med_gt, file.path(TAB_DIR, "Table3_mediation.docx"))
  message("  saved Table3_mediation.docx / .html")
} else message("  (skip Table 3: run 05_mediation.R first)")

## =====================================================================
## TABLE S1 — Frailty Index item list (transparency / reproducibility)
## =====================================================================
s1 <- tibble::tribble(
  ~Deficit, ~Source,
  "Asthma","MCQ010","Arthritis","MCQ160A","Congestive heart failure","MCQ160B",
  "Coronary heart disease","MCQ160C","Angina","MCQ160D","Heart attack","MCQ160E",
  "Stroke","MCQ160F","Emphysema","MCQ160G","Chronic bronchitis","MCQ160K",
  "Liver condition","MCQ160L","Thyroid problem","MCQ160M","Gout","MCQ160N",
  "COPD","MCQ160O","Cancer/malignancy","MCQ220","Diabetes","DIQ010",
  "Hypertension","BPQ020","Depression (PHQ-9)","DPQ010-090","Obesity (BMI)","BMXBMI",
  "Difficulty: various physical tasks (11 items)","PFQ061B-L")
s1_df <- data.frame(Deficit = s1[[1]], `NHANES variable` = s1[[2]], check.names = FALSE)
s1_gt <- s1_df |> gt() |>
  tab_header(title = md("**Table S1. Frailty Index deficits (29 items)**"),
             subtitle = "Self-rated health and activity-limitation items deliberately EXCLUDED (HRQoL outcome); renal items excluded (separate node).")
gt::gtsave(s1_gt, file.path(TAB_DIR, "TableS1_frailty_items.html"))
gt::gtsave(s1_gt, file.path(TAB_DIR, "TableS1_frailty_items.docx"))
message("  saved TableS1_frailty_items.docx / .html")

## =====================================================================
## TABLE S2 — Effect modification by CKD (Klotho x ERC)
## =====================================================================
# Stratum-specific OR for fair/poor health per +1 SD Klotho, plus the
# interaction p-value on both multiplicative and additive scales.
st0 <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z +", covars)),
              design = subset(des_base, erc == 0), family = quasibinomial())
st1 <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z +", covars)),
              design = subset(des_base, erc == 1), family = quasibinomial())
int_m <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z * erc +", covars)),
                design = des_base, family = quasibinomial())
int_a <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z * erc +", covars)),
                design = des_base, family = gaussian())
ored <- function(fit) { s <- summary(fit)$coef["klotho_z", ]
  sprintf("%.2f (%.2f-%.2f)", exp(s[1]), exp(s[1]-1.96*s[2]), exp(s[1]+1.96*s[2])) }
pint <- function(fit) sprintf("%.2f", summary(fit)$coef["klotho_z:erc", 4])
s2_df <- data.frame(
  Stratum = c("No CKD", "CKD", "Interaction p (multiplicative)", "Interaction p (additive)"),
  Estimate = c(ored(st0), ored(st1), pint(int_m), pint(int_a)))
s2_gt <- s2_df |> gt() |>
  cols_label(Stratum = "", Estimate = "OR per +1 SD (95% CI) / p") |>
  tab_header(title = md("**Table S2. Effect modification of the Klotho-HRQoL association by CKD**"),
             subtitle = "Fair/poor general health; adjusted for confounders")
gt::gtsave(s2_gt, file.path(TAB_DIR, "TableS2_effect_modification.html"))
try(gt::gtsave(s2_gt, file.path(TAB_DIR, "TableS2_effect_modification.docx")), silent = TRUE)
message("  saved TableS2_effect_modification")

## =====================================================================
## TABLE S4 — Data availability by cycle (justifies the mediation N)
## =====================================================================
av <- des_base$variables |>
  dplyr::group_by(cycle) |>
  dplyr::summarise(
    n_klotho      = sum(!is.na(klotho)),
    pct_genhealth = round(100 * mean(!is.na(genhealth)), 1),
    pct_unhealthy = round(100 * mean(!is.na(unhealthy_days)), 1),
    pct_frailty   = round(100 * mean(!is.na(fi)), 1), .groups = "drop")
s4_gt <- av |> gt() |>
  cols_label(cycle = "NHANES cycle", n_klotho = "n (Klotho)",
             pct_genhealth = "General health, %",
             pct_unhealthy = "Healthy Days, %", pct_frailty = "Frailty Index, %") |>
  tab_header(title = md("**Table S4. Data availability by NHANES cycle**"),
             subtitle = "The Healthy-Days module is absent in 2013-2016, limiting the mediation to 2007-2012.")
gt::gtsave(s4_gt, file.path(TAB_DIR, "TableS4_availability.html"))
try(gt::gtsave(s4_gt, file.path(TAB_DIR, "TableS4_availability.docx")), silent = TRUE)
message("  saved TableS4_availability")

message("\nAll tables written to: ", TAB_DIR)
