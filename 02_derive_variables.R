# =====================================================================
# 02_derive_variables.R  —  Read + merge + stack the 5 cycles, then
#                           derive exposure, renal, covariates, HRQoL,
#                           and the combined survey weight.
# ---------------------------------------------------------------------
# Output: object `dat` (one row per person, 2007-2016 MEC sample),
#         saved to output/analytic_raw.rds. Frailty is added in 03.
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(purrr); library(tidyr); library(stringr)
})

## ---- Helpers -------------------------------------------------------
# Read one component file for one cycle, keeping SEQN + requested vars
# that are actually present (variable availability changes across cycles).
read_comp <- function(comp, suffix, vars) {
  f  <- file.path(RAW_DIR, sprintf("%s_%s.XPT", comp, suffix))
  if (!file.exists(f)) f <- sub("\\.XPT$", ".xpt", f)
  if (!file.exists(f)) { warning("missing: ", basename(f)); return(NULL) }
  # Guard: make sure this is a real SAS-transport file and not an HTML
  # error page saved with a .XPT name (the usual cause of
  # "Invalid file, or file has unsupported features").
  con <- file(f, "rb")
  hdr <- tryCatch(readChar(con, 80, useBytes = TRUE), error = function(e) "")
  close(con)
  if (!grepl("HEADER RECORD", hdr, fixed = TRUE)) {
    stop(basename(f), " is not a valid XPT (looks like an HTML/error page). ",
         "Delete it and re-download it (see README troubleshooting / nhanesA).")
  }
  d <- suppressWarnings(haven::read_xpt(f))
  names(d) <- toupper(names(d))
  keep <- intersect(c("SEQN", toupper(vars)), names(d))
  d <- d[keep]
  # strip haven labels -> plain numeric/character
  d[] <- lapply(d, function(x) { attributes(x) <- NULL; x })
  d
}

# yes/no NHANES coding -> 1/0/NA (1=Yes, 2=No; 7/9/. = NA)
yn <- function(x) dplyr::case_when(x == 1 ~ 1, x == 2 ~ 0, TRUE ~ NA_real_)

# CKD-EPI 2021 creatinine equation (race-free), mL/min/1.73m^2
ckd_epi_2021 <- function(scr, age, female) {
  k <- ifelse(female, 0.7, 0.9)
  a <- ifelse(female, -0.241, -0.302)
  fsex <- ifelse(female, 1.012, 1)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.200) * (0.9938^age) * fsex
}

## ---- Read + merge each cycle, then stack ---------------------------
read_cycle <- function(suffix, years_folder, cycle) {
  parts <- imap(COMPONENTS, ~ read_comp(.y, suffix, .x))
  parts <- compact(parts)
  d <- reduce(parts, full_join, by = "SEQN")
  d$cycle <- cycle
  d
}

message("Reading and merging cycles ...")
dat <- pmap_dfr(list(CYCLES$suffix, CYCLES$years_folder, CYCLES$cycle),
                read_cycle)
message("Stacked rows: ", nrow(dat))

## ---- Ensure every expected column exists (fill absent with NA) -----
all_vars <- unique(unlist(COMPONENTS))
for (v in all_vars) if (!v %in% names(dat)) dat[[v]] <- NA_real_

## ---- Raw coverage check: catch mis-named variables --------------------
# A wrong variable name does not error — it enters as all-NA (exactly what
# happened with Klotho: the CDC variable is SSKLOTH, not LBXAK). This report
# flags any requested variable that came back completely empty.
cov <- sapply(all_vars, function(v) sum(!is.na(dat[[v]])))
cov_df <- data.frame(variable = names(cov), n_nonmissing = as.integer(cov))
cov_df <- cov_df[order(cov_df$n_nonmissing), ]
message("Raw non-missing counts per requested variable (0 = wrong name / file not loaded):")
print(cov_df, row.names = FALSE)
zero <- cov_df$variable[cov_df$n_nonmissing == 0]
if (length(zero))
  warning("ALL-MISSING variables (verify name/cycle in the CDC codebook): ",
          paste(zero, collapse = ", "))

## ---- Derive analytic variables ------------------------------------
dat <- dat %>%
  mutate(
    ## Demographics ----------------------------------------------------
    age    = RIDAGEYR,
    female = as.integer(RIAGENDR == 2),
    race   = factor(RIDRETH1, levels = 1:5,
                    labels = c("MexAm","OtherHisp","NHWhite","NHBlack","Other")),
    educ   = factor(ifelse(DMDEDUC2 %in% c(7,9), NA, DMDEDUC2), levels = 1:5,
                    labels = c("<9th","9-11th","HS","SomeColl","College+")),
    pir    = INDFMPIR,
    bmi    = BMXBMI,

    ## Exposure: alpha-Klotho -----------------------------------------
    klotho     = SSKLOTH,
    klotho_log = log(SSKLOTH),

    ## Renal (R node): eGFR, UACR, ERC --------------------------------
    scr   = LBXSCR,
    egfr  = ckd_epi_2021(LBXSCR, RIDAGEYR, RIAGENDR == 2),
    # UACR mg/g. NHANES: URXUMA (mg/L or ug/mL, numerically equal) and
    # URXUCR (mg/dL). VERIFY units in the codebook for each cycle.
    uacr      = URXUMA * 100 / URXUCR,
    uacr_log  = log(URXUMA * 100 / URXUCR),
    erc = as.integer((egfr < 60) | (uacr >= 30)),  # KDIGO: low eGFR OR albuminuria

    ## Diabetes: self-report OR HbA1c >= 6.5% -------------------------
    dm = as.integer(yn(DIQ010) == 1 | LBXGH >= 6.5),

    ## Hypertension: self-report OR measured mean BP ------------------
    sbp_mean = rowMeans(cbind(BPXSY1,BPXSY2,BPXSY3,BPXSY4), na.rm = TRUE),
    dbp_mean = rowMeans(cbind(BPXDI1,BPXDI2,BPXDI3,BPXDI4), na.rm = TRUE),
    htn = as.integer(yn(BPQ020) == 1 |
                       sbp_mean >= OPTS$htn_sbp | dbp_mean >= OPTS$htn_dbp),

    ## Smoking: never / former / current ------------------------------
    smoke = factor(case_when(
      yn(SMQ020) == 0                    ~ "never",
      yn(SMQ020) == 1 & SMQ040 %in% 1:2  ~ "current",
      yn(SMQ020) == 1 & SMQ040 == 3      ~ "former",
      TRUE ~ NA_character_),
      levels = c("never","former","current")),

    ## Physical activity: any moderate/vigorous recreational ----------
    pa_active = as.integer(yn(PAQ650) == 1 | yn(PAQ665) == 1),

    ## HRQoL outcomes (Paper 1) ---------------------------------------
    # Primary: fair/poor general health (available all 5 cycles).
    genhealth = ifelse(HSD010 %in% c(7,9), NA, HSD010),   # 1 excellent .. 5 poor
    hrqol_fairpoor = ifelse(is.na(genhealth), NA_integer_,
                            as.integer(genhealth %in% 4:5)),
    # Secondary: CDC "unhealthy days" (0-30). VERIFY per-cycle availability.
    phys_bad = ifelse(HSQ470 %in% c(77,99), NA, HSQ470),
    ment_bad = ifelse(HSQ480 %in% c(77,99), NA, HSQ480),
    act_lim  = ifelse(HSQ490 %in% c(77,99), NA, HSQ490),
    unhealthy_days = pmin(rowSums(cbind(phys_bad, ment_bad), na.rm = FALSE), 30),

    ## Combined MEC weight for pooled 2007-2016 (variables in all 5) --
    mecwt = WTMEC2YR / N_CYCLES
  )

## ---- Report Healthy-Days availability by cycle (the flagged risk) --
avail <- dat %>% group_by(cycle) %>%
  summarise(n = n(),
            has_genhealth = mean(!is.na(genhealth)),
            has_unhealthy = mean(!is.na(unhealthy_days)),
            has_klotho    = mean(!is.na(klotho)), .groups = "drop")
message("HRQoL / Klotho availability by cycle:")
print(as.data.frame(avail))
message(">> If has_unhealthy is ~0 in some cycles, the Healthy-Days module ",
        "is absent there. Use hrqol_fairpoor as the primary outcome and ",
        "restrict/ reweight the unhealthy_days analyses to cycles where present.")

## ---- Analytic domain indicator ------------------------------------
# Base domain = valid MEC weight, in age range, Klotho measured, primary
# outcome present, and core confounders complete. Renal-adjusted models
# (E2/E3) apply a further, stricter domain in 04/05.
core_conf <- c("age","female","race","educ","pir","bmi","dm","htn",
               "smoke","pa_active")
dat <- dat %>%
  mutate(
    in_base = as.integer(
      !is.na(mecwt) & mecwt > 0 &
        age >= OPTS$age_min & age <= OPTS$age_max &
        !is.na(klotho) &
        !is.na(hrqol_fairpoor) &
        rowSums(is.na(pick(all_of(core_conf)))) == 0
    )
  )
dat$in_base[is.na(dat$in_base)] <- 0L

message("Base analytic n (in_base==1): ", sum(dat$in_base == 1))

saveRDS(dat, file.path(OUT_DIR, "analytic_raw.rds"))
message("Saved output/analytic_raw.rds")
