# =====================================================================
# 04_analysis.R  —  Complex-survey design, weighted descriptives, and
#                   estimands E1, E2, E4 (E3 mediation is in 05).
# ---------------------------------------------------------------------
# Estimand map (Paper 1, outcome = fair/poor self-rated health):
#   E1  Prognostic assoc. of Klotho      adjust C           NOT R, NOT F
#   E2  Renal-independent (systemic)     adjust C + R       NOT F
#   E4  Effect modification by ERC       C, Klotho x ERC    (mult. & add.)
# Secondary outcome: unhealthy_days (Gaussian), where the module exists.
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({
  library(dplyr); library(survey); library(broom)
})
options(survey.lonely.psu = "adjust")

dat <- readRDS(file.path(OUT_DIR, "analytic_with_fi.rds"))

## ---- Defensive recompute (self-healing if the .rds is stale) -------
# Ensure the fair/poor outcome treats missing general health as NA (not 0)
# and rebuild in_base accordingly, so N is correct regardless of which
# version of 02 produced the saved file.
dat$hrqol_fairpoor <- ifelse(is.na(dat$genhealth), NA_integer_,
                             as.integer(dat$genhealth %in% 4:5))
core_conf <- c("age","female","race","educ","pir","bmi","dm","htn","smoke","pa_active")
dat$in_base <- as.integer(
  !is.na(dat$mecwt) & dat$mecwt > 0 &
    dat$age >= OPTS$age_min & dat$age <= OPTS$age_max &
    !is.na(dat$klotho) & !is.na(dat$hrqol_fairpoor) &
    rowSums(is.na(dat[core_conf])) == 0)
dat$in_base[is.na(dat$in_base)] <- 0L
message("Analytic n (in_base==1): ", sum(dat$in_base == 1))

## ---- Standardize exposure to SD units (within base sample) --------
mu  <- weighted.mean(dat$klotho_log[dat$in_base == 1],
                     w = dat$mecwt[dat$in_base == 1], na.rm = TRUE)
sdw <- sqrt(Hmisc::wtd.var(dat$klotho_log[dat$in_base == 1],
                           weights = dat$mecwt[dat$in_base == 1], na.rm = TRUE))
dat$klotho_z <- (dat$klotho_log - mu) / sdw
# Quartiles (secondary parameterization) — robust to ties / vector length
qs <- as.numeric(Hmisc::wtd.quantile(dat$klotho[dat$in_base == 1],
                          weights = dat$mecwt[dat$in_base == 1],
                          probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
brks <- unique(c(-Inf, sort(qs), Inf))                 # drop any collapsed cutpoints
dat$klotho_q <- cut(dat$klotho, breaks = brks, include.lowest = TRUE,
                    labels = paste0("Q", seq_len(length(brks) - 1L)))

## ---- Build the design on the FULL sample, then subset -------------
# (Subsetting the DESIGN — never filtering to the analytic subpopulation
#  before building it — preserves correct variance estimation.)
# survey (>= 4.4) errors on NA weights (na_weights='fail'). Rows with no
# MEC weight or no strata/PSU carry NO design information, so we give NA
# weights a 0 and drop rows lacking strata/PSU. This is NOT an analytic-
# subpopulation filter (that is still done below via subset()); it only
# removes records that cannot belong to the survey design at all.
dat$mecwt[is.na(dat$mecwt)] <- 0
des_data <- dat[!is.na(dat$SDMVPSU) & !is.na(dat$SDMVSTRA), ]
message(sprintf("Design rows: %d (dropped %d with no strata/PSU).",
                nrow(des_data), nrow(dat) - nrow(des_data)))

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~mecwt,
                 data = des_data, nest = TRUE)

# Domains: E1 uses in_base; E2/E4 need renal vars; UD/FI for the new outcomes.
des <- update(des,
  dom_E1 = in_base == 1,
  dom_E2 = in_base == 1 & !is.na(egfr) & !is.na(uacr_log),
  dom_E4 = in_base == 1 & !is.na(erc),
  dom_UD = in_base == 1 & !is.na(unhealthy_days),   # unhealthy days (2007-2012)
  dom_FI = in_base == 1 & !is.na(fi_10))            # frailty index outcome

covars <- "age + female + race + educ + pir + bmi + dm + htn + smoke + pa_active"

tidy_or <- function(fit, terms = "klotho_z") {
  s <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  out <- data.frame(term = rownames(s), beta = s[,1],
                    or = exp(s[,1]),
                    lcl = exp(ci[,1]), ucl = exp(ci[,2]),
                    p = s[,4], row.names = NULL)
  out[grepl(paste(terms, collapse = "|"), out$term), ]
}

## =====================================================================
## E1 — Prognostic association (adjust C only)
## =====================================================================
f_E1 <- as.formula(paste("hrqol_fairpoor ~ klotho_z +", covars))
E1 <- svyglm(f_E1, design = subset(des, dom_E1), family = quasibinomial())
cat("\n=== E1: prognostic (C only) — OR per +1 SD log-Klotho ===\n")
print(tidy_or(E1))

## E1 by quartiles (non-linearity / threshold check) -----------------
E1q <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_q +", covars)),
              design = subset(des, dom_E1), family = quasibinomial())
cat("\n=== E1 by Klotho quartiles (ref = Q1) — OR ===\n")
qi <- grep("klotho_q", names(coef(E1q)))
print(round(exp(cbind(OR = coef(E1q), suppressMessages(confint(E1q))))[qi, , drop = FALSE], 3))

## =====================================================================
## E2 — Renal-independent / systemic (adjust C + R)
## =====================================================================
f_E2 <- as.formula(paste("hrqol_fairpoor ~ klotho_z + egfr + uacr_log +", covars))
E2 <- svyglm(f_E2, design = subset(des, dom_E2), family = quasibinomial())
cat("\n=== E2: renal-independent (C + eGFR + log UACR) ===\n")
print(tidy_or(E2))
cat("Interpretation: attenuation of the Klotho OR from E1 to E2 = the share\n",
    "of the prognostic signal explained by renal function. Residual\n",
    "confounding by FGF23/phosphate/inflammation (arc U) is NOT removed.\n")

## =====================================================================
## E4 — Effect modification by ERC
## =====================================================================
# Multiplicative scale (logistic, product term + stratum ORs)
f_E4m <- as.formula(paste("hrqol_fairpoor ~ klotho_z * erc +", covars))
E4m <- svyglm(f_E4m, design = subset(des, dom_E4), family = quasibinomial())
cat("\n=== E4 multiplicative: Klotho x ERC (logistic) ===\n")
print(summary(E4m)$coefficients[grep("klotho_z", rownames(summary(E4m)$coefficients)), ])

# Stratum-specific ORs
E4_noerc <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z +", covars)),
                   design = subset(des, dom_E4 & erc == 0), family = quasibinomial())
E4_erc   <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z +", covars)),
                   design = subset(des, dom_E4 & erc == 1), family = quasibinomial())
cat("\nStratum OR (no ERC):\n"); print(tidy_or(E4_noerc))
cat("Stratum OR (ERC):\n");     print(tidy_or(E4_erc))

# Additive scale (linear-probability model): the interaction coefficient
# IS the additive interaction on the probability scale. Public-health
# relevance is additive, so report this too.
f_E4a <- as.formula(paste("hrqol_fairpoor ~ klotho_z * erc +", covars))
E4a <- svyglm(f_E4a, design = subset(des, dom_E4), family = gaussian())
cat("\n=== E4 additive: Klotho x ERC (linear-probability model) ===\n")
print(summary(E4a)$coefficients[grep("klotho_z", rownames(summary(E4a)$coefficients)), ])

## =====================================================================
## SECONDARY / ALTERNATIVE OUTCOMES
## =====================================================================
tidy_beta <- function(fit) round(summary(fit)$coefficients["klotho_z", ], 4)

## Frailty as an OUTCOME (the most promising path: Klotho -> frailty) --
# fi_10 = frailty index per 0.1 units. Negative beta = more Klotho, less frailty.
FI1 <- svyglm(as.formula(paste("fi_10 ~ klotho_z +", covars)),
              design = subset(des, dom_FI), family = gaussian())
FI2 <- svyglm(as.formula(paste("fi_10 ~ klotho_z + egfr + uacr_log +", covars)),
              design = subset(des, dom_FI), family = gaussian())
cat("\n=== Frailty (FI x10) ~ Klotho  [beta = FI*10 change per +1 SD] ===\n")
cat("C only:      "); print(tidy_beta(FI1))
cat("C + renal:   "); print(tidy_beta(FI2))

## Unhealthy days as a continuous HRQoL outcome (2007-2012 only) ------
UD1 <- svyglm(as.formula(paste("unhealthy_days ~ klotho_z +", covars)),
              design = subset(des, dom_UD), family = gaussian())
UD2 <- svyglm(as.formula(paste("unhealthy_days ~ klotho_z + egfr + uacr_log +", covars)),
              design = subset(des, dom_UD), family = gaussian())
cat("\n=== Unhealthy days ~ Klotho  [beta = days per +1 SD] (2007-2012) ===\n")
cat("C only:      "); print(tidy_beta(UD1))
cat("C + renal:   "); print(tidy_beta(UD2))

## Effective sample sizes per analysis --------------------------------
cat("\n=== N per analysis ===\n")
cat("E1/E4 (fair-poor): ", sum(des$variables$dom_E1, na.rm = TRUE),
    " | E2 (+renal): ",     sum(des$variables$dom_E2, na.rm = TRUE),
    " | frailty: ",         sum(des$variables$dom_FI, na.rm = TRUE),
    " | unhealthy days: ",  sum(des$variables$dom_UD, na.rm = TRUE), "\n")

## ---- Save fitted objects ------------------------------------------
saveRDS(list(E1 = E1, E1q = E1q, E2 = E2, E4m = E4m, E4a = E4a,
             E4_noerc = E4_noerc, E4_erc = E4_erc,
             FI1 = FI1, FI2 = FI2, UD1 = UD1, UD2 = UD2,
             klotho_center = mu, klotho_sd = sdw),
        file.path(OUT_DIR, "models_paper1.rds"))
message("\nSaved output/models_paper1.rds")

## ---- Weighted Table 1 by ERC (descriptives) -----------------------
tab_des <- subset(des, dom_E4)
cat("\n=== Weighted means/props by ERC (Table 1 skeleton) ===\n")
print(svyby(~ klotho + fi + age + bmi + pir,     ~ erc, tab_des, svymean, na.rm = TRUE))
print(svyby(~ hrqol_fairpoor + frail + dm + htn, ~ erc, tab_des, svymean, na.rm = TRUE))
print(svytable(~ erc, tab_des))
