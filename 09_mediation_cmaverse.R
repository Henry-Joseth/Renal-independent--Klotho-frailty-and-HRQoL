# =====================================================================
# 09_mediation_cmaverse.R  —  Sensitivity mediation with the BINARY
#   HRQoL outcome (fair/poor health, all 5 cycles), complementing the
#   continuous-outcome mediation in 05.
#
# METHOD NOTE (pragmatic path): CMAverse::cmest does not natively ingest
# a complex survey design (strata + PSU). The PRIMARY result here is the
# regression-based natural-effects estimator of VanderWeele — the same
# estimator cmest(model="rb") implements — but fit with survey weights
# and a design-based subbootstrap that DOES respect strata/PSU. CMAverse
# is then run as a cross-check (its bootstrap resamples rows and ignores
# the clustering, so treat its CIs as approximate).
#
#   source("R/09_mediation_cmaverse.R")
#   NBOOT = 200 to test; raise to 1000 for the final run.
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({ library(survey); library(dplyr) })
options(survey.lonely.psu = "adjust")
set.seed(OPTS$seed)
NBOOT <- 200                                            # <- raise to 1000 for final
TAB_DIR <- file.path(OUT_DIR, "tables"); dir.create(TAB_DIR, showWarnings = FALSE)

## ---- Load + defensive recompute ------------------------------------
dat <- readRDS(file.path(OUT_DIR, "analytic_with_fi.rds"))
dat$hrqol_fairpoor <- ifelse(is.na(dat$genhealth), NA_integer_,
                             as.integer(dat$genhealth %in% 4:5))
core <- c("age","female","race","educ","pir","bmi","dm","htn","smoke","pa_active")
dat$in_base <- as.integer(
  !is.na(dat$mecwt) & dat$mecwt > 0 & dat$age >= OPTS$age_min & dat$age <= OPTS$age_max &
    !is.na(dat$klotho) & !is.na(dat$hrqol_fairpoor) & rowSums(is.na(dat[core])) == 0)
dat$in_base[is.na(dat$in_base)] <- 0L
m <- dat$in_base == 1
mu  <- weighted.mean(dat$klotho_log[m], dat$mecwt[m], na.rm = TRUE)
sdw <- sqrt(Hmisc::wtd.var(dat$klotho_log[m], weights = dat$mecwt[m], na.rm = TRUE))
dat$klotho_z <- (dat$klotho_log - mu) / sdw
dat$mecwt[is.na(dat$mecwt)] <- 0

C <- c("age","female","race","educ","pir","bmi","dm","htn","smoke","pa_active")
R <- c("egfr","uacr_log")
conf <- paste(c(C, R), collapse = " + ")

dat$dom_med <- as.integer(dat$in_base == 1 & !is.na(dat$fi_10) & !is.na(dat$hrqol_fairpoor) &
                            rowSums(is.na(dat[c(C, R)])) == 0)
dat$dom_med[is.na(dat$dom_med)] <- 0L
message("Binary-outcome mediation domain n = ", sum(dat$dom_med == 1))

## ---- Weighted VanderWeele natural effects (OR scale) ---------------
# Mediator (linear): M = b0 + b1 A + b'C ; residual variance sig2
# Outcome (logistic): logit P(Y) = t0 + t1 A + t2 M + t3 A*M + t'C
# log OR^NDE = t1 + t3(E[M|a*] + t2 sig2) + 0.5 t3^2 sig2      (a=1, a*=0)
# log OR^NIE = t2 b1 + t3 b1
mediation_or <- function(des) {
  med <- svyglm(as.formula(paste("fi_10 ~ klotho_z +", conf)), design = des)
  out <- svyglm(as.formula(paste("hrqol_fairpoor ~ klotho_z * fi_10 +", conf)),
                design = des, family = quasibinomial())
  b1   <- coef(med)["klotho_z"]
  Abar <- mean(med$model$klotho_z)
  mA0  <- mean(as.numeric(fitted(med))) - b1 * Abar          # E[M | A=0], averaged over C
  sig2 <- summary(med)$dispersion[1]
  t1 <- coef(out)["klotho_z"]; t2 <- coef(out)["fi_10"]; t3 <- coef(out)["klotho_z:fi_10"]
  logNDE <- t1 + t3 * (mA0 + t2 * sig2) + 0.5 * t3^2 * sig2
  logNIE <- t2 * b1 + t3 * b1
  orNDE <- exp(logNDE); orNIE <- exp(logNIE); orTE <- orNDE * orNIE
  pm <- (orNDE * (orNIE - 1)) / (orNDE * orNIE - 1)
  c(OR_NDE = unname(orNDE), OR_NIE = unname(orNIE),
    OR_TE = unname(orTE), prop_mediated = unname(pm))
}

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~mecwt,
                 data = dat[!is.na(dat$SDMVPSU) & !is.na(dat$SDMVSTRA), ], nest = TRUE)
des <- update(des, dom_med = dom_med == 1)
des_med <- subset(des, dom_med)

point <- mediation_or(des_med)
cat("\n=== Binary-outcome mediation, weighted point estimates (OR scale) ===\n")
print(round(point, 4))

cat("\nDesign-based subbootstrap (", NBOOT, " reps)...\n", sep = "")
rep_des <- as.svrepdesign(des_med, type = "subbootstrap", replicates = NBOOT)
W <- weights(rep_des, type = "analysis")
boot <- t(apply(W, 2, function(wc) {
  d_i <- rep_des$variables; d_i$.w <- wc
  di <- svydesign(ids = ~1, weights = ~.w, data = d_i)
  tryCatch(mediation_or(di), error = function(e) rep(NA_real_, 4))
}))
ci  <- apply(boot, 2, quantile, probs = c(.025, .975), na.rm = TRUE)
res <- data.frame(estimate = round(point, 3),
                  lcl = round(ci[1, ], 3), ucl = round(ci[2, ], 3))
cat("\n=== Binary-outcome mediation with 95% subbootstrap CIs (OR scale) ===\n")
print(res)
write.csv(res, file.path(TAB_DIR, "TableS3_mediation_binary.csv"))

## ---- CMAverse cross-check (approximate: ignores strata/PSU) --------
if (requireNamespace("CMAverse", quietly = TRUE)) {
  cm <- tryCatch({
    d_cm <- des_med$variables
    args <- list(data = d_cm, model = "rb", outcome = "hrqol_fairpoor",
                 exposure = "klotho_z", mediator = "fi_10", EMint = TRUE,
                 basec = c(C, R), yreg = "logistic", mreg = list("linear"),
                 astar = 0, a = 1, mval = list(median(d_cm$fi_10, na.rm = TRUE)),
                 estimation = "paramfunc", inference = "bootstrap", nboot = NBOOT)
    if ("weights" %in% names(formals(CMAverse::cmest))) args$weights <- d_cm$mecwt
    do.call(CMAverse::cmest, args)
  }, error = function(e) { message("CMAverse cross-check skipped: ", conditionMessage(e)); NULL })
  if (!is.null(cm)) { cat("\n--- CMAverse cmest(model='rb') cross-check ---\n"); print(summary(cm)) }
} else {
  message("CMAverse not installed; skipping cross-check. install.packages('CMAverse') to enable.")
}

message("\nSaved output/tables/TableS3_mediation_binary.csv")
cat("\nNote: OR^NIE > or < 1 indicates the frailty-mediated part of the Klotho-HRQoL\n",
    "odds ratio. Proportion mediated is on the OR scale and is interpretable only\n",
    "with the usual caution when the outcome is not rare (~18-22% here).\n")
