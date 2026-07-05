# =====================================================================
# 05_mediation.R  —  E3: mediation of the Klotho -> HRQoL association
#                    by frailty, under the complex survey design.
# ---------------------------------------------------------------------
# Regression-based natural (in)direct effects with an exposure-mediator
# interaction (VanderWeele), for a CONTINUOUS outcome (linear-linear),
# so the closed-form NDE/NIE are exact. CIs come from a survey
# subbootstrap that resamples PSUs within strata (design-consistent).
#
#   A = klotho_z   (exposure, +1 SD log-Klotho)
#   M = fi_10      (mediator, frailty per 0.1 FI)
#   Y = unhealthy_days  (continuous HRQoL; see note on binary outcome)
#   confounders for all three links: C + R  (age..pa_active, egfr, uacr_log)
#
# WARNINGS baked into the design (from the DAG discussion):
#   * Do NOT enter inflammation (CRP) as a plain covariate here: it is a
#     mediator-outcome confounder AFFECTED BY the exposure. If you want to
#     account for it, use interventional (in)direct effects (g-formula),
#     not this regression form. Left out by default.
#   * The FI must exclude HRQoL-overlapping items (handled in 03).
#   * A and M are measured contemporaneously -> the exposure->mediator
#     link is cross-sectional. Report NDE/NIE as a statistical
#     decomposition, not proof of a temporal mechanism.
#
# Binary outcome (fair/poor health): linear-linear closed forms do not
# apply. Either (a) use unhealthy_days as here, or (b) use CMAverse::cmest
# / regmedint with a bootstrap over the same replicate weights.
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({ library(dplyr); library(survey) })
options(survey.lonely.psu = "adjust")
set.seed(OPTS$seed)

dat <- readRDS(file.path(OUT_DIR, "analytic_with_fi.rds"))

# recompute klotho_z consistently with 04
mu  <- weighted.mean(dat$klotho_log[dat$in_base == 1], dat$mecwt[dat$in_base == 1], na.rm = TRUE)
sdw <- sqrt(Hmisc::wtd.var(dat$klotho_log[dat$in_base == 1],
                           weights = dat$mecwt[dat$in_base == 1], na.rm = TRUE))
dat$klotho_z <- (dat$klotho_log - mu) / sdw

C  <- c("age","female","race","educ","pir","bmi","dm","htn","smoke","pa_active")
R  <- c("egfr","uacr_log")
conf_rhs <- paste(c(C, R), collapse = " + ")

# Mediation analytic domain: complete on A, M, Y, C, R
dat$dom_med <- as.integer(
  dat$in_base == 1 &
    !is.na(dat$klotho_z) & !is.na(dat$fi_10) & !is.na(dat$unhealthy_days) &
    rowSums(is.na(dat[c(C, R)])) == 0
)
dat$dom_med[is.na(dat$dom_med)] <- 0L
message("Mediation domain n = ", sum(dat$dom_med == 1),
        " (needs the Healthy-Days module; 0 => use fair/poor + CMAverse).")

## ---- Point estimate of NDE / NIE from a fitted survey design ------
# Linear mediator:  M = b0 + bA*A + bC*C
# Linear outcome:   Y = t0 + tA*A + tM*M + tAM*A*M + tC*C
# Natural direct (change A: a* -> a) and indirect effects (VanderWeele,
# continuous M and Y), evaluated at mediator confounders' means:
#   NDE = (tA + tAM*(b0 + bA*a* + bC*Cbar)) * (a - a*)
#   NIE = (tM*bA + tAM*bA*a) * (a - a*)
# Here a = 1, a* = 0 (a +1 SD contrast in log-Klotho).
mediation_effects <- function(design) {
  med_fit <- svyglm(as.formula(paste("fi_10 ~ klotho_z +", conf_rhs)),
                    design = design)
  out_fit <- svyglm(as.formula(paste("unhealthy_days ~ klotho_z * fi_10 +", conf_rhs)),
                    design = design)
  bA <- coef(med_fit)["klotho_z"]
  b0 <- coef(med_fit)["(Intercept)"]
  # mean of mediator confounders' linear-predictor contribution, obtained
  # from fitted values (Gaussian => fitted == linear predictor):
  A_mean <- mean(med_fit$model$klotho_z)
  Cbar_contrib <- mean(as.numeric(fitted(med_fit))) - (b0 + bA * A_mean)
  tA  <- coef(out_fit)["klotho_z"]
  tM  <- coef(out_fit)["fi_10"]
  tAM <- coef(out_fit)["klotho_z:fi_10"]
  a <- 1; astar <- 0
  m_at_astar <- b0 + bA * astar + Cbar_contrib
  NDE <- (tA + tAM * m_at_astar) * (a - astar)
  NIE <- (tM * bA + tAM * bA * a) * (a - astar)
  TE  <- NDE + NIE
  PM  <- NIE / TE
  c(NDE = unname(NDE), NIE = unname(NIE), TE = unname(TE), prop_mediated = unname(PM))
}

## ---- Survey subbootstrap for CIs ----------------------------------
dat$mecwt[is.na(dat$mecwt)] <- 0                       # survey >= 4.4: no NA weights
des_data <- dat[!is.na(dat$SDMVPSU) & !is.na(dat$SDMVSTRA), ]
des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~mecwt,
                 data = des_data, nest = TRUE)
des <- update(des, dom_med = dom_med == 1)
des_med <- subset(des, dom_med)

point <- mediation_effects(des_med)
cat("\n=== E3 point estimates (units: unhealthy days per +1 SD log-Klotho) ===\n")
print(round(point, 4))

cat("\nBootstrapping CIs (", OPTS$boot_reps, " replicates)...\n", sep = "")
rep_des <- as.svrepdesign(des_med, type = "subbootstrap", replicates = OPTS$boot_reps)
# refit on each replicate weight column
W <- weights(rep_des, type = "analysis")
boot <- t(apply(W, 2, function(wcol) {
  d_i <- rep_des$variables
  d_i$.w <- wcol
  di <- svydesign(ids = ~1, weights = ~.w, data = d_i)  # replicate: weights carry design
  tryCatch(mediation_effects(di), error = function(e) rep(NA_real_, 4))
}))
ci <- apply(boot, 2, quantile, probs = c(.025, .975), na.rm = TRUE)
cat("\n=== E3 with 95% subbootstrap CIs ===\n")
res <- data.frame(estimate = point, lcl = ci[1, ], ucl = ci[2, ])
print(round(res, 4))

saveRDS(list(point = point, boot = boot, ci = ci, res = res),
        file.path(OUT_DIR, "mediation_E3.rds"))
message("\nSaved output/mediation_E3.rds")

cat("\nReporting note: NIE = the part of Klotho's association with worse\n",
    "HRQoL that runs through higher frailty; NDE = the remainder. Frame as\n",
    "a cross-sectional decomposition. For the binary fair/poor outcome,\n",
    "re-run with CMAverse::cmest(model='rb', ...) over the same replicate\n",
    "weights, or use interventional effects if CRP/FGF23 are added.\n")

