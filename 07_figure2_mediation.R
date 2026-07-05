# =====================================================================
# 07_figure2_mediation.R  —  Figure 2: mediation path diagram
# ---------------------------------------------------------------------
# Draws alpha-Klotho -> Frailty -> HRQoL with the a-path, b-path and the
# NDE/NIE/total-effect decomposition, exported as PDF (vector) and PNG
# (300 dpi) into output/figures/. Coefficients are computed live so the
# figure always matches the analysis.
#
#   source("R/07_figure2_mediation.R")
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({ library(survey); library(dplyr) })
options(survey.lonely.psu = "adjust")
FIG_DIR <- file.path(OUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)

## ---- Rebuild design + refit the mediation models -------------------
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

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~mecwt,
                 data = dat[!is.na(dat$SDMVPSU) & !is.na(dat$SDMVSTRA), ], nest = TRUE)
des_med <- subset(des, in_base == 1 & !is.na(fi_10) & !is.na(unhealthy_days) &
                    !is.na(egfr) & !is.na(uacr_log))
conf <- "age + female + race + educ + pir + bmi + dm + htn + smoke + pa_active + egfr + uacr_log"

med_fit <- svyglm(as.formula(paste("fi_10 ~ klotho_z +", conf)), design = des_med)
out_fit <- svyglm(as.formula(paste("unhealthy_days ~ klotho_z + fi_10 +", conf)), design = des_med)
a_path <- coef(med_fit)["klotho_z"]                    # Klotho -> frailty
b_path <- coef(out_fit)["fi_10"]                       # frailty -> unhealthy days
cprime <- coef(out_fit)["klotho_z"]                    # direct

E3 <- readRDS(file.path(OUT_DIR, "mediation_E3.rds"))$res
fmt  <- function(x) sprintf("%.2f", x)
fmtci<- function(r) sprintf("%.2f (%.2f, %.2f)", r["estimate"], r["lcl"], r["ucl"])
lab_a   <- sprintf("a = %.3f per SD", a_path)
lab_b <- sprintf("b = %.2f per 0.1 FI", b_path)
lab_nie <- paste0("Indirect (NIE): ", fmtci(E3["NIE", ]))
lab_nde <- paste0("Direct (NDE): ", fmtci(E3["NDE", ]), "  ns")
lab_te  <- paste0("Total effect: ", fmtci(E3["TE", ]))

## ---- Draw ----------------------------------------------------------
draw_fig <- function() {
  teal <- "#0F6E56"; gray <- "#888780"; purple <- "#534AB7"; coral <- "#993C1D"
  par(mar = c(1, 1, 1, 1))
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "", asp = 1)
  boxf <- function(xc, yc, w, h, title, sub, col) {
    rect(xc-w/2, yc-h/2, xc+w/2, yc+h/2, border = col, lwd = 2,
         col = paste0(col, "18"), xpd = NA)
    text(xc, yc+2.5, title, font = 2, col = col, cex = 0.95)
    text(xc, yc-3.5, sub, col = col, cex = 0.7)
  }
  # boxes: A left, M top-center, Y right
  Ax <- 16; Ay <- 30; Mx <- 50; My <- 82; Yx <- 84; Yy <- 30
  # arrows (indirect route, teal)
  arrows(Ax+6, Ay+9, Mx-14, My-9, col = teal, lwd = 3, length = 0.12)
  arrows(Mx+14, My-9, Yx-6, Yy+9, col = teal, lwd = 3, length = 0.12)
  # direct route (gray dashed)
  arrows(Ax+13, Ay-3, Yx-13, Yy-3, col = gray, lwd = 2.5, lty = 2, length = 0.12)
  # boxes on top
  boxf(Ax, Ay, 26, 16, "alpha-Klotho", "exposure (+1 SD log)", purple)
  boxf(Mx, My, 26, 16, "Frailty index", "mediator", teal)
  boxf(Yx, Yy, 26, 16, "Unhealthy days", "HRQoL outcome", coral)
  # path labels
  text(28, 62, lab_a, col = teal, cex = 0.72, font = 2)
  text(72, 62, lab_b, col = teal, cex = 0.72, font = 2)
  text(50, 40, lab_nie, col = teal, cex = 0.82, font = 2)
  text(50, 22, lab_nde, col = gray, cex = 0.8)
  text(50, 8,  lab_te, col = "black", cex = 0.85, font = 2)
}

pdf(file.path(FIG_DIR, "Figure2_mediation.pdf"), width = 7.2, height = 4.2)
draw_fig(); dev.off()
png(file.path(FIG_DIR, "Figure2_mediation.png"), width = 7.2, height = 4.2,
    units = "in", res = 300)
draw_fig(); dev.off()

message("Figure 2 saved to ", FIG_DIR, " (PDF + PNG).")
message(sprintf("  a-path=%.3f  b-path=%.3f  direct(c')=%.3f", a_path, b_path, cprime))
message("  ", lab_nie, " | ", lab_nde, " | ", lab_te)
