# =====================================================================
# 08_sensitivity.R  —  Two light sensitivity runs
#   PART 1: "strict" Frailty Index (no PFQ items) vs full FI, for the
#           Klotho -> frailty association (feeds Table S3).
#   PART 2: restricted-spline shape of Klotho -> frailty, with 95% band
#           and a non-linearity test (Figure 3).
#
#   source("R/08_sensitivity.R")
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({ library(survey); library(dplyr); library(splines) })
options(survey.lonely.psu = "adjust")
FIG_DIR <- file.path(OUT_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
TAB_DIR <- file.path(OUT_DIR, "tables");  dir.create(TAB_DIR, showWarnings = FALSE)

## ---- Load + defensive recompute (as in 04) -------------------------
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

## =====================================================================
## PART 1 — Strict Frailty Index (drop the 11 PFQ physical-function items)
## =====================================================================
yn01 <- function(x) dplyr::case_when(x == 1 ~ 1, x == 2 ~ 0, TRUE ~ NA_real_)
phq_total <- function(df) {
  it <- paste0("DPQ0", c("10","20","30","40","50","60","70","80","90"))
  mm <- as.matrix(df[intersect(it, names(df))]); mm[mm %in% c(7,9)] <- NA
  rowSums(mm, na.rm = FALSE)
}
cond <- c(MCQ010=NA,MCQ160A=NA,MCQ160B=NA,MCQ160C=NA,MCQ160D=NA,MCQ160E=NA,
          MCQ160F=NA,MCQ160G=NA,MCQ160K=NA,MCQ160L=NA,MCQ160M=NA,MCQ160N=NA,
          MCQ160O=NA,MCQ220=NA)
defs <- lapply(names(cond), function(v) if (v %in% names(dat)) yn01(dat[[v]]) else NA_real_)
defm <- do.call(cbind, defs)
defm <- cbind(defm,
              d_dm  = yn01(dat$DIQ010),
              d_htn = yn01(dat$BPQ020),
              d_dep = { p <- phq_total(dat); dplyr::case_when(p>=10~1,p>=5~0.5,p>=0~0,TRUE~NA_real_) },
              d_obe = dplyr::case_when(dat$bmi>=30~1, dat$bmi>=25~0.5, dat$bmi>=0~0, TRUE~NA_real_))
n_it   <- ncol(defm)                                   # 18 items, no PFQ
n_ok   <- rowSums(!is.na(defm))
fi_str <- ifelse(n_ok/n_it >= OPTS$fi_min_items_frac, rowSums(defm, na.rm=TRUE)/n_ok, NA_real_)
dat$fi10_strict <- fi_str * 10
message(sprintf("Strict FI: %d items; computable for %d in-base persons.",
                n_it, sum(!is.na(dat$fi10_strict[m]))))

## ---- Design + compare full vs strict FI (Klotho -> frailty) --------
des <- svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~mecwt,
                 data = dat[!is.na(dat$SDMVPSU) & !is.na(dat$SDMVSTRA), ], nest=TRUE)
C  <- "age + female + race + educ + pir + bmi + dm + htn + smoke + pa_active"
R  <- "egfr + uacr_log"
b1 <- function(fit) round(summary(fit)$coef["klotho_z", c(1,2,4)], 4)

fitC  <- function(y) svyglm(as.formula(paste(y,"~ klotho_z +",C)),
                            design = subset(des, in_base==1 & !is.na(get(y))))
fitCR <- function(y) svyglm(as.formula(paste(y,"~ klotho_z +",C,"+",R)),
                            design = subset(des, in_base==1 & !is.na(get(y))))
s3 <- rbind(
  `Full FI (+C)`         = b1(fitC("fi_10")),
  `Full FI (+C+renal)`   = b1(fitCR("fi_10")),
  `Strict FI (+C)`       = b1(fitC("fi10_strict")),
  `Strict FI (+C+renal)` = b1(fitCR("fi10_strict")))
colnames(s3) <- c("beta","SE","p")
cat("\n=== Table S3 (part): Klotho -> frailty, full vs strict FI ===\n"); print(s3)
write.csv(s3, file.path(TAB_DIR, "TableS3_strictFI.csv"))

## =====================================================================
## PART 2 — Spline shape Klotho -> frailty (Figure 3)
## =====================================================================
des_fi <- subset(des, in_base==1 & !is.na(fi_10))
spl  <- svyglm(as.formula(paste("fi_10 ~ ns(klotho_z, df=3) +", C)), design = des_fi)
quad <- svyglm(as.formula(paste("fi_10 ~ klotho_z + I(klotho_z^2) +", C)), design = des_fi)
p_nl <- summary(quad)$coef["I(klotho_z^2)", "Pr(>|t|)"]
cat(sprintf("\nNon-linearity test (quadratic term) p = %.3f\n", p_nl))

## ---- prediction grid in pg/mL, confounders at reference ------------
d0 <- des_fi$variables
kg  <- seq(quantile(d0$klotho, .02, na.rm=TRUE), quantile(d0$klotho, .98, na.rm=TRUE), length.out = 120)
kz  <- (log(kg) - mu) / sdw
refval <- function(v) { x <- d0[[v]]
  if (is.factor(x)) factor(levels(x)[1], levels = levels(x))
  else if (all(na.omit(x) %in% c(0,1))) 0 else weighted.mean(x, d0$mecwt, na.rm=TRUE) }
Cv <- c("age","female","race","educ","pir","bmi","dm","htn","smoke","pa_active")
newd <- data.frame(klotho_z = kz)
for (v in Cv) newd[[v]] <- rep(refval(v), length(kz))

pr  <- predict(spl, newdata = newd, se.fit = TRUE)
est <- as.numeric(pr) / 10                             # back to FI 0-1 scale
se  <- as.numeric(SE(pr)) / 10
lo  <- est - 1.96*se; hi <- est + 1.96*se

## ---- draw Figure 3 -------------------------------------------------
draw_fig3 <- function() {
  teal <- "#0F6E56"
  par(mar = c(4.2, 4.4, 2.2, 1))
  plot(kg, est, type="n", ylim = range(c(lo,hi)),
       xlab = "Serum α-Klotho (pg/mL)", ylab = "Adjusted frailty index",
       main = sprintf("Klotho–frailty dose-response (non-linearity p = %.2f)", p_nl),
       cex.main = 0.95, font.main = 1)
  polygon(c(kg, rev(kg)), c(lo, rev(hi)), col = paste0(teal,"22"), border = NA)
  lines(kg, est, col = teal, lwd = 2.5)
  rug(sample(d0$klotho[!is.na(d0$klotho)], min(1500, sum(!is.na(d0$klotho)))),
      col = "#88878055", ticksize = 0.02)
  abline(h = pretty(range(c(lo,hi))), col = "#00000010")
}
pdf(file.path(FIG_DIR, "Figure3_klotho_frailty_spline.pdf"), width = 6.6, height = 4.4)
draw_fig3(); dev.off()
png(file.path(FIG_DIR, "Figure3_klotho_frailty_spline.png"), width = 6.6, height = 4.4,
    units = "in", res = 300)
draw_fig3(); dev.off()
message("Figure 3 saved to ", FIG_DIR, " (PDF + PNG).")
message("Table S3 (strict FI) saved to ", TAB_DIR, "/TableS3_strictFI.csv")
