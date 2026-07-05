# =====================================================================
# 03_frailty_index.R  —  Rockwood-style deficit-accumulation Frailty
#                        Index (FI), harmonized across 2007-2016.
# ---------------------------------------------------------------------
# Design decisions baked in (per the DAG discussion):
#   * The FI EXCLUDES self-rated general health (HSD010) and the
#     activity-limitation Healthy-Days item, because those ARE the HRQoL
#     outcome. Mixing them would make the mediation partly tautological.
#   * Kidney/renal items are EXCLUDED from the FI so the mediator does
#     not contaminate the separate renal (R) node.
#   * PFQ physical-function items are included by default but can be
#     dropped (OPTS$fi_include_pfq = FALSE) for a "strict" sensitivity FI.
#
# Rules (Searle et al. 2008):
#   - each deficit coded 0 (absent) .. 1 (present), graded allowed;
#   - a person needs >= fi_min_items_frac of items non-missing;
#   - FI = mean of available deficits (0-1). frail = FI >= 0.25.
# =====================================================================

source("R/00_config.R")
suppressPackageStartupMessages({ library(dplyr); library(purrr) })

dat <- readRDS(file.path(OUT_DIR, "analytic_raw.rds"))

## ---- deficit recoders ---------------------------------------------
yn01 <- function(x) dplyr::case_when(x == 1 ~ 1, x == 2 ~ 0, TRUE ~ NA_real_)

# PFQ061*: 1 no difficulty, 2 some, 3 much, 4 unable, 7/9 NA
pf_grade <- function(x) dplyr::case_when(
  x == 1 ~ 0, x == 2 ~ 0.5, x %in% c(3,4) ~ 1, TRUE ~ NA_real_)

# PHQ-9 total (0-27); graded deficit 0 / 0.5 / 1
phq_total <- function(df) {
  items <- paste0("DPQ0", c("10","20","30","40","50","60","70","80","90"))
  m <- as.matrix(df[intersect(items, names(df))])
  m[m %in% c(7,9)] <- NA
  rowSums(m, na.rm = FALSE)
}

## ---- assemble deficit matrix --------------------------------------
d <- dat

# 1) Self-reported chronic conditions (MCQ + DIQ + BPQ) --------------
cond_yn <- list(
  asthma      = "MCQ010",  arthritis = "MCQ160A", chf   = "MCQ160B",
  chd         = "MCQ160C", angina    = "MCQ160D", mi    = "MCQ160E",
  stroke      = "MCQ160F", emphysema = "MCQ160G", chronbronch = "MCQ160K",
  liver       = "MCQ160L", thyroid   = "MCQ160M", gout  = "MCQ160N",
  copd        = "MCQ160O", cancer    = "MCQ220"
)
def <- imap_dfc(cond_yn, function(var, nm) {
  tibble::tibble(!!paste0("d_", nm) := if (var %in% names(d)) yn01(d[[var]]) else NA_real_)
})
def$d_diabetes <- yn01(d$DIQ010)
def$d_htn      <- yn01(d$BPQ020)

# 2) Depression (PHQ-9) ---------------------------------------------
phq <- phq_total(d)
def$d_depression <- dplyr::case_when(phq >= 10 ~ 1, phq >= 5 ~ 0.5,
                                     phq >= 0 ~ 0, TRUE ~ NA_real_)

# 3) Anthropometry ---------------------------------------------------
def$d_obese <- dplyr::case_when(d$bmi >= 30 ~ 1, d$bmi >= 25 ~ 0.5,
                                d$bmi >= 0  ~ 0, TRUE ~ NA_real_)

# 4) Physical function (optional; PFQ061*) ---------------------------
pfq_vars <- grep("^PFQ061", names(d), value = TRUE)
if (OPTS$fi_include_pfq && length(pfq_vars) > 0) {
  pf <- as.data.frame(lapply(d[pfq_vars], pf_grade))
  names(pf) <- paste0("d_", tolower(pfq_vars))
  def <- dplyr::bind_cols(def, pf)
}

## ---- NOTE on deliberate exclusions --------------------------------
# NOT included as deficits (by design):
#   HSD010  (self-rated health)      -> HRQoL outcome
#   HSQ490  (activity-limited days)  -> HRQoL outcome
#   any KIQ renal items              -> separate R node
# Do not add them here.

## ---- compute FI with missingness rule -----------------------------
def_mat  <- as.matrix(def)
n_items  <- ncol(def_mat)
n_ok     <- rowSums(!is.na(def_mat))
frac_ok  <- n_ok / n_items
fi_sum   <- rowSums(def_mat, na.rm = TRUE)
fi       <- ifelse(frac_ok >= OPTS$fi_min_items_frac, fi_sum / n_ok, NA_real_)

dat$fi          <- fi
dat$fi_n_items  <- n_items
dat$fi_n_ok     <- n_ok
dat$frail       <- as.integer(fi >= 0.25)
# scaled version (per 0.1 FI) for interpretable model coefficients
dat$fi_10       <- fi * 10

message(sprintf("Frailty Index: %d candidate deficits (PFQ %s).",
                n_items, ifelse(OPTS$fi_include_pfq, "included", "excluded")))
message(sprintf("FI computable for %d of %d in-base persons; median FI = %.3f.",
                sum(!is.na(dat$fi[dat$in_base == 1])),
                sum(dat$in_base == 1),
                median(dat$fi[dat$in_base == 1], na.rm = TRUE)))

saveRDS(dat, file.path(OUT_DIR, "analytic_with_fi.rds"))
message("Saved output/analytic_with_fi.rds")
