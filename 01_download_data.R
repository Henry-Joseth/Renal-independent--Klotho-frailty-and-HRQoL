# =====================================================================
# 01_download_data.R  —  Download the NHANES .XPT files (via nhanesA)
# ---------------------------------------------------------------------
# Direct download.file() from wwwn.cdc.gov returns HTML error pages on
# many networks (and during CDC outages), so this script uses the
# nhanesA package, which fetches reliably. It pulls EVERY component x
# cycle listed in 00_config.R and writes each as data/raw/<TABLE>.XPT,
# exactly what 02_derive_variables.R expects.
#
#   source("R/01_download_data.R")
#
# Already-downloaded, valid files are skipped, so re-running is cheap and
# safe (it will NOT re-download or clobber good files).
# =====================================================================

source("R/00_config.R")
if (!requireNamespace("nhanesA", quietly = TRUE))
  install.packages("nhanesA", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({ library(nhanesA); library(haven) })

# A valid SAS-transport file begins with "HEADER RECORD" and is > 1 KB.
valid_xpt <- function(p) {
  if (!file.exists(p) || file.info(p)$size < 1024) return(FALSE)
  con <- file(p, "rb"); on.exit(close(con))
  grepl("HEADER RECORD",
        tryCatch(readChar(con, 80, useBytes = TRUE), error = function(e) ""),
        fixed = TRUE)
}

# Fallback for the Klotho surplus file if nhanesA can't fetch it: the CDC
# "Public" data path (occasionally works when nhanesA doesn't).
sskl_direct <- function(suf, dest) {
  yy  <- substr(CYCLES$years_folder[CYCLES$suffix == suf], 1, 4)
  for (u in c(
    sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/SSKL_%s.xpt", yy, suf),
    sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/SSKL_%s.XPT", yy, suf))) {
    ok <- tryCatch({ download.file(u, dest, mode = "wb", method = "libcurl", quiet = TRUE)
                     valid_xpt(dest) }, error = function(e) FALSE)
    if (isTRUE(ok)) return(TRUE)
    if (file.exists(dest)) unlink(dest)
  }
  FALSE
}

# Full list of tables: COMPONENT_SUFFIX for every cycle.
grid <- do.call(rbind, lapply(seq_len(nrow(CYCLES)), function(i)
  data.frame(comp = names(COMPONENTS), suf = CYCLES$suffix[i],
             stringsAsFactors = FALSE)))
grid$name <- paste0(grid$comp, "_", grid$suf)

options(timeout = max(600, getOption("timeout")))
message("Fetching ", nrow(grid), " NHANES tables via nhanesA into ", RAW_DIR, " ...")

fetch_one <- function(comp, suf, name) {
  dest <- file.path(RAW_DIR, paste0(name, ".XPT"))
  if (valid_xpt(dest)) return("exists")
  # nhanesA returns RAW coded values (translated = FALSE) with correct
  # upper-case variable names; NULL/empty if the table is absent that cycle.
  d <- tryCatch(nhanesA::nhanes(name, translated = FALSE), error = function(e) NULL)
  if (!is.null(d) && nrow(d)) {
    haven::write_xpt(d, dest)
    if (valid_xpt(dest)) { message("  ok:   ", name, "  (", nrow(d), "x", ncol(d), ")"); return("ok") }
  }
  if (comp == "SSKL" && sskl_direct(suf, dest)) { message("  ok:   ", name, " (direct)"); return("ok") }
  message("  n/a:  ", name); "fail"
}

res <- mapply(fetch_one, grid$comp, grid$suf, grid$name)
n_ok <- sum(res %in% c("ok", "exists"))
message(sprintf("Download step done: %d/%d tables available. Files in data/raw: %d",
                n_ok, length(res), length(list.files(RAW_DIR))))
if (any(res == "fail"))
  message("Not downloaded (may not exist for that cycle — the coverage ",
          "check in 02 will confirm if it matters): ",
          paste(grid$name[res == "fail"], collapse = ", "))

