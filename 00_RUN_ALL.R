#===============================================================================
# 00_RUN_ALL.R  --  Pathways to Desistance (ICPSR 29961)
# Runs every analysis behind the manuscript, in order, and then checks your
# results against the numbers in the manuscript.
#
# Before you start
#   1. This folder (the one holding 00_RUN_ALL.R and R/) sits inside your ICPSR_29961
#      folder, next to DS0001 ... DS0011. Nothing else needs to be moved.
#   2. Install the packages once (skip any you already have):
#        install.packages(c("lavaan", "dpm", "panelr", "mice", "sandwich",
#                           "lmtest", "fixest", "ggplot2"))
#   3. In RStudio, open this file and click Source (Ctrl+Shift+S).
#
# What runs                                    Writes to (inside ICPSR_29961)
#   01_data_prep.R                 annual panel      panel_desistance.rds, panel_6month.rds
#   03_robustness_checks.R         earlier models    output_checks/
#   04_simulation_lag_artifacts.R  simulations 1-2   output_checks/
#   05_riclpm_main.R               main RI-CLPMs     output_main/
#   06_simulation_common_state.R   simulation 3      output_main/
#   09_multivariate_riclpm.R       several measures  output_main/ (M1-M3)
#   07_figures.R                   Fig. 1 and 2      figures/
#   08_verify_against_reference.R  comparison        verification_report.txt
#
# Your older panel_desistance.rds is copied to panel_desistance_before_2026-10-12.rds
# before 01 replaces it. Your output/ and output_r/ folders are not touched.
#
# Time: roughly 1 to 2 hours, most of it in 05 and 09. Leave RStudio
# open; the console shows each step as it starts and finishes. If one step fails,
# the steps that do not depend on it still run, and the errors are listed at the end.
#===============================================================================

#---- Where things are ----------------------------------------------------------
pkg_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
if (is.na(pkg_dir) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  pkg_dir <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA)
if (is.na(pkg_dir) || !file.exists(file.path(pkg_dir, "R", "05_riclpm_main.R")))
  pkg_dir <- "path/to/ICPSR_29961/pathways-riclpm-replication"   # <- edit if the folder is not found automatically
icpsr_dir <- dirname(pkg_dir)
if (!dir.exists(file.path(icpsr_dir, "DS0001")))
  stop("Could not find DS0001 in ", icpsr_dir, ". Put the package folder inside ICPSR_29961 or edit pkg_dir.")
setwd(icpsr_dir)
Sys.setenv(REF_DIR = file.path(pkg_dir, "reference_outputs"), RUN_PARTS = "EFGHIJKL", REPS = "50")

need <- c("lavaan", "dpm", "panelr", "mice", "sandwich", "lmtest", "fixest", "ggplot2")
missing_pkgs <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs)) stop("Please install first: install.packages(c(\"",
                               paste(missing_pkgs, collapse = "\", \""), "\"))")

cat("Package folder:", pkg_dir, "\nWorking folder:", getwd(), "\n")
cat(R.version.string, "| lavaan", as.character(packageVersion("lavaan")),
    "| dpm", as.character(packageVersion("dpm")), "\n")

#---- Keep a copy of the older panel file -----------------------------------------
if (file.exists("panel_desistance.rds") && !file.exists("panel_desistance_before_2026-10-12.rds")) {
  file.copy("panel_desistance.rds", "panel_desistance_before_2026-10-12.rds")
  cat("Backed up your older panel_desistance.rds\n")
}

#---- Run the steps ----------------------------------------------------------------
errors <- list(); times <- list(); t_all <- Sys.time()
run_step <- function(f, needs = character(0)) {
  if (any(needs %in% names(errors))) {
    errors[[f]] <<- paste("skipped because", paste(intersect(needs, names(errors)), collapse = ", "), "failed")
    cat("\n==== SKIPPED", f, "(", errors[[f]], ")\n"); return(invisible(FALSE))
  }
  cat("\n==================================================================\n")
  cat("==== ", f, " started at ", format(Sys.time(), "%H:%M"), "\n", sep = "")
  cat("==================================================================\n")
  t0 <- Sys.time()
  ok <- tryCatch({ source(file.path(pkg_dir, "R", f), local = new.env()); TRUE },
                 error = function(e) { errors[[f]] <<- conditionMessage(e); FALSE })
  while (sink.number() > 0) sink()                       # close any log left open by an error
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  times[[f]] <<- mins
  cat("==== ", f, if (ok) " finished" else " FAILED", " after ", mins, " minutes\n", sep = "")
  if (!ok) cat("     Error:", errors[[f]], "\n")
  invisible(ok)
}

run_step("01_data_prep.R")
run_step("03_robustness_checks.R",        needs = "01_data_prep.R")
run_step("04_simulation_lag_artifacts.R", needs = "03_robustness_checks.R")
run_step("05_riclpm_main.R",              needs = "01_data_prep.R")
run_step("06_simulation_common_state.R")
run_step("09_multivariate_riclpm.R",      needs = "01_data_prep.R")
run_step("07_figures.R",                  needs = "05_riclpm_main.R")
run_step("08_verify_against_reference.R")

#---- Summary ------------------------------------------------------------------------
cat("\n\n================ SUMMARY ================\n")
for (f in names(times)) cat(sprintf("  %-32s %6.1f min  %s\n", f, times[[f]], if (f %in% names(errors)) "FAILED" else "ok"))
cat(sprintf("  Total: %.0f minutes\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
if (length(errors)) { cat("\nProblems:\n"); for (f in names(errors)) cat("  ", f, ":", errors[[f]], "\n") } else
  cat("\nAll steps finished. Open verification_report.txt in", getwd(), "\n")
cat("If a step failed, the run_log files in output_checks/ and output_main/ show where it stopped.\n")
