#===============================================================================
# 09_multivariate_riclpm.R  --  Pathways to Desistance (ICPSR 29961)
# Multivariate RI-CLPMs: several correlates and offending in one model
#
# The main analysis (05) pairs one correlate with offending at a time, so a
# correlate's cross-lag can carry the effect of another correlate that moves with
# it. Here each model holds a whole set of correlates and offending together:
# every within-person deviation at age a-1 predicts every deviation at age a.
# A forward effect is then the effect of last year's measure on this year's
# offending, holding last year's offending and last year's other measures constant.
#
#   M1 Cognitive shifts         expectations, moral disengagement, social costs,
#                               certainty, personal rewards, social rewards, identity
#   M2 Psychosocial maturation  impulse control, suppression of aggression,
#                               consideration of others, future orientation,
#                               resistance to peer influence (the overall PSMI score
#                               is left out because it contains the identity subscale)
#   M3 Social bonds             employment, school, living with a partner, parent
#   M4 Across domains           the measures with a prospective effect in the
#                               pairwise models: moral disengagement, certainty,
#                               personal rewards, social rewards, impulse control,
#                               suppression of aggression, peer antisocial influence
#
# Specification (same as 05 except where noted)
#   - ages 16-23, interviews placed at age; random intercept for every variable
#   - autoregressive paths free at every age
#   - cross-lags ONE value for all ages (BIC favored this over age bands in all 18
#     pairwise models, J1); comparable to the "one cross-lag for all ages" columns of
#     J1_age_bands_vs_invariant.csv
#   - same-age residual covariances free at every age
#   - scales and offending per pooled SD; bonds on their 0-1 scale
#   - FIML, robust (MLR) standard errors
#
# Requires: panel_desistance.rds from 01_data_prep.R; output_main/J1 from 05 (for
#           the side-by-side comparison; the models run without it)
# Packages: lavaan
# Output  : output_main/M1_multivariate_cross_lags.csv  forward and reverse, with pairwise
#           output_main/M2_multivariate_fit.csv  convergence, size, log-likelihood (no CFI/RMSEA)
#           output_main/M3_multivariate_all_lagged_paths.csv  every lagged path
#           output_main/run_log_multivariate.txt
# Run     : from the ICPSR_29961 folder, source("09_multivariate_riclpm.R")  (about 30 to 60 minutes)
#===============================================================================

suppressPackageStartupMessages(library(lavaan))
data_dir <- "."
out_dir  <- file.path(data_dir, "output_main")
dir.create(out_dir, showWarnings = FALSE)

panel <- readRDS(file.path(data_dir, "panel_desistance.rds"))
if (is.null(panel$YRS)) panel$YRS <- panel$year
panel$AGE <- panel$BL_AGE + panel$YRS
AGES <- 16:23
ids  <- sort(unique(panel$CASEID))
bonds <- c("WORK_PROP", "SCHOOL_PROP", "PARTNER", "PARENT")
lab_map <- c(WORK_PROP = "Employment", SCHOOL_PROP = "School", PARTNER = "Living with a partner",
             PARENT = "Parent", EXPECT = "Expectations for the future", MORDIS = "Moral disengagement",
             SOCCST = "Social costs of punishment", PUNYOU = "Certainty of punishment",
             PERRWD = "Personal rewards of crime", SOCRWD = "Social rewards of crime",
             PSMIID = "Identity (PSMI)", IMPULS = "Impulse control", SUPAGG = "Suppression of aggression",
             CONSID = "Consideration of others", FUTURE = "Future orientation",
             PEERPR = "Resistance to peer influence", PSMI = "Psychosocial maturity (PSMI)",
             PRINFL = "Peer antisocial influence", OFF = "Offending")

models <- list(
  M1_cognitive  = c("EXPECT", "MORDIS", "SOCCST", "PUNYOU", "PERRWD", "SOCRWD", "PSMIID"),
  M2_maturation = c("IMPULS", "SUPAGG", "CONSID", "FUTURE", "PEERPR"),
  M3_bonds      = c("WORK_PROP", "SCHOOL_PROP", "PARTNER", "PARENT"),
  M4_across     = c("MORDIS", "PUNYOU", "PERRWD", "SOCRWD", "IMPULS", "SUPAGG", "PRINFL"))
only <- Sys.getenv("MV_MODELS", "")                  # optional: e.g. MV_MODELS=M3_bonds
if (nzchar(only)) models <- models[intersect(names(models), strsplit(only, ",")[[1]])]

zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
panel$OFF <- zs(panel$OFFVAR)

make_wide <- function(vars) {
  w <- data.frame(CASEID = ids)
  for (v in c(vars, "OFF")) {
    x <- if (v %in% bonds || v == "OFF") panel[[v]] else zs(panel[[v]])
    for (a in AGES) {
      r <- which(panel$AGE == a); m <- match(ids, panel$CASEID[r])
      w[[paste0(v, "_", a)]] <- x[r][m]
    }
  }
  w[rowSums(!is.na(w[, paste0("OFF_", AGES)])) > 0, ]   # same people as 05: offending seen at least once
}

mv_syntax <- function(vars) {
  V <- c(vars, "OFF"); s <- character(0)
  for (v in V) {
    s <- c(s, sprintf("RI_%s =~ %s", v, paste0("1*", v, "_", AGES, collapse = " + ")))
    for (a in AGES) s <- c(s, sprintf("w_%s_%d =~ 1*%s_%d", v, a, v, a), sprintf("%s_%d ~~ 0*%s_%d", v, a, v, a))
  }
  for (i in 2:length(AGES)) {
    a <- AGES[i]; p <- AGES[i - 1]
    for (v in V) {
      rhs <- c(sprintf("ar_%s_%d*w_%s_%d", v, a, v, p),
               sapply(setdiff(V, v), function(u) sprintf("c_%s__%s*w_%s_%d", u, v, u, p)))   # c_FROM__TO
      s <- c(s, sprintf("w_%s_%d ~ %s", v, a, paste(rhs, collapse = " + ")))
    }
  }
  for (a in AGES) {                                    # same-age (residual) covariances, free at each age
    s <- c(s, sprintf("w_%s_%d ~~ w_%s_%d", V, a, V, a))
    pr <- combn(V, 2)
    s <- c(s, sprintf("w_%s_%d ~~ w_%s_%d", pr[1, ], a, pr[2, ], a))
  }
  RI <- paste0("RI_", V)
  s <- c(s, sprintf("%s ~~ %s", RI, RI))
  pr <- combn(RI, 2); s <- c(s, sprintf("%s ~~ %s", pr[1, ], pr[2, ]))
  for (r in RI) s <- c(s, sprintf("%s ~~ 0*w_%s_%d", r, V, AGES[1]))
  paste(s, collapse = "\n")
}

sink(file.path(out_dir, "run_log_multivariate.txt"), split = TRUE)
cat("Multivariate RI-CLPMs --", format(Sys.time()), "\n")
fwd <- list(); fit_tab <- list(); all_paths <- list()
for (mn in names(models)) {
  vars <- models[[mn]]
  w <- make_wide(vars)
  cat("\n=====", mn, ":", paste(vars, collapse = ", "), "| N =", nrow(w), "=====\n")
  t0 <- Sys.time()
  # h1 and baseline models are skipped: with 40 to 64 observed variables they take far
  # longer than the model itself, and they are needed only for CFI and RMSEA
  fit <- tryCatch(lavaan(mv_syntax(vars), data = w, missing = "ml", estimator = "MLR",
                         meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE,
                         h1 = FALSE, baseline = FALSE),
                  error = function(e) { cat("  model failed:", conditionMessage(e), "\n"); NULL })
  cat("  time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  if (is.null(fit)) next
  ok <- lavInspect(fit, "converged") && isTRUE(suppressWarnings(lavInspect(fit, "post.check")))
  fit_tab[[mn]] <- data.frame(model = mn, n = nrow(w), variables = length(vars) + 1,
                              parameters = lavInspect(fit, "npar"), converged_proper = ok,
                              loglik = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit))
  print(fit_tab[[mn]], row.names = FALSE)
  pe <- parameterEstimates(fit)
  cl <- pe[grepl("^c_", pe$label) & pe$op == "~", ]
  cl <- cl[!duplicated(cl$label), ]
  parts <- do.call(rbind, strsplit(sub("^c_", "", cl$label), "__", fixed = TRUE))
  paths <- data.frame(model = mn, from = parts[, 1], to = parts[, 2], est = cl$est,
                      lo = cl$ci.lower, hi = cl$ci.upper, p = cl$pvalue)
  all_paths[[mn]] <- paths
  for (v in vars) {
    f <- paths[paths$from == v & paths$to == "OFF", ]; r <- paths[paths$from == "OFF" & paths$to == v, ]
    fwd[[paste(mn, v)]] <- data.frame(model = mn, variable = v, predictor = unname(lab_map[v]),
      mv_X_to_off = f$est, mv_lo = f$lo, mv_hi = f$hi, mv_p = f$p,
      mv_off_to_X = r$est, mv_lo_r = r$lo, mv_hi_r = r$hi, mv_p_r = r$p)
  }
}

save_tab <- function(t, name, note = NULL) {
  num <- sapply(t, is.numeric); t[num] <- lapply(t[num], round, 3)
  write.csv(t, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat("\n=====", name, "=====\n"); if (!is.null(note)) cat(note, "\n"); print(t, row.names = FALSE)
}
M1 <- do.call(rbind, fwd)
j1_path <- file.path(out_dir, "J1_age_bands_vs_invariant.csv")
if (file.exists(j1_path)) {                      # pairwise models with one cross-lag for all ages (05, part J)
  j1 <- read.csv(j1_path)
  m <- match(M1$variable, j1$variable)
  M1$pair_X_to_off <- j1$inv_X_to_off[m]; M1$pair_lo <- j1$inv_lo[m]; M1$pair_hi <- j1$inv_hi[m]
  M1$pair_off_to_X <- j1$inv_off_to_X[m]; M1$pair_lo_r <- j1$inv_lo_r[m]; M1$pair_hi_r <- j1$inv_hi_r[m]
}
save_tab(M1, "M1_multivariate_cross_lags",
         "Forward (measure -> next year's offending) and reverse cross-lags, one value for all ages; mv_ = multivariate model, pair_ = pairwise model (J1).")
save_tab(do.call(rbind, fit_tab), "M2_multivariate_fit")
save_tab(do.call(rbind, all_paths), "M3_multivariate_all_lagged_paths",
         "Every lagged path between different variables (from -> to), one value for all ages.")
sink()
