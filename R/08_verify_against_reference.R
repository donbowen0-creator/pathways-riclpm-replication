#===============================================================================
# 08_verify_against_reference.R  --  Pathways to Desistance (ICPSR 29961)
#
# Checks your run against the tables behind the submitted manuscript.
#
# Part 1. Every CSV in output_main/ and output_checks/ is compared, row by row,
#         with the reference copy in reference_outputs/. Each file gets a verdict:
#           MATCH   all numbers within .0015 (or 0.1% for large values such as AIC)
#           CLOSE   all numbers within .01 (expected if a package version differs)
#           DIFFERS at least one number off by more than .01 -- look at the file
#           MISSING the file was not produced by your run
# Part 2. About 85 numbers quoted in the manuscript text are recomputed from YOUR
#         tables and compared with the printed value (PASS / FAIL).
# Part 3. Package versions, so the Method section can cite the versions you used.
#
# Requires: output_main/ and output_checks/ from scripts 03 to 06 and 09
# Packages: none beyond base R
# Output  : verification_report.txt (also printed to the console)
# Run     : from the ICPSR_29961 folder, source("<package folder>/R/08_verify_against_reference.R")
#           The reference folder is found automatically when run through 00_RUN_ALL.R;
#           otherwise set ref_dir below.
#===============================================================================

ref_dir <- Sys.getenv("REF_DIR", "pathways-riclpm-replication/reference_outputs")
report  <- "verification_report.txt"
TOL_MATCH <- 0.0015; TOL_CLOSE <- 0.01

if (!dir.exists(ref_dir)) stop("Reference folder not found: ", ref_dir,
                               "\nSet REF_DIR or edit ref_dir at the top of this script.")
sink(report, split = TRUE)
cat("Verification of the desistance analysis\n")
cat("Run on", format(Sys.time(), "%Y-%m-%d %H:%M"), "in", getwd(), "\n")
cat("Reference tables:", normalizePath(ref_dir), "\n\n")

#---- Part 1: file-by-file comparison --------------------------------------------
interval_pat <- "^-?[0-9.]+ \\(-?[0-9.]+ to -?[0-9.]+\\)$"
is_interval_col <- function(x) is.character(x) && all(grepl(interval_pat, x[!is.na(x) & x != ""]))
nums_of <- function(x) lapply(regmatches(x, gregexpr("-?[0-9]*\\.?[0-9]+", x)), as.numeric)

compare_file <- function(sub, f) {
  ref <- read.csv(file.path(ref_dir, sub, f), stringsAsFactors = FALSE)
  mine_path <- file.path(sub, f)
  if (!file.exists(mine_path)) return(data.frame(file = file.path(sub, f), verdict = "MISSING",
                                                 max_diff = NA, rows_unmatched = NA, worst = ""))
  mine <- read.csv(mine_path, stringsAsFactors = FALSE)
  keys <- names(ref)[sapply(names(ref), function(v) (is.character(ref[[v]]) && !is_interval_col(ref[[v]])) ||
                                                      v %in% c("partial_pre", "direction", "timing", "model", "band", "lag"))]
  keys <- intersect(keys, names(mine))
  mk <- function(d) { k <- if (length(keys)) do.call(paste, c(d[keys], sep = " | ")) else rep("", nrow(d))
                      paste(k, ave(seq_along(k), k, FUN = seq_along), sep = " #") }
  rk <- mk(ref); mkk <- mk(mine)
  m <- match(rk, mkk)
  unmatched <- sum(is.na(m))
  worst <- ""; maxd <- 0
  for (v in setdiff(intersect(names(ref), names(mine)), keys)) {
    a <- ref[[v]]; b <- mine[[v]][m]
    if (is.logical(a) || is.logical(b)) {
      bad <- which(!is.na(m) & as.character(a) != as.character(b))
      if (length(bad)) { maxd <- max(maxd, 1); worst <- paste0(v, " differs at ", rk[bad[1]]) }
      next
    }
    if (is.character(a)) { A <- nums_of(a); B <- nums_of(as.character(b))
      d <- mapply(function(p, q) if (length(p) == length(q) && length(p)) max(abs(p - q)) else NA, A, B)
      scale <- sapply(A, function(p) if (length(p)) max(abs(p)) else 0)
    } else { d <- abs(as.numeric(a) - as.numeric(b)); scale <- abs(as.numeric(a)) }
    d_adj <- ifelse(scale > 10, d / scale * TOL_MATCH / 0.001, d)   # large values: 0.1% counts as .0015
    if (any(!is.na(d_adj))) {
      i <- which.max(d_adj)
      if (d_adj[i] > maxd) { maxd <- d_adj[i]; worst <- paste0(v, " at ", rk[i], " (reference ", a[i], ", yours ", b[i], ")") }
    }
  }
  verdict <- if (unmatched > 0) "DIFFERS" else if (maxd <= TOL_MATCH) "MATCH" else if (maxd <= TOL_CLOSE) "CLOSE" else "DIFFERS"
  data.frame(file = file.path(sub, f), verdict = verdict, max_diff = round(maxd, 4),
             rows_unmatched = unmatched, worst = if (verdict == "MATCH") "" else worst)
}

res <- do.call(rbind, lapply(c("output_main", "output_checks"), function(sub)
  do.call(rbind, lapply(list.files(file.path(ref_dir, sub), pattern = "\\.csv$"), function(f) compare_file(sub, f)))))
cat("PART 1. Tables compared with the reference run\n")
cat(sprintf("  %-55s %-8s %9s %10s\n", "file", "verdict", "max diff", "unmatched"))
for (i in seq_len(nrow(res))) {
  cat(sprintf("  %-55s %-8s %9s %10s\n", res$file[i], res$verdict[i],
              ifelse(is.na(res$max_diff[i]), "", format(res$max_diff[i])),
              ifelse(is.na(res$rows_unmatched[i]), "", res$rows_unmatched[i])))
  if (nzchar(res$worst[i])) cat("      largest difference:", res$worst[i], "\n")
}
cat("\n  Summary:", paste(names(table(res$verdict)), table(res$verdict), sep = " = ", collapse = ", "), "\n\n")

#---- Part 2: numbers quoted in the manuscript ------------------------------------
rd <- function(sub, f) { p <- file.path(sub, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }
E0 <- rd("output_main", "E0_fit.csv");  E1 <- rd("output_main", "E1_cross_lags_by_band.csv")
E3 <- rd("output_main", "E3_summary_by_band.csv"); G1 <- rd("output_main", "G1_cross_lags_concurrent.csv")
G2 <- rd("output_main", "G2_concurrent_paths.csv"); H2 <- rd("output_main", "H2_childhood_on_stable_offending.csv")
I1 <- rd("output_main", "I1_dynamic_panel_ml.csv"); J1 <- rd("output_main", "J1_age_bands_vs_invariant.csv")
L1 <- rd("output_main", "L1_cross_lags_secure_time_both.csv"); F2 <- rd("output_main", "F2_secure_time_effect.csv")
S4 <- rd("output_main", "S4_common_state_simulation.csv"); MV <- rd("output_main", "M1_multivariate_cross_lags.csv")
A1 <- rd("output_checks", "A1_positive_control.csv"); S1 <- rd("output_checks", "S1_simulation_lag_artifacts.csv")
S2 <- rd("output_checks", "S2_calibrated_null_cross_lags.csv"); S3 <- rd("output_checks", "S3_hybrid_vs_riclpm_calibrated.csv")

scales <- c("EXPECT", "MORDIS", "SOCCST", "PUNYOU", "PERRWD", "SOCRWD", "PSMIID",
            "IMPULS", "SUPAGG", "CONSID", "FUTURE", "PEERPR", "PSMI")          # 13 attitude and maturity scales
risk   <- c("MORDIS", "PERRWD", "SOCRWD", "PRINFL")                             # theory predicts more offending
all14  <- c(scales, "PRINFL")
bonds4 <- c("WORK_PROP", "SCHOOL_PROP", "PARTNER", "PARENT")
e3  <- function(v, col) E3[E3$variable == v & E3$band == "all", col]
e1  <- function(v, band, col, d = E1) d[d$model == unique(d$model)[1] & d$variable == v & d$band == band, col]
ex0 <- function(lo, hi) (lo > 0 & hi > 0) | (lo < 0 & hi < 0)
first_num <- function(s) as.numeric(sub("^(-?[0-9.]+).*$", "\\1", s))
flip <- function(v) ifelse(v %in% risk, 1, -1)
i1 <- function(pp) { d <- I1[I1$direction == "X_to_off" & as.logical(I1$partial_pre) == pp & I1$variable %in% all14, ]
                     d$b <- d$est * flip(d$variable); d$x0 <- ex0(d$lo, d$hi); d }
j1 <- function() { d <- J1[J1$variable %in% all14, ]; d$b <- d$inv_X_to_off * flip(d$variable)
                   d$x0 <- ex0(d$inv_lo, d$inv_hi); d }

checks <- list(
  list("Analytic sample, N = 1,349",                       function() E0$n[1], 1349, 0),
  list("Scale models: robust CFI, lowest (.967)",          function() min(E0$CFI[!E0$variable %in% bonds4]), .967, 3),
  list("Scale models: robust CFI, highest (.993)",         function() max(E0$CFI[!E0$variable %in% bonds4]), .993, 3),
  list("Scale models: robust RMSEA, lowest (.021)",        function() min(E0$RMSEA[!E0$variable %in% bonds4]), .021, 3),
  list("Scale models: robust RMSEA, highest (.044)",       function() max(E0$RMSEA[!E0$variable %in% bonds4]), .044, 3),
  list("Scale models: SRMR, lowest (.034)",                function() min(E0$SRMR[!E0$variable %in% bonds4]), .034, 3),
  list("Scale models: SRMR, highest (.051)",               function() max(E0$SRMR[!E0$variable %in% bonds4]), .051, 3),
  list("Bond models: robust CFI, lowest (.938)",           function() min(E0$CFI[E0$variable %in% bonds4]), .938, 3),
  list("Bond models: robust CFI, highest (1.000)",         function() max(E0$CFI[E0$variable %in% bonds4]), 1.000, 3),
  list("Bond models: robust RMSEA, lowest (.000)",         function() min(E0$RMSEA[E0$variable %in% bonds4]), .000, 3),
  list("Bond models: robust RMSEA, highest (.059)",        function() max(E0$RMSEA[E0$variable %in% bonds4]), .059, 3),
  list("Bond models: SRMR, lowest (.035)",                 function() min(E0$SRMR[E0$variable %in% bonds4]), .035, 3),
  list("Bond models: SRMR, highest (.059)",                function() max(E0$SRMR[E0$variable %in% bonds4]), .059, 3),
  list("Scale stability, lowest (.13)",                    function() min(E3$stability_X[E3$band == "all" & E3$variable %in% all14]), .13, 2),
  list("Scale stability, highest (.24)",                   function() max(E3$stability_X[E3$band == "all" & E3$variable %in% all14]), .24, 2),
  list("Offending stability, lowest (.15)",                function() min(E3$stability_off[E3$band == "all"]), .15, 2),
  list("Offending stability, highest (.18)",               function() max(E3$stability_off[E3$band == "all"]), .18, 2),
  list("Parenthood stability (.76)",                       function() e3("PARENT", "stability_X"), .76, 2),
  list("RI correlation, peer influence (.70)",             function() e3("PRINFL", "RI_corr"), .70, 2),
  list("RI correlation, suppression of aggression (-.51)",  function() e3("SUPAGG", "RI_corr"), -.51, 2),
  list("Within-person share, lowest scale (.47)",          function() min(E3$share_within_X[E3$band == "all" & E3$variable %in% all14]), .47, 2),
  list("Within-person share, highest scale (.75)",         function() max(E3$share_within_X[E3$band == "all" & E3$variable %in% all14]), .75, 2),
  list("Same-age r, peer influence (.36)",                 function() e3("PRINFL", "same_age_r"), .36, 2),
  list("Same-age r, impulse control (-.18)",               function() e3("IMPULS", "same_age_r"), -.18, 2),
  list("Same-age r, suppression of aggression (-.17)",     function() e3("SUPAGG", "same_age_r"), -.17, 2),
  list("Same-age r, personal rewards (.15)",               function() e3("PERRWD", "same_age_r"), .15, 2),
  list("Same-age r, moral disengagement (.14)",            function() e3("MORDIS", "same_age_r"), .14, 2),
  list("Same-age r, social rewards (.13)",                 function() e3("SOCRWD", "same_age_r"), .13, 2),
  list("Same-age r, expectations (-.12)",                  function() e3("EXPECT", "same_age_r"), -.12, 2),
  list("Forward p < .05 among 39 scale estimates (6)",     function() sum(E1$model == "main" & E1$variable %in% scales & E1$p < .05), 6, 0),
  list("Forward, moral disengagement 16-18 (.088)",        function() e1("MORDIS", "16-18", "X_to_off"), .088, 3),
  list("Forward, moral disengagement 21-23 (.076)",        function() e1("MORDIS", "21-23", "X_to_off"), .076, 3),
  list("Forward, personal rewards 16-18 (.122)",           function() e1("PERRWD", "16-18", "X_to_off"), .122, 3),
  list("Forward, social rewards 16-18 (.078)",             function() e1("SOCRWD", "16-18", "X_to_off"), .078, 3),
  list("Forward, impulse control 21-23 (-.065)",           function() e1("IMPULS", "21-23", "X_to_off"), -.065, 3),
  list("Forward, certainty 18-21 (.046)",                  function() e1("PUNYOU", "18-21", "X_to_off"), .046, 3),
  list("Forward, certainty 18-21, p (.049)",               function() e1("PUNYOU", "18-21", "p"), .049, 3),
  list("Forward, peer influence 18-21 (.073)",             function() e1("PRINFL", "18-21", "X_to_off"), .073, 3),
  list("Reverse p < .05 at 16-18 among 13 scales (8)",     function() sum(E1$model == "main" & E1$variable %in% scales & E1$band == "16-18" & E1$p_r < .05), 8, 0),
  list("Reverse, offending to employment 21-23 (-.044)",   function() e1("WORK_PROP", "21-23", "off_to_X"), -.044, 3),
  list("Reverse, offending to partner 21-23 (-.051)",      function() e1("PARTNER", "21-23", "off_to_X"), -.051, 3),
  list("Secure time both: employment 21-23 (-.027)",       function() e1("WORK_PROP", "21-23", "off_to_X", L1), -.027, 3),
  list("Secure time both: school 16-18 (-.015)",           function() e1("SCHOOL_PROP", "16-18", "off_to_X", L1), -.015, 3),
  list("Secure time both: partner 21-23 (-.039)",          function() e1("PARTNER", "21-23", "off_to_X", L1), -.039, 3),
  list("Secure-time effect on offending, smallest (-.23)", function() max(F2$secure_time_effect), -.23, 2),
  list("Secure-time effect on offending, largest (-.40)",  function() min(F2$secure_time_effect), -.40, 2),
  list("Same-age path model: personal rewards 16-18 (.082)", function() e1("PERRWD", "16-18", "X_to_off", G1), .082, 3),
  list("Same-age directed path, peer influence (.33)",     function() G2$same_age_path[G2$variable == "PRINFL" & G2$direction == "xy"], .33, 2),
  list("Same-age directed path, impulse control (-.20)",   function() G2$same_age_path[G2$variable == "IMPULS" & G2$direction == "xy"], -.20, 2),
  list("Childhood: no early behavior problems (-.277)",    function() H2$beta[H2$factor == "cp_noearly"], -.277, 3),
  list("Stable offending R-squared (.19)",                 function() H2$R2_RIy[1], .19, 2),
  list("BIC favors one cross-lag for all 18 pairs (18)",   function() sum(J1$dBIC_bands_vs_invariant > 0), 18, 0),
  list("BIC difference, smallest (4.5)",                   function() min(J1$dBIC_bands_vs_invariant), 4.5, 1),
  list("BIC difference, largest (27.4)",                   function() max(J1$dBIC_bands_vs_invariant), 27.4, 1),
  list("AIC favors bands (7 pairs)",                       function() sum(J1$dAIC_bands_vs_invariant < 0), 7, 0),
  list("One cross-lag: personal rewards (.056)",           function() J1$inv_X_to_off[J1$variable == "PERRWD"], .056, 3),
  list("One cross-lag: peer influence (.062)",             function() J1$inv_X_to_off[J1$variable == "PRINFL"], .062, 3),
  list("One cross-lag: suppression of aggression (-.037)", function() J1$inv_X_to_off[J1$variable == "SUPAGG"], -.037, 3),
  list("One cross-lag: reverse significant, 13 scales (9)", function() { d <- J1[J1$variable %in% scales, ]; sum(ex0(d$inv_lo_r, d$inv_hi_r)) }, 9, 0),
  list("Lag-only dynamic panel: wrong side of zero (14)",  function() sum(i1(FALSE)$b < 0), 14, 0),
  list("Lag-only dynamic panel: CI excludes zero (8)",     function() sum(i1(FALSE)$x0), 8, 0),
  list("Same-year link allowed: theory direction (10)",    function() sum(i1(TRUE)$b > 0.0005), 10, 0),
  list("Same-year link allowed: CI excludes zero (7)",     function() sum(i1(TRUE)$x0), 7, 0),
  list("Same-year link allowed: personal rewards (.098)",  function() I1$est[I1$variable == "PERRWD" & I1$direction == "X_to_off" & as.logical(I1$partial_pre)], .098, 3),
  list("One cross-lag RI-CLPM: expected direction (11)",   function() sum(j1()$b > 0), 11, 0),
  list("One cross-lag RI-CLPM: CI excludes zero (3)",      function() sum(j1()$x0), 3, 0),
  list("Earlier hybrid model, peer influence OR (0.89)",   function() A1$OR_within[A1$timing == "Previous year"][1], 0.89, 2),
  list("Simulation 1: hybrid OR, swings do not persist (0.56)", function() first_num(S1$hybrid_OR_with_prior_offending[2]), 0.56, 2),
  list("Simulation 1: RI-CLPM, swings do not persist (-.15)",  function() first_num(S1$riclpm_cross_lag[2]), -.15, 2),
  list("Simulation 1: hybrid OR, swings persist (0.84)",   function() first_num(S1$hybrid_OR_with_prior_offending[3]), 0.84, 2),
  list("Simulation 1: RI-CLPM, swings persist (.10)",      function() first_num(S1$riclpm_cross_lag[3]), .10, 2),
  list("Simulation 2: null cross-lags, lowest (-.009)",    function() min(S2$null_mean), -.009, 3),
  list("Simulation 2: null cross-lags, highest (.007)",    function() max(S2$null_mean), .007, 3),
  list("Simulation 2: hybrid OR, peer influence (0.77)",   function() S3$hybrid_OR_mean[S3$predictor == "Peer antisocial influence" & S3$true_lagged_effect == 0], 0.77, 2),
  list("Simulation 3: spurious cross-lag, share .15 (.037)", function() first_num(S4$riclpm_cross_lag[S4$persistence == 0.5 & S4$share_x == 0.15 & S4$share_y == 0.15]), .037, 3),
  list("Simulation 3: spurious cross-lag, share .30 (.084)", function() first_num(S4$riclpm_cross_lag[S4$persistence == 0.5 & S4$share_x == 0.30 & S4$share_y == 0.30]), .084, 3),
  list("Multivariate, cognitive: personal rewards forward (.050)", function() MV$mv_X_to_off[MV$model == "M1_cognitive" & MV$variable == "PERRWD"], .050, 3),
  list("Multivariate, across domains: personal rewards forward (.048)", function() MV$mv_X_to_off[MV$model == "M4_across" & MV$variable == "PERRWD"], .048, 3),
  list("Multivariate, across domains: peer influence forward (.061)", function() MV$mv_X_to_off[MV$model == "M4_across" & MV$variable == "PRINFL"], .061, 3),
  list("Multivariate, maturation: suppression of aggression forward (-.029)", function() MV$mv_X_to_off[MV$model == "M2_maturation" & MV$variable == "SUPAGG"], -.029, 3),
  list("Multivariate, forward CI excludes zero, all models (3)", function() sum(ex0(MV$mv_lo, MV$mv_hi)), 3, 0),
  list("Multivariate, domain models: reverse paths kept by scales (3)", function() { d <- MV[MV$model %in% c("M1_cognitive", "M2_maturation"), ]; sum(ex0(d$mv_lo_r, d$mv_hi_r)) }, 3, 0),
  list("Multivariate, reverse: expectations (-.033)", function() MV$mv_off_to_X[MV$model == "M1_cognitive" & MV$variable == "EXPECT"], -.033, 3),
  list("Multivariate, bonds: offending to employment (-.016)", function() MV$mv_off_to_X[MV$model == "M3_bonds" & MV$variable == "WORK_PROP"], -.016, 3),
  list("Multivariate, bonds: offending to school (-.007)", function() MV$mv_off_to_X[MV$model == "M3_bonds" & MV$variable == "SCHOOL_PROP"], -.007, 3)
)

rhu <- function(x, d) sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d   # round half up, as printed
cat("PART 2. Numbers quoted in the manuscript, recomputed from your tables\n")
n_fail <- 0
for (ck in checks) {
  val <- tryCatch(ck[[2]](), error = function(e) NA)
  val <- if (length(val) == 0) NA else val[1]
  ok  <- !is.na(val) && rhu(val, ck[[4]]) == rhu(ck[[3]], ck[[4]])
  if (!ok) n_fail <- n_fail + 1
  cat(sprintf("  %-4s %-58s printed %8s   yours %8s\n", if (ok) "PASS" else "FAIL", ck[[1]],
              format(ck[[3]]), if (is.na(val)) "missing" else format(round(val, ck[[4]] + 1))))
}
cat(sprintf("\n  %d of %d checks passed.\n", length(checks) - n_fail, length(checks)))
if (n_fail > 0) cat("  A FAIL means the printed number in the manuscript should be updated from your tables,\n",
                    " or the difference should be explained.\n")

#---- Part 3: versions --------------------------------------------------------------
cat("\nPART 3. Versions used in this run (the manuscript cites R 4.6.1 and lavaan 0.7-2)\n")
cat("  ", R.version.string, "\n")
for (p in c("lavaan", "dpm", "panelr", "mice", "fixest", "sandwich", "lmtest", "ggplot2")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
  cat(sprintf("   %-9s %s\n", p, v))
}
cat("\n  If your R or lavaan version differs, the Method sentence and the R Core Team reference\n",
    " should cite your versions.\n")
sink()
cat("\nReport saved to", normalizePath(report), "\n")
