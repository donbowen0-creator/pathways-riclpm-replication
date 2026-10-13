#===============================================================================
# 05_riclpm_main.R  --  Pathways to Desistance (ICPSR 29961)
# Main analysis for the RI-CLPM version of the desistance paper (Draft 1, 2026-10-03)
#
# Each model pairs ONE predictor (X) with offending variety (Y). Interviews are
# placed at the youth's age (baseline age + years since baseline), ages 16-23,
# so every lag is one year. Random intercepts hold each person's stable level of
# X and Y; the cross-lags are among within-person deviations from those levels.
#
# E. Main RI-CLPMs: 14 attitude and maturity scales, 4 social bonds.
#    Free model (fit; per-lag standardized estimates) and an age-band model
#    (cross-lags equal within ages 16-18, 18-21, 21-23; Chen & Sutton 2025).
# F. Street-time adjustment: share of each recall period spent in a secure
#    facility as a time-varying covariate of offending.
# G. Concurrent-path check (Speyer et al. 2025; Usami 2021): a same-age directed
#    path X -> Y (and, separately, Y -> X) replaces the same-age residual covariance.
# H. Childhood promotive factors and the random intercepts: do stable childhood
#    factors explain the stable, between-person part of offending and of each
#    predictor?
# I. Maximum-likelihood dynamic panel models (Allison, Williams & Moral-Benito
#    2017; dpm package): person fixed effects, lagged offending, lagged predictor
#    treated as predetermined (offending may feed back into later X).
#
# Requires: panel_desistance.rds from 01_data_prep.R (2026-10-02 version)
# Packages: install.packages(c("lavaan", "dpm", "panelr"))
# Output  : output_main/  (CSV tables + run_log.txt)
# Run     : Rscript 05_riclpm_main.R           (about 15-30 minutes)
#===============================================================================

suppressPackageStartupMessages({ library(lavaan); library(dpm); library(panelr) })

data_dir <- "."
out_dir  <- file.path(data_dir, "output_main")
dir.create(out_dir, showWarnings = FALSE)
set.seed(20261003)

panel <- readRDS(file.path(data_dir, "panel_desistance.rds"))
if (is.null(panel$YRS)) panel$YRS <- panel$year
panel$AGE <- panel$BL_AGE + panel$YRS
AGES <- 16:23
ids  <- sort(unique(panel$CASEID))

lab_map <- c(WORK_PROP = "Employment (share of period)", SCHOOL_PROP = "School (share of period)",
             PARTNER = "Living with partner", PARENT = "Parent",
             EXPECT = "Expectations for the future", MORDIS = "Moral disengagement",
             SOCCST = "Social costs of punishment", PUNYOU = "Certainty of punishment",
             PERRWD = "Personal rewards of crime", SOCRWD = "Social rewards of crime",
             PSMIID = "Identity (PSMI)",
             IMPULS = "Impulse control", SUPAGG = "Suppression of aggression",
             CONSID = "Consideration of others", FUTURE = "Future orientation",
             PEERPR = "Resistance to peer influence", PSMI = "Psychosocial maturity (PSMI)",
             PRINFL = "Peer antisocial influence")
domain_map <- c(EXPECT = "Cognitive", MORDIS = "Cognitive", SOCCST = "Cognitive", PUNYOU = "Cognitive",
                PERRWD = "Cognitive", SOCRWD = "Cognitive", PSMIID = "Cognitive",
                IMPULS = "Maturation", SUPAGG = "Maturation", CONSID = "Maturation",
                FUTURE = "Maturation", PEERPR = "Maturation", PSMI = "Maturation",
                PRINFL = "Positive control",
                WORK_PROP = "Bonds", SCHOOL_PROP = "Bonds", PARTNER = "Bonds", PARENT = "Bonds")
lab <- function(v) unname(lab_map[v])

scales <- c("EXPECT", "MORDIS", "SOCCST", "PUNYOU", "PERRWD", "SOCRWD", "PSMIID",
            "IMPULS", "SUPAGG", "CONSID", "FUTURE", "PEERPR", "PSMI", "PRINFL")
bonds  <- c("WORK_PROP", "SCHOOL_PROP", "PARTNER", "PARENT")
all_x  <- c(scales, bonds)
args <- commandArgs(trailingOnly = TRUE)        # optional: run only the named predictors (quick test)
if (length(args) > 0) all_x <- intersect(all_x, args)

save_tab <- function(t, name, note = NULL) {
  num <- sapply(t, is.numeric); t[num] <- lapply(t[num], round, 3)
  write.csv(t, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat("\n=====", name, "=====\n"); if (!is.null(note)) cat(note, "\n")
  print(t, row.names = FALSE)
}
# Which parts to run (default all). Example: RUN_PARTS=JK Rscript 05_riclpm_main.R
RUN <- strsplit(Sys.getenv("RUN_PARTS", "EFGHIJKL"), "")[[1]]
sink(file.path(out_dir, paste0("run_log", if (length(RUN) < 8) paste0("_", paste(RUN, collapse = "")) else "", ".txt")), split = TRUE)
cat("Pathways RI-CLPM main analysis --", format(Sys.time()), "\n")

#-------------------------------------------------------------------------------
# Data: one wide row per person, x16..x23, y16..y23, s16..s23 (street-time control)
#-------------------------------------------------------------------------------
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
panel$Y_z <- zs(panel$OFFVAR)             # offending variety, per pooled SD

make_wide <- function(v) {
  x <- if (v %in% bonds) panel[[v]] else zs(panel[[v]])   # scales per pooled SD; bonds on 0-1 scale
  w <- data.frame(CASEID = ids)
  for (a in AGES) {
    r <- which(panel$AGE == a); m <- match(ids, panel$CASEID[r])
    w[[paste0("x", a)]] <- x[r][m]
    w[[paste0("y", a)]] <- panel$Y_z[r][m]
    w[[paste0("s", a)]] <- panel$PROP_SECURE[r][m]
  }
  cp <- panel[panel$year == 0, ][match(ids, panel$CASEID[panel$year == 0]), ]
  w$cp_noearly <- zs(cp$CP_NO_EARLY_PROBLEMS); w$cp_twobio <- cp$CP_TWO_BIO_PARENTS
  w$cp_pareduc <- zs(cp$CP_PARENT_EDUC);       w$cp_grades <- zs(cp$CP_SCHOOL_GRADES)
  w$cp_iq      <- zs(cp$CP_IQ);                w$cp_noparrest <- cp$CP_NO_PARENT_ARREST
  w$male <- as.numeric(cp$BL_SEX == 1); w$phoenix <- as.numeric(cp$BL_SITE == 2)
  w$black <- as.numeric(cp$BL_RACE == 2); w$hispanic <- as.numeric(cp$BL_RACE == 3)
  w$other <- as.numeric(cp$BL_RACE == 4); w$bl_age <- cp$BL_AGE
  keep <- rowSums(!is.na(w[, c(paste0("x", AGES), paste0("y", AGES))])) > 0
  w[keep, ]
}

band_of <- function(a) ifelse(a <= 18, "16-18", ifelse(a <= 21, "18-21", "21-23"))   # by the lag's END age
cp_vars   <- c("cp_noearly", "cp_twobio", "cp_pareduc", "cp_grades", "cp_iq", "cp_noparrest")
ctrl_vars <- c("male", "phoenix", "black", "hispanic", "other", "bl_age")

#-------------------------------------------------------------------------------
# RI-CLPM syntax
#   bands  : cross-lags equal within age bands
#   conc   : "none" | "xy" (same-age X -> Y) | "yx" (same-age Y -> X)
#   street : offending regressed on same-age share of time in secure care
#   ri_on  : time-invariant predictors of both random intercepts
#-------------------------------------------------------------------------------
riclpm_syntax <- function(bands = TRUE, conc = "none", street = FALSE, ri_on = NULL) {
  s <- c(paste("RIx =~", paste0("1*x", AGES, collapse = " + ")),
         paste("RIy =~", paste0("1*y", AGES, collapse = " + ")))
  for (a in AGES) s <- c(s, sprintf("wx%d =~ 1*x%d", a, a), sprintf("wy%d =~ 1*y%d", a, a),
                         sprintf("x%d ~~ 0*x%d", a, a), sprintf("y%d ~~ 0*y%d", a, a))
  for (i in 2:length(AGES)) {
    a <- AGES[i]; p <- AGES[i - 1]
    tag <- if (bands) gsub("-", "_", band_of(a)) else a
    xs <- sprintf("wx%d ~ arx%d*wx%d + cyx%s*wy%d", a, a, p, tag, p)
    ys <- sprintf("wy%d ~ ary%d*wy%d + cxy%s*wx%d", a, a, p, tag, p)
    if (conc == "xy") ys <- paste0(ys, sprintf(" + cc*wx%d", a))
    if (conc == "yx") xs <- paste0(xs, sprintf(" + cc*wy%d", a))
    s <- c(s, xs, ys)
  }
  cov_ages <- if (conc == "none") AGES else AGES[1]
  s <- c(s, sprintf("wx%d ~~ cov%d*wy%d", cov_ages, cov_ages, cov_ages),
         sprintf("wx%d ~~ wx%d", AGES, AGES), sprintf("wy%d ~~ wy%d", AGES, AGES),
         "RIx ~~ RIx", "RIy ~~ RIy", "RIx ~~ RIy",
         sprintf("RIx ~~ 0*wx%d", AGES[1]), sprintf("RIx ~~ 0*wy%d", AGES[1]),
         sprintf("RIy ~~ 0*wx%d", AGES[1]), sprintf("RIy ~~ 0*wy%d", AGES[1]))
  if (!isFALSE(street)) {
    # Share of the recall period in secure care: its own random intercept (stable
    # differences in confinement, correlated with stable X and stable offending)
    # plus a same-age effect on offending, equal across ages.
    s <- c(s, paste("RIs =~", paste0("1*s", AGES, collapse = " + ")),
           sprintf("y%d ~ st*s%d", AGES, AGES), sprintf("s%d ~~ s%d", AGES, AGES),
           "RIs ~~ RIs", "RIs ~~ RIx + RIy",
           sprintf("RIs ~~ 0*wx%d", AGES[1]), sprintf("RIs ~~ 0*wy%d", AGES[1]))
    # street = "both": secure time also enters the predictor's equation, so time spent
    # locked up cannot pass for an effect of offending on work, school or partnering
    if (identical(street, "both")) s <- c(s, sprintf("x%d ~ sx*s%d", AGES, AGES))
  }
  if (!is.null(ri_on)) s <- c(s, paste("RIx ~", paste(ri_on, collapse = " + ")),
                              paste("RIy ~", paste(ri_on, collapse = " + ")))
  paste(s, collapse = "\n")
}

# FIML for missing ages. Models with time-invariant covariates (section H) keep the
# covariates fixed and use missing = "ml.x" so the few people missing a childhood
# item stay in the model.
fit_ri <- function(w, ..., ri_on = NULL) {
  miss <- if (is.null(ri_on)) "ml" else "ml.x"
  tryCatch(lavaan(riclpm_syntax(..., ri_on = ri_on), data = w, missing = miss, estimator = "MLR",
                  meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE, fixed.x = TRUE),
           error = function(e) { cat("  model failed:", conditionMessage(e), "\n"); NULL })
}

band_rows <- function(fit, v, model) {
  if (is.null(fit)) return(NULL)
  pe <- parameterEstimates(fit)
  do.call(rbind, lapply(c("16_18", "18_21", "21_23"), function(b) {
    xy <- pe[pe$label == paste0("cxy", b), ][1, ]; yx <- pe[pe$label == paste0("cyx", b), ][1, ]
    data.frame(model = model, variable = v, predictor = lab(v), domain = domain_map[[v]],
               band = sub("_", "-", b),
               X_to_off = xy$est, lo = xy$ci.lower, hi = xy$ci.upper, p = xy$pvalue,
               off_to_X = yx$est, lo_r = yx$ci.lower, hi_r = yx$ci.upper, p_r = yx$pvalue)
  }))
}

if (any(c("E", "F", "G", "H") %in% RUN)) {
#===============================================================================
# E + F + G. RI-CLPMs for every predictor
#===============================================================================
E_fit <- list(); E_band <- list(); E_lag <- list(); F_band <- list(); F_st <- list(); G_band <- list(); G_cc <- list()
H_ri <- list(); wides <- list()

for (v in all_x) {
  cat("\n---", lab(v), "---\n")
  w <- make_wide(v); wides[[v]] <- w

  # E. free model (fit + per-lag standardized) and band model
  f_free <- fit_ri(w, bands = FALSE)
  f_band <- fit_ri(w, bands = TRUE)
  if (!is.null(f_free)) {
    fm <- fitMeasures(f_free, c("cfi.robust", "rmsea.robust", "srmr"))
    E_fit[[v]] <- data.frame(variable = v, predictor = lab(v), domain = domain_map[[v]],
                             n = lavInspect(f_free, "nobs"), converged = lavInspect(f_free, "converged"),
                             CFI = fm[[1]], RMSEA = fm[[2]], SRMR = fm[[3]])
    ss <- standardizedSolution(f_free)
    pe <- parameterEstimates(f_free)
    ri_r <- ss[ss$lhs == "RIx" & ss$rhs == "RIy" & ss$op == "~~", ]
    vr  <- function(n) pe$est[pe$lhs == n & pe$rhs == n & pe$op == "~~"]
    E_lag[[v]] <- do.call(rbind, lapply(AGES[-1], function(a) {
      g <- function(l, r, op) ss[ss$lhs == l & ss$rhs == r & ss$op == op, ]
      xy <- g(paste0("wy", a), paste0("wx", a - 1), "~"); yx <- g(paste0("wx", a), paste0("wy", a - 1), "~")
      cw <- g(paste0("wx", a), paste0("wy", a), "~~")
      sx <- g(paste0("wx", a), paste0("wx", a - 1), "~"); sy <- g(paste0("wy", a), paste0("wy", a - 1), "~")
      data.frame(variable = v, predictor = lab(v), domain = domain_map[[v]],
                 lag = paste0(a - 1, "->", a), band = band_of(a),
                 stability_X = sx$est.std, stability_off = sy$est.std,
                 same_age_r = cw$est.std, same_age_lo = cw$ci.lower, same_age_hi = cw$ci.upper,
                 X_to_off = xy$est.std, X_to_off_p = xy$pvalue,
                 off_to_X = yx$est.std, off_to_X_p = yx$pvalue,
                 RI_corr = ri_r$est.std,
                 share_within_X = 1 - vr("RIx") / (vr("RIx") + mean(sapply(AGES, function(k) vr(paste0("wx", k))))),
                 share_within_off = 1 - vr("RIy") / (vr("RIy") + mean(sapply(AGES, function(k) vr(paste0("wy", k))))))
    }))
  }
  E_band[[v]] <- band_rows(f_band, v, "main")

  # F. street-time adjustment
  f_st <- fit_ri(w, bands = TRUE, street = TRUE)
  F_band[[v]] <- band_rows(f_st, v, "street_time")
  if (!is.null(f_st)) {
    r <- parameterEstimates(f_st); r <- r[r$label == "st", ][1, ]
    F_st[[v]] <- data.frame(variable = v, predictor = lab(v), secure_time_effect = r$est, lo = r$ci.lower, hi = r$ci.upper)
  }

  # G. concurrent-path checks
  for (cc in c("xy", "yx")) {
    f_cc <- fit_ri(w, bands = TRUE, conc = cc)
    G_band[[paste(v, cc)]] <- band_rows(f_cc, v, paste0("concurrent_", cc))
    if (!is.null(f_cc)) {
      pe <- parameterEstimates(f_cc); r <- pe[pe$label == "cc", ][1, ]
      G_cc[[paste(v, cc)]] <- data.frame(variable = v, predictor = lab(v), direction = cc,
                                         same_age_path = r$est, lo = r$ci.lower, hi = r$ci.upper, p = r$pvalue)
    }
  }

  # H. childhood promotive factors (and baseline controls) -> both random intercepts
  f_h <- fit_ri(w, bands = TRUE, ri_on = c(cp_vars, ctrl_vars))
  if (!is.null(f_h)) {
    ss <- standardizedSolution(f_h)
    r2 <- lavInspect(f_h, "r2")
    h <- ss[ss$op == "~" & ss$lhs %in% c("RIx", "RIy") & ss$rhs %in% cp_vars, ]
    H_ri[[v]] <- data.frame(variable = v, predictor = lab(v), intercept = h$lhs, factor = h$rhs,
                            beta = h$est.std, lo = h$ci.lower, hi = h$ci.upper, p = h$pvalue,
                            R2 = ifelse(h$lhs == "RIx", r2[["RIx"]], r2[["RIy"]]))
    # did adding them change the cross-lags? (they shouldn't: the cross-lags are within-person)
    E_band[[paste(v, "H")]] <- band_rows(f_h, v, "with_childhood_on_intercepts")
  }
}

E_lag_tab <- do.call(rbind, E_lag)
save_tab(do.call(rbind, E_fit), "E0_fit", "Free RI-CLPM fit (robust CFI/RMSEA, SRMR)")
save_tab(do.call(rbind, E_band), "E1_cross_lags_by_band",
         paste("Within-person cross-lags, equal within age bands (by the lag's end age).",
               "Scales: SD of offending per SD of X. Bonds: SD of offending per 0-1 change in X.",
               "off_to_X: SD of X (bonds: 0-1 units) per SD of offending."))
save_tab(E_lag_tab, "E2_each_lag_standardized",
         "Free model, fully standardized: stability, same-age within-person correlation, cross-lags")
E3 <- do.call(rbind, lapply(split(E_lag_tab, E_lag_tab$variable), function(t) {
  do.call(rbind, lapply(c("all", "16-18", "18-21", "21-23"), function(b) {
    u <- if (b == "all") t else t[t$band == b, ]
    data.frame(variable = t$variable[1], predictor = t$predictor[1], domain = t$domain[1], band = b,
               lags = nrow(u), stability_X = mean(u$stability_X), stability_off = mean(u$stability_off),
               same_age_r = mean(u$same_age_r), same_age_sig = sum(u$same_age_lo > 0 | u$same_age_hi < 0),
               X_to_off = mean(u$X_to_off), off_to_X = mean(u$off_to_X),
               ratio_same_to_cross = abs(mean(u$same_age_r)) / max(abs(mean(u$X_to_off)), 1e-6),
               RI_corr = t$RI_corr[1], share_within_X = t$share_within_X[1], share_within_off = t$share_within_off[1])
  }))
}))
save_tab(E3, "E3_summary_by_band", "Averages of the free model's standardized estimates")
save_tab(do.call(rbind, F_band), "F1_cross_lags_street_time",
         "Band model with offending regressed on same-age share of time in secure care")
save_tab(do.call(rbind, F_st), "F2_secure_time_effect",
         "Same-age effect of a full recall period in secure settings on offending (SD units)")
save_tab(do.call(rbind, G_band), "G1_cross_lags_concurrent", "Band model with a same-age directed path")
save_tab(do.call(rbind, G_cc), "G2_concurrent_paths",
         "Same-age path (xy: X -> offending; yx: offending -> X), equal across ages, unstandardized")
save_tab(do.call(rbind, H_ri), "H1_childhood_on_intercepts",
         "Standardized effects of childhood promotive factors on the random intercepts (controls: sex, site, race/ethnicity, baseline age)")

}  # end E-H
if ("H" %in% RUN) {
#===============================================================================
# H2. Offending alone: random-intercept autoregressive model, childhood factors
#     predicting the stable offending level (one clean estimate for the paper)
#===============================================================================
w <- if (exists("wides") && length(wides)) wides[[1]] else make_wide(all_x[1])
s <- c(paste("RIy =~", paste0("1*y", AGES, collapse = " + ")))
for (a in AGES) s <- c(s, sprintf("wy%d =~ 1*y%d", a, a), sprintf("y%d ~~ 0*y%d", a, a), sprintf("wy%d ~~ wy%d", a, a))
for (i in 2:length(AGES)) s <- c(s, sprintf("wy%d ~ ary%d*wy%d", AGES[i], AGES[i], AGES[i - 1]))
s <- c(s, sprintf("RIy ~~ 0*wy%d", AGES[1]), "RIy ~~ RIy",
       paste("RIy ~", paste(c(cp_vars, ctrl_vars), collapse = " + ")))
f_h2 <- lavaan(paste(s, collapse = "\n"), data = w, missing = "ml.x", estimator = "MLR",
               meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE, fixed.x = TRUE)
ss <- standardizedSolution(f_h2)
h2 <- ss[ss$op == "~" & ss$lhs == "RIy" & ss$rhs %in% c(cp_vars, ctrl_vars), ]
h2 <- data.frame(factor = h2$rhs, beta = h2$est.std, lo = h2$ci.lower, hi = h2$ci.upper, p = h2$pvalue)
h2$R2_RIy <- lavInspect(f_h2, "r2")[["RIy"]]
save_tab(h2, "H2_childhood_on_stable_offending",
         "Random-intercept model of offending, ages 16-23: standardized effects on the stable offending level")

}  # end H2
if ("I" %in% RUN) {
#===============================================================================
# I. Maximum-likelihood dynamic panel models (dpm), ages 16-23
#    Y(a) = b1*Y(a-1) + b2*X(a-1) + alpha_i + e(a); X predetermined.
#    partial.pre = TRUE also lets X(a-1) correlate with the same-period error.
#===============================================================================
long <- panel[panel$AGE %in% AGES, ]
grid <- expand.grid(CASEID = ids, AGE = AGES)
long <- merge(grid, long, by = c("CASEID", "AGE"), all.x = TRUE)
for (v in scales) long[[v]] <- zs(panel[[v]])[match(paste(long$CASEID, long$AGE), paste(panel$CASEID, panel$AGE))]
long$Y <- panel$Y_z[match(paste(long$CASEID, long$AGE), paste(panel$CASEID, panel$AGE))]
pd <- panel_data(long[, c("CASEID", "AGE", "Y", all_x)], id = CASEID, wave = AGE)

# dpm() builds the Allison et al. model; it is fitted here with lavaan directly
# (print.only = TRUE gives the syntax), which avoids a dpm/lavaan version clash.
wide_dpm <- function(dv, xv) {
  w <- data.frame(CASEID = ids)
  for (k in seq_along(AGES)) {
    r <- long[long$AGE == AGES[k], ]; m <- match(ids, r$CASEID)
    w[[paste0(dv, "_", k)]] <- r[[dv]][m]; w[[paste0(xv, "_", k)]] <- r[[xv]][m]
  }
  w
}
I_tab <- list()
for (v in all_x) {
  for (pp in c(FALSE, TRUE)) {
    for (dir in c("X_to_off", "off_to_X")) {
      dv <- if (dir == "X_to_off") "Y" else v; xv <- if (dir == "X_to_off") v else "Y"
      syn <- paste(capture.output(dpm(as.formula(paste(dv, "~ pre(lag(", xv, "))")), data = pd,
                                      partial.pre = pp, print.only = TRUE)), collapse = "\n")
      wd <- wide_dpm(dv, xv)
      m <- tryCatch(sem(syn, data = wd, missing = "ml", estimator = "MLR", fixed.x = FALSE),
                    error = function(e) tryCatch(sem(syn, data = wd, missing = "ml", estimator = "MLR",
                                                     fixed.x = FALSE, h1 = FALSE, baseline = FALSE),
                                                 error = function(e2) { cat("  dpm failed:", v, dir, conditionMessage(e2), "\n"); NULL }))
      if (is.null(m)) next
      pe <- parameterEstimates(m); r <- pe[pe$label == "en1", ][1, ]
      I_tab[[paste(v, pp, dir)]] <- data.frame(variable = v, predictor = lab(v), direction = dir,
                                               partial_pre = pp, est = r$est, lo = r$ci.lower, hi = r$ci.upper,
                                               p = r$pvalue, converged = lavInspect(m, "converged"))
    }
  }
}
save_tab(do.call(rbind, I_tab), "I1_dynamic_panel_ml",
         paste("Allison-style ML dynamic panel models with person fixed effects, ages 16-23.",
               "Scales per pooled SD; bonds 0-1. Lagged-predictor coefficient, one effect for all ages."))

}  # end I
#===============================================================================
# J. Do the age bands earn their place? Age-invariant vs age-band vs free
#    cross-lags, compared by AIC and BIC (lower = better; differences > 10 are
#    large). Chen et al. (2026) kept age-invariant cross-lags in Pathways.
# K. Random-slope check (Andersen 2022; Sorjonen et al. 2026): adds a person-
#    specific linear slope over age for both X and offending, so slow personal
#    trends cannot pass for cross-lags (an RI-CLPM with structured residuals
#    around individual trajectories, as in the LCM-SR / ALT-SR family).
#===============================================================================
if (any(c("J", "K") %in% RUN)) {
  invariant_syntax <- function() gsub("cxy(16_18|18_21|21_23)", "cxyall",
                         gsub("cyx(16_18|18_21|21_23)", "cyxall", riclpm_syntax(bands = TRUE)))
  slope_syntax <- function() {
    s <- riclpm_syntax(bands = TRUE)
    load <- paste0(AGES - AGES[1], "*x", AGES, collapse = " + ")
    s <- paste(s,
      paste("Sx =~", load), paste("Sy =~", gsub("\\*x", "*y", load)),
      "Sx ~~ Sx", "Sy ~~ Sy", "Sx ~~ Sy + RIx + RIy", "Sy ~~ RIx + RIy",
      sprintf("Sx ~~ 0*wx%d", AGES[1]), sprintf("Sx ~~ 0*wy%d", AGES[1]),
      sprintf("Sy ~~ 0*wx%d", AGES[1]), sprintf("Sy ~~ 0*wy%d", AGES[1]), sep = "\n")
    s
  }
  fit_syn <- function(syn, w) tryCatch(lavaan(syn, data = w, missing = "ml", estimator = "MLR",
                         meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE),
                         error = function(e) { cat("  model failed:", conditionMessage(e), "\n"); NULL })
  J_tab <- list(); K_band <- list(); K_fit <- list()
  for (v in all_x) {
    cat("\n--- J/K:", lab(v), "---\n")
    w <- make_wide(v)
    if ("J" %in% RUN) {
      fits <- list(invariant = fit_syn(invariant_syntax(), w), bands = fit_syn(riclpm_syntax(bands = TRUE), w),
                   free = fit_syn(riclpm_syntax(bands = FALSE), w))
      ic <- sapply(fits, function(f) if (is.null(f)) c(NA, NA) else c(AIC(f), BIC(f)))
      J_tab[[v]] <- data.frame(variable = v, predictor = lab(v),
                               dAIC_bands_vs_invariant = ic[1, "bands"] - ic[1, "invariant"],
                               dBIC_bands_vs_invariant = ic[2, "bands"] - ic[2, "invariant"],
                               dAIC_free_vs_bands = ic[1, "free"] - ic[1, "bands"],
                               dBIC_free_vs_bands = ic[2, "free"] - ic[2, "bands"])
      if (!is.null(fits$invariant)) {
        pe <- parameterEstimates(fits$invariant)
        xy <- pe[pe$label == "cxyall", ][1, ]; yx <- pe[pe$label == "cyxall", ][1, ]
        J_tab[[v]] <- cbind(J_tab[[v]], inv_X_to_off = xy$est, inv_lo = xy$ci.lower, inv_hi = xy$ci.upper,
                            inv_off_to_X = yx$est, inv_lo_r = yx$ci.lower, inv_hi_r = yx$ci.upper)
      }
    }
    if ("K" %in% RUN) {
      fk <- fit_syn(slope_syntax(), w)
      K_band[[v]] <- band_rows(fk, v, "random_slopes")
      if (!is.null(fk)) K_fit[[v]] <- data.frame(variable = v, predictor = lab(v),
                                               converged = lavInspect(fk, "converged"),
                                               slope_var_X = parameterEstimates(fk)$est[parameterEstimates(fk)$lhs == "Sx" & parameterEstimates(fk)$rhs == "Sx"][1],
                                               slope_var_off = parameterEstimates(fk)$est[parameterEstimates(fk)$lhs == "Sy" & parameterEstimates(fk)$rhs == "Sy"][1])
    }
  }
  if ("J" %in% RUN) save_tab(do.call(rbind, J_tab), "J1_age_bands_vs_invariant",
     "Information criteria differences (negative = the first model fits better) and the age-invariant cross-lags")
  if ("K" %in% RUN) {
    save_tab(do.call(rbind, K_band), "K1_cross_lags_random_slopes", "Band model plus person-specific linear slopes for X and offending")
    save_tab(do.call(rbind, K_fit), "K2_random_slope_variances", "Slope variances (per year of age)")
  }
}

#===============================================================================
# L. Incapacitation check (added after the Draft 1 review). Employment and school
#    count community time only, and partner status is 0 while a youth lives in an
#    institution, so offending -> confinement -> fewer bonds next year is partly
#    built in. Here secure time enters BOTH equations (offending and the predictor).
#===============================================================================
if ("L" %in% RUN) {
  L_band <- list(); L_sx <- list()
  for (v in all_x) {
    cat("\n--- L:", lab(v), "---\n")
    w <- make_wide(v)
    fl <- fit_ri(w, bands = TRUE, street = "both")
    L_band[[v]] <- band_rows(fl, v, "secure_time_both_equations")
    if (!is.null(fl)) {
      pe <- parameterEstimates(fl); r <- pe[pe$label == "sx", ][1, ]
      L_sx[[v]] <- data.frame(variable = v, predictor = lab(v), secure_time_on_X = r$est, lo = r$ci.lower, hi = r$ci.upper)
    }
  }
  save_tab(do.call(rbind, L_band), "L1_cross_lags_secure_time_both",
           "Band model with time in secure settings in both the offending and the predictor equations")
  save_tab(do.call(rbind, L_sx), "L2_secure_time_on_predictor",
           "Same-age effect of a full recall period in secure settings on the predictor (bonds: 0-1 units; scales: SD)")
}

cat("\nDone:", format(Sys.time()), "\n")
sink()
