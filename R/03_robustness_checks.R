#===============================================================================
# 03_robustness_checks.R  --  Pathways to Desistance (ICPSR 29961)
# Checks requested after the literature review (2026-10-02)
#
# A. POSITIVE CONTROL. Peer antisocial influence and PSMI identity are added to
#    the main lagged hybrid model. Forney & Ward (2019) found both predict later
#    offending within-person. If our model picks them up, the nulls for the
#    other predictors are not just a design that can't find lagged effects.
# B. WITHIN-PERSON VARIATION. How much each predictor actually changes inside a
#    person (Lanfear & Oliveira 2026: too little change makes lagged effects
#    hard to find).
# C. FIXED-EFFECTS POISSON on the number of offense types, same period vs
#    previous period, on observed data (Loughran et al. 2016; Thomas et al. 2023).
# D. RANDOM-INTERCEPT CROSS-LAGGED PANEL MODELS by age, 16-23
#    (Chen & Sutton 2025; Forney & Ward 2019).
#
# Requires panel_desistance.rds from 01_data_prep.R (2026-10-02 version, which
# adds PRINFL, PRBEHV, PSMIID and SOCRWD).
# Packages: install.packages(c("mice", "sandwich", "lmtest", "fixest", "lavaan"))
# Output  : output_checks/ (CSV tables, run_log.txt, PNG chart)
# Run     : Rscript 03_robustness_checks.R        (M = 20 imputations)
#           Rscript 03_robustness_checks.R 5      (quick look, M = 5)
#===============================================================================

library(sandwich)
library(lmtest)
library(mice)
library(splines)
library(fixest)
library(lavaan)

data_dir <- "."
out_dir  <- file.path(data_dir, "output_checks")
dir.create(out_dir, showWarnings = FALSE)
M <- 20
args <- commandArgs(trailingOnly = TRUE); if (length(args) > 0) M <- as.integer(args[1])
set.seed(20261002)

panel <- readRDS(file.path(data_dir, "panel_desistance.rds"))
stopifnot(all(c("PRINFL", "PSMIID", "SOCRWD") %in% names(panel)))   # needs the updated 01_data_prep.R
if (is.null(panel$YRS)) panel$YRS <- panel$year
panel   <- panel[order(panel$CASEID, panel$year), ]
T_MAX   <- max(panel$year)
yrs_map <- tapply(panel$YRS, panel$year, function(x) x[1])

#-------------------------------------------------------------------------------
# 0. Variable sets (same as 02_run_all_analyses.R, plus the two check variables)
#-------------------------------------------------------------------------------
bonds      <- c("WORK_PROP", "SCHOOL_PROP", "PARTNER", "PARENT")
cognitive  <- c("EXPECT", "MORDIS", "SOCCST", "PUNYOU", "PERRWD")
maturation <- c("IMPULS", "SUPAGG", "CONSID", "FUTURE", "PEERPR", "PSMI")
checks     <- c("PRINFL", "PSMIID")        # peer antisocial influence; identity
childhood  <- c("CP_NO_EARLY_PROBLEMS", "CP_TWO_BIO_PARENTS", "CP_PARENT_EDUC",
                "CP_SCHOOL_GRADES", "CP_IQ", "CP_NO_PARENT_ARREST")
baseline   <- c("BL_AGE", "BL_SEX", "BL_SITE", "BL_RACE", "BL_PRIOR_PETITIONS")
tv_all     <- c("OFFVAR", "PROP_SECURE", bonds, cognitive, maturation, checks)
followup_only <- c("PROP_SECURE", bonds)
scales     <- c(cognitive, maturation, checks, "CP_NO_EARLY_PROBLEMS", "CP_PARENT_EDUC",
                "CP_SCHOOL_GRADES", "CP_IQ")
zc <- sapply(scales, function(v) c(m = mean(panel[[v]], na.rm = TRUE), s = sd(panel[[v]], na.rm = TRUE)))

lab_map <- c(WORK_PROP = "Employment (share of period)", SCHOOL_PROP = "School (share of period)",
             PARTNER = "Living with partner", PARENT = "Parent",
             EXPECT = "Expectations", MORDIS = "Moral disengagement",
             SOCCST = "Social costs of punishment", PUNYOU = "Certainty of punishment",
             PERRWD = "Personal rewards of crime", SOCRWD = "Social rewards of crime",
             IMPULS = "Impulse control", SUPAGG = "Suppression of aggression",
             CONSID = "Consideration of others", FUTURE = "Future orientation",
             PEERPR = "Resistance to peer influence", PSMI = "Psychosocial maturity",
             PRINFL = "Peer antisocial influence", PSMIID = "Identity (PSMI)")
lab <- function(v) unname(lab_map[sub("_z$", "", v)])

save_tab <- function(t, name, note = NULL) {
  num <- sapply(t, is.numeric); t[num] <- lapply(t[num], round, 3)
  write.csv(t, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat("\n=====", name, "=====\n"); if (!is.null(note)) cat(note, "\n")
  print(t, row.names = FALSE)
}

sink(file.path(out_dir, "run_log.txt"), split = TRUE)
cat("Pathways robustness checks -- R run,", format(Sys.time()), " M =", M, "\n")

#===============================================================================
# B. Within-person variation (observed data, follow-up Years 1-7)
#===============================================================================
pv <- panel[panel$year >= 1, ]
wvar <- do.call(rbind, lapply(c(bonds, cognitive, maturation, checks), function(v) {
  x <- pv[[v]]; id <- pv$CASEID
  pm  <- ave(x, id, FUN = function(z) mean(z, na.rm = TRUE))
  dev <- x - pm
  chg <- ave(x, id, FUN = function(z) c(NA, diff(z)))           # year-to-year change
  ever <- tapply(x, id, function(z) { z <- z[!is.na(z)]; length(z) >= 2 && diff(range(z)) > 0 })
  sd_tot <- sd(x, na.rm = TRUE)
  data.frame(variable = v, label = lab(v),
             sd_total   = sd_tot,
             sd_between = sd(tapply(x, id, mean, na.rm = TRUE), na.rm = TRUE),
             sd_within  = sd(dev, na.rm = TRUE),
             within_share_of_variance = var(dev, na.rm = TRUE) / var(x, na.rm = TRUE),
             within_sd_in_pooled_sd   = sd(dev, na.rm = TRUE) / sd_tot,
             median_yearly_change_in_pooled_sd = median(abs(chg), na.rm = TRUE) / sd_tot,
             pct_people_who_change = 100 * mean(ever, na.rm = TRUE))
}))
save_tab(wvar, "B_within_person_variation",
         paste("Observed data, Years 1-7. within_sd_in_pooled_sd: a typical within-person swing,",
               "in the per-1-SD units the ORs use. Binary bonds: pct_people_who_change = share",
               "who ever switch status."))

#===============================================================================
# A. Positive control (multiple imputation, same design as the main model)
#===============================================================================
ids  <- sort(unique(panel$CASEID))
wide <- data.frame(CASEID = ids)
for (v in tv_all) for (y in 0:T_MAX) {
  if (v %in% followup_only && y == 0) next
  x <- panel[panel$year == y, c("CASEID", v)]
  wide[[paste0(v, "_Y", y)]] <- x[[v]][match(ids, x$CASEID)]
}
b0 <- panel[panel$year == 0, c("CASEID", childhood, baseline)]
wide <- cbind(wide, b0[match(ids, b0$CASEID), -1])
wide$FEMALE  <- as.numeric(wide$BL_SEX == 2)
wide$PHOENIX <- as.numeric(wide$BL_SITE == 2)
wide$BLACK   <- as.numeric(wide$BL_RACE == 2)
wide$HISP    <- as.numeric(wide$BL_RACE == 3)
wide$OTHER   <- as.numeric(wide$BL_RACE == 4)
impcols <- setdiff(names(wide), c("CASEID", "BL_SEX", "BL_SITE", "BL_RACE"))

pred <- quickpred(wide[, impcols], mincor = 0.20,
                  include = c("FEMALE", "PHOENIX", "BLACK", "HISP", "OTHER", "BL_AGE"))
t0  <- Sys.time()
imp <- mice(wide[, impcols], m = M, method = "pmm", predictorMatrix = pred, maxit = 10, printFlag = FALSE)
cat("\nImputation done: M =", M, "in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")

obs_out <- !is.na(panel$OFFVAR); names(obs_out) <- paste(panel$CASEID, panel$year)
tvz    <- c(bonds, paste0(c(cognitive, maturation, checks), "_z"))
lagged <- paste0(tvz, "_L1")

build_long <- function(w) {
  long <- do.call(rbind, lapply(0:T_MAX, function(y) {
    d <- data.frame(CASEID = ids, year = y)
    for (v in tv_all) {
      col <- paste0(v, "_Y", y)
      d[[v]] <- if (col %in% names(w)) w[[col]] else NA_real_
    }
    d
  }))
  long <- merge(long, cbind(CASEID = ids, wide[, c("BL_SEX", "BL_SITE", "BL_RACE")],
                            w[, c(childhood, "BL_AGE", "BL_PRIOR_PETITIONS")]), by = "CASEID")
  long <- long[order(long$CASEID, long$year), ]
  long$outcome_observed <- obs_out[paste(long$CASEID, long$year)]
  long$outcome_observed[is.na(long$outcome_observed)] <- FALSE
  long$ANY_OFFEND   <- as.numeric(long$OFFVAR > 0)
  long$AGE_T        <- long$BL_AGE + yrs_map[as.character(long$year)]
  long$SEX_FEMALE   <- as.numeric(long$BL_SEX == 2)
  long$SITE_PHOENIX <- as.numeric(long$BL_SITE == 2)
  long$RACE_BLACK   <- as.numeric(long$BL_RACE == 2)
  long$RACE_HISP    <- as.numeric(long$BL_RACE == 3)
  long$RACE_OTHER   <- as.numeric(long$BL_RACE == 4)
  for (v in scales) long[[paste0(v, "_z")]] <- (long[[v]] - zc["m", v]) / zc["s", v]
  lag1 <- function(x) ave(x, long$CASEID, FUN = function(z) c(NA, head(z, -1)))
  long$OFFEND_L1      <- lag1(long$ANY_OFFEND)
  long$PROP_SECURE_L1 <- lag1(long$PROP_SECURE)
  for (v in tvz) long[[paste0(v, "_L1")]] <- lag1(long[[v]])
  long
}
build_analytic <- function(long) {
  d <- long[long$year >= 2 & long$outcome_observed, ]
  for (v in lagged) { pm <- ave(d[[v]], d$CASEID, FUN = mean)
    d[[paste0("BW_", v)]] <- pm; d[[paste0("WI_", v)]] <- d[[v]] - pm }
  for (v in tvz) { pm <- ave(d[[v]], d$CASEID, FUN = mean)
    d[[paste0("BW_", v, "_T0")]] <- pm; d[[paste0("WI_", v, "_T0")]] <- d[[v]] - pm }
  d
}
d_list <- lapply(1:M, function(i) build_analytic(build_long(complete(imp, i))))
cat("Analytic person-years:", nrow(d_list[[1]]), "from", length(unique(d_list[[1]]$CASEID)), "people\n")

hyb <- function(vars) c(paste0("WI_", vars, "_L1"), paste0("BW_", vars, "_L1"))
T0  <- function(vars) c(paste0("WI_", vars, "_T0"), paste0("BW_", vars, "_T0"))
main_tv  <- c(bonds, paste0(c(cognitive, maturation), "_z"))
check_tv <- paste0(checks, "_z")
T_child  <- c("CP_NO_EARLY_PROBLEMS_z", "CP_TWO_BIO_PARENTS", "CP_PARENT_EDUC_z",
              "CP_SCHOOL_GRADES_z", "CP_IQ_z", "CP_NO_PARENT_ARREST")
controls <- c("ns(AGE_T, df = 3)", "SEX_FEMALE", "SITE_PHOENIX", "RACE_BLACK", "RACE_HISP",
              "RACE_OTHER", "BL_PRIOR_PETITIONS", "OFFEND_L1", "PROP_SECURE", "PROP_SECURE_L1",
              "factor(year)")
fml <- function(y, terms) as.formula(paste(y, "~", paste(terms, collapse = " + ")))
fit_one <- function(d, terms) {
  m <- glm(fml("ANY_OFFEND", terms), data = d, family = binomial())
  list(b = coef(m), V = vcovCL(m, cluster = d$CASEID, type = "HC1"))
}
run_pooled <- function(terms) {
  fits <- lapply(d_list, fit_one, terms = terms)
  B <- sapply(fits, `[[`, "b"); W <- Reduce(`+`, lapply(fits, `[[`, "V")) / length(fits)
  list(b = rowMeans(B), V = W + (1 + 1 / length(fits)) * cov(t(B)))
}
or_ci <- function(pool, term) {
  b <- pool$b[term]; se <- sqrt(pool$V[term, term])
  c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}

P_prev_full  <- run_pooled(c(hyb(main_tv), hyb(check_tv), T_child, controls))   # main model + checks
P_same_full  <- run_pooled(c(T0(c(main_tv, check_tv)), T_child, controls))
P_prev_alone <- run_pooled(c(hyb(check_tv), T_child, controls))                  # checks on their own
P_same_alone <- run_pooled(c(T0(check_tv), T_child, controls))
P_prev_main  <- run_pooled(c(hyb(main_tv), T_child, controls))                   # main model as in 02

row_or <- function(pred, model, timing, x)
  data.frame(predictor = pred, model = model, timing = timing,
             OR_within = unname(x[1]), lo = unname(x[2]), hi = unname(x[3]), p = unname(x[4]))
posc <- do.call(rbind, lapply(check_tv, function(v) rbind(
  row_or(lab(v), "Alone (+ childhood, controls)", "Same period",   or_ci(P_same_alone, paste0("WI_", v, "_T0"))),
  row_or(lab(v), "Alone (+ childhood, controls)", "Previous year", or_ci(P_prev_alone, paste0("WI_", v, "_L1"))),
  row_or(lab(v), "With all four domains",         "Same period",   or_ci(P_same_full,  paste0("WI_", v, "_T0"))),
  row_or(lab(v), "With all four domains",         "Previous year", or_ci(P_prev_full,  paste0("WI_", v, "_L1"))))))
save_tab(posc, "A1_positive_control",
         "Within-person ORs per 1 SD. Previous year = measured a year before the offending period; controls prior offending.")

# Did adding the checks change the main within-person estimates?
shift <- do.call(rbind, lapply(main_tv, function(v) {
  a <- or_ci(P_prev_main, paste0("WI_", v, "_L1")); b <- or_ci(P_prev_full, paste0("WI_", v, "_L1"))
  data.frame(predictor = lab(v), OR_main_model = a[1], lo = a[2], hi = a[3],
             OR_with_checks = b[1], lo_checks = b[2], hi_checks = b[3])
}))
save_tab(shift, "A2_main_estimates_with_checks_added",
         "Previous-year within-person ORs, main model vs main model + peer influence and identity")

# Per typical within-person change (rescales the per-1-SD ORs by the within SD)
wsd <- setNames(wvar$within_sd_in_pooled_sd, wvar$variable)
resc <- do.call(rbind, lapply(paste0(c(cognitive, maturation, checks), "_z"), function(v) {
  s <- or_ci(P_same_full, paste0("WI_", v, "_T0")); p <- or_ci(P_prev_full, paste0("WI_", v, "_L1"))
  k <- wsd[sub("_z$", "", v)]
  data.frame(predictor = lab(v), within_sd = k,
             same_OR_per_SD = s[1], prev_OR_per_SD = p[1],
             same_OR_per_within_SD = s[1]^k, same_lo = s[2]^k, same_hi = s[3]^k,
             prev_OR_per_within_SD = p[1]^k, prev_lo = p[2]^k, prev_hi = p[3]^k)
}))
save_tab(resc, "B2_ORs_per_typical_within_change",
         "Same-period and previous-year within-person ORs rescaled to a typical (1 within-SD) change")

# Is there enough lagged within-person variation in the analytic rows?
d1 <- d_list[[1]]
lagvar <- data.frame(predictor = sapply(tvz, lab),
                     sd_within_same = sapply(tvz, function(v) sd(d1[[paste0("WI_", v, "_T0")]])),
                     sd_within_prev = sapply(tvz, function(v) sd(d1[[paste0("WI_", v, "_L1")]])),
                     row.names = NULL)
save_tab(lagvar, "B3_within_sd_same_vs_previous_rows",
         "Within-person SD of each predictor in the analytic rows (imputation 1): same-period vs previous-year versions")

#===============================================================================
# C. Fixed-effects Poisson, same period vs previous period (observed data)
#===============================================================================
# Outcome: number of offense types in the recall period (variety proportion x 22
# items, rounded; the proportion is stored to 2 decimals). Person and year fixed
# effects; SEs clustered on person. No prior-offending term (person fixed effects
# with a lagged outcome are biased). People who never offend drop out.
ob <- panel
ob$OFFCNT <- round(ob$OFFVAR * 22)
fe_vars <- c(main_tv, check_tv, "SOCRWD_z")
lag1o <- function(x) ave(x, ob$CASEID, FUN = function(z) c(NA, head(z, -1)))
for (v in c(fe_vars, "PROP_SECURE")) ob[[paste0(v, "_L1")]] <- lag1o(ob[[v]])
ob <- ob[ob$year >= 2, ]

fe_pair <- function(vars, name) {
  same_t <- c(vars, "PROP_SECURE"); prev_t <- c(paste0(vars, "_L1"), "PROP_SECURE", "PROP_SECURE_L1")
  dd <- ob[complete.cases(ob[, c("OFFCNT", same_t, prev_t)]), ]          # same rows for both
  f_s <- as.formula(paste("OFFCNT ~", paste(same_t, collapse = " + "), "| CASEID + year"))
  f_p <- as.formula(paste("OFFCNT ~", paste(prev_t, collapse = " + "), "| CASEID + year"))
  ms <- fepois(f_s, data = dd, cluster = ~CASEID); mp <- fepois(f_p, data = dd, cluster = ~CASEID)
  cat("\n", name, ": ", nobs(ms), " person-years, ", length(unique(dd$CASEID[obs(ms)])),
      " people with any offending used\n", sep = "")
  get <- function(m, term) { b <- coef(m)[term]; se <- se(m)[term]
    c(IRR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se))) }
  t <- do.call(rbind, lapply(vars, function(v) {
    s <- get(ms, v); p <- get(mp, paste0(v, "_L1"))
    data.frame(predictor = lab(v), same_IRR = s[1], same_lo = s[2], same_hi = s[3], same_p = s[4],
               prev_IRR = p[1], prev_lo = p[2], prev_hi = p[3], prev_p = p[4])
  }))
  save_tab(t, name, "Fixed-effects Poisson IRRs (person + year FE). IRR < 1 = fewer offense types.")
  t
}
fe_lough <- fe_pair(c("PUNYOU_z", "PERRWD_z", "SOCRWD_z", "SOCCST_z"), "C1_fepois_loughran_set")
fe_full  <- fe_pair(c(main_tv, check_tv), "C2_fepois_full_battery")

#===============================================================================
# D. Random-intercept cross-lagged panel models by age (16-23)
#===============================================================================
# Each person's interviews are placed at their age (baseline age + years), so
# time points are ages 16 ... 23 and lags are one year. X = the predictor (per
# pooled SD); Y = offending variety (per pooled SD). FIML for missing ages, MLR.
# Cross-lag X(a) -> Y(a+1) = within-person change in X at age a predicting
# within-person offending in the following year. Bands follow Chen & Sutton (2025).
AGES  <- 16:23
pa    <- panel; pa$AGE <- pa$BL_AGE + pa$YRS
pa$Y_z <- (pa$OFFVAR - mean(pa$OFFVAR, na.rm = TRUE)) / sd(pa$OFFVAR, na.rm = TRUE)
band_of <- function(a) ifelse(a <= 18, "16-18", ifelse(a <= 21, "18-21", "21-23"))   # by the lag's END age

riclpm_syntax <- function(bands = FALSE) {
  s <- c(paste("RIx =~", paste0("1*x", AGES, collapse = " + ")),
         paste("RIy =~", paste0("1*y", AGES, collapse = " + ")))
  for (a in AGES) s <- c(s, sprintf("wx%d =~ 1*x%d", a, a), sprintf("wy%d =~ 1*y%d", a, a),
                         sprintf("x%d ~~ 0*x%d", a, a), sprintf("y%d ~~ 0*y%d", a, a))
  for (i in 2:length(AGES)) {
    a <- AGES[i]; p <- AGES[i - 1]
    tag <- if (bands) gsub("-", "_", band_of(a)) else a
    s <- c(s, sprintf("wx%d ~ arx%d*wx%d + cyx%s*wy%d", a, a, p, tag, p),
              sprintf("wy%d ~ ary%d*wy%d + cxy%s*wx%d", a, a, p, tag, p))
  }
  s <- c(s, sprintf("wx%d ~~ cov%d*wy%d", AGES, AGES, AGES),
         sprintf("wx%d ~~ wx%d", AGES, AGES), sprintf("wy%d ~~ wy%d", AGES, AGES),
         "RIx ~~ RIx", "RIy ~~ RIy", "RIx ~~ RIy",
         sprintf("RIx ~~ 0*wx%d", AGES[1]), sprintf("RIx ~~ 0*wy%d", AGES[1]),
         sprintf("RIy ~~ 0*wx%d", AGES[1]), sprintf("RIy ~~ 0*wy%d", AGES[1]))
  paste(s, collapse = "\n")
}

ri_vars <- c("MORDIS", "PERRWD", "PUNYOU", "EXPECT", "IMPULS", "SUPAGG", "FUTURE", "PEERPR",
             "PSMI", "CONSID", "SOCCST", "PRINFL", "PSMIID")
ri_lag <- list(); ri_band <- list(); ri_fit <- list()
for (v in ri_vars) {
  pa$X_z <- (pa[[v]] - zc["m", v]) / zc["s", v]
  w <- data.frame(CASEID = ids)
  for (a in AGES) {
    r <- pa[pa$AGE == a, ]
    w[[paste0("x", a)]] <- r$X_z[match(ids, r$CASEID)]
    w[[paste0("y", a)]] <- r$Y_z[match(ids, r$CASEID)]
  }
  w <- w[rowSums(!is.na(w[, -1])) > 0, ]
  fit_free <- lavaan(riclpm_syntax(FALSE), data = w, missing = "ml", estimator = "MLR",
                     meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE)
  fit_band <- lavaan(riclpm_syntax(TRUE), data = w, missing = "ml", estimator = "MLR",
                     meanstructure = TRUE, int.ov.free = TRUE, int.lv.free = FALSE)
  fm <- fitMeasures(fit_free, c("cfi.robust", "rmsea.robust", "srmr"))
  ri_fit[[v]] <- data.frame(predictor = lab(v), n = lavInspect(fit_free, "nobs"),
                            converged = lavInspect(fit_free, "converged"),
                            CFI = fm[1], RMSEA = fm[2], SRMR = fm[3])
  ss <- standardizedSolution(fit_free)
  same <- ss[ss$op == "~~" & grepl("^wx", ss$lhs) & grepl("^wy", ss$rhs), ]
  ri_lag[[v]] <- do.call(rbind, lapply(AGES[-1], function(a) {
    xy <- ss[ss$lhs == paste0("wy", a) & ss$rhs == paste0("wx", a - 1) & ss$op == "~", ]
    yx <- ss[ss$lhs == paste0("wx", a) & ss$rhs == paste0("wy", a - 1) & ss$op == "~", ]
    cw <- same[same$lhs == paste0("wx", a), ]
    sx <- ss[ss$lhs == paste0("wx", a) & ss$rhs == paste0("wx", a - 1) & ss$op == "~", ]
    sy <- ss[ss$lhs == paste0("wy", a) & ss$rhs == paste0("wy", a - 1) & ss$op == "~", ]
    data.frame(predictor = lab(v), lag = paste0(a - 1, "->", a),
               stability_X = sx$est.std, stability_offending = sy$est.std,
               same_age_r = cw$est.std, same_age_p = cw$pvalue,
               X_to_offending = xy$est.std, X_to_off_p = xy$pvalue,
               offending_to_X = yx$est.std, off_to_X_p = yx$pvalue)
  }))
  pe <- parameterEstimates(fit_band)
  ri_band[[v]] <- do.call(rbind, lapply(c("16_18", "18_21", "21_23"), function(b) {
    xy <- pe[pe$label == paste0("cxy", b), ][1, ]; yx <- pe[pe$label == paste0("cyx", b), ][1, ]
    data.frame(predictor = lab(v), band = sub("_", "-", b),
               X_to_offending = xy$est, lo = xy$ci.lower, hi = xy$ci.upper, p = xy$pvalue,
               offending_to_X = yx$est, lo_r = yx$ci.lower, hi_r = yx$ci.upper, p_r = yx$pvalue)
  }))
}
save_tab(do.call(rbind, ri_fit), "D0_riclpm_fit", "RI-CLPM fit (free cross-lags), ages 16-23")
save_tab(do.call(rbind, ri_band), "D1_riclpm_by_age_band",
         paste("Within-person cross-lagged effects, constrained equal within age bands (by the lag's end age).",
               "Units: SD of offending variety per SD of X (and the reverse)."))
ri_lag_tab <- do.call(rbind, ri_lag)
save_tab(ri_lag_tab, "D2_riclpm_each_lag",
         "Free model, standardized: same-age within-person correlation, X -> next-year offending, offending -> next-year X")

# Summary: average same-age within-person correlation vs average cross-lag (standardized)
ri_sum <- do.call(rbind, lapply(split(ri_lag_tab, ri_lag_tab$predictor), function(t)
  data.frame(predictor = t$predictor[1], mean_stability_X = mean(t$stability_X),
             mean_same_age_r = mean(t$same_age_r),
             same_age_sig = sum(t$same_age_p < .05), mean_X_to_off = mean(t$X_to_offending),
             X_to_off_sig = sum(t$X_to_off_p < .05), mean_off_to_X = mean(t$offending_to_X), mean_stability_offending = mean(t$stability_offending),
             off_to_X_sig = sum(t$off_to_X_p < .05), lags = nrow(t))))
save_tab(ri_sum[order(-abs(ri_sum$mean_same_age_r)), ], "D3_riclpm_summary",
         "Averages over the 7 one-year lags; *_sig = number of lags with p < .05")

#===============================================================================
# Chart: same-age within-person correlation vs cross-lag to next-year offending (RI-CLPM)
#===============================================================================
local({
  BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#0b0b0b"; MUTED <- "#52514e"
  t <- ri_sum[order(ri_sum$mean_same_age_r), ]; k <- nrow(t); y <- seq_len(k)
  png(file.path(out_dir, "fig_riclpm_same_age_vs_cross_lag.png"), width = 1250, height = 240 + 58 * k,
      res = 120, bg = "#fcfcfb")
  par(mar = c(4.2, 15.5, 5.6, 1.5), fg = "#cfcdc8", col.axis = MUTED)
  xr <- range(c(t$mean_same_age_r, t$mean_X_to_off, 0)) + c(-0.03, 0.03)
  plot(NA, xlim = xr, ylim = c(0.5, k + 0.5), yaxt = "n", xlab = "", ylab = "", bty = "l")
  abline(v = 0, col = MUTED, lty = 2)
  axis(2, at = y, labels = t$predictor, las = 1, cex.axis = 0.85, tick = FALSE, col.axis = INK)
  segments(t$mean_X_to_off, y, t$mean_same_age_r, y, col = "#cfcdc8", lwd = 2)
  points(t$mean_same_age_r, y, pch = 19, col = ORANGE, cex = 1.3)
  points(t$mean_X_to_off, y, pch = 19, col = BLUE, cex = 1.3)
  mtext("Standardized within-person estimate (average over ages 16-23)", side = 1, line = 2.6, cex = 0.8, col = MUTED)
  mtext("RI-CLPM: same-age association vs prediction of next year's offending", side = 3, line = 3.8,
        adj = 0, at = par("usr")[1], cex = 1.05, font = 2, col = INK)
  legend(x = par("usr")[1], y = k + 0.5 + k * 0.12, xpd = NA, bty = "n", cex = 0.85, horiz = TRUE,
         legend = c("Same age (within-person correlation)", "Predicts next year's offending"),
         col = c(ORANGE, BLUE), pch = 19, text.col = INK)
  dev.off()
})

cat("\nDone:", format(Sys.time()), "\n")
sink()
