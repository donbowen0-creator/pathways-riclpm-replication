#===============================================================================
# 01_data_prep.R  --  Pathways to Desistance (ICPSR 29961)
#
# Builds the 8-wave annual person-year panel used for the desistance paper.
# Predictors are organised into four domains drawn from established desistance
# theory:
#
#   1. Social bonds            (age-graded informal social control; Sampson & Laub)
#   2. Cognitive shifts        (cognitive transformation / perceived costs; Giordano,
#                               Paternoster & Bushway, rational-choice perceptions)
#   3. Maturation              (psychosocial maturity; Steinberg, Cauffman, Monahan)
#   4. Childhood promotive factors (time-invariant, measured at baseline)
#
# No Criminogenic Embeddedness Model (CEM) measures or composites are built here.
# Every item is kept as its own variable -- nothing is averaged into an index --
# so each one can be reported with its own effect size.
#
# Output: panel_desistance.rds (one row per person-year, years 0-7)
#         panel_6month.rds    (one row per person-wave, waves 0-6, 6 months apart)
# Run this once; re-run only if you change the variable lists in Section 3.
#===============================================================================

# ---- 0. Folder that holds the DS0001 ... DS0011 sub-folders --------------------
data_dir <- "."   # <- CHANGE to the ICPSR_29961 folder if you run from elsewhere

# ---- 1. Helpers ----------------------------------------------------------------

# ICPSR stores categorical items as R factors with labels like
# "(01) Mostly A s". This turns them back into their numeric codes.
decode_factor_column <- function(x) {
  if (!is.factor(x)) return(x)
  codes <- suppressWarnings(as.numeric(sub("^\\((-?[0-9]+)\\).*$", "\\1", levels(x))))
  codes[as.integer(x)]
}
decode_all_factors <- function(df) {
  for (col in names(df)) if (is.factor(df[[col]])) df[[col]] <- decode_factor_column(df[[col]])
  df
}

# In the ICPSR R files, missing-data codes are already NA. This is a safety net
# in case any negative sentinel codes (-9, -99, -101, -201 ...) survive.
recode_missing <- function(x, thresh = -8) {
  x[!is.na(x) & x <= thresh] <- NA
  x
}

zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Grab a variable if the wave has it, otherwise return NA (with a warning)
get_var <- function(df, name) {
  if (name %in% names(df)) return(recode_missing(df[[name]]))
  warning("Missing variable: ", name)
  rep(NA_real_, nrow(df))
}

# ---- 2. Load all 11 waves -------------------------------------------------------

wave_prefixes <- c("S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","SA")

parts <- vector("list", 11)
for (i in 1:11) {
  fname <- file.path(data_dir, sprintf("DS%04d", i), sprintf("29961-%04d-Data.rda", i))
  if (!file.exists(fname)) fname <- file.path(data_dir, sprintf("29961-%04d-Data.rda", i))
  nm <- load(fname)
  parts[[i]] <- decode_all_factors(get(nm))
  rm(list = nm)
  cat("Loaded wave", i - 1, "--", nrow(parts[[i]]), "rows,", ncol(parts[[i]]), "cols\n")
}

# ---- 3. Variable lists (edit here to change the measures) ----------------------

# Cognitive shifts -- scale scores, kept in their natural direction
cognitive_stems <- c(
  EXPECT = "EXPECT",   # Expectations to have work, family, and stay out of trouble with the law
  MORDIS = "MORDIS",   # Moral disengagement (higher = MORE disengagement)
  SOCCST = "SOCCST",   # Perceived social costs of punishment
  PUNYOU = "PUNYOU",   # Perceived certainty of punishment for self (0-10)
  PERRWD = "PERRWD"    # Perceived personal rewards of crime (higher = MORE reward)
)

# Maturation -- psychosocial maturity scale scores
# (TEMPER is left out on purpose: it combines IMPULS and SUPAGG, so including it
#  would count the same items twice)
maturation_stems <- c(
  IMPULS = "IMPULS",   # WAI impulse control
  SUPAGG = "SUPAGG",   # WAI suppression of aggression
  CONSID = "CONSID",   # WAI consideration of others
  FUTURE = "FUTURE",   # Future orientation
  PEERPR = "PEERPR",   # Resistance to peer influence
  PSMI   = "PSMI"      # Psychosocial maturity inventory (responsibility), overall
)

# Check variables -- NOT part of the four domains. Used only by
# 03_robustness_checks.R (positive control and Loughran-style replication).
check_stems <- c(
  PRINFL = "PRINFL",   # Peer antisocial influence (Forney & Ward 2019 found a lagged effect)
  PRBEHV = "PRBEHV",   # Peer antisocial behavior
  PSMIID = "PSMIID",   # PSMI identity subscale (Forney & Ward's identity measure)
  STLRWD = "STLRWD",   # Social rewards of crime - stealing
  FGTRWD = "FGTRWD",   # Social rewards of crime - fighting
  ROBRWD = "ROBRWD"    # Social rewards of crime - robbery
)

# Social bonds -- built from the follow-up calendars in Section 4 (not asked at
# baseline in this form, so they are NA for Year 0):
#   WORK_PROP   share of the recall period employed in the community
#   SCHOOL_PROP share of the recall period enrolled in a community school
#   PARTNER     lived with a significant other (spouse, partner, parent of child)
#   PARENT      has at least one living child

# ---- 4. Extract per-wave variables ---------------------------------------------

wave_list <- vector("list", 11)
for (wi in 0:10) {
  p  <- wave_prefixes[wi + 1]
  df <- parts[[wi + 1]]
  out <- data.frame(CASEID = df$CASEID, wave = wi)

  for (s in names(cognitive_stems))  out[[s]] <- get_var(df, paste0(p, cognitive_stems[[s]]))
  for (s in names(maturation_stems)) out[[s]] <- get_var(df, paste0(p, maturation_stems[[s]]))
  for (s in names(check_stems))      out[[s]] <- get_var(df, paste0(p, check_stems[[s]]))
  out$SOCRWD <- rowMeans(out[, c("STLRWD", "FGTRWD", "ROBRWD")], na.rm = TRUE)   # social rewards, mean of 3
  out$SOCRWD[is.nan(out$SOCRWD)] <- NA

  # Outcome: variety of self-reported offending in the recall period
  # (baseline uses the 6-month-before-baseline version)
  out$OFFVAR <- get_var(df, paste0(p, if (wi == 0) "SRO6" else "SROPRV"))

  if (wi >= 1) {
    nmonths <- get_var(df, paste0(p, "NMONTHS"))          # length of recall period
    interviewed <- !is.na(nmonths)                         # NA = missed interview

    # Work: weeks worked in the community / weeks in the recall period
    wk <- get_var(df, paste0(p, "JOBCAL_NWEEKSCU"))
    out$WORK_PROP <- pmin(wk / (nmonths * 4.345), 1)

    # School: months in a community school / months in the recall period.
    # ICPSR skips this item when the youth was not enrolled; in the R file the
    # skip shows up as NA, so for interviewed youth NA is set to 0.
    sch <- get_var(df, paste0(p, "SCHCAL_MAINFREECOUNT"))
    sch[is.na(sch) & interviewed] <- 0
    out$SCHOOL_PROP <- pmin(sch / nmonths, 1)

    # Partner: the household question is skipped when the youth's main location
    # was an institution, so for interviewed youth NA is set to 0 (no partner
    # in the household).
    prt <- get_var(df, paste0(p, "HOMEQ_SIGOTH"))
    prt[is.na(prt) & interviewed] <- 0
    out$PARTNER <- prt

    # Parenthood: any living child
    out$PARENT <- as.numeric(get_var(df, paste0(p, "CHLD_COUNT")) > 0)

    # Time in locked settings during the recall period (exposure control)
    out$PROP_SECURE <- get_var(df, paste0(p, "PROPTIMESECURESETTINGS"))
  } else {
    out$WORK_PROP <- out$SCHOOL_PROP <- out$PARTNER <- out$PARENT <- out$PROP_SECURE <- NA
  }

  wave_list[[wi + 1]] <- out
}

# ---- 5. Childhood promotive factors and baseline covariates -------------------

bl <- parts[[1]]
grades <- recode_missing(bl$S0DEM29)
grades[!is.na(grades) & grades > 8] <- NA          # 9-11 = invalid / no grades / other

cov <- data.frame(
  CASEID = bl$CASEID,

  # Childhood promotive factors (all coded so that higher = more promotive)
  CP_NO_EARLY_PROBLEMS = 5 - recode_missing(bl$S0NEARPRO),        # 0-5; problems BEFORE age 11, reversed
  CP_TWO_BIO_PARENTS   = as.numeric(recode_missing(bl$S0NBPHOM) == 2),
  CP_PARENT_EDUC       = recode_missing(bl$S0PAEDUC),             # parents' education level
  CP_SCHOOL_GRADES     = 9 - grades,                              # 1 = mostly below Ds ... 8 = mostly As
  CP_IQ                = recode_missing(bl$S0FSIQ),               # WASI full-scale IQ
  CP_NO_PARENT_ARREST  = as.numeric(recode_missing(bl$S0ARRBTH) == 0),

  # Baseline covariates
  BL_AGE             = bl$S0AGE,
  BL_SEX             = bl$S0SGEND,          # 1 = Male, 2 = Female
  BL_SITE            = bl$S0SITE,           # 1 = Philadelphia, 2 = Phoenix
  BL_RACE            = bl$S0ETHN_R,         # 1 White, 2 Black, 3 Hispanic, 4 Other
  BL_PRIOR_PETITIONS = recode_missing(bl$PRI_B4_BL_EVER)
)

# ---- 6. Build the annual panel (Year 0-7) -------------------------------------
# Waves 1, 3 and 5 (the 6-month interviews) are dropped so that each row is one
# year. Note: Years 1-3 come from 6-month recall interviews, Years 4-7 from
# 12-month recall interviews; WORK_PROP and SCHOOL_PROP are proportions so they
# are comparable across the two.

annual_map <- c(`0` = 0, `2` = 1, `4` = 2, `6` = 3, `7` = 4, `8` = 5, `9` = 6, `10` = 7)

panel <- do.call(rbind, lapply(names(annual_map), function(w) {
  d <- wave_list[[as.integer(w) + 1]]
  d$year <- annual_map[[w]]
  d
}))
panel <- merge(panel, cov, by = "CASEID", all.x = TRUE)
panel <- panel[order(panel$CASEID, panel$year), ]
panel$YRS <- panel$year                      # years since baseline (used for age)

# Second panel at 6-month spacing (waves 0-6: baseline to 36 months), used to
# test whether a shorter lag reveals within-person effects a 1-year lag misses.
panel6 <- do.call(rbind, lapply(0:6, function(w) {
  d <- wave_list[[w + 1]]
  d$year <- w                                # time index = wave number
  d
}))
panel6 <- merge(panel6, cov, by = "CASEID", all.x = TRUE)
panel6 <- panel6[order(panel6$CASEID, panel6$year), ]
panel6$YRS <- panel6$year / 2                # years since baseline

# ---- 7. Standardised versions of the continuous scales ------------------------
# z-scored across all pooled person-years, so effects read as "per 1 SD".
# Binary bonds (PARTNER, PARENT) and the proportions are left on their own scale.

for (s in c(names(cognitive_stems), names(maturation_stems), "PRINFL", "PRBEHV", "PSMIID", "SOCRWD")) {
  panel[[paste0(s, "_z")]] <- zscore(panel[[s]])
}
for (s in c("CP_NO_EARLY_PROBLEMS", "CP_PARENT_EDUC", "CP_SCHOOL_GRADES", "CP_IQ")) {
  panel[[paste0(s, "_z")]] <- zscore(panel[[s]])
}

# ---- 8. Save ------------------------------------------------------------------

saveRDS(panel, file.path(data_dir, "panel_desistance.rds"))
cat("\nSaved panel_desistance.rds --", nrow(panel), "person-years,",
    length(unique(panel$CASEID)), "people,", ncol(panel), "columns\n")
saveRDS(panel6, file.path(data_dir, "panel_6month.rds"))
cat("Saved panel_6month.rds --", nrow(panel6), "person-waves (baseline to 36 months, 6-month spacing)\n")

cat("\nNon-missing values by year:\n")
print(aggregate(cbind(OFFVAR, WORK_PROP, PARTNER, EXPECT, IMPULS) ~ year, data = panel,
                FUN = function(x) sum(!is.na(x)), na.action = na.pass))
