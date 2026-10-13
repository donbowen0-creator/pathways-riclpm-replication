#===============================================================================
# 06_simulation_common_state.R  --  added 2026-10-03 after the novelty check
#
# Sorjonen et al. (2026, PLOS ONE) show that RI-CLPM cross-lags can be spurious
# when both variables share a slow-moving, auto-correlated "state" (for example
# a stretch of bad circumstances that raises both antisocial attitudes and
# offending for a few years). They recommend checking with person-mean-centred
# scores. This script asks how large a spurious cross-lag each estimator gives
# under that process, with NO causal effect in either direction and NO
# same-period causal link, when the shared state is tuned to Pathways:
#   - same-age within-person correlation about 0.15 (all of it from the shared state)
#   - within-person stability of X and Y about 0.2 (as in E3_summary_by_band.csv)
#
#   C(t)  = rho * C(t-1) + noise                         shared state, persistence rho
#   x(t)  = stable_x + sqrt(s) * C(t) + sqrt(1-s) * ux(t)  ux own AR(1), phi_own
#   y(t)  = stable_y + sqrt(s) * C(t) + sqrt(1-s) * uy(t)  uy own AR(1), phi_own
#   share = 0.15 or 0.30 of within-person variance from the shared state (equal in X and
#   offending, plus two conditions where it differs: 0.30 vs 0.10 and 0.10 vs 0.30);
#   phi_own chosen so total within-person stability is about 0.20 where possible.
#   Compare each scenario's simulated same-age correlation and stability with the
#   observed values (same-age r about 0.14-0.36, stability about 0.13-0.24) to see
#   which shared-state processes the data allow.
#
# Estimators (all should give 0):
#   1. RI-CLPM cross-lags x(t-1) -> y(t) and y(t-1) -> x(t), standardized (as in 05: E2/E3)
#   2. Hybrid (person-mean-centred) lagged regression of y(t) on x(t-1), with
#      and without prior y  -- Sorjonen et al.'s recommended check
#   3. ML dynamic panel model, standard and with partial predetermination
# Packages: lavaan, dpm, panelr
# Output  : output_main/S4_common_state_simulation.csv
#===============================================================================
suppressPackageStartupMessages({ library(lavaan); library(dpm); library(panelr) })
set.seed(20261003)
N <- 1300; A <- 1:8; REPS <- as.integer(Sys.getenv("REPS", "50")); S_SHARE <- 0.15; TARGET_STAB <- 0.20

simulate <- function(rho, share = S_SHARE, share_y = share) {
  phi_x <- max(0, (TARGET_STAB - share * rho) / (1 - share))
  phi_y <- max(0, (TARGET_STAB - share_y * rho) / (1 - share_y))
  ar <- function(p) { m <- matrix(0, N, 8); m[, 1] <- rnorm(N)
    for (t in 2:8) m[, t] <- p * m[, t - 1] + rnorm(N, sd = sqrt(1 - p^2)); m }
  C <- ar(rho); ux <- ar(phi_x); uy <- ar(phi_y)
  a <- rnorm(N); b <- 0.5 * a + rnorm(N, sd = 0.9)          # correlated stable levels
  list(x = a + sqrt(share) * C + sqrt(1 - share) * ux,
       y = b + sqrt(share_y) * C + sqrt(1 - share_y) * uy)
}

syn_ri <- local({
  s <- c(paste("RIx =~", paste0("1*x", A, collapse = " + ")), paste("RIy =~", paste0("1*y", A, collapse = " + ")))
  for (a in A) s <- c(s, sprintf("wx%d =~ 1*x%d", a, a), sprintf("wy%d =~ 1*y%d", a, a),
                      sprintf("x%d ~~ 0*x%d", a, a), sprintf("y%d ~~ 0*y%d", a, a))
  for (a in A[-1]) s <- c(s, sprintf("wx%d ~ wx%d + wy%d", a, a - 1, a - 1),
                          sprintf("wy%d ~ wy%d + wx%d", a, a - 1, a - 1))
  s <- c(s, sprintf("wx%d ~~ wy%d", A, A), sprintf("wx%d ~~ wx%d", A, A), sprintf("wy%d ~~ wy%d", A, A),
         "RIx ~~ RIx", "RIy ~~ RIy", "RIx ~~ RIy", "RIx ~~ 0*wx1", "RIx ~~ 0*wy1", "RIy ~~ 0*wx1", "RIy ~~ 0*wy1")
  paste(s, collapse = "\n")
})
riclpm <- function(s) {
  w <- data.frame(s$x / sd(s$x), s$y / sd(s$y)); names(w) <- c(paste0("x", A), paste0("y", A))
  f <- lavaan(syn_ri, data = w, meanstructure = TRUE, int.ov.free = TRUE)
  ss <- standardizedSolution(f)
  c(cross = mean(ss$est.std[ss$op == "~" & grepl("^wy", ss$lhs) & grepl("^wx", ss$rhs)]),
    reverse = mean(ss$est.std[ss$op == "~" & grepl("^wx", ss$lhs) & grepl("^wy", ss$rhs)]),
    same_age = mean(ss$est.std[ss$op == "~~" & grepl("^wx", ss$lhs) & grepl("^wy", ss$rhs)]),
    stab_x = mean(ss$est.std[ss$op == "~" & grepl("^wx", ss$lhs) & grepl("^wx", ss$rhs)]))
}
hybrid <- function(s, ldv) {   # linear, standardized units (SD of y per SD of x)
  x <- s$x / sd(s$x); y <- s$y / sd(s$y)
  d <- do.call(rbind, lapply(3:8, function(t) data.frame(id = 1:N, t = t, y = y[, t], y_l1 = y[, t - 1], x_l1 = x[, t - 1])))
  d$BW <- ave(d$x_l1, d$id); d$WI <- d$x_l1 - d$BW
  f <- if (ldv) y ~ WI + BW + y_l1 + factor(t) else y ~ WI + BW + factor(t)
  unname(coef(lm(f, data = d))["WI"])
}
# dpm syntax built once on a template panel, then fitted directly in lavaan (as in 05)
tmpl <- panel_data(data.frame(id = rep(1:3, each = 8), wave = rep(A, 3), y = rnorm(24), x = rnorm(24)), id = id, wave = wave)
syn_dpm <- lapply(c(std = FALSE, pp = TRUE), function(pp)
  paste(capture.output(dpm(y ~ pre(lag(x)), data = tmpl, partial.pre = pp, print.only = TRUE)), collapse = "\n"))
dpm_fit <- function(s, which) {
  w <- data.frame(s$y / sd(s$y), s$x / sd(s$x)); names(w) <- c(paste0("y_", A), paste0("x_", A))
  f <- tryCatch(sem(syn_dpm[[which]], data = w, fixed.x = FALSE), error = function(e) NULL)
  if (is.null(f)) return(NA)
  pe <- parameterEstimates(f); pe$est[pe$label == "en1"][1]
}

scen <- data.frame(rho = c(0, 0.5, 0.5, 0.8, 0.8, 0.95, 0.95, 0.5, 0.5),
                   share = c(0, 0.15, 0.30, 0.15, 0.30, 0.15, 0.30, 0.30, 0.10),
                   share_y = c(0, 0.15, 0.30, 0.15, 0.30, 0.15, 0.30, 0.10, 0.30))
scen$label <- ifelse(scen$share == 0, "No shared state",
                     sprintf("Shared state: persistence %.2f, share in X %.2f, share in offending %.2f",
                             scen$rho, scen$share, scen$share_y))
q <- function(v) sprintf("%.3f (%.3f to %.3f)", mean(v, na.rm = TRUE), quantile(v, .025, na.rm = TRUE), quantile(v, .975, na.rm = TRUE))
res <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) {
  r <- replicate(REPS, {
    s <- simulate(scen$rho[i], share = scen$share[i], share_y = scen$share_y[i])
    ri <- riclpm(s)
    c(ri, hyb_ldv = hybrid(s, TRUE), hyb_noldv = hybrid(s, FALSE),
      dpm_std = dpm_fit(s, "std"), dpm_pp = dpm_fit(s, "pp"))
  })
  data.frame(scenario = scen$label[i], persistence = scen$rho[i], share_x = scen$share[i], share_y = scen$share_y[i],
             sim_same_age_r = q(r["same_age", ]), sim_stability_x = q(r["stab_x", ]),
             riclpm_cross_lag = q(r["cross", ]), riclpm_reverse_lag = q(r["reverse", ]),
             hybrid_with_prior_y = q(r["hyb_ldv", ]), hybrid_no_prior_y = q(r["hyb_noldv", ]),
             dpm_standard = q(r["dpm_std", ]), dpm_partial_pre = q(r["dpm_pp", ]))
}))
cat("No causal effect in either direction in any scenario; every estimate should be 0.\n",
    REPS, "replications; mean (2.5% to 97.5%). Standardized units.\n\n")
print(res, row.names = FALSE)
dir.create("output_main", showWarnings = FALSE)
write.csv(res, file.path("output_main", "S4_common_state_simulation.csv"), row.names = FALSE)
