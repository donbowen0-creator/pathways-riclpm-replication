#===============================================================================
# 04_simulation_lag_artifacts.R
# Can a lagged within-person estimate be non-zero when the predictor has NO
# effect on later offending, only a same-period link? (2026-10-02)
#
# Simulated panel shaped like Pathways: 1,300 people, 8 yearly waves.
#   predictor  x = stable level + within-person swing u (persistence phi)
#   offending  y* = stable level + gamma * u (same period only) + e (own persistence 0.3)
#   NO lagged effect of x on y in any scenario.
# Each scenario is fitted three ways:
#   1. hybrid logit, previous-period x, prior-offending control (our main model)
#   2. hybrid logit, previous-period x, no prior-offending control
#   3. RI-CLPM cross-lag x(t-1) -> y(t) on the continuous outcome (standardized)
# A correct estimator would give OR = 1 and cross-lag = 0 throughout.
# Packages: lavaan
#===============================================================================
library(lavaan)
set.seed(20261002)
N <- 1300; A <- 1:8; REPS <- 50

simulate <- function(gamma, phi) {
  a <- rnorm(N); b <- 0.5 * a + rnorm(N, sd = 0.9)
  u <- matrix(0, N, 8); e <- matrix(0, N, 8); u[, 1] <- rnorm(N); e[, 1] <- rnorm(N)
  for (t in 2:8) {
    u[, t] <- phi * u[, t - 1] + rnorm(N, sd = sqrt(1 - phi^2))
    e[, t] <- 0.3 * e[, t - 1] + rnorm(N, sd = sqrt(1 - 0.09))
  }
  list(x = a + u, ystar = b + gamma * u + e)
}
hybrid <- function(s, ldv) {
  y <- (s$ystar > 0.8) * 1                                   # about 1 in 3 person-years offend
  d <- do.call(rbind, lapply(3:8, function(t) data.frame(id = 1:N, t = t, y = y[, t], y_l1 = y[, t - 1],
                                                         x_l1 = s$x[, t - 1])))
  d$BW <- ave(d$x_l1, d$id); d$WI <- d$x_l1 - d$BW
  f <- if (ldv) y ~ WI + BW + y_l1 + factor(t) else y ~ WI + BW + factor(t)
  unname(exp(coef(glm(f, data = d, family = binomial()))["WI"] * sd(s$x)))
}
syn <- local({
  s <- c(paste("RIx =~", paste0("1*x", A, collapse = " + ")), paste("RIy =~", paste0("1*y", A, collapse = " + ")))
  for (a in A) s <- c(s, sprintf("wx%d =~ 1*x%d", a, a), sprintf("wy%d =~ 1*y%d", a, a),
                      sprintf("x%d ~~ 0*x%d", a, a), sprintf("y%d ~~ 0*y%d", a, a))
  for (a in A[-1]) s <- c(s, sprintf("wx%d ~ wx%d + wy%d", a, a - 1, a - 1),
                          sprintf("wy%d ~ wy%d + cxy*wx%d", a, a - 1, a - 1))
  s <- c(s, sprintf("wx%d ~~ wy%d", A, A), sprintf("wx%d ~~ wx%d", A, A), sprintf("wy%d ~~ wy%d", A, A),
         "RIx ~~ RIx", "RIy ~~ RIy", "RIx ~~ RIy", "RIx ~~ 0*wx1", "RIx ~~ 0*wy1", "RIy ~~ 0*wx1", "RIy ~~ 0*wy1")
  paste(s, collapse = "\n")
})
riclpm <- function(s) {
  w <- data.frame(s$x / sd(s$x), s$ystar / sd(s$ystar)); names(w) <- c(paste0("x", A), paste0("y", A))
  ss <- standardizedSolution(lavaan(syn, data = w, meanstructure = TRUE, int.ov.free = TRUE))
  mean(ss$est.std[ss$op == "~" & grepl("^wy", ss$lhs) & grepl("^wx", ss$rhs)])   # standardized, as in 03 (D3)
}

scen <- data.frame(label = c("No same-period link", "Same-period link, swings don't persist",
                             "Same-period link, swings persist"),
                   gamma = c(0, 0.6, 0.6), phi = c(0.5, 0, 0.5))
res <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) {
  r <- replicate(REPS, { s <- simulate(scen$gamma[i], scen$phi[i])
    c(hyb_ldv = hybrid(s, TRUE), hyb_noldv = hybrid(s, FALSE), ri = riclpm(s)) })
  q <- function(v) sprintf("%.3f (%.3f to %.3f)", mean(v), quantile(v, .025), quantile(v, .975))
  data.frame(scenario = scen$label[i], phi = scen$phi[i],
             hybrid_OR_with_prior_offending = q(r["hyb_ldv", ]),
             hybrid_OR_no_prior_offending   = q(r["hyb_noldv", ]),
             riclpm_cross_lag               = q(r["ri", ]))
}))
cat("True lagged effect = 0 in every scenario (OR 1, cross-lag 0).", REPS, "replications; mean (2.5% to 97.5%).\n\n")
print(res, row.names = FALSE)
dir.create("output_checks", showWarnings = FALSE)
write.csv(res, file.path("output_checks", "S1_simulation_lag_artifacts.csv"), row.names = FALSE)

#===============================================================================
# Part 2. Calibrated to Pathways: how big a cross-lag can each predictor's
# same-age link produce on its own? Uses the RI-CLPM estimates from
# 03_robustness_checks.R (D3_riclpm_summary.csv): the same-age within-person
# correlation, the predictor's within-person stability and offending's stability.
# Again NO lagged effect is simulated. If the real cross-lag is well outside the
# simulated range, the same-age link alone can't explain it.
#===============================================================================
obs <- read.csv(file.path("output_checks", "D3_riclpm_summary.csv"))
stab_y <- mean(obs$mean_stability_offending)
cal <- do.call(rbind, lapply(seq_len(nrow(obs)), function(i) {
  r <- obs$mean_same_age_r[i]; g <- r / sqrt(1 - r^2)          # gives that same-age correlation
  phx <- max(obs$mean_stability_X[i], 0)
  sims <- replicate(30, {
    a <- rnorm(N); b <- 0.5 * a + rnorm(N, sd = 0.9)
    u <- matrix(0, N, 8); e <- matrix(0, N, 8); u[, 1] <- rnorm(N); e[, 1] <- rnorm(N)
    for (t in 2:8) {
      u[, t] <- phx * u[, t - 1] + rnorm(N, sd = sqrt(1 - phx^2))
      e[, t] <- stab_y * e[, t - 1] + rnorm(N, sd = sqrt(1 - stab_y^2))
    }
    riclpm(list(x = a + u, ystar = b + g * u + e))
  })
  data.frame(predictor = obs$predictor[i], same_age_r = r, stability_X = phx,
             observed_mean_cross_lag = obs$mean_X_to_off[i],
             null_mean = mean(sims), null_lo = quantile(sims, .025), null_hi = quantile(sims, .975),
             outside_null_range = obs$mean_X_to_off[i] < quantile(sims, .025) |
                                  obs$mean_X_to_off[i] > quantile(sims, .975), row.names = NULL)
}))
num <- sapply(cal, is.numeric); cal[num] <- lapply(cal[num], round, 3)
cat("\nCalibrated null (no lagged effect), offending stability =", round(stab_y, 3), "\n")
print(cal, row.names = FALSE)
write.csv(cal, file.path("output_checks", "S2_calibrated_null_cross_lags.csv"), row.names = FALSE)

#===============================================================================
# Part 3. Under Pathways-like dynamics, what does our main (hybrid) model show
# (a) when there is NO lagged effect, and (b) when the lagged effect equals the
# RI-CLPM estimate? Five key predictors, 30 replications each.
#===============================================================================
simulate_cal <- function(r, phx, beta) {
  g <- r / sqrt(1 - r^2)
  a <- rnorm(N); b <- 0.5 * a + rnorm(N, sd = 0.9)
  u <- matrix(0, N, 8); e <- matrix(0, N, 8); u[, 1] <- rnorm(N); e[, 1] <- rnorm(N)
  for (t in 2:8) {
    u[, t] <- phx * u[, t - 1] + rnorm(N, sd = sqrt(1 - phx^2))
    e[, t] <- stab_y * e[, t - 1] + rnorm(N, sd = sqrt(1 - stab_y^2))
  }
  ys <- b + g * u + e
  ys[, 2:8] <- ys[, 2:8] + beta * sqrt(g^2 + 1) * u[, 1:7]   # standardized within-person cross-lag ~ beta
  list(x = a + u, ystar = ys)
}
key <- c("Peer antisocial influence", "Moral disengagement", "Personal rewards of crime",
         "Impulse control", "Certainty of punishment")
p3 <- do.call(rbind, lapply(key, function(k) {
  o <- obs[obs$predictor == k, ]
  do.call(rbind, lapply(c(0, o$mean_X_to_off), function(beta) {
    r <- replicate(30, { s <- simulate_cal(o$mean_same_age_r, max(o$mean_stability_X, 0), beta)
      c(hyb = hybrid(s, TRUE), ri = riclpm(s)) })
    data.frame(predictor = k, true_lagged_effect = beta,
               hybrid_OR_mean = mean(r["hyb", ]), hybrid_lo = quantile(r["hyb", ], .025),
               hybrid_hi = quantile(r["hyb", ], .975),
               riclpm_mean = mean(r["ri", ]), riclpm_lo = quantile(r["ri", ], .025),
               riclpm_hi = quantile(r["ri", ], .975), row.names = NULL)
  }))
}))
num <- sapply(p3, is.numeric); p3[num] <- lapply(p3[num], round, 3)
cat("\nPart 3: hybrid model and RI-CLPM under Pathways-like dynamics\n")
print(p3, row.names = FALSE)
write.csv(p3, file.path("output_checks", "S3_hybrid_vs_riclpm_calibrated.csv"), row.names = FALSE)

#===============================================================================
# Part 4. Chart: what our main model would show with no lagged effect, what it
# would show if the RI-CLPM lagged effects are real, and what it actually showed
#===============================================================================
a1 <- read.csv(file.path("output_checks", "A1_positive_control.csv"))
a2 <- read.csv(file.path("output_checks", "A2_main_estimates_with_checks_added.csv"))
actual <- c(a2$OR_with_checks[match(key[-1], a2$predictor)],
            a1$OR_within[a1$predictor == key[1] & a1$model == "With all four domains" & a1$timing == "Previous year"])
names(actual) <- c(key[-1], key[1])
local({
  GREY <- "#8f8d88"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#0b0b0b"; MUTED <- "#52514e"
  k <- length(key); y <- rev(seq_len(k))
  nul <- p3$hybrid_OR_mean[p3$true_lagged_effect == 0]
  alt <- p3$hybrid_OR_mean[p3$true_lagged_effect != 0]
  act <- actual[key]
  png(file.path("output_checks", "fig_main_model_masks_lagged_effects.png"), width = 1250, height = 560,
      res = 120, bg = "#fcfcfb")
  par(mar = c(4.2, 14, 6.4, 1.5), fg = "#cfcdc8", col.axis = MUTED)
  plot(NA, xlim = c(0.72, 1.2), ylim = c(0.5, k + 0.5), log = "x", yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "l")
  xt <- c(0.75, 0.8, 0.9, 1, 1.1, 1.2); axis(1, at = xt, labels = xt, cex.axis = 0.85)
  abline(v = xt, col = "#ebe9e4", lwd = 0.8); abline(v = 1, col = MUTED, lty = 2)
  axis(2, at = y, labels = key, las = 1, cex.axis = 0.85, tick = FALSE, col.axis = INK)
  segments(pmin(nul, alt, act), y, pmax(nul, alt, act), y, col = "#dcdad5", lwd = 2)
  points(nul, y, pch = 19, col = GREY, cex = 1.4)
  points(alt, y, pch = 19, col = BLUE, cex = 1.4)
  points(act, y, pch = 18, col = ORANGE, cex = 2)
  mtext("Previous-year within-person odds ratio from our main (hybrid) model", side = 1, line = 2.6, cex = 0.8, col = MUTED)
  mtext("Our main model hides small lagged effects", side = 3, line = 4.6, adj = 0, at = 10^par("usr")[1],
        cex = 1.05, font = 2, col = INK)
  legend(x = 10^par("usr")[1], y = k + 0.5 + k * 0.33, xpd = NA, bty = "n", cex = 0.82, ncol = 3,
         legend = c("Simulated: no lagged effect", "Simulated: RI-CLPM lagged effect is real", "Actual Pathways result"),
         col = c(GREY, BLUE, ORANGE), pch = c(19, 19, 18), pt.cex = c(1.2, 1.2, 1.7), text.col = INK)
  dev.off()
})
