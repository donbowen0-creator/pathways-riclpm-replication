#===============================================================================
# 07_figures.R  --  Pathways to Desistance (ICPSR 29961)
# Figures 1 and 2 of the manuscript.
#
# Fig. 1  Same-age associations and next-year effects, 14 scales
#         (05_riclpm_main.R -> output_main/E3_summary_by_band.csv, band "all")
# Fig. 2  Next-year effect of each scale on offending from three estimators:
#         lag-only dynamic panel, dynamic panel with the same-year link allowed
#         (partial predetermination), and the age-invariant RI-CLPM
#         (output_main/I1_dynamic_panel_ml.csv and J1_age_bands_vs_invariant.csv)
#
# Requires: the output_main/ folder from 05_riclpm_main.R
# Packages: install.packages("ggplot2")
# Output  : figures/Fig1.tif, figures/Fig1.eps, figures/Fig2.tif, figures/Fig2.eps
#           (174 mm wide; TIFF at 600 dpi; grayscale; Arial/Helvetica lettering)
# Run     : from the ICPSR_29961 folder, source("07_figures.R")
#===============================================================================

suppressPackageStartupMessages(library(ggplot2))

in_dir  <- "output_main"
out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)

short <- c(EXPECT = "Expectations for the future", MORDIS = "Moral disengagement",
           SOCCST = "Social costs of punishment", PUNYOU = "Certainty of punishment",
           PERRWD = "Personal rewards of crime", SOCRWD = "Social rewards of crime",
           PSMIID = "Identity (PSMI)", IMPULS = "Impulse control",
           SUPAGG = "Suppression of aggression", CONSID = "Consideration of others",
           FUTURE = "Future orientation", PEERPR = "Resistance to peer influence",
           PSMI = "Psychosocial maturity (PSMI)", PRINFL = "Peer antisocial influence (control)")

theme_pub <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
        axis.text.y = element_text(colour = "black", size = 8),
        axis.text.x = element_text(colour = "grey25", size = 7.5),
        axis.title.x = element_text(size = 8, margin = margin(t = 6)),
        legend.position = "top", legend.justification = "left",
        legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(4, 12, 4, 4))

save_fig <- function(p, name, height_mm) {
  ggsave(file.path(out_dir, paste0(name, ".tif")), p, width = 174, height = height_mm,
         units = "mm", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  eps_ok <- tryCatch({
    ggsave(file.path(out_dir, paste0(name, ".eps")), p, width = 174, height = height_mm,
           units = "mm", device = cairo_ps, bg = "white"); TRUE
  }, error = function(e) FALSE)
  if (!eps_ok) ggsave(file.path(out_dir, paste0(name, ".eps")), p, width = 174, height = height_mm,
                      units = "mm", device = "eps", bg = "white")
  cat("Saved", file.path(out_dir, paste0(name, ".tif")), "and .eps\n")
}

#---- Fig. 1 ---------------------------------------------------------------------
e3 <- read.csv(file.path(in_dir, "E3_summary_by_band.csv"))
e3 <- e3[e3$band == "all" & e3$variable %in% names(short), ]
e3$label <- short[e3$variable]
e3 <- e3[order(abs(e3$same_age_r)), ]                     # largest at the top
e3$label <- factor(e3$label, levels = e3$label)

long1 <- rbind(
  data.frame(label = e3$label, est = e3$same_age_r, what = "Same age"),
  data.frame(label = e3$label, est = e3$X_to_off,   what = "Measure to next year's offending"),
  data.frame(label = e3$label, est = e3$off_to_X,   what = "Offending to next year's measure"))
long1$what <- factor(long1$what, levels = c("Same age", "Measure to next year's offending",
                                            "Offending to next year's measure"))
span <- data.frame(label = e3$label,
                   lo = pmin(e3$same_age_r, e3$X_to_off, e3$off_to_X),
                   hi = pmax(e3$same_age_r, e3$X_to_off, e3$off_to_X))
ratio <- data.frame(label = e3$label,
                    txt = ifelse(abs(e3$same_age_r) >= 0.12,
                                 sprintf("%.1f x", abs(e3$same_age_r) / pmax(abs(e3$X_to_off), 0.001)), ""))

p1 <- ggplot(long1, aes(x = est, y = label)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_segment(data = span, aes(x = lo, xend = hi, y = label, yend = label),
               inherit.aes = FALSE, colour = "grey80", linewidth = 1.6, lineend = "round") +
  # the two next-year estimates are nudged up and down slightly so that equal values stay visible
  geom_point(aes(shape = what, fill = what, size = what,
                 y = as.numeric(label) + c(0, 0.17, -0.17)[as.numeric(what)]),
             colour = "black", stroke = 0.5) +
  geom_text(data = ratio, aes(x = 0.46, y = label, label = txt), inherit.aes = FALSE,
            hjust = 1, size = 2.7) +
  annotate("text", x = 0.46, y = nrow(e3) + 0.9, label = "Same age / next year",
           hjust = 1, size = 2.6, colour = "grey25") +
  scale_shape_manual(values = c(21, 21, 21)) +
  scale_fill_manual(values = c("black", "grey60", "white")) +
  scale_size_manual(values = c(2.4, 2.0, 2.0)) +
  scale_x_continuous(limits = c(-0.22, 0.47), breaks = seq(-0.2, 0.4, 0.1),
                     labels = function(b) sub("^(-?)0\\.", "\\1.", sprintf("%.1f", b))) +
  coord_cartesian(clip = "off") +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE),
         size = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Standardized within-person estimate, ages 16 to 23\n(left of zero = lower offending)",
       y = NULL) +
  theme_pub
save_fig(p1, "Fig1", height_mm = 112)

#---- Fig. 2 ---------------------------------------------------------------------
i1 <- read.csv(file.path(in_dir, "I1_dynamic_panel_ml.csv"))
j1 <- read.csv(file.path(in_dir, "J1_age_bands_vs_invariant.csv"))
risk <- c("MORDIS", "PERRWD", "SOCRWD", "PRINFL")             # theory predicts more offending
rows <- list()
for (v in names(short)) {
  s <- if (v %in% risk) 1 else -1                              # flip protective measures
  d <- i1[i1$variable == v & i1$direction == "X_to_off", ]
  lag <- d[!as.logical(d$partial_pre), ]; pp <- d[as.logical(d$partial_pre), ]
  j <- j1[j1$variable == v, ]
  fl <- function(e, lo, hi) c(e * s, min(lo * s, hi * s), max(lo * s, hi * s))
  rows[[v]] <- rbind(
    data.frame(variable = v, est = "Lag-only dynamic panel", t(fl(lag$est, lag$lo, lag$hi))),
    data.frame(variable = v, est = "Dynamic panel, same-year link allowed", t(fl(pp$est, pp$lo, pp$hi))),
    data.frame(variable = v, est = "RI-CLPM, one effect for all ages", t(fl(j$inv_X_to_off, j$inv_lo, j$inv_hi))))
}
f2 <- do.call(rbind, rows); names(f2)[3:5] <- c("b", "lo", "hi")
ord <- f2[f2$est == "Dynamic panel, same-year link allowed", ]
f2$label <- factor(short[f2$variable], levels = short[ord$variable[order(ord$b)]])
f2$est <- factor(f2$est, levels = c("Lag-only dynamic panel", "Dynamic panel, same-year link allowed",
                                    "RI-CLPM, one effect for all ages"))
dodge <- position_dodge(width = 0.7)
p2 <- ggplot(f2, aes(x = b, y = label, group = est)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = est), position = dodge, linewidth = 0.45) +
  geom_point(aes(shape = est, fill = est), position = dodge, size = 1.9, colour = "black", stroke = 0.45) +
  annotate("text", x = -0.005, y = length(levels(f2$label)) + 0.75, label = "Opposite to theory",
           hjust = 1, size = 2.6, colour = "grey25") +
  annotate("text", x = 0.005, y = length(levels(f2$label)) + 0.75, label = "Direction theory predicts",
           hjust = 0, size = 2.6, colour = "grey25") +
  scale_shape_manual(values = c(21, 21, 21)) +
  scale_fill_manual(values = c("grey60", "black", "white")) +
  scale_colour_manual(values = c("grey55", "black", "grey55")) +
  scale_x_continuous(breaks = seq(-0.15, 0.15, 0.05),
                     labels = function(b) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", b))) +
  coord_cartesian(clip = "off") +
  guides(colour = "none", shape = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Effect of last year's measure on this year's offending (SD units, 95% CI)\nProtective measures reversed, so right of zero is the direction theory predicts",
       y = NULL) +
  theme_pub
save_fig(p2, "Fig2", height_mm = 150)

cat("\nFigures written to", normalizePath(out_dir), "\n")
