rm(list = ls())

## =========================================================
## LRR analysis: model + EMMs + a/b letters + publication plot
## =========================================================

## --- Libraries (keep only what is actually used)
library(dplyr)      # data wrangling
library(ggplot2)    # plotting
library(glmmTMB)    # mixed models
library(DHARMa)     # diagnostics
library(emmeans)    # marginal means + contrasts
library(multcomp)   # cld() for compact letter display
library(multcompView)# cld() for compact letter display

## =========================================================
## Load and prepare data
## =========================================================
data <- read.csv2("data/fireflies2024.csv", header = TRUE, sep = ";", dec = ".")

# Remove biodynamic treatment (not part of this comparison)
data <- subset(data, type != "biodyn")

# Quick check: total males (optional; keep if you like reporting it)
sum(data$splend)

# Convert to factors used in model / plotting
data <- data %>%
  mutate(
    region     = factor(region),
    locality   = factor(locality),
    side       = factor(side),
    lux_factor = factor(lux),
    type       = factor(type)
  )

## =========================================================
## Compute LRR (log-response ratio) relative to baseline (0.1 lux)
## Baseline is computed per locality.
## eps avoids division by zero and log(0) issues.
## =========================================================
eps <- 0.5

data <- data %>%
  group_by(locality) %>%
  mutate(
    baseline = mean(splend[lux == 0.1], na.rm = TRUE),
    rel_abund = (splend + eps) / (baseline + eps),
    LRR = log(rel_abund)
  ) %>%
  ungroup()

# Optional sanity check: baseline values per locality
# data %>% group_by(locality) %>% summarise(baseline = mean(splend[lux == 0.1], na.rm = TRUE))

## =========================================================
## Model: LRR ~ lux * type with locality random intercept
## =========================================================
mod1 <- glmmTMB(LRR ~ lux_factor * type + (1 | locality), data = data)

## --- Diagnostics (keep during development; comment out for final run if desired)
sim_res <- simulateResiduals(mod1)
plot(sim_res)
plotResiduals(sim_res, data$lux_factor)

## --- Model output (useful for Results)
summary(mod1)

# Pairwise comparisons among all lux_factor × type combinations (Tukey-adjusted)
post_type <- emmeans(mod1, pairwise ~ lux_factor * type)
post_type

## =========================================================
## EMMs + compact letters (a/b) within each lux level
## =========================================================

# EMMs for: type differences within each lux_factor
emm <- emmeans(mod1, ~ type | lux_factor)

# Table for plotting: EMM + 95% CI
emm_df <- as.data.frame(emm) %>%
  transmute(
    lux_factor,
    type,
    emmean,
    lower.CL = asymp.LCL,
    upper.CL = asymp.UCL
  )

# Compact letter display (Tukey-adjusted) within each lux_factor
cld_df <- as.data.frame(multcomp::cld(emm, adjust = "tukey", Letters = "abcd")) %>%
  transmute(
    lux_factor,
    type,
    group = gsub("\\s+", "", .group)   # remove spaces in letter groupings
  )

# Merge letters with EMMs and set y position for labels per lux level
emm_lab <- left_join(emm_df, cld_df, by = c("lux_factor", "type")) %>%
  group_by(lux_factor) %>%
  mutate(label_y = max(upper.CL, na.rm = TRUE) + 0.25) %>%  # headroom for a/b letters
  ungroup()

## =========================================================
## Plot: raw data (subtle) + EMM ± CI (dominant) + a/b letters
## =========================================================
pd <- position_dodge(width = 0.5)

lamp_cols <- c(
  LED    = "#0070c0",  # blue
  sodium = "#ffc000"   # orange
)

p_final_ab <- ggplot() +
  # Raw data as background texture (same shape, coloured by type)
  geom_point(
    data = data,
    aes(x = lux_factor, y = LRR, color = type),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.5),
    alpha = 0.25,
    size = 1.8,
    shape = 16
  ) +
  # Model-based means
  geom_point(
    data = emm_df,
    aes(x = lux_factor, y = emmean, color = type),
    position = pd,
    size = 3,
    shape = 16
  ) +
  # Model-based 95% confidence intervals (no caps for consistent style)
  geom_linerange(
    data = emm_df,
    aes(x = lux_factor, ymin = lower.CL, ymax = upper.CL, color = type),
    position = pd,
    linewidth = 0.9
  ) +
  # a/b letters (Tukey-adjusted), one per type within each lux level
  geom_text(
    data = emm_lab,
    aes(x = lux_factor, y = label_y, label = group, color = type),
    position = pd,
    show.legend = FALSE,
    fontface = "bold",
    size = 5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = lamp_cols, name = "Type") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(
    x = "Light intensity (lux)",
    y = "Log-response ratio (LRR)"
  ) +
  theme_classic(base_size = 12)

p_final_ab

## =========================================================
## Export figure (journal-ready TIFF)
## =========================================================
ggsave("figures/Fig4.tiff", plot = p_final_ab,
       width = 7, height = 5, units = "in",
       dpi = 600, compression = "lzw")
