rm(list = ls())
## =========================================================
##  Effect of the light introduction (Biodynamic) at pristine 
##  nature (2025)
## =========================================================

## =========================================================
## Load required libraries
## =========================================================
library(dplyr)       # data manipulation
library(ggplot2)     # plotting
library(gridExtra)   # arranging plots (not used here, but kept)
library(lmerTest)    # mixed-model utilities
library(sjPlot)      # mixed-model plotting helpers
library(glmmTMB)     # GLMM fitting
library(emmeans)     # estimated marginal means & contrasts
library(DHARMa)      # model diagnostics
library(ggbreak)     # add scale break

## =========================================================
## Load data
## =========================================================
data <- read.csv2(
  "data/fireflies2025.csv",
  header = TRUE, sep = ";", dec = "."
)

## =========================================================
## Factor preparation
## =========================================================
data$line     <- as.factor(data$line)
data$trap     <- as.factor(data$trap)
data$locality <- as.factor(data$locality)

# IMPORTANT: set timing order BEFORE model fitting
# This defines reference levels and ensures consistent EMM ordering
data$timing <- factor(
  data$timing,
  levels = c("control1", "light", "control2")
)

## =========================================================
## Colour grouping for plotting
## (collapse both controls into one colour)
## =========================================================
data <- data %>%
  mutate(
    treat_col = ifelse(grepl("^control", timing), "control", "light"),
    treat_col = factor(treat_col, levels = c("control", "light"))
  )

# Publication-safe colour palette
treat_cols <- c(
  control = "#0070C0",  # blue (control before + after)
  light   = "#FFC000"   # orange (light treatment)
)
data$luxF <- as.factor(data$lux)
summary(data)

data$abund[data$lux >=0.1 & data$timing== "light"]


sum(data$abund, na.rm = TRUE)
sum(data$abund[data$timing== "light"], na.rm = TRUE)
sum(data$abund[data$timing== "control1"], na.rm = TRUE)
sum(data$abund[data$timing== "control2"], na.rm = TRUE)
## =========================================================
## Adding a distance as a factor
## =========================================================
data$distance_f <- as.factor(data$distance)

## =========================================================
## Model: male abundance ~ timing
## =========================================================
mod1 <- glmmTMB(
  abund ~ timing + locality+ (1 | line) +(1 | trap),
  ziformula = ~1,          # constant zero-inflation
  family = poisson,
  data = data,
  
)


## =========================================================
## Model diagnostics
## =========================================================
sim_res_f <- simulateResiduals(mod1)
plot(sim_res_f)             # residual diagnostics
testZeroInflation(sim_res_f)


## =========================================================
## Model summary and post hoc contrasts
## =========================================================
summary(mod1)

performance::icc(mod1)

# Pairwise comparisons among timing levels (Tukey-adjusted)
post_type <- emmeans(mod1, pairwise ~ timing)
post_type

## =========================================================
## Figure 5: Male abundance across temporal treatments
## =========================================================

## --- Estimated marginal means on response scale (counts)
emm_resp <- emmeans(mod1, ~ timing, type = "response")
emm_resp

emm_df <- as.data.frame(emm_resp) %>%
  transmute(
    timing,
    mean  = rate,        # back-transformed mean count
    lower = asymp.LCL,   # lower 95% CI
    upper = asymp.UCL    # upper 95% CI
  )


##Prepare a tabel 3
coef_tab <- summary(mod1)$coefficients$cond

ci_tab <- data.frame(
  parameter = rownames(coef_tab),
  estimate = coef_tab[, "Estimate"],
  se = coef_tab[, "Std. Error"],
  lower_95 = coef_tab[, "Estimate"] - 1.96 * coef_tab[, "Std. Error"],
  upper_95 = coef_tab[, "Estimate"] + 1.96 * coef_tab[, "Std. Error"],
  z = coef_tab[, "z value"],
  p = coef_tab[, "Pr(>|z|)"]
) %>%
  mutate(
    across(c(estimate, se, lower_95, upper_95, z), ~ round(.x, 3)),
    p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

ci_tab
## --- Compact letter display (Tukey-adjusted)
## Letters are computed on the link (log) scale,
## but can be safely displayed on the response scale
emm_link <- emmeans(mod1, ~ timing)
cld_df <- as.data.frame(
  multcomp::cld(emm_link, adjust = "tukey", Letters = "abcd")
) %>%
  transmute(
    timing,
    group = gsub("\\s+", "", .group)   # clean letter formatting
  )

## --- Combine EMMs, letters, and colour grouping
plot_df <- emm_df %>%
  left_join(cld_df, by = "timing") %>%
  mutate(
    treat_col = ifelse(grepl("^control", timing), "control", "light"),
    treat_col = factor(treat_col, levels = c("control", "light")),
    label_y   = upper + 0.35            # vertical offset for letters
  )

## --- Determine y-axis headroom dynamically
y_top <- max(c(data$abund, plot_df$upper), na.rm = TRUE)

## =========================================================
## Prepare shape grouping
## =========================================================
## Distinguish only the illuminated treatment from both controls
data$shape_grp <- ifelse(data$timing == "light", "light", "control")
plot_df$shape_grp <- ifelse(plot_df$timing == "light", "light", "control")


## =========================================================
## Base plot
## =========================================================
p <- ggplot() +
  
  # Raw observations (background texture only)
  geom_jitter(
    data = data,
    aes(
      timing,
      abund,
      color = treat_col,
      shape = shape_grp
    ),
    width = 0.08, height = 0,
    alpha = 0.12,
    size = 1.6
  ) +
  
  # Model-based 95% confidence intervals
  geom_linerange(
    data = plot_df,
    aes(
      timing,
      ymin = lower,
      ymax = upper,
      color = treat_col
    ),
    linewidth = 1.2
  ) +
  
  # Model-based means
  geom_point(
    data = plot_df,
    aes(
      timing,
      mean,
      color = treat_col,
      shape = shape_grp
    ),
    size = 4.2
  ) +
  
  # Significance letters (a/b)
  geom_text(
    data = plot_df,
    aes(timing, label_y, label = group),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  
  # Colour scale
  scale_color_manual(
    values = treat_cols,
    name = NULL
  ) +
  
  # Shape scale: light differs from both controls
  scale_shape_manual(
    values = c(
      control = 16,  # filled circle
      light   = 17   # filled triangle
    ),
    name = NULL
  ) +
  
  # Initial y-limits (expanded later by axis break)
  coord_cartesian(
    ylim = c(0, y_top + 2.0),
    clip = "off"
  ) +
  
  # Axis labels
  labs(
    x = "Timing of the experiment",
    y = "Abundance (count/trap night)"
  ) +
  
  # Theme
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

p


## =========================================================
## Broken y-axis to handle single high outlier
## =========================================================
p <- p +
  scale_y_break(
    c(6, 12),
    ticklabels = c(12, 13),
    space = 0.1
  )

p


## =========================================================
## Export figure (journal-ready TIFF)
## =========================================================
ggsave(
  "figures/Fig9.tiff",
  plot = p,
  width = 5,
  height = 3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)