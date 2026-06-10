rm(list = ls())

## =========================================================
## Fireflies 2024
## Alternative exposure metric: photon flux
##
## Main response:
## LRR = log-response ratio of male abundance relative to
## the locality-specific 0.1 lux baseline.
##
## Random effects:
## locality + locality-specific trap line
## =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)

## =========================================================
## User settings
## =========================================================

data_file <- "data/fireflies2024.csv"
out_dir <- "figures"

reference_lux <- 0.1
eps <- 0.5

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

lamp_cols <- c(
  HPS = "#ffc000",
  LED = "#0070c0"
)

## =========================================================
## Load and clean data
## =========================================================

fireflies_raw <- read.csv2(
  data_file,
  header = TRUE,
  sep = ";",
  dec = "."
)

fireflies <- fireflies_raw %>%
  filter(type != "biodyn") %>%
  mutate(
    lux = as.numeric(lux),
    splend = as.numeric(splend),
    
    ## Recode consistently.
    type = case_when(
      type %in% c("sodium", "hps", "HPS") ~ "HPS",
      type %in% c("LED", "led") ~ "LED",
      TRUE ~ as.character(type)
    ),
    type = factor(type, levels = c("HPS", "LED")),
    
    region = factor(region),
    locality = factor(locality),
    side = factor(side)
  ) %>%
  filter(
    !is.na(lux),
    !is.na(splend),
    !is.na(locality),
    !is.na(side),
    !is.na(type)
  ) %>%
  droplevels() %>%
  mutate(
    ## Locality-specific trap line.
    ## Equivalent to locality:side.
    line_id = interaction(locality, side, drop = TRUE)
  )

if (any(fireflies$lux <= 0, na.rm = TRUE)) {
  stop("lux contains zero or negative values.")
}

cat("Total males:", sum(fireflies$splend, na.rm = TRUE), "\n")
cat("Number of observations:", nrow(fireflies), "\n\n")

cat("Observations by lamp type and lux:\n")
print(xtabs(~ type + lux, data = fireflies))
cat("\n")

cat("Observations by locality-side line and lux:\n")
print(
  fireflies %>%
    group_by(locality, side, line_id, lux) %>%
    summarise(
      n = n(),
      males = sum(splend, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(locality, side, lux)
)
cat("\n")

## =========================================================
## Photon flux lookup table
## =========================================================
## These values must match the units used in your manuscript.

photon_table <- tibble(
  lux = c(0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8),
  HPS = c(
    0.00065936,
    0.001318719,
    0.002637439,
    0.006593597,
    0.013187194,
    0.026374387,
    0.052748775,
    0.105497549
  ),
  LED = c(
    0.000842739,
    0.001685479,
    0.003370957,
    0.008427394,
    0.016854787,
    0.033709574,
    0.067419149,
    0.134838297
  )
)

photon_long <- photon_table %>%
  pivot_longer(
    cols = c("HPS", "LED"),
    names_to = "type",
    values_to = "photon_flux"
  ) %>%
  mutate(type = factor(type, levels = levels(fireflies$type)))

fireflies <- fireflies %>%
  left_join(photon_long, by = c("lux", "type"))

if (any(is.na(fireflies$photon_flux))) {
  print(
    fireflies %>%
      filter(is.na(photon_flux)) %>%
      distinct(lux, type)
  )
  stop("Some photon flux values were not matched. Check lux and type coding.")
}

if (any(fireflies$photon_flux <= 0, na.rm = TRUE)) {
  stop("photon_flux contains zero or negative values.")
}

## =========================================================
## Derived variables
## =========================================================

fireflies <- fireflies %>%
  group_by(locality) %>%
  mutate(
    baseline = mean(splend[lux == reference_lux], na.rm = TRUE),
    rel_abund = (splend + eps) / (baseline + eps),
    LRR = log(rel_abund),
    log2_lux = log2(lux),
    log2_photon_flux = log2(photon_flux),
    found = as.integer(splend > 0),
    lux_factor = factor(
      lux,
      levels = c(0.1, 0.5, 1),
      labels = c("0.1", "0.5", "1")
    )
  ) %>%
  ungroup()

if (any(is.na(fireflies$baseline))) {
  print(
    fireflies %>%
      filter(is.na(baseline)) %>%
      distinct(locality)
  )
  stop("Some localities do not have a valid 0.1 lux baseline.")
}

cat("Baseline per locality:\n")
print(
  fireflies %>%
    distinct(locality, baseline) %>%
    arrange(locality)
)
cat("\n")

cat("Photon flux values used in the analysis:\n")
print(
  fireflies %>%
    distinct(type, lux, photon_flux, log2_photon_flux) %>%
    arrange(type, lux)
)
cat("\n")

## =========================================================
## Main inferential dataset
## =========================================================
## 0.1 lux is excluded because it defines the LRR reference.

fireflies_model <- fireflies %>%
  filter(lux != reference_lux) %>%
  droplevels()

cat("Number of observations used in the LRR model:", nrow(fireflies_model), "\n\n")

cat("Model observations by lamp type and lux:\n")
print(xtabs(~ type + lux, data = fireflies_model))
cat("\n")

## =========================================================
## Model selection: LRR ~ photon flux
## =========================================================
## Same random-effect structure in all models:
## (1 | locality) + (1 | line_id)

mod_full <- glmmTMB(
  LRR ~ log2_photon_flux * type + (1 | locality) + (1 | line_id),
  data = fireflies_model
)

mod_add <- glmmTMB(
  LRR ~ log2_photon_flux + type + (1 | locality) + (1 | line_id),
  data = fireflies_model
)

mod_flux_only <- glmmTMB(
  LRR ~ log2_photon_flux + (1 | locality) + (1 | line_id),
  data = fireflies_model
)

## Final model.
## Keep the additive model if the interaction is not supported.
mod_final <- mod_add

cat("=====================================================\n")
cat("MODEL COMPARISON: LRR PHOTON-FLUX MODELS\n")
cat("=====================================================\n")

cat("AIC comparison:\n")
print(AIC(mod_flux_only, mod_add, mod_full))
cat("\n")

cat("Likelihood-ratio test: additive vs interaction model\n")
print(anova(mod_add, mod_full))
cat("\n")

cat("Likelihood-ratio test: flux-only vs additive model\n")
print(anova(mod_flux_only, mod_add))
cat("\n")

## =========================================================
## Final model summary
## =========================================================

cat("=====================================================\n")
cat("FINAL PHOTON-FLUX MODEL SUMMARY\n")
cat("=====================================================\n")
print(summary(mod_final))
cat("\n")

## =========================================================
## Fixed-effect confidence intervals and effect sizes
## =========================================================

coef_tab <- summary(mod_final)$coefficients$cond

ci_tab <- data.frame(
  term = rownames(coef_tab),
  estimate = coef_tab[, "Estimate"],
  se = coef_tab[, "Std. Error"],
  lower_95 = coef_tab[, "Estimate"] - 1.96 * coef_tab[, "Std. Error"],
  upper_95 = coef_tab[, "Estimate"] + 1.96 * coef_tab[, "Std. Error"],
  row.names = NULL
) %>%
  mutate(
    ratio = exp(estimate),
    ratio_low = exp(lower_95),
    ratio_high = exp(upper_95),
    percent_change = (ratio - 1) * 100
  )

cat("95% Wald CI for fixed effects:\n")
print(ci_tab)
cat("\n")

b <- fixef(mod_final)$cond

beta_flux <- unname(b["log2_photon_flux"])
flux_ratio <- exp(beta_flux)

cat("Effect of photon flux:\n")
cat("- beta for log2(photon flux):", beta_flux, "\n")
cat("- multiplicative change in relative abundance for each doubling of photon flux:",
    flux_ratio, "\n")
cat("- percent change for each doubling of photon flux:",
    round((flux_ratio - 1) * 100, 1), "%\n\n")

if ("typeLED" %in% names(b)) {
  beta_LED <- unname(b["typeLED"])
  LED_ratio <- exp(beta_LED)
  
  cat("Effect of lamp type, LED vs HPS at the same photon flux:\n")
  cat("- beta for typeLED:", beta_LED, "\n")
  cat("- multiplicative ratio in relative abundance:", LED_ratio, "\n")
  cat("- percent difference of LED relative to HPS:",
      round((LED_ratio - 1) * 100, 1), "%\n\n")
}

## =========================================================
## Diagnostics
## =========================================================

sim_final <- simulateResiduals(mod_final, n = 1000)

plot(sim_final)

cat("DHARMa diagnostics for final photon-flux model:\n")
print(testUniformity(sim_final))
print(testDispersion(sim_final))
print(testOutliers(sim_final))
print(testQuantiles(sim_final))
cat("\n")

plotResiduals(sim_final, form = fireflies_model$log2_photon_flux)
plotResiduals(sim_final, form = fireflies_model$type)
plotResiduals(sim_final, form = fireflies_model$locality)
plotResiduals(sim_final, form = fireflies_model$line_id)

## =========================================================
## Sensitivity to pseudo-count eps
## =========================================================

run_eps_sensitivity <- function(eps_value, data_input) {
  
  tmp <- data_input %>%
    group_by(locality) %>%
    mutate(
      baseline_tmp = mean(splend[lux == reference_lux], na.rm = TRUE),
      LRR_tmp = log((splend + eps_value) / (baseline_tmp + eps_value)),
      log2_photon_flux_tmp = log2(photon_flux)
    ) %>%
    ungroup() %>%
    filter(lux != reference_lux) %>%
    droplevels()
  
  m <- glmmTMB(
    LRR_tmp ~ log2_photon_flux_tmp + type + (1 | locality) + (1 | line_id),
    data = tmp
  )
  
  ct <- summary(m)$coefficients$cond
  
  data.frame(
    eps = eps_value,
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"],
    row.names = NULL
  ) %>%
    mutate(
      ratio = exp(estimate),
      percent_change = (ratio - 1) * 100
    )
}

eps_sens <- bind_rows(
  run_eps_sensitivity(0.1, fireflies),
  run_eps_sensitivity(0.5, fireflies),
  run_eps_sensitivity(1.0, fireflies)
)

cat("Sensitivity of fixed effects to pseudo-count eps:\n")
print(eps_sens)
cat("\n")

## =========================================================
## Optional analysis A: binary capture/success model
## =========================================================
## This model asks whether at least one male was captured in a trap.
## I would not include lamp type here as a main model, because type is
## confounded with locality and the response is not baseline-standardized.

mod_prob_flux <- glmmTMB(
  found ~ log2_photon_flux + (1 | locality) + (1 | line_id),
  family = binomial(link = "logit"),
  data = fireflies
)

cat("=====================================================\n")
cat("OPTIONAL MODEL: BINARY CAPTURE PROBABILITY\n")
cat("=====================================================\n")
print(summary(mod_prob_flux))
cat("\n")

sim_prob_flux <- simulateResiduals(mod_prob_flux, n = 1000)
plot(sim_prob_flux)

cat("DHARMa diagnostics for binary photon-flux model:\n")
print(testUniformity(sim_prob_flux))
print(testDispersion(sim_prob_flux))
print(testOutliers(sim_prob_flux))
cat("\n")

## Optional exploratory version with lamp type.
## Do not use as main inference unless you explicitly frame it as exploratory.
mod_prob_flux_type <- glmmTMB(
  found ~ log2_photon_flux + type + (1 | locality) + (1 | line_id),
  family = binomial(link = "logit"),
  data = fireflies
)

cat("Exploratory binary model with lamp type:\n")
print(summary(mod_prob_flux_type))
cat("\n")

## =========================================================
## Publication-style plot
## =========================================================
## Because the model uses only 0.5 and 1 lux, the plot is shown
## at the observed flux levels rather than as a long extrapolated curve.

pred_grid <- expand.grid(
  lux = sort(unique(fireflies_model$lux)),
  type = levels(fireflies_model$type)
) %>%
  left_join(photon_long, by = c("lux", "type")) %>%
  mutate(
    log2_photon_flux = log2(photon_flux),
    type = factor(type, levels = levels(fireflies_model$type)),
    
    ## Required by the model formula, ignored because re.form = NA.
    locality = fireflies_model$locality[1],
    line_id = fireflies_model$line_id[1],
    
    lux_factor = factor(
      lux,
      levels = c(0.1, 0.5, 1),
      labels = c("0.1", "0.5", "1")
    )
  )

pred <- predict(
  mod_final,
  newdata = pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA
)

pred_grid <- pred_grid %>%
  mutate(
    fit_LRR = pred$fit,
    ci_low_LRR = pred$fit - 1.96 * pred$se.fit,
    ci_high_LRR = pred$fit + 1.96 * pred$se.fit
  )

pd <- position_dodge(width = 0.5)

p_flux <- ggplot() +
  
  ## Reference 0.1 lux observations, shown only as visual reference.
  geom_point(
    data = fireflies %>% filter(lux == reference_lux),
    aes(x = lux_factor, y = LRR),
    position = position_jitter(
      width = 0.10,
      height = 0,
      seed = 1
    ),
    shape = 21,
    fill = "white",
    colour = "grey45",
    alpha = 0.65,
    size = 2.0,
    stroke = 0.6
  ) +
  
  ## Raw observations used in the model.
  geom_point(
    data = fireflies_model,
    aes(x = lux_factor, y = LRR, colour = type),
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.5,
      seed = 1
    ),
    alpha = 0.25,
    size = 2.0,
    shape = 16
  ) +
  
  ## Model-based 95% CI.
  geom_linerange(
    data = pred_grid,
    aes(
      x = lux_factor,
      ymin = ci_low_LRR,
      ymax = ci_high_LRR,
      colour = type,
      group = type
    ),
    position = pd,
    linewidth = 1.0
  ) +
  
  ## Model-based predictions.
  geom_point(
    data = pred_grid,
    aes(
      x = lux_factor,
      y = fit_LRR,
      colour = type,
      group = type
    ),
    position = pd,
    size = 4.2,
    shape = 16
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_colour_manual(
    values = lamp_cols,
    breaks = c("HPS", "LED"),
    name = "Lamp type"
  ) +
  
  scale_x_discrete(
    limits = c("0.1", "0.5", "1")
  ) +
  
  labs(
    x = "Illuminance level used in the experiment (lux)",
    y = "Log-response ratio (LRR)"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

print(p_flux)

ggsave(
  file.path(out_dir, "Fig_photonflux_LRR_lineRE.tiff"),
  plot = p_flux,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)