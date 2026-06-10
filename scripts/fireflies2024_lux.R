rm(list = ls())

## =========================================================
## Fireflies 2024
##
## Main question:
## How does illuminance and lamp type affect male firefly abundance?
##
## Response:
## LRR = log-response ratio of male abundance relative to the
## locality-specific 0.1 lx reference level.
##
## Important:
## - The main inferential model excludes 0.1 lx, because this level
##   defines the reference used to calculate LRR.
## - The final figure includes 0.1 lx only as visual reference.
## - Trap line is included as a locality-specific random effect.
## =========================================================


## =========================================================
## 1. Libraries
## =========================================================

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)


## =========================================================
## 2. User settings
## =========================================================

data_file <- "data/fireflies2024.csv"
out_dir   <- "figures"

reference_lux <- 0.1
eps <- 0.5

lamp_cols <- c(
  LED = "#0070c0",
  HPS = "#ffc000"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## =========================================================
## 3. Load and clean data
## =========================================================

fireflies_raw <- read.csv2(
  data_file,
  header = TRUE,
  sep = ";",
  dec = "."
)

## Remove biodynamic treatment.
## This analysis compares only HPS and LED.
fireflies <- fireflies_raw %>%
  filter(type != "biodyn") %>%
  mutate(
    lux    = as.numeric(lux),
    splend = as.numeric(splend),
    
    ## Recode lamp type for plotting and interpretation.
    ## HPS is used instead of "sodium" for consistency with the manuscript.
    type = case_when(
      type %in% c("sodium", "hps", "HPS") ~ "HPS",
      type %in% c("LED", "led")           ~ "LED",
      TRUE                                ~ as.character(type)
    ),
    
    ## Set HPS as the reference level in models.
    type = factor(type, levels = c("HPS", "LED")),
    
    region   = factor(region),
    locality = factor(locality),
    side     = factor(side)
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
    ## This is equivalent to locality:side.
    line_id = interaction(locality, side, drop = TRUE)
  )

## log2(lux) requires strictly positive values.
if (any(fireflies$lux <= 0, na.rm = TRUE)) {
  stop("lux contains zero or negative values. log2(lux) requires strictly positive values.")
}

## Basic data checks.
cat("=====================================================\n")
cat("DATA CHECKS\n")
cat("=====================================================\n")
cat("Total males:", sum(fireflies$splend, na.rm = TRUE), "\n")
cat("Number of observations:", nrow(fireflies), "\n\n")

cat("Observations by lamp type and illuminance:\n")
print(xtabs(~ type + lux, data = fireflies))
cat("\n")

cat("Observations by locality-side line and illuminance:\n")
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
## 4. Compute LRR relative to 0.1 lx
## =========================================================
## For each locality, the mean abundance at 0.1 lx is used as the
## local reference value.
##
## eps is added to avoid division by zero and log(0).

fireflies <- fireflies %>%
  group_by(locality) %>%
  mutate(
    baseline = mean(splend[lux == reference_lux], na.rm = TRUE),
    rel_abund = (splend + eps) / (baseline + eps),
    LRR = log(rel_abund),
    log2_lux = log2(lux),
    lux_factor = factor(
      lux,
      levels = c(0.1, 0.5, 1),
      labels = c("0.1", "0.5", "1")
    )
  ) %>%
  ungroup()

## Check that every locality has a usable 0.1 lx baseline.
if (any(is.na(fireflies$baseline))) {
  bad_localities <- fireflies %>%
    filter(is.na(baseline)) %>%
    distinct(locality)
  
  print(bad_localities)
  stop("Some localities do not have a valid 0.1 lx baseline.")
}

cat("Baseline per locality:\n")
baseline_check <- fireflies %>%
  distinct(locality, baseline) %>%
  arrange(locality)

print(baseline_check)
cat("\n")

cat("Baseline observations by locality and lamp type:\n")
baseline_type_check <- fireflies %>%
  filter(lux == reference_lux) %>%
  group_by(locality, type) %>%
  summarise(
    n = n(),
    baseline_type_mean = mean(splend, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(locality, type)

print(baseline_type_check)
cat("\n")

cat("Lux values and their log2-transformed values:\n")
lux_log2_table <- fireflies %>%
  distinct(lux, log2_lux) %>%
  arrange(lux)

print(lux_log2_table)
cat("\n")


## =========================================================
## 5. Main inferential dataset
## =========================================================
## The 0.1 lx level is excluded from the main model because it defines
## the LRR reference. The model therefore evaluates the response at
## 0.5 and 1 lx relative to the 0.1 lx baseline.

fireflies_model <- fireflies %>%
  filter(lux != reference_lux) %>%
  droplevels()

cat("Number of observations used in the main model:", nrow(fireflies_model), "\n\n")

cat("Model observations by lamp type and illuminance:\n")
print(xtabs(~ type + lux, data = fireflies_model))
cat("\n")


## =========================================================
## 6. Main model selection
## =========================================================
## Candidate models use the same random-effect structure:
##
##   (1 | locality) + (1 | line_id)
##
## where line_id is the locality-by-side trap-line identity.
##
## mod_full:
##   Tests whether the slope of illuminance differs between HPS and LED.
##
## mod_add:
##   Additive model: independent effects of illuminance and lamp type.
##
## mod_lux_only:
##   Illuminance-only model, used to check whether lamp type improves
##   the model.

mod_full <- glmmTMB(
  LRR ~ log2_lux * type + (1 | locality) + (1 | line_id),
  data = fireflies_model
)

mod_add <- glmmTMB(
  LRR ~ log2_lux + type + (1 | locality) + (1 | line_id),
  data = fireflies_model
)

mod_lux_only <- glmmTMB(
  LRR ~ log2_lux + (1 | locality) + (1 | line_id),
  data = fireflies_model
)

## Final model:
## The additive model is used when the interaction is not supported.
mod_final <- mod_add

cat("=====================================================\n")
cat("MODEL COMPARISON\n")
cat("=====================================================\n")

cat("AIC comparison:\n")
print(AIC(mod_lux_only, mod_add, mod_full))
cat("\n")

cat("Likelihood-ratio test: additive vs interaction model\n")
print(anova(mod_add, mod_full))
cat("\n")

cat("Likelihood-ratio test: lux-only vs additive model\n")
print(anova(mod_lux_only, mod_add))
cat("\n")


## =========================================================
## 7. Final model summary
## =========================================================

cat("=====================================================\n")
cat("FINAL MODEL SUMMARY\n")
cat("=====================================================\n")
print(summary(mod_final))
cat("\n")


## =========================================================
## 8. Fixed-effect confidence intervals and effect sizes
## =========================================================
## Effects are converted from the LRR scale to multiplicative ratios.
##
## Example:
## ratio = exp(beta)
##
## ratio < 1 means lower relative abundance.
## ratio = 0.57 means abundance is 57% of the comparison level.

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

## Interpret illuminance effect.
b <- fixef(mod_final)$cond

beta_lux <- unname(b["log2_lux"])
lux_ratio <- exp(beta_lux)

cat("Effect of illuminance:\n")
cat("- beta for log2(lux):", beta_lux, "\n")
cat("- multiplicative change in relative abundance for each doubling of lux:",
    lux_ratio, "\n")
cat("- percent change for each doubling of lux:",
    round((lux_ratio - 1) * 100, 1), "%\n\n")

## Interpret LED effect relative to HPS.
if ("typeLED" %in% names(b)) {
  beta_LED <- unname(b["typeLED"])
  LED_ratio <- exp(beta_LED)
  
  cat("Effect of lamp type, LED vs HPS at the same log2(lux):\n")
  cat("- beta for typeLED:", beta_LED, "\n")
  cat("- multiplicative ratio in relative abundance:", LED_ratio, "\n")
  cat("- percent difference of LED relative to HPS:",
      round((LED_ratio - 1) * 100, 1), "%\n\n")
}


## =========================================================
## 9. Sensitivity to pseudo-count eps
## =========================================================
## This is not the main analysis.
## It checks whether conclusions depend strongly on the chosen eps.
##
## The same random-effect structure is used as in the final model.

run_eps_sensitivity <- function(eps_value, data_input) {
  
  tmp <- data_input %>%
    group_by(locality) %>%
    mutate(
      baseline_tmp = mean(splend[lux == reference_lux], na.rm = TRUE),
      LRR_tmp = log((splend + eps_value) / (baseline_tmp + eps_value)),
      log2_lux_tmp = log2(lux)
    ) %>%
    ungroup() %>%
    filter(lux != reference_lux) %>%
    droplevels()
  
  m <- glmmTMB(
    LRR_tmp ~ log2_lux_tmp + type + (1 | locality) + (1 | line_id),
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
## 10. Diagnostics for the final model
## =========================================================
## DHARMa diagnostics are simulation-based and are used here to check:
## - residual uniformity
## - dispersion
## - outliers
## - residual patterns against predictors and grouping variables

sim_final <- simulateResiduals(mod_final, n = 1000)

plot(sim_final)

cat("DHARMa diagnostics for final model:\n")
print(testUniformity(sim_final))
print(testDispersion(sim_final))
print(testOutliers(sim_final))
print(testQuantiles(sim_final))
cat("\n")

## Residual checks against key predictors and grouping variables.
plotResiduals(sim_final, form = fireflies_model$log2_lux)
plotResiduals(sim_final, form = fireflies_model$type)
plotResiduals(sim_final, form = fireflies_model$locality)
plotResiduals(sim_final, form = fireflies_model$line_id)


## =========================================================
## 11. Prediction grid for the final publication plot
## =========================================================
## Predictions are based on the final inferential model:
##
##   LRR ~ log2_lux + type + (1 | locality) + (1 | line_id)
##
## The 0.1 lx level is not predicted because it defines the
## reference used to calculate LRR and was excluded from the
## final model.
##
## The 0.1 lx observations are shown only as visual reference.
## Population-level predictions are shown, with random effects excluded.

pred_grid <- expand.grid(
  lux = sort(unique(fireflies_model$lux)),
  type = levels(fireflies_model$type)
) %>%
  mutate(
    log2_lux = log2(lux),
    type = factor(type, levels = levels(fireflies_model$type)),
    
    ## Required by the model formula.
    ## Ignored in prediction because re.form = NA.
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


## =========================================================
## 12. Final publication plot
## =========================================================
## Figure style:
## - grey open points = 0.1 lx reference observations
## - small coloured points = raw observations used in the final model
## - large coloured points = population-level predictions from the final model
## - vertical lines = 95% confidence intervals
## - dashed horizontal line = no change relative to 0.1 lx baseline

pd <- position_dodge(width = 0.5)

p_2024_LRR <- ggplot() +
  
  ## Reference observations at 0.1 lx.
  ## These points are shown only to make the reference level visible.
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
  
  ## Raw observations used in the final model.
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
  
  ## Model-based 95% confidence intervals from the final model.
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
  
  ## Model-based predictions from the final model.
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
  
  ## Reference line: no change relative to 0.1 lx.
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_colour_manual(
    values = lamp_cols,
    breaks = c("LED", "HPS"),
    name = "Lamp type"
  ) +
  
  scale_x_discrete(
    limits = c("0.1", "0.5", "1")
  ) +
  
  scale_y_continuous(
    breaks = c(-4, -2, 0, 2),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  labs(
    x = "Light intensity (lux)",
    y = "Log-response ratio (LRR)"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

print(p_2024_LRR)


## =========================================================
## 13. Export figure
## =========================================================

ggsave(
  file.path(out_dir, "Fig5.tiff"),
  plot = p_2024_LRR,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)