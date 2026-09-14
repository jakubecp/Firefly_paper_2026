rm(list = ls())

## =========================================================
## Fireflies 2023
## Main analysis = photopic illuminance in lux
## Response = whether at least one male found the female
## Trap line is included as a design-level random effect
## =========================================================

library(dplyr)
library(glmmTMB)
library(DHARMa)
library(ggplot2)

## ---------------------------------------------------------
## Load data
## ---------------------------------------------------------
data <- read.csv2("data/fireflies2023.csv", sep = ";", dec = ".")

## ---------------------------------------------------------
## Basic cleaning
## ---------------------------------------------------------
data$date0 <- as.Date(data$date, format = "%m.%d.%Y")

if (any(is.na(data$date0))) {
  stop("Some dates could not be parsed. Check the format in data$date.")
}

## Keep only localities A and B
data2 <- data %>%
  filter(loc %in% c("A", "B")) %>%
  filter(!is.na(males), !is.na(lux_f), !is.na(date0), !is.na(side)) %>%
  mutate(
    loc = factor(loc, levels = c("A", "B")),
    side = factor(side),
    date0 = as.Date(date0)
  ) %>%
  droplevels()

## Check that lux can be log-transformed
if (any(data2$lux_f <= 0, na.rm = TRUE)) {
  stop("lux_f contains zero or negative values. log2(lux_f) requires strictly positive values.")
}

## ---------------------------------------------------------
## Derived variables
## ---------------------------------------------------------
date_levels <- sort(unique(data2$date0))

data2 <- data2 %>%
  mutate(
    date_f = factor(date0, levels = date_levels),
    found = as.integer(males > 0),
    log2_lux = log2(lux_f),
    line_id = interaction(loc, side, drop = TRUE)
  )

## ---------------------------------------------------------
## Basic checks
## ---------------------------------------------------------
cat("Number of observations:", nrow(data2), "\n")
cat("Total males:", sum(data2$males, na.rm = TRUE), "\n")
cat("Number of successful observations (found = 1):",
    sum(data2$found, na.rm = TRUE), "\n\n")

cat("Successes by lux and locality:\n")
print(
  data2 %>%
    group_by(loc, lux_f) %>%
    summarise(
      n = n(),
      successes = sum(found),
      failures = n - successes,
      prop_success = mean(found),
      .groups = "drop"
    )
)
cat("\n")

## =========================================================
## MAIN MODEL
## Probability that at least one male found the female
## Predictor: photopic illuminance in lux
## Random effects: sampling date and trap line
## =========================================================

mod_prob <- glmmTMB(
  found ~ log2_lux + loc + (1 | date_f) + (1 | line_id),
  family = binomial(link = "logit"),
  data = data2
)

cat("=====================================================\n")
cat("MAIN MODEL: probability of finding the female\n")
cat("Predictor: photopic illuminance in lux\n")
cat("Trap line included as a design-level random intercept\n")
cat("=====================================================\n")

mod_sum <- summary(mod_prob)

cat("Fixed effects:\n")
print(mod_sum$coefficients$cond)
cat("\n")

## ---------------------------------------------------------
## Confidence intervals for fixed effects
## ---------------------------------------------------------
coef_tab <- mod_sum$coefficients$cond

ci_fixed <- data.frame(
  term = rownames(coef_tab),
  estimate = coef_tab[, "Estimate"],
  se = coef_tab[, "Std. Error"],
  lower_95 = coef_tab[, "Estimate"] - 1.96 * coef_tab[, "Std. Error"],
  upper_95 = coef_tab[, "Estimate"] + 1.96 * coef_tab[, "Std. Error"],
  row.names = NULL
)

cat("95% Wald CI for fixed effects:\n")
print(ci_fixed)
cat("\n")

## ---------------------------------------------------------
## Effect size: odds ratio per doubling of lux
## ---------------------------------------------------------
lux_effect <- ci_fixed %>%
  filter(term == "log2_lux") %>%
  mutate(
    OR_per_doubling_lux = exp(estimate),
    CI_low = exp(lower_95),
    CI_high = exp(upper_95)
  ) %>%
  select(term, OR_per_doubling_lux, CI_low, CI_high)

cat("Odds ratio for each doubling of lux:\n")
print(lux_effect)
cat("\n")

## =========================================================
## Diagnostics
## =========================================================

sim_prob <- simulateResiduals(mod_prob, n = 1000)

plot(sim_prob)
print(testUniformity(sim_prob))
print(testDispersion(sim_prob))
print(testOutliers(sim_prob))

## Temporal autocorrelation check
## Residuals are aggregated by date first, then tested in chronological order
sim_prob_date <- recalculateResiduals(sim_prob, group = data2$date_f)

cat("Temporal autocorrelation test, residuals aggregated by date:\n")
print(testTemporalAutocorrelation(
  sim_prob_date,
  time = as.numeric(date_levels)
))
cat("\n")

## Optional residual check against the main continuous predictor
plotResiduals(sim_prob, form = data2$lux_f)

## =========================================================
## Estimated lux thresholds
## Lux at which P(found) drops below selected probabilities
## Population-level estimates; random effects excluded
## =========================================================

coef_prob <- fixef(mod_prob)$cond

lux_at_p <- function(p, intercept, slope) {
  if (is.na(slope) || slope >= 0) return(NA_real_)
  2 ^ ((qlogis(p) - intercept) / slope)
}

int_A <- unname(coef_prob["(Intercept)"])
slp <- unname(coef_prob["log2_lux"])

locB_effect <- if ("locB" %in% names(coef_prob)) {
  unname(coef_prob["locB"])
} else {
  0
}

int_B <- int_A + locB_effect

lux_range <- range(data2$lux_f, na.rm = TRUE)

thresholds <- data.frame(
  locality = c("A", "B"),
  lux_at_P_0.05 = c(
    lux_at_p(0.05, int_A, slp),
    lux_at_p(0.05, int_B, slp)
  ),
  lux_at_P_0.01 = c(
    lux_at_p(0.01, int_A, slp),
    lux_at_p(0.01, int_B, slp)
  )
)

thresholds <- thresholds %>%
  mutate(
    P_0.05_inside_observed_range =
      lux_at_P_0.05 >= lux_range[1] & lux_at_P_0.05 <= lux_range[2],
    P_0.01_inside_observed_range =
      lux_at_P_0.01 >= lux_range[1] & lux_at_P_0.01 <= lux_range[2]
  )

cat("Estimated thresholds from the lux model:\n")
print(thresholds)
cat("\n")

cat("Observed lux range:", lux_range[1], "to", lux_range[2], "\n\n")

## =========================================================
## Graph
## Observed proportions + fitted probabilities
## x-axis = lux
## Trap line is included in the model but not displayed
## Population-level predictions, random effects excluded
## =========================================================

obs_prob <- data2 %>%
  group_by(loc, lux_f) %>%
  summarise(
    n = n(),
    successes = sum(found),
    prop = mean(found),
    .groups = "drop"
  )

## Exact binomial CI for observed proportions
obs_prob$ci_low <- NA_real_
obs_prob$ci_high <- NA_real_

for (i in seq_len(nrow(obs_prob))) {
  bt <- binom.test(obs_prob$successes[i], obs_prob$n[i])
  obs_prob$ci_low[i] <- bt$conf.int[1]
  obs_prob$ci_high[i] <- bt$conf.int[2]
}

## Prediction grid
pred_grid <- expand.grid(
  lux_f = exp(seq(
    log(min(data2$lux_f)),
    log(max(data2$lux_f)),
    length.out = 200
  )),
  loc = levels(data2$loc)
)

pred_grid <- pred_grid %>%
  mutate(
    log2_lux = log2(lux_f),
    loc = factor(loc, levels = levels(data2$loc)),
    date_f = data2$date_f[1],
    line_id = data2$line_id[1]
  )

pred_link <- predict(
  mod_prob,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

pred_grid$fit <- plogis(pred_link$fit)
pred_grid$ci_low <- plogis(pred_link$fit - 1.96 * pred_link$se.fit)
pred_grid$ci_high <- plogis(pred_link$fit + 1.96 * pred_link$se.fit)

p1 <- ggplot() +
  geom_ribbon(
    data = pred_grid,
    aes(x = lux_f, ymin = ci_low, ymax = ci_high),
    alpha = 0.2
  ) +
  geom_line(
    data = pred_grid,
    aes(x = lux_f, y = fit),
    linewidth = 1
  ) +
  geom_point(
    data = obs_prob,
    aes(x = lux_f, y = prop),
    size = 2
  ) +
  geom_linerange(
    data = obs_prob,
    aes(x = lux_f, ymin = ci_low, ymax = ci_high)
  ) +
  scale_x_log10(
    labels = scales::label_number()
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(
    ~ loc,
    labeller = as_labeller(c(
      A = "Locality A",
      B = "Locality B"
    ))
  ) +
  labs(
    x = "Illuminance (lux)",
    y = "Probability of success"
  ) +
  theme_bw()

print(p1)

dir.create("figures", showWarnings = FALSE)

ggsave(
  "figures/Fig7.tiff",
  plot = p1,
  width = 5,
  height = 3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

## =========================================================
## Notes:
## - The response is binary: whether at least one male found the female.
## - The lux effect is interpreted per doubling of illuminance.
## - Trap line is included as a design-level random effect.
## - Thresholds are model-based estimates and should not be overinterpreted
##   if they fall outside the observed lux range.
## =========================================================