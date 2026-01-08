rm(list = ls())

## =========================================================
## Male abundance vs light level 2023
## =========================================================

library(dplyr)      # data wrangling
library(glmmTMB)    # GLMM fitting
library(DHARMa)     # residual diagnostics for GLMMs

## ---------------------------------------------------------
## Load data
## ---------------------------------------------------------
data <- read.csv2("data/fireflies2023.csv", sep = ";", dec = ".")

# Parse date/time (NOTE: this is not used in the model below; keep only if needed later)
data$date <- as.POSIXct(data$date, format = "%m.%d.%Y")

# Convert grouping variables to factors
data$loc  <- factor(data$loc)
data$side <- factor(data$side)

## ---------------------------------------------------------
## Filter data: keep only localities A and B
## ---------------------------------------------------------
data2 <- subset(data, loc != "C")

## ---------------------------------------------------------
## Optional: create trap-level ID (loc + row)
## (Currently NOT used in the fitted model; keep only if you plan
##  to use (1 | id2) later.)
## ---------------------------------------------------------
data2$row <- factor(data2$row)
data2$id2 <- factor(paste(data2$loc, data2$row))

# Quick check: total number of males in the analysed subset
sum(data2$males, na.rm = TRUE)

## ---------------------------------------------------------
## Model: males ~ lux_f with locality random intercept
## ---------------------------------------------------------
mod1 <- glmmTMB(
  males ~ lux_f + (1 | loc),
  family = poisson,
  data = data2
)

## ---------------------------------------------------------
## Diagnostics (DHARMa)
## ---------------------------------------------------------
sim_res_f <- simulateResiduals(mod1)
plot(sim_res_f)                 # overall residual diagnostics
testZeroInflation(sim_res_f)    # test for excess zeros vs model expectation

## ---------------------------------------------------------
## Model fit and summary output
## ---------------------------------------------------------
summary(mod1)

