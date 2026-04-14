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

foton <- tibble(lux = c(0.05,0.1,0.2,0.5,1,2,4,8), fotonflowHPS = c(0.00065936, 
                                                                     0.001318719,
                                                                     0.002637439,
                                                                     0.006593597,
                                                                     0.013187194,
                                                                     0.026374387,
                                                                     0.052748775,
                                                                     0.105497549
                                                                     ))

data2$fotonHPS[data2$lux_f == 0.05] <- foton$fotonflowHPS[1] 
data2$fotonHPS[data2$lux_f == 0.1] <- foton$fotonflowHPS[2] 
data2$fotonHPS[data2$lux_f == 0.2] <- foton$fotonflowHPS[3] 
data2$fotonHPS[data2$lux_f == 0.5] <- foton$fotonflowHPS[4] 
data2$fotonHPS[data2$lux_f == 1] <- foton$fotonflowHPS[5] 
data2$fotonHPS[data2$lux_f == 2] <- foton$fotonflowHPS[6] 
data2$fotonHPS[data2$lux_f == 4] <- foton$fotonflowHPS[7] 
data2$fotonHPS[data2$lux_f == 8] <- foton$fotonflowHPS[8] 



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

