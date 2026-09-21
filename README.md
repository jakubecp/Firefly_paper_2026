# Haunted by the light: Effects of street lamps on green ghost firefly *Lamprohiza splendidula*

This repository contains data and R scripts used to analyse the effects of artificial light at night (ALAN) on the courtship behavior of the green ghost firefly *Lamprohiza splendidula*.

The study evaluates the impact of three commonly used street lighting technologies—high-pressure sodium lamps (HPS), warm white LEDs, and biodynamic lamps in PC-amber mode—on male mate-searching activity in both urban and previously unilluminated (“pristine”) habitats. In addition, short-term persistence effects following experimental illumination are assessed.

## Associated manuscript

Haunted by the light: Effects of street lamps on green ghost firefly *Lamprohiza splendidula*  
Martin Novák, Pavel Jakubec, Jana Svobodová, Petr Chajma, Tomáš Kadlec  
Faculty of Environmental Sciences, Czech University of Life Sciences Prague

## Repository structure

data/
  fireflies2023.csv   # HPS lamps – initial experiment
  
  fireflies2024.csv   # HPS vs warm LED – main experiment
  
  fireflies2025.csv   # Biodynamic lamps – pristine habitats
  
  weather_hourly.xlsx # hourly weather data for all years Pisel (2023)
  
  female_observations_2023.xlsx   # Observations of the females in 2023

scripts/
  fireflies2023_lux.R     # Analysis of the 2023 data
  
  fireflies2024_lux.R     # Analysis of the 2024 data
  
  biodynamic_lux.R        # Analysis of the 2025 data
  
  00_run_all.R        # Master script running the full analysis

figures/
  Fig4.tiff
  Fig5.tiff
  ...

renv.lock             # Locked R package versions
.Rprofile             # Project-specific R environment

## Data description

The datasets consist of pitfall-trap counts of male Lamprohiza splendidula attracted to artificial glow lures under different light intensities and lamp types.

2023: Initial experiment under isolated HPS lamps  
2024: Expanded experiment comparing HPS and warm LED lamps across multiple localities  
2025: Experimental introduction of biodynamic lighting into previously dark natural habitats  

Each row represents a sampling unit (typically a trap-night, or—in the 2024 dataset—the content of one trap over a longer sampling interval).

## Software and reproducibility

R version 4.5.2 (2025-10-31, ucrt)  
Windows 10 x64 (build 19045)

Exact package versions are recorded in renv.lock.

Key packages:
glmmTMB, lmerTest, emmeans, DHARMa, dplyr, ggplot2, sjPlot, ggbreak, multcomp, multcompView

## Reproducing the analysis

install.packages("renv")
renv::restore()
source("scripts/00_run_all.R")

## Data availability

All data are included in this repository.  
Archived releases will be deposited with a DOI (e.g. Zenodo).

## Funding

This work was supported by the Technology Agency of the Czech Republic  
(grant no. SS06010373).

## References

Pisel T (2023). openmeteo: Retrieve Weather Data from the Open-Meteo API. R package version 0.2.4, https://github.com/cran/openmeteo. (accessed 21 Sep. 2026) 
