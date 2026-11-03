# MT4608 Sampling Theory Group Project
## Euan Allison, Atton Cheung, Sean Gallacher, Ben Porksen-manning, Sofia Soderberg

Estimate of weekday (Mon-Fri, 10:00-16:00) crossings of Granny Clark's Wynd, St Andrews, from a simple random sample of 21 of 120 fifteen-minute intervals.

## Structure
- `data/main_survey.csv` - main survey counts (pedestrians, cyclists, cars) per sampled interval.
- `data/pilot_survey.csv` - pilot counts in 5-minute bins (3 sessions).
- `R/analysis.R` - loads data, computes the estimates, standard errors, t and bootstrap CIs.
- `report.Rmd` - the report; sources `R/analysis.R`.

## Reproduce
From the project root:
```r
source("R/analysis.R")
rmarkdown::render("report.Rmd")        # requires rmarkdown and LaTeX
```
or `Rscript R/analysis.R` to print the estimates. The bootstrap uses `set.seed(2026)`.
