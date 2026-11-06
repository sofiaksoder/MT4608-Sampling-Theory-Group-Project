# MT4608 Sampling Theory Group Project
## Euan Allison, Atton Cheung, Sean Gallacher, Ben Porksen-manning, Sofia Soderberg

Estimate of weekday (Mon-Fri, 10:00-16:00) crossings of Granny Clark's Wynd, St Andrews, from a simple random sample of 21 of 120 fifteen-minute intervals.

## Structure
- `data/main_survey.csv` - main survey counts (pedestrians, cyclists, cars) per sampled interval.
- `data/pilot_survey.csv` - pilot counts in 5-minute bins (3 sessions).
- `analysis_srs.R` - loads data, computes the estimates, standard errors, t and bootstrap CIs.
