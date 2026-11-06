# SRS analysis of the main survey.
# Frame: N = 120 fifteen-minute intervals (10:00-16:00, Mon-Fri).
# Sample: n intervals drawn without replacement; each interval is one element.
# Run from the project root: source("analysis_srs.R") or Rscript analysis_srs.R
# Defines: N, n, f, main, S2, pilot_se15, target_se, n_needed,
#          res_people, res_cars, pilot_hours, main_hours.

set.seed(2026)
N <- 120
B <- 20000
days <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
datafile <- function(f) if (file.exists(file.path("data", f))) file.path("data", f) else f

# main survey ------------------------------------------------------------------
main <- read.csv(datafile("main_survey.csv"), stringsAsFactors = FALSE)
hm <- do.call(rbind, strsplit(main$time, ":"))
main$time <- as.numeric(hm[, 1]) + as.numeric(hm[, 2]) / 60
main$people <- main$pedestrians + main$cyclists
day_idx <- match(main$day, days) - 1
stopifnot(!anyNA(day_idx))
main$interval <- day_idx * 24 + (main$time - 10) * 4 + 1
stopifnot(all(main$interval == round(main$interval)),
          all(main$interval %in% 1:N),
          anyDuplicated(main$interval) == 0)   # sampled without replacement
n <- nrow(main)
f <- n / N

# SRS estimator ----------------------------------------------------------------
# Bootstrap means are rescaled about ybar so their variance equals the design
# variance (1 - n/N) s^2 / n, i.e. includes the finite population correction.
srs_estimate <- function(y, N, B = 20000) {
  n <- length(y); ybar <- mean(y); s2 <- var(y); fpc <- 1 - n / N
  T_hat <- N * ybar
  se <- N * sqrt(fpc * s2 / n)
  tq <- qt(0.975, n - 1)
  bm <- replicate(B, mean(sample(y, n, replace = TRUE)))
  boot <- N * (ybar + sqrt(fpc * n / (n - 1)) * (bm - ybar))
  b <- quantile(boot, c(0.025, 0.975))
  list(mean = ybar, s2 = s2, T = T_hat, se = se, rel = se / T_hat,
       t_lo = T_hat - tq * se, t_hi = T_hat + tq * se,
       b_lo = unname(b[1]), b_hi = unname(b[2]))
}
res_people <- srs_estimate(main$people, N, B)
res_cars <- srs_estimate(main$cars, N, B)

# pilot: variance and required sample size ---------------------------------------
pilot <- read.csv(datafile("pilot_survey.csv"), stringsAsFactors = FALSE)
pilot$people <- pilot$pedestrians + pilot$cyclists
pilot$block15 <- ceiling(pilot$bin / 3)
pilot15 <- aggregate(cbind(people, cars) ~ session + block15, data = pilot, FUN = sum)
pilot_people <- pilot15$people

S2 <- var(pilot_people)
pilot_n <- length(pilot_people)
pilot_se15 <- N * sqrt((1 - pilot_n / N) * S2 / pilot_n)
target_se <- round(pilot_se15 / 2, -1)                  # half the pilot SE
n_needed <- N^2 * S2 / (target_se^2 + N * S2)           # SRS sample size for target SE

pilot_hours <- 3 * 0.5 * 2     # 3 sessions x 30 min x 2 observers
main_hours <- n * 0.25

if (sys.nframe() == 0) {
  cat("n =", n, " f =", round(f, 3), " pilot S2 =", round(S2, 1),
      " target SE =", target_se, " n needed =", round(n_needed, 1), "\n")
  print(sapply(list(people = res_people, cars = res_cars), unlist))
}