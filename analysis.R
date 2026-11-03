# Analysis for the Granny Clark's Wynd crossing survey.
# Run from the project root (source("R/analysis.R") or Rscript R/analysis.R).
# Defines: N, n, main, pilot_walk, pilot_cars, res_people, res_cars.

set.seed(2026)          # for reproducible bootstrap
N <- 120                # 15-min intervals, 10:00-16:00, Mon-Fri
B <- 20000              # bootstrap resamples

# ---- Main survey -----------------------------------------------------------
main <- read.csv("data/main_survey.csv", stringsAsFactors = FALSE)
hm <- do.call(rbind, strsplit(main$time, ":"))
main$time   <- as.numeric(hm[, 1]) + as.numeric(hm[, 2]) / 60   # decimal hours
main$people <- main$pedestrians + main$cyclists                 # people = pedestrians + cyclists
day_idx <- match(main$day, c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")) - 1
main$interval <- day_idx * 24 + (main$time - 10) * 4 + 1        # interval number in 1..N
n <- nrow(main)

# ---- Pilot survey: sum 5-min bins to 15-min blocks ---------------------------
pilot <- read.csv("data/pilot_survey.csv", stringsAsFactors = FALSE)
pilot$block15 <- ceiling(pilot$bin / 3)
pilot15 <- aggregate(cbind(pedestrians, cars) ~ session + block15, data = pilot, FUN = sum)
pilot_walk <- pilot15$pedestrians       # pedestrians only (cyclists/dogs excluded)
pilot_cars <- pilot15$cars

# ---- Sample size from pilot ------------------------------------------------
S2 <- var(pilot_walk)
target_se <- 475
n_needed <- N^2 * S2 / (target_se^2 + N * S2)

# ---- Estimation ------------------------------------------------------------
estimate <- function(y, N, B = 20000) {
  n <- length(y); s2 <- var(y)
  T_hat <- N * mean(y)
  se <- N * sqrt((1 - n / N) * s2 / n)
  tq <- qt(0.975, n - 1)
  boot <- replicate(B, N * mean(sample(y, n, replace = TRUE)))
  b <- quantile(boot, c(0.025, 0.975))
  list(T = T_hat, se = se, rel = se / T_hat,
       t_lo = T_hat - tq * se, t_hi = T_hat + tq * se,
       b_lo = unname(b[1]), b_hi = unname(b[2]))
}
res_people <- estimate(main$people, N, B)
res_cars   <- estimate(main$cars,   N, B)

if (sys.nframe() == 0) {   # printed only when run as a script
  cat("n =", n, " pilot S2 =", round(S2, 1), " n needed =", round(n_needed, 1), "\n")
  print(sapply(list(people = res_people, cars = res_cars), unlist))
}
