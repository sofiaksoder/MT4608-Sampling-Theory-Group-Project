# Run from the project root (source("R/analysis.R") or Rscript R/analysis.R).
# Defines: N, n, f, main, block, nh, pilot_walk, pilot_cars, res_people, res_cars,
#          ps_people, len_tab.

set.seed(2026)          # for reproducible bootstrap
N <- 120                # 15-min intervals, 10:00-16:00, Mon-Fri
B <- 20000              # bootstrap resamples

# main survey ----------
main <- read.csv("data/main_survey.csv", stringsAsFactors = FALSE)
hm <- do.call(rbind, strsplit(main$time, ":"))
main$time   <- as.numeric(hm[, 1]) + as.numeric(hm[, 2]) / 60   # decimal hours
main$people <- main$pedestrians + main$cyclists                 # people = pedestrians + cyclists
day_idx <- match(main$day, c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")) - 1
main$interval <- day_idx * 24 + (main$time - 10) * 4 + 1        # interval number in 1..N
n <- nrow(main)
f <- n / N
stopifnot(anyDuplicated(main$interval) == 0)   # sampled without replacement

# pilot survey: sum 5-min bins to 15-min blocks ----------
pilot <- read.csv("data/pilot_survey.csv", stringsAsFactors = FALSE)
pilot$block15 <- ceiling(pilot$bin / 3)
pilot15 <- aggregate(cbind(pedestrians, cars) ~ session + block15, data = pilot, FUN = sum)
pilot_walk <- pilot15$pedestrians
pilot_cars <- pilot15$cars

# sample size from pilot ----------
S2 <- var(pilot_walk)
target_se <- 475
n_needed <- N^2 * S2 / (target_se^2 + N * S2)

# estimation ----------
# Bootstrap resample means are rescaled about ybar so that their variance equals
# the design variance (1 - n/N) s^2 / n, i.e. the finite population correction
# is included (a plain bootstrap would overstate the SE by about 10% here).
estimate <- function(y, N, B = 20000) {
  n <- length(y); s2 <- var(y); ybar <- mean(y)
  T_hat <- N * ybar
  se <- N * sqrt((1 - n / N) * s2 / n)
  tq <- qt(0.975, n - 1)
  bm <- replicate(B, mean(sample(y, n, replace = TRUE)))
  boot <- N * (ybar + sqrt((1 - n / N) * n / (n - 1)) * (bm - ybar))
  b <- quantile(boot, c(0.025, 0.975))
  list(T = T_hat, se = se, rel = se / T_hat,
       t_lo = T_hat - tq * se, t_hi = T_hat + tq * se,
       b_lo = unname(b[1]), b_hi = unname(b[2]))
}
res_people <- estimate(main$people, N, B)
res_cars   <- estimate(main$cars,   N, B)

# post-stratification by two-hour block (40 intervals per block) ----------
block <- cut(main$time, c(10, 12, 14, 16), right = FALSE, labels = c("10-12", "12-14", "14-16"))
Nh <- rep(N / 3, 3)
nh <- as.vector(table(block))
post_strat <- function(y) {
  m <- tapply(y, block, mean); v <- tapply(y, block, var)
  c(T = sum(Nh * m), se = sqrt(sum(Nh^2 * (1 - nh / Nh) * v / nh)))
}
ps_people <- post_strat(main$people)

# pilot: SE of the estimated total for different interval lengths ----------
# each session is 30 min = 6 five-minute bins; lengths of 5, 10, 15, 30 min.
pilot_se <- function(len_bins) {
  blk <- (pilot$bin - 1) %/% len_bins
  y <- as.vector(tapply(pilot$pedestrians, list(pilot$session, blk), sum))
  Nlen <- 360 / len_bins          # intervals of this length in the 30-hour frame
  nlen <- length(y)
  c(len = 5 * len_bins, n = nlen, S2 = var(y),
    SE = Nlen * sqrt((1 - nlen / Nlen) * var(y) / nlen))
}
len_tab <- as.data.frame(t(sapply(c(1, 2, 3, 6), pilot_se)))

if (sys.nframe() == 0) {
  cat("n =", n, " pilot S2 =", round(S2, 1), " n needed =", round(n_needed, 1), "\n")
  print(sapply(list(people = res_people, cars = res_cars), unlist))
  cat("post-stratified people:", round(ps_people), "\n")
  print(len_tab)
}
