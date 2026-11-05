# Run from the project root: source("analysis.R") or Rscript analysis.R
# Data files may sit in the project root or in data/.
# Defines: N, n, f, main, block, nh, nd, pilot, pilot_people, pilot_cars,
#          len_tab, pilot_se15, target_se, S2, n_needed, res_people, res_cars,
#          ps_block, ps_day, day_tab, member_tab, pilot_hours, main_hours.

set.seed(2026)                       # reproducible bootstrap
N <- 120                             # 15-min intervals, 10:00-16:00, Mon-Fri
B <- 20000                           # bootstrap resamples
days <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
datafile <- function(f) if (file.exists(file.path("data", f))) file.path("data", f) else f

# main survey ----------------------------------------------------------------
main <- read.csv(datafile("main_survey.csv"), stringsAsFactors = FALSE)
hm <- do.call(rbind, strsplit(main$time, ":"))
main$time <- as.numeric(hm[, 1]) + as.numeric(hm[, 2]) / 60       # decimal hours
main$people <- main$pedestrians + main$cyclists                    # people = pedestrians + cyclists
day_idx <- match(main$day, days) - 1
stopifnot(!anyNA(day_idx))
main$interval <- day_idx * 24 + (main$time - 10) * 4 + 1           # interval number in 1..N
stopifnot(all(main$interval == round(main$interval)),
          all(main$interval %in% 1:N),
          anyDuplicated(main$interval) == 0)                       # sampled without replacement
main$day <- factor(main$day, levels = days)
n <- nrow(main)
f <- n / N

# pilot survey: people = pedestrians + cyclists, 5-min bins -> 15-min blocks ---
pilot <- read.csv(datafile("pilot_survey.csv"), stringsAsFactors = FALSE)
pilot$people <- pilot$pedestrians + pilot$cyclists
pilot$block15 <- ceiling(pilot$bin / 3)
pilot15 <- aggregate(cbind(people, cars) ~ session + block15, data = pilot, FUN = sum)
pilot_people <- pilot15$people
pilot_cars <- pilot15$cars

# pilot: effect of interval length (sessions are 30 min = 6 five-minute bins) --
bins_per_frame <- N * 3              # 5-min bins in the frame
len_row <- function(len_bins) {
  blk <- (pilot$bin - 1) %/% len_bins
  y <- as.vector(tapply(pilot$people, list(pilot$session, blk), sum))
  Nlen <- bins_per_frame / len_bins; nlen <- length(y); s2 <- var(y)
  c(len = 5 * len_bins, n = nlen, S2 = s2,
    SE = Nlen * sqrt((1 - nlen / Nlen) * s2 / nlen), Nlen = Nlen)
}
len_tab <- as.data.frame(t(sapply(c(1, 2, 3, 6), len_row)))

# target SE = half the pilot SE at 15-min intervals (nearest 10) ---------------
pilot_se15 <- len_tab$SE[len_tab$len == 15]
target_se <- round(pilot_se15 / 2, -1)
len_tab$n_req <- with(len_tab, Nlen^2 * S2 / (target_se^2 + Nlen * S2))   # intervals needed
len_tab$hours <- len_tab$n_req * len_tab$len / 60                          # counting hours needed
S2 <- var(pilot_people)
n_needed <- N^2 * S2 / (target_se^2 + N * S2)

pilot_hours <- 3 * 0.5 * 2           # 3 sessions x 30 min x 2 observers
main_hours <- n * 0.25

# estimation -------------------------------------------------------------------
# Bootstrap resample means are rescaled about ybar so their variance equals the
# design variance (1 - n/N) s^2 / n, i.e. the finite population correction.
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
res_cars <- estimate(main$cars, N, B)

# post-stratification (equal stratum sizes Ng) ---------------------------------
post_strat <- function(y, g, Ng) {
  m <- tapply(y, g, mean); v <- tapply(y, g, var); k <- as.vector(table(g))
  c(T = sum(Ng * m), se = sqrt(sum(Ng^2 * (1 - k / Ng) * v / k)))
}
block <- cut(main$time, c(10, 12, 14, 16), right = FALSE, labels = c("10-12", "12-14", "14-16"))
nh <- as.vector(table(block))                 # intervals per two-hour block
nd <- as.vector(table(main$day))              # intervals per day
ps_block <- post_strat(main$people, block, N / 3)      # 40 intervals per block
ps_day <- post_strat(main$people, main$day, N / 5)     # 24 intervals per day

# who counted what ---------------------------------------------------------------
day_tab <- do.call(rbind, lapply(days, function(d) {
  s <- main[main$day == d, ]; tb <- table(s$member)
  data.frame(Day = d, n = nrow(s), Mean = round(mean(s$people), 1),
             Observers = paste(names(tb), tb, collapse = ", "))
}))
member_tab <- aggregate(people ~ member, data = main, FUN = function(x) c(n = length(x), mean = mean(x)))

if (sys.nframe() == 0) {
  cat("n =", n, " pilot S2 =", round(S2, 1), " target SE =", target_se,
      " n needed =", round(n_needed, 1), "\n")
  print(sapply(list(people = res_people, cars = res_cars), unlist))
  cat("post-stratified (block):", round(ps_block), "  (day):", round(ps_day), "\n")
  print(len_tab); print(day_tab)
}
