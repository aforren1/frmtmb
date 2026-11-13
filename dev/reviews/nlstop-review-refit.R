# Review of lane nlstop: does prior_conflict_escape() run on refit() (the
# parametric bootstrap's engine, start = NULL with a warm template), and
# what does it cost there? A bernoulli GLMM, 5000 rows, 100 groups, 5
# slopes at 2 under normal(0, 0.25). Arms interleaved: "on" the lane,
# "off" the escape as the identity. 3 rounds of 3 refits per arm.
# Usage: Rscript dev/reviews/nlstop-review-refit.R
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
orig <- ns$prior_conflict_escape
cnt <- new.env(); cnt$k <- 0L
set_arm <- function(arm) {
  utils::assignInNamespace("prior_conflict_escape",
    if (arm == "on") function(...) { cnt$k <- cnt$k + 1L; orig(...) } else
      function(obj, opt, ...) opt, "frmtmb")
}
set.seed(14); n <- 5000; ng <- 100
d <- data.frame(g = factor(sample(ng, n, TRUE)))
for (k in 1:5) d[[paste0("x", k)]] <- rnorm(n)
u <- rnorm(ng, 0, 0.7)
d$y <- rbinom(n, 1, plogis(-0.5 + 2 * (d$x1 + d$x2 + d$x3 + d$x4 + d$x5) +
                            u[d$g]))
set_arm("on")
fit <- frm(y ~ x1 + x2 + x3 + x4 + x5 + (1 | g), family = bernoulli(),
           data = d, prior = set_prior("normal(0, 0.25)", class = "b"))
sims <- simulate(fit, nsim = 3, seed = 1)
cat("ROWHEAD round arm secs calls ll\n")
for (rd in 1:3) for (arm in if (rd %% 2) c("on", "off") else c("off", "on")) {
  set_arm(arm)
  cnt$k <- 0L
  t0 <- proc.time()[["elapsed"]]
  ll <- vapply(seq_len(ncol(sims)), function(i) {
    as.numeric(logLik_penalized(suppressWarnings(refit(fit, sims[[i]]))))
  }, 0)
  cat(sprintf("ROW %d %s %.2f %d %s\n", rd, arm,
              proc.time()[["elapsed"]] - t0, cnt$k,
              paste(sprintf("%.6f", ll), collapse = ",")))
  flush.console()
}
cat("DONE\n")
