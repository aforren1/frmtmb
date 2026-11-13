# Lane nlstop, fix round: the review's cost cases for the prior restart
# (dev/reviews/nlstop-review-prior.R a4 and a4g, -refit.R), arms rotated
# per round in one process: "on" the lane's build, "off" and "off2" (the
# control) with prior_conflict_escape() as the identity, 0.81.0's path.
# Min of 3 rounds; also whether the restart ran (a counter on the escape)
# and the penalized log-likelihood of every arm.
# Usage: NLSTOP_ARM=after Rscript dev/nlstop-fix-cost.R
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
orig <- ns$prior_conflict_escape
cnt <- new.env(); cnt$k <- 0L; cnt$won <- 0L
set_arm <- function(arm) {
  utils::assignInNamespace("prior_conflict_escape",
    if (arm == "on") function(...) {
      cnt$k <- cnt$k + 1L
      o <- orig(...)
      if (!is.null(o$prior_escape)) cnt$won <- cnt$won + 1L
      o
    } else function(obj, opt, ...) opt, "frmtmb")
}
set.seed(14); n <- 20000; ng <- 200
d <- data.frame(g = factor(sample(ng, n, TRUE)))
for (k in 1:5) d[[paste0("x", k)]] <- rnorm(n)
u <- rnorm(ng, 0, 0.7)
d$y <- rbinom(n, 1, plogis(-0.5 + 2 * (d$x1 + d$x2 + d$x3 + d$x4 + d$x5) +
                            u[d$g]))
d$yg <- 1 + 5 * d$x1 + u[d$g] + rnorm(n)
set.seed(14); n2 <- 5000; ng2 <- 100
d2 <- data.frame(g = factor(sample(ng2, n2, TRUE)))
for (k in 1:5) d2[[paste0("x", k)]] <- rnorm(n2)
u2 <- rnorm(ng2, 0, 0.7)
d2$y <- rbinom(n2, 1, plogis(-0.5 + 2 * (d2$x1 + d2$x2 + d2$x3 + d2$x4 +
                                          d2$x5) + u2[d2$g]))
pr <- set_prior("normal(0, 0.25)", class = "b")
set_arm("on")
fit2 <- frm(y ~ x1 + x2 + x3 + x4 + x5 + (1 | g), family = bernoulli(),
            data = d2, prior = pr)
sims <- simulate(fit2, nsim = 3, seed = 1)
cases <- list(
  a4 = function() frm(y ~ x1 + x2 + x3 + x4 + x5 + (1 | g),
                      family = bernoulli(), data = d, prior = pr),
  a4g = function() frm(yg ~ x1 + x2 + (1 | g), data = d,
                       prior = set_prior("normal(0, 1)", class = "b")),
  refit = function() {
    lapply(seq_len(ncol(sims)), function(i) {
      suppressWarnings(refit(fit2, sims[[i]]))
    })
  })
llp_of <- function(r) {
  if (inherits(r, "frmtmb_fit")) r <- list(r)
  paste(sprintf("%.6f", vapply(r, function(f) {
    as.numeric(logLik_penalized(f))
  }, 0)), collapse = ",")
}
arms <- c("off", "on", "off2")
cat("ROWHEAD case round arm secs calls won llp\n")
X <- list()
for (cn in names(cases)) for (rd in 1:3) {
  for (a in arms[(seq_along(arms) + rd - 2L) %% 3L + 1L]) {
    set_arm(sub("2$", "", a))
    cnt$k <- 0L; cnt$won <- 0L
    t0 <- proc.time()[["elapsed"]]
    r <- suppressMessages(cases[[cn]]())
    s <- proc.time()[["elapsed"]] - t0
    cat(sprintf("ROW %s %d %s %.2f %d %d %s\n", cn, rd, a, s, cnt$k,
                cnt$won, llp_of(r)))
    X[[length(X) + 1L]] <- data.frame(case = cn, arm = a, s = s)
    flush.console()
  }
}
set_arm("on")
X <- do.call(rbind, X)
for (cn in names(cases)) {
  m <- tapply(X$s[X$case == cn], X$arm[X$case == cn], min)
  cat(sprintf("%s: min s off %.2f on %.2f off2 %.2f | on/off %.3f control off2/off %.3f\n",
              cn, m[["off"]], m[["on"]], m[["off2"]], m[["on"]] / m[["off"]],
              m[["off2"]] / m[["off"]]))
}
