# Lane gpstop: the designs of dev/review-arstep-cost.R (the review of lane
# arstep), in ONE process on the lane library, the lane against the mutant
# "norelin" (gp_body_relin_applies() returns FALSE: the 0.81.0 route),
# interleaved. None of these fits holds a nonlinear body, so the two arms
# run the same code but for that one test.
#   noar       gp_marginal_ar() refuses every AR term: the 0.79.1 route for
#              ar() (the working start), second descent still on
#   nosecond   gp_second_descent() returns its input: the 0.79.1 route off
#              the gaussian
# Designs: (a) n distinct positions with ar(time); (b) panels: S subjects by
# T times, x = time, ar(time, gr = subject), so rows outnumber positions;
# (c) off the gaussian, poisson and bernoulli y ~ gp(x) at n 120 and 240.
# Usage: Rscript review-arstep-cost.R <design ids, comma> <reps>
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
a <- commandArgs(TRUE)
ids <- strsplit(a[1], ",")[[1]]
reps <- as.integer(a[2])
rv_libs("lane")
ns <- asNamespace("frmtmb")
orig_ar <- ns$gp_marginal_ar
orig_sd <- ns$gp_second_descent
orig_ra <- ns$gp_body_relin_applies
set_arm <- function(arm) {
  rv_put("gp_marginal_ar", orig_ar)
  rv_put("gp_second_descent", orig_sd)
  rv_put("gp_body_relin_applies", orig_ra)
  if (arm == "norelin") rv_put("gp_body_relin_applies", function(...) FALSE)
  if (arm == "noar") rv_put("gp_marginal_ar", function(frame) {
    if (length(frame[["autocor"]] %||% list())) FALSE else NULL
  })
  if (arm == "nosecond") rv_put("gp_second_descent", function(obj, opt, ...) {
    if (exists(".frmtmb_gp_alt", envir = obj$env, inherits = FALSE)) {
      rm(".frmtmb_gp_alt", envir = obj$env)
    }
    opt
  })
}
designs <- list(
  a500 = function() {
    set.seed(9500); n <- 500
    x <- sort(runif(n, 0, 10))
    list(d = data.frame(x = x, t = seq_len(n),
                        y = 2 * (x > 5) - 1 + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.1))),
         f = y ~ gp(x) + ar(time = t), fam = gaussian(), other = "norelin")
  },
  a300 = function() {
    set.seed(9300); n <- 300
    x <- sort(runif(n, 0, 10))
    list(d = data.frame(x = x, t = seq_len(n),
                        y = sin(x) + as.numeric(arima.sim(list(ar = 0.6), n, sd = 0.1))),
         f = y ~ gp(x) + ar(time = t), fam = gaussian(), other = "norelin")
  },
  p600 = function() {
    set.seed(9600); S <- 12; T <- 50
    d <- expand.grid(t = seq_len(T), s = factor(seq_len(S)))
    d$x <- d$t / 5
    d$y <- sin(d$x) + rnorm(S, 0, 0.3)[d$s] +
      unlist(lapply(seq_len(S), function(i) as.numeric(arima.sim(list(ar = 0.5), T, sd = 0.2))))
    list(d = d, f = y ~ gp(x) + ar(time = t, gr = s), fam = gaussian(), other = "norelin")
  },
  p1000 = function() {
    set.seed(9601); S <- 25; T <- 40
    d <- expand.grid(t = seq_len(T), s = factor(seq_len(S)))
    d$x <- d$t / 4
    d$y <- sin(d$x) + unlist(lapply(seq_len(S), function(i)
      as.numeric(arima.sim(list(ar = 0.5), T, sd = 0.2))))
    list(d = d, f = y ~ gp(x) + ar(time = t, gr = s), fam = gaussian(), other = "norelin")
  },
  pois120 = function() {
    set.seed(9120); n <- 120; x <- sort(runif(n, 0, 10))
    list(d = data.frame(x = x, y = rpois(n, exp(1 + sin(3 * x) / 2))), f = y ~ gp(x),
         fam = poisson(), other = "norelin")
  },
  pois240 = function() {
    set.seed(9240); n <- 240; x <- sort(runif(n, 0, 10))
    list(d = data.frame(x = x, y = rpois(n, exp(1 + (2 * (x > 5) - 1) / 2))), f = y ~ gp(x),
         fam = poisson(), other = "norelin")
  },
  bern240 = function() {
    set.seed(9241); n <- 240; x <- sort(runif(n, 0, 10))
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(2 * sin(3 * x) + 0.5 * x - 2.5))),
         f = y ~ gp(x), fam = bernoulli(), other = "norelin")
  })
for (id in ids) {
  m <- designs[[id]]()
  for (r in seq_len(reps)) for (arm in if (r %% 2) c("lane", m$other) else c(m$other, "lane")) {
    set_arm(arm)
    invisible(gc(reset = TRUE))
    t0 <- proc.time()[[3]]
    q <- rv_quiet(frm(m$f, data = m$d, family = m$fam))
    el <- proc.time()[[3]] - t0
    g <- gc()
    fit <- q$value
    if (inherits(fit, "error")) {
      cat(sprintf("COST %s rep %d %-8s ERROR %s\n", id, r, arm, conditionMessage(fit)))
      next
    }
    cat(sprintf("COST %s n %d rep %d %-8s wall %.2f ll %.6f conv %d w %d marg %s second %s Rmax_MB %.0f\n",
                id, nrow(m$d), r, arm, el, as.numeric(logLik(fit)), fit$opt$convergence,
                length(q$warnings), !is.null(fit$obj$env$.frmtmb_gp_marginal),
                if (is.null(fit$opt$gp_second)) "-" else paste(signif(fit$opt$gp_second, 3), collapse = ","),
                sum(g[, ncol(g)])))
    flush.console()
  }
}
set_arm("lane")
cat("DONE\n")
