# Review of lane nlstart: cost of the probe on fresh designs, probe off /
# on / off again in rotating order within one process, 3 rounds, min per
# arm. Also an nlf() body with frm_ode() (frmtmb.ode) from its default start.
# Usage: REV_ARM=after Rscript dev/review-nlstart-cost.R
source("dev/review-nlstart-common.R")
time1 <- function(expr_fn, probe) {
  t0 <- proc.time()[["elapsed"]]
  r <- withr::with_options(list(frmtmb.nl_sign_probe = probe), tt(expr_fn()))
  list(t = proc.time()[["elapsed"]] - t0, ll = llv(r),
       multi = any(attr(r, "mc") == "frmtmb_nl_multistart"))
}
bench <- function(label, expr_fn, rounds = 3) {
  arms <- c("off", "on", "off2")
  tm <- list(off = c(), on = c(), off2 = c()); ll <- list()
  for (r in seq_len(rounds)) {
    ord <- arms[(seq_along(arms) + r - 2) %% 3 + 1]
    for (a in ord) {
      x <- time1(expr_fn, a == "on")
      tm[[a]] <- c(tm[[a]], x$t); ll[[a]] <- x$ll
      if (a == "on") multi <- x$multi
    }
  }
  m <- vapply(tm, min, 0)
  cat(sprintf("%-26s min s: off %.2f on %.2f off2 %.2f | on/off %.2f control %.2f | ll off %.6f on %.6f multi %s\n",
              label, m[["off"]], m[["on"]], m[["off2"]],
              m[["on"]] / mean(c(m[["off"]], m[["off2"]])),
              m[["off2"]] / m[["off"]], ll$off, ll$on, multi))
  flush.console()
}
set.seed(8101)
n <- as.integer(Sys.getenv("REV_N", "20000"))
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                g = factor(sample(1:30, n, TRUE)))
d$ug <- rnorm(30, 0, 0.4)[d$g]
d$y_mm <- 3 * d$xp / (1 + d$xp) + rnorm(n, 0, 0.3)
d$y_l3 <- rpois(n, 3 / (1 + exp((0.3 - d$x) / 0.7)))
d$y_s4 <- 0.5 + 2.5 / (1 + exp(-d$x / 0.6)) + rnorm(n, 0, 0.3)
d$y_5 <- exp(0.5 * d$x) - (1 - exp(0.5 * d$z)) + 0.3 * d$w + rnorm(n, 0, 0.3)
d$y_re <- (1.5 + d$ug) * exp(0.4 * d$x) + rnorm(n, 0, 0.3)
d$y_be <- rbinom(n, 1, plogis(1.5 * exp(-0.6 * d$x) - 2))
d$y_lin <- 1 + 0.5 * d$x + rnorm(n)
bench("mm gaussian", function() frm(bf(y_mm ~ Vm * xp / (K + xp),
  Vm + K ~ 1, nl = TRUE), data = d))
bench("logis3 poisson", function() frm(bf(y_l3 ~ Asym / (1 +
  exp(-k * (x - xmid))), Asym + xmid + k ~ 1, nl = TRUE),
  family = poisson("identity"), data = d))
bench("sig4 gaussian", function() frm(bf(y_s4 ~ lo + (hi - lo) /
  (1 + exp(-k * (x - xmid))), lo + hi + xmid + k ~ 1, nl = TRUE),
  data = d))
bench("five gaussian", function() frm(bf(y_5 ~ b1 * exp(b2 * x) +
  b3 * (1 - exp(-b4 * z)) + b5 * w, b1 + b2 + b3 + b4 + b5 ~ 1, nl = TRUE),
  data = d))
bench("exp RE gaussian", function() frm(bf(y_re ~ b1 * exp(b2 * x),
  b1 ~ 1 + (1 | g), b2 ~ 1, nl = TRUE), data = d))
bench("exp bernoulli", function() frm(bf(y_be ~ b1 * exp(b2 * x),
  b1 + b2 ~ 1, nl = TRUE), family = bernoulli(), data = d))
bench("a + b*x (not probed)", function() frm(bf(y_lin ~ a + b * x,
  a + b ~ 1, nl = TRUE), data = d))

if (requireNamespace("frmtmb.ode", quietly = TRUE) &&
      requireNamespace("RTMBode", quietly = TRUE)) {
  suppressPackageStartupMessages(library(frmtmb.ode))
  cat("frmtmb.ode from", dirname(find.package("frmtmb.ode")), "\n")
  pk_dyn <- function(t, y, p) {
    list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
  }
  set.seed(2026)
  tq <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  id <- factor(rep(seq_len(8), each = length(tq)))
  time <- rep(tq, 8)
  ka <- exp(stats::rnorm(8, 0, 0.3))[as.integer(id)]
  ke <- exp(log(0.2) + stats::rnorm(8, 0, 0.25))[as.integer(id)]
  mu <- 100 * ka / (10 * (ka - ke)) * (exp(-ke * time) - exp(-ka * time))
  dd <- data.frame(id = id, time = time, dose = 100,
                   conc = mu + stats::rnorm(length(mu), 0, 0.3))
  bench("ode nlf, default start", function() frm(bf(conc ~ pk / 10, nl = TRUE) +
    nlf(pk ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                     parms = list(exp(lka), exp(lke)),
                     group = id, output = 2L)) +
    lf(lka ~ 1, lke ~ 1) + gaussian(), data = dd), rounds = 2)
} else cat("frmtmb.ode or RTMBode not available\n")
cat("DONE\n")
