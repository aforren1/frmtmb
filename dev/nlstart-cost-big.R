# Wall time of the sign probe on larger fits, in one process: each design
# fitted with the probe off, on and off again, in rotating order, over 3
# rounds; the minimum per arm is reported. off2 / off is the control.
# Seed 7301.
source("dev/nlstart-common.R")
set.seed(7301)
n <- 20000
g <- factor(sample(1:200, n, TRUE))
x <- round(rnorm(n), 3)
u <- rnorm(200, 0, 0.3)
d1 <- data.frame(x = x, y = 1.5 * exp(0.6 * x) + rnorm(n, 0, 0.5))
d2 <- data.frame(x = x, g = g,
                 y = rpois(n, exp(log(1.5) + u[as.integer(g)]) *
                             exp(0.4 * x)))
d3 <- data.frame(x = x, y = 2 * (1 - exp(-0.8 * x)) + rnorm(n, 0, 0.5))
designs <- list(
  gauss_exp_2e4 = list(f = bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
                       fam = gaussian(), d = d1),
  gauss_sat_2e4 = list(f = bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1,
                              nl = TRUE), fam = gaussian(), d = d3),
  pois_re_2e4 = list(f = bf(y ~ b1 * exp(b2 * x), b1 ~ 1 + (1 | g), b2 ~ 1,
                            nl = TRUE), fam = poisson("identity"), d = d2))
fit_t <- function(ds, probe) {
  old <- options(frmtmb.nl_sign_probe = probe)
  on.exit(options(old))
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(ds$f, family = ds$fam, data = ds$d))
  c(t = proc.time()[["elapsed"]] - t0,
    ll = if (iserr(r)) NA else as.numeric(logLik(r)))
}
for (nm in names(designs)) {
  ds <- designs[[nm]]
  res <- list(off = NULL, on = NULL, off2 = NULL)
  arms <- c(off = FALSE, on = TRUE, off2 = FALSE)
  for (round in 0:2) {
    ord <- names(arms)[(seq_along(arms) + round) %% 3L + 1L]
    for (a in ord) res[[a]] <- rbind(res[[a]], fit_t(ds, arms[[a]]))
  }
  m <- vapply(res, function(z) min(z[, "t"]), 0)
  cat(sprintf(paste0("BIG %s min s: off %.2f on %.2f off2 %.2f | on/off %.3f",
                     " control %.3f | ll off %.6f on %.6f\n"), nm, m[["off"]],
              m[["on"]], m[["off2"]], m[["on"]] / m[["off"]],
              m[["off2"]] / m[["off"]], res$off[1, "ll"], res$on[1, "ll"]))
  flush.console()
}
