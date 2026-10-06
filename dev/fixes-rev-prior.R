# Reviewer of lane fixes, claim 4b: per-threshold Intercept coef rows.
# default_prior() rows against brms 2.23.0's, and the penalized
# objective against a hand computation WITH a predictor, so the
# centering offset enters, for both orders of a class row and a coef row.
#   Rscript dev/fixes-rev-prior.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(9090)
n <- 250
d <- data.frame(x = rnorm(n, 1, 1), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -0.2) + (u > 0.8) + (u > 1.8) + (u > 2.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
d$y2 <- 1L + (u + rnorm(n) > 1)
rows_f <- function(dp) {
  r <- dp[dp$class %in% c("Intercept", "delta"), ]
  paste(r$class, r$group, r$coef, if (!is.null(r$resp)) r$resp, sep = "|")
}
rows_b <- function(dp) {
  r <- dp[dp$class %in% c("Intercept", "delta"), ]
  paste(r$class, r$group, r$coef, if (nzchar(paste(r$resp, collapse = "")))
    r$resp, sep = "|")
}
cases <- list(
  cum = list(y ~ x, "cumulative", list()),
  hurdle = list(yh ~ x, "hurdle_cumulative", list()),
  acat_cs = list(y ~ x + cs(z), "acat", list()),
  cr_gr = list(yg | thres(gr = g) ~ x, "cratio", list()),
  cum_eq = list(y ~ x, "cumulative", list(threshold = "equidistant")),
  cum_eq_gr = list(yg | thres(gr = g) ~ x, "cumulative",
                   list(threshold = "equidistant")),
  sr_stz = list(y ~ x, "sratio", list(threshold = "sum_to_zero")),
  binary_cum = list(y2 ~ x, "cumulative", list()))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  ff <- do.call(get(cs[[2]], asNamespace("frmtmb")), cs[[3]])
  fb <- do.call(get(cs[[2]], asNamespace("brms")), cs[[3]])
  a <- tryCatch(rows_f(default_prior(cs[[1]], data = d, family = ff)),
                error = function(e) paste("ERROR", conditionMessage(e)))
  b <- tryCatch(rows_b(brms::default_prior(cs[[1]], data = d, family = fb)),
                error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-11s identical %s\n  frm : %s\n  brms: %s\n", nm,
              identical(a, b), paste(a, collapse = "  "),
              paste(b, collapse = "  ")))
}
# multivariate: two ordinal responses
d$w <- 1L + (u + rnorm(n) > 0.5) + (u + rnorm(n) > 1.5)
a <- tryCatch(rows_f(default_prior(bf(y ~ x) + bf(w ~ x), data = d,
                                   family = cumulative())),
              error = function(e) paste("ERROR", conditionMessage(e)))
b <- tryCatch(rows_b(brms::default_prior(brms::bf(y ~ x) + brms::bf(w ~ x),
                                         data = d,
                                         family = brms::cumulative())),
              error = function(e) paste("ERROR", conditionMessage(e)))
cat(sprintf("%-11s identical %s\n  frm : %s\n  brms: %s\n", "mv_cum",
            identical(a, b), paste(a, collapse = "  "),
            paste(b, collapse = "  ")))

# the penalized objective: fit with a prior, evaluated against the
# prior-free objective at the same point, against a hand log prior on
# the CENTERED thresholds tau - mean(x) * b_x with the ordered map's
# log-Jacobian sum(log increments)
hand <- function(fit, dens) {
  raw <- fit$estimates$tau_raw
  tau <- frmtmb:::ord_threshold_values(family(fit), raw)
  bx <- fit$estimates$beta[1]
  tc <- tau - mean(d$x) * bx
  lj <- if (family(fit)$family %in% c("cumulative")) sum(raw[-1]) else 0
  sum(vapply(seq_along(tc), function(k) dens(k, tc[k]), 0)) + lj
}
pen <- function(fit, f0) {
  p <- fit$opt$par
  -(fit$obj$fn(p) - f0$obj$fn(p))
}
for (fam in c("cumulative", "sratio")) {
  famf <- get(fam, asNamespace("frmtmb"))()
  f0 <- frm(y ~ x, family = famf, data = d)
  pc <- set_prior("normal(0, 3)", class = "Intercept")
  p2 <- set_prior("normal(0.7, 0.4)", class = "Intercept", coef = "2")
  p4 <- set_prior("student_t(3, 2, 1)", class = "Intercept", coef = "4")
  dn <- function(k, v) {
    if (k == 2) stats::dnorm(v, 0.7, 0.4, log = TRUE) else
      if (k == 4) stats::dt(v - 2, 3, log = TRUE) else
        stats::dnorm(v, 0, 3, log = TRUE)
  }
  dn_only2 <- function(k, v) if (k == 2) stats::dnorm(v, 0.7, 0.4, log = TRUE)
    else 0
  for (lab in c("class+coef2+coef4", "coef2+class+coef4", "coef4+coef2+class",
                "coef2 only")) {
    pl <- switch(lab, "class+coef2+coef4" = pc + p2 + p4,
                 "coef2+class+coef4" = p2 + pc + p4,
                 "coef4+coef2+class" = p4 + p2 + pc,
                 "coef2 only" = p2)
    fit <- tryCatch(frm(y ~ x, family = famf, data = d, prior = pl),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      cat(fam, lab, "ERROR", conditionMessage(fit), "\n"); next
    }
    h <- hand(fit, if (lab == "coef2 only") dn_only2 else dn)
    cat(sprintf("%-10s %-18s penalty %.10f hand %.10f rel %.3g\n", fam, lab,
                pen(fit, f0), h, abs(pen(fit, f0) - h) / abs(h)))
  }
}
# a coef row for a threshold of one level of grouped thresholds, with a
# predictor: brms does not center grouped thresholds, says the resolver
fg0 <- frm(yg | thres(gr = g) ~ x, family = cumulative(), data = d)
fg <- frm(yg | thres(gr = g) ~ x, family = cumulative(), data = d,
          prior = set_prior("normal(0.5, 0.4)", class = "Intercept",
                            group = "b", coef = "2"))
raw <- fg$estimates$tau_raw
tau <- frmtmb:::ord_threshold_values(family(fg), raw)
h <- stats::dnorm(tau[6], 0.5, 0.4, log = TRUE) + sum(raw[6:7])
cat(sprintf("grouped b coef2: penalty %.10f hand(uncentered) %.10f\n",
            pen(fg, fg0), h))
# prior_summary() and print of a fit carrying a coef row
fit <- frm(y ~ x, family = cumulative(), data = d,
           prior = set_prior("normal(0, 3)", class = "Intercept") +
             set_prior("normal(0.7, 0.4)", class = "Intercept", coef = "2"))
print(tryCatch(prior_summary(fit), error = function(e) conditionMessage(e)))
# simulate from the prior, and get_prior(route = "sample")
cat("simulate(prior draws):", tryCatch({
  s <- simulate(fit, nsim = 1, newparams = "prior"); "ok"
}, error = function(e) conditionMessage(e)), "\n")
