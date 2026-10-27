# Reviewer, fix round 2: the b/Intercept branch of resolve_priorlist()
# changed. frm() MAP fits with b, Intercept, sigma and dpar priors on
# models WITHOUT a cov = TRUE term must be unchanged against rellib-r13.
# Prints the objective at the optimum and every estimate to 12 digits;
# the two logs are diffed minus their first line.
# Usage: Rscript dev/review-acprior-r2-objs.R before|after
# Log: dev/review-acprior-log/r2-objs-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(14)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = rep(1:12, each = 10),
                t = rep(1:10, 12))
d$y <- 1 + 0.5 * d$x + rnorm(12)[d$g] * 0.5 + rnorm(n, 0, exp(0.2 * d$z))
d$y2 <- 2 - 0.4 * d$x + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
d$yb <- rbinom(n, 1, plogis(0.2 + d$x))
cases <- list(
  gauss = list(bf(y ~ x + z), gaussian(),
               set_prior("normal(0, 1)", class = "b") +
                 set_prior("student_t(3, 1, 2)", class = "Intercept") +
                 set_prior("normal(1, 0.5)", class = "sigma")),
  re = list(bf(y ~ x + (1 + x | g)), gaussian(),
            set_prior("normal(0, 2)", class = "b", coef = "x") +
              set_prior("normal(0, 1)", class = "sd") +
              set_prior("lkj(2)", class = "cor")),
  dist = list(bf(y ~ x, sigma ~ z), gaussian(),
              set_prior("normal(0, 1)", class = "b", dpar = "sigma") +
                set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma")),
  pois = list(bf(cnt ~ x), poisson(),
              set_prior("normal(0, 0.5)", class = "b") +
                set_prior("normal(0, 2)", class = "Intercept")),
  logit = list(bf(yb ~ x + z), bernoulli(),
               set_prior("student_t(3, 0, 1)", class = "b")),
  mv = list(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), NULL,
            set_prior("normal(0, 1)", class = "b", resp = "y2") +
              set_prior("normal(1, 1)", class = "sigma", resp = "y") +
              set_prior("lkj(2)", class = "rescor")),
  arfalse = list(bf(y ~ x + ar(t, g)), gaussian(),
                 set_prior("normal(0, 1)", class = "b") +
                   set_prior("normal(1, 1)", class = "sigma", lb = 0.2)),
  cosy = list(bf(y ~ x + cosy(t, g)), gaussian(),
              set_prior("normal(1, 1)", class = "sigma") +
                set_prior("normal(0, 1)", class = "Intercept"))
)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  r <- tryCatch({
    fit <- if (is.null(cs[[2]])) {
      suppressWarnings(frm(cs[[1]], data = d, prior = cs[[3]]))
    } else suppressWarnings(frm(cs[[1]], data = d, family = cs[[2]],
                                prior = cs[[3]]))
    est <- unlist(fit$estimates)
    sprintf("obj %.12f | %s", fit$obj$fn(fit$obj$env$last.par.best),
            paste(sprintf("%.12f", est), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-8s %s\n", nm, r))
}
