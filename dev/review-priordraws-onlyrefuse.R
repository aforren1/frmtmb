# Reviewer: which models sample_prior = "only" refuses, brms 2.23.0
# (stancode, no compile) against frm_sample() on the lane build.
#   Rscript dev/review-priordraws-onlyrefuse.R
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$yc <- rpois(60, 3)
d$yo <- pmin(d$yc, 4) + 1
d$t <- rep(1:12, 5)
d$yb <- rbinom(60, 1, 0.5)
Pb <- function(...) brms::set_prior(...)
Pf <- function(...) frmtmb::set_prior(...)
nb <- "normal(0, 1)"
cases <- list(
  smooth_bprior = list(y ~ s(x), "gaussian", list(list(nb, "b"))),
  smooth_default = list(y ~ s(x), "gaussian", list()),
  sigma_x_default = list(list(y ~ x, sigma ~ z), "gaussian",
                         list(list(nb, "b"))),
  sigma_x_proper = list(list(y ~ x, sigma ~ z), "gaussian",
                        list(list(nb, "b"), list(nb, "b", "sigma"))),
  cumulative = list(yo ~ x, "cumulative", list(list(nb, "b"))),
  ar_default = list(y ~ x + ar(time = t, gr = h), "gaussian",
                    list(list(nb, "b"))),
  ar_cov_default = list(y ~ x + ar(time = t, gr = h, cov = TRUE), "gaussian",
                        list(list(nb, "b"))),
  cosy_default = list(y ~ x + cosy(time = t, gr = h), "gaussian",
                      list(list(nb, "b"))),
  cor_re = list(y ~ x + (1 + x | g), "gaussian", list(list(nb, "b"))),
  zip = list(yc ~ x, "zero_inflated_poisson", list(list(nb, "b"))),
  mixture = list(y ~ x, "mixture2", list(list(nb, "b", "mu1"),
                                         list(nb, "b", "mu2"),
                                         list(nb, "Intercept", "mu1"),
                                         list("normal(1, 1)", "Intercept",
                                              "mu2"))),
  bernoulli = list(yb ~ x + (1 | g), "bernoulli", list(list(nb, "b"))),
  student_re = list(y ~ x + (1 | gr(g, dist = "student")), "gaussian",
                    list(list(nb, "b"))),
  mo = list(y ~ mo(m), "gaussian", list(list(nb, "b"))),
  noint = list(y ~ 0 + x, "gaussian", list(list(nb, "b"))),
  int_as_b = list(list(y ~ 0 + Intercept + x), "gaussian", list(list(nb, "b"))),
  offset_only = list(y ~ 1 + offset(x), "gaussian", list())
)
fam_of <- function(s, pkg) {
  if (s == "mixture2") {
    return(if (pkg == "brms") brms::mixture(brms::gaussian(), brms::gaussian())
           else frmtmb::mixture(gaussian(), gaussian()))
  }
  if (pkg == "brms") get(s, asNamespace("brms"))() else {
    f <- tryCatch(get(s, asNamespace("frmtmb")), error = function(e) get(s))
    f()
  }
}
mk <- function(pl, P) {
  out <- NULL
  for (p in pl) {
    one <- P(p[[1]], class = p[[2]],
             dpar = if (length(p) > 2) p[[3]] else "")
    out <- if (is.null(out)) one else out + one
  }
  out
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n=====", nm, "=====\n")
  rb <- tryCatch({
    brms::stancode(do.call(brms::bf, if (is.list(cs[[1]])) cs[[1]] else list(cs[[1]])), data = d, family = fam_of(cs[[2]], "brms"),
                   prior = mk(cs[[3]], Pb), sample_prior = "only")
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("brms:", rb, "\n")
  rf <- tryCatch({
    ds <- suppressWarnings(suppressMessages(frm_sample(
      do.call(frmtmb::bf, if (is.list(cs[[1]])) cs[[1]] else list(cs[[1]])),
      data = d, family = fam_of(cs[[2]], "frm"), prior = mk(cs[[3]], Pf),
      sample_prior = "only", chains = 1, iter = 60, seed = 1, refresh = 0)))
    paste("OK", ndraws(ds))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("frm: ", rf, "\n")
}
