# Two more naming cases against brms 2.23.0, after a review of the plan:
# smooths in two predictors (sds_ order and prefix) and residual
# correlation priors on both responses of a multivariate model.
#   Rscript dev/priordraws-names2.R
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
loadNamespace("brms")
cat("lib:", find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60),
                h = factor(rep(1:5, each = 12)), t = rep(1:12, 5))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$w <- 0.3 * d$z + rnorm(60)
q <- function(e) suppressWarnings(suppressMessages(e))
cases <- list(
  smooths = list(
    brms = list(brms::bf(y ~ s(x) + s(z), sigma ~ s(x))),
    frm = list(frmtmb::bf(y ~ s(x) + s(z), sigma ~ s(x)), NULL)),
  mv_ar = list(
    brms = list(brms::bf(y ~ x + ar(time = t, gr = h)) +
                  brms::bf(w ~ z + ar(time = t, gr = h)) +
                  brms::set_rescor(FALSE),
                prior = c(brms::prior(normal(0, 0.5), class = ar,
                                      resp = y),
                          brms::prior(normal(0, 0.4), class = ar,
                                      resp = w))),
    frm = list(frmtmb::bf(y ~ x + ar(time = t, gr = h)) +
                 frmtmb::bf(w ~ z + ar(time = t, gr = h)) +
                 frmtmb::set_rescor(FALSE),
               frmtmb::set_prior("normal(0, 0.5)", class = "ar",
                                 resp = "y") +
                 frmtmb::set_prior("normal(0, 0.4)", class = "ar",
                                   resp = "w"))))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fb <- q(do.call(brms::brm, c(cs$brms, list(data = d, sample_prior = "yes",
                                             chains = 1, iter = 200,
                                             seed = 1, refresh = 0))))
  ds <- q(frm_sample(cs$frm[[1]], data = d, family = gaussian(),
                     prior = cs$frm[[2]], sample_prior = "yes", chains = 1,
                     iter = 200, seed = 1, refresh = 0))
  bn <- names(brms::prior_draws(fb))
  fn <- names(prior_draws(ds))
  cat(sprintf("%-8s %s\n  frmtmb: %s\n  brms:   %s\n", nm,
              if (identical(fn, bn)) "SAME" else "DIFF",
              paste(fn, collapse = " "), paste(bn, collapse = " ")))
}
