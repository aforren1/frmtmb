# brms's Stan log density against frmtmb's joint density at frmtmb's
# estimates, for the gp(by = ) Hilbert-space forms (check C of
# dev/brms-likelihood-tests.md). Flat priors on both sides. Prints the
# measured constant and the largest inner-parameter gradient.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-gpby/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("lib:", find.package("frmtmb"), "\n")
env <- testthat::test_env("frmtmb")
for (h in list.files("tests/testthat", "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
src <- parse("tests/testthat/test-gp-by.R")
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
}
options(frmtmb.brms_lp_report = TRUE)
d <- env$gpby_data()
d$z <- { set.seed(77); stats::runif(nrow(d), 0, 3) }
# z must carry signal, or its own length scale runs off to where the
# spectral density underflows and brms's zgp is not defined
d$y <- d$y + 0.8 * sin(2 * d$z)
forms <- c("y ~ gp(x, by = f, k = 8)",
           "y ~ gp(x, by = w, k = 8)",
           "y ~ gp(x, by = f, k = 8, cmc = FALSE)",
           "y ~ gp(x, by = f, k = 8, gr = FALSE)",
           "y ~ gp(x, z, by = f, k = 5, iso = FALSE)",
           "y ~ gp(x, z, by = f, k = 5)")
for (form in forms) {
  cat("==", form, "\n")
  fit <- frm(bf(stats::as.formula(form)), data = d)
  r <- try(local({
    environment(env$brms_lp_check) <- env
    env$brms_lp_check(brms::bf(stats::as.formula(form)), gaussian(), d,
                      fit, joint = TRUE)
  }))
  if (!inherits(r, "try-error")) {
    cat(sprintf("RES %s | const %.3e | max_grad %.3e | ours %.6f\n", form,
                r$measured_const, r$max_grad, r$ours))
  }
}
