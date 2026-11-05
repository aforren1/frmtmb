# nl + ar(cov = FALSE): fitted() against brms's posterior_epred at the same
# values (spec ac_ar of dev/nlfix-l1-specs.R, seed 1582). Which reads the
# ARMA shift: the body alone, or the body plus the shift? And on the same
# model written with a linear predictor (log-linear body is not linear, so a
# plain y ~ x with nl spelled b0 + b1 * x, seed 1582).
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
d <- E$eq_dat(seed = 1582)
d$y <- 2 + 0.8 * d$x + rnorm(nrow(d), 0, 0.5)
fl <- frm(bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d)
fn <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
                              autocor = ~ ar(time = t, gr = g, cov = FALSE)),
                           data = d,
                           start = list(beta = c("b0_(Intercept)" = 2,
                                                 "b1_(Intercept)" = 0.5))))
cat("logLik linear, nl:", logLik(fl), logLik(fn), "\n")
b <- fn$estimates$beta
body <- b[[1]] + b[[2]] * d$x
fv <- fitted(fn)[, "Estimate"]
fv2 <- fitted(fl)[, "Estimate"]
cat("max|frm nl fitted - frm linear fitted|", max(abs(fv - fv2)), "\n")
cat("max|frm nl fitted - body|", max(abs(fv - body)), "\n")
# brms at the same values
fx <- E$brms_fixed_fit
print(args(fx))
bl <- brms::bf(y ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE))
bn <- brms::bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE,
               autocor = ~ ar(time = t, gr = g, cov = FALSE))
# brms's own R code path, read rather than recalled
cat(capture.output(print(body(brms:::predictor.bprepnl))), sep = "\n")
cat(capture.output(print(body(brms:::predictor.bprepl))), sep = "\n")
