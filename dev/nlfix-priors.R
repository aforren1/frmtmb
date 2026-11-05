# The prior table and a threshold prior on an nl ordinal fit, both routes,
# against the linear model and brms's default_prior(). Seed 49.
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
d <- E$eq_dat("ordinal", seed = 49)
cat("== frm default_prior, linear\n")
print(default_prior(bf(y ~ x), data = d, family = cumulative()))
cat("== frm default_prior, nl\n")
print(default_prior(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                    family = cumulative()))
cat("== frm default_prior, nl, disc\n")
print(default_prior(bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE), data = d,
                    family = cumulative()))
cat("== sample route, nl\n")
print(default_prior(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                    family = cumulative(), route = "sample"))
cat("== brms\n")
print(brms::default_prior(brms::bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), data = d,
                          family = brms::cumulative())[, c("prior", "class",
                                                           "coef", "dpar",
                                                           "nlpar")])
pr <- set_prior("normal(0, 1)", class = "Intercept")
fl <- suppressWarnings(frm(bf(y ~ x), family = cumulative(), data = d,
                           prior = pr))
cat("linear with Intercept prior:", logLik(fl), "\n")
fn <- tryCatch(suppressWarnings(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
                                    family = cumulative(), data = d,
                                    prior = pr)),
               error = function(e) e)
if (inherits(fn, "error")) {
  cat("nl with Intercept prior ERROR:", conditionMessage(fn), "\n")
} else {
  cat("nl with Intercept prior:", logLik(fn), "\n")
  print(prior_summary(fn))
}
fn2 <- tryCatch(suppressWarnings(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE),
                                     family = cumulative(), data = d,
                                     prior = set_prior("normal(0, 1)", class = "b", nlpar = "b1"))),
                error = function(e) e)
if (inherits(fn2, "error")) cat("ERR", conditionMessage(fn2), "\n") else
  cat("nl with b1 prior:", logLik(fn2), "\n")
