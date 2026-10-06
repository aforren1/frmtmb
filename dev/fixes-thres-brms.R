# Lane fixes, item 4a: brms 2.23.0's posterior_linpred(incl_thres = TRUE)
# on the ordinal families at frmtmb's estimates (Fixed_param draws):
# its array shape, its dimnames, and its values against
# disc * (thres - mu) (cumulative, sratio) or disc * (mu - thres)
# (cratio, acat) from frmtmb's estimates.
#   Rscript dev/fixes-thres-brms.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
# the helpers read frmtmb internals, as they do under test_file()
henv <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = henv)
}
brms_fixed_fit <- henv$brms_fixed_fit
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
# group b never reaches the top category, so it has one threshold fewer
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
shapes <- list(
  cum = list(bf(y ~ x), cumulative(), brms::bf(y ~ x), brms::cumulative()),
  sr_cs = list(bf(y ~ x + cs(z)), sratio(), brms::bf(y ~ x + cs(z)),
               brms::sratio()),
  cr_disc = list(bf(y ~ x, disc ~ 0 + z), cratio(),
                 brms::bf(y ~ x, disc ~ 0 + z), brms::cratio()),
  acat_eq = list(bf(y ~ x), acat(threshold = "equidistant"),
                 brms::bf(y ~ x), brms::acat(threshold = "equidistant")),
  sr_stz = list(bf(y ~ x), sratio(threshold = "sum_to_zero"),
                 brms::bf(y ~ x), brms::sratio(threshold = "sum_to_zero")),
  cum_stz_probit = list(bf(y ~ x), cumulative("probit", threshold = "sum_to_zero"),
                 brms::bf(y ~ x), brms::cumulative("probit", threshold = "sum_to_zero")),
  cum_gr = list(bf(yg | thres(gr = g) ~ x), cumulative(),
                brms::bf(yg | thres(gr = g) ~ x), brms::cumulative()),
  hurdle = list(bf(yh ~ x), hurdle_cumulative(), brms::bf(yh ~ x),
                brms::hurdle_cumulative()))
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("==", nm, "\n")
  fit <- tryCatch(suppressWarnings(frm(s[[1]], family = s[[2]], data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat("  frm ERROR:", conditionMessage(fit), "\n"); next
  }
  bb <- tryCatch(brms_fixed_fit(s[[3]], s[[4]], d, fit, ndraws = 4),
                 error = function(e) e)
  if (inherits(bb, "error")) {
    cat("  brms fixed fit ERROR:", conditionMessage(bb), "\n"); next
  }
  pl <- tryCatch(brms::posterior_linpred(bb, incl_thres = TRUE),
                 error = function(e) e)
  if (inherits(pl, "error")) {
    cat("  brms ERROR:", conditionMessage(pl), "\n"); next
  }
  cat("  dim:", dim(pl), "\n  dimnames:\n")
  str(dimnames(pl))
  cat("  draw 1, rows 1-2:\n"); print(pl[1, 1:2, ])
  # the brms value of the plain linear predictor and its disc, for the
  # by-hand reconstruction
  mu <- brms::posterior_linpred(bb)
  cat("  draws identical:", isTRUE(all.equal(pl[1, , ], pl[4, , ])), "\n")
  e <- variables(bb)
  cat("  brms variables (b_Intercept*, delta, disc*):",
      grep("^b_Intercept|^delta|^disc|^b_disc|^Intercept", e, value = TRUE),
      "\n")
  saveRDS(list(pl = pl, mu = mu), file.path("dev/fixes-log",
                                            paste0("thres-", nm, ".rds")))
}
