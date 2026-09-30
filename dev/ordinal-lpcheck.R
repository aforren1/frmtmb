# Log-density identity against brms 2.23.0's compiled Stan program for
# the ordinal families with disc, each threshold structure, and acat
# off the logit. Runs the
# test suite's own harness (brms_lp_check(), helper-brms.R) outside
# testthat so every row prints its measured constant and gradient.
# Seed 20260930 for the data. Output: dev/ordinal-log-lpcheck.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordinal/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("frmtmb from", find.package("frmtmb"), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R")) {
  sys.source(file.path("tests/testthat", h), envir = env)
}
options(frmtmb.brms_lp_report = TRUE)
lp_tight <- frmtmb_control(grad_tol = 1e-6, restarts = 3)

ord_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  d$yh <- ifelse(runif(n) < plogis(-1 + 0.3 * d$x), 0L, d$y)
  d
}
d <- ord_data(20260930)

rows <- list(
  list("cumulative logit, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(cumulative()), quote(brms::cumulative())),
  list("cumulative probit, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(cumulative("probit")), quote(brms::cumulative("probit"))),
  list("sratio logit, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(sratio()), quote(brms::sratio())),
  list("cratio cloglog, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(cratio("cloglog")), quote(brms::cratio("cloglog"))),
  list("cratio probit, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(cratio("probit")), quote(brms::cratio("probit"))),
  list("acat logit, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(acat()), quote(brms::acat())),
  list("sratio, cs(z) and disc ~ 0 + x", list(y ~ x + cs(z), disc ~ 0 + x),
       quote(sratio()), quote(brms::sratio())),
  list("acat, cs(z) and disc ~ 0 + x", list(y ~ x + cs(z), disc ~ 0 + x),
       quote(acat()), quote(brms::acat())),
  list("cumulative, disc ~ 0 + z, link_disc softplus",
       list(y ~ x, disc ~ 0 + z),
       quote(cumulative(link_disc = "softplus")),
       quote(brms::cumulative(link_disc = "softplus"))),
  list("cumulative equidistant", list(y ~ x),
       quote(cumulative(threshold = "equidistant")),
       quote(brms::cumulative(threshold = "equidistant"))),
  list("cumulative probit equidistant, disc ~ 0 + z",
       list(y ~ x, disc ~ 0 + z),
       quote(cumulative("probit", threshold = "equidistant")),
       quote(brms::cumulative("probit", threshold = "equidistant"))),
  list("sratio equidistant, cs(z)", list(y ~ x + cs(z)),
       quote(sratio(threshold = "equidistant")),
       quote(brms::sratio(threshold = "equidistant"))),
  list("cratio equidistant", list(y ~ x),
       quote(cratio(threshold = "equidistant")),
       quote(brms::cratio(threshold = "equidistant"))),
  list("acat equidistant, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(acat(threshold = "equidistant")),
       quote(brms::acat(threshold = "equidistant"))),
  list("cumulative sum_to_zero", list(y ~ x),
       quote(cumulative(threshold = "sum_to_zero")),
       quote(brms::cumulative(threshold = "sum_to_zero"))),
  list("cumulative sum_to_zero, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(cumulative(threshold = "sum_to_zero")),
       quote(brms::cumulative(threshold = "sum_to_zero"))),
  list("cumulative probit sum_to_zero", list(y ~ x),
       quote(cumulative("probit", threshold = "sum_to_zero")),
       quote(brms::cumulative("probit", threshold = "sum_to_zero"))),
  list("sratio sum_to_zero, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(sratio(threshold = "sum_to_zero")),
       quote(brms::sratio(threshold = "sum_to_zero"))),
  list("cratio probit sum_to_zero", list(y ~ x),
       quote(cratio("probit", threshold = "sum_to_zero")),
       quote(brms::cratio("probit", threshold = "sum_to_zero"))),
  list("acat sum_to_zero, cs(z)", list(y ~ x + cs(z)),
       quote(acat(threshold = "sum_to_zero")),
       quote(brms::acat(threshold = "sum_to_zero"))),
  list("cumulative thres(gr = g) equidistant", list(y | thres(gr = g) ~ x),
       quote(cumulative(threshold = "equidistant")),
       quote(brms::cumulative(threshold = "equidistant"))),
  list("sratio thres(gr = g) sum_to_zero, disc ~ 0 + z",
       list(y | thres(gr = g) ~ x, disc ~ 0 + z),
       quote(sratio(threshold = "sum_to_zero")),
       quote(brms::sratio(threshold = "sum_to_zero"))),
  list("acat thres(gr = g), disc ~ 0 + z",
       list(y | thres(gr = g) ~ x, disc ~ 0 + z),
       quote(acat()), quote(brms::acat())),
  list("acat probit, disc ~ 0 + z", list(y ~ x, disc ~ 0 + z),
       quote(acat("probit")), quote(brms::acat("probit"))),
  list("acat cloglog, cs(z)", list(y ~ x + cs(z)),
       quote(acat("cloglog")), quote(brms::acat("cloglog"))),
  list("acat cauchit equidistant", list(y ~ x),
       quote(acat("cauchit", threshold = "equidistant")),
       quote(brms::acat("cauchit", threshold = "equidistant"))),
  list("acat probit_approx sum_to_zero", list(y ~ x),
       quote(acat("probit_approx", threshold = "sum_to_zero")),
       quote(brms::acat("probit_approx", threshold = "sum_to_zero"))),
  # brms's softit helper divides a vector by a vector with `/`, which
  # the installed stanc refuses, so softit is checked on the R side
  # alone (test-ordinal-disc-thres.R)
  list("acat softit thres(gr = g), disc ~ 0 + z",
       list(y | thres(gr = g) ~ x, disc ~ 0 + z),
       quote(acat("softit")), quote(brms::acat("softit"))),
  list("acat probit thres(gr = g), disc ~ 0 + z",
       list(y | thres(gr = g) ~ x, disc ~ 0 + z),
       quote(acat("probit")), quote(brms::acat("probit"))),
  list("hurdle_cumulative probit equidistant", list(yh ~ x),
       quote(hurdle_cumulative("probit", threshold = "equidistant")),
       quote(brms::hurdle_cumulative("probit", threshold = "equidistant"))),
  list("hurdle_cumulative sum_to_zero, hu ~ x", list(yh ~ x, hu ~ x),
       quote(hurdle_cumulative(threshold = "sum_to_zero")),
       quote(brms::hurdle_cumulative(threshold = "sum_to_zero")))
)
res <- data.frame(row = character(), logLik = numeric(), const = numeric(),
                  rel = numeric(), grad = numeric(), stringsAsFactors = FALSE)
for (r in rows) {
  cat("\n==", r[[1]], "==\n")
  out <- tryCatch({
    fit <- frm(do.call(bf, r[[2]]), family = eval(r[[3]]), data = d,
               control = lp_tight)
    bf_brms <- do.call(brms::bf, r[[2]])
    chk <- env$brms_lp_check(bf_brms, eval(r[[4]]), d, fit)
    c(logLik = chk$ours, const = chk$measured_const,
      rel = abs(chk$measured_const) / max(1, abs(chk$ours)),
      grad = chk$max_grad)
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    c(logLik = NA, const = NA, rel = NA, grad = NA)
  })
  res[nrow(res) + 1L, ] <- list(r[[1]], out[["logLik"]], out[["const"]],
                                out[["rel"]], out[["grad"]])
}
cat("\n== SUMMARY (seed 20260930, n = 300) ==\n")
print(format(res, digits = 4), row.names = FALSE)
