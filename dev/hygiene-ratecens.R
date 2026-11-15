# Item 1519: rate() with cens() or trunc() against brms 2.23.0's
# compiled program, through tests/testthat/helper-brms.R's
# brms_lp_check() (check A: log_prob at frmtmb's estimates equals
# logLik; check B: brms's gradient vanishes there). Seed 1519.
#
#   Rscript dev/hygiene-ratecens.R <mine|base> > dev/hygiene-out/ratecens.log
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "mine")) "C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-hygiene"
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = file.path(WT, "dev", "stan-cache"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
H <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(WT, "tests", "testthat", "helper-brms.R"), envir = H)
brms_lp_check <- H[["brms_lp_check"]]
options(frmtmb.brms_lp_report = TRUE)
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "brms",
    as.character(packageVersion("brms")), "\n")

set.seed(1519)
n <- 300
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
dr$y <- rpois(n, exp(0.2 + 0.5 * dr$x) * dr$time)
dr$yn <- rnbinom(n, mu = exp(0.2 + 0.5 * dr$x) * dr$time,
                 size = 2 * dr$time)
# right (1) and left (-1) censoring at the observed count
dr$cc <- rep(c(0, 0, 1, -1), length.out = n)
# a truncation window that holds every count
dr$ub <- pmax(dr$y, dr$yn) + 3L

run <- function(label, expr) {
  cat("\n== ", label, "\n", sep = "")
  r <- tryCatch(
    with_reporter("check", expr),
    error = function(e) cat("ERROR", conditionMessage(e), "\n"))
  invisible(r)
}
row <- function(label, fit_f, bform, fam, d = dr) {
  run(label, {
    fit <- frm(fit_f, data = d, family = fam)
    r <- brms_lp_check(bform, fam, d, fit)
    cat(sprintf("RESULT %s const %.3g grad %.3g logLik %.6f\n", label,
                r$measured_const, r$max_grad, r$ours))
  })
}
row("poisson rate+cens", y | rate(time) + cens(cc) ~ x,
    brms::bf(y | rate(time) + cens(cc) ~ x), poisson())
row("poisson(sqrt) rate+cens", y | rate(time) + cens(cc) ~ x,
    brms::bf(y | rate(time) + cens(cc) ~ x), poisson("sqrt"))
row("negbinomial rate+cens", yn | rate(time) + cens(cc) ~ x,
    brms::bf(yn | rate(time) + cens(cc) ~ x), brms::negbinomial())
row("geometric rate+cens", y | rate(time) + cens(cc) ~ x,
    brms::bf(y | rate(time) + cens(cc) ~ x), brms::geometric())
row("poisson rate+trunc(lb=0,ub)", y | rate(time) + trunc(lb = 0, ub = ub) ~ x,
    brms::bf(y | rate(time) + trunc(lb = 0, ub = ub) ~ x), poisson())
row("poisson rate+trunc(lb=1)",
    y | rate(time) + trunc(lb = 1) ~ x,
    brms::bf(y | rate(time) + trunc(lb = 1) ~ x), poisson(),
    d = dr[dr$y >= 1, ])
row("negbinomial rate+trunc(lb=0,ub)",
    yn | rate(time) + trunc(lb = 0, ub = ub) ~ x,
    brms::bf(yn | rate(time) + trunc(lb = 0, ub = ub) ~ x),
    brms::negbinomial())
