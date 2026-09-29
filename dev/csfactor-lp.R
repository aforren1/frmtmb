# Log-density identity against brms for cs() on a FACTOR, one row per
# sequential ordinal family, plus a character column and a factor with
# numeric-looking levels. Reuses tests/testthat/helper-brms.R's
# brms_lp_check(): rstan::log_prob() at frmtmb's estimates against
# logLik(fit), and the gradient at that point.
#
#   Rscript dev/csfactor-lp.R <lib> > dev/csfactor-log/lp.txt
args <- commandArgs(TRUE)
LIB <- if (length(args)) args[[1L]] else
  "C:/Users/adf44/source/r/wt-csfactor-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-csfactor"
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = file.path(WT, "dev", "csfactor-stan-cache"))
# R reads only <HOME>/.R/Makevars.win and HOME here is not Documents, so
# rstan compiles without -std=gnu++17 and reports "invalid connection"
# rather than a compiler error (dev/release/run-tests.R, lines 19-27)
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
dir.create(Sys.getenv("FRMTMB_STAN_CACHE"), showWarnings = FALSE,
           recursive = TRUE)
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
# the helper reaches internals (ord_tau_from_raw and the rest), which
# testthat::test_env("frmtmb") supplies inside the suite; here the
# namespace has to be the helper environment's parent
H <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(WT, "tests", "testthat", "helper-brms.R"), envir = H)
brms_lp_check <- H[["brms_lp_check"]]
options(frmtmb.brms_lp_report = TRUE)
cat("frmtmb", as.character(packageVersion("frmtmb")),
    "brms", as.character(packageVersion("brms")), "\n")

set.seed(405)
n <- 400
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + 0.4 * x + eff)
p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(x = x, fc = fc,
                y = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$fch <- as.character(fc)
d$fnum <- factor(as.integer(fc))

row <- function(lab, bform, fform, bfam, ffam) {
  cat("\n--", lab, "--\n")
  r <- tryCatch({
    fit <- frm(fform, family = ffam, data = d)
    cat("bcs names:", paste(grep("^bcs", variables(fit), value = TRUE),
                            collapse = " "), "\n")
    res <- brms_lp_check(bform, bfam, d, fit)
    sprintf("PASS const %.3e grad %.3e logLik %.10g",
            res$measured_const, res$max_grad, res$ours)
  }, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(r, "\n")
}

row("sratio, cs(fc) factor", brms::bf(y ~ x + cs(fc)), bf(y ~ x + cs(fc)),
    brms::sratio(), sratio())
row("cratio, cs(fc) factor", brms::bf(y ~ x + cs(fc)), bf(y ~ x + cs(fc)),
    brms::cratio(), cratio())
row("acat, cs(fc) factor", brms::bf(y ~ x + cs(fc)), bf(y ~ x + cs(fc)),
    brms::acat(), acat())
row("sratio, cs(fch) character", brms::bf(y ~ x + cs(fch)),
    bf(y ~ x + cs(fch)), brms::sratio(), sratio())
row("sratio, cs(fnum) numeric-looking levels", brms::bf(y ~ x + cs(fnum)),
    bf(y ~ x + cs(fnum)), brms::sratio(), sratio())
row("sratio, cs(fc) + cs(x) two terms", brms::bf(y ~ cs(fc) + cs(x)),
    bf(y ~ cs(fc) + cs(x)), brms::sratio(), sratio())
cat("\ndone\n")
