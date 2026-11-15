# Review of lane priors2, item 6: the open chart where the box is not
# the plain case: a user lb below the support, a closed upper end beside
# an open lower one, and a badly scaled column that autoscale rescales.
# Rscript priors2-review-chart.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(seq_len(n / 12), each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(n / 12)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
d$zs <- d$z * 1e-4
one <- function(lab, f, pr, ...) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(frm(f, data = d, prior = pr, ...),
                                    warning = function(x) {
                                      w <<- c(w, conditionMessage(x))
                                      invokeRestart("muffleWarning")
                                    }),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%-40s ERROR %s", lab, substr(conditionMessage(r), 1, 150)),
        "\n")
    return(invisible())
  }
  b <- r$estimates$beta
  cat(sprintf("%-40s code %d obj %.12g slope %.10g", lab,
              r$opt$convergence, r$opt$objective, b[length(b)]), "\n")
  if (length(w)) cat("    warn:", substr(w, 1, 110), sep = "\n    ")
}
ln <- function(...) set_prior("lognormal(0, 1)", class = "b", ...)
one("lb = -1 (below support)", yneg ~ z, ln(coef = "z", lb = -1))
one("lb = 0", yneg ~ z, ln(coef = "z", lb = 0))
one("lb = 0, ub = 5 (closed upper)", yneg ~ z, ln(coef = "z", lb = 0, ub = 5))
one("lb = 0, ub = 0.004 (upper binds)", yneg ~ z,
    ln(coef = "z", lb = 0, ub = 0.004))
one("lb = 1e-8, ub = 0.004", yneg ~ z, ln(coef = "z", lb = 1e-8, ub = 0.004))
# autoscale: the column is 1e-4 of z, so the slope is 1e4 times larger;
# lognormal(9, 1) puts the prior on that scale
lz <- function(...) set_prior("lognormal(9, 1)", class = "b", coef = "zs",
                              ...)
one("zs lb = 0, autoscale default", yneg ~ zs, lz(lb = 0))
one("zs lb = 0, autoscale FALSE", yneg ~ zs, lz(lb = 0),
    control = frmtmb_control(autoscale = FALSE))
one("zs lb = 1e-4, autoscale default", yneg ~ zs, lz(lb = 1e-4))
one("zs lb = 1e-4, autoscale FALSE", yneg ~ zs, lz(lb = 1e-4),
    control = frmtmb_control(autoscale = FALSE))
# the class-wide b bound with the open end on two slopes at once
d$yneg2 <- d$yneg - 0.6 * d$x
one("two slopes lb = 0 class-wide", yneg2 ~ z + x,
    set_prior("lognormal(0, 1)", class = "b", lb = 0))
one("two slopes lb = 1e-8 class-wide", yneg2 ~ z + x,
    set_prior("lognormal(0, 1)", class = "b", lb = 1e-8))
