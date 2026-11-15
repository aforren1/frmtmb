# Lane priors2: the four interior maxima of dev/priors2-survey.R that the
# open chart moved, against a tight optimum: nlminb with rel.tol 1e-15
# from the box's own point on the lane build's objective (the objective
# is the same function on both builds). Prints, per case, the distance
# of each build's point to the tight one, relative to max(1, |p|), and
# the objective gaps.
# Rscript dev/priors2-survey-tight.R > dev/priors2-out/survey-tight.txt
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
rd <- function(f) {
  x <- strsplit(readLines(f)[-1L], "\t", fixed = TRUE)
  stats::setNames(lapply(x, function(s) {
    as.numeric(strsplit(s[5L], ",", fixed = TRUE)[[1L]])
  }), vapply(x, `[`, "", 1L))
}
pb <- rd("dev/priors2-out/survey-before.tsv")
pa <- rd("dev/priors2-out/survey-after.tsv")
set.seed(1)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:15, 10)),
                t = factor(rep(1:5, 30)))
d$y <- 1 + 0.4 * d$x - 0.3 * d$z + rnorm(15)[d$g] * 0.5 + rnorm(n, sd = 0.7)
cases <- list(
  b_lognormal_in = set_prior("lognormal(0, 1)", class = "b", coef = "x",
                             lb = 0),
  b_gamma2_in = set_prior("gamma(2, 2)", class = "b", coef = "x", lb = 0),
  b_beta_in = set_prior("beta(2, 2)", class = "b", coef = "x", lb = 0,
                        ub = 1),
  b_invgamma_in = set_prior("inv_gamma(3, 1)", class = "b", coef = "x",
                            lb = 0))
for (nm in names(cases)) {
  f <- suppressWarnings(frm(y ~ x + z, data = d, prior = cases[[nm]]))
  o <- stats::nlminb(pb[[nm]], f$obj$fn, f$obj$gr,
                     control = list(rel.tol = 1e-15, x.tol = 1e-15,
                                    eval.max = 1000, iter.max = 1000))
  rel <- function(p) max(abs(p - o$par) / pmax(1, abs(o$par)))
  cat(sprintf(paste("%-15s tight obj %.17g; box (0.82.0) %.3g off,",
                    "objective +%.3g; open chart %.3g off, objective",
                    "+%.3g\n"), nm, o$objective, rel(pb[[nm]]),
              f$obj$fn(pb[[nm]]) - o$objective, rel(pa[[nm]]),
              f$obj$fn(pa[[nm]]) - o$objective))
}
