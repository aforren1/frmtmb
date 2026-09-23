# Lane wt-reunc, punch round 1: the minors that are claims about
# BEHAVIOR, checked on the build rather than taken from the round's
# notes. brms is not needed for any of these.
#
#   Rscript dev/reunc-minors.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(3)
G <- 8; m <- 6
d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                h = factor(rep(1:4, length.out = G * m)),
                x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.8)[d$g] + rnorm(G * m, 0, 0.7)
fit <- frm(bf(y ~ x + (1 | g) + (1 | h)), data = d)

cat("\nminor 2: re_formula = ~1 reads two ways\n")
s <- function(rf) {
  set.seed(1)
  if (missing(rf)) simulate(fit, nsim = 3) else
    simulate(fit, nsim = 3, re_formula = rf)
}
base <- s()
cat("  simulate(~1)       identical to NULL: ", identical(s(~1), base),
    "\n", sep = "")
cat("  simulate(~0)       identical to NULL: ", identical(s(~0), base),
    "\n", sep = "")
cat("  simulate(~(1 | g)) identical to NULL: ",
    identical(s(~(1 | g)), base), "\n", sep = "")
cat("  simulate(NA)       identical to NULL: ", identical(s(NA), base),
    "\n", sep = "")
p <- function(rf) {
  set.seed(1)
  if (missing(rf)) predict(fit, ndraws = 20) else
    predict(fit, ndraws = 20, re_formula = rf)
}
cat("  predict(~1)        identical to NA:   ",
    identical(p(~1), p(NA)), "\n", sep = "")
cat("  predict(~1)        identical to NULL: ",
    identical(p(~1), p()), "\n", sep = "")
cat("  fitted(~1) identical to NA: ",
    identical(fitted(fit, re_formula = ~1), fitted(fit, re_formula = NA)),
    "\n", sep = "")

cat("\nminor 3: a b-perturbed fit does not inherit the cache\n")
# joint_precision() memoizes into fit$cache. A list copy shares the
# environment, so a perturbed fit would have answered with the
# original's matrix. Check that the copy's cache is a DIFFERENT
# environment and that writing to it leaves the original alone.
jp <- frmtmb:::joint_precision(fit)
cat("  original memoized: ", length(ls(fit$cache)) > 0, "\n", sep = "")
g <- fit
g$estimates[["b"]][1] <- g$estimates[["b"]][1] + 1
g$cache <- new.env(parent = emptyenv())
cat("  copy's cache is a different environment: ",
    !identical(g$cache, fit$cache), "\n", sep = "")
cat("  copy's cache is empty: ", length(ls(g$cache)) == 0, "\n", sep = "")
cat("  fit_set_outer() does the same: ",
    length(ls(frmtmb:::fit_set_outer(
      fit, frmtmb:::fit_outer_vector(fit, frmtmb:::outer_par_map(fit)))$cache
    )) == 0, "\n", sep = "")

cat("\nminor 1: an unmatched term, including a nested one\n")
for (rf in list(~(1 | nosuch), ~(1 | g/h), ~(1 | g:h))) {
  e <- tryCatch({
    predict(fit, ndraws = 5, re_formula = rf)
    "accepted"
  }, frmtmb_error = function(e) paste("refused:", conditionMessage(e)))
  cat("  ", deparse(rf), " -> ", substr(e, 1, 60), "\n", sep = "")
}

cat("\nminor 5: a quadrature fit is silent on the scalar route\n")
set.seed(10)
db <- data.frame(g = factor(rep(1:10, each = 6)), x = rnorm(60))
db$y <- rbinom(60, 1, plogis(db$x + rnorm(10)[db$g]))
fq <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = db,
                           quadrature = TRUE))
fl <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = db)
nd <- db[1:2, c("x", "g")]
sc <- function(f, rf) {
  a <- if (missing(rf)) fitted(f, newdata = nd) else
    fitted(f, newdata = nd, re_formula = rf)
  round(as.numeric(a[, "Est.Error"]), 5)
}
cat("  quadrature, known level: ", paste(sc(fq), collapse = " "),
    "; re_formula = NA: ", paste(sc(fq, NA), collapse = " "), "\n", sep = "")
cat("  Laplace,    known level: ", paste(sc(fl), collapse = " "),
    "; re_formula = NA: ", paste(sc(fl, NA), collapse = " "), "\n", sep = "")
cat("  scalar route warns on the quadrature fit: ",
    tryCatch({ fitted(fq, newdata = nd); FALSE },
             warning = function(w) TRUE), "\n", sep = "")
