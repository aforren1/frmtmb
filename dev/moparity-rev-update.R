# Reviewer: update(prior = ) merging on the lane build, across the
# spellings frm() takes for a prior.
# Usage: Rscript dev/moparity-rev-update.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(2)
n <- 80
dd <- data.frame(x = rnorm(n), z = rnorm(n))
dd$y <- 1 + 0.6 * dd$x + 0.2 * dd$z + rnorm(n)
p1 <- set_prior("normal(0, 0.05)", class = "b", coef = "x") +
  set_prior("normal(0, 10)", class = "Intercept")
fit <- frm(bf(y ~ x + z), family = gaussian(), data = dd, prior = p1)
bx <- function(f) if ("x" %in% rownames(fixef(f))) format(fixef(f)["x", 1], digits = 5) else "(no x)"
rows <- function(f) {
  if (is.null(f$prior)) return("<none>")
  paste(vapply(unclass(f$prior), function(s) {
    paste0(s$class, "/", s$coef, ":",
           s$prior %||% frmtmb:::prior_dist_string(s$dist))
  }, ""), collapse = "; ")
}
cat("base fit: b_x", bx(fit), "| prior", rows(fit), "\n")
try1 <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR", substr(conditionMessage(r), 1, 140), "\n")
  } else {
    cat(lab, ": b_x", bx(r), "| prior", rows(r), "\n")
  }
}
try1("set_prior sigma", update(fit, prior = set_prior("exponential(2)",
                                                       class = "sigma")))
try1("set_prior new b x", update(fit, prior = set_prior("normal(0, 5)",
                                                         class = "b",
                                                         coef = "x")))
try1("prior = NULL", update(fit, prior = NULL))
try1("prior() spelling", update(fit, prior = prior(exponential(2),
                                                   class = sigma)))
if (requireNamespace("brms", quietly = TRUE)) {
  bp <- brms::prior(exponential(2), class = sigma)
  try1("brms::prior() object", update(fit, prior = bp))
}
try1("unmatched new row", update(fit, prior = set_prior("normal(0, 1)",
                                                         class = "b",
                                                         coef = "nosuch")))
try1("formula drop x + sigma", update(fit, . ~ . - x,
                                      prior = set_prior("exponential(2)",
                                                        class = "sigma")))
try1("no prior arg, formula drop x", update(fit, . ~ . - x))
# a stored brmsprior
if (requireNamespace("brms", quietly = TRUE)) {
  bp2 <- brms::prior(normal(0, 0.05), class = b, coef = x)
  fitb <- frm(bf(y ~ x + z), family = gaussian(), data = dd, prior = bp2)
  cat("stored brmsprior: class of fit$prior", class(fitb$prior)[1], "\n")
  try1("  then set_prior sigma", update(fitb, prior = set_prior(
    "exponential(2)", class = "sigma")))
  try1("  then brms::prior sigma", update(fitb, prior = brms::prior(
    exponential(2), class = sigma)))
}
