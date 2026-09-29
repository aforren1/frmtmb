# Reviewer re-check, item (d): is the rewritten refusal message true of
# the code? It names "the one shape where a formula reaches inside a
# term". Two shapes reach the refusal; this asks which.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
tryv <- function(x) tryCatch(x, error = function(e) e)

set.seed(5)
n <- 300L
d <- data.frame(x = stats::runif(n),
                f = factor(rep(c("a", "b"), length.out = n)),
                h = factor(rep(1:6, length.out = n)),
                g = factor(rep(1:10, each = 30L)))
d$y <- sin(2 * pi * d$x) + stats::rnorm(10, 0, 0.5)[d$g] * d$x +
  stats::rnorm(n, 0, 0.3)

cat("== shape 1: a WHOLE term dropped, nothing reached inside ==\n")
f1 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | f) +
                               (1 | h)), family = gaussian(), data = d))
r <- tryv(frm_linpred(f1, re_formula = ~ (1 | f)))
cat("re_formula = ~(1|f) on a fit with (1|f) + (1|h):\n  ",
    if (inherits(r, "error")) conditionMessage(r) else "NO ERROR", "\n\n")

cat("== shape 2: a formula that DOES reach inside a term ==\n")
f2 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 + x | h)),
                          family = gaussian(), data = d))
r <- tryv(frm_linpred(f2, re_formula = ~ (1 | h)))
cat("re_formula = ~(1|h) on a fit with (1 + x | h):\n  ",
    if (inherits(r, "error")) conditionMessage(r) else "NO ERROR", "\n\n")

cat("== the same two shapes with NO smooth, which must NOT refuse ==\n")
f3 <- suppressWarnings(frm(bf(y ~ s(x) + (1 | f) + (1 | h)),
                          family = gaussian(), data = d))
r <- tryv(frm_linpred(f3, re_formula = ~ (1 | f)))
cat("  ~(1|f) on s(x) + (1|f) + (1|h):",
    if (inherits(r, "error")) paste("ERROR:",
                                    substr(conditionMessage(r), 1, 70)) else
      "OK", "\n")
f4 <- suppressWarnings(frm(bf(y ~ s(x) + (1 + x | h)), family = gaussian(),
                          data = d))
r <- tryv(frm_linpred(f4, re_formula = ~ (1 | h)))
cat("  ~(1|h) on s(x) + (1 + x | h):",
    if (inherits(r, "error")) paste("ERROR:",
                                    substr(conditionMessage(r), 1, 70)) else
      "OK", "\n")
