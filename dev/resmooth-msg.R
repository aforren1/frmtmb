# Lane wt-resmooth, wording fix. The refusal message on the case that
# fires most, a formula naming ONE of two bar terms beside a factor
# smooth, and on the case that names one term's columns partly. Both must
# read true of the code, and both must still carry the fragment the two
# pinning tests match.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
set.seed(5)
n <- 300
d <- data.frame(x = runif(n), f = factor(rep(c("a", "b"), length.out = n)),
                h = factor(rep(1:4, length.out = n)),
                g = factor(rep(1:10, each = 30L)))
d$y <- sin(2 * pi * d$x) + rnorm(10, 0, 0.5)[d$g] * d$x +
  rnorm(4, 0, 0.4)[d$h] + rnorm(2, 0, 0.3)[d$f] + rnorm(n, 0, 0.3)
fit <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | f) +
                               (1 | h)), data = d))
show <- function(lab, expr) {
  m <- tryCatch({ force(expr); "OK (accepted)" },
                error = function(e) conditionMessage(e))
  cat("==", lab, "\n")
  cat(strwrap(m, 74, prefix = "   "), sep = "\n")
  cat("   carries 'cannot name':", grepl("cannot name", m, fixed = TRUE),
      "\n")
}
# the whole (1 | h) term is dropped: reaches inside nothing
show("re_formula = ~(1 | f)", frm_linpred(fit, re_formula = ~ (1 | f)))
# and a term whose columns are named in part
fit2 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) +
                                (1 + x | h)), data = d))
show("re_formula = ~(1 | h) on a (1 + x | h) fit",
     frm_linpred(fit2, re_formula = ~ (1 | h)))
# the shapes that must still be ACCEPTED
show("re_formula = ~(1 | f) + (1 | h)",
     frm_linpred(fit, re_formula = ~ (1 | f) + (1 | h)))
show("re_formula = NA", frm_linpred(fit, re_formula = NA))
show("re_formula = NULL", frm_linpred(fit, re_formula = NULL))
