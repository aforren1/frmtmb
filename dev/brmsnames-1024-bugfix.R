## The two answers the parser's alias got wrong beyond names, per arm:
## drop1() and step() on y ~ exp(x) + z (seed 3, n = 200), and
## fitted(re_formula = ~ (1 + exp(x) | g)) on
## y ~ exp(x) + (1 + exp(x) | g) (seed 77, n = 300). The constructions
## are the review's (dev/reviews/2026-10-24-brmsnames.md, finding 1),
## and are the tests' in tests/testthat/test-brms-names-specials.R.
##
##   Rscript dev/brmsnames-1024-bugfix.R base|lane > dev/brmsnames-1024-log/bugfix-<arm>.txt
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(3)
n <- 200
d <- data.frame(x = runif(n), z = rnorm(n))
d$y <- 1 + 0.8 * exp(d$x) + 0.05 * d$z + rnorm(n)
fit <- q(frm(y ~ exp(x) + z, data = d))
print(q(drop1(fit)))
cat("AIC of y ~ z refitted:", format(AIC(q(frm(y ~ z, data = d))),
                                      digits = 8), "\n")
cat("AIC of y ~ exp(x) refitted:", format(AIC(q(frm(y ~ exp(x), data = d))),
                                          digits = 8), "\n")
s <- q(step(fit, trace = 0))
cat("step() keeps:", rownames(fixef(s)), "\n\n")

set.seed(77)
n <- 300
d <- data.frame(x = runif(n, 0.2, 2), g = factor(rep(1:20, each = 15)))
d$y <- 1 + 0.2 * exp(d$x) + rnorm(20, 0, 0.4)[d$g] +
  rnorm(20, 0, 0.5)[d$g] * exp(d$x) + rnorm(n, 0, 0.3)
d$xx <- exp(d$x)
nd <- d[c(1, 20, 40), ]
fit <- q(frm(y ~ exp(x) + (1 + exp(x) | g), data = d))
full <- fitted(fit, newdata = nd)[, 1]
same <- tryCatch(fitted(fit, newdata = nd,
                        re_formula = ~ (1 + exp(x) | g))[, 1],
                 error = function(e) conditionMessage(e))
ref <- fitted(q(frm(y ~ xx + (1 + xx | g), data = d)), newdata = nd)[, 1]
cat("fitted(), re_formula = NULL:            ", format(full, digits = 10),
    "\n")
cat("fitted(), re_formula = ~ (1 + exp(x) | g):",
    if (is.character(same)) same else format(same, digits = 10), "\n")
cat("the same model on a column xx = exp(x): ", format(ref, digits = 10),
    "\n")
