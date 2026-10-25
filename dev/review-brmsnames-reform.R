## predict(re_formula =) naming the fit's own transformed slope term.
## Seed 77.  Rscript dev/review-brmsnames-reform.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(77)
n <- 300
d <- data.frame(x = runif(n, 0.2, 2), g = factor(rep(1:20, each = 15)))
d$y <- 1 + 0.2 * exp(d$x) + rnorm(20, 0, 0.4)[d$g] +
  rnorm(20, 0, 0.5)[d$g] * exp(d$x) + rnorm(n, 0, 0.3)
d$xx <- exp(d$x)
nd <- d[c(1, 20, 40), ]
for (f in list(y ~ exp(x) + (1 + exp(x) | g), y ~ log(x) + (1 + log(x) | g),
               y ~ xx + (1 + xx | g))) {
  fit <- q(frm(f, data = d))
  rhs <- f[[3]][[3]]
  full <- fitted(fit, newdata = nd)[, 1]
  same <- tryCatch(fitted(fit, newdata = nd,
                           re_formula = stats::as.formula(call("~", rhs)))[, 1],
                   error = function(e) conditionMessage(e))
  icpt <- fitted(fit, newdata = nd, re_formula = ~ (1 | g))[, 1]
  cat(deparse1(f), "\n  full    ", format(full, digits = 10),
      "\n  same re ", if (is.character(same)) same else
        format(same, digits = 10),
      "\n  (1|g)   ", format(icpt, digits = 10), "\n")
  cat("  ranef cols:", colnames(ranef(fit)$g), "\n")
}
