# Reviewer, re-check: step()'s "object 'change' not found" on a plain
# frmtmb fit, per arm; stepAIC() on data without NA; BIC through
# insight::get_loglikelihood(); frm_bootstrap() row counts.
#   Rscript dev/aterms2-rev-22-step.R <base|lane>
# Seed 2101 data (no NA). Log: dev/aterms2-rev-log-22-step.txt
arm <- commandArgs(TRUE)[[1]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
tr <- function(e) tryCatch(q(e), error = function(e) paste("ERROR:", conditionMessage(e)))
set.seed(2101)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                s = rep(c(TRUE, TRUE, FALSE, FALSE), length.out = n))
d$y <- 1 + d$x + 0.5 * d$z + rnorm(n)
g <- q(frm(y ~ x + z + w, data = d[d$s, ]))
r <- tr(stats::step(g, trace = 0))
cat(arm, "step(g):", if (is.character(r)) r else deparse1(formula(r)), "\n")
r <- tr(MASS::stepAIC(g, trace = 0))
cat(arm, "stepAIC(g):", if (is.character(r)) r else deparse1(formula(r)), "\n")
if (arm == "lane") {
  f <- q(frm(y | subset(s) ~ x + z + w, data = d))
  r <- tr(MASS::stepAIC(f, trace = 0))
  cat(arm, "stepAIC(f):", if (is.character(r)) r else deparse1(formula(r)), "\n")
  r1 <- tr(MASS::stepAIC(f, trace = 0, k = log(nobs(f))))
  r2 <- tr(MASS::stepAIC(g, trace = 0, k = log(nobs(g))))
  cat(arm, "stepAIC(k = log(nobs)) f:", if (is.character(r1)) r1 else deparse1(formula(r1)),
      " g:", if (is.character(r2)) r2 else deparse1(formula(r2)), "\n")
  cat(arm, "BIC(f)", BIC(f), " BIC(insight::get_loglikelihood(f))",
      BIC(insight::get_loglikelihood(f)), " BIC(g)", BIC(g),
      " BIC(insight::get_loglikelihood(g))", BIC(insight::get_loglikelihood(g)), "\n")
  b <- tr(frm_bootstrap(f, nsim = 3, seed = 1))
  bg <- tr(frm_bootstrap(g, nsim = 3, seed = 1))
  cat(arm, "frm_bootstrap identical f vs g:", identical(b, bg), " class", class(b)[1], " dim",
      paste(dim(b), collapse = "x"), "\n")
  d1 <- as.data.frame(q(drop1(f, test = "Chisq")))
  d2 <- as.data.frame(q(drop1(g, test = "Chisq")))
  cat(arm, "drop1 identical f vs g:", identical(d1, d2), "\n")
  print(d1)
  cat(arm, "drop1 numeric columns identical:",
      identical(unname(as.matrix(d1)), unname(as.matrix(d2))),
      " all.equal:", paste(all.equal(d1, d2), collapse = "; "), "\n")
}
