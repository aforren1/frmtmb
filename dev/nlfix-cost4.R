# The saddle table's rows that fell with the stop at two agreeing restarts: each
# restart's value, code and warnings. Data of dev/nlfix-saddle.R, seeds 9202, 9205.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
ns <- asNamespace("frmtmb")
ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
each <- function(lab, f, fam, d) {
  cat("---", lab, "\n")
  at0 <- suppressWarnings(frm(f, family = fam, data = d, control = ctl))
  sd <- ns$nl_saddle_find(at0, use_curvature = FALSE)$saddle
  fit1 <- suppressWarnings(suppressMessages(frm(f, family = fam, data = d)))
  cat(sprintf("frm() whole: ll %.4f code %d\n", as.numeric(logLik(fit1)), fit1$opt$convergence))
  for (st in ns$nl_saddle_starts(at0$frame, sd, NULL)) {
    w <- character(0)
    r <- withCallingHandlers(suppressMessages(
      frm(f, family = fam, data = d, start = st)),
      warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 30)); invokeRestart("muffleWarning") })
    cat(sprintf("start %s: ll %.4f code %d beta %s flat %s warn[%s]\n",
                paste(round(unlist(st$beta), 3), collapse = ","), as.numeric(logLik(r)),
                r$opt$convergence, paste(round(r$estimates$beta, 3), collapse = ","),
                ns$nl_flat_warning(lapply(w, simpleWarning)), paste(unique(w), collapse = " | ")))
  }
}
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
set.seed(9205)
d$yg <- 1.5 * exp(0.6 * d$x) + rnorm(n, 0, 0.5)
d$yc <- rpois(n, 1.5 * exp(0.6 * d$x))
d$yzi <- ifelse(runif(n) < 0.25, 0L, rpois(n, 1.5 * exp(0.6 * d$x)))
d$yhu <- ifelse(runif(n) < 0.25, 0L, 1L + rpois(n, 1.5 * exp(0.6 * d$x)))
d$ypo <- rgamma(n, 3, 3 / (1.5 * exp(0.6 * d$x)))
d$ybi <- rbinom(n, 1, plogis(1.5 * exp(0.6 * d$x) - 2))
f <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
each("cratio", f, cratio(), d)
each("sratio", f, sratio(), d)
each("cumulative", f, cumulative(), d)
each("acat", f, acat(), d)
fb <- bf(ybi ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
each("bernoulli", fb, bernoulli(), d)
