# Lane postfix: the nonlinear-intercept scan (nl_intercept_scan(),
# R/fit.R) on nonlinear models whose default start is outside the
# family's domain, and on models whose default start works.
#   Rscript dev/postfix-nlstart.R <arm> [lib]
# For each model: the fit from the default start, and the fit from a
# start placed by hand near the truth; the two must reach the same
# log-likelihood. Seeds in line.
a <- commandArgs(TRUE)
arm <- a[1]
lib <- if (length(a) > 1L) a[2] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("ARM", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
fit_one <- function(expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, substr(conditionMessage(x), 1, 70))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    if (grepl("scan", conditionMessage(m))) {
      w <<- c(w, paste("[message]", substr(conditionMessage(m), 120, 330)))
    }
    invokeRestart("muffleMessage")
  }),
  error = function(e) structure(conditionMessage(e), class = "err"))
  list(r = r, w = w)
}
report <- function(label, dflt, hand) {
  d <- dflt$r
  h <- hand$r
  if (inherits(d, "err")) {
    cat(sprintf("%-34s default: ERROR %s\n", label, substr(d, 1, 90)))
  } else {
    cat(sprintf("%-34s default: conv %d logLik %.10f\n", label,
                d$opt$convergence, as.numeric(logLik(d))))
  }
  if (inherits(h, "err")) {
    cat(sprintf("%-34s hand:    ERROR %s\n", "", substr(h, 1, 90)))
  } else {
    cat(sprintf("%-34s hand:    conv %d logLik %.10f\n", "",
                h$opt$convergence, as.numeric(logLik(h))))
  }
  if (!inherits(d, "err") && !inherits(h, "err")) {
    cat(sprintf("%-34s logLik(default) - logLik(hand) = %.3e\n", "",
                as.numeric(logLik(d)) - as.numeric(logLik(h))))
  }
  if (length(dflt$w)) {
    cat("    default warns:", unique(dflt$w), sep = "\n      ")
  }
  cat("\n")
}

## 1. the ledger row: update(fit2, formula. = bf(count ~ a + b, nl))
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
fit2 <- brms_fixture(2)
d2 <- brms_fixture_data(2)
report("1 brmsfit-methods:955 update",
       fit_one(update(fit2, formula. = bf(count ~ a + b, nl = TRUE))),
       fit_one(update(fit2, formula. = bf(count ~ a + b, nl = TRUE),
                      start = list(beta = c(mean(d2$count) / 2, 0,
                                            mean(d2$count) / 2, 0)))))

## 2. a + b, identified (a ~ 1, b ~ 0 + x), Gamma("identity")
set.seed(31)
n <- 200
x <- runif(n, 0, 2)
dd <- data.frame(x = x, y = rgamma(n, shape = 5, rate = 5 / (5 + 2 * x)))
f2 <- bf(y ~ a + b, a ~ 1, b ~ 0 + x, nl = TRUE)
report("2 a + b, Gamma identity",
       fit_one(frm(f2, family = Gamma("identity"), data = dd)),
       fit_one(frm(f2, family = Gamma("identity"), data = dd,
                   start = list(beta = c(4, 1)))))

## 3. b1 * exp(b2 * x), Gamma("identity")
set.seed(32)
dd3 <- data.frame(x = runif(n, 0, 2))
dd3$y <- rgamma(n, shape = 4, rate = 4 / (3 * exp(0.5 * dd3$x)))
f3 <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
report("3 b1 exp(b2 x), Gamma identity",
       fit_one(frm(f3, family = Gamma("identity"), data = dd3)),
       fit_one(frm(f3, family = Gamma("identity"), data = dd3,
                   start = list(beta = c(3, 0.5)))))

## 4. logistic growth, gaussian: the scale parameter divides
set.seed(33)
dd4 <- data.frame(x = runif(n, 0, 20))
dd4$y <- 30 / (1 + exp((10 - dd4$x) / 2)) + rnorm(n)
f4 <- bf(y ~ asym / (1 + exp((xmid - x) / scal)), asym + xmid + scal ~ 1,
         nl = TRUE)
report("4 logistic growth, gaussian",
       fit_one(frm(f4, data = dd4)),
       fit_one(frm(f4, data = dd4, start = list(beta = c(30, 10, 2)))))

## 5. log(a) + b x on poisson's log link: log(0) at the default start
set.seed(34)
dd5 <- data.frame(x = rnorm(n))
dd5$y <- rpois(n, exp(log(4) + 0.3 * dd5$x))
f5 <- bf(y ~ log(a) + b * x, a + b ~ 1, nl = TRUE)
report("5 log(a) + b x, poisson",
       fit_one(frm(f5, family = poisson(), data = dd5)),
       fit_one(frm(f5, family = poisson(), data = dd5,
                   start = list(beta = c(4, 0.3)))))

## 7. brms_nonlinear's insurance-loss curve without its priors (the
## vignette's fit_loss; vigport's workaround row "no prior, no start")
data("loss", package = "brms")
f7 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
         ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE)
report("7 brms_nonlinear fit_loss",
       fit_one(frm(f7, data = loss, family = gaussian())),
       fit_one(frm(f7, data = loss, family = gaussian(),
                   start = list(beta = c(5000, 1, 45)))))

## 6. a body whose default start works: brms_nonlinear's fit1, which
## the scan must not move (it runs only after a failure)
set.seed(35)
dd6 <- data.frame(x = rnorm(100))
dd6$y <- rnorm(100, 2 * exp(0.5 * dd6$x), 0.5)
f6 <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
r6 <- fit_one(frm(f6, data = dd6))
report("6 b1 exp(b2 x), gaussian (works)", r6,
       fit_one(frm(f6, data = dd6, start = list(beta = c(2, 0.5)))))
if (!inherits(r6$r, "err")) {
  cat("6 estimates, every digit:", sprintf("%.17g", r6$r$opt$par), "\n")
}
