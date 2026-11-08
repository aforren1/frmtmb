# Review of lane nlstart: the sign probe under a prior, on this lane and on
# a trial merge with lane loglik (rev-nlstart-lib/merge-loglik), where
# logLik() drops the log prior and logLik_penalized() keeps it.
# Usage: REV_ARM=after|<lib> Rscript dev/review-nlstart-loglik.R
source("dev/review-nlstart-common.R")
has_pen <- exists("logLik_penalized", asNamespace("frmtmb"))
cat("logLik_penalized:", has_pen, "\n")
n <- 300
set.seed(101)
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- 1.5 * exp(0.6 * d$x) + 0.3 * d$z + rnorm(n, 0, 0.5)
fm <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
pr <- set_prior("normal(0, 2)", nlpar = "b1")
for (reml in c(FALSE, TRUE)) {
  f <- tt(frm(fm, data = d, prior = pr, REML = reml))
  cat(sprintf("REML %s: %s\n", reml, if (iserr(f)) paste("ERROR:", f$msg) else
    sprintf("-objective %.6f logLik %s penalized %s conv %s", -f$opt$objective,
            tryCatch(format(as.numeric(logLik(f)), digits = 10),
                     error = function(e) paste("ERR", conditionMessage(e))),
            if (has_pen) format(as.numeric(logLik_penalized(f)), digits = 10)
            else "-", f$opt$convergence)))
  for (m in attr(f, "m")) cat("   MSG:", substr(m, 1, 400), "\n")
  off <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                             tt(frm(fm, data = d, prior = pr, REML = reml)))
  if (!iserr(off)) {
    cat(sprintf("   probe off: -objective %.6f\n", -off$opt$objective))
  } else cat("   probe off ERROR:", off$msg, "\n")
}
