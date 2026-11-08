# Rebuild the brief's table: logLik and AIC of y ~ x + (1 | g) under
# priors on b of shrinking scale (seed 4101). Arg 1: "before" drops the
# lane library, so the base build answers.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "after"
libs <- c("C:/Users/adf44/source/r/wt-loglik-lib",
          "C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "before") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")),
    " version:", format(packageVersion("frmtmb")), "\n")
set.seed(4101)
n <- 200; ng <- 20
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = n / ng)))
d$y <- 1 + 0.5 * d$x + rnorm(ng)[d$g] + rnorm(n)
fits <- list(none = frm(y ~ x + (1 | g), data = d))
for (s in c(100, 10, 1, 0.1)) {
  p <- set_prior(paste0("normal(0, ", s, ")"), class = "b")
  fits[[paste0("normal(0, ", s, ")")]] <- frm(y ~ x + (1 | g), data = d,
                                              prior = p)
}
for (nm in names(fits)) {
  f <- fits[[nm]]
  ll <- logLik(f)
  aic <- tryCatch(format(AIC(f), digits = 10),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  pen <- if (exists("logLik_penalized", asNamespace("frmtmb")))
    format(as.numeric(frmtmb::logLik_penalized(f)), digits = 10) else "NA"
  cat(sprintf("%-16s logLik %.10f  AIC %s  b_x %.6f  pen %s\n", nm,
              as.numeric(ll), aic, f$opt$par[2], pen))
}
# the no-prior fit's logLik at each MAP estimate: the data part the split
# must reproduce
nll0 <- fits$none$obj$fn
for (nm in names(fits)[-1]) {
  f <- fits[[nm]]
  v <- -nll0(f$opt$par)
  cat(sprintf("%-16s data logLik at MAP via no-prior tape %.10f\n", nm, v))
}
