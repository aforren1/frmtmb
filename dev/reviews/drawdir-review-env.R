# Reviewer of lane drawdir: what the libraries hold.
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("frmtmb", "frmtmb.sample", "tmbstan", "rstan", "lme4", "loo",
            "brms", "ordinal", "MASS")) {
  f <- tryCatch(find.package(p), error = function(e) NA)
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  cat(sprintf("%-15s %-10s %s\n", p, v, f))
}
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(g = factor(rep(1:15, each = 8)), x = rnorm(120))
d$y <- 1 + 0.5 * d$x + rnorm(15, 0, 0.8)[d$g] + rnorm(15, 0, 0.1)[d$g] * d$x +
  rnorm(120)
fit <- frm(bf(y ~ x + (1 + x | g)), data = d)
print(fit$opt$par)
print(names(fit$obj$env$parameters))
print(class(fit$prior)); print(fit$prior)
cat("fn at opt:", fit$obj$fn(fit$opt$par), "\n")
print(ls(fit$cache))
