# Control for dev/nb2dm-priors.R: the same tables for negbinomial() and
# multinomial(), the neighbors, so a row that differs from brms there
# is not this lane's.
#   Rscript dev/nb2dm-priors-ctl.R [lib] > dev/nb2dm-out/priors-ctl.txt 2>&1
a <- commandArgs(trailingOnly = TRUE)
lib <- if (length(a)) a[1] else "C:/Users/adf44/source/r/wt-nb2dm-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r10",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "\n")
set.seed(5)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(6, 10))
d$y <- rnbinom(n, mu = exp(1 + 0.3 * d$x), size = 2)
Y <- t(sapply(seq_len(n), function(i) rmultinom(1, 12, c(0.2, 0.3, 0.5))))
colnames(Y) <- c("a", "b", "c")
d$Y <- Y
d$tr <- rowSums(Y)
print(as.data.frame(frmtmb::default_prior(y ~ x + (1 | g), d,
                                          family = frmtmb::negbinomial(),
                                          route = "sample"))[, 1:6])
print(as.data.frame(frmtmb::default_prior(Y | trials(tr) ~ x + (1 | g), d,
                                          family = frmtmb::multinomial(),
                                          route = "sample"))[, 1:6])
