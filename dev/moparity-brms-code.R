# What brms 2.23.0 writes for mo(x) * f and (mo(x) | g): Stan code,
# Stan data and default priors. No compile.
# Usage: Rscript dev/moparity-brms-code.R > dev/moparity-log/brms-code.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 120
dd <- data.frame(
  x = sample(0:3, n, TRUE),
  f = factor(sample(c("a", "b", "c"), n, TRUE)),
  f2 = factor(sample(c("a", "b"), n, TRUE)),
  g = factor(sample(paste0("g", 1:8), n, TRUE)),
  z = rnorm(n)
)
dd$y <- rnorm(n)
show <- function(fo, ...) {
  cat("\n==========", deparse1(fo$formula), "\n")
  tryCatch(show0(fo, ...), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    print(conditionCall(e))
  })
}
show0 <- function(fo, ...) {
  print(default_prior(fo, data = dd, ...))
  cat(stancode(fo, data = dd, ...))
  sd <- standata(fo, data = dd, ...)
  str(unclass(sd)[setdiff(names(sd), c("Y"))], vec.len = 6)
}
show(bf(y ~ mo(x) * f))
show(bf(y ~ mo(x) * f2))
show(bf(y ~ mo(x):f))
show(bf(y ~ mo(x):f2))
show(bf(y ~ f + mo(x):f))
show(bf(y ~ (mo(x) | g)))
show(bf(y ~ mo(x) + (mo(x) | g)))
show(bf(y ~ mo(x) + (1 + mo(x) | g)))
show(bf(y ~ mo(x) * z + (mo(x) | g)))
show(bf(y ~ mo(x) + (mo(x) || g)))
show(bf(y ~ mo(x) + (mo(x) + z | g)))
