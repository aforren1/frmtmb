# What brms 2.23.0 generates for gp(x, by = ) in its several forms.
# No compile: stancode(), standata(), default_prior() only.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 40
dd <- data.frame(x = runif(n, 0, 5), z = runif(n, 0, 3),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 w = runif(n, 0.5, 2))
dd$y <- sin(dd$x) + rnorm(n, 0, 0.3)
show <- function(form, ...) {
  cat("\n==========", deparse(form), "\n")
  sc <- stancode(form, data = dd, ...)
  print(sc)
  sd <- standata(form, data = dd, ...)
  str(sd[grepl("gp|Igp|Cgp|Jgp", names(sd))], give.attr = FALSE,
      vec.len = 4)
  print(default_prior(form, data = dd, ...))
}
show(y ~ gp(x, by = f))
show(y ~ gp(x, by = f, gr = FALSE))
show(y ~ gp(x, by = w))
show(y ~ gp(x, z, by = f))
show(y ~ gp(x, z, by = f, iso = FALSE))
show(y ~ gp(x, by = f, k = 8))
show(y ~ gp(x, by = f, cmc = FALSE))
show(y ~ gp(x, z))
