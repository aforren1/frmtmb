# Reviewer, re-check (m8): the quadrature refusal for a dim-1 group-level mo() block.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib", "C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1); n <- 200
d <- data.frame(x = sample(0:3, n, TRUE), g = factor(sample(1:20, n, TRUE)))
d$y <- rbinom(n, 1, plogis(-0.5 + 0.4 * d$x + rnorm(20, 0, 0.5)[d$g]))
for (fo in list(y ~ mo(x) + (0 + mo(x) | g), y ~ mo(x) + (1 | g))) {
  r <- tryCatch(frm(bf(fo), family = bernoulli(), data = d, quadrature = TRUE), error = function(e) e)
  cat(deparse1(fo), ":", if (inherits(r, "error")) conditionMessage(r) else paste("fit, logLik", format(as.numeric(logLik(r)), digits = 8)), "\n")
}
