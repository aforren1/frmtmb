# First smoke of ordinal mixtures on the lane build. Seed 20261005.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261005)
n <- 600
x <- rnorm(n)
z <- rnorm(n)
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.5 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
d <- data.frame(y, x, z)
print(table(d$y))
tryf <- function(label, expr) {
  cat("\n=====", label, "\n")
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
    print(r)
  invisible(r)
}
f1 <- tryf("none", frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative()),
                       data = d))
if (!inherits(f1, "error")) {
  print(logLik(f1))
  print(f1$estimates)
  tryf("fixef", fixef(f1))
  tryf("fitted", head(fitted(f1)))
  tryf("predict", head(predict(f1)))
  tryf("simulate", table(simulate(f1, nsim = 1, seed = 1)[[1]]))
  tryf("mixture_probs", head(mixture_probs(f1)))
  tryf("variables", variables(f1))
  tryf("default_prior", default_prior(bf(y ~ x), data = d,
                                      family = mixture(cumulative(),
                                                       cumulative())))
}
f2 <- tryf("shared", frm(bf(y ~ x), family = mixture(cumulative(),
                                                      cumulative(),
                                                      order = "mu"),
                         data = d))
if (!inherits(f2, "error")) {
  print(logLik(f2))
  print(f2$estimates)
}
