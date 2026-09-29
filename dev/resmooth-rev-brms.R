# Reviewer, claims 2 and 3: brms's own rule, on a fit this review made.
# StanHeaders 2.32.10 sits in the worker's library and is first on the
# path, so rstan 2.32.7 compiles without R_MAKEVARS_USER.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms:", as.character(packageVersion("brms")),
    "| StanHeaders:", as.character(packageVersion("StanHeaders")),
    "| rstan:", as.character(packageVersion("rstan")), "\n")

set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)

cache <- "dev/resmooth-rev-brms-fits.rds"
if (file.exists(cache)) {
  fits <- readRDS(cache)
} else {
  fits <- list()
  fits$fs <- brm(y ~ s(x, g, bs = "fs", k = 5), data = d, chains = 2,
                 iter = 1000, seed = 7, refresh = 0)
  fits$byg <- brm(y ~ s(x, by = g), data = d, chains = 2, iter = 1000,
                  seed = 7, refresh = 0)
  saveRDS(fits, cache)
}

for (nm in names(fits)) {
  fit <- fits[[nm]]
  a <- posterior_epred(fit, re_formula = NA)
  b <- posterior_epred(fit, re_formula = NULL)
  cat(sprintf("%-5s identical(NA, NULL) = %-5s max|NA-NULL| = %.3e\n",
              nm, identical(a, b), max(abs(a - b))))
  # a one-sided formula on a fit with no bar term
  p <- tryCatch(posterior_epred(fit, re_formula = ~ (1 | g)),
                error = function(e) conditionMessage(e))
  cat("  re_formula = ~(1|g):",
      if (is.character(p)) paste("ERROR:", substr(p, 1, 110)) else
        sprintf("OK max|-NULL| = %.3e", max(abs(p - b))), "\n")
  # a NEW level of g
  nd <- d[1:5, ]
  nd$g <- factor("99", levels = c(levels(d$g), "99"))
  for (rf in list(NA, NULL)) {
    for (anl in c(FALSE, TRUE)) {
      r <- tryCatch(posterior_epred(fit, newdata = nd, re_formula = rf,
                                    allow_new_levels = anl),
                    error = function(e) conditionMessage(e))
      cat(sprintf("  new level g=99 re_formula=%-4s allow_new_levels=%-5s %s\n",
                  if (is.null(rf)) "NULL" else "NA", anl,
                  if (is.character(r)) paste("ERROR:", substr(r, 1, 95)) else
                    paste("OK", sprintf("%.4f", colMeans(r)[1]))))
    }
  }
  # the grouping column missing
  nd2 <- nd[, setdiff(names(nd), "g"), drop = FALSE]
  for (anl in c(FALSE, TRUE)) {
    r <- tryCatch(posterior_epred(fit, newdata = nd2, re_formula = NA,
                                  allow_new_levels = anl),
                  error = function(e) conditionMessage(e))
    cat(sprintf("  no g column, NA, allow_new_levels=%-5s %s\n", anl,
                if (is.character(r)) paste("ERROR:", substr(r, 1, 95)) else
                  "OK"))
  }
}
cat("DONE\n")
