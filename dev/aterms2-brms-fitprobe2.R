# second half of aterms2-brms-fitprobe.R with idx values that match; seed 12
# What brms 2.23.0 does after the fit for subset(), mi(idx =) and rate()
# models. Seed 11. Output: dev/aterms2-log-brms-fitprobe.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
try_show <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.array(r) || is.matrix(r)) {
    cat("dim:", paste(dim(r), collapse = " x "), "\n")
    cat("dimnames[[length]]:", paste(utils::head(dimnames(r)[[length(dim(r))]]),
                                    collapse = " "), "\n")
  } else if (is.list(r) && !is.data.frame(r)) {
    str(r, max.level = 1)
  } else {
    print(utils::head(r))
  }
  invisible(r)
}
# mi idx
set.seed(12)
dm <- data.frame(g1 = sample(seq(1, 59, 2), 60, TRUE), g2 = 1:60, s = rep(c(TRUE, FALSE), 30))
dm$x <- rnorm(60)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + rnorm(60, sd = 0.5)
dm$x[c(3, 7, 8)] <- NA
bm <- bf(y ~ mi(x, idx = g1)) + bf(x | mi() + index(g2) + subset(s) ~ 1) +
  set_rescor(FALSE)
fm <- brm(bm, dm, chains = 1, iter = 300, refresh = 0, seed = 12)
try_show("mi idx fitted y", fitted(fm, resp = "y"))
try_show("mi idx fitted x", fitted(fm, resp = "x"))
try_show("mi idx log_lik", log_lik(fm))
try_show("mi idx nobs", nobs(fm))
try_show("mi idx fitted newdata y", fitted(fm, newdata = dm[1:10, ], resp = "y"))
try_show("mi idx variables", grep("Ymi", variables(fm), value = TRUE))

# rate
dr <- data.frame(y = rpois(30, 3), x = rnorm(30), time = rep(1:3, 10))
fr <- brm(y | rate(time) ~ x, dr, poisson(), chains = 1, iter = 300,
          refresh = 0, seed = 13)
try_show("rate fitted", fitted(fr))
try_show("rate fitted/denom check",
         range(fitted(fr)[, 1] / (fitted(fr, dpar = "mu")[, 1] * dr$time)))
try_show("rate linpred", posterior_linpred(fr)[1, 1:3])
try_show("rate newdata without time", fitted(fr, newdata = dr[1:3, "x", drop = FALSE]))
try_show("rate newdata", fitted(fr, newdata = dr[1:3, ]))
try_show("rate conditional_effects",
         conditional_effects(fr)$x[1:2, c("x", "time", "estimate__")])
