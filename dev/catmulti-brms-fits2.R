# brms 2.23.0, lane catmulti, second batch: the linear-scale shapes of a
# mixture (several mu predictors outside the categorical families), a
# two-category categorical (one predictor), and fitted()/predict() of a
# multinomial at newdata with and without trials.
# Output: dev/catmulti-out/brms-fits2.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
out <- "dev/catmulti-out"
source("dev/catmulti-data.R")
sec <- function(t) cat("\n========", t, "========\n")
shape <- function(nm, x) {
  cat(nm, ": dim", paste(dim(x), collapse = " x "), "\n")
  dn <- dimnames(x)
  for (k in seq_along(dn)) {
    cat("  dimnames[[", k, "]]: ",
        if (is.null(dn[[k]])) "NULL" else
          paste(utils::head(dn[[k]], 6), collapse = ", "), "\n", sep = "")
  }
}
try_shape <- function(nm, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) cat(nm, ":", r, "\n") else shape(nm, r)
}
ctl <- list(chains = 2, iter = 1000, warmup = 500, refresh = 0,
            backend = "rstan", seed = 2)
sec("mixture of two gaussians")
set.seed(21)
dmx <- data.frame(x = rnorm(300))
z <- rbinom(300, 1, 0.4)
dmx$y <- ifelse(z == 1, -2 + dmx$x, 2 - 0.5 * dmx$x) + rnorm(300, 0, 0.7)
fmx <- do.call(brm, c(list(bf(y ~ x), data = dmx,
                           family = mixture(gaussian, gaussian),
                           prior = c(prior(normal(-2, 1), Intercept,
                                           dpar = mu1),
                                     prior(normal(2, 1), Intercept,
                                           dpar = mu2))), ctl))
saveRDS(fmx, file.path(out, "brms-mix.rds"))
try_shape("posterior_linpred", posterior_linpred(fmx))
try_shape("fitted(scale = 'linear')", fitted(fmx, scale = "linear"))
try_shape("posterior_linpred(transform = TRUE)",
          posterior_linpred(fmx, transform = TRUE))
sec("categorical, two categories")
dc2 <- catmulti_cat_data(14)
dc2$f <- factor(ifelse(dc2$f == "a", "a", "b"))
fc2 <- do.call(brm, c(list(bf(f ~ x), data = dc2, family = categorical()),
                      ctl))
try_shape("posterior_linpred", posterior_linpred(fc2))
try_shape("fitted(scale = 'linear')", fitted(fc2, scale = "linear"))
sec("multinomial newdata")
fm <- readRDS(file.path(out, "brms-mult-bare.rds"))
nd <- data.frame(x = c(-1, 0, 1), N = c(5L, 10L, 20L))
print(fitted(fm, newdata = nd))
try_shape("fitted(newdata without trials)",
          fitted(fm, newdata = data.frame(x = c(-1, 0, 1))))
sec("multinomial refcat y2 fixef")
fr <- readRDS(file.path(out, "brms-mult-refy2.rds"))
print(fixef(fr))
