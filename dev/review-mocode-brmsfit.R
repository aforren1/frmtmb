## Reviewer probe: brms 2.23.0 fits (rstan), for the value
## conditional_effects() holds a mo() variable at, and the posterior at
## an empty category inside the range.
##   Rscript dev/review-mocode-brmsfit.R > dev/review-mocode-log/brmsfit.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
fz <- q(brm(y ~ mo(k) + z, data = d, chains = 2, iter = 2000, seed = 3,
            refresh = 0, backend = "rstan"))
ce <- q(conditional_effects(fz, effects = "z"))[[1L]]
cat("ce(z): held k =", format(unique(ce$k)), "; estimate__ at first z (",
    format(ce$z[1L]), ") =", format(ce$estimate__[1L]), "\n")
ce2 <- q(conditional_effects(fz, effects = "z:k"))[[1L]]
cat("ce(z:k): values of k =", format(unique(ce2$k)), "\n")
ce3 <- q(conditional_effects(fz, effects = "k"))[[1L]]
cat("ce(k): grid of k =", format(ce3$k), "\n")
cat("        estimate__ =", format(signif(ce3$estimate__, 6)), "\n")
r <- tryCatch(q(posterior_epred(fz, newdata = data.frame(k = 2.5, z = 0))),
              error = function(e) e)
cat("posterior_epred(k = 2.5):",
    if (inherits(r, "error")) paste("ERROR:", conditionMessage(r)) else
      format(mean(r)), "\n")
r <- tryCatch(q(posterior_epred(fz, newdata = data.frame(k = 0.6, z = 0))),
              error = function(e) e)
cat("posterior_epred(k = 0.6):",
    if (inherits(r, "error")) paste("ERROR:", conditionMessage(r)) else
      format(mean(r)), "\n")
r <- tryCatch(q(posterior_epred(fz, newdata = data.frame(k = 4.4, z = 0))),
              error = function(e) e)
cat("posterior_epred(k = 4.4):",
    if (inherits(r, "error")) paste("ERROR:", conditionMessage(r)) else
      format(mean(r)), "\n")
if (requireNamespace("emmeans", quietly = TRUE)) {
  r <- tryCatch(q(summary(emmeans::emmeans(fz, ~ z))), error = function(e) e)
  cat("emmeans(~ z):",
      if (inherits(r, "error")) paste("ERROR:", conditionMessage(r)) else
        paste(format(r$emmean), "at k =", format(r$k %||% NA)), "\n")
  rg <- tryCatch(q(emmeans::ref_grid(fz)), error = function(e) e)
  cat("ref_grid k:", if (inherits(rg, "error")) conditionMessage(rg) else
        format(rg@grid$k), "\n")
}

set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(brm(y ~ mo(x), data = dg, chains = 2, iter = 2000, seed = 3,
            refresh = 0, backend = "rstan"))
ep <- q(posterior_epred(fg, newdata = data.frame(x = 1:5)))
cat("gap fit x = 1:5: mean =", format(signif(colMeans(ep), 6)), "\n")
cat("                 sd   =", format(signif(apply(ep, 2, sd), 6)), "\n")
