# Reviewer of lane ndaterm: what brms 2.23.0 does with newdata that
# lacks an addition-term variable, per post-processing function.
# Log: dev/ndaterm-rev-out/brms.txt (redirected by the caller). Seed 7.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
tr <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) {
    structure(paste("ERROR:", conditionMessage(e)), class = "err")
  })
  if (inherits(r, "err")) {
    cat(sprintf("%-46s %s\n", label, gsub("\n", " ", substr(r, 1, 200))))
  } else {
    v <- if (is.data.frame(r)) r$estimate__ else if (is.matrix(r) && "Estimate" %in% colnames(r)) r[, "Estimate"]
         else if (is.matrix(r)) colMeans(r) else r
    cat(sprintf("%-46s OK %s\n", label,
                paste(format(head(as.numeric(v), 4), digits = 5),
                      collapse = " ")))
  }
  invisible(r)
}
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE),
                w = runif(n, 0.5, 2), s = runif(n, 0.2, 1))
d$y <- rbinom(n, d$N, plogis(0.3 + 0.5 * d$x))
d$g <- 1 + d$x + rnorm(n, 0, d$s)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
d$cc <- ifelse(d$cnt > 4, "right", "none")
d$cnt2 <- pmin(d$cnt, 4)
fit <- function(f, fam) {
  suppressMessages(brm(f, family = fam, data = d, chains = 1, iter = 400,
                       warmup = 200, refresh = 0, seed = 7))
}
fb <- fit(bf(y | trials(N) + weights(w) ~ x), binomial())
nd0 <- data.frame(x = 0)
cat("\n== binomial y | trials(N) + weights(w), newdata x only\n")
tr("posterior_predict", posterior_predict(fb, newdata = nd0))
tr("posterior_epred", posterior_epred(fb, newdata = nd0))
tr("fitted", fitted(fb, newdata = nd0))
tr("predict", predict(fb, newdata = nd0))
tr("posterior_linpred", posterior_linpred(fb, newdata = nd0))
tr("log_lik (with y, N)", log_lik(fb, newdata = data.frame(x = 0, y = 3,
                                                          N = 10)))
tr("predict, N given, w absent", predict(fb, newdata = data.frame(
  x = 0, N = 1000)))
ce <- tr("CE default", conditional_effects(fb, "x")[[1]])
if (!inherits(ce, "err")) cat("CE default grid N:", unique(ce$N), "\n")
fg <- fit(bf(g | se(s) ~ x), gaussian())
cat("\n== gaussian g | se(s), newdata x only\n")
tr("posterior_predict", posterior_predict(fg, newdata = nd0))
tr("posterior_epred", posterior_epred(fg, newdata = nd0))
tr("fitted", fitted(fg, newdata = nd0))
ce <- tr("CE default", conditional_effects(fg, "x")[[1]])
if (!inherits(ce, "err")) cat("CE default grid s:", unique(ce$s), "\n")
ce <- tr("CE predict", conditional_effects(fg, "x", method =
                                             "posterior_predict")[[1]])
fc <- fit(bf(cnt2 | cens(cc) ~ x), poisson())
cat("\n== poisson cnt2 | cens(cc), newdata x only\n")
tr("posterior_predict", posterior_predict(fc, newdata = nd0))
tr("fitted", fitted(fc, newdata = nd0))
cat("done\n")
