# Review, cross-lane: mocode's display and standard-error changes on a
# model with acprior's and loglik's residual correlation terms:
# conditional_effects(), emmeans, marginaleffects and fitted(se) of
# y ~ x + mo(k) + ar(t, g, cov = TRUE) with k from 2, and frm() with a
# class-sigma prior (the innovation reading in the MAP penalty).
#   Rscript dev/review-rel076-ce.R > dev/review-rel076-out/ce.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20261029)
d <- expand.grid(t = 1:5, g = 1:24)
d$x <- stats::rnorm(nrow(d))
d$k <- sample(2:5, nrow(d), replace = TRUE)
e <- unlist(lapply(1:24, function(g) {
  as.numeric(stats::arima.sim(list(ar = 0.6), 5, sd = 0.6))
}))
d$y <- 1 + 0.5 * d$x + c(0, 0.5, 0.8, 1.4)[d$k - 1] + e
fo <- y ~ x + mo(k) + ar(t, g, cov = TRUE)
f <- frm(fo, data = d)
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n"); print(r)
}
show("conditional_effects(f, 'k') grid", {
  ce <- conditional_effects(f, "k")[[1]]
  ce[, c("k", "estimate__", "se__")]
})
show("held value of k in conditional_effects(f, 'x')", {
  unique(conditional_effects(f, "x")[[1]]$k)
})
show("fitted(newdata, se) at k = 2:5", {
  nd <- data.frame(k = 2:5, x = 0, t = 1, g = 1)
  fitted(f, newdata = nd)
})
show("emmeans(~ x) without at (refused: k held at its mean)", {
  summary(emmeans::emmeans(f, ~ x))
})
show("emmeans(~ k, at = list(k = 2:5))", {
  summary(emmeans::emmeans(f, ~ k, at = list(k = 2:5)))
})
show("marginaleffects::predictions at k = 2:5", {
  p <- marginaleffects::predictions(f, newdata = data.frame(
    k = 2:5, x = 0, t = 1, g = 1))
  data.frame(k = 2:5, est = p$estimate, se = p$std.error)
})
# the MAP penalty under a class-sigma prior: brms's density at sigma * c
pr <- set_prior("normal(0, 0.3)", class = "sigma")
fm <- frm(fo, data = d, prior = pr)
th <- fm$estimates[["thetaac"]]
ar <- th / sqrt(1 + th^2)
s <- sigma(fm)
cat("-- frm(prior = normal(0, 0.3) on sigma): converged",
    isTRUE(fm$opt$convergence == 0), " sigma (marginal)", format(s),
    " ar", format(ar), " innovation sd", format(s * sqrt(1 - ar^2)), "\n")
cat("   ML: sigma", format(sigma(f)), "\n")
