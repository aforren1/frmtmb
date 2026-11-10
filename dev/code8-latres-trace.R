# The cancelling-roots arma fit: which tier decides its SEs on the lane.
source("dev/code8-common.R")
ns <- asNamespace("frmtmb")
trace("se_ridge_ritz", tracer = quote(cat("  [ritz in] S:", S, " ev:",
  format(e$values, digits = 3), "\n")), exit = quote({
  rv <- returnValue(); cat("  [ritz] values:", format(rv$values, digits = 3),
                           " noise:", format(rv$noise[1], digits = 3), "\n")
  print(round(rv$vectors, 4)) }), where = ns, print = FALSE)
trace("se_tier3", tracer = quote(cat("  [tier3 in] diagH", format(diag(H),
  digits = 3), "\n")), exit = quote({ rv <- returnValue()
  cat("  [se_tier3] lost:", paste(names(rv$lost), rv$lost), "\n") }),
  where = ns, print = FALSE)
trace("cov_from_hessian", exit = quote({ rv <- returnValue()
  cat("  [cov_from_hessian] tier", rv$tier, "lost", names(rv$lost), "\n") }),
  where = ns, print = FALSE)
set.seed(2L)
da <- expand.grid(t = 1:6, g = factor(1:10))
da$x <- stats::rnorm(60L)
da$y <- stats::rpois(60L, exp(0.5 + 0.3 * da$x +
                                rep(stats::rnorm(10L, 0, 0.5), each = 6L)))
r <- rfit(bf(y ~ x + arma(t, g, cov = TRUE)), family = poisson(), data = da)
h <- fns$fit_outer_hessian(r$f)
print(signif(h$H, 3)); print(signif(h$E, 3))
