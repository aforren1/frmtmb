# Test edits for test-autocor.R: compare gls()/lme() sigma, the marginal
# residual sd, with frmtmb's sigma (the innovation sd) times the
# marginal factor.
f <- "tests/testthat/test-autocor.R"
x <- paste(readLines(f), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop("not once: ", old, " (", n, ")")
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("# The standard errors of ac_pars(), by the delta method",
"# The marginal residual sd over sigma (the innovation sd under ar(),
# ma() and arma(), as in brms): the square root of the process variance
# at unit innovation variance, summed from the psi weights of
# stats::ARMAtoMA(), so independent of frmtmb's own solve. gls() and
# lme() report the marginal sd as their sigma.
ac_marg <- function(fit) {
  ac <- fit$frame$autocor[[1L]]
  nat <- ac_pars(fit)
  ar <- nat[seq_len(ac$p)]
  ma <- nat[ac$p + seq_len(ac$q)]
  psi <- c(1, stats::ARMAtoMA(ar = ar, ma = ma, lag.max = 5000))
  sqrt(sum(psi^2))
}

# The standard errors of ac_pars(), by the delta method")
rep1("    expect_within_se(log(sigma(f)), log(g$sigma),",
     "    expect_within_se(log(sigma(f) * ac_marg(f)), log(g$sigma),")
rep1("    # the sigma dpar is log-linked, so exp() of its coefficients are the
    # two residual SDs varIdent reports as sigma times its ratios
    expect_equal(sort(exp(unname(f$estimates$betad))),",
"    # the sigma dpar is log-linked, so exp() of its coefficients are the
    # two innovation SDs; times the marginal factor they are the two
    # residual SDs varIdent reports as sigma times its ratios
    expect_equal(sort(exp(unname(f$estimates$betad)) * ac_marg(f)),")
rep1("    expect_equal(sigma(f), g$sigma, tolerance = 1e-5)",
     "    expect_equal(sigma(f) * ac_marg(f), g$sigma, tolerance = 1e-5)")
rep1("  expect_equal(stats::sd(as.vector(r)), sigma(f), tolerance = 0.05)",
"  # the marginal sd, sigma (the innovation sd) times the factor
  expect_equal(stats::sd(as.vector(r)), sigma(f) * ac_marg(f),
               tolerance = 0.05)")
rep1("  # pearson divides by the marginal SD, which is sigma (R is unit
  # diagonal), so it stays the plain standardization
  expect_equal(residuals(f, type = \"pearson\")[, \"Estimate\"],
               residuals(f)[, \"Estimate\"] / sigma(f), tolerance = 1e-12)",
"  # pearson divides by the marginal SD of a row, sigma (the innovation
  # SD) times the marginal factor
  expect_equal(residuals(f, type = \"pearson\")[, \"Estimate\"],
               residuals(f)[, \"Estimate\"] / (sigma(f) * ac_marg(f)),
               tolerance = 1e-12)")
writeLines(x, f)
