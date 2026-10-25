p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/tests/testthat/test-weights-unbounded.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
  # the block adds nothing a row can use, so the fit is the one without
  # it, to the 1e-4 sd its unused direction keeps
  r0 <- wu_fit(bf(y | weights(w15) ~ x), family = G, data = d)
  expect_equal(unname(coef(r$fit)), unname(coef(r0$fit)), tolerance = 1e-3)",
"  # the block adds nothing a row can use, so the maximum is the one
  # without it, to the 1e-4 sd its unused direction keeps; its own sd
  # is flat there, which the fit says in its warnings (wu_fit() keeps
  # them)
  r0 <- wu_fit(bf(y | weights(w15) ~ x), family = G, data = d)
  ll <- as.numeric(logLik(r$fit))
  expect_true(is.finite(ll))
  expect_equal(ll, as.numeric(logLik(r0$fit)), tolerance = 1e-6)
  expect_equal(unname(fixef(r$fit)[, 1]), unname(fixef(r0$fit)[, 1]),
               tolerance = 1e-4)")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
