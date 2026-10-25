p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/tests/testthat/test-weights-unbounded.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  d <- wu_data()
  d$ly <- log(d$y)
  d$w15 <- 1.5
  d$one <- 1
  d$one2 <- 1
  # the registered full-rank block is refused as the built-in one is
  e_us <- tryCatch(frm(bf(ly | weights(w15) ~ x + (1 | obs)), data = d),
                   error = conditionMessage)
  e_my <- tryCatch(frm(bf(ly | weights(w15) ~ x + myiid(1 | obs)),
                       data = d), error = conditionMessage)",
"  d <- wu_data()
  d$w15 <- 1.5
  d$one <- 1
  d$one2 <- 1
  G <- Gamma(link = \"log\")
  # the registered full-rank block is refused as the built-in one is
  e_us <- tryCatch(frm(bf(y | weights(w15) ~ x + (1 | obs)), family = G,
                       data = d), error = conditionMessage)
  e_my <- tryCatch(frm(bf(y | weights(w15) ~ x + myiid(1 | obs)),
                       family = G, data = d), error = conditionMessage)")
rep1("  r <- wu_fit(bf(ly | weights(w15) ~ x + rankone(0 + one + one2 | obs)),
              data = d)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
  r0 <- wu_fit(bf(ly | weights(w15) ~ x), data = d)
  expect_equal(sigma(r$fit), sigma(r0$fit), tolerance = 1e-3)",
"  r <- wu_fit(bf(y | weights(w15) ~ x + rankone(0 + one + one2 | obs)),
              family = G, data = d)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
  # the block adds nothing a row can use, so the fit is the one without
  # it, to the 1e-4 sd its unused direction keeps
  r0 <- wu_fit(bf(y | weights(w15) ~ x), family = G, data = d)
  expect_equal(unname(coef(r$fit)), unname(coef(r0$fit)), tolerance = 1e-3)")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
