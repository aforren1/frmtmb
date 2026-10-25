p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/tests/testthat/test-fcor.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  nu <- summary(fit)$spec_pars[\"nu\", \"Estimate\"]\n", "")
rep1("  expect_equal(unname(rowMeans(sims)), unname(frm_linpred(fit)),
               tolerance = 6 * sqrt(max(diag(S)) / 4000) /
                 mean(abs(frm_linpred(fit))))",
"  expect_lt(max(abs(rowMeans(sims) - frm_linpred(fit))),
            6 * sqrt(max(diag(S)) / 4000))")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
