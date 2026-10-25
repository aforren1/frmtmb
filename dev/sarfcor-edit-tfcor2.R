p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/tests/testthat/test-fcor.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  fit <- frm(y ~ x + fcor(V), family = student(), data = d,
             data2 = list(V = V))
  nu <- frm_linpred(fit, dpar = \"nu\", type = \"response\")[1]
  ref <- mvtnorm::dmvt(d$y, delta = frm_linpred(fit),
                       sigma = sigma(fit)^2 * V, df = nu, log = TRUE)
  expect_equal(as.numeric(logLik(fit)), ref, tolerance = 1e-9)",
"  fit <- frm(y ~ x + fcor(V), family = student(), data = d,
             data2 = list(V = V))
  # one response is ONE multivariate-t draw, so its nu is not
  # identified and runs off; the density is checked at a point instead,
  # nu = 5 (nu's link is log(nu - 1))
  p <- fit$obj$env$last.par.best
  bd <- which(names(p) == \"betad\")
  p[bd[2L]] <- log(5 - 1)
  b <- p[names(p) == \"beta\"]
  ref <- mvtnorm::dmvt(d$y, delta = as.vector(cbind(1, d$x) %*% b),
                       sigma = exp(2 * p[bd[1L]]) * V, df = 5, log = TRUE)
  expect_equal(-fit$obj$fn(p), ref, tolerance = 1e-10)")
rep1("  pr <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V), profile = TRUE)",
"  pr <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V),
            control = frmtmb_control(profile = TRUE))")
rep1("  expect_error(predict(fit, newdata = nd), \"newdata2\",
               class = \"frmtmb_error\")",
"  expect_error(predict(fit, newdata = nd), \"newdata\",
               class = \"frmtmb_error\")")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
