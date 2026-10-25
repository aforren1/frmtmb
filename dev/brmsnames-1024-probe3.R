## brms 2.23.0's default_prior() coef for the renamed columns the suite
## diff lists (dev/brmsnames-1024-suitediff.R), beside frmtmb's.
.libPaths(c("C:/Users/adf44/source/r/wt-brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
set.seed(3)
n <- 100
d <- data.frame(x = runif(n), x2 = runif(n, 1, 2), kx = 2,
                f = factor(sample(c("a", "c-d", "e"), n, TRUE)))
d$c3 <- cut(d$x, 5)
d$y <- rnorm(n)
fs <- c("y ~ splines::ns(x2, 3)", "y ~ c3", "y ~ f", "y ~ I(x * kx)",
        "y ~ x:log(x2)", "y ~ log(abs(x) + 1)", "y ~ I(x^2) + x2")
for (s in fs) {
  b <- suppressMessages(brms::default_prior(brms::bf(as.formula(s)),
                                            data = d))
  f <- frmtmb::default_prior(frmtmb::bf(as.formula(s)), data = d,
                             family = gaussian())
  cb <- b$coef[b$class == "b" & nzchar(b$coef)]
  cf <- f$coef[f$class == "b" & nzchar(f$coef)]
  cat(sprintf("%-26s brms: %s\n%-26s lane: %s  identical=%s\n", s,
              paste(cb, collapse = " "), "", paste(cf, collapse = " "),
              identical(cb, cf)))
}
