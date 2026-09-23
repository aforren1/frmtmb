# brms-parity status of the drmTMB features recorded as frmtmb candidates:
# is each one a brms feature? Read from the installed brms namespace.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("brms")
cat("brms", format(packageVersion("brms")), "\n")
for (f in c("hurdle_negbinomial", "zero_one_inflated_beta", "fcor", "gr",
            "mi", "icc", "heritability", "variance_decomposition")) {
  cat(sprintf("%-24s exported: %s\n", f,
              f %in% getNamespaceExports("brms")))
}
cat("gr() arguments:", paste(names(formals(brms::gr)), collapse = ", "),
    "\n")
cat("fcor() arguments:", paste(names(formals(brms::fcor)), collapse = ", "),
    "\n")
# negbinomial truncation: brms builds Stan code with a CDF term.
d <- data.frame(y = c(1, 2, 3, 5, 8), x = rnorm(5))
sc <- brms::make_stancode(brms::bf(y | trunc(lb = 1) ~ x), data = d,
                          family = brms::negbinomial())
cat("negbinomial trunc(lb = 1) Stan code uses a CDF:",
    any(grepl("neg_binomial_2_l(c|cc)df", sc)), "\n")
# brms, like drmTMB's sd(g) ~ w, requires by to be constant within g.
sc <- brms::make_stancode(brms::bf(y ~ x + (1 | gr(g, by = f))),
                          data = data.frame(y = rnorm(20), x = rnorm(20),
                                            g = rep(1:10, each = 2),
                                            f = rep(c("a", "b"), each = 10)))
cat("gr(by =) declares one SD set per level of by (Nby_1 in Stan code):",
    any(grepl("Nby_1", sc, fixed = TRUE)), "\n")
# What fcor(M) means: the Stan lines that build the residual covariance.
M <- diag(5)
sc <- brms::make_stancode(brms::bf(y ~ x + fcor(M)), data = d,
                          data2 = list(M = M))
cat(grep("fcor|Mfcor|sigma", strsplit(sc, "\n")[[1]], value = TRUE),
    sep = "\n")
