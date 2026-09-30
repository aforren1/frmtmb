# Priors, names and other methods on the new terms. Seed 26.
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
show <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  print(r)
  invisible(r)
}
set.seed(26)
n <- 120
dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                 s = rep(c(TRUE, FALSE), n / 2), w = rnorm(n))
dm$x <- rnorm(n)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w + rnorm(n, sd = 0.5)
dm$x[c(3, 7)] <- NA
bm <- bf(y ~ mi(x, idx = g1) + w) + bf(x | mi() + index(g2) + subset(s) ~ 1)
show("default_prior", default_prior(bm, data = dm, family = gaussian()))
fm <- frm(bm, data = dm, family = gaussian())
show("variables", variables(fm))
show("prior on mixidxEQg1", fixef(frm(bm, data = dm, family = gaussian(),
                                      prior = set_prior("normal(0, 0.01)",
                                                        class = "b",
                                                        coef = "mixidxEQg1",
                                                        resp = "y"))))
show("hypothesis", hypothesis(fm, "y_mixidxEQg1 > 0"))
show("VarCorr", VarCorr(fm))
show("ce", names(conditional_effects(fm, resp = "y")))
dr <- data.frame(y = rpois(50, 3), x = rnorm(50), t = runif(50, 1, 2))
show("default_prior rate", default_prior(y | rate(t) ~ x, data = dr,
                                         family = poisson()))
show("frm_compat rate", frm_compat("rate()", "poisson"))
show("frm_compat subset", frm_compat("subset()"))
show("anova rate vs none", anova(frm(y | rate(t) ~ x, data = dr, family = poisson()),
                                 frm(y | rate(t) ~ 1, data = dr, family = poisson())))
show("bootstrap rate", {
  b <- frm_bootstrap(frm(y | rate(t) ~ x, data = dr, family = poisson()),
                     nboot = 5, seed = 1)
  class(b)
})
show("emmeans rate", {
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::emmeans(frm(y | rate(t) ~ x, data = dr, family = poisson()),
                     ~ 1)
  }
})
