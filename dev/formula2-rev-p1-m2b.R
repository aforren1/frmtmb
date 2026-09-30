# Reviewer, punch round 1: the cnms of fit2 in formula2-rev-p1-m12.R and a
# non-vacuous check that re_formula selected something.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(51)
n <- 360
d <- data.frame(f = factor(sample(c("a", "b", "c"), n, TRUE)),
                x = rnorm(n), h = factor(rep(1:18, each = 20)))
d$y <- rnorm(n, as.numeric(d$f) + d$x, exp(0.3 * rnorm(18)[d$h]))
fit2 <- suppressWarnings(frm(bf(y ~ x + (0 + f | h)) +
                lf(sigma ~ 1 + (0 + f | h), cmc = FALSE), data = d))
for (b in fit2$frame$re_blocks) for (cp in b$components)
  cat(cp$lp_key, ":", cp$cnms, "\n")
nd <- d[1:10, ]
for (dp in c("mu", "sigma")) {
  a <- fitted(fit2, newdata = nd, dpar = dp, re_formula = ~ (0 + f | h))
  b <- fitted(fit2, newdata = nd, dpar = dp)
  z <- fitted(fit2, newdata = nd, dpar = dp, re_formula = NA)
  cat(dp, "term == full:", isTRUE(all.equal(a, b)), " full != NA:",
      !isTRUE(all.equal(b, z)), "\n")
}
