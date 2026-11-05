# The merge of this round with predread: S2. Arg: the merged library.
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(9202); n <- 300
d <- data.frame(x = round(rnorm(n), 3), g = factor(rep(1:30, each = 10)), t = rep(1:10, 30))
d$y <- ordered(cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4))
d$yn <- 1 + 0.5 * d$x + rnorm(n)
fn <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d); fl <- frm(y ~ x, family = cumulative(), data = d)
nd <- d[1:3, ]
cat("frm_linpred(dpar = mu, type = response): body", format(frm_linpred(fn, newdata = nd, dpar = "mu", type = "response"), digits = 5), "| linear", format(frm_linpred(fl, newdata = nd, dpar = "mu", type = "response"), digits = 5), "\n")
cat("plogis of the link scale:", format(plogis(frm_linpred(fn, newdata = nd, dpar = "mu", type = "link")), digits = 5), "\n")
x <- fitted(fn, newdata = nd, dpar = "mu"); y <- fitted(fl, newdata = nd, dpar = "mu"); cat("fitted(dpar = mu): body", format(x[, 1], digits = 5), "| linear", format(y[, 1], digits = 5), "\n")
c1 <- tryCatch(as.data.frame(frm_curve(fn, newdata = nd, dpar = "mu", transform = TRUE))$.estimate, error = function(e) conditionMessage(e)); c2 <- tryCatch(as.data.frame(frm_curve(fl, newdata = nd, dpar = "mu", transform = TRUE))$.estimate, error = function(e) conditionMessage(e))
cat("frm_curve(transform): body", format(c1, digits = 5), "| linear", format(c2, digits = 5), "\n")
fm <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = mixture(cumulative(), cumulative()), data = d)
cat("ordinal mixture body, dpar = mu1, response:", format(frm_linpred(fm, newdata = nd, dpar = "mu1", type = "response"), digits = 5), " link:", format(frm_linpred(fm, newdata = nd, dpar = "mu1", type = "link"), digits = 5), "\n")
set.seed(4)
fa <- suppressWarnings(frm(bf(yn ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d))
fb <- suppressWarnings(frm(bf(yn ~ x, autocor = ~ ar(time = t, gr = g, cov = FALSE)), data = d))
s1 <- rowMeans(as.matrix(simulate(fa, nsim = 2000, newdata = d[1:20, ], seed = 2))); s2 <- rowMeans(as.matrix(simulate(fb, nsim = 2000, newdata = d[1:20, ], seed = 2)))
cat("simulate(newdata), body against linear, max row-mean difference:", max(abs(s1 - s2)), "\n")
