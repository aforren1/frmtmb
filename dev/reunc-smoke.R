# Lane wt-reunc: smoke test of partial re_formula and the group-effect
# draw. Rscript dev/reunc-smoke.R [lib]
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(3)
G <- 15; m <- 8
d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                h = factor(rep(1:5, times = G * m / 5)), x = rnorm(G * m))
ug <- rnorm(G, 0, 0.8); sg <- rnorm(G, 0, 0.4); uh <- rnorm(5, 0, 0.6)
d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + uh[d$h] + rnorm(G * m)
fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d)
fe <- fixef(fit)[, "Estimate"]
re <- ranef(fit)
rg <- re$g
print(dim(re$g))
lin <- function(ig, sl, ih) {
  fe[1] + fe[2] * d$x + ig * rg[as.character(d$g), "Intercept"] +
    sl * rg[as.character(d$g), "x"] * d$x +
    ih * re$h[as.character(d$h), "Intercept"]
}
chk <- function(lab, got, want) {
  cat(sprintf("%-34s max abs diff %.3e\n", lab, max(abs(got - want))))
}
chk("NULL", frm_linpred(fit), lin(1, 1, 1))
chk("~(1 + x | g)", frm_linpred(fit, re_formula = ~(1 + x | g)),
    lin(1, 1, 0))
chk("~(1 | g)", frm_linpred(fit, re_formula = ~(1 | g)), lin(1, 0, 0))
chk("~(0 + x | g)", frm_linpred(fit, re_formula = ~(0 + x | g)),
    lin(0, 1, 0))
chk("~(1 | h)", frm_linpred(fit, re_formula = ~(1 | h)), lin(0, 0, 1))
chk("~(1 | g) + (1 | h)", frm_linpred(fit, re_formula = ~(1 | g) + (1 | h)),
    lin(1, 0, 1))
chk("~1", frm_linpred(fit, re_formula = ~1), lin(0, 0, 0))
chk("full formula identical to NULL",
    frm_linpred(fit, re_formula = ~(1 + x | g) + (1 | h)),
    frm_linpred(fit))
cat("identical full:", identical(frm_linpred(fit,
    re_formula = ~(1 + x | g) + (1 | h)), frm_linpred(fit)), "\n")
nd <- d[c(1, 2, 7, 13), ]
nd_noh <- nd; nd_noh$h <- NULL
chk("newdata ~(1 | g), no h column",
    frm_linpred(fit, newdata = nd_noh, re_formula = ~(1 | g)),
    lin(1, 0, 0)[c(1, 2, 7, 13)])
e <- tryCatch(frm_linpred(fit, re_formula = ~(1 | nosuch)),
              error = function(e) e)
cat("nosuch class:", class(e)[1:2], "\n  ", conditionMessage(e), "\n")
e <- tryCatch(fitted(fit, re_formula = ~(1 + z | g)), error = function(e) e)
cat("(1+z|g):", conditionMessage(e), "\n")
se0 <- frm_linpred(fit, se.fit = TRUE)$se.fit
se1 <- frm_linpred(fit, se.fit = TRUE, re_formula = ~(1 | g))$se.fit
seN <- frm_linpred(fit, se.fit = TRUE, re_formula = NA)$se.fit
cat("se NULL / (1|g) / NA row1:", se0[1], se1[1], seN[1], "\n")
f1 <- fitted(fit, re_formula = ~(1 | h))
chk("fitted ~(1|h) estimate", f1[, 1], lin(0, 0, 1))
# predict: draws
set.seed(9)
# propagate_error = FALSE holds the group effects too since the
# 2026-09-23 rename, so this arm is now the noise-only predictive
# distribution and its correlations are ~0 by construction
P <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000,
             propagate_error = FALSE)
print(round(cor(P), 3))
set.seed(9)
P2 <- predict(fit, newdata = nd, ndraws = 4000, re_formula = ~(1 | g))
print(P2)
