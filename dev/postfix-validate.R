# Lane postfix: the pooled methods of a frm_multiple() result against
# mice::pool() and emmeans's own mira pooling of glm()/lm() fits of the
# same imputations.
#   Rscript dev/postfix-validate.R <lib>
# Output: dev/postfix-log/validate.txt (by redirection). Seeds in line.
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(mice)
  library(emmeans)
})
cat("frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
pr <- function(label, v) cat(sprintf("%-58s %.3e\n", label, v))

## 1. Poisson, where the ML covariance of glm() is the one frm() has
set.seed(11)
n <- 120
d <- data.frame(x1 = rnorm(n), x2 = factor(sample(c("a", "b", "c"), n,
                                                  TRUE)))
d$y <- rpois(n, exp(0.4 + 0.5 * d$x1 + c(a = 0, b = 0.3, c = -0.2)[d$x2]))
d$x1[sample(n, 30)] <- NA
imp <- mice(d, m = 5, print = FALSE, seed = 12)
fm <- suppressMessages(frm_multiple(y ~ x1 + x2, family = poisson(),
                                    data = imp))
mi <- with(imp, glm(y ~ x1 + x2, family = poisson()))
pl <- summary(pool(mi))
fx <- fixef(fm)
cat("\n## poisson, m = 5, n = 120 (seed 11, mice seed 12)\n")
pr("fixef Estimate vs pool()", rel(fx[, "Estimate"], pl$estimate))
pr("fixef Est.Error vs pool() std.error", rel(fx[, "Est.Error"],
                                              pl$std.error))
pd <- pool(mi)$pooled
pr("pooled df vs pool() df", rel(fm$pooled$df, pd$df))
# the total covariance by hand from glm(): W + (1 + 1/m) B
gl <- getfit(mi)
W <- Reduce(`+`, lapply(gl, vcov)) / 5
Q <- sapply(gl, coef)
Tv <- W + (1 + 1 / 5) * cov(t(Q))
V <- vcov(fm)
pr("vcov() vs glm() total covariance (all entries)",
   max(abs(V - Tv)) / max(abs(Tv)))
pr("diag(vcov()) vs pool() t", rel(diag(V), pd$t))
# fitted() on newdata against Rubin over predict.glm(type = "link"),
# the scale fitted() pools on since fix round 1, carried back by exp()
nd <- data.frame(x1 = c(-1, 0, 1.5), x2 = factor(c("a", "b", "c"),
                                                 levels = c("a", "b", "c")))
pg <- lapply(gl, function(g) predict(g, nd, type = "link",
                                     se.fit = TRUE))
qe <- sapply(pg, `[[`, "fit")
ue <- sapply(pg, function(p) p$se.fit^2)
pp <- lapply(seq_len(nrow(nd)), function(i) {
  # n - k is the complete-data df, df.residual() of the fit
  pool.scalar(qe[i, ], ue[i, ], n = 120, k = 4)
})
qbar <- vapply(pp, `[[`, 0, "qbar")
tse <- sqrt(vapply(pp, `[[`, 0, "t"))
fi <- fitted(fm, newdata = nd)
pr("fitted(newdata) Estimate vs exp(mice::pool.scalar() of eta)",
   rel(fi[, "Estimate"], exp(qbar)))
pr("fitted(newdata) Est.Error vs exp(qbar) x pooled se of eta",
   rel(fi[, "Est.Error"], exp(qbar) * tse))
pr("fitted(newdata) Q97.5 vs exp(qbar + t(df) se), mice's df",
   rel(fi[, "Q97.5"], exp(qbar + qt(0.975, vapply(pp, `[[`, 0, "df")) *
                              tse)))
cat("  pooled lower bounds (all > 0):",
    format(fi[, "Q2.5"], digits = 4), "\n")
# emmeans of the pooled fit against emmeans of the mira
e1 <- as.data.frame(emmeans(fm, ~ x2))
e2 <- as.data.frame(emmeans(mi, ~ x2))
pr("emmeans(~ x2) emmean vs emmeans(mira)", rel(e1$emmean, e2$emmean))
pr("emmeans(~ x2) SE vs emmeans(mira)", rel(e1$SE, e2$SE))
pr("emmeans(~ x2) df vs emmeans(mira)", rel(e1$df, e2$df))
c1 <- as.data.frame(pairs(emmeans(fm, ~ x2)))
c2 <- as.data.frame(pairs(emmeans(mi, ~ x2)))
pr("pairs() estimate vs emmeans(mira)", rel(c1$estimate, c2$estimate))
pr("pairs() SE vs emmeans(mira)", rel(c1$SE, c2$SE))
pr("pairs() df vs emmeans(mira)", rel(c1$df, c2$df))
print(cbind(frm = e1[, c("emmean", "SE", "df")],
            mira = e2[, c("emmean", "SE", "df")]))

## 2. gaussian on nhanes: lm()'s covariance rescaled to the ML one
cat("\n## gaussian, nhanes, m = 5 (mice seed 1)\n")
imp2 <- mice(nhanes, m = 5, print = FALSE, seed = 1)
fg <- suppressMessages(frm_multiple(bmi ~ age * chl, data = imp2))
ml <- getfit(with(imp2, lm(bmi ~ age * chl)))
nn <- nrow(nhanes)
pp <- 4
Wml <- Reduce(`+`, lapply(ml, function(f) vcov(f) * (nn - pp) / nn)) / 5
Qg <- sapply(ml, coef)
Tg <- Wml + (1 + 1 / 5) * cov(t(Qg))
pr("fixef Estimate vs lm()", rel(fixef(fg)[, "Estimate"], rowMeans(Qg)))
pr("vcov() vs lm() total covariance at the ML variance",
   max(abs(vcov(fg) - Tg)) / max(abs(Tg)))
pr("nobs()", abs(nobs(fg) - 25))

## 3. a mixed model: internal identities and Rubin by hand
cat("\n## gaussian with (1 + x | g), 20 groups, hot-deck imputations,",
    "m = 5 (seed 21)\n")
set.seed(21)
nm <- 300
g <- factor(rep(1:20, length.out = nm))
u0 <- rnorm(20, 0, 0.8)
u1 <- rnorm(20, 0, 0.4)
x <- rnorm(nm)
y <- 1 + u0[g] + (0.5 + u1[g]) * x + rnorm(nm)
x[sample(nm, 60)] <- NA
dl <- lapply(1:5, function(i) {
  xi <- x
  xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
  data.frame(y = y, x = xi, g = g)
})
fmm <- suppressWarnings(suppressMessages(
  frm_multiple(y ~ x + (1 + x | g), data = dl)))
re <- ranef(fmm, condVar = TRUE)
per <- lapply(fmm$fits, ranef, condVar = TRUE)
E <- sapply(per, function(r) r$g[, 1])
S <- sapply(per, function(r) attr(r$g, "condSD")[, 1])
pr("ranef() mode vs mean of the modes", rel(re$g[, 1], rowMeans(E)))
pr("ranef(condVar) condSD vs Rubin by hand",
   rel(attr(re$g, "condSD")[, 1],
       sqrt(rowMeans(S^2) + (1 + 1 / 5) * apply(E, 1, var))))
cf <- coef(fmm)
pr("coef() = fixef() + ranef() (mean of sums)",
   rel(cf$g[["Intercept"]], fixef(fmm)["Intercept", "Estimate"] +
         re$g[, 1]))
vc <- VarCorr(fmm)
pv <- fmm$pooled_varcorr
pr("VarCorr() sds vs the pooled_varcorr table (log scale)",
   rel(vc$g$sd[, "Estimate"], pv$estimate[pv$type == "sd"]))
pr("VarCorr() sd 2.5% vs pooled_varcorr lwr (same t quantile)",
   rel(vc$g$sd[, "Q2.5"], pv$lwr[pv$type == "sd"]))
pr("VarCorr() cor vs the pooled_varcorr table (Fisher z)",
   rel(vc$g$cor[1, "Estimate", 2], pv$estimate[pv$type == "cor"]))
pr("VarCorr() cor 97.5% vs pooled_varcorr upr",
   rel(vc$g$cor[1, "Q97.5", 2], pv$upr[pv$type == "cor"]))
# Rubin by hand on the log scale over each fit's VarCorr()
per_vc <- lapply(fmm$fits, VarCorr)
L <- sapply(per_vc, function(v) log(v$g$sd[, "Estimate"]))
Ul <- sapply(per_vc, function(v) (v$g$sd[, "Est.Error"] /
                                    v$g$sd[, "Estimate"])^2)
pr("VarCorr() sd vs Rubin by hand on the log scale",
   rel(vc$g$sd[, "Estimate"], exp(rowMeans(L))))
pr("VarCorr() sd Est.Error vs the delta method by hand",
   rel(vc$g$sd[, "Est.Error"],
       exp(rowMeans(L)) * sqrt(rowMeans(Ul) + 1.2 * apply(L, 1, var))))
print(vc$g$sd)
print(vc$g$cor[1, , 2])
print(pv)
