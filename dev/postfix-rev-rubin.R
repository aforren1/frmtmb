# Reviewer of lane postfix: Rubin's rules by hand from x$fits.
# Usage: Rscript dev/postfix-rev-rubin.R
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(mice)})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

## 1. mixed model, poisson, mice imputations of a covariate
set.seed(41)
n <- 200; G <- 15
g <- factor(sample(seq_len(G), n, TRUE))
u <- rnorm(G, 0, 0.5)
x <- rnorm(n); z <- rnorm(n)
y <- rpois(n, exp(0.3 + 0.4 * x - 0.3 * z + u[g]))
x[sample(n, 50)] <- NA
d <- data.frame(y, x, z, g)
imp <- mice(d, m = 5, seed = 42, print = FALSE,
            method = c("", "pmm", "", ""))
fm <- frm_multiple(y ~ x + z + (1 | g), family = poisson(), data = imp)
m <- fm$m
fx <- fixef(fm)
## vcov by hand
Q <- sapply(fm$fits, function(f) fixef(f)[, "Estimate"])
W <- Reduce(`+`, lapply(fm$fits, vcov)) / m
Tv <- W + (1 + 1/m) * cov(t(Q))
cat("mixed pois: vcov vs hand             ", rel(vcov(fm), Tv), "\n")
cat("mixed pois: sqrt diag vcov vs fixef  ",
    rel(sqrt(diag(vcov(fm))), fx[, "Est.Error"]), "\n")
## Barnard-Rubin df by hand vs fixef's df column (if any) and mice formula
dfcom <- df.residual(fm$fits[[1]])
B <- apply(Q, 1, var); Ubar <- diag(W); Tt <- Ubar + (1 + 1/m) * B
lam <- (1 + 1/m) * B / Tt
dfold <- (m - 1) / lam^2
dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
dfBR <- dfold * dfobs / (dfold + dfobs)
cat("dfcom", dfcom, "\n")
print(fm$pooled[, intersect(c("term", "estimate", "df", "fmi"), names(fm$pooled))])
cat("BR df by hand:", format(dfBR, digits = 6), "\n")
## fitted in sample: Rubin per row on the response scale
fi <- fitted(fm)
per <- lapply(fm$fits, fitted)
Qf <- sapply(per, function(p) p[, "Estimate"])
Uf <- sapply(per, function(p) p[, "Est.Error"]^2)
Tf <- rowMeans(Uf) + (1 + 1/m) * apply(Qf, 1, var)
cat("fitted Estimate vs hand              ", rel(fi[, "Estimate"], rowMeans(Qf)), "\n")
cat("fitted Est.Error vs hand             ", rel(fi[, "Est.Error"], sqrt(Tf)), "\n")
lamf <- (1 + 1/m) * apply(Qf, 1, var) / Tf
dfo <- (m - 1) / lamf^2
dfb <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lamf)
dff <- dfo * dfb / (dfo + dfb)
lo <- rowMeans(Qf) + qt(0.025, dff) * sqrt(Tf)
cat("fitted Q2.5 vs hand t(BR df)         ", rel(fi[, "Q2.5"], lo), "\n")
## one fit's fitted(): symmetric on the response scale?
p1 <- per[[1]]
cat("one fit: Q2.5 = Est - 1.96 SE ?      ",
    max(abs(p1[, "Q2.5"] - (p1[, "Estimate"] + qnorm(0.025) * p1[, "Est.Error"]))), "\n")
cat("one fit: any Q2.5 < 0 (poisson mean)  ", any(p1[, "Q2.5"] < 0), "\n")
cat("pooled: any Q2.5 < 0                  ", any(fi[, "Q2.5"] < 0), "\n")
## fitted(newdata) vs conditional_effects of the pooled fit at the same rows
ce <- conditional_effects(fm, "x")$x
nd <- ce[, c("x", "z"), drop = FALSE]
nd$g <- NA
fn <- fitted(fm, newdata = ce[, c("x", "z")], re_formula = NA)
cat("ce vs fitted(newdata, re_formula=NA): estimate rel ",
    rel(fn[, "Estimate"], ce$estimate__), "\n")
cat("ce vs fitted: lower rel ", rel(fn[, "Q2.5"], ce$lower__),
    " upper rel ", rel(fn[, "Q97.5"], ce$upper__), "\n")
print(head(cbind(ce[, c("x", "estimate__", "lower__", "upper__")],
                 fn[, c("Estimate", "Q2.5", "Q97.5")]), 3))
## ranef/coef shapes
cat("class coef(fit1):", class(coef(fm$fits[[1]])), " class $g:",
    class(coef(fm$fits[[1]])$g), "\n")
str(coef(fm$fits[[1]])$g[1:2, ])
re <- ranef(fm, condVar = TRUE)
E <- sapply(fm$fits, function(f) ranef(f)$g[, 1])
S <- sapply(fm$fits, function(f) attr(ranef(f, condVar = TRUE)$g, "condSD")[, 1])
cat("ranef mode vs hand ", rel(re$g[, 1], rowMeans(E)), "\n")
cat("ranef condSD vs hand ", rel(attr(re$g, "condSD")[, 1],
    sqrt(rowMeans(S^2) + (1 + 1/m) * apply(E, 1, var))), "\n")
cf <- coef(fm)
cat("coef = mean coef ", rel(cf$g[, 1], rowMeans(sapply(fm$fits, function(f) coef(f)$g[, 1]))), "\n")
## VarCorr: log-sd pooling by hand
vcs <- lapply(fm$fits, VarCorr)
es <- sapply(vcs, function(v) v$g$sd[, "Estimate"])
ss <- sapply(vcs, function(v) v$g$sd[, "Est.Error"])
qv <- log(es); uv <- (ss / es)^2
Tq <- mean(uv) + (1 + 1/m) * var(qv)
vc <- VarCorr(fm)
cat("VarCorr sd vs exp(mean log)", rel(vc$g$sd[, "Estimate"], exp(mean(qv))),
    " EstErr vs delta", rel(vc$g$sd[, "Est.Error"], exp(mean(qv)) * sqrt(Tq)), "\n")
print(vc$g$sd); print(es)

## 2. bernoulli: fitted response-scale pooling outside [0, 1]?
set.seed(43)
n <- 40
x <- rnorm(n); y <- rbinom(n, 1, plogis(2 + 2 * x))
x[sample(n, 10)] <- NA
imp2 <- mice(data.frame(y, x), m = 5, seed = 44, print = FALSE)
fb <- frm_multiple(y ~ x, family = bernoulli(), data = imp2)
fbn <- fitted(fb, newdata = data.frame(x = c(-2, 0, 1.5, 2.5)))
print(fbn)
fb1 <- fitted(fb$fits[[1]], newdata = data.frame(x = c(-2, 0, 1.5, 2.5)))
print(fb1)
ceb <- conditional_effects(fb, "x")$x; ceb <- ceb[c(1, 50, 100), ]; print(fitted(fb, newdata = data.frame(x = ceb$x)))
print(ceb[, c("x", "estimate__", "lower__", "upper__")])

## 3. small-sample df: gaussian n = 15 against mice::pool(lm)
set.seed(45)
n <- 15
x <- rnorm(n); y <- 1 + x + rnorm(n)
x[sample(n, 5)] <- NA
imp3 <- mice(data.frame(y, x), m = 5, seed = 46, print = FALSE)
fg <- frm_multiple(y ~ x, family = gaussian(), data = imp3)
pl <- pool(with(imp3, lm(y ~ x)))$pooled
cat("dfcom frm", df.residual(fg$fits[[1]]), " mice dfcom", pl$dfcom[1], "\n")
print(fg$pooled[, intersect(c("term", "estimate", "se", "df"), names(fg$pooled))])
print(pl[, c("term", "estimate", "df")])
fgn <- fitted(fg, newdata = data.frame(x = 0))
print(fgn)
