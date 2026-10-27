## Fix round 1 check: mo() prediction standard errors against exact
## references. When the simplex is interior, y ~ mo(k) (+ z) on four
## categories is a reparametrization of y ~ factor(k) (+ z), and
## mo(k) * z of factor(k) * z, so the ML fitted values and their
## delta-method SEs equal lm()'s with the ML sigma (sqrt(RSS / n)).
## On the gap data (1, 2, 5) the observed categories are three cells of
## a model that fits them exactly, so their SE is sigma_ML / sqrt(n_x).
##   Rscript dev/review-mocode-se.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", arm, "lib:", dirname(system.file(package = "frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
ml_se <- function(m, nd) {
  p <- predict(m, newdata = nd, se.fit = TRUE)
  n <- length(residuals(m))
  list(fit = unname(p$fit),
       se = unname(p$se.fit) * sqrt(m$df.residual / n))
}
cmp <- function(lab, f, r) {
  cat(sprintf("%s\n  frmtmb est %s\n  lm     est %s\n  frmtmb se  %s\n",
              lab, paste(format(f[, 1], digits = 9), collapse = " "),
              paste(format(r$fit, digits = 9), collapse = " "),
              paste(format(f[, 2], digits = 7), collapse = " ")))
  cat(sprintf("  exact  se  %s\n  ratio      %s\n",
              paste(format(r$se, digits = 7), collapse = " "),
              paste(format(f[, 2] / r$se, digits = 6), collapse = " ")))
}

## A. seed 5, y ~ mo(k) + z
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
nd <- data.frame(k = 1:4, z = 0)
fa <- q(frm(y ~ mo(k) + z, data = d))
cmp("A. y ~ mo(k) + z, seed 5, k = 1..4, z = 0",
    q(fitted(fa, newdata = nd))[, 1:2],
    ml_se(lm(y ~ factor(k) + z, data = d), nd))
fa0 <- q(frm(y ~ mo(k), data = d))
cmp("A0. y ~ mo(k), seed 5, k = 1..4",
    q(fitted(fa0, newdata = data.frame(k = 1:4)))[, 1:2],
    ml_se(lm(y ~ factor(k), data = d), data.frame(k = 1:4)))
cat("  in-sample fitted() se, first 4 rows (k =",
    paste(d$k[1:4], collapse = " "), "):",
    paste(format(q(fitted(fa0))[1:4, 2], digits = 7), collapse = " "), "\n")
pl <- q(frm_linpred(fa0, newdata = data.frame(k = 1:4), se.fit = TRUE))
cat("  frm_linpred(se.fit) se:",
    paste(format(unname(pl$se.fit), digits = 7), collapse = " "), "\n")

## C. mo(k) * z with increasing slopes, seed 11
set.seed(11)
dc <- data.frame(k = sample(1:4, n, replace = TRUE))
dc$z <- rnorm(n)
dc$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[dc$k] +
  dc$z * c(0.2, 1, 1.5, 1.7)[dc$k] + rnorm(n)
ndc <- data.frame(k = 1:4, z = 1)
fc <- q(frm(y ~ mo(k) * z, data = dc))
cmp("C. y ~ mo(k) * z, seed 11, k = 1..4, z = 1",
    q(fitted(fc, newdata = ndc))[, 1:2],
    ml_se(lm(y ~ factor(k) * z, data = dc), ndc))

## B. the gap case, seed 7
set.seed(7)
dg <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
dg$y <- dg$x + rnorm(200)
fg <- q(frm(y ~ mo(x), data = dg))
lg <- lm(y ~ factor(x), data = dg)
r <- ml_se(lg, data.frame(x = c(1L, 2L, 5L)))
fgv <- q(fitted(fg, newdata = data.frame(x = 1:5)))
cmp("B. gap 1, 2, 5, seed 7, observed x = 1, 2, 5",
    fgv[c(1, 2, 5), 1:2], r)
cat("  x = 3, 4 Est.Error:", format(fgv[3:4, 2]), "\n")
