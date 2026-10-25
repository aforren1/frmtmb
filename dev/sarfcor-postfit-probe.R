# The rest of the post-fit methods on fcor() and sar() fits: the
# re_formula mean, conditional_effects(), the autocor= and acformula()
# spellings, update(), AIC(), print(). Seed 9.
# Run: Rscript dev/sarfcor-postfit-probe.R > dev/sarfcor-out/postfit-probe.txt
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.spatial"),
    "\n")
set.seed(9)
g <- expand.grid(r = 1:6, c = 1:6)
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
V <- 0.5^abs(outer(1:n, 1:n, "-"))
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, 6)))
d$y <- as.vector(solve(diag(n) - 0.5 * W,
                       1 + 0.5 * d$x + rnorm(6)[d$h] + rnorm(n)))
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("---", lab, "\n")
  print(r)
}
fl <- frm(y ~ x + (1 | h) + sar(W), data = d, data2 = list(W = W))
rho <- unname(frmtmb:::autocor_natural(fl$estimates$thetaac,
                                       fl$frame$autocor[[1]]))
Ai <- solve(diag(n) - rho * W)
show("lag: fitted(re_formula = NA) equals A^-1 X b",
     max(abs(fitted(fl, re_formula = NA)[, 1] -
               as.vector(Ai %*% frm_linpred(fl, re_formula = NA)))))
show("lag: fitted() equals A^-1 (X b + Z u)",
     max(abs(fitted(fl)[, 1] - as.vector(Ai %*% frm_linpred(fl)))))
show("lag: conditional_effects() (the term dropped, as brms's incl_autocor = FALSE)",
     head(conditional_effects(fl)[[1]][, c("x", "estimate__")], 3))
show("autocor = ~ sar(W) is the same model",
     as.numeric(logLik(frm(bf(y ~ x + (1 | h), autocor = ~ sar(W)), data = d,
                           data2 = list(W = W)))) -
       as.numeric(logLik(fl)))
show("acformula(~ fcor(V)) is the same model as + fcor(V)",
     as.numeric(logLik(frm(bf(y ~ x) + acformula(~ fcor(V)), data = d,
                           data2 = list(V = V)))) -
       as.numeric(logLik(frm(y ~ x + fcor(V), data = d,
                             data2 = list(V = V)))))
show("update() of a sar fit keeps the term",
     as.numeric(logLik(update(fl, . ~ . - x))) <
       as.numeric(logLik(fl)))
show("AIC of the sar fit", AIC(fl))
print(fl)
ff <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
print(summary(ff))
show("fcor term object", fcor(V))
show("sar term object", sar(W, type = "error"))
show("fcor in a nonlinear body",
     frm(bf(y ~ a * x, a ~ 1 + fcor(V), nl = TRUE), data = d,
         data2 = list(V = V)))
show("sar crossed", frm(y ~ x * sar(W), data = d, data2 = list(W = W)))
