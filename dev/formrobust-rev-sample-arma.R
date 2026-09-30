# Reviewer: frmtmb.sample's per-draw fill of a missing newdata response
# under ar(cov = FALSE), against brms's .predictor_arma() run on the
# SAME posterior draws. Data seed 31, sampler seed 3, 1 chain, 4000
# iter (the lane's formrobust-repro3c.R construction). RNG seed 2.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/formrobust-rev-stan-cache")
suppressMessages(library(frmtmb.sample))
pred_arma <- get(".predictor_arma", asNamespace("brms"))
set.seed(31)
G <- 30; Tn <- 8
d <- expand.grid(t = 1:Tn, g = factor(1:G))
d$x <- rnorm(nrow(d))
e <- as.vector(apply(matrix(rnorm(G * Tn), Tn, G), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$y <- 1 + 0.5 * d$x + e
ds <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ x + ar(t, g, p = 1)), family = gaussian(), data = d,
             chains = 1, iter = 4000, refresh = 0, seed = 3)))
nd <- d[d$g == "1", ]
nd$y[nd$t >= 5] <- NA
dr <- posterior::as_draws_df(ds)
cat("draw columns:", head(names(dr), 12), "\n")
S <- nrow(dr)
pick <- function(pat) {
  nm <- grep(pat, names(dr), value = TRUE)
  cat("  ", pat, "->", nm, "\n"); as.numeric(dr[[nm[1]]])
}
b0 <- pick("^b_Intercept$"); bx <- pick("^b_x$")
sig <- pick("^sigma$")
cat("names(ds):", names(ds), "\n")
fit0 <- ds[["fit"]]
ac <- fit0$frame$autocor[[1]]
th <- pick("^thetaac_1$")
ar1 <- vapply(th, function(v) frmtmb:::autocor_cond_coefs(v, ac)$ar, 0)
cat("ar posterior mean", mean(ar1), " sigma mean", mean(sig), "\n")
eta <- outer(b0, rep(1, nrow(nd))) + outer(bx, nd$x)
prep <- structure(list(family = list(fun = "gaussian"),
                       dpars = list(sigma = sig), ndraws = S,
                       data = list()), class = "brmsprep")
set.seed(2)
shifted <- pred_arma(eta, ar = matrix(ar1, S, 1), Y = nd$y,
                     J_lag = c(rep(1, nrow(nd) - 1), 0), fprep = prep)
yb <- shifted + matrix(rnorm(S * nrow(nd)), S) * sig
set.seed(2)
pf <- posterior_predict(ds, newdata = nd)
ef <- posterior_epred(ds, newdata = nd)
rows <- 4:8
cat("S =", S, "\n")
cat("predict sd, brms-on-same-draws:", format(apply(yb, 2, sd)[rows],
                                              digits = 4), "\n")
cat("predict sd, frmtmb.sample     :", format(apply(pf, 2, sd)[rows],
                                              digits = 4), "\n")
cat("ratio frmtmb/brms             :",
    format(apply(pf, 2, sd)[rows] / apply(yb, 2, sd)[rows], digits = 4),
    "\n")
cat("predict mean diff / (sd/sqrt(S)):",
    format((colMeans(pf) - colMeans(yb))[rows] /
             (apply(yb, 2, sd)[rows] / sqrt(S)), digits = 3), "\n")
cat("epred sd, brms-on-same-draws  :", format(apply(shifted, 2, sd)[rows],
                                              digits = 4), "\n")
cat("epred sd, frmtmb.sample       :", format(apply(ef, 2, sd)[rows],
                                              digits = 4), "\n")
ks <- vapply(rows, function(j) suppressWarnings(
  ks.test(pf[, j], yb[, j])$p.value), 0)
cat("KS p, predict rows 4..8:", format(ks, digits = 3), "\n")
# the mutant: the expected fill in the draws, which the lane's test
# does not catch; its row-6 sd shows what the instrument can see
src <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-mut/S01_draws_fill_expected.R"
source(src, local = new.env())
set.seed(2)
pm <- posterior_predict(ds, newdata = nd)
cat("MUTANT expected fill, predict sd:", format(apply(pm, 2, sd)[rows],
                                                digits = 4), "\n")
