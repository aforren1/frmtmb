# Reviewer of lane ordmix: odds and ends on the lane build.
# (a) an ordinal-mixture response in a multivariate model (brms builds it)
# (b) cens() and trunc() on an ordinal mixture (brms refuses both)
# (c) theta1 ~ z alone with three components (brms refuses)
# (d) incl_thres on frmtmb.sample's posterior_linpred() is lane fixes';
#     here: what posterior_linpred() does on an ordinal mixture
# Seed 20261095.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261095)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), ce = sample(0:1, n, TRUE))
lat <- ifelse(rbinom(n, 1, 0.4) == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
  rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
d$g <- d$z + rnorm(n)
tr <- function(lab, expr) {
  r <- tryCatch({
    f <- suppressWarnings(expr)
    sprintf("FITS logLik %.6f", as.numeric(logLik(f)))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-26s %s\n", lab, substr(r, 1, 300)))
}
tr("(a) mvbf ordinal mixture",
   frm(mvbf(bf(y ~ x) + mixture(cumulative(), cumulative()),
            bf(g ~ z) + gaussian()), data = d))
tr("(a2) mvbf plain cumulative",
   frm(mvbf(bf(y ~ x) + cumulative(), bf(g ~ z) + gaussian()), data = d))
tr("(b) cens()",
   frm(bf(y | cens(ce) ~ x), family = mixture(cumulative(), sratio()),
       data = d))
tr("(b2) trunc()",
   frm(bf(y | trunc(ub = 4) ~ x), family = mixture(cumulative(), sratio()),
       data = d))
tr("(c) theta1 ~ z, K = 3",
   frm(bf(y ~ x, theta1 ~ z), family = mixture(cumulative(), sratio(),
                                               acat()), data = d))
tr("(e) mixture of 3 cumulatives",
   frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                   cumulative()), data = d))
f <- frm(bf(y ~ x), family = mixture(cumulative(), sratio()), data = d)
for (dp in c("mu1", "mu2", "theta1")) {
  r <- tryCatch(head(as.numeric(predict(f, dpar = dp)), 2),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("(d) predict(dpar =", dp, "):", substr(paste(r, collapse = " "), 1,
                                              200), "\n")
}
cat("(f) print:\n")
print(f)
cat("(g) mixture_probs head:\n")
print(head(mixture_probs(f), 3))

## (h) the multivariate fit is the sum of its two univariate fits: the
## responses share no parameter and no residual correlation
fm <- suppressWarnings(frm(mvbf(bf(y ~ x) + mixture(cumulative(), cumulative()),
                                bf(g ~ z) + gaussian()), data = d))
f_y <- suppressWarnings(frm(bf(y ~ x), family = mixture(cumulative(),
                                                        cumulative()),
                            data = d))
f_g <- frm(bf(g ~ z), family = gaussian(), data = d)
cat(sprintf("(h) mv logLik %.10f; univariate sum %.10f; diff %.3g\n",
            as.numeric(logLik(fm)),
            as.numeric(logLik(f_y)) + as.numeric(logLik(f_g)),
            as.numeric(logLik(fm)) - as.numeric(logLik(f_y)) -
              as.numeric(logLik(f_g))))
cat("(h) mv fixef rows:", paste(rownames(fixef(fm)), collapse = " "), "\n")
P1 <- fitted(fm, resp = "y")
P2 <- fitted(f_y)
cat(sprintf("(h) fitted(resp = 'y') vs univariate: max|diff| %.3g\n",
            max(abs(P1[, "Estimate", ] - P2[, "Estimate", ]))))
