# logistic_normal() and dirichlet() inside mvbf() on lane catmulti's
# build: the fit against the two univariate fits it factors into (no
# residual correlation, so the log-likelihood is their sum), and every
# post-fit route that reads the logistic-normal correlations.
# Arg 1: "base" or "lane". Output: dev/catmulti-out/mvln-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/wt-catmulti-lib"
base <- "C:/Users/adf44/source/r/rellib-r9"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(lib, base, user) else c(base, user))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
source("dev/catmulti-data.R")
tr <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}
sec <- function(t) cat("\n========", t, "========\n")
dd <- catmulti_mv_data(13)
# a second composition for two logistic_normal() responses
set.seed(131)
G <- matrix(stats::rgamma(nrow(dd) * 3, shape = 4), nrow(dd), 3)
W <- G / rowSums(G)
colnames(W) <- c("u", "v", "w")
dd$W <- W
ctl <- frmtmb_control(grad_tol = 1e-6, restarts = 3)
sec("mv logistic_normal + gaussian")
fm <- tr(frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                  bf(w ~ x, family = gaussian())), data = dd, control = ctl))
print(if (is.character(fm)) fm else "fitted")
if (!is.character(fm)) {
  f1 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd, control = ctl)
  f2 <- frm(bf(w ~ x), family = gaussian(), data = dd, control = ctl)
  cat(sprintf("logLik mv %.10f  sum of univariate %.10f  rel diff %.3g\n",
              as.numeric(logLik(fm)),
              as.numeric(logLik(f1)) + as.numeric(logLik(f2)),
              abs(as.numeric(logLik(fm)) - as.numeric(logLik(f1)) -
                    as.numeric(logLik(f2))) / abs(as.numeric(logLik(fm)))))
  print(summary(fm)$spec_pars)
  print(summary(f1)$spec_pars)
  print(tr(variables(fm)))
  print(tr(hypothesis(fm, "lncor_Y__b__c > 0", class = NULL)))
  print(tr(hypothesis(f1, "lncor__b__c > 0", class = NULL)))
  print(tr(default_prior(fm)))
  fp <- tr(frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                    bf(w ~ x, family = gaussian())), data = dd,
               prior = set_prior("lkj(2)", class = "lncor", resp = "Y")))
  print(tr(summary(fp)$spec_pars))
  fpw <- tr(frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                     bf(w ~ x, family = gaussian())), data = dd,
                prior = set_prior("lkj(2)", class = "lncor", resp = "w")))
  print(if (is.character(fpw)) fpw else "no refusal for resp = w")
  set.seed(1)
  s <- tr(simulate(fm, nsim = 2))
  print(tr(str(s)))
  print(tr(dim(predict(fm, resp = "Y"))))
  print(tr(confint(fm)))
  ds <- tr(suppressWarnings(frm_sample(fm, chains = 1, iter = 200,
                                       refresh = 0, seed = 2)))
  if (!is.character(ds)) {
    print(grep("lncor", variables(ds), value = TRUE))
    sm <- tr(summary(ds)); print(class(sm)); print(tr(names(sm)))
    ds1 <- tr(suppressWarnings(frm_sample(f1, chains = 1, iter = 200, refresh = 0, seed = 2)))
    sm1 <- tr(summary(ds1)); print(class(sm1)); print(tr(names(sm1)))
    print(tr(posterior_summary(ds)[grep("lncor", variables(ds)), , drop = FALSE]))
  } else print(ds)
}
sec("two logistic_normal responses")
f3 <- tr(frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                  bf(W ~ x, family = logistic_normal())), data = dd,
             control = ctl))
if (!is.character(f3)) {
  print(summary(f3)$spec_pars)
  fw <- frm(bf(W ~ x), family = logistic_normal(), data = dd, control = ctl)
  f1 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd, control = ctl)
  cat(sprintf("rel diff logLik %.3g\n",
              abs(as.numeric(logLik(f3)) - as.numeric(logLik(f1)) -
                    as.numeric(logLik(fw))) / abs(as.numeric(logLik(f3)))))
  print(tr(default_prior(f3)))
} else print(f3)
sec("mv dirichlet + gaussian against the brms fit")
fd <- tr(frm(mvbf(bf(Y ~ x, family = dirichlet()),
                  bf(w ~ x, family = gaussian())), data = dd, control = ctl))
if (!is.character(fd)) {
  bd <- readRDS("dev/catmulti-out/brms-mv-dir.rds")
  fb <- brms::fixef(bd)
  ff <- fixef(fd)
  nm <- rownames(fb)
  cat("names identical:", identical(rownames(ff), nm), "\n")
  print(round(cbind(frmtmb = ff[nm, "Estimate"], brms = fb[nm, "Estimate"],
                    z = (ff[nm, "Estimate"] - fb[nm, "Estimate"]) /
                      fb[nm, "Est.Error"],
                    se_ratio = ff[nm, "Est.Error"] / fb[nm, "Est.Error"]), 4))
  sb <- summary(bd)$spec_pars
  print(sb)
  print(summary(fd)$spec_pars)
  print(dim(fitted(fd)))
  print(dimnames(fitted(fd))[[3]])
}
