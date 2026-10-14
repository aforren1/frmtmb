# Reviewer of lane postfix: guard-absent cases of the refusals.
# Usage: Rscript dev/postfix-rev-refuse.R <lib or "base">
a <- commandArgs(TRUE)
lib <- if (length(a) && a[1] != "base") a[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
try1 <- function(lbl, expr) {
  r <- tryCatch(withCallingHandlers({
    v <- expr; paste("OK", class(v)[1])
  }, message = function(m) {
    cat("  [message]", conditionMessage(m), "\n"); invokeRestart("muffleMessage")
  }, warning = function(w) {
    cat("  [warning]", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-50s %s\n", lbl, substr(r, 1, 300)))
}
set.seed(1)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- 1 + dd$x + rnorm(60)
dd$o <- cut(dd$y, c(-Inf, 0, 1, 2, Inf), ordered_result = TRUE)
f <- frm(bf(y ~ x), family = gaussian(), data = dd)
try1("summary(f, waic = FALSE)", summary(f, waic = FALSE))
try1("summary(f, waic = TRUE)", summary(f, waic = TRUE))
try1("summary(f, wai = TRUE)", summary(f, wai = TRUE))
try1("summary(f, waic = NA)", summary(f, waic = NA))
try1("summary(f, waic = 'no')", summary(f, waic = "no"))
try1("frm(o ~ x, cumulative(), threshold = 'flexible')",
     frm(bf(o ~ x), family = cumulative(), data = dd, threshold = "flexible"))
try1("frm(o ~ x, cumulative(), thres = 'equidistant')",
     frm(bf(o ~ x), family = cumulative(), data = dd, thres = "equidistant"))
try1("frm(o ~ x, cumulative(), thresh = 'equidistant')",
     frm(bf(o ~ x), family = cumulative(), data = dd, thresh = "equidistant"))
try1("frm(o ~ x, cumulative(), th = 'equidistant')",
     frm(bf(o ~ x), family = cumulative(), data = dd, th = "equidistant"))
try1("frm(o ~ x, cumulative(), threshold = NULL)",
     frm(bf(o ~ x), family = cumulative(), data = dd, threshold = NULL))
try1("frm(bf(o ~ x) + cumulative(), threshold = 'eq')",
     frm(bf(o ~ x) + cumulative(), data = dd, threshold = "equidistant"))
try1("do.call(frm, threshold)",
     do.call(frm, list(bf(o ~ x), family = cumulative(), data = dd,
                       threshold = "equidistant")))
fam <- sratio()
try1("frm(family = fam (a variable), threshold =)",
     frm(bf(o ~ x), family = fam, data = dd, threshold = "equidistant"))
## fitted/residuals ndraws
try1("fitted(f, ndraws = 10)", fitted(f, ndraws = 10))
try1("fitted(f, ndraws = 0)", fitted(f, ndraws = 0))
try1("fitted(f, ndraws = 1.5)", fitted(f, ndraws = 1.5))
try1("fitted(f, NULL, NULL, 'response', NULL, NULL, NULL, 10)",
     fitted(f, NULL, NULL, "response", NULL, NULL, NULL, 10))
try1("residuals(f, ndraws = 10)", residuals(f, ndraws = 10))
try1("predict(f, ndraws = 10)", predict(f, ndraws = 10))
try1("posterior_epred(f, ndraws = 10)", posterior_epred(f, ndraws = 10))
try1("fitted(f, draw_ids = 1)", fitted(f, draw_ids = 1))
## core log_lik with newdata on an ML fit
nd <- dd[1:5, ]
try1("log_lik(f, newdata = nd)", log_lik(f, newdata = nd))
## frm_multiple: a factor level absent from one imputation's rows but
## kept in its levels
if (exists("frm_multiple")) {
  imps <- lapply(1:3, function(i) {
    d <- dd
    d$f <- factor(rep(c("a", "b", "c"), 20), levels = c("a", "b", "c"))
    d$y <- d$y + rnorm(60, 0, 0.1)
    d
  })
  imps[[2]]$f[imps[[2]]$f == "c"] <- "b"   # level c empty, still a level
  try1("frm_multiple, level 'c' empty in imputation 2",
       frm_multiple(bf(y ~ x + f), family = gaussian(), data = imps))
  f2 <- frm(bf(y ~ x + f), family = gaussian(), data = imps[[2]])
  cat("one fit on imputation 2: fixef rows:",
      paste(rownames(fixef(f2)), collapse = ","), "\n")
  print(fixef(f2))
  ## a group level empty in imputation 2 (factor level kept)
  imps2 <- imps
  imps2[[2]]$g[imps2[[2]]$g == "6"] <- "5"
  fm2 <- frm_multiple(bf(y ~ x + (1 | g)), family = gaussian(), data = imps2)
  try1("ranef(), group level 6 empty in imp 2", ranef(fm2))
  try1("fitted(newdata), group level 6 empty in imp 2",
       fitted(fm2, newdata = dd[1:3, ]))
  try1("fitted(newdata, re_formula = NA), level empty in imp 2",
       fitted(fm2, newdata = dd[1:3, ], re_formula = NA))
  try1("VarCorr(), level empty in imp 2", VarCorr(fm2))
  try1("vcov(), level empty in imp 2", vcov(fm2))
  fmg <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps)
  try1("summary(fm, waic = FALSE)", summary(fmg, waic = FALSE))
  try1("summary(fm, waic = TRUE)", summary(fmg, waic = TRUE))
  try1("fitted(fm, scale = 'linear')", fitted(fmg, scale = "linear"))
  try1("fitted(fm, probs = c(.1,.5,.9))", fitted(fmg, probs = c(.1, .5, .9)))
  try1("vcov(fm, pars = 'x')", vcov(fmg, pars = "x"))
  try1("predict(fm, newdata = nd)", predict(fmg, newdata = nd))
  try1("AIC(fm)", AIC(fmg))
  try1("BIC(fm)", BIC(fmg))
  try1("nobs(fm)", nobs(fmg))
  try1("family(fm)", family(fmg))
  try1("formula(fm)", formula(fmg))
  try1("residuals(fm)", residuals(fmg))
  try1("fitted(fm, newdata=nd, allow_new_levels)", fitted(fmg, newdata = nd))
  try1("frm_multiple(..., thres = 1)",
       frm_multiple(bf(y ~ x), family = gaussian(), data = imps, thres = 1))
  ## categorical: fitted in [0,1]
  imo <- lapply(1:3, function(i) {
    d <- dd; d$o <- cut(d$y + rnorm(60, 0, .3), c(-Inf, 0, 1, 2, Inf),
                        ordered_result = TRUE); d })
  fmo <- frm_multiple(bf(o ~ x), family = cumulative(), data = imo)
  fo <- fitted(fmo, newdata = data.frame(x = c(-3, 0, 3)))
  cat("cumulative pooled fitted dim", dim(fo), " range of bounds",
      range(fo[, 3:4, ]), "\n")
}
