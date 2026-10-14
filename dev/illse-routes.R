# Lane illse: every reader of the covariance on one ill-conditioned fit
# (gaussian degree-8 raw polynomial with (1 | g), seed 3 of
# dev/illse-designs.R, and the same without the group term), against
# the same quantity on the poly(x, 8) fit of the same likelihood: the
# ratio of each SE, and the warnings at the call.
#   ILLSE_ARM=base|lane Rscript dev/illse-routes.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
suppressMessages({library(emmeans); library(marginaleffects)})
for (re in c(FALSE, TRUE)) {
  ds <- des_poly(8, re)(3)
  d <- ds$data
  raw <- cap(frm(ds$formula, data = d))
  f <- raw$v
  g <- suppressWarnings(frm(ds$poly_formula, data = d))
  cat(sprintf("\n== %s | code %d | lost %s | fit warnings %d\n", ds$lab,
              f$opt$convergence, lost_of(f), length(raw$w)))
  nd <- data.frame(x = c(1.1, 1.5, 1.9), g = d$g[1])
  q <- function(lab, expr) {
    a <- cap(eval(expr, list(m = f)))
    b <- cap(eval(expr, list(m = g)))
    v <- if (inherits(a$v, "illse_err") || inherits(b$v, "illse_err")) {
      paste("ERR", substr(c(a$v, b$v)[1], 1, 70))
    } else paste(sprintf("%.5f", a$v / b$v), collapse = " ")
    cat(sprintf("  %-34s %s | warnings %d\n", lab, v, length(a$w)))
  }
  q("fitted(re_formula = NA)", quote(
    fitted(m, newdata = nd, re_formula = NA)[, "Est.Error"]))
  if (re) q("fitted() at a group", quote(fitted(m, newdata = nd)[, "Est.Error"]))
  q("frm_linpred(se.fit)", quote(frm_linpred(m, newdata = nd,
                                            re_formula = NA,
                                            se.fit = TRUE)$se.fit))
  q("frm_lp_basis() through V_root", quote({
    b <- frm_lp_basis(m, newdata = nd, re_formula = NA)
    if (is.null(b$V_root)) sqrt(rowSums((b$A %*% b$V) * b$A)) else
      sqrt(rowSums((b$A %*% b$V_root)^2) + rowSums((b$A %*% b$V_rest) * b$A))
  }))
  q("predict() Est.Error (4000 draws)", quote({
    set.seed(1)
    predict(m, newdata = nd, re_formula = NA, ndraws = 4000)[, "Est.Error"]
  }))
  q("conditional_effects() band", quote({
    ce <- conditional_effects(m, effects = "x", re_formula = NA)[[1]]
    ce$se__[c(10, 50, 90)]
  }))
  q("emmeans() at x", quote(summary(emmeans::emmeans(
    m, ~ x, at = list(x = nd$x)))$SE))
  q("emtrends(x)", quote(summary(emmeans::emtrends(
    m, ~ x, var = "x", at = list(x = nd$x)))$SE))
  q("sigma's SE (summary)", quote(summary(m)$spec_pars[, "Est.Error"]))
  if (re) q("VarCorr() sd's SE", quote(VarCorr(m)$g$sd[, "Est.Error"]))
  if (re) q("ranef(condVar) condSD, levels 1-3", quote(
    attr(ranef(m, condVar = TRUE)$g, "condSD")[1:3, 1]))
  q("marginaleffects predictions()", quote(
    marginaleffects::predictions(m, newdata = nd,
                                 re.form = NA)$std.error[1:3]))
}
