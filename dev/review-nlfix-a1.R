# Area 1: every reader of an ordinal nl fit against the linear spelling.
source("dev/review-nlfix-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(4101)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:15, each = 20)))
d$w <- rep(1:20, 15)
lat <- 0.8 * d$x + 0.5 * d$x * d$z + rlogis(n)
d$y <- cut(lat, c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
nd <- d[1:7, ]
nd$x <- nd$x + 1
cmp <- function(lab, a, b) {
  if (iserr(a) || iserr(b)) {
    cat(sprintf("  %-22s nl=%s | lin=%s\n", lab,
        if (iserr(a)) paste("ERR", substr(a$msg, 1, 200)) else "ok",
        if (iserr(b)) paste("ERR", substr(b$msg, 1, 200)) else "ok"))
  } else {
    a <- suppressWarnings(as.numeric(unlist(a))); b <- suppressWarnings(as.numeric(unlist(b)))
    if (length(a) != length(b)) {
      cat(sprintf("  %-22s lengths %d %d\n", lab, length(a), length(b)))
      return(invisible())
    }
    cat(sprintf("  %-22s max|nl-lin| %.3e  (range nl %.3g..%.3g)\n", lab,
                max(abs(a - b), na.rm = TRUE), min(a, na.rm = TRUE), max(a, na.rm = TRUE)))
  }
}
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  F <- get(fam)
  cat("\n=====", fam, "\n")
  fn <- tt(frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = F(), data = d))
  fl <- tt(frm(y ~ 0 + x, family = F(), data = d))
  cat("  logLik nl lin gap", as.numeric(logLik(fn)), as.numeric(logLik(fl)),
      as.numeric(logLik(fn)) - as.numeric(logLik(fl)), "\n")
  cat("  estimates nl:", format(unlist(fn$estimates), digits = 8), "\n")
  rd <- list(
    fitted = function(f) fitted(f),
    fitted_nd = function(f) fitted(f, newdata = nd),
    predict = function(f) { set.seed(5); predict(f) },
    predict_nd = function(f) { set.seed(5); predict(f, newdata = nd) },
    linpred = function(f) frm_linpred(f),
    linpred_nd = function(f) frm_linpred(f, newdata = nd),
    residuals = function(f) residuals(f, type = "deviance"),
    resid_osa = function(f) { set.seed(1); residuals(f, type = "osa") },
    resid_ord = function(f) residuals(f, type = "ordinary"),
    loglik_pw = function(f) { x <- tryCatch(frmtmb:::pointwise_loglik(f), error = function(e) NULL); if (is.null(x)) stop("no pointwise") else x },
    sim = function(f) { set.seed(1); s <- simulate(f, nsim = 4); as.numeric(unlist(lapply(s, as.integer))) },
    curve = function(f) { x <- frm_curve(f, newdata = nd); as.data.frame(x)$estimate },
    curve_ci = function(f) { x <- as.data.frame(frm_curve(f, newdata = nd)); x$conf.high %||% x$upper %||% x[[ncol(x)]] },
    vcov = function(f) vcov(f),
    emmeans = function(f) { e <- emmeans::ref_grid(f, at = list(x = c(0, 1))); e@bhat },
    emmeans2 = function(f) { e <- emmeans::emmeans(f, ~ x, at = list(x = c(0, 1))); as.data.frame(e)$emmean },
    me = function(f) { m <- marginaleffects::predictions(f, newdata = nd); m$estimate },
    me_se = function(f) { m <- marginaleffects::predictions(f, newdata = nd); m$std.error },
    me_slopes = function(f) { m <- marginaleffects::avg_slopes(f, variables = "x"); m$estimate },
    confint = function(f) confint(f),
    update = function(f) { u <- update(f, data = d[1:250, ]); as.numeric(logLik(u)) },
    refit = function(f) { u <- update(f); as.numeric(logLik(u)) },
    condeff = function(f) { ce <- conditional_effects(f, effects = "x"); nrow(ce[[1]]) },
    condeff_pred = function(f) { ce <- conditional_effects(f, effects = "x", method = "predict"); nrow(ce[[1]]) },
    pp_check_data = function(f) { s <- posterior_predict(f, nsim = 2); dim(s) }
  )
  for (nm in names(rd)) {
    a <- tt(rd[[nm]](fn)); b <- tt(rd[[nm]](fl))
    cmp(nm, a, b)
    if (length(attr(a, "w"))) cat("     nl warn:", substr(attr(a, "w")[1], 1, 180), "\n")
  }
}
