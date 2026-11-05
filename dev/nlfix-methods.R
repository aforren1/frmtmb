# Post-fit methods on an nl ordinal fit, against the same model written
# linearly (y ~ x + z): every method must give the linear fit's values, up
# to the coefficient names. Seed 49 (dev/equiv-specs.R eq_dat).
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("dev/equiv-specs.R", envir = E)
d <- E$eq_dat("ordinal", seed = 49)
set.seed(49)
d$y0 <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)),
                levels = 0:4)
nd <- d[c(3, 7, 11, 19), ]
rd <- function(a, b) {
  if (is.null(a) || is.null(b)) return("NULL value")
  a <- suppressWarnings(as.numeric(as.matrix(a)))
  b <- suppressWarnings(as.numeric(as.matrix(b)))
  if (length(a) != length(b)) return(sprintf("length %d vs %d", length(a),
                                             length(b)))
  if (anyNA(a) || anyNA(b)) return(sprintf("NA %d %d", sum(is.na(a)),
                                           sum(is.na(b))))
  sprintf("%.2e", max(abs(a - b)) / max(1, abs(b)))
}
try1 <- function(lab, fa, fb) {
  a <- tryCatch(suppressWarnings(suppressMessages(fa())),
                error = function(e) e)
  b <- tryCatch(suppressWarnings(suppressMessages(fb())),
                error = function(e) e)
  ea <- inherits(a, "error")
  eb <- inherits(b, "error")
  if (ea || eb) {
    cat(sprintf("  %-28s nl: %s | linear: %s\n", lab,
                if (ea) paste("ERR", substr(conditionMessage(a), 1, 100)) else "ok",
                if (eb) paste("ERR", substr(conditionMessage(b), 1, 60)) else "ok"))
  } else {
    cat(sprintf("  %-28s rel diff %s\n", lab, rd(a, b)))
  }
}
seeded <- function(f) function() {
  set.seed(1)
  f()
}
for (fam in c("cumulative", "sratio", "cratio", "acat", "hurdle_cumulative")) {
  Y <- if (fam == "hurdle_cumulative") "y0" else "y"
  F <- get(fam, envir = asNamespace("frmtmb"))
  ff <- function(s) stats::as.formula(gsub("Y", Y, s, fixed = TRUE))
  fl <- suppressWarnings(frm(bf(ff("Y ~ x + z")), family = F(), data = d))
  fn <- suppressWarnings(frm(bf(ff("Y ~ bx * x + bz * z"), bx + bz ~ 1,
                                nl = TRUE), family = F(), data = d))
  cat("\n==", fam, " logLik nl", format(as.numeric(logLik(fn)), digits = 10),
      "linear", format(as.numeric(logLik(fl)), digits = 10), "\n")
  try1("fitted(response)", function() fitted(fn), function() fitted(fl))
  try1("fitted(newdata)", function() fitted(fn, newdata = nd),
       function() fitted(fl, newdata = nd))
  try1("fitted(linear) estimate", function() fitted(fn, scale = "linear")[, 1],
       function() fitted(fl, scale = "linear")[, 1])
  try1("predict(seed)", seeded(function() predict(fn)),
       seeded(function() predict(fl)))
  try1("predict(newdata,seed)", seeded(function() predict(fn, newdata = nd)),
       seeded(function() predict(fl, newdata = nd)))
  try1("frm_linpred", function() frm_linpred(fn),
       function() frm_linpred(fl))
  try1("frm_linpred(newdata)", function() frm_linpred(fn, newdata = nd),
       function() frm_linpred(fl, newdata = nd))
  try1("frm_curve(newdata)",
       function() frm_curve(fn, newdata = nd)[["estimate"]],
       function() frm_curve(fl, newdata = nd)[["estimate"]])
  try1("simulate", seeded(function() simulate(fn, nsim = 1)),
       seeded(function() simulate(fl, nsim = 1)))
  try1("simulate(newdata)",
       seeded(function() simulate(fn, nsim = 1, newdata = nd)),
       seeded(function() simulate(fl, nsim = 1, newdata = nd)))
  try1("residuals(osa)", seeded(function() residuals(fn, type = "osa")),
       seeded(function() residuals(fl, type = "osa")))
  try1("log_lik", function() log_lik(fn), function() log_lik(fl))
  try1("loo", function() loo(fn)$estimates, function() loo(fl)$estimates)
  try1("vcov", function() {
    v <- vcov(fn)
    v[grep("^b[xz]_", rownames(v)), grep("^b[xz]_", colnames(v))]
  }, function() {
    v <- vcov(fl)
    v[c("x", "z"), c("x", "z")]
  })
  try1("confint(wald) slopes",
       function() confint(fn, parm = c("bx_(Intercept)", "bz_(Intercept)")),
       function() confint(fl, parm = c("x", "z")))
  try1("conditional_effects", function() {
    conditional_effects(fn, effects = "x")[[1]]$estimate__
  }, function() conditional_effects(fl, effects = "x")[[1]]$estimate__)
  try1("ce(categorical)", function() {
    conditional_effects(fn, effects = "x", categorical = TRUE)[[1]]$estimate__
  }, function() {
    conditional_effects(fl, effects = "x", categorical = TRUE)[[1]]$estimate__
  })
  try1("update(newdata)", function() logLik(update(fn, newdata = d[-1, ])),
       function() logLik(update(fl, newdata = d[-1, ])))
  try1("get_prior rows less 1", function() {
    nrow(get_prior(bf(ff("Y ~ bx * x + bz * z"), bx + bz ~ 1, nl = TRUE),
                   data = d, family = F())) - 1L
  }, function() {
    nrow(get_prior(bf(ff("Y ~ x + z")), data = d, family = F()))
  })
  try1("hypothesis", function() {
    hypothesis(fn, "bx_Intercept > 0")$hypothesis$Estimate
  }, function() hypothesis(fl, "x > 0")$hypothesis$Estimate)
}
cat("\n== na.exclude on an nl gaussian and an nl cumulative\n")
dn <- d
dn$x[c(2, 5)] <- NA
ga <- frm(bf(eta ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = gaussian(),
          data = dn, na.action = na.exclude)
print(dim(fitted(ga)))
print(length(frm_linpred(ga)))
oa <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = dn,
          na.action = na.exclude)
print(dim(fitted(oa)))
