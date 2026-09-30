# Reviewer: the new exports, autocor(), bf(autocor =), acformula(),
# 0 + intercept, drop_unused_levels (claim 6). Data seed 61.
# REVLIB="" for the base arm.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
msgs <- character(0)
withCallingHandlers(library(frmtmb), message = function(m) {
  msgs <<- c(msgs, conditionMessage(m))
}, packageStartupMessage = function(m) {
  msgs <<- c(msgs, conditionMessage(m))
})
cat("frmtmb from", find.package("frmtmb"), "\n")
cat("brms namespace loaded after library(frmtmb):",
    isNamespaceLoaded("brms"), "\n")
cat("conflicts with stats:", intersect(ls("package:frmtmb"),
                                       ls("package:stats")), "\n")
tr <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) cat(sprintf("[%s] ERROR: %s\n", label, r))
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, unique(w)))
  invisible(if (inherits(r, "err")) NULL else r)
}
same <- function(label, a, b) cat(sprintf("[%s] identical: %s\n", label,
                                          identical(a, b)))
cat("\n== ar masking ==\n")
r <- tr("ar(lh)", ar(lh))
cat("class of ar(lh):", class(r), "\n")
r2 <- stats::ar(lh)
cat("stats::ar(lh) order:", r2$order, "\n")
cat("ar found on search path from:", environmentName(environment(ar)), "\n")
cat("\n== autocor() without brms loaded ==\n")
set.seed(61)
d <- expand.grid(t = 1:6, g = factor(1:15))
d$x <- rnorm(nrow(d)); d$y <- rnorm(nrow(d)); d$y2 <- rnorm(nrow(d))
fit <- frm(bf(y ~ x + ar(t, g)), data = d)
cat("brms loaded before autocor():", isNamespaceLoaded("brms"), "\n")
r <- tr("autocor(fit)", autocor(fit))
cat("value is NULL:", is.null(r), "\n")
tr("autocor(fit, resp = 'zz')", autocor(fit, resp = "zz"))
tr("autocor(fit, foo = 1)", autocor(fit, foo = 1))
tr("autocor on a non-fit", autocor(1))
loadNamespace("brms")
r <- tr("autocor(fit) after loadNamespace(brms)", autocor(fit))
r <- tr("brms::autocor(fit)", brms::autocor(fit))
cat("brms::autocor(fit) is NULL:", is.null(r), "\n")
suppressMessages(library(brms, warn.conflicts = FALSE))
r <- tr("autocor(fit) with brms attached", autocor(fit))
cat("  NULL:", is.null(r), "; ar now from:",
    environmentName(environment(ar)), "\n")
tr("frmtmb::bf(y ~ x) + brms ar()", frmtmb::bf(y ~ x) + ar(t, g))
tr("frmtmb::bf(y ~ x) + frmtmb::ar()", frmtmb::bf(y ~ x) + frmtmb::ar(t, g))
detach("package:brms")
bf <- frmtmb::bf

cat("\n== bf(autocor =), acformula() ==\n")
f1 <- frm(bf(y ~ x, autocor = ~ ar(t, g)), data = d)
same("bf(autocor =) logLik", logLik(f1), logLik(fit))
f2 <- frm(bf(y ~ x) + acformula(~ ar(t, g)), data = d)
same("acformula logLik", logLik(f2), logLik(fit))
tr("acformula(~ x)", acformula(~ x))
tr("acformula(~ ar(t, g) + x)", acformula(~ ar(t, g) + x))
tr("mv + acformula without resp", bf(y ~ x) + bf(y2 ~ x) +
     acformula(~ ar(t, g)))
fm <- tr("mv + acformula(resp = 'y')",
         frm(bf(y ~ x) + bf(y2 ~ x) + acformula(~ ar(t, g), resp = "y"),
             data = d))
fm2 <- tr("mv with ar in y's formula", frm(bf(y ~ x + ar(t, g)) +
                                              bf(y2 ~ x), data = d))
if (!is.null(fm) && !is.null(fm2)) same("mv acformula logLik",
                                        logLik(fm), logLik(fm2))
tr("bf(mvbind(y, y2) ~ x, autocor = ~ ar(t, g))", {
  fmb <- frm(bf(mvbind(y, y2) ~ x, autocor = ~ ar(t, g)), data = d)
  fmc <- frm(bf(mvbind(y, y2) ~ x + ar(t, g)), data = d)
  same("mvbind autocor logLik", logLik(fmb), logLik(fmc)) })
tr("bf(y ~ x, autocor = ~ ar(t, g)) twice via bf update",
   print(bf(bf(y ~ x), autocor = ~ ar(t, g))$formula))
tr("bf(y ~ x + ar(t, g), autocor = ~ ma(t, g))",
   frm(bf(y ~ x + ar(t, g), autocor = ~ ma(t, g)), data = d))
tr("brms: bf(y ~ x + ar(t, g), autocor = ~ ma(t, g))",
   print(brms::bf(y ~ x + ar(t, g), autocor = ~ ma(t, g))))
tr("update(fit, autocor = ~ ma(t, g))", update(fit, autocor = ~ ma(t, g)))

cat("\n== 0 + intercept ==\n")
fi <- tr("0 + intercept + x", frm(bf(y ~ 0 + intercept + x), data = d))
f0 <- frm(bf(y ~ x), data = d)
if (!is.null(fi)) {
  same("0 + intercept logLik vs y ~ x", logLik(fi), logLik(f0))
  cat("  coef names:", rownames(fixef(fi)), "\n")
  tr("predict newdata 0 + intercept", predict(fi, newdata = d[1:3, ]))
  tr("CE 0 + intercept", cat("  CE displays:",
                             names(conditional_effects(fi)), "\n"))
}
dd <- d; dd$intercept <- 2
tr("0 + intercept, data intercept = 2", frm(bf(y ~ 0 + intercept + x),
                                            data = dd))
dd$intercept <- rnorm(nrow(dd))
tr("y ~ x + intercept (covariate, with intercept)",
   print(rownames(fixef(frm(bf(y ~ x + intercept), data = dd)))))
tr("brms standata y ~ x + intercept (covariate)",
   print(colnames(brms::standata(y ~ x + intercept, data = dd)$X)))

cat("\n== drop_unused_levels ==\n")
d$xf <- factor(sample(c("a", "b"), nrow(d), TRUE), levels = c("a", "b", "c"))
d$yy <- rnorm(nrow(d), as.numeric(d$xf))
fT <- frm(bf(yy ~ xf), data = d)
fF <- tr("drop_unused_levels = FALSE", frm(bf(yy ~ xf), data = d,
                                           drop_unused_levels = FALSE))
if (!is.null(fF)) {
  cat("  coefs TRUE:", rownames(fixef(fT)), " FALSE:", rownames(fixef(fF)),
      "\n")
  same("logLik TRUE vs FALSE", logLik(fT), logLik(fF))
  nd <- data.frame(xf = factor(c("a", "c"), levels = c("a", "b", "c")))
  tr("fitted(FALSE) on level c", print(fitted(fF, newdata = nd)))
  tr("fitted(TRUE) on level c", print(fitted(fT, newdata = nd)))
  s <- brms::standata(yy ~ xf, data = d, drop_unused_levels = FALSE)
  cat("  brms X columns:", colnames(s$X), " xfc all zero:",
      all(s$X[, "xfc"] == 0), "\n")
  tr("emmeans(FALSE)", print(emmeans::emmeans(fF, ~ xf)))
  tr("CE(FALSE)", print(conditional_effects(fF)[[1]][, c("xf",
                                                         "estimate__")]))
  tr("influence(FALSE)", invisible(influence(fF)))
}
tr("drop_unused_levels = NA", frm(bf(yy ~ xf), data = d,
                                  drop_unused_levels = NA))
