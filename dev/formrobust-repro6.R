# Item 6: 0 + intercept, drop_unused_levels, autocorrelation term
# objects, autocor(). Seed 61. FORMROBUST_LIB="" for the before arm.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(lab, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  cat(sprintf("%-34s %s%s\n", lab, paste(format(r, digits = 10),
                                         collapse = " "),
              if (length(w)) paste0("  [warning: ", paste(w, collapse = "; "),
                                    "]") else ""))
}
set.seed(61)
d <- data.frame(x = rnorm(40), g = gl(4, 10))
d$y <- 1 + 0.5 * d$x + rnorm(40)
tr("0 + intercept + x coefs",
   rownames(fixef(frm(bf(y ~ 0 + intercept + x), data = d))))
tr("same logLik as y ~ x",
   identical(logLik(frm(bf(y ~ 0 + intercept + x), data = suppressWarnings(d)))[1],
             logLik(frm(bf(y ~ x), data = d))[1]))
f1 <- suppressWarnings(frm(bf(y ~ 0 + intercept + x), data = d))
tr("predict(newdata) without column", dim(fitted(f1, newdata = d[1:3, "x",
                                                            drop = FALSE])))
tr("intercept column not ones", {
  d2 <- d; d2$intercept <- 2
  frm(bf(y ~ 0 + intercept + x), data = d2)
})
dd <- data.frame(y = rnorm(10), x = factor(c("a", "b"),
                                           levels = c("a", "b", "c")))
tr("drop_unused_levels = FALSE X",
   colnames(frm(bf(y ~ x), data = dd, drop_unused_levels = FALSE,
                dry_run = "frame")$linpreds[[1]]$X))
tr("drop_unused_levels = TRUE X",
   colnames(frm(bf(y ~ x), data = dd, dry_run = "frame")$linpreds[[1]]$X))
tr("ar() object", class(ar(week, subj)))
tr("bf + arma(x)", bf(y ~ 1) + arma(x))
tr("bf + acformula", deparse1((bf(y ~ x) + acformula(~ ar(t, g)))$formula))
tr("bf(autocor =)", deparse1(bf(y ~ x, autocor = ~ arma(t, g))$formula))
dg <- expand.grid(t = 1:5, g = factor(1:12))
dg$y <- rnorm(60)
fa <- frm(bf(y ~ 1 + ar(t, g)), data = dg)
tr("autocor(fit)", is.null(autocor(fa)))
tr("autocor(fit, resp = 'zz')", autocor(fa, resp = "zz"))
tr("acformula fit identical",
   identical(logLik(frm(bf(y ~ 1) + acformula(~ ar(t, g)), data = dg))[1],
             logLik(fa)[1]))
