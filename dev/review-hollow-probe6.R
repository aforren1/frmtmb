# Reviewer probe 6 (fix round): try to break the fixed
# draws_newdata_is_fitted() again
.libPaths(c("C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb.sample"), "\n")
isf <- frmtmb.sample:::draws_newdata_is_fitted
mk <- function(fit) structure(list(stanfit = NULL, draws = matrix(0, 1, 1),
                                   fit = fit), class = "frmtmb_draws")
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(6)
n <- 40
dd <- data.frame(x = rnorm(n), z = rnorm(n), t2 = runif(n),
                 w = runif(n, 0.5, 2), s = rep(c(TRUE, TRUE, FALSE, TRUE), 10),
                 sdx = runif(n, 0.1, 0.3), age = runif(n, 1, 3),
                 cc = rep(c(0, 0, 0, 2), 10), u = runif(n))
dd$y <- 1 + 0.5 * dd$x + rnorm(n)
dd$y2 <- dd$y + 1
dd$hi <- dd$y + 2
y2e <- dd$y2  # a bound found in the formula environment, not the data

bump <- function(d, v, i = 4) {
  if (is.logical(d[[v]])) d[[v]][i] <- !d[[v]][i] else d[[v]][i] <- d[[v]][i] + 0.5
  d
}
row <- function(lbl, fit, changed) {
  if (inherits(fit, "error")) {
    cat(sprintf("%-26s FIT ERROR: %s\n", lbl, conditionMessage(fit)))
    return(invisible())
  }
  ds <- mk(fit)
  ch <- vapply(changed, function(v) isf(ds, bump(dd, v)), NA)
  cat(sprintf("%-26s frame=[%s] raw=[%s] | data %s, fit$data %s | changed: %s\n",
              lbl, paste(names(fit$data), collapse = ","),
              paste(names(fit$frame$raw_vars), collapse = ","),
              isf(ds, dd), isf(ds, fit$data),
              paste(sprintf("%s=%s", changed, ifelse(ch, "ACCEPTED", "refused")),
                    collapse = " ")))
}
tryfit <- function(expr) tryCatch(q(expr), error = function(e) e)

row("nl body only (t2)",
    tryfit(frm(bf(y ~ a * exp(b * t2), a ~ 1, b ~ 1, nl = TRUE), data = dd)),
    c("t2", "u"))
row("dpar sigma ~ z",
    tryfit(frm(bf(y ~ x, sigma ~ z), data = dd)), c("z", "x", "u"))
row("offset only (age)",
    tryfit(frm(y ~ x + offset(log(age)), data = dd)), c("age", "u"))
row("weights(w)",
    tryfit(frm(y | weights(w) ~ x, data = dd)), c("w", "u"))
row("mvbf subset(s)",
    tryfit(frm(bf(y | subset(s) ~ x) + bf(z ~ x) + set_rescor(FALSE),
               data = dd)), c("s", "x", "u"))
row("me(x, sdx)",
    tryfit(frm(y ~ me(x, sdx), data = dd)), c("sdx", "x", "u"))
row("cens(cc, y2) column",
    tryfit(frm(y | cens(cc, y2) ~ x, data = dd)), c("y2", "cc", "u"))
row("cens(cc, 5) number",
    tryfit(frm(y | cens(cc, 5) ~ x, data = dd)), c("y2", "cc", "u"))
row("trunc(ub = hi) column",
    tryfit(frm(y | trunc(ub = hi) ~ x, data = dd)), c("hi", "u"))
row("trunc(ub = 10) number",
    tryfit(frm(y | trunc(ub = 10) ~ x, data = dd)), c("hi", "u"))
fe <- tryfit(frm(y | cens(cc, y2e) ~ x, data = dd))
row("cens(cc, y2e) env bound", fe, c("y2", "u"))
if (!inherits(fe, "error")) {
  dse <- mk(fe)
  d_ok <- cbind(dd, y2e = y2e)
  d_mv <- d_ok; d_mv$y2e[4] <- d_mv$y2e[4] + 3
  cat("  env bound given as a column: equal", isf(dse, d_ok),
      " moved", isf(dse, d_mv), "\n")
}
# the identical() shortcut on a cens model: the frame object, and the
# frame with a moved bound beside it
fc <- q(frm(y | cens(cc, y2) ~ x, data = dd))
dsc <- mk(fc)
fr <- fc$data
fr2 <- fr; fr2$y2 <- dd$y2; fr2$y2[4] <- fr2$y2[4] + 3
fr3 <- fr; attr(fr3, "terms") <- NULL
cat("cens shortcut: frame object", isf(dsc, fr), " frame + moved y2",
    isf(dsc, fr2), " frame without terms attr", isf(dsc, fr3), "\n")
