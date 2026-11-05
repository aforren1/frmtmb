# Short check, item 7: the hook failure of test-nl-family.R on the merge. Seed 11001, the
# test's own data and formula. Rscript dev/review-predread-r1-hook.R <core lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(11001)
n <- 300
dd <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
dd$yo <- ordered(cut(1.5 * exp(0.6 * dd$x) + 0.7 * dd$z + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
fm <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) + nlf(disc ~ exp(d1 * d2 * z)) + lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) + lf(d2 ~ 1)
ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
w <- character()
fit <- withCallingHandlers(frm(fm, family = cumulative(), data = dd, control = ctl),
  warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage"))
cat("warnings:", length(w), "\n")
cat(substr(w, 1, 160), sep = "\n")
lp <- fit$frame$linpreds[["yo.disc"]]
cat("disc lp X:", class(lp$X), " dim", dim(lp$X), " nl_body", !is.null(lp$nl_body), " constant", !is.null(lp$constant), " n_param_cols", lp$n_param_cols, "\n")
r <- withCallingHandlers(tryCatch(frmtmb:::ord_fit_check(fit, "yo"), error = function(e) paste("ERROR:", conditionMessage(e))),
  error = function(e) {
    cs <- sys.calls()
    cat("CALLS:", paste(tail(vapply(cs, function(cl) substr(deparse(cl)[1], 1, 70), ""), 8), collapse = "\n  "), "\n")
  })
print(r)
print(body(frmtmb:::ord_fit_check)[[5]])
