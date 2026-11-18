# rel084-kt, merge step 8 (drawdir x gpselimit): on gpselimit's 13 white
# rows (gp() fits at their basis's white limit; rows 1-5, 9, 10, 14, 21, 22,
# 25-27 of the sweep, data from rv_sweep_row() of the review lib, the same
# generator as gs_sweep_row()), run predict() on the merged build and record
# the cached draw check (fire, excess, rel, n, dir, dir_par, dir_n,
# warn_only) and every condition class raised by the fit and by predict().
# Usage: Rscript dev/rel084-kt/probe-step8-drawdir-gpselimit.R
source("C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/gpselimit-review-lib.R")
rv_libs("lane")
ns <- asNamespace("frmtmb")
cls_of <- function(expr) {
  said <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      said <<- c(said, paste0("w:", paste(class(w), collapse = "/")))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      said <<- c(said, paste0("m:", paste(class(m), collapse = "/")))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      said <<- c(said, paste0("E:", paste(class(e), collapse = "/"), ": ",
                              conditionMessage(e)))
      NULL
    })
  list(v = v, said = said)
}
fmt <- function(x) if (is.null(x)) "NULL" else paste(format(x, digits = 6), collapse = ",")
for (row in c(1:5, 9, 10, 14, 21, 22, 25:27)) {
  s <- rv_sweep_row(row)
  f <- cls_of(frm(bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian()))
  fit <- f$v
  if (is.null(fit)) {
    cat(sprintf("ROW S%02d fit ERROR | %s\n", row, paste(f$said, collapse = " ; ")))
    next
  }
  wh <- vapply(rv_hsgp(fit), function(bk) {
    ns$gp_basis_white(bk, fit$estimates$theta[bk$theta_idx])
  }, NA)
  p1 <- cls_of({ set.seed(5); predict(fit, ndraws = 200, summary = FALSE) })
  chk <- fit$cache[["draw_check"]]
  p2 <- cls_of(predict(fit))
  cat(sprintf(paste0("ROW S%02d ll %.6f white %s | fit conds [%s]\n",
                     "    predict(ndraws=200) ok %s conds [%s]\n",
                     "    predict() ok %s conds [%s]\n",
                     "    draw_check: fire %s warn_only %s dir %s dir_par %s ",
                     "dir_n %s | excess %s rel %s n %s\n"),
              row, as.numeric(logLik(fit)),
              paste(as.integer(wh), collapse = ""),
              paste(f$said, collapse = " ; "),
              !is.null(p1$v), paste(p1$said, collapse = " ; "),
              !is.null(p2$v), paste(p2$said, collapse = " ; "),
              fmt(chk$fire), fmt(chk$warn_only), fmt(chk$dir),
              fmt(chk$dir_par), fmt(chk$dir_n), fmt(chk$excess),
              fmt(chk$rel), fmt(chk$n)))
  cat("    opt$par names:", paste(names(fit$opt$par), collapse = " "), "\n")
}
cat("DONE\n")
