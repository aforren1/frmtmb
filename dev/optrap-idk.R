# Lane optrap, item 3: test-id-kron.R's merged-dpar fit, printed to full
# precision, with what could differ between processes: the BLAS, the
# thread settings, and the first evaluations of the objective.
#   Rscript dev/optrap-idk.R <lib or ""> <label>
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "| lib", dirname(find.package("frmtmb")), "\n")
ex <- parse("tests/testthat/test-id-kron.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, globalenv())
  }
}
d <- kron_traits()
form <- bf(y1 ~ 1 + (1 | q | gr(id, cov = A)),
           sigma ~ 1 + (1 | q | gr(id, cov = A))) + gaussian()
u <- frm(form, data = d$wide, data2 = list(A = d$A),
         control = frmtmb_control(check_olre = "ignore"),
         dry_run = "objective")
# the objective and gradient at the start, bit for bit
f0 <- u$obj$fn(u$obj$par)
g0 <- u$obj$gr(u$obj$par)
cs <- character(0)
f <- withCallingHandlers(
  frm(form, data = d$wide, data2 = list(A = d$A),
      control = frmtmb_control(check_olre = "ignore")),
  warning = function(w) {
    cs <<- c(cs, paste("W", substr(conditionMessage(w), 1, 50)))
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cs <<- c(cs, paste("M", substr(conditionMessage(m), 1, 50)))
    invokeRestart("muffleMessage")
  })
o <- f$opt
cat(sprintf("%s | f0 %a | g0 %s | conv %s obj %.12g iter %s | %s\n",
            a[2], f0, paste(sprintf("%a", g0), collapse = ","),
            o$convergence, o$objective, o$iterations,
            paste(cs, collapse = "; ")))
cat(a[2], "| par", sprintf("%.10g", o$par), "\n")
cat(a[2], "| grad", sprintf("%.3g", f$obj$gr(o$par)), "\n")
cat(a[2], "| boundary stop", isTRUE(f$cache$se_boundary_stop), "kind",
    paste(f$cache$se_boundary_kind, collapse = ","), "\n")
lost <- tryCatch(frmtmb:::sdr_of(f)$se_lost, error = function(e) NULL)
cat(a[2], "| se lost:", paste(names(lost), lost, sep = ":", collapse = ","),
    "\n")
# where the two core types first disagree: the objective and gradient at
# fixed points off the start
set.seed(1)
for (k in 1:4) {
  p <- u$obj$par + stats::rnorm(length(u$obj$par), 0, 0.3)
  cat(a[2], "| probe", k, sprintf("%a", u$obj$fn(p)),
      paste(sprintf("%a", u$obj$gr(p)), collapse = ","), "\n")
}
# the same with the random effects held (no inner solve): the joint
# negative log-likelihood and its gradient
ff <- u$obj$env$f
pj <- u$obj$env$last.par
for (k in 1:2) {
  q <- pj + stats::rnorm(length(pj), 0, 0.3)
  cat(a[2], "| joint", k, sprintf("%a", ff(q, order = 0)),
      sprintf("%a", sum(ff(q, order = 1))), "\n")
}
