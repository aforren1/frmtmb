# Review of lane adhess, item 5: the fits whose pins moved, against a
# Richardson Jacobian of the AD gradient at the covariance's own point.
#   REV_ARM=lane|base Rscript dev/adhess-rev-pins.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
tt <- "C:/Users/adf44/source/r/frmtmb-wt-adhess/tests/testthat/"
ex <- parse(paste0(tt, "test-se-check.R"))
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, globalenv())
}
src <- readLines(paste0(tt, "helper-mo-flat.R"), warn = FALSE)
ok <- try(eval(parse(text = src), globalenv()), silent = TRUE)
show <- function(nm, fit) {
  at <- the_point(fit)
  V <- suppressWarnings(vcov(fit, full = TRUE))
  Hr <- rich_hess(fit$obj, at)
  Vr <- tryCatch(solve(Hr), error = function(e) NULL)
  if (is.null(Vr)) {
    st <- ns$obj_state_save(fit$obj); Hx <- fit$obj$he(at); ns$obj_state_restore(fit$obj, st)
    D <- sqrt(abs(diag(Hr)))
    cat(sprintf("%-10s Richardson Hessian singular; max |sym(AD H) - Rich| / DD %.2e; diag(V) %s
", nm, max(abs((Hx + t(Hx)) / 2 - Hr) / outer(D, D), na.rm = TRUE), paste(signif(diag(V), 3), collapse = "/")))
    return(invisible())
  }
  e <- cov_err(V, Vr)
  st <- ns$obj_state_save(fit$obj)
  eo <-cov_err(solve(ns$fd_hessian(at, fit$obj$gr)$H), Vr)
  ns$obj_state_restore(fit$obj, st)
  cat(sprintf("%-10s np %d re %d units %s | vcov vs Richardson se %.2e cor %.2e | optimHess route se %.2e\n",
              nm, length(at), length(fit$obj$env$random),
              !is.null(fit$par_units), e[["se"]], e[["cor"]], eo[["se"]]))
}
d <- se_mo_data(7)
d$ls <- d$ls / 16
show("healthy", frm(ls ~ mo(income), data = d))
if (exists("mo_flat_data")) {
  d2 <- mo_flat_data()
  short <- suppressWarnings(frm(ls ~ mo(inc) + age, data = d2,
    control = frmtmb_control(optCtrl = list(iter.max = 3, eval.max = 5),
                             restarts = 0)))
  show("short", short)
}
ss <- lme4::sleepstudy
show("sleep_ms", frm(Reaction ~ Days + (Days | Subject), data = ss))
ss$Reaction <- ss$Reaction / 32
show("sleep_32", frm(Reaction ~ Days + (Days | Subject), data = ss))
