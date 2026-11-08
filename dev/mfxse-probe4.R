# Probe: the exact gp() row at an unseen position, seed 8150: its kriging
# variance, the order of get_vcov() and get_predict() calls, and the row
# check's verdict.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
ns <- asNamespace("frmtmb")
set.seed(8150)
m <- 90
dg <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
dg$y <- sin(dg$x) + 0.3 * dg$z + rnorm(m, sd = 0.3)
fit <- quiet(frm(y ~ z + gp(x), data = dg))
ndg <- data.frame(x = c(dg$x[5], 3.33, 7, 12), z = 3)
cat("theta:", fit$estimates$theta, "\n")
for (lp in fit$frame$linpreds) {
  ed <- ns$lp_eta_design(fit, lp, ndg, TRUE, FALSE)
  for (sp in ed$sm_parts) cat("extra_var:", sp$extra_var, "\n")
}
l <- frm_linpred(fit, newdata = ndg, se.fit = TRUE)
cat("linpred se:", l$se.fit, "\n")
ord <- character(0)
for (nm in c("get_vcov", "get_predict")) {
  local({
    f0 <- get(paste0(nm, ".frmtmb_fit"), envir = ns)
    nm0 <- nm
    registerS3method(nm, "frmtmb_fit", function(model, ...) {
      ord <<- c(ord, nm0)
      f0(model, ...)
    }, envir = asNamespace("marginaleffects"))
  })
}
p <- withCallingHandlers(predictions(fit, newdata = ndg),
                         warning = function(w) {
                           cat("WARN:", conditionMessage(w), "\n")
                           invokeRestart("muffleWarning")
                         })
cat("order:", paste(head(rle(ord)$values, 6), head(rle(ord)$lengths, 6)),
    "\n")
cat("mfx se:", p$std.error, "\n")
ft <- fitted(fit, newdata = ndg)
cat("fitted se:", ft[, "Est.Error"], "\n")
