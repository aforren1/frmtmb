# Does the frame's factor order still meet TMB's compressed-sweep bug if b is not put
# first? dev/ncgp-grcheck.R and REML at low noise on five seeds, with gp_tape_order()
# the identity. dev/ncgp-remlarm2.R reads it on a hand tape over a triangular factor.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
utils::assignInNamespace("gp_tape_order", function(frame, tpl) tpl, ns = "frmtmb")
ex <- parse("dev/ncgp-truth.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]])[1] %in% c("marginal", "closed_max", "REL")) eval(e)
}
for (seed in 1:5) {
  set.seed(seed)
  x <- round(stats::runif(60, 0, 6), 3)
  y <- sin(2 * x) + 1e-3 * stats::rnorm(60)
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ gp(x)), data = data.frame(x = x, y = y), family = gaussian(), REML = TRUE),
    warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 50)); invokeRestart("muffleWarning") })
  p <- fit$opt$par
  g <- fit$obj$gr(p)
  fd <- vapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, 1e-5)
    (fit$obj$fn(p + e) - fit$obj$fn(p - e)) / 2e-5
  }, 0)
  cat(sprintf("NCGPNO seed %d conv %d iter %d sigma %.4g | max|gr| %.2e | max|gr - fd| %.2e | %s\n", seed,
              fit$opt$convergence, fit$opt$iterations, exp(fit$estimates$betad[1]), max(abs(g)),
              max(abs(g - fd)), paste(unique(w), collapse = "; ")))
}
