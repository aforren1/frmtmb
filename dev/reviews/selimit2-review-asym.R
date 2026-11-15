# Review of lane selimit2, user question (b): asym_laplace() without random
# effects. Are its point estimates quantile regression's (quantreg::rq),
# and what do rq's standard errors (nid, ker, iid) read beside frmtmb's?
.libPaths(c("C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(quantreg))
for (s in 1:4) {
  set.seed(s)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rexp(n) * (1 + 0.3 * abs(d$x))
  for (tau in c(0.3, 0.5)) {
    ws <- character(0)
    f <- withCallingHandlers(
      frm(bf(y ~ x, quantile = tau), data = d, family = asym_laplace()),
      warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
      }, message = function(m) invokeRestart("muffleMessage"))
    b <- f$estimates$beta
    r <- quantreg::rq(y ~ x, tau = tau, data = d)
    chk <- function(bb) sum((d$y - bb[1] - bb[2] * d$x) *
                              (tau - (d$y - bb[1] - bb[2] * d$x < 0)))
    se <- function(m) tryCatch(
      summary(r, se = m)$coefficients[, 2], error = function(e) c(NA, NA))
    cat(sprintf(paste0("seed %d tau %.1f | code %d | frm %.5f %.5f | rq %.5f %.5f",
                       " | check loss frm-rq %.3g | nwarn %d\n"),
                s, tau, f$opt$convergence, b[1], b[2], coef(r)[1], coef(r)[2],
                chk(b) - chk(coef(r)), length(ws)))
    cat("   rq se nid", format(se("nid"), digits = 4), "| ker",
        format(se("ker"), digits = 4), "| iid", format(se("iid"), digits = 4),
        "| frm vcov", format(suppressWarnings(sqrt(diag(vcov(f)))[1:2]),
                             digits = 4), "\n")
  }
}
