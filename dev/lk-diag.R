suppressMessages(library(frmtmb))
options(digits = 8)
dpars_at <- function(fam, etas) {
  out <- list()
  for (nm in names(etas)) {
    out[[nm]] <- fam$links[[nm]]$linkinv(etas[[nm]])
    out[[paste0(".eta_", nm)]] <- etas[[nm]]
  }
  out
}
tape_nll <- function(fam, y, aterms, base_eta, sweep) {
  function(p) {
    etas <- base_eta
    etas[[sweep]] <- p[1] + 0 * base_eta[[sweep]]
    -sum(fam$lpdf(y, dpars_at(fam, etas), aterms))
  }
}
cat("== Beta: which link fails at |eta| = 30 ==\n")
for (lk in c("probit", "probit_approx", "cauchit", "softit")) {
  for (e0 in c(-30, 30, -20, 20)) {
    f <- tape_nll(Beta(lk), c(0.2, 0.5, 0.9), list(),
                  list(mu = 0, phi = log(5)), "mu")
    tp <- RTMB::MakeTape(f, e0)
    cat(sprintf("%-14s eta %4d  val %12.4g  grad %12.4g  logodds %10.4g\n",
                lk, e0, tp(e0), as.numeric(tp$jacobian(e0)),
                { lo <- frmtmb:::frmtmb_links[[lk]]$logit_eta
                  if (is.null(lo)) NA else lo(e0) }))
  }
}
cat("\n== negbinomial positive-mean links: gradient vs central difference ==\n")
for (lk in c("softplus", "squareplus", "sqrt")) {
  for (e0 in c(-30, 30)) {
    f <- tape_nll(negbinomial(lk), c(0, 1, 3), list(),
                  list(mu = 0, shape = log(2)), "mu")
    tp <- RTMB::MakeTape(f, e0)
    g <- as.numeric(tp$jacobian(e0))
    for (h in c(1e-4 * 30, 1e-5, 1e-6)) {
      fd <- (f(e0 + h) - f(e0 - h)) / (2 * h)
      cat(sprintf("%-11s eta %4d h %8.1e  AD %14.8g  fd %14.8g  rel %9.3g\n",
                  lk, e0, h, g, fd, abs(g - fd) / max(1, abs(fd))))
    }
  }
}
