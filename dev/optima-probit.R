# Lane optima, item 2: the probit's log-odds, q(eta) = log(Phi(eta)) -
# log(Phi(-eta)), against Rmpfr at 256 bits, value and derivative,
# through RTMB's tape as a fit reads them. Two forms:
#   old: log(RTMB::pnorm(eta)) - log(RTMB::pnorm(-eta))   (0.68.1)
#   new: RTMB::pnorm(eta, log.p = TRUE) - RTMB::pnorm(-eta, log.p = TRUE)
# and, when the installed frmtmb is given, the link's own logit_eta.
#   Rscript dev/optima-probit.R base|lane
args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
old <- function(e) log(RTMB::pnorm(e)) - log(RTMB::pnorm(-e))
new <- function(e) {
  RTMB::pnorm(e, log.p = TRUE) - RTMB::pnorm(-e, log.p = TRUE)
}
link <- frmtmb:::frmtmb_links$probit$logit_eta
forms <- list(old = old, new = new, installed = link)

# exact reference: q and dq/deta = phi/Phi(eta) + phi/Phi(-eta)
prec <- 256
ref <- function(eta) {
  x <- mpfr(eta, prec)
  P <- pnorm(x)
  M <- pnorm(-x)
  ph <- dnorm(x)
  list(q = as.numeric(log(P) - log(M)),
       dq = as.numeric(ph / P + ph / M))
}

relerr <- function(a, b) {
  ifelse(is.finite(a), abs(a - b) / pmax(abs(b), .Machine$double.xmin), Inf)
}
report <- function(eta, lab) {
  R <- ref(eta)
  cat("==", lab, ":", length(eta), "points on [", min(eta), ",", max(eta),
      "]\n")
  for (nm in names(forms)) {
    f <- forms[[nm]]
    v <- RTMB::MakeTape(function(e) f(e), eta)(eta)
    g <- diag(RTMB::MakeTape(function(e) f(e), eta)$jacobian(eta))
    ev <- relerr(v, R$q)
    eg <- relerr(g, R$dq)
    cat(sprintf(paste0("  %-9s value: non-finite %4d, max rel err %.3g,",
                       " median %.3g | derivative: non-finite %4d, max rel",
                       " err %.3g, median %.3g\n"),
                nm, sum(!is.finite(v)), max(ev[is.finite(ev)], -Inf),
                stats::median(ev[is.finite(ev)]), sum(!is.finite(g)),
                max(eg[is.finite(eg)], -Inf),
                stats::median(eg[is.finite(eg)])))
  }
}
report(seq(-38, 38, length.out = 2001), "the old form's finite range")
report(seq(-200, 200, length.out = 4001), "wide")
report(c(-1e4, -1e3, -500, -100, -60, -40, -38.5, 38.5, 40, 60, 100, 500,
         1e3, 1e4), "far points")
eta <- seq(-38, 38, length.out = 2001)
vo <- RTMB::MakeTape(old, eta)(eta)
vn <- RTMB::MakeTape(new, eta)(eta)
cat("old vs new on [-38, 38]:", sum(vo != vn), "of", length(eta),
    "values differ, max relative", format(max(relerr(vo, vn)), digits = 3),
    "\n")
