# Reviewer's check of lane bin2port's proposed fix for the cumulative
# tied-threshold NaN (dev/bin2port-findings.md, log_lik:462): patch
# ord_cumulative_logpmf() IN THIS SESSION ONLY (assignInNamespace; no
# library on disk changes) and compare with the unpatched function.
#   Rscript dev/bin2port-rev-cumfix.R      seed 20261014
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("ord_cumulative_logpmf", ns)
src <- deparse(orig)
pat <- "a <- q(disc * (tau[ym] - eta))"
stopifnot(sum(grepl(pat, src, fixed = TRUE)) == 1L)
src2 <- sub(pat, "a <- q(disc * (tau[ym] - eta)) + (i1 + iK)", src,
            fixed = TRUE)
patched <- eval(parse(text = src2))
environment(patched) <- ns
use <- function(f) {
  unlockBinding("ord_cumulative_logpmf", ns)
  assign("ord_cumulative_logpmf", f, envir = ns)
  lockBinding("ord_cumulative_logpmf", ns)
}
exact <- function(y, mu, tau, link = "logit") {
  Fc <- switch(link, logit = stats::plogis, probit = stats::pnorm,
               cloglog = function(x) 1 - exp(-exp(x)))
  vapply(seq_along(mu), function(s) {
    F <- c(0, Fc(tau - mu[s]), 1)
    log(F[y[s] + 1L] - F[y[s]])
  }, 0)
}
dens <- function(fam, y, mu, tau, hu = NULL) {
  lp <- fam[["ord_lpdf_make"]](identity)
  dp <- list(mu = mu, disc = 1)
  if (!is.null(hu)) dp$hu <- hu
  as.numeric(lp(y, dp, list(), list(tau_raw = tau)))
}
set.seed(20261014)
cat("== 1. tied thresholds, patched vs exact (logit, probit, cloglog)\n")
use(patched)
mu <- c(-0.7, 0, 1.3)
for (link in c("logit", "probit", "cloglog")) {
  fam <- frmtmb::cumulative(link = link)
  for (tau in list(c(0, 0, 0), c(-1, -1, 1), c(-1, 1, 1))) {
    for (y in 1:4) {
      a <- dens(fam, rep(y, 3), mu, tau)
      e <- exact(rep(y, 3), mu, tau, link)
      ok <- identical(is.finite(a), is.finite(e)) &&
        isTRUE(all.equal(a[is.finite(e)], e[is.finite(e)]))
      if (!ok || y %in% c(1, 4)) {
        cat(sprintf("  %-7s tau (%s) y=%d  patched %s  exact %s  %s\n", link,
                    paste(tau, collapse = ","), y,
                    paste(format(a, digits = 6), collapse = " "),
                    paste(format(e, digits = 6), collapse = " "),
                    if (ok) "ok" else "MISMATCH"))
      }
    }
  }
}
cat("== 2. hurdle_cumulative at tau = (0, 0, 1), y = 0, 1, 4\n")
hfam <- frmtmb::hurdle_cumulative()
for (f in list(orig = orig, patched = patched)) {
  use(f)
  cat("  ", paste(format(dens(hfam, c(0, 1, 4), mu, c(0, 0, 1),
                               hu = rep(0.2, 3)), digits = 6),
                  collapse = " "), "\n")
}
cat("  exact y=1,4 parts: log(0.8) +",
    paste(format(exact(c(1, 4), mu[2:3], c(0, 0, 1)), digits = 6),
          collapse = " "), "; y=0: log(0.2) =", format(log(0.2)), "\n")
cat("== 3. distinct thresholds, 2000 random rows per link: patched",
    "identical to unpatched on doubles?\n")
for (link in c("logit", "probit", "cloglog", "cauchit")) {
  fam <- frmtmb::cumulative(link = link)
  y <- sample(1:5, 2000, TRUE)
  m <- rnorm(2000, 0, 3)
  tau <- sort(rnorm(4, 0, 2))
  use(orig); a0 <- dens(fam, y, m, tau)
  use(patched); a1 <- dens(fam, y, m, tau)
  cat(sprintf("  %-8s identical %s  (NaN in unpatched: %d)\n", link,
              identical(a0, a1), sum(is.nan(a0))))
}
cat("== 4. a fit, unpatched and patched: same tape?\n")
n <- 400
x <- rnorm(n)
z <- 0.8 * x + stats::rlogis(n)
d <- data.frame(x = x, y = factor(cut(z, c(-Inf, -1, 0, 1.2, Inf),
                                      labels = FALSE), ordered = TRUE))
fits <- list()
for (nm in c("orig", "patched")) {
  use(if (nm == "orig") orig else patched)
  fits[[nm]] <- frm(bf(y ~ x), family = cumulative(), data = d)
}
p0 <- fits$orig$opt$par %||% fits$orig$tmb$par
cat("  logLik identical:", identical(as.numeric(logLik(fits$orig)),
                                      as.numeric(logLik(fits$patched))),
    " coef identical:", identical(coef(fits$orig), coef(fits$patched)),
    " vcov identical:", identical(vcov(fits$orig), vcov(fits$patched)),
    "\n")
obj <- function(fit) fit$tmb %||% fit$obj
o0 <- obj(fits$orig); o1 <- obj(fits$patched)
if (!is.null(o0) && is.function(o0$fn)) {
  pp <- o0$par + rnorm(length(o0$par), 0, 0.1)
  cat("  objective at a perturbed point identical:",
      identical(o0$fn(pp), o1$fn(pp)), " gradient identical:",
      identical(o0$gr(pp), o1$gr(pp)), "\n")
} else {
  cat("  (no objective handle found on the fit; names:",
      paste(names(fits$orig), collapse = " "), ")\n")
}
use(orig)
