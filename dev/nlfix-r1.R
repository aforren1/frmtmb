# R1: the two designs of the short check, with its seed, plus mvbf and the disc
# body. Prints the value from the default start against every coefficient moved.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
source("dev/review-nlfix-common.R", local = TRUE, echo = FALSE)
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2", .libPaths()))
show <- function(lab, expr) {
  msg <- character(0); warn <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(r, "error")) { cat(lab, ": ERROR", conditionMessage(r), "\n"); return(invisible(r)) }
  cat(sprintf("%-22s ll %.4f code %d\n", lab, as.numeric(logLik(r)), r$opt$convergence))
  for (m in msg) cat("   MSG ", substr(m, 1, 300), "\n")
  for (w in warn) cat("   WARN", substr(w, 1, 300), "\n")
  invisible(r)
}
set.seed(11001)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3), g = factor(rep(1:15, each = 20)), h = factor(rep(1:30, 10)))
u <- rnorm(15, 0, 0.5)
lat <- 1.5 * exp(0.6 * d$x) + 0.7 * d$z
d$yo <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
d$yn <- 1 + 0.7 * d$x + 0.6 * d$z + rnorm(n, 0, exp(-0.5 + 0.4 * d$z))
fm <- bf(yn ~ 1) + nlf(mu ~ a + b1 * c1 * x) + nlf(sigma ~ exp(s0 + s1 * s2 * z)) + lf(a ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) + lf(s0 ~ 1) + lf(s1 ~ 1) + lf(s2 ~ 1)
r0 <- show("gaussian default", frm(fm, data = d))
cat("  sigma range", diff(range(frm_linpred(r0, dpar = "sigma"))), "\n")
print(round(unlist(r0$estimates[c("beta", "betad")]), 3))
show("gaussian linear", frm(bf(yn ~ x, sigma ~ z), data = d))
set.seed(11001); invisible(runif(1))
d$ymix <- ifelse(runif(n) < 0.4, -2 + 0.8 * d$x, 2 - 0.6 * d$x) + rnorm(n, 0, 0.5)
fm2 <- bf(ymix ~ 1) + nlf(mu1 ~ a1 + b1 * c1 * x) + nlf(mu2 ~ a2 + b2 * c2 * x) + lf(a1 ~ 1) + lf(b1 ~ 1) + lf(c1 ~ 1) + lf(a2 ~ 1) + lf(b2 ~ 1) + lf(c2 ~ 1)
fam <- mixture(gaussian(), gaussian())
r2 <- show("mixture default", frm(fm2, family = fam, data = d))
cat("  mu1 range", diff(range(frm_linpred(r2, dpar = "mu1"))), " mu2 range", diff(range(frm_linpred(r2, dpar = "mu2"))), "\n")
show("mixture linear", frm(bf(ymix ~ 1, mu1 ~ x, mu2 ~ x), family = fam, data = d))
set.seed(11001)
d$y1 <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
d$y2 <- ordered(cut(1.2 * exp(-0.5 * d$z) + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
f1 <- bf(y1 ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
f2 <- bf(y2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE)
s1 <- show("ord y1", frm(f1, family = cumulative(), data = d))
s2 <- show("ord y2", frm(f2, family = cumulative(), data = d))
cat("  sum", as.numeric(logLik(s1)) + as.numeric(logLik(s2)), "\n")
show("mvbf default", frm(mvbf(f1, f2, rescor = FALSE), family = cumulative(), data = d))
fm4 <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) + nlf(disc ~ exp(d1 * d2 * z)) + lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) + lf(d2 ~ 1)
show("disc body", frm(fm4, family = cumulative(), data = d))
