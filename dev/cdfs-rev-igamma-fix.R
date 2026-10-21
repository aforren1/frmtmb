# Reviewer: a candidate fix for the log_igamma_two() NaN, applied in
# this session only (assignInNamespace), to show the failing fits are
# that defect: cap the exponent of the quadrature integrand, whose
# window crosses the mode when x is below it; that arm is weighted 0
# there, so any finite value serves.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
patched <- ns$igamma_gl_tails
b <- paste(deparse(body(patched), width.cutoff = 500L), collapse = "
")
b <- gsub("exp(lf(lo + (x - lo) * gl40$x[k]) - lfx)",
          "exp(ad_min(lf(lo + (x - lo) * gl40$x[k]) - lfx, 600))", b, fixed = TRUE)
b <- gsub("exp(lf(x + wup * gl40$x[k]) - lfx)",
          "exp(ad_min(lf(x + wup * gl40$x[k]) - lfx, 600))", b, fixed = TRUE)
stopifnot(lengths(regmatches(b, gregexpr("ad_min(lf", b, fixed = TRUE))) == 2)
body(patched) <- parse(text = b)[[1]]
environment(patched) <- ns
utils::assignInNamespace("igamma_gl_tails", patched, "frmtmb")
cat("== patched: NaN map and accuracy vs stats::pgamma\n")
worst <- 0
for (a in c(1.3e4, 1e5, 1e6)) {
  r <- c(0.01, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.5)
  p <- ns$log_igamma_two(r * a, a)
  rF <- stats::pgamma(r * a, a, log.p = TRUE)
  rS <- stats::pgamma(r * a, a, lower.tail = FALSE, log.p = TRUE)
  e <- pmax(abs(p$lP - rF) / pmax(1, abs(rF)), abs(p$lQ - rS) / pmax(1, abs(rS)))
  cat(sprintf("a %g: NaN %d, worst rel log err %.2e\n", a,
              sum(is.nan(p$lP) | is.nan(p$lQ)), max(e)))
}
cat("== patched: gradient on the tape at a = 1e5, x = 0.7 a (d/d log a, d/d log x)\n")
f <- function(p) {
  g <- ns$log_igamma_two(exp(p[2]), exp(p[1]))
  g$lP + g$lQ
}
tp <- RTMB::MakeTape(f, c(log(1e5), log(7e4)))
print(tp$jacobian(c(log(1e5), log(7e4))))
cat("== patched: the two fits that failed\n")
set.seed(3)
d <- data.frame(x = rnorm(200))
d$y <- rgamma(200, shape = 1e5, rate = 1e5 / exp(1 + 0.1 * d$x))
qq <- quantile(d$y, 0.8)
d$cc <- as.integer(d$y > qq)
d$y[d$cc == 1] <- qq
fit <- frm(y | cens(cc) ~ x, family = Gamma(link = "log"), data = d)
cat(sprintf("Gamma shape 1e5: par %s logLik %.8f conv %d\n",
            paste(signif(fit$opt$par, 8), collapse = " "),
            as.numeric(logLik(fit)), fit$opt$convergence))
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
d$y <- rpois(n, exp(11 + 0.3 * d$x))
d$cc <- 0
d$cc[d$y > 75000] <- 1
d$y[d$cc == 1] <- 75000
r <- tryCatch({
  fit <- frm(y | cens(cc) ~ x, family = poisson(), data = d)
  sprintf("poisson large counts: par %s logLik %.6f conv %d",
          paste(signif(fit$opt$par, 10), collapse = " "),
          as.numeric(logLik(fit)), fit$opt$convergence)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, "\n")
