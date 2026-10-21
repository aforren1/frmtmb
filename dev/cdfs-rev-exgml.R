# Reviewer, fix round 1: the censored exgaussian case "u2" of
# cdfs-rev-brmsfit-r1.R (same seed and draws), an independent ML fit on
# brms::dexgaussian() / pexgaussian(), against frmtmb's.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
options(frmtmb.notices = FALSE)
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-cdfs/dev/cdfs-rev-brmsfit-r1.R")
stop_at <- grep("^summ_frm <- function", src) - 1
eval(parse(text = src[grep("^set.seed", src):stop_at]))
for (nm in c("exgaussian_u2", "exgaussian_u25")) {
  d <- cases[[nm]]$d
  nll <- function(p) {
    mu <- p[1] + p[2] * d$x; s <- exp(p[3]); b <- exp(p[4])
    -sum(ifelse(d$cc == 0, dexgaussian(d$y, mu, s, b, log = TRUE),
                pexgaussian(d$y, mu, s, b, lower.tail = FALSE, log.p = TRUE)))
  }
  best <- NULL
  for (lb in c(-6, -3, -1.5, -0.5)) {
    o <- nlminb(c(1, 0.4, log(0.5), lb), nll)
    if (is.null(best) || o$objective < best$objective) best <- o
    cat(sprintf("%s start log beta %g: logLik %.6f beta %.5f sigma %.5f\n", nm, lb, -o$objective,
                exp(o$par[4]), exp(o$par[3])))
  }
  fit <- frm(y | cens(cc, y2) ~ x, family = exgaussian(), data = d)
  cat(sprintf("%s frmtmb: logLik %.6f | best independent %.6f\n", nm, as.numeric(logLik(fit)), -best$objective))
}
cat("== profile of frmtmb's exact objective in log beta, case u2\n")
d <- cases[["exgaussian_u2"]]$d
fit <- frm(y | cens(cc, y2) ~ x, family = exgaussian(), data = d)
ob <- fit$obj
cat("par names:", names(fit$opt$par), " opt:", signif(fit$opt$par, 6), "\n")
j <- length(fit$opt$par)
for (lb in c(-12, -7, -4, -3, -2, -1.7, -1.5, -1)) {
  o <- nlminb(fit$opt$par[-j], function(q) ob$fn(c(q, lb)))
  cat(sprintf("beta %.5f: profile logLik %.6f\n", exp(lb), -o$objective))
}
