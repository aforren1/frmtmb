# Lane arstep: the student step rows. Default fit, then the review's 15
# reference starts with where each ends, and the working marginal's local
# maxima. Data as dev/review-ncgp-r3-ng.R (seed 74000 + 31 seed + n).
# Args: arm case n seed [curve]
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; cs <- args[2]; n <- as.integer(args[3]); seed <- as.integer(args[4])
cn <- if (length(args) > 4) args[5] else "step"
rev_libs(arm)
ns <- asNamespace("frmtmb")
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x, two = function(x) sin(x) + 0.4 * sin(10 * x), smooth = function(x) sin(x), step = function(x) 2 * (x > 5) - 1)
set.seed(74000 + seed * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
eta <- curves[[cn]](x)
m <- switch(cs,
  student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
  gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
fitit <- function(start = NULL) { a <- list(m$fml, m$d, family = m$fam); if (!is.null(start)) a$start <- start; do.call(rev_fit, a) }
show <- function(lbl, r) {
  e <- r$fit$estimates
  cat(sprintf("D3 %-22s ll %.4f conv %d w %d | sdgp %.4g lscale %.4g | betad %s | beta %s %s\n", lbl,
              as.numeric(logLik(r$fit)), r$fit$opt$convergence, length(r$warn), exp(e$theta[1]), exp(e$theta[2]),
              paste(signif(e$betad, 5), collapse = ","), paste(signif(e$beta, 5), collapse = ","),
              if (!is.null(e$thetaac)) paste0("ac ", signif(e$thetaac, 5)) else ""))
}
r0 <- fitit(); show("default", r0)
w <- ns$gp_start_working(r0$fit$frame)
cat(sprintf("D3 working start sdgp %.4g lscale %.4g; alts: %s\n", exp(w$theta[1]), exp(w$theta[2]),
            paste(vapply(w$alts, function(a) sprintf("(%.4g, %.4g)", exp(a[1]), exp(a[2])), ""), collapse = " ")))
for (rho0 in c(0.1, 0.3, 1, 3, 8)) for (sd0 in c(0.3, 1.5, 6)) {
  r2 <- fitit(list(theta = c(log(sd0), log(rho0))))
  if (!inherits(r2$fit, "error")) show(sprintf("start sd %g lscale %g", sd0, rho0), r2)
}
cat("DONE\n")
