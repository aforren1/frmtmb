# Reviewer of lane surface, claim 9: an ML fit whose true variance is 0.
# Does the floor create an optimum at the floor? One process, two arms
# (floor; then sd2_floored() replaced by exp(2 * log_sd), as 0.68.1).
# Per structure: the fitted log sd, logLik, convergence, identical()
# between arms, and the marginal objective along log sd at the
# estimates of the other parameters.
#   Rscript dev/surface-rev-zerovar.R > dev/surface-rev-out/zerovar.txt
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages(library(frmtmb))
set.seed(42)
n <- 120
d <- data.frame(x = round(runif(n, 0, 6), 2), y = rnorm(n),
                g = factor(rep(1:20, each = 6)),
                tt = rep(c(0, 0.3, 0.4, 1.1, 1.5, 2.7), 20),
                p1 = runif(n), p2 = runif(n))
d$tim <- num_factor(d$tt)
pc <- cbind(round(runif(25), 2) * 10, round(runif(25), 2) * 10)
d$pos <- num_factor(rep(pc[, 1], length.out = n), rep(pc[, 2], length.out = n))
d$one <- factor(rep(1, n))
d$sp <- factor(rep(1:20, each = 6))
A <- diag(20); dimnames(A) <- list(levels(d$sp), levels(d$sp))
A[A == 0] <- 0.3
ff <- list(
  gp = bf(y ~ gp(x)),
  ou = bf(y ~ 1 + ou(tim + 0 | g)),
  homcs = bf(y ~ 1 + homcs(tim + 0 | g)),
  exp = bf(y ~ 1 + exp(pos + 0 | one)),
  gr_cov = bf(y ~ 1 + (1 | gr(sp, cov = A))))
fit_arm <- function(nm) {
  f <- tryCatch(suppressWarnings(suppressMessages(
    frm(ff[[nm]], family = gaussian(), data = d, data2 = list(A = A)))),
    error = function(e) e)
  f
}
prof <- function(f, at = c(-20, -100, -300, -330, -372, -400, -1000)) {
  p <- f$obj$env$last.par.best
  op <- f$obj$par
  par <- p[names(p) != "b"]
  i <- which(names(par) == "theta")[1]
  par <- par[names(par) %in% names(op)]
  vapply(at, function(v) {
    q <- par; q[i] <- v
    tryCatch(f$obj$fn(q), error = function(e) NA_real_)
  }, 0)
}
res <- list()
for (arm in c("floor", "nofloor")) {
  if (arm == "nofloor") {
    assignInNamespace("sd2_floored", function(log_sd) exp(2 * log_sd),
                      "frmtmb")
  }
  for (nm in names(ff)) res[[arm]][[nm]] <- fit_arm(nm)
}
for (nm in names(ff)) {
  a <- res$floor[[nm]]; b <- res$nofloor[[nm]]
  if (inherits(a, "error") || inherits(b, "error")) {
    cat(sprintf("%-7s ERROR floor: %s | nofloor: %s\n", nm,
                if (inherits(a, "error")) conditionMessage(a) else "ok",
                if (inherits(b, "error")) conditionMessage(b) else "ok"))
    next
  }
  cat(sprintf(paste("%-7s log sd %.4f | logLik %.10f | conv %s |",
                    "identical est %s logLik %s\n"), nm,
              a$estimates$theta[1], as.numeric(logLik(a)),
              format(a$opt$convergence %||% NA),
              identical(a$estimates, b$estimates),
              identical(as.numeric(logLik(a)), as.numeric(logLik(b)))))
  cat("   objective (nll) at log sd -20 -100 -300 -330 -372 -400 -1000\n")
  cat("     floor  :", format(prof(a), digits = 12), "\n")
  cat("     nofloor:", format(prof(b), digits = 12), "\n")
}
