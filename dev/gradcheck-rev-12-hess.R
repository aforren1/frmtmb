## Reviewer, claim 1(b): WHERE the stage-2 Hessian comes from, and
## whether it is the Hessian of the MARGINAL objective on a Laplace
## model. grad_headroom() prefers fit$obj$he(p) and only falls back to
## stats::optimHess() on the marginal gradient. If he() on a random
## effects object returned anything other than the marginal Hessian of
## the outer parameters, every GLMM headroom would be wrong.
##
## Also the direct cost of the stage, timed in blocks past 1.2 s with a
## minimum over rounds, against one whole fit timed the same way.
## usage: Rscript gradcheck-rev-12-hess.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x

## a Laplace GLMM, small enough to difference by hand
set.seed(8501)
ng <- 50
dd <- data.frame(g = factor(rep(seq_len(ng), 20)))
dd$x <- rnorm(nrow(dd))
re <- rnorm(ng, 0, 0.7)
dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
f <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = poisson(),
                          data = dd))
p <- f$opt$par
np <- length(p)
cat("outer np =", np, "  random effects =",
    length(f$obj$env$random %||% integer(0)), "\n")
he <- tryCatch(f$obj$he(p), error = function(e) conditionMessage(e))
if (is.character(he)) {
  cat("obj$he() on the Laplace object ERRORS:", he, "\n")
} else {
  cat("obj$he() returns dim", paste(dim(he), collapse = "x"),
      " (outer np is", np, ")\n")
}
Hfd <- stats::optimHess(p, function(q) as.numeric(f$obj$fn(q)),
                        function(q) drop(f$obj$gr(q)))
cat("optimHess(marginal gradient) dim", paste(dim(Hfd), collapse = "x"), "\n")
## an independent central difference of the marginal gradient, so the
## comparison does not rest on optimHess's own step rule
eps <- 1e-5
Hcd <- matrix(NA_real_, np, np)
for (j in seq_len(np)) {
  qp <- p; qp[j] <- qp[j] + eps
  qm <- p; qm[j] <- qm[j] - eps
  Hcd[, j] <- (drop(f$obj$gr(qp)) - drop(f$obj$gr(qm))) / (2 * eps)
}
Hcd <- (Hcd + t(Hcd)) / 2
Hfs <- (Hfd + t(Hfd)) / 2
cat("max relative gap optimHess vs central difference:",
    format(max(abs(Hfs - Hcd)) / max(abs(Hcd)), digits = 4), "\n")
if (!is.character(he) && identical(dim(he), c(np, np))) {
  hs <- (he + t(he)) / 2
  cat("max relative gap he() vs central difference:",
      format(max(abs(hs - Hcd)) / max(abs(Hcd)), digits = 4),
      "  <- must be small, or he() is NOT the marginal Hessian\n")
}
g <- drop(f$obj$gr(p))
hr <- frmtmb:::grad_headroom(f, g, rep(FALSE, np))
chd <- chol(Hcd)
cat("headroom from grad_headroom():",
    format(hr, digits = 10), "\n")
cat("headroom from the central-difference Hessian:",
    format(0.5 * sum(backsolve(chd, g, transpose = TRUE)^2), digits = 10),
    "\n\n")

## ---- the direct cost of the stage, on the two models the brief names
bench <- function(tag, dd, form, fam, extra = list()) {
  ctl <- do.call(frmtmb_control, c(list(restarts = 0), extra))
  f <- suppressWarnings(frm(form, family = fam, data = dd, control = ctl))
  p <- f$opt$par; np <- length(p)
  g <- drop(f$obj$gr(p))
  act <- rep(FALSE, np)
  blk <- function(fn) {
    k <- 1L
    repeat {
      t0 <- proc.time()[["elapsed"]]
      for (i in seq_len(k)) fn()
      el <- proc.time()[["elapsed"]] - t0
      if (el > 1.2 || k > 4096L) return(el / k)
      k <- k * 2L
    }
  }
  rounds <- 3L
  hd <- fitt <- grt <- rep(NA_real_, rounds)
  for (r in seq_len(rounds)) {
    hd[r] <- blk(function() frmtmb:::grad_headroom(f, g, act))
    grt[r] <- blk(function() f$obj$gr(p))
    t0 <- proc.time()[["elapsed"]]
    invisible(suppressWarnings(frm(form, family = fam, data = dd,
                                  control = ctl)))
    fitt[r] <- proc.time()[["elapsed"]] - t0
  }
  cat(sprintf(paste0("%-16s np=%-3d n=%-6d one grad %.5f s  ",
                     "grad_headroom %.5f s  whole fit %.4f s  ",
                     "stage2/fit %.4f\n"),
              tag, np, nrow(dd), min(grt), min(hd), min(fitt),
              min(hd) / min(fitt)))
}
set.seed(8502)
ng <- 5000
d5 <- data.frame(g = factor(rep(seq_len(ng), 4)))
d5$x <- rnorm(nrow(d5))
re <- rnorm(ng, 0, 0.7)
d5$y <- rpois(nrow(d5), exp(0.3 + 0.5 * d5$x + re[d5$g]))
bench("glmm 5000 re", d5, bf(y ~ x + (1 | g)), poisson())

set.seed(8503)
d6 <- data.frame(g = factor(rep(seq_len(ng), 4)))
d6$x <- rnorm(nrow(d6))
rs <- rnorm(ng, 0, 0.4)
d6$y <- rpois(nrow(d6), exp(0.3 + (0.5 + rs[d6$g]) * d6$x + re[d6$g]))
bench("glmm 5000 slope", d6, bf(y ~ x + (x | g)), poisson())

set.seed(8504)
n <- 20000
d7 <- data.frame(x = sort(runif(n, 0, 8)), z = rnorm(n))
d7$y <- rnorm(n, sin(d7$x) + 0.4 * d7$z, 0.5)
bench("smooth n=20000", d7, bf(y ~ s(x) + z), gaussian())

set.seed(8505)
n <- 700
d8 <- data.frame(t = sort(runif(n, 0, 12)))
d8$y <- rnorm(n, sin(d8$t) + 0.2 * d8$t, 0.3)
bench("gp k=40 n=700", d8, bf(y ~ gp(t, k = 40)), gaussian())

## a many-parameter LMM, the worker's worst case
set.seed(8506)
ng <- 200
d9 <- data.frame(g = factor(rep(seq_len(ng), 20)))
for (j in 1:6) d9[[paste0("x", j)]] <- rnorm(nrow(d9))
b <- rnorm(6, 0, 0.4)
d9$y <- rnorm(nrow(d9),
              as.matrix(d9[, paste0("x", 1:6)]) %*% b +
                rnorm(ng, 0, 0.6)[d9$g], 1)
bench("LMM us 6x6", d9,
      bf(y ~ x1 + x2 + x3 + x4 + x5 + x6 +
           (1 + x1 + x2 + x3 + x4 + x5 | g)), gaussian())
cat("DONE hess\n")
