# Lane selimit2: the three expectations of test-gp-noncentered.R that
# failed once inside R CMD check --as-cran (lines 481, 527, 528) and pass
# one file per process and in one session. Prints the measured values
# against the tests' bounds, to see whether they sit near a bound.
#   Rscript dev/selimit2-gpnc.R <lib> <tag>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(4402)
n <- 40L
x <- sort(round(stats::runif(n, 0, 6), 2))
d <- data.frame(x = x, y = stats::rpois(n, exp(1.5 + sin(2 * x))))
fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = poisson()))
sm <- frmtmb:::smooth_b_idx(fit)
ch <- frmtmb:::fd_eta_chain(fit, NULL, NULL, FALSE, TRUE, sm)
f <- function(g) frmtmb:::fitted_point(g, NULL, NA)
fd <- as.vector(frmtmb:::fit_fd_se(fit, f, b_idx = sm, b_chain = ch$chain))
de <- frm_linpred(fit, type = "response", se.fit = TRUE)$se.fit
r481 <- max(abs(fd / de - 1))
set.seed(20261102)
n <- 40L
x <- round(stats::runif(n, 0, 6), 3)
y <- sin(2 * x) + 1e-2 * stats::rnorm(n)
u <- sort(unique(x))
J <- match(x, u)
D2 <- outer(u, u, "-")^2
stub <- list(re_blocks = list(list(covstruct = "gp")))
ff <- function(p) {
  "[<-" <- RTMB::ADoverload("[<-")
  L <- t(chol(exp(-D2 / (2 * exp(2 * p$lrho))) + diag(1e-12, length(u))))
  eta <- p$mu + exp(p$lsd) * as.vector(L %*% p$b)[J]
  -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) -
    sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
}
pars <- list(mu = 0, ls = log(0.01), lsd = 1.08, lrho = 0.3,
             b = numeric(length(u)))
at <- c(ls = log(0.01), lsd = 1.08396, lrho = 0.27)
plain <- RTMB::MakeADFun(ff, pars, random = c("mu", "b"), silent = TRUE)
g0 <- plain$gr(at)
v0 <- plain$fn(at)
fast <- frmtmb:::gp_tape_scope(stub, RTMB::MakeADFun(
  ff, frmtmb:::gp_tape_order(stub, pars), random = c("mu", "b"),
  silent = TRUE))
r527 <- abs(fast$fn(at) - v0) / abs(v0)
r528 <- max(abs(fast$gr(at) - g0)) / max(abs(g0))
cat(sprintf("%s | 481: %.3e (bound 1e-3) | 527: %.3e (bound 1e-8) | 528: %.3e (bound 1e-5) | code %d\n",
            a[2], r481, r527, r528, fit$opt$convergence))
