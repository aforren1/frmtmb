# p069's own data: the importance fit against my own exact marginal (product
# Gauss-Hermite around each group's mode), at 200, 512 and 2000 draws and over
# draw seeds. The correlation's chart sits at the boundary, so the reference is
# profiled over the sds and the means with the correlation at its limit, and
# with the correlation free on a bounded scale.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/equiv-specs.R", envir = E)
d <- E$eq_dat("poisson", seed = 69)
dl <- split(d, d$g)
g1 <- gh_rule(31L)
grid <- as.matrix(expand.grid(seq_len(31L), seq_len(31L)))
xn <- matrix(g1$x[grid], nrow(grid), 2L)
lw <- rowSums(matrix(log(g1$w)[grid], nrow(grid), 2L)) + rowSums(xn^2)
ll <- function(p, dg, U) {
  eta <- outer(p[1] + p[2] * dg$x, rep(1, nrow(U))) + cbind(1, dg$x) %*% t(U)
  colSums(matrix(dpois(rep(dg$y, nrow(U)), exp(as.vector(eta)), log = TRUE),
                 nrow(dg), nrow(U)))
}
nll <- function(p) {
  # p: b0, b1, log sd1, log sd2, atanh(cor)
  s <- exp(p[3:4])
  r <- tanh(p[5])
  S <- diag(s) %*% matrix(c(1, r, r, 1), 2) %*% diag(s) + diag(1e-12, 2)
  ch <- chol(S)
  Pi <- chol2inv(ch)
  ld <- 2 * sum(log(diag(ch)))
  lpr <- function(U) -0.5 * rowSums((U %*% Pi) * U) - log(2 * pi) - 0.5 * ld
  tot <- 0
  for (dg in dl) {
    neg <- function(u) {
      U <- matrix(u, 1L, 2L)
      -(ll(p, dg, U) + lpr(U))
    }
    o <- optim(c(0, 0), neg, method = "BFGS", control = list(reltol = 1e-11))
    h <- optimHess(o$par, neg)
    Lc <- backsolve(chol(h), diag(2))
    U <- matrix(o$par, nrow(xn), 2L, byrow = TRUE) + sqrt(2) * xn %*% t(Lc)
    a <- ll(p, dg, U) + lpr(U) + lw
    tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(2) + sum(log(abs(diag(Lc)))))
  }
  tot
}
lap <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)), family = poisson(), data = d))
v <- frmtmb:::varcorr_matrices(lap)[[1L]]
st <- c(lap$opt$par[1:2], log(sqrt(diag(v))), atanh(min(0.999, v[1, 2] / sqrt(v[1, 1] * v[2, 2]))))
ref <- nlminb(unname(st), nll, lower = c(-Inf, -Inf, -Inf, -Inf, -6),
              upper = c(Inf, Inf, Inf, Inf, 6), control = list(iter.max = 80))
cat("exact: par", format(ref$par, digits = 6), "ll", format(-ref$objective, digits = 9),
    "cor", tanh(ref$par[5]), "\n")
cat("Laplace ll", format(as.numeric(logLik(lap)), digits = 9), "\n")
for (N in c(200L, 512L, 2000L)) {
  for (s in c(1L, 2L, 3L)) {
    w <- character(0)
    f <- withCallingHandlers(
      frm(bf(y ~ x + (1 + x | g)), family = poisson(), data = d, importance = N,
          control = frmtmb_control(importance_seed = s)),
      warning = function(x) {
        w <<- c(w, substr(conditionMessage(x), 1, 60))
        invokeRestart("muffleWarning")
      })
    vv <- frmtmb:::varcorr_matrices(f)[[1L]]
    cat(sprintf("N %d seed %d ll %.5f mcse %.4f err/mcse %.2f rounds %d max move %.3g ess %.3f sd %.4f %.4f cor %.4f warns %s\n",
                N, s, as.numeric(logLik(f)), f$importance$mcse,
                (as.numeric(logLik(f)) + ref$objective) / f$importance$mcse,
                f$importance$rounds, max(f$importance$moves), f$importance$ess_min,
                sqrt(vv[1, 1]), sqrt(vv[2, 2]), vv[1, 2] / sqrt(vv[1, 1] * vv[2, 2]),
                paste(w, collapse = "; ")))
  }
}
