# Where would a non-centered exact gp() spend its time on the tape? Hand
# closures, no frmtmb: y ~ N(mu + f[J], sigma), n positions, smooth truth,
# seed 20261103. Arms: centered at 1e-6, non-centered at 1e-12 with RTMB's
# defaults, with TapeConfig(matmul = "atomic"), and the closed-form marginal.
#
# Usage: Rscript dev/ncgp-proto.R <n> <arm> [nrep]
# Prints NCGPP lines: tape build seconds, nodes per tape, peak working set,
# per-evaluation minimum of fn, gr at theta away from the last inner mode.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
library(RTMB)
args <- commandArgs(TRUE)
n <- as.integer(args[1])
arm <- args[2]
nrep <- if (length(args) > 2) as.integer(args[3]) else 5L

set.seed(20261103)
x <- sort(runif(n, 0, 10))
y <- sin(x) + 0.3 * rnorm(n)
D2 <- outer(x, x, "-")^2

if (grepl("_atomic$", arm)) TapeConfig(matmul = "atomic")
if (grepl("_vec$", arm)) TapeConfig(matmul = "atomic", vectorize = "enable")
if (grepl("^cmp", arm)) { TMB::config(tmbad.sparse_hessian_compress = 1L, DLL = "RTMB"); arm <- sub("^cmp_", "", arm) }

kern <- function(rho, nug) {
  K <- exp(-D2 / (2 * rho^2))
  K + RTMB::diag(nug, n)
}

closures <- list(
  centered = function(p) {
    s <- exp(p$lsd); r <- exp(p$lrho)
    K <- s^2 * kern(r, 1e-6)
    nll <- -sum(RTMB::dmvnorm(p$b, 0, K, log = TRUE))
    nll - sum(RTMB::dnorm(y, p$mu + p$b, exp(p$ls), log = TRUE))
  },
  ncp = function(p) {
    s <- exp(p$lsd); r <- exp(p$lrho)
    L <- t(chol(kern(r, 1e-12)))
    f <- s * as.vector(L %*% p$b)
    nll <- -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE))
    nll - sum(RTMB::dnorm(y, p$mu + f, exp(p$ls), log = TRUE))
  },
  marginal = function(p) {
    s <- exp(p$lsd); r <- exp(p$lrho)
    V <- s^2 * kern(r, 1e-12) + RTMB::diag(exp(2 * p$ls), n)
    -sum(RTMB::dmvnorm(y - p$mu, 0, V, log = TRUE))
  }
)
closures$ncp_atomic <- closures$ncp
closures$ncp_vec <- closures$ncp
closures$centered_atomic <- closures$centered

par0 <- list(mu = 0, lsd = 0, lrho = 0, ls = 0)
rnd <- NULL
if (!arm %in% "marginal") {
  par0$b <- numeric(n)
  rnd <- "b"
}
gc(reset = TRUE)
t0 <- proc.time()[[3]]
obj <- MakeADFun(closures[[arm]], par0, random = rnd, silent = TRUE)
tb <- proc.time()[[3]] - t0
nodes <- function(nm) {
  tp <- tryCatch(GetTape(obj, nm), error = function(e) NULL)
  if (is.null(tp)) return(NA_real_)
  nrow(tp$data.frame())
}
nn <- c(ADFun = nodes("ADFun"), ADGrad = nodes("ADGrad"), ADHess = nodes("ADHess"))
g1 <- gc()
th <- obj$par
minelapsed <- function(f, rounds = nrep) {
  out <- Inf
  for (i in seq_len(rounds)) {
    t0 <- Sys.time()
    k <- 0L
    while (TRUE) {
      f()
      k <- k + 1L
      el <- as.numeric(Sys.time() - t0, units = "secs")
      if (el > 1.2) break
    }
    out <- min(out, el / k)
  }
  out
}
# moved theta, so the inner solve happens: each eval starts from the last
# inner mode, as during an optimization
i <- 0L
wob <- function() {
  i <<- i + 1L
  th + 0.01 * sin(i) * seq_along(th) / length(th)
}
tf <- minelapsed(function() obj$fn(wob()))
tg <- minelapsed(function() obj$gr(wob()))
t0 <- proc.time()[[3]]
op <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 300))
topt <- proc.time()[[3]] - t0
pk <- sum(gc()[, 6])
cat(sprintf(paste0("NCGPP n=%d arm=%s build=%.2fs nodes(F,G,H)=%s fn=%.4fs ",
                   "gr=%.4fs opt=%.2fs iter=%d feval=%d geval=%d conv=%d ",
                   "obj=%.4f peakMb=%.0f\n"),
            n, arm, tb, paste(nn, collapse = ","), tf, tg, topt,
            op$iterations, op$evaluations[[1]], op$evaluations[[2]],
            op$convergence, op$objective, pk))
