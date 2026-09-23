# Two open questions after the first pass of punch-round-1 fixes.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
op <- options(digits = 12)

cat("\n=== Q1: is `alpha ~ 0 + x` at ITS OWN optimum?\n")
set.seed(7)
n <- 300
x <- rnorm(n)
y <- 1 + 0.5 * x + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5
dd <- data.frame(y = y, x = x)
fo <- bf(y ~ x, sigma ~ 1, alpha ~ 0 + x)
def <- suppressWarnings(frm(fo, family = skew_normal(), data = dd))
cat(sprintf("  default          ll %16.7f  b %10.4f  esc %s\n",
            as.numeric(logLik(def)), def$estimates$betad[[2]],
            if (is.null(def$opt[["stationary_escape"]])) "-" else "fired"))
best <- as.numeric(logLik(def)); bb <- def$estimates$betad[[2]]
for (b0 in c(-20, -8, -4, -2, -1, -0.25, 0.25, 1, 2, 4, 8, 20)) {
  f <- try(suppressWarnings(frm(fo, family = skew_normal(), data = dd,
                                start = list(betad = c(0, b0)))), silent = TRUE)
  if (inherits(f, "try-error")) next
  ll <- as.numeric(logLik(f))
  cat(sprintf("  start b=%-6g     ll %16.7f  b %10.4f\n", b0, ll,
              f$estimates$betad[[2]]))
  if (ll > best) { best <- ll; bb <- f$estimates$betad[[2]] }
}
cat(sprintf("  BEST over starts ll %16.7f  b %10.4f   default short by %.6f\n",
            best, bb, best - as.numeric(logLik(def))))

cat("\n=== Q2: what causes the M3 convergence-1 warning?\n")
set.seed(8)
nb <- 200
xb <- rnorm(nb, 0, 1e4)
yb <- 0.001 * xb + RTMBdist::rskewnorm2(nb, 0, 1.5, 4)
d3 <- data.frame(y = yb, x = xb)
f0 <- bf(y ~ x, sigma ~ 1, alpha ~ 1)
r <- residuals(stats::lm(y ~ x, data = d3))
sk <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
m3r <- sk(r); m3y <- sk(d3$y)
cat(sprintf("  start LANE  sigma %10.5f alpha %8.4f\n",
            log(stats::sd(r)), 2 * sign(m3r) + 0.5 * m3r))
cat(sprintf("  start BASE  sigma %10.5f alpha %8.4f\n",
            log(stats::sd(d3$y)), 2 * sign(m3y) + 0.5 * m3y))
probe <- function(lab, ...) {
  w <- 0L
  f <- withCallingHandlers(frm(f0, family = skew_normal(), data = d3, ...),
                           warning = function(x) {
                             w <<- w + 1L; invokeRestart("muffleWarning")
                           })
  g <- max(abs(f$obj$gr(f$opt$par)))
  u <- f$par_units
  cat(sprintf("  %-34s ll %15.8f conv %d warn %d max|grad| %9.3g units %s\n",
              lab, as.numeric(logLik(f)), f$opt$convergence, w, g,
              if (is.null(u)) "none" else paste(format(range(u), digits = 3),
                                                collapse = "..")))
  invisible(f)
}
probe("default (lane start, autoscale auto)")
probe("lane start, autoscale FALSE",
      control = frmtmb_control(autoscale = FALSE))
probe("BASE start, autoscale auto",
      start = list(betad = c(log(stats::sd(d3$y)),
                             2 * sign(m3y) + 0.5 * m3y)))
probe("BASE start, autoscale FALSE",
      start = list(betad = c(log(stats::sd(d3$y)), 2 * sign(m3y) + 0.5 * m3y)),
      control = frmtmb_control(autoscale = FALSE))
probe("lane start, more restarts",
      control = frmtmb_control(restarts = 5, grad_tol = 1e-6))
options(op)
