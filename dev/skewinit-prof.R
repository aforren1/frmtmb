# The one job in the no-move battery where the fixed build comes out
# BELOW the base build: skew_normal with a mu covariate under
# profile = TRUE. What is the optimum, and does either build reach it?
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

set.seed(11)
n <- 300
x <- rnorm(n)
z <- runif(n)
g <- factor(rep(1:20, length.out = n))
gg <- rnorm(20, 0, 0.7)[as.integer(g)]
eta <- 0.4 + 0.8 * x
dd <- data.frame(
  x = x, z = z, g = g,
  ygau = eta + gg + rnorm(n, 0, 0.9),
  ypos = rgamma(n, shape = 3, scale = exp(eta) / 3),
  ycnt = rpois(n, exp(0.5 + 0.4 * x)),
  ybin = rbinom(n, 1, plogis(eta)),
  ybeta = rbeta(n, 2 * plogis(eta) * 5, 2 * (1 - plogis(eta)) * 5),
  ysn0 = 2 + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5,
  ysn1 = eta + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5
)
fo <- bf(ysn1 ~ x, sigma ~ 1, alpha ~ 1)
r <- residuals(stats::lm(ysn1 ~ x, data = dd))
op <- options(digits = 12)

one <- function(lab, ...) {
  f <- suppressWarnings(frm(fo, family = skew_normal(), data = dd, ...))
  e <- f$opt[["stationary_escape"]]
  cat(sprintf("%-26s ll %18.10f  alpha %12.5f  max|grad| %10.4g  esc %s\n",
              lab, as.numeric(logLik(f)),
              f$estimates$betad[["alpha_(Intercept)"]],
              max(abs(f$obj$gr(f$opt$par))),
              if (is.null(e)) "-" else format(e[["gain"]], digits = 4)))
  invisible(f)
}
cat("sd(y)", stats::sd(dd$ysn1), " sd(resid)", stats::sd(r),
    " skew(y)", skew(dd$ysn1), " skew(resid)", skew(r), "\n\n")
cat("== ML\n")
one("default")
one("start raw sd(y), a=+2", start = list(betad = c(log(stats::sd(dd$ysn1)), 2)))
one("start resid sd, a=+2", start = list(betad = c(log(stats::sd(r)), 2)))
one("start resid sd, a=-2", start = list(betad = c(log(stats::sd(r)), -2)))
cat("\n== profile = TRUE\n")
pf <- frmtmb_control(profile = TRUE)
one("default", control = pf)
one("start raw sd(y), a=+2", control = pf,
    start = list(betad = c(log(stats::sd(dd$ysn1)), 2)))
one("start resid sd, a=+2", control = pf,
    start = list(betad = c(log(stats::sd(r)), 2)))
one("start resid sd, a=-2", control = pf,
    start = list(betad = c(log(stats::sd(r)), -2)))
cat("\n== profile = TRUE, tighter optimizer\n")
pf2 <- frmtmb_control(profile = TRUE, restarts = 5, grad_tol = 1e-6)
one("default, 5 restarts", control = pf2)
cat("\n== sn::selm\n")
m <- sn::selm(ysn1 ~ x, family = "SN", data = dd)
cat("sn logLik", as.numeric(m@logL), " alpha",
    as.numeric(m@param$dp[["alpha"]]), "\n")
options(op)
