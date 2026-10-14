# Review round 2: nan_unless_pos() on the tape at small, denormal and
# huge x; the wiener and frechet cases of fixes 1 and 2; a fit inside
# the support under the guard; ad_tree_extreme() gradients at ties.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("frmtmb from", find.package("frmtmb"), "\n")
g <- frmtmb:::nan_unless_pos

cat("\n=== A. nan_unless_pos on the tape: value, gradient, Hessian\n")
xs <- c(-1, -1e-300, -0, 0, 4.9e-324, 1e-310, 2.3e-308, 1e-300, 1e-20, 1,
        1e300, 8.9e307, 9e307, 1e308)
tp <- RTMB::MakeTape(function(p) g(p), 1)
for (x in xs) {
  v <- tp(x)
  j <- tp$jacobian(x)
  h <- tryCatch(RTMB::MakeTape(function(p) tp$jacfun()(p), 1)$jacobian(x),
                error = function(e) NA)
  cat(sprintf("x = %-10g double %-4s tape %-4s grad %-5s hess %-5s\n", x,
              format(g(x)), format(v), format(j), format(h)))
}
# inside a sum with other terms, as a density uses it: the gradient of
# p^2 + guard(p) must be 2p exactly
tp2 <- RTMB::MakeTape(function(p) p^2 + g(p), 1)
for (x in c(1e-300, 4.9e-324, 0.5)) {
  cat(sprintf("d/dp [p^2 + guard(p)] at %g: %s (want %g)\n", x,
              format(tp2$jacobian(x)), 2 * x))
}
# can u / u fold? x a tape variable, u built from it, divided by itself
tp3 <- RTMB::MakeTape(function(p) { u <- p + abs(p); u / u }, 1)
cat("u / u on the tape at p = -3:", tp3(-3), "(NaN means no fold)\n")

cat("\n=== B. wiener identity links: objective outside the support\n")
set.seed(20261016)
d <- frmtmb.eam:::ddm_simulate(60, mu = 0.8, bs = 1.3, ndt = 0.3)
fid <- frm(bf(rt | dec(upper) ~ 1),
           family = wiener(link_bs = "identity", link_ndt = "identity",
                           link_bias = "identity"), data = d)
par <- fid$opt$par
cat("identity fit logLik", as.numeric(logLik(fid)), "(was -14.53702)\n")
for (j in 1:3) for (v in c(1.2, -0.5, 1e-300)) {
  p2 <- par
  p2[which(names(par) == "betad")[j]] <- v
  cat(sprintf("betad[%d] = %g: objective %s\n", j, v, format(fid$obj$fn(p2))))
}

cat("\n=== C. frechet(link_nu = 'identity') at nu = 0.8 on the tape\n")
fam <- frechet(link_nu = "identity")
y <- c(0.5, 1, 3)
f <- function(p) {
  fam$lpdf(y, list(mu = exp(p[1]), nu = p[2], .eta_mu = p[1], .eta_nu = p[2]),
           list())
}
tpf <- RTMB::MakeTape(f, c(log(2), 3))
cat("nu = 0.8:", tpf(c(log(2), 0.8)), "; nu = 1:", tpf(c(log(2), 1)),
    "; nu = 1 + 1e-12:", tpf(c(log(2), 1 + 1e-12)), "; nu = 3:",
    tpf(c(log(2), 3)), "\n")

cat("\n=== D. a fit inside the support is unchanged by the guard\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(623)
y <- rg(60, 0, 1, -0.6)
fi <- frm(y ~ 1, family = gen_extreme_value(link_xi = "identity"),
          data = data.frame(y = y))
cat(sprintf("GEV identity xi, seed 623: logLik %.10f (round 1: -57.9967) xi %.6f conv %d\n",
            as.numeric(logLik(fi)), summary(fi)$spec_pars["xi", 1],
            fi$opt$convergence))
set.seed(801)
nu0 <- 1.15
dd <- data.frame(y = 2 / gamma(1 - 1 / nu0) * rexp(200)^(-1 / nu0))
ff <- frm(y ~ 1, family = frechet(link_nu = "identity"), data = dd)
cat(sprintf("frechet identity nu, seed 801: logLik %.10f (round 1: -150.5084) nu %.6f\n",
            as.numeric(logLik(ff)), ff$opt$par[2]))

cat("\n=== E. gev_scale_xi() at ties of the min and the max\n")
sx <- frmtmb:::gev_scale_xi
for (yy in list(c(-1, -1, 0.3, 2, 2), c(-1, 0.3, 2), c(-1, -1, -1, 2, 2, 2, 2))) {
  fn <- function(p) sx(p[3], yy, p[1], exp(p[2]))
  p0 <- c(0.1, log(1.2), 0.4)
  tpx <- RTMB::MakeTape(fn, p0)
  ja <- tpx$jacobian(p0)
  h <- 1e-6
  jn <- vapply(1:3, function(k) {
    e <- replace(numeric(3), k, h)
    (fn(p0 + e) - fn(p0 - e)) / (2 * h)
  }, 0)
  cat(sprintf("y = (%s): tape grad %s | central diff %s | max rel %.2e\n",
              paste(yy, collapse = ", "), paste(signif(ja, 8), collapse = " "),
              paste(signif(jn, 8), collapse = " "),
              max(abs(ja - jn) / pmax(1e-12, abs(jn)))))
}
# and a tie between a row and the -1 hold: x such that -1 / max(x) = -1
yy <- c(-1, 0.3, 1.1)
fn <- function(p) sx(p[3], yy, p[1], exp(p[2]))
p0 <- c(0.1, log(1), 0.4)
cat("hold tie (max x = 1): value", fn(p0), " tape grad",
    signif(RTMB::MakeTape(fn, p0)$jacobian(p0), 6), "\n")
