# Reviewer probe: frechet() under link_nu = "identity" on the tape,
# fitted() as the mean, and recovery with nu near 1.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

cat("\n=== A. nu <= 1 under link_nu = 'identity'\n")
fam <- frechet(link_nu = "identity")
y <- c(0.5, 1, 3)
pl <- fam$lpdf(y, list(mu = 2, nu = 0.8, .eta_mu = log(2), .eta_nu = 0.8),
               list())
cat("plain doubles, nu = 0.8:", pl, "\n")
f <- function(p) {
  fam$lpdf(y, list(mu = exp(p[1]), nu = p[2], .eta_mu = p[1], .eta_nu = p[2]),
           list())
}
tp <- RTMB::MakeTape(f, c(log(2), 3))
cat("on the tape, nu = 0.8:", tp(c(log(2), 0.8)), "\n")
cat("on the tape, nu = 3  :", tp(c(log(2), 3)), "\n")
cat("plain, nu = 3        :",
    fam$lpdf(y, list(mu = 2, nu = 3, .eta_mu = log(2), .eta_nu = 3), list()),
    "\n")
# a fit whose optimizer could step below 1
set.seed(801)
n <- 200
nu0 <- 1.15
dd <- data.frame(y = 2 / gamma(1 - 1 / nu0) * rexp(n)^(-1 / nu0))
r <- tryCatch(withCallingHandlers(
  frm(y ~ 1, family = frechet(link_nu = "identity"), data = dd),
  warning = function(w) { cat("WARN:", substr(conditionMessage(w), 1, 150), "\n"); invokeRestart("muffleWarning") }),
  error = function(e) e)
if (inherits(r, "error")) cat("ERROR", conditionMessage(r), "\n") else {
  cat("identity fit nu0 1.15: conv", r$opt$convergence, "logLik",
      as.numeric(logLik(r)), "par", r$opt$par, "\n")
}
r2 <- frm(y ~ 1, family = frechet(), data = dd)
cat("logm1 fit   nu0 1.15: conv", r2$opt$convergence, "logLik",
    as.numeric(logLik(r2)), "nu",
    frm_linpred(r2, dpar = "nu", type = "response")[1], "\n")

cat("\n=== B. fitted() is the mean\n")
set.seed(802)
dd <- data.frame(x = rnorm(300))
nu <- 3
dd$y <- exp(1 + 0.3 * dd$x) / gamma(1 - 1 / nu) * rexp(300)^(-1 / nu)
fit <- frm(y ~ x, family = frechet(), data = dd)
fm <- fitted(fit)[, "Estimate"]
mu <- frm_linpred(fit, dpar = "mu", type = "response")
nuh <- frm_linpred(fit, dpar = "nu", type = "response")
qm <- vapply(1:3, function(i) {
  s <- mu[i] / gamma(1 - 1 / nuh[i])
  integrate(function(u) s * u^(-1 / nuh[i]) * exp(-u), 0, Inf,
            rel.tol = 1e-12)$value
}, 0)
cat("fitted - mu, max abs rel:", max(abs(fm - mu) / mu),
    "; mu vs quadrature mean of the fitted density, rows 1:3, max rel:",
    max(abs(mu[1:3] - qm) / qm), "\n")
cat("predict() mean of draws vs fitted, rel:",
    abs(mean(simulate(fit, nsim = 200, seed = 1)[[1]]) - fm[1]) / fm[1],
    "(one row, one sim column; noisy)\n")

cat("\n=== C. recovery with nu near 1, 60 seeds each, n = 300\n")
for (nu0 in c(1.1, 1.3)) {
  res <- t(vapply(1:60, function(s) {
    set.seed(9000 + s)
    d <- data.frame(x = rnorm(300))
    d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(300)^(-1 / nu0)
    w <- 0
    f <- withCallingHandlers(frm(y ~ x, family = frechet(), data = d),
                             warning = function(e) { w <<- w + 1; invokeRestart("muffleWarning") })
    cf <- fixef(f)
    sm <- summary(f)$spec_pars
    c(conv = f$opt$convergence, warn = w, b = cf["x", 1], se = cf["x", 2],
      nu = sm["nu", 1], nuse = sm["nu", 2], a = cf["Intercept", 1],
      ase = cf["Intercept", 2])
  }, numeric(8)))
  cat(sprintf("nu0 %.1f: conv0 %d/60, warned %d; x: mean %.4f empSD %.4f meanSE %.4f cover %.3f; nu: mean %.4f empSD %.4f meanSE %.4f cover %.3f; Intercept: mean %.4f empSD %.4f meanSE %.4f cover %.3f\n",
              nu0, sum(res[, "conv"] == 0), sum(res[, "warn"] > 0),
              mean(res[, "b"]), sd(res[, "b"]), mean(res[, "se"]),
              mean(abs(res[, "b"] - 0.3) < 1.96 * res[, "se"]),
              mean(res[, "nu"]), sd(res[, "nu"]), mean(res[, "nuse"]),
              mean(abs(res[, "nu"] - nu0) < 1.96 * res[, "nuse"]),
              mean(res[, "a"]), sd(res[, "a"]), mean(res[, "ase"]),
              mean(abs(res[, "a"] - 0.5) < 1.96 * res[, "ase"])))
}
