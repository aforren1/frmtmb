# Log density, log CDF and log survivor of frechet, gen_extreme_value and
# discrete_weibull against three references on a grid with edge values:
#   * brms 2.23.0's Stan functions (its generated functions block, and
#     Stan's built-in frechet), exposed with rstan::expose_stan_functions();
#   * evd (dgev, dfrechet);
#   * the same formulas in 256-bit arithmetic (Rmpfr), the arbiter where
#     two double-precision references disagree.
# frmtmb is evaluated off the tape (plain numbers) and on the tape (RTMB,
# with the linear predictors build_objective() stores beside each dpar).
# The error measure is |a - b| / max(1, |b|): relative for large values,
# absolute below 1.
# Output: dev/famuni-out/density-ref.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(rstan))

# brms's own function blocks, read off its generated programs
fb <- function(code) {
  lines <- strsplit(code, "\n")[[1]]
  i0 <- grep("^functions [{]", lines)
  i1 <- grep("^data [{]", lines) - 2L
  paste(lines[(i0 + 1L):i1], collapse = "\n")
}
set.seed(1)
d0 <- data.frame(y = rexp(20) + 1, yi = rpois(20, 3), x = rnorm(20))
gev_code <- suppressWarnings(as.character(brms::stancode(
  brms::bf(y ~ x), data = d0, family = brms::gen_extreme_value())))
dw_code <- as.character(brms::stancode(
  brms::bf(yi ~ x), data = d0, family = brms::brmsfamily("discrete_weibull")))
# scale_xi() is overloaded, which Rcpp cannot expose, and is not a density
gev_fun <- sub("(?s)[/][*] scale auxiliary parameter xi.*$", "", fb(gev_code),
               perl = TRUE)
fun <- paste0("functions {\n", gev_fun, "\n", fb(dw_code), "
  real fr_lpdf(real y, real nu, real s) { return frechet_lpdf(y | nu, s); }
  real fr_lcdf(real y, real nu, real s) { return frechet_lcdf(y | nu, s); }
  real fr_lccdf(real y, real nu, real s) { return frechet_lccdf(y | nu, s); }
  real gev_lpdf_w(real y, real mu, real s, real xi) {
    return gen_extreme_value_lpdf(y | mu, s, xi); }
  real gev_lcdf_w(real y, real mu, real s, real xi) {
    return gen_extreme_value_lcdf(y | mu, s, xi); }
  real gev_lccdf_w(real y, real mu, real s, real xi) {
    return gen_extreme_value_lccdf(y | mu, s, xi); }
  real dw_lpmf_w(int y, real mu, real k) {
    return discrete_weibull_lpmf(y | mu, k); }
  real dw_lcdf_w(int y, real mu, real k) {
    return discrete_weibull_lcdf(y | mu, k); }
  real dw_lccdf_w(int y, real mu, real k) {
    return discrete_weibull_lccdf(y | mu, k); }
}
model {}
")
env <- new.env()
# the C++ compile needs pkgbuild's toolchain check first in this process
stopifnot(pkgbuild::has_build_tools())
expose_stan_functions(stanc(model_code = fun), env = env)
# frmtmb is attached only after the compile, which failed with it attached
suppressMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("frmtmb from", find.package("frmtmb"), "\n")
cat("brms", format(packageVersion("brms")), "rstan",
    format(packageVersion("rstan")), "evd", format(packageVersion("evd")),
    "Rmpfr", format(packageVersion("Rmpfr")), "\n")

PB <- 256
mp <- function(x) mpfr(x, PB)
err <- function(a, b) abs(a - b) / pmax(1, abs(b))
# a one-line summary of a against the reference b: the worst error over
# the points where both are finite, and where they disagree on finiteness
cmp <- function(label, a, b) {
  ok <- is.finite(a) & is.finite(b)
  bad <- xor(is.finite(a), is.finite(b))
  cat(sprintf("  %-42s n = %3d  max err = %.2e  finite in one only: %d\n",
              label, sum(ok), if (any(ok)) max(err(a[ok], b[ok])) else NA,
              sum(bad)))
  invisible(which(bad))
}

on_tape <- function(fam, y, eta, what = "lpdf") {
  dp <- names(eta)
  len <- vapply(eta, length, 1L)
  f <- function(p) {
    pars <- list()
    at <- 0L
    for (nm in dp) {
      e <- p[at + seq_len(len[[nm]])]
      at <- at + len[[nm]]
      pars[[nm]] <- fam$links[[nm]]$linkinv(e)
      pars[[paste0(".eta_", nm)]] <- e
    }
    fam[[what]](y, pars, list())
  }
  tp <- RTMB::MakeTape(f, unlist(eta, use.names = FALSE))
  out <- tp(unlist(eta, use.names = FALSE))
  if (what == "lcdf") log(out) else out
}
off_tape <- function(fam, y, eta, what = "lpdf") {
  pars <- lapply(names(eta), function(nm) fam$links[[nm]]$linkinv(eta[[nm]]))
  names(pars) <- names(eta)
  out <- fam[[what]](y, pars, list())
  if (what == "lcdf") log(out) else out
}

# --------------------------------------------------------------- frechet
cat("\n=== frechet: Stan's frechet(nu, mu / gamma(1 - 1/nu))\n")
g <- expand.grid(y = c(1e-3, 0.05, 0.5, 1, 2.5, 10, 100, 1e4),
                 mu = c(0.5, 1, 4, 50), nu = c(1.05, 1.5, 3, 12, 60))
cat("grid:", nrow(g), "points; y in", range(g$y), "mu in", range(g$mu),
    "nu in", range(g$nu), "\n")
fam <- frechet()
eta <- list(mu = log(g$mu), nu = log(g$nu - 1))
nu_d <- 1 + exp(eta$nu)
s <- g$mu / gamma(1 - 1 / nu_d)
# reference in 256 bits, from the same double scale Stan is handed
zm <- log(mp(g$y)) - log(mp(s))
u <- exp(-mp(nu_d) * zm)
ref <- list(lpdf = as.numeric(log(mp(nu_d)) - log(mp(s)) -
                                (1 + mp(nu_d)) * zm - u),
            lcdf = as.numeric(-u),
            lccdf = as.numeric(log(-expm1(-u))))
for (w in names(ref)) {
  cat(w, "\n")
  stan <- mapply(get(paste0("fr_", w), env), g$y, nu_d, s)
  cmp("frmtmb off tape vs 256-bit", off_tape(fam, g$y, eta, w), ref[[w]])
  cmp("frmtmb on tape vs 256-bit", on_tape(fam, g$y, eta, w), ref[[w]])
  cmp("brms/Stan vs 256-bit", stan, ref[[w]])
  cmp("frmtmb on tape vs brms/Stan", on_tape(fam, g$y, eta, w), stan)
}
ev <- log(evd::dfrechet(g$y, loc = 0, scale = s, shape = nu_d))
cat("lpdf\n")
cmp("evd::dfrechet vs 256-bit", ev, ref$lpdf)
cmp("frmtmb off tape vs evd::dfrechet", off_tape(fam, g$y, eta), ev)

# ---------------------------------------------- generalized extreme value
cat("\n=== gen_extreme_value: brms's Stan functions\n")
g <- expand.grid(y = c(-30, -3, -0.5, 0, 0.7, 2, 8, 50),
                 mu = c(-1, 0.5), sigma = c(0.3, 1, 4),
                 xi = c(-0.9, -0.4, -1e-7, 0, 1e-9, 1e-4, 0.01, 0.3, 1.2))
fam <- gen_extreme_value()
eta <- list(mu = g$mu, sigma = log(g$sigma), xi = log1p(g$xi))
# every reference is handed the xi and sigma the tape sees, which are
# expm1(log1p(xi)) and exp(log(sigma)), not always the grid values to
# the last bit
xi_d <- expm1(eta$xi)
sg_d <- exp(eta$sigma)
inside <- 1 + xi_d * (g$y - g$mu) / sg_d > 0
cat("grid:", nrow(g), "points,", sum(inside), "inside the support; xi in",
    paste(unique(g$xi), collapse = ", "), "\n")
xm <- (mp(g$y) - mp(g$mu)) / mp(sg_d)
xim <- mp(xi_d)
gum <- xi_d == 0
lt <- log1p(xim * xm)
hm <- lt / xim
hm[gum] <- xm[gum]
lt[gum] <- mp(0)
ref <- list(lpdf = as.numeric(-log(mp(sg_d)) - lt - hm - exp(-hm)),
            lcdf = as.numeric(-exp(-hm)),
            lccdf = as.numeric(log(-expm1(-exp(-hm)))))
for (w in names(ref)) {
  cat(w, "(inside the support)\n")
  stan <- mapply(get(paste0("gev_", w, "_w"), env), g$y, g$mu, sg_d, xi_d)
  ft <- on_tape(fam, g$y, eta, w)
  cmp("frmtmb off tape vs 256-bit", off_tape(fam, g$y, eta, w)[inside],
      ref[[w]][inside])
  cmp("frmtmb on tape vs 256-bit", ft[inside], ref[[w]][inside])
  cmp("brms/Stan vs 256-bit", stan[inside], ref[[w]][inside])
  big <- inside & abs(xi_d) >= 1e-4
  cmp("brms/Stan vs 256-bit, |xi| >= 1e-4", stan[big], ref[[w]][big])
  cmp("frmtmb on tape vs brms/Stan, |xi| >= 1e-4", ft[big], stan[big])
}
cat("outside the support (", sum(!inside), "points):\n")
lp_out <- on_tape(fam, g$y[!inside], lapply(eta, `[`, !inside))
t_out <- (1 + xi_d * (g$y - g$mu) / sg_d)[!inside]
cat("  density exp(lpdf) is exactly 0 at", sum(exp(lp_out) == 0), "of",
    length(lp_out), "points; the largest log density there is",
    format(max(lp_out), digits = 4), "at t =",
    format(t_out[which.max(lp_out)], digits = 4),
    "(brms/Stan: NaN, which rejects the value)
")
cf_out <- exp(on_tape(fam, g$y[!inside], lapply(eta, `[`, !inside), "lcdf"))
want <- ifelse(xi_d[!inside] > 0, 0, 1)
cat("  CDF: max |F - (0 below the lower end, 1 above the upper end)| =",
    format(max(abs(cf_out - want)), digits = 3), "\n")
ev <- mapply(function(y, m, s, k) {
  evd::dgev(y, loc = m, scale = s, shape = k, log = TRUE)
}, g$y, g$mu, sg_d, xi_d)
cat("lpdf\n")
cmp("evd::dgev vs 256-bit (inside)", ev[inside], ref$lpdf[inside])
cmp("frmtmb off tape vs evd::dgev (inside)",
    off_tape(fam, g$y, eta)[inside], ev[inside])

cat("\n--- d lpdf / d eta_xi across the series/quotient blend:\n")
cat("    tape gradient against a 256-bit central difference (h = 1e-30)\n")
xs <- c(-0.05, -0.016, -0.012, -0.01, -0.005, -1e-3, -1e-6, -1e-12, 0,
        1e-12, 1e-6, 1e-3, 0.005, 0.01, 0.012, 0.016, 0.05)
ref_lpdf <- function(yv, m, s, e) {
  x <- (mp(yv) - mp(m)) / mp(s)
  k <- expm1(e)
  l <- log1p(k * x)
  -log(mp(s)) - l - l / k - exp(-l / k)
}
tape_at <- function(yv, x0, order) {
  f <- function(p) {
    fam$lpdf(yv, list(mu = 0.2, sigma = 1.3, xi = expm1(p), .eta_xi = p),
             list())
  }
  tp <- RTMB::MakeTape(f, log1p(x0))
  if (order == 1) tp$jacobian(log1p(x0))[1, 1] else
    tp$jacfun()$jacobian(log1p(x0))[1, 1]
}
for (yv in c(-1.5, 0.4, 3)) {
  gt <- vapply(xs, function(x0) tape_at(yv, x0, 1), 0)
  h <- mp(1e-30)
  # at e = 0 the central difference straddles the Gumbel point, where
  # the 256-bit quotient is finite on both sides
  gn <- vapply(xs, function(x0) {
    e <- mp(log1p(x0))
    as.numeric((ref_lpdf(yv, 0.2, 1.3, e + h) -
                  ref_lpdf(yv, 0.2, 1.3, e - h)) / (2 * h))
  }, 0)
  cat(sprintf("  y = %4.1f: max err over %d xi = %.2e\n", yv, length(xs),
              max(err(gt, gn))))
}
cat("--- second derivative (tape Hessian vs 256-bit second difference)\n")
for (yv in c(-1.5, 0.4, 3)) {
  ht <- vapply(xs, function(x0) tape_at(yv, x0, 2), 0)
  h <- mp(1e-25)
  hn <- vapply(xs, function(x0) {
    e <- mp(log1p(x0))
    as.numeric((ref_lpdf(yv, 0.2, 1.3, e + h) -
                  2 * ref_lpdf(yv, 0.2, 1.3, e) +
                  ref_lpdf(yv, 0.2, 1.3, e - h)) / h^2)
  }, 0)
  ok <- xs != 0
  cat(sprintf("  y = %4.1f: max err over %d xi (xi = 0 left out) = %.2e\n",
              yv, sum(ok), max(err(ht[ok], hn[ok]))))
}

# ------------------------------------------------------- discrete weibull
cat("\n=== discrete_weibull: brms's Stan functions\n")
g <- expand.grid(y = c(0, 1, 2, 5, 20, 200),
                 mu = c(1e-6, 0.05, 0.5, 0.9, 0.999),
                 k = c(0.2, 0.7, 1, 2.5, 6))
cat("grid:", nrow(g), "points\n")
fam <- discrete_weibull()
eta <- list(mu = qlogis(g$mu), shape = log(g$k))
mu_d <- plogis(eta$mu)
k_d <- exp(eta$shape)
Lm <- log(mp(mu_d))
a <- Lm * mp(g$y)^mp(k_d)
b <- Lm * (mp(g$y) + 1)^mp(k_d)
# log(exp(a) - exp(b)) written as a + log1p(-exp(b - a)), an identity: at
# shape 6 and y = 200 exp(a) is below even 256-bit MPFR's exponent range
# (it underflowed to -Inf at 6 of the 150 points written the other way)
ref <- list(lpdf = as.numeric(a + log1p(-exp(b - a))),
            lcdf = as.numeric(log(-expm1(b))),
            lccdf = as.numeric(b))
for (w in names(ref)) {
  cat(w, "\n")
  sf <- get(paste0("dw_", if (w == "lpdf") "lpmf" else w, "_w"), env)
  stan <- mapply(sf, as.integer(g$y), mu_d, k_d)
  cmp("frmtmb off tape vs 256-bit", off_tape(fam, g$y, eta, w), ref[[w]])
  cmp("frmtmb on tape vs 256-bit", on_tape(fam, g$y, eta, w), ref[[w]])
  bad <- cmp("brms/Stan vs 256-bit", stan, ref[[w]])
  if (length(bad)) {
    cat("    where brms/Stan is not finite and the reference is:",
        paste(sprintf("(y=%g, mu=%g, k=%g: %.4g)", g$y[bad], g$mu[bad],
                      g$k[bad], ref[[w]][bad])[seq_len(min(4, length(bad)))],
              collapse = " "), "\n")
  }
  cmp("frmtmb on tape vs brms/Stan", on_tape(fam, g$y, eta, w), stan)
}
cat("on the tape at eta_mu = 40 (mu rounds to 1 in double), shape 1:\n")
ft <- on_tape(fam, c(0, 3, 50), list(mu = rep(40, 3), shape = rep(0, 3)))
L40 <- -log1p(exp(mp(-40)))
rf <- as.numeric(c(log(-expm1(L40)), 3 * L40 + log(-expm1(L40)),
                   50 * L40 + log(-expm1(L40))))
cmp("frmtmb on tape vs 256-bit", ft, rf)
cat("  brms/Stan at mu = plogis(40):",
    format(mapply(env$dw_lpmf_w, c(0L, 3L, 50L), plogis(40), 1),
           digits = 4), "\n")
s1 <- vapply(c(0.3, 0.9), function(m) {
  sum(exp(off_tape(fam, 0:20000, list(mu = rep(qlogis(m), 20001),
                                      shape = rep(log(0.8), 20001)))))
}, 0)
cat("sum of the pmf over 0..20000 at shape 0.8, mu = 0.3 and 0.9:",
    format(s1, digits = 15), "\n")
