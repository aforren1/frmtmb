# The extreme-value and discrete Weibull families: frechet(),
# gen_extreme_value() and discrete_weibull(), brms 2.23.0's three
# univariate families that frm() did not take before. They live in a
# file of their own so that families.R, which other work also edits,
# carries only their registry lines.

#' Extreme-value and discrete Weibull families
#'
#' `frechet()`, `gen_extreme_value()` and `discrete_weibull()` are
#' brms's families of the same names, with brms's parameters, links and
#' defaults.
#'
#' @section Frechet:
#' `frechet()` is for a strictly positive response with a heavy right
#' tail. Its parameters are `mu` and `nu`. As in brms, `mu` is the MEAN,
#' not the scale: the density is Stan's `frechet(nu, s)` with scale
#' `s = mu / gamma(1 - 1 / nu)`, which needs `nu > 1`. The default link
#' of `nu` is `logm1`, `log(nu - 1)`, which keeps it there. Under
#' `link_nu = "identity"` a `nu` at or below 1 has no mean, and the
#' density is `NaN` there, as brms's is, in the fitted objective as well
#' as on doubles.
#' The variance is finite only for `nu > 2`; `residuals(type =
#' "pearson")` divides by `Inf` where it is not.
#'
#' @section Generalized extreme value:
#' `gen_extreme_value()` has a location `mu`, a scale `sigma` and a
#' shape `xi`. The density is
#' `exp(-t^(-1/xi)) t^(-1 - 1/xi) / sigma` with
#' `t = 1 + xi (y - mu) / sigma`, and its limit at `xi = 0` is the
#' Gumbel density, which the family evaluates exactly at `xi = 0` and
#' to full precision near it. So `bf(y ~ x, xi = 0)` fits a Gumbel
#' model, as in brms.
#'
#' The support depends on the parameters: `y > mu - sigma / xi` for
#' `xi > 0` and `y < mu - sigma / xi` for `xi < 0`. Outside the
#' support the density is 0. Its log falls by `1e8` per unit of
#' `t` below 0 rather than jumping to `-Inf`, so a row a millionth
#' outside has a log density below -100 and the optimizer meets a slope
#' that turns it back, not a wall. The start, `xi = 0`, is inside the
#' support for every data set. For `xi > -1` the density goes to 0 at
#' an end of the support, so a maximum has every row inside it. That is
#' the region brms's `scale_xi()` maps its `tmp_xi` onto.
#'
#' In a [mixture()] a component that does not reach a row adds nothing
#' to that row's density, which is what a mixture density is. brms's
#' density is `NaN` there, which rejects the parameter value, so a brms
#' mixture holds each component's support over every row.
#'
#' The default link of `xi` is `log1p`, `log(1 + xi)`, brms's link
#' for a predicted `xi`, which keeps `xi > -1`. `link_xi = "identity"`
#' is the other choice brms allows. Under it, `xi <= -1` has no
#' likelihood, so every link holds `xi` above -1. Below -1 the
#' likelihood is unbounded and has no maximum.
#'
#' brms parameterizes an `xi` without a formula differently. It samples
#' an unbounded `tmp_xi` and scales it into the support the current `mu`
#' and `sigma` leave, and it reports `xi` from it after the fit. That is
#' a change of variable, so the maximum of the likelihood is the same
#' point. A fit does not use it, because its range reaches below
#' `xi = -1`. `frmtmb.sample::frm_sample()` does, because there the
#' support's edge is a wall the sampler meets; it holds the range above
#' -1 (see [frmtmb-gev-xi-chart]).
#'
#' A fit whose `xi` stops within 1e-3 of -1 with the likelihood not
#' stationary in `xi` there warns: its log-likelihood is a supremum, not a
#' maximum. A maximum inside that band gets a message. For
#' `-1 < xi < -0.5` the maximum exists but is not regular (Smith, 1985),
#' so its standard errors do not have their usual meaning; outside the
#' band no message says so.
#'
#' The mean, which [fitted()] returns, is
#' `mu + sigma (gamma(1 - xi) - 1) / xi`, and `mu + sigma * 0.5772...`
#' (Euler's constant) at `xi = 0`. It is infinite for `xi >= 1`, and
#' `fitted()` returns `Inf` there. brms's `posterior_epred()` uses the
#' same formula and returns `NaN` at `xi = 0` and `xi = 1` and a finite,
#' wrong number for `xi > 1`. The variance is finite only for
#' `xi < 1 / 2`.
#'
#' brms 2.23.0 marks `gen_extreme_value` deprecated and warns when it
#' builds the Stan code. frmtmb does not repeat the warning: it is about
#' brms's plans for its own code, not about the model.
#'
#' @section Discrete Weibull:
#' `discrete_weibull()` is for a count `0, 1, 2, ...`. Its parameters
#' are `mu` in `(0, 1)` and a positive `shape`, with
#' `P(Y >= y) = mu^(y^shape)`, so `P(Y = y) = mu^(y^shape) -
#' mu^((y + 1)^shape)`. `mu` is not the mean: it is the probability of a
#' count of at least 1. The mean, which [fitted()] returns, is the sum
#' of `mu^(y^shape)` over `y >= 1`. frmtmb sums it to a relative
#' precision of about 1e-12, the tail beyond 1000 terms through its
#' integral. brms's `posterior_epred()` stops the sum at the first term
#' below 0.001, so the two differ where the tail is long.
#'
#' `cens()` and `trunc()` follow brms's discrete conventions (see
#' [frmtmb_family()], Censoring a discrete response): a right-censored
#' count at `k` is `P(Y > k)`, and `trunc(lb = )` reads its bound as
#' inclusive.
#'
#' @section Addition terms:
#' All three take `weights()`, `cens()` and `trunc()`. brms also takes
#' `mi()` on a `frechet` or `gen_extreme_value` response; frmtmb takes
#' `mi()` on a gaussian or student response only.
#'
#' @references
#' Smith, R. L. (1985). Maximum likelihood estimation in a class of
#' nonregular cases. *Biometrika*, 72(1), 67-90.
#'
#' @inheritParams frmtmb-families
#' @param link_nu Link for `nu`: `"logm1"` (the default) or `"identity"`.
#' @param link_xi Link for `xi`: `"log1p"` (the default) or `"identity"`.
#' @param link_sigma Link for the GEV's `sigma`: `"log"` (the default),
#'   `"identity"`, `"softplus"` or `"squareplus"`.
#' @param link_shape Link for the discrete Weibull's `shape`: `"log"` (the
#'   default), `"identity"`, `"softplus"` or `"squareplus"`.
#' @return A `frmtmb_family`.
#' @examples
#' set.seed(1)
#' n <- 200
#' dd <- data.frame(x = rnorm(n))
#'
#' # Frechet: mu is the mean
#' nu <- 4
#' mu <- exp(1 + 0.3 * dd$x)
#' dd$yf <- mu / gamma(1 - 1 / nu) * rexp(n)^(-1 / nu)
#' fit_f <- frm(yf ~ x, family = frechet(), data = dd)
#' fixef(fit_f)
#'
#' # generalized extreme value
#' xi <- 0.2
#' dd$yg <- 1 + 0.5 * dd$x + 0.8 * (rexp(n)^(-xi) - 1) / xi
#' fit_g <- frm(bf(yg ~ x, sigma ~ 1), family = gen_extreme_value(),
#'              data = dd)
#' fixef(fit_g)
#'
#' # discrete Weibull counts
#' q <- plogis(1 + 0.4 * dd$x)
#' dd$yd <- ceiling((log1p(-runif(n)) / log(q))^(1 / 1.3) - 1)
#' fit_d <- frm(yd ~ x, family = discrete_weibull(), data = dd)
#' head(fitted(fit_d))
#' @name frmtmb-extreme-families
NULL

#' @rdname frmtmb-extreme-families
#' @export
frechet <- function(link = "log", link_nu = "logm1") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["frechet"]], "log")
  fam_frechet(link, link_nu)
}

#' @rdname frmtmb-extreme-families
#' @export
gen_extreme_value <- function(link = "identity", link_sigma = "log",
                              link_xi = "log1p") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["gen_extreme_value"]], "identity")
  fam_gen_extreme_value(link, link_sigma, link_xi)
}

#' @rdname frmtmb-extreme-families
#' @export
discrete_weibull <- function(link = "logit", link_shape = "log") {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["discrete_weibull"]], "logit")
  fam_discrete_weibull(link, link_shape)
}

# brms:::links_dpars("xi"): the shape of the GEV must stay above -1 for
# its likelihood to have a maximum, which log1p gives.
dpar_links_xi <- c("log1p", "identity")

#' `exp(-w)`, held at `exp(700)` where it would overflow.
#'
#' The Frechet and GEV densities subtract `exp(-w)` and their CDFs are
#' `exp(-exp(-w))`. Far into the lower tail `exp(-w)` overflows, the log
#' density is `-Inf` and its gradient `Inf * 0`, a `NaN`, which in a
#' `mixture()` reaches the whole row through the log-sum, where that
#' component's weight is 0. Past `w = -700` the density is below
#' `exp(-1e304)` either way, so holding it there changes no value a
#' double can carry.
#'
#' @noRd
exp_neg_capped <- function(w) exp(-ad_max(w, -700))

# ---------------------------------------------------------------- frechet

#' 0 where `x > 0` and `NaN` where `x <= 0`, on the RTMB tape as on
#' doubles: a term that takes a parameter out of the likelihood where
#' the density does not exist. Its derivatives are 0 where it is 0, with
#' three exceptions at the ends of the doubles (the review's round 2,
#' dev/famuni-rev2-out/guard.txt and hess.txt): the gradient is NaN for
#' `0 < x < 2.3e-308`, where `x` is subnormal; the Hessian is NaN for
#' `0 < x <= 1e-155`; and the value itself is NaN for `x >= 9e307`, where
#' `u` overflows. The GEV and the Frechet cannot reach these: `1 + xi`
#' and `nu - 1` are 0 or at least 1.1e-16 in double precision. frmtmb.eam's
#' wiener guard reaches them only with a parameter below 1e-154.
#'
#' With `u = x + |x|` it is `u / u - 1`. It must not have a constant 0
#' as an operand: RTMB folds `0 * x` and `0 / x` to the constant 0 when
#' it records the tape, so `0 * log(x)` is NaN on doubles and 0 in every
#' fit (dev/rtmb-pitfalls.md, item 22). It raises no warning on doubles,
#' where `log(x) - log(x)` would.
#'
#' @noRd
nan_unless_pos <- function(x) {
  u <- x + abs(x)
  u / u - 1
}

#' The log of the Frechet scale, `log(mu) - lgamma(1 - 1 / nu)`, and a
#' term that is `NaN` where `nu <= 1` under a link that allows it. There
#' `gamma(1 - 1 / nu)` is negative or infinite, and `lgamma()` would
#' return the log of its absolute value, a finite and wrong scale.
#'
#' @noRd
frechet_log_scale <- function(dpars, lk, lk_nu) {
  nu <- dpars[["nu"]]
  out <- dpar_log(dpars, "mu", lk) - lgamma(1 - 1 / nu)
  if (identical(lk_nu$name, "logm1")) return(out)
  out + nan_unless_pos(nu - 1)
}

#' Frechet family, dpars `mu` (the mean) and `nu` (the shape, above 1).
#'
#' @noRd
fam_frechet <- function(link = "log", link_nu = "logm1") {
  lk <- mu_link(link, "frechet")
  lk_nu <- dpar_link(link_nu, "nu", "frechet", dpar_links_above1)
  # z = log(y / scale): the density, the CDF and the log survivor are
  # each a function of nu * z, which keeps them finite far into both
  # tails where (y / scale)^(-nu) itself would overflow or underflow
  logz <- function(q, dpars) log(q) - frechet_log_scale(dpars, lk, lk_nu)
  frmtmb_family(
    "frechet",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = c("mu", "nu"),
    links = list(mu = lk, nu = lk_nu),
    lpdf = function(y, dpars, aterms) {
      ov <- osa_unwrap(y)
      nu <- dpars[["nu"]]
      ls <- frechet_log_scale(dpars, lk, lk_nu)
      z <- log(if (is.null(ov)) y else ov$y) - ls
      ll <- log(nu) - ls - (1 + nu) * z - exp_neg_capped(nu * z)
      if (is.null(ov)) ll else ll * ov$keep
    },
    lcdf = function(q, dpars, aterms) {
      exp(-exp_neg_capped(dpars[["nu"]] * logz(q, dpars)))
    },
    lccdf = function(q, dpars, aterms) {
      log(-expm1(-exp_neg_capped(dpars[["nu"]] * logz(q, dpars))))
    },
    valid_y = positive_y("frechet"),
    init_dpars = list(
      # at nu = 3 the mean is 1.2 times the median; the median is the
      # start that a heavy right tail cannot drag away
      mu = function(y, aterms) 1.2 * stats::median(y),
      nu = function(y, aterms) 3
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) {
        nu <- dpars[["nu"]]
        g1 <- gamma(1 - 1 / nu)
        g2 <- suppressWarnings(gamma(1 - 2 / nu))
        sc <- dpars[["mu"]] / g1
        ifelse(nu > 2, sc^2 * (g2 - g1^2), Inf)
      },
      # Y = scale * E^(-1/nu) for E ~ Exp(1), so the partial mean over
      # (lb, ub) is scale * gamma(1 - 1/nu) * (P(a, e_lb) - P(a, e_ub))
      # with a = 1 - 1/nu, e_q = (q / scale)^(-nu), and P the
      # regularized lower incomplete gamma; scale * gamma(a) is mu
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        nu <- dpars[["nu"]]
        sc <- dpars[["mu"]] / gamma(1 - 1 / nu)
        e_lb <- (pmax(lb, 0) / sc)^(-nu)
        e_ub <- (ub / sc)^(-nu)
        a <- 1 - 1 / nu
        dpars[["mu"]] * (stats::pgamma(e_lb, a) - stats::pgamma(e_ub, a)) /
          (exp(-e_ub) - exp(-e_lb))
      }
    ),
    sim = function(dpars, aterms, n) {
      nu <- dpars[["nu"]]
      dpars[["mu"]] / gamma(1 - 1 / nu) * stats::rexp(n)^(-1 / nu)
    }
  )
}

# -------------------------------------------------- generalized extreme value

#' `log1p(z) / z`, finite at `z = 0` with all the derivatives a fit
#' takes of it.
#'
#' The GEV density divides `log(1 + xi x)` by `xi`, which is `0 / 0` at
#' the Gumbel point `xi = 0`, and a tape cannot branch on `xi` there.
#' So the ratio is a blend: its Taylor series where `|z| < 0.016`, the
#' quotient itself where `|z| > 0.016`, joined by the weight
#' `b = exp(-(z / 0.01)^8)`. The series to `z^11` is exact to 1e-19
#' where `b` is not negligible, and `b` is below 1e-17 wherever the
#' series is not exact. The quotient is evaluated at `z + 1e-50 b`,
#' which is `z` to the last bit wherever the quotient's weight `1 - b`
#' is not negligible, and which keeps it and its derivatives finite at
#' `z = 0`, where `1 - b` is exactly 0 and a `NaN` would survive the
#' product. Both clamps to `[-1, 1]` keep the polynomial and the weight
#' from overflowing at a far `z`, where `b` is exactly 0.
#'
#' @noRd
log1p_ratio <- function(z) {
  zc <- ad_min(ad_max(z, -1), 1)
  u <- (zc / 0.01)^2
  b <- exp(-(u * u)^2)
  s <- 0
  for (k in 11:0) s <- (-1)^k / (k + 1) + zc * s
  ze <- z + 1e-50 * b
  b * s + (1 - b) * (log1p(ze) / ze)
}

#' The pieces of the GEV density shared by its density, CDF and log
#' survivor: `l = log(t)` and `h = log(t) / xi`, with `t = 1 + xi x`,
#' so that `exp(-h)` is `t^(-1/xi)`, and `out`, how far a row lies
#' outside the support, where `t <= 0`.
#'
#' Outside the support the density is 0 and the distribution function 0
#' below the lower end (`xi > 0`) or 1 above the upper end (`xi < 0`).
#' `t` is held at 1e-15 there, which keeps `h` finite with the sign of
#' `x` and at least the edge's size, so `exp(-exp(-h))` is that 0 or 1
#' to within 1e-13 at any `xi` up to 10, as a `trunc()` bound or a
#' censoring point beyond the support needs.
#'
#' @noRd
gev_parts <- function(q, dpars) {
  x <- (q - dpars[["mu"]]) / dpars[["sigma"]]
  z <- dpars[["xi"]] * x
  # how far outside, taken on its own: written as zc - z it is 0 to the
  # bit, but the density scales it by 1e8, and the tape then sums -1e8
  # and +1e8 into the gradient of z beside its O(1) terms, which cost
  # those terms eight digits (dev/famuni-gradprobe3.R)
  out <- ad_pos(-1 + 1e-15 - z)
  # z + out rounds to -1 itself once |z| is large (at z = -27 the 1e-15
  # is below the spacing of 26), and log1p(-1) is -Inf: the second clamp
  # starts at -1 and is exact there (dev/famuni-cdfprobe.R)
  zc <- ad_max(z + out, -1 + 1e-15)
  list(l = log1p(zc), h = x * log1p_ratio(zc), out = out)
}

# zeta(2), ..., zeta(22), from 256-bit Rmpfr::zeta(): the coefficients
# of log gamma(1 - x) = gamma_E x + sum_{k >= 2} zeta(k) x^k / k
gev_zeta <- c(1.6449340668482264, 1.2020569031595942, 1.0823232337111381,
  1.03692775514337, 1.0173430619844492, 1.0083492773819229,
  1.0040773561979444, 1.0020083928260821, 1.000994575127818,
  1.0004941886041194, 1.000246086553308, 1.0001227133475785,
  1.0000612481350588, 1.000030588236307, 1.0000152822594086,
  1.0000076371976379, 1.000003817293265, 1.0000019082127165,
  1.0000009539620338, 1.0000004769329869, 1.0000002384505027)

#' The two ratios the GEV mean and variance are made of, for
#' `|xi| <= 0.05`, from the series of `log gamma(1 - x)`:
#' `a = log gamma(1 - xi) / xi` and
#' `d = (log gamma(1 - 2 xi) - 2 log gamma(1 - xi)) / xi^2`, in which
#' the terms that cancel have been cancelled by hand.
#'
#' `gamma()` is exact to a few eps relative to 1 near 1, so
#' `gamma(1 - xi) - 1` keeps only `eps / xi` of its digits and the
#' variance's `gamma(1 - 2 xi) - gamma(1 - xi)^2` only `eps / xi^2`: 8.8e-5
#' relative at `xi = 1e-6` (dev/famuni-gevmean-mp.R). The series to
#' `k = 22` is exact to 1e-19 in this range.
#'
#' @noRd
gev_lgamma_ratios <- function(xi) {
  k <- seq_along(gev_zeta) + 1L
  a <- -digamma(1) + as.vector(outer(xi, k - 1L, `^`) %*% (gev_zeta / k))
  d <- as.vector(outer(xi, k - 2L, `^`) %*% (gev_zeta * (2^k - 2) / k))
  list(a = a, d = d)
}

#' `expm1(u) / u`, 1 at `u = 0`.
#'
#' @noRd
expm1_ratio <- function(u) ifelse(u == 0, 1, expm1(u) / ifelse(u == 0, 1, u))

#' The GEV mean `mu + sigma (gamma(1 - xi) - 1) / xi`, off the tape;
#' `mu + sigma * 0.5772...` (Euler's constant) at `xi = 0` and `Inf` for
#' `xi >= 1`, where the mean does not exist.
#'
#' @noRd
gev_mean <- function(mu, sigma, xi) {
  small <- abs(xi) <= 0.05
  xs <- ifelse(small | xi >= 1, 0.5, xi)
  s <- gev_lgamma_ratios(ifelse(small, xi, 0))
  # (gamma(1 - xi) - 1) / xi = a * expm1(a xi) / (a xi)
  r <- ifelse(small, s$a * expm1_ratio(s$a * xi),
              expm1(lgamma(1 - xs)) / xs)
  ifelse(xi >= 1, Inf, mu + sigma * r)
}

#' The GEV variance, `sigma^2 (g2 - g1^2) / xi^2` with
#' `g_k = gamma(1 - k xi)`, off the tape; `sigma^2 pi^2 / 6` at `xi = 0`
#' and `Inf` for `xi >= 1/2`.
#'
#' @noRd
gev_var <- function(sigma, xi) {
  small <- abs(xi) <= 0.05
  xs <- ifelse(small | xi >= 0.5, 0.25, xi)
  s <- gev_lgamma_ratios(ifelse(small, xi, 0))
  # g2 - g1^2 = g1^2 expm1(D), D = d xi^2, and g1^2 = exp(2 a xi)
  gs <- exp(2 * s$a * xi) * s$d * expm1_ratio(s$d * xi^2)
  g <- exp(2 * lgamma(1 - xs)) *
    expm1(lgamma(1 - 2 * xs) - 2 * lgamma(1 - xs)) / xs^2
  ifelse(xi >= 0.5, Inf, sigma^2 * ifelse(small, gs, g))
}

#' The GEV partial mean over `(lb, ub)`, divided by its mass, off the
#' tape, by adaptive quadrature over the quantile scale: the integral of
#' `Q(p)` over `(F(lb), F(ub))`. Quadrature rather than the incomplete
#' gamma form because that form divides by `xi` and needs a positive
#' shape `1 - xi`, so it fails at the Gumbel point and for `xi >= 1`,
#' where a truncated mean still exists.
#'
#' @noRd
gev_trunc_mean <- function(mu, sigma, xi, lb, ub) {
  n <- max(length(mu), length(sigma), length(xi), length(lb), length(ub))
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)
  lb <- rep_len(lb, n)
  ub <- rep_len(ub, n)
  cdf <- function(q, m, s, k) {
    x <- (q - m) / s
    t <- pmax(1 + k * x, 0)
    if (k == 0) exp(-exp(-x)) else exp(-t^(-1 / k))
  }
  qf <- function(p, m, s, k) {
    e <- -log(p)
    if (k == 0) m - s * log(e) else m + s * expm1(-k * log(e)) / k
  }
  vapply(seq_len(n), function(i) {
    pl <- cdf(lb[i], mu[i], sigma[i], xi[i])
    pu <- cdf(ub[i], mu[i], sigma[i], xi[i])
    if (!(pu > pl)) return(NaN)
    v <- tryCatch(stats::integrate(qf, pl, pu, m = mu[i], s = sigma[i],
                                   k = xi[i], rel.tol = 1e-10)$value,
                  error = function(e) NaN)
    v / (pu - pl)
  }, numeric(1))
}

#' Draws from the GEV by its quantile function at a uniform.
#'
#' @noRd
rgev <- function(n, mu, sigma, xi) {
  e <- stats::rexp(n)
  xi <- rep_len(xi, n)
  mu + sigma * ifelse(xi == 0, -log(e),
                      expm1(-xi * log(e)) / ifelse(xi == 0, 1, xi))
}

#' Warn when a fitted `xi` has run to -1, its edge.
#'
#' For `xi <= -1` the GEV likelihood has no maximum, and both links stop
#' the optimizer above -1: `log1p` by construction, `identity` by the
#' `NaN` it gives there. A fit that stops on the edge is the supremum of
#' a likelihood that rises toward it, not an interior maximum, and the
#' optimizer's own warning there, "false convergence (8)", does not say
#' which parameter or why. Measured on the review's cases
#' (dev/reviews/2026-10-13-famuni.md, fix 3 and fix 4): seeds 621 and 622
#' under the identity link end at xi = -1.0000, and seed 1402 under
#' `log1p` at -0.9999998.
#'
#' A fitted `xi` within 1e-3 of -1 is not by itself an edge: the review's
#' 720 correct-model fits held one true maximum at -0.999092 (identity
#' link, n = 1000, seed 39), and the first version of this check called
#' it a supremum. So a scalar `xi` is judged by `d logLik / d xi` at the
#' stop, from the objective's gradient carried to the `xi` scale. Where
#' it is not near 0 the edge, not a maximum, stopped the fit, and the
#' warning says so: negative, the likelihood still rises toward -1;
#' large and positive, the fit sits in the corner where moving `xi` up
#' would push the largest row out of the support, which the identity
#' link reaches (3 of the 720 fits, at xi + 1 below 1e-15, d logLik /
#' d xi 154 to 238). Near 0 the fit is at a maximum inside the band and
#' gets a message instead. On the 720 fits (dev/famuni-r2-xiedge.R) the
#' edge fits have `|d logLik / d xi|` of 0.88 or more and the one
#' interior maximum 3.3e-4; 0.1 separates them with room on both sides.
#' A predicted `xi` has no one direction to read, and gets a warning
#' that names the distance only.
#'
#' @noRd
gev_fit_check <- function(fit, resp) {
  multi <- length(fit$spec$responses) > 1L
  xi <- tryCatch(
    frm_linpred(fit, dpar = "xi", type = "response",
                resp = if (multi) resp),
    error = function(e) NULL)
  if (!is.numeric(xi) || !length(xi) || !any(xi < -1 + 1e-3, na.rm = TRUE)) {
    return(invisible(NULL))
  }
  of <- if (multi) paste0("of response ", resp, " ")
  at <- format(min(xi, na.rm = TRUE), digits = 7)
  lp <- fit$frame[["linpreds"]][[linpred_key(resp, "xi")]]
  if (!is.null(lp[["constant"]])) return(invisible(NULL))
  dll <- if (!dpar_has_predictor(fit$spec, lp)) gev_xi_slope(fit, lp)
  if (is.null(dll) || !is.finite(dll)) {
    frm_warning("gen_extreme_value: xi ", of, "reaches ", at, ", within ",
                "1e-3 of the edge at -1. Below -1 the likelihood has no ",
                "maximum, and near it the standard errors do not have ",
                "their usual meaning. Check whether the fit stopped on ",
                "that edge", call. = FALSE, class = "frmtmb_gev_xi_edge")
  } else if (abs(dll) > 0.1) {
    frm_warning("gen_extreme_value: xi ", of, "ran to ", at, ", the edge ",
                "at -1, and the likelihood is not stationary in xi there ",
                "(d logLik / d xi = ", format(dll, digits = 3), "): the ",
                "edge, not a maximum, stopped the fit. Below -1 ",
                "the likelihood has no maximum; on the edge the ",
                "log-likelihood is a supremum and the standard errors do ",
                "not have their usual meaning. A ",
                "shorter upper tail than the GEV allows is the usual ",
                "cause; compare the fit with a bounded family",
                call. = FALSE, class = "frmtmb_gev_xi_edge")
  } else {
    frm_message("gen_extreme_value: xi ", of, "is ", at, ", within 1e-3 ",
                "of -1, at a maximum (d logLik / d xi = ",
                format(dll, digits = 3), "). Below -0.5 the maximum is ",
                "not regular, so the standard errors do not have their ",
                "usual meaning", class = "frmtmb_gev_xi_near_edge")
  }
  invisible(NULL)
}

#' `d logLik / d xi` at the fit's stop for a scalar `xi`: the objective's
#' gradient in the `xi` coefficient, divided by the link's `dxi / deta`.
#' `NULL` where it cannot be read.
#'
#' @noRd
gev_xi_slope <- function(fit, lp) {
  p <- fit$opt$par
  g <- tryCatch(drop(fit$obj$gr(p)), error = function(e) NULL)
  if (is.null(g) || length(g) != length(p)) return(NULL)
  fx <- if (identical(lp[["par"]], "betad")) fit$frame[["betad_fixed_idx"]]
  free <- setdiff(seq_along(fit$frame[["par_template"]][[lp[["par"]]]]), fx)
  j <- which(names(p) == lp[["par"]])[match(lp[["idx"]], free)]
  if (length(j) != 1L || is.na(j)) return(NULL)
  eta <- fit$estimates[[lp[["par"]]]][lp[["idx"]]]
  # the objective is -logLik, judged in the natural units the optimizer
  # used, as grad_verdict() judges it
  u <- rep_len(fit$par_units %||% 1, length(p))
  -g[j] * u[j] / lp[["link"]]$mu_eta(eta)
}
#' Generalized extreme value family, dpars `mu` (location), `sigma`
#' (scale) and `xi` (shape).
#'
#' @noRd
fam_gen_extreme_value <- function(link = "identity", link_sigma = "log",
                                  link_xi = "log1p") {
  fam_gev_impl(link, link_sigma, link_xi, chart = "xi")
}

#' The GEV family in one of two charts for `xi`. `"xi"` is the family
#' `gen_extreme_value()` returns. `"tmp"` is brms's: the `xi` dpar holds
#' `tmp_xi` on an identity link, and the density maps it into the support
#' the current `mu` and `sigma` leave (`gev_scale_xi()`), REPORTing the
#' `xi` it used as `frmtmb_gev_xi`. Only the sampler tapes it
#' (`gev_xi_chart_frame()`); it has no CDF, so a `cens()` or `trunc()`
#' response keeps the `"xi"` chart.
#'
#' @noRd
fam_gev_impl <- function(link, link_sigma, link_xi, chart) {
  lk <- mu_link(link, "gen_extreme_value")
  lk_sigma <- dpar_link(link_sigma, "sigma", "gen_extreme_value",
                        dpar_links_positive)
  lk_xi <- dpar_link(link_xi, "xi", "gen_extreme_value", dpar_links_xi)
  # log1p holds xi above -1 by construction. Under the identity link
  # nothing does, and below -1 the likelihood is unbounded, so the soft
  # support, which caps an edge row at t = 1e-15, would give a finite
  # and wrong value there: seeds 621 and 622 ran to xi = -1.007 with a
  # row outside the support (dev/reviews/2026-10-13-famuni.md, fix 3)
  xi_guard <- if (identical(lk_xi$name, "log1p")) {
    function(dpars) 0
  } else {
    function(dpars) nan_unless_pos(1 + dpars[["xi"]])
  }
  core_lpdf <- function(y, dpars) {
    ov <- osa_unwrap(y)
    g <- gev_parts(if (is.null(ov)) y else ov$y, dpars)
    # a row outside the support has density 0, which is -Inf on the
    # log scale and a NaN in every gradient that passes through it, so
    # in a mixture() one component that does not reach a row poisoned
    # the row. The log density falls off as -1e8 per unit of t below
    # 0 instead: below exp(-100) a millionth outside the support, with
    # an ordinary gradient that points back in. It is exactly 0 where
    # t > 1e-15, which the clamp in gev_parts() starts at.
    ll <- -dpar_log(dpars, "sigma", lk_sigma) - g[["l"]] - g[["h"]] -
      exp_neg_capped(g[["h"]]) - 1e8 * g[["out"]]
    if (is.null(ov)) ll else ll * ov$keep
  }
  lpdf <- function(y, dpars, aterms) core_lpdf(y, dpars) + xi_guard(dpars)
  lcdf <- function(q, dpars, aterms) {
    exp(-exp_neg_capped(gev_parts(q, dpars)[["h"]])) + xi_guard(dpars)
  }
  lccdf <- function(q, dpars, aterms) {
    log(-expm1(-exp_neg_capped(gev_parts(q, dpars)[["h"]]))) +
      xi_guard(dpars)
  }
  links <- list(mu = lk, sigma = lk_sigma, xi = lk_xi)
  if (identical(chart, "tmp")) {
    # no guard here: gev_scale_xi() keeps every row inside the support
    # and holds its range above -1, where brms's goes on
    lpdf <- function(y, dpars, aterms) {
      frmtmb_gev_xi <- gev_scale_xi(dpars[["xi"]], y, dpars[["mu"]],
                                    dpars[["sigma"]])
      RTMB::REPORT(frmtmb_gev_xi)
      dpars[["xi"]] <- frmtmb_gev_xi
      core_lpdf(y, dpars)
    }
    lcdf <- NULL
    lccdf <- NULL
    links$xi <- get_link("identity", dpar = "xi")
  }
  frmtmb_family(
    "gen_extreme_value",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = c("mu", "sigma", "xi"),
    links = links,
    lpdf = lpdf,
    lcdf = lcdf,
    lccdf = lccdf,
    init_dpars = list(
      # the Gumbel moments: xi = 0 has unbounded support, so every data
      # set is inside it at the start
      mu = function(y, aterms) {
        mean(y) + digamma(1) * stats::sd(y) * sqrt(6) / pi
      },
      sigma = function(y, aterms) stats::sd(y) * sqrt(6) / pi,
      xi = function(y, aterms) 0
    ),
    type = "continuous",
    post = list(
      fit_check = gev_fit_check,
      # the sampler's chart for a scalar xi (gev_xi_chart_frame())
      tmp_xi_family = function() {
        fam_gev_impl(link, link_sigma, link_xi, chart = "tmp")
      },
      mean_fn = function(dpars, aterms) {
        gev_mean(dpars[["mu"]], dpars[["sigma"]], dpars[["xi"]])
      },
      var_fn = function(dpars, aterms) {
        n <- max(length(dpars[["mu"]]), length(dpars[["xi"]]))
        rep_len(gev_var(dpars[["sigma"]], dpars[["xi"]]), n)
      },
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        gev_trunc_mean(dpars[["mu"]], dpars[["sigma"]], dpars[["xi"]],
                       lb, ub)
      }
    ),
    sim = function(dpars, aterms, n) {
      rgev(n, dpars[["mu"]], dpars[["sigma"]], dpars[["xi"]])
    }
  )
}

#' The largest or smallest element of `v`, on the tape: a pairwise
#' tree of the tape-safe `ad_max()` or `ad_min()`, `log2(n)` vectorized
#' steps rather than `n` scalar ones. An odd length repeats its last
#' element, which changes neither.
#'
#' @noRd
ad_tree_extreme <- function(v, f) {
  n <- length(v)
  while (n > 1L) {
    if (n %% 2L) v <- v[c(seq_len(n), n)]
    i <- seq.int(1L, length(v), by = 2L)
    v <- f(v[i], v[i + 1L])
    n <- length(v)
  }
  v
}

#' brms's `scale_xi()`: the GEV shape from an unbounded `tmp_xi`,
#' `lb + (ub - lb) plogis(tmp_xi)`, where `lb` and `ub` are the two values
#' `-1 / min(x)` and `-1 / max(x)` for `x = (y - mu) / sigma`, in order:
#' the interval of `xi` that keeps every observed row inside the support.
#'
#' With one change: `lb` is held at -1 or above. brms's interval reaches
#' below -1, where the likelihood is unbounded; frmtmb's `xi` is above -1
#' on every route, and the fit's `log1p` link, which the draws are carried
#' back to, cannot represent a value below it. On the review's data at
#' xi = -0.85 (seed 1402, n = 40) brms's 95 percent interval reaches
#' -1.238; a first version without the hold gave a NaN `xi` for every draw
#' below -1 (dev/famuni-out/tmpxi-sample.txt). `ub` is held at `lb` or
#' above for the case where every row lies above `mu`, where the two
#' values are both below 0.
#'
#' @noRd
gev_scale_xi <- function(tmp, y, mu, sigma) {
  x <- (y - mu) / sigma
  b1 <- -1 / ad_tree_extreme(x, ad_min)
  b2 <- -1 / ad_tree_extreme(x, ad_max)
  lb <- ad_max(ad_min(b1, b2), -1)
  ub <- ad_max(ad_max(b1, b2), lb)
  lb + (ub - lb) / (1 + exp(-tmp))
}

#' The sampler's chart for a scalar GEV shape
#'
#' A `gen_extreme_value()` `xi` without a formula is sampled as brms
#' samples it, as `tmp_xi` on the whole line, mapped into the support the
#' current `mu` and `sigma` leave by brms's `scale_xi()`. Sampled on its
#' `log1p` scale, the support's edge is a wall in the sampled space: on
#' the review's data (n = 40, `y ~ 1`, 4 chains x 2000) 3388 of 4000
#' transitions diverged at xi = -0.85, and none did through `tmp_xi`
#' (dev/famuni-out/tmpxi.txt).
#'
#' `gev_xi_chart_terms()` lists each such `xi`: its response, linear
#' predictor key, parameter component and index, the fit's link, and
#' `blocked`, the reasons the chart cannot be used (a `cens()` or
#' `trunc()` term, whose CDF the chart does not carry, a second GEV
#' response, which would share the REPORTed name, or a [mixture()]
#' component, whose `xi1`, `xi2`, ... the chart does not map).
#' `gev_xi_chart_frame()`
#' returns the frame with each unblocked one in the `tmp_xi` chart, and
#' `gev_xi_chart_start()` carries a parameter list's `xi` coefficients
#' to `tmp_xi` at the fit's `mu` and `sigma`. `dpar_has_predictor()`
#' tells `xi ~ 1`, which brms samples on `log1p`, from no formula.
#'
#' @param fit A `frmtmb_fit`.
#' @param frame A model frame.
#' @param terms What `gev_xi_chart_terms()` returned.
#' @param pars A parameter list in the frame's template.
#' @param spec,lp A model spec and one of its frame's linear predictors.
#' @keywords internal
#' @name frmtmb-gev-xi-chart
#' @aliases dpar_has_predictor
#' @usage dpar_has_predictor(spec, lp)
#' @rawNamespace export(gev_xi_chart_terms, gev_xi_chart_frame,
#'   gev_xi_chart_start, dpar_has_predictor)
NULL

#' @rdname frmtmb-gev-xi-chart
gev_xi_chart_terms <- function(fit) {
  out <- list()
  for (rs in fit$spec$responses) {
    fname <- rs$family[["family"]]
    # a mixture's GEV components carry xi1, xi2, ...; the chart covers a
    # plain GEV response only, and a mixture is listed so that the
    # sampler says so (a two-GEV mixture at xi -0.3 had 84 of 4000
    # transitions divergent, dev/famuni-rev2-out/mixsample.txt)
    mix <- !identical(fname, "gen_extreme_value")
    if (mix && !(startsWith(fname, "mixture(") &&
                   grepl("gen_extreme_value", fname, fixed = TRUE))) next
    r <- rs$resp_name
    dps <- if (mix) grep("^xi[0-9]+$", rs$family[["dpars"]], value = TRUE)
           else "xi"
    for (dp in dps) {
      key <- linpred_key(r, dp)
      lp <- fit$frame[["linpreds"]][[key]]
      if (is.null(lp) || !is.null(lp[["constant"]]) ||
            !is.null(lp[["equate"]]) || !is.null(lp[["nl_body"]]) ||
            dpar_has_predictor(fit$spec, lp)) next
      av <- fit$frame[["aterm_values"]][[r]] %||% list()
      blocked <- c(if (mix) "a mixture()",
                   if (!is.null(av[["cens"]])) "cens()",
                   if (!is.null(av[["trunc_lb"]]) ||
                         !is.null(av[["trunc_ub"]])) "trunc()")
      out[[length(out) + 1L]] <- list(resp = r, dpar = dp, key = key,
                                      par = lp[["par"]], idx = lp[["idx"]],
                                      link = lp[["link"]], blocked = blocked)
    }
  }
  # the REPORTed name is one per model, so two charted responses would
  # share it; a mixture's components are not charted anyway
  plain <- vapply(out, function(tm) !"a mixture()" %in% tm$blocked, NA)
  if (sum(plain) > 1L) {
    out[plain] <- lapply(out[plain], function(tm) {
      tm$blocked <- c(tm$blocked, "a second gen_extreme_value() response")
      tm
    })
  }
  out
}

#' @rdname frmtmb-gev-xi-chart
gev_xi_chart_frame <- function(frame, terms) {
  for (tm in terms) {
    if (length(tm$blocked)) next
    fam <- frame[["spec"]]$responses[[tm$resp]]$family
    tmp_fam <- fam[["post"]][["tmp_xi_family"]]()
    frame[["spec"]]$responses[[tm$resp]]$family <- tmp_fam
    frame[["linpreds"]][[tm$key]][["link"]] <-
      get_link("identity", dpar = "xi")
  }
  frame
}

#' @rdname frmtmb-gev-xi-chart
gev_xi_chart_start <- function(fit, terms, pars) {
  multi <- length(fit$spec$responses) > 1L
  for (tm in terms) {
    if (length(tm$blocked)) next
    rs <- if (multi) tm$resp
    lpr <- function(dp) {
      frm_linpred(fit, dpar = dp, type = "response", resp = rs)
    }
    x <- (fit$frame[["y"]][[tm$resp]] - lpr("mu")) / lpr("sigma")
    b <- -1 / range(x)
    xi <- lpr("xi")[1L]
    q <- (xi - min(b)) / (max(b) - min(b))
    # the maximum is inside the support, so q is inside (0, 1); held off
    # its ends so a start on the edge has a finite tmp_xi
    pars[[tm$par]][tm$idx] <- stats::qlogis(min(max(q, 1e-8), 1 - 1e-8))
  }
  pars
}

# ------------------------------------------------------- discrete weibull

#' `q^k` for an integer `q`, with `0^k = 0` and `q^k = 0` for `q < 0`:
#' the survivor `mu^(q^k)` is then 1 there, as it is below the support.
#' Written as `exp(k log(max(q, 1)))` times `min(max(q, 0), 1)`, which
#' is 1 at every integer `q >= 1` and 0 below, because `0^k`
#' differentiates to `0 * log(0)`, a `NaN` in the gradient with respect
#' to `k`. The clamps are the tape-safe ones, so this also holds when
#' `q` is the response on the one-step tape, where it cannot be
#' compared.
#'
#' @noRd
dw_pow <- function(q, k) {
  ad_max(ad_min(q, 1), 0) * exp(k * log(ad_max(q, 1)))
}

#' The mean and the second moment of the discrete Weibull, off the tape:
#' `E[Y] = sum_{y >= 1} S(y)` and `E[Y^2] = sum_{y >= 1} (2y - 1) S(y)`
#' with `S(y) = mu^(y^k)`. The first `m` terms are summed, and the rest
#' through the Euler-Maclaurin form, the integral from `m + 1` plus half
#' the first omitted term less a twelfth of its slope, which with
#' `m = 1000` is exact to about 1e-12 relative. The integrals are upper
#' incomplete gammas, `int_c^Inf t^j exp(L t^k) dt =
#' Gamma((j + 1) / k, -L c^k) / (k (-L)^((j + 1) / k))`, taken on the
#' log scale because `gamma((j + 1) / k)` overflows for a small shape.
#'
#' @noRd
dw_moments <- function(mu, shape, m = 1000L) {
  n <- max(length(mu), length(shape))
  mu <- rep_len(mu, n)
  shape <- rep_len(shape, n)
  y <- seq_len(m)
  out <- vapply(seq_len(n), function(i) {
    L <- log(mu[i])
    k <- shape[i]
    if (!is.finite(L) || !is.finite(k) || k <= 0 || L > 0) {
      return(c(NaN, NaN))
    }
    if (L == 0) return(c(Inf, Inf))
    S <- exp(L * y^k)
    m1 <- sum(S)
    m2 <- sum((2 * y - 1) * S)
    c0 <- m + 1
    f0 <- exp(L * c0^k)
    if (f0 > 0) {
      x0 <- -L * c0^k
      tail_int <- function(j) {
        a <- (j + 1) / k
        exp(lgamma(a) + stats::pgamma(x0, a, lower.tail = FALSE,
                                      log.p = TRUE) - log(k) - a * log(-L))
      }
      # S'(t) = S(t) L k t^(k - 1)
      d0 <- f0 * L * k * c0^(k - 1)
      m1 <- m1 + tail_int(0) + f0 / 2 - d0 / 12
      m2 <- m2 + 2 * tail_int(1) - tail_int(0) +
        (2 * c0 - 1) * f0 / 2 - (2 * f0 + (2 * c0 - 1) * d0) / 12
    }
    c(m1, m2)
  }, numeric(2))
  list(mean = out[1L, ], m2 = out[2L, ])
}

#' The discrete Weibull partial mean over `lb..ub` (both included),
#' divided by its mass, off the tape. A finite window is summed outright;
#' an open upper end uses the moments above, less the terms below `lb`.
#'
#' @noRd
dw_trunc_mean <- function(mu, shape, lb, ub) {
  n <- max(length(mu), length(shape), length(lb), length(ub))
  mu <- rep_len(mu, n)
  shape <- rep_len(shape, n)
  lb <- rep_len(pmax(lb, 0), n)
  ub <- rep_len(ub, n)
  vapply(seq_len(n), function(i) {
    L <- log(mu[i])
    k <- shape[i]
    S <- function(y) ifelse(y <= 0, 1, exp(L * pmax(y, 0)^k))
    mass <- S(lb[i]) - S(ub[i] + 1)
    if (is.finite(ub[i]) && ub[i] - lb[i] <= 1e6) {
      yy <- seq(lb[i], ub[i])
      return(sum(yy * (S(yy) - S(yy + 1))) / mass)
    }
    # sum_{y >= lb} y P(y) = lb S(lb) + sum_{y > lb} S(y)
    full <- dw_moments(mu[i], k)$mean
    below <- if (lb[i] >= 1) sum(S(seq_len(lb[i]))) else 0
    (lb[i] * S(lb[i]) + full - below) / mass
  }, numeric(1))
}

#' Discrete Weibull family, dpars `mu` (in `(0, 1)`, `P(Y >= 1)`) and
#' `shape`.
#'
#' @noRd
fam_discrete_weibull <- function(link = "logit", link_shape = "log") {
  lk <- mu_link(link, "discrete_weibull")
  lk_shape <- dpar_link(link_shape, "shape", "discrete_weibull",
                        dpar_links_positive)
  frmtmb_family(
    "discrete_weibull",
    accepts_aterms = c("weights", "cens", "trunc"),
    dpars = c("mu", "shape"),
    links = list(mu = lk, shape = lk_shape),
    lpdf = function(y, dpars, aterms) {
      # log(S(y) - S(y + 1)) = a + log(1 - exp(b - a)) with a = L y^k and
      # b = L (y + 1)^k, which keeps the difference of two numbers near
      # 1 off the tape's arithmetic; L = log(mu) comes from the log odds
      ov <- osa_unwrap(y)
      yy <- if (is.null(ov)) y else ov$y
      L <- dpar_log(dpars, "mu", lk)
      k <- dpars[["shape"]]
      a <- L * dw_pow(yy, k)
      ll <- a + log(-expm1(L * (dw_pow(yy + 1, k) - dw_pow(yy, k))))
      if (is.null(ov)) ll else ll * ov$keep
    },
    lcdf = function(q, dpars, aterms) {
      -expm1(dpar_log(dpars, "mu", lk) * dw_pow(q + 1, dpars[["shape"]]))
    },
    lccdf = function(q, dpars, aterms) {
      dpar_log(dpars, "mu", lk) * dw_pow(q + 1, dpars[["shape"]])
    },
    valid_y = count_y("discrete_weibull"),
    init_dpars = list(
      # at shape 1 the family is the geometric with P(Y >= y) = mu^y,
      # whose mean is mu / (1 - mu)
      mu = function(y, aterms) {
        m <- max(mean(y), 0.05)
        m / (1 + m)
      },
      shape = function(y, aterms) 1
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        dw_moments(dpars[["mu"]], dpars[["shape"]])$mean
      },
      var_fn = function(dpars, aterms) {
        mm <- dw_moments(dpars[["mu"]], dpars[["shape"]])
        mm$m2 - mm$mean^2
      },
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        dw_trunc_mean(dpars[["mu"]], dpars[["shape"]], lb, ub)
      },
      # mu is P(Y >= 1), not a location, so a mixture component starts
      # it from the quantile it is placed at: at shape 1 the median m
      # has mu^(m + 1) = 1/2
      mixture_mu_start = function(y, aterms, p) {
        q <- stats::quantile(y, p, names = FALSE)
        min(max(0.5^(1 / (q + 1)), 0.02), 0.98)
      }
    ),
    sim = function(dpars, aterms, n) {
      u <- stats::runif(n)
      k <- dpars[["shape"]]
      pmax(ceiling((log1p(-u) / log(dpars[["mu"]]))^(1 / k) - 1), 0)
    }
  )
}
