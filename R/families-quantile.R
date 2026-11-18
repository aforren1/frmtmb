# quantile_elf(): quantile regression through the extended log-F (ELF)
# density of Fasiolo et al. (2021), a smoothed asymmetric Laplace. It
# lives in a file of its own so that families.R, which other work also
# edits, carries only its registry line.

#' Smoothed quantile regression: the extended log-F family
#'
#' `quantile_elf()` fits quantile regression with the extended log-F
#' (ELF) density of Fasiolo et al. (2021). qgam fits the same density.
#' It is a smoothed [asym_laplace()]. Its log density has curvature in
#' the location everywhere. So a fit has Hessian standard errors, and
#' the Laplace approximation of a model with a group-level term in `mu`
#' has a maximum. `asym_laplace()` has neither, because its log density
#' is piecewise linear in the location.
#'
#' @section The density:
#' The parameters are `mu`, `sigma` and `quantile`, as in
#' `asym_laplace()`. A fixed smoothing `lambda` (the argument
#' `smoothing`) is added. Let `e = y - mu`, `p = quantile`,
#' `a = lambda (1 - p)` and `b = lambda p`. Then the log density is
#'
#' `(1 - p) e / sigma - lambda log(1 + exp(e / (lambda sigma))) -
#' log(lambda sigma) - log B(a, b)`.
#'
#' This is the ELF density of Fasiolo et al. (2021) with their
#' `lambda`. It is also the density of qgam's `elf()` family, whose
#' argument `co` is `lambda sigma`. The response is
#' `y = mu + lambda sigma logit(Z)` with `Z ~ Beta(a, b)`. So:
#'
#' - the distribution function is the incomplete beta
#'   `pbeta(plogis(e / (lambda sigma)), a, b)`;
#' - the mean is `mu + lambda sigma (digamma(a) - digamma(b))`;
#' - the variance is `(lambda sigma)^2 (trigamma(a) + trigamma(b))`.
#'
#' As `lambda` goes to 0, the density goes to the asymmetric Laplace
#' density of `asym_laplace()` with the same `mu`, `sigma` and
#' `quantile`.
#'
#' Hold `quantile` at the quantile you want, as in
#' `bf(y ~ x, quantile = 0.9)`. If you do not, `quantile` is estimated,
#' as in `asym_laplace()`. Then `mu` is no longer a quantile of
#' interest.
#'
#' @section What mu estimates:
#' At its maximum, the likelihood sets the mean of
#' `plogis(e / (lambda sigma))` to `1 - p`. So `mu` is the `p` quantile
#' of `y - lambda sigma L`, where `L` is a standard logistic variable.
#' This is the quantile of the response smoothed by a logistic kernel
#' of scale `lambda sigma`. It is not the quantile of `y` itself.
#'
#' The difference is 0 for a residual density that is symmetric about
#' the quantile, such as the median of a symmetric law. It grows with
#' the square of `lambda`. It does not shrink as the sample grows.
#'
#' Let `f` be the residual density and `q` the quantile. For a kernel
#' width that is small against the scale on which `f` changes, the
#' difference is about `-(pi^2 / 6) (lambda sigma)^2 f'(q) / f(q)`.
#' This is an approximation only. At `smoothing = 1` it can be off by a
#' factor of 2. Near an end of the support it can have the wrong sign.
#'
#' The table gives the bias of `mu` as `n` goes to infinity. The units
#' are the residual sd, and the standard error of the quantile estimate
#' at `n = 1000`. It was computed by numerical integration
#' (`dev/squant-pop.R` in the source repository), not from the
#' approximation.
#'
#' \tabular{lrrrr}{
#'   \strong{residuals, quantile} \tab \strong{0.25: sd} \tab
#'     \strong{0.25: SE} \tab \strong{1: sd} \tab \strong{1: SE} \cr
#'   gaussian, 0.5 \tab 0 \tab 0 \tab 0 \tab 0 \cr
#'   gaussian, 0.75 \tab 0.0068 \tab 0.17 \tab 0.080 \tab 2.1 \cr
#'   gaussian, 0.9 \tab 0.0040 \tab 0.08 \tab 0.057 \tab 1.2 \cr
#'   t with 3 df, 0.75 \tab 0.0102 \tab 0.35 \tab 0.103 \tab 3.4 \cr
#'   lognormal, 0.5 \tab 0.0137 \tab 0.81 \tab 0.083 \tab 4.5 \cr
#'   lognormal, 0.9 \tab 0.0069 \tab 0.08 \tab 0.109 \tab 1.2
#' }
#'
#' @section Choosing the smoothing:
#' `smoothing` is a fixed constant of the family, like `k` of [huber()].
#' It is not estimated. Maximum likelihood would choose the `lambda`
#' that fits the shape of the residuals. A large `lambda` moves `mu`
#' away from the quantile.
#'
#' qgam chooses the absolute width `lambda sigma` from a rule that
#' balances the bias against the variance of the estimate (its argument
#' `err`). With this rule, the width shrinks as `n^(-1/3)`. qgam then
#' chooses `sigma`, its "learning rate", to calibrate its posterior
#' credible intervals. It does not use the likelihood for either. Here
#' `sigma` is a parameter of the density, estimated by maximum
#' likelihood, and `lambda` is relative to it. The default is fixed. It
#' does not change with the data. Use these steps:
#'
#' 1. Keep the default `0.25` for a model with fixed effects only and a
#'    quantile between 0.1 and 0.9.
#' 2. Use a smaller value, such as `0.1`, for a large sample. There a
#'    bias of a fixed size is a larger part of the standard error.
#' 3. Use a larger value, such as `1`, for a model with a group-level
#'    term in `mu` (see "Group-level terms"). Accept a larger bias away
#'    from the median.
#' 4. To see the effect of the choice, fit two values and compare `mu`.
#'
#' @section Standard errors:
#' The default `summary()`, `fixef()`, `confint()` and `hypothesis()`
#' show the Hessian standard errors of the ELF model. The residuals of
#' a quantile regression do not usually follow the ELF density. Then
#' these standard errors are not those of the quantile estimate. Away
#' from the median, their 95% intervals cover only about 0.7 to 0.8
#' (table below). [frm_bootstrap()] draws from the fitted ELF, so it has
#' the same problem.
#'
#' Use the sandwich covariance [vcov_cluster()] with one cluster per
#' row. It is the quantile regression's own covariance, with a kernel
#' estimate of the residual density. Do these steps:
#'
#' 1. Add a column that numbers the rows, `dd$row <- seq_len(nrow(dd))`,
#'    before the fit.
#' 2. Make the covariance: `V <- vcov_cluster(fit, ~ row, "CR0",
#'    full = TRUE)`.
#' 3. Give it to the inference methods: `summary(fit, vcov = V)` and
#'    `confint(fit, vcov = V)`.
#'
#' `vcov(fit, cluster = ~ row)` gives the coefficient block of the
#' sandwich, with the default small-sample correction of
#' [vcov_cluster()].
#'
#' The table gives the ratio of the mean standard error to the Monte
#' Carlo sd of the slope. It also gives the coverage of 95% Wald
#' intervals for the true quantile slope. The settings are
#' `smoothing = 0.25`, 400 rows, and 300 data sets per row
#' (`dev/squant-se.R`). The bootstrap had 200 draws on 3 data sets of
#' the second row.
#'
#' \tabular{lrrrrr}{
#'   \strong{residuals, quantile} \tab \strong{Hessian} \tab
#'     \strong{coverage} \tab \strong{sandwich} \tab \strong{coverage}
#'     \tab \strong{bootstrap} \cr
#'   gaussian, 0.5 \tab 0.829 \tab 0.893 \tab 0.977 \tab 0.943 \tab \cr
#'   gaussian, 0.9 \tab 0.592 \tab 0.777 \tab 1.002 \tab 0.913 \tab
#'     0.39 to 0.44 \cr
#'   sd growing with x, 0.9 \tab 0.591 \tab 0.767 \tab 0.995 \tab 0.907
#'     \tab \cr
#'   lognormal, 0.5 \tab 1.028 \tab 0.953 \tab 0.985 \tab 0.950 \tab \cr
#'   lognormal, 0.9 \tab 0.505 \tab 0.700 \tab 1.012 \tab 0.913 \tab
#' }
#'
#' With 100 rows and gaussian residuals, the sandwich coverage is 0.93
#' at the median and 0.84 at the 0.9 quantile. A larger smoothing helps
#' there (0.91 at `smoothing = 1`), because the kernel estimate of the
#' density is less noisy. `quantreg::summary.rq(se = "nid")` covered
#' 0.94 to 0.96 in every row.
#'
#' @section Group-level terms:
#' With `(1 | g)` in `mu`, the Laplace approximation falls as the group
#' sd grows without end, so it has a highest point. That point can be
#' the boundary `sd = 0`. For `asym_laplace()` it rises without end
#' instead. The test
#' designs had 80 rows in groups of 2 or 4, or 100 rows in groups of 10,
#' with 10 data sets each. On them, every fit stopped at code 0 with no
#' warning.
#'
#' But at a small smoothing the approximation is poor. Its profile in
#' the group sd is a sawtooth with several local maxima. The cause is
#' the log determinant of the inner Hessian, which jumps as a group's
#' mode crosses a data point. So a fit can stop on a local maximum, with
#' code 0 and no warning.
#'
#' Compared with `quadrature = TRUE` at the median, the log group sd
#' differed by up to 0.87 at `smoothing = 0.1`, 0.59 at `0.25`, 0.31 at
#' `0.5` and 0.18 at `1`. There were 40 such fits and 10 at the 0.9
#' quantile. Of these, at least 7, 2, 1 and 0 stopped below a higher
#' point of a 10-point profile (`dev/squant-re.R`). On one data set at
#' the 0.9 quantile, a finer profile was 3.42 higher than the fit at
#' `smoothing = 0.1`, and 1.76 higher at `0.25`. At the 0.9 quantile
#' with groups of 4, the log sd differed from quadrature by 0.28 even at
#' `smoothing = 1`. Quadrature fitted all 200 of these models.
#'
#' An exact `gp()` term is slow at a small smoothing. At the default
#' `0.25`, `y ~ gp(x) + x2` warned that the Laplace approximation is
#' not reliable, and took 129 s at 40 rows and 501 s at 80. At
#' `smoothing = 1` it took 0.6 s and 3.1 s. The time is in the inner
#' Newton search (`dev/reviews/squant-short-gp.R`).
#'
#' Use these steps for a model with a group-level term:
#'
#' 1. Use `smoothing = 1`.
#' 2. For the estimates, use `quadrature = TRUE` where it applies (one
#'    scalar grouping factor).
#' 3. For standard errors, know what each route gives:
#'    - The Laplace fit: `vcov_cluster(fit, ~ g, "CR0", full = TRUE)`
#'      is the cluster sandwich, robust to the ELF's misspecification.
#'      Its coverage was not measured for this family.
#'    - The quadrature fit: [vcov_cluster()] refuses it, so only the
#'      Hessian standard errors of the ELF model are available. They are
#'      not calibrated for the quantile estimate.
#'
#' @section How to choose between asym_laplace() and quantile_elf():
#' Both families estimate a quantile. Use the steps below to choose.
#'
#' 1. If the model has a group-level term (`(1 | g)`, `s()`, `gp()`)
#'    in `mu`, and you fit it with the Laplace approximation, use
#'    `quantile_elf()`. `asym_laplace()` has no maximum there, and the
#'    fit is refused.
#' 2. If you need standard errors or [confint()] from the fit itself,
#'    use `quantile_elf()` with the sandwich covariance:
#'    `confint(fit, vcov = vcov_cluster(fit, ~ row, "CR0", full = TRUE))`
#'    (see "Standard errors"). `asym_laplace()` has no Hessian standard
#'    error.
#' 3. If you need the exact sample quantile regression estimate, as
#'    `quantreg::rq()` gives it, and only fixed effects, use
#'    `asym_laplace()`. Or use `quantile_elf()` with a small `smoothing`.
#' 4. If you will sample the model with frmtmb.sample's `frm_sample()`
#'    to match a brms model, use `asym_laplace()`. brms has no ELF
#'    family.
#'
#' @references
#' Fasiolo, M., Wood, S. N., Zaffran, M., Nedellec, R. and Goude, Y.
#' (2021). Fast calibrated additive quantile regression. *Journal of
#' the American Statistical Association* 116, 1402-1412.
#' \doi{10.1080/01621459.2020.1725521}
#'
#' Fasiolo, M., Wood, S. N., Zaffran, M., Nedellec, R. and Goude, Y.
#' (2021). qgam: Bayesian nonparametric quantile regression modeling in
#' R. *Journal of Statistical Software* 100(9), 1-31.
#' \doi{10.18637/jss.v100.i09}
#'
#' Koenker, R. and Bassett, G. (1978). Regression quantiles.
#' *Econometrica* 46, 33-50.
#'
#' @param link The link of `mu`: the set of `asym_laplace()`.
#' @param link_sigma The link of `sigma`: `"log"` (default),
#'   `"identity"`, `"softplus"` or `"squareplus"`.
#' @param link_quantile The link of `quantile`: `"logit"` (default) or
#'   `"identity"`.
#' @param smoothing The smoothing `lambda`, one finite number above 0
#'   (default `0.25`). It is the width of the smoothing kernel in units
#'   of `sigma`. Smaller values give a smaller bias and less curvature.
#' @return A `frmtmb_family`.
#' @seealso [asym_laplace()] for the unsmoothed density, which brms has;
#'   [vcov_cluster()] for the sandwich covariance.
#' @examples
#' set.seed(1)
#' n <- 300
#' dd <- data.frame(x = runif(n), row = seq_len(n))
#' dd$y <- 1 + 2 * dd$x + (0.5 + dd$x) * rnorm(n)
#' fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(),
#'            data = dd)
#'
#' # the default summary shows the ELF model's Hessian standard errors
#' summary(fit)
#'
#' # the sandwich standard errors of the quantile estimate
#' V <- vcov_cluster(fit, ~ row, "CR0", full = TRUE)
#' summary(fit, vcov = V)
#' confint(fit, vcov = V)
#' sqrt(diag(vcov(fit, cluster = ~ row)))
#' @export
quantile_elf <- function(link = "identity", link_sigma = "log",
                         link_quantile = "logit", smoothing = 0.25) {
  link <- link_arg_value(substitute(link), link,
                         brms_mu_links[["quantile_elf"]], "identity")
  fam_quantile_elf(link, link_sigma, link_quantile, smoothing)
}

#' The ELF family object. `smoothing` is the `lambda` of Fasiolo et al.
#' (2021), relative to `sigma`.
#'
#' @noRd
fam_quantile_elf <- function(link = "identity", link_sigma = "log",
                             link_quantile = "logit", smoothing = 0.25) {
  lk_sigma <- dpar_link(
    link_sigma, "sigma", "quantile_elf", dpar_links_positive)
  lk_quantile <- dpar_link(
    link_quantile, "quantile", "quantile_elf", dpar_links_unit)
  if (!is.numeric(smoothing) || length(smoothing) != 1L ||
        !is.finite(smoothing) || smoothing <= 0) {
    frm_stop("quantile_elf(smoothing =): the smoothing must be one finite ",
             "number above 0, not ", arg_desc(smoothing), call. = FALSE)
  }
  lam <- as.numeric(smoothing)
  frmtmb_family(
    "quantile_elf",
    accepts_aterms = c("weights", "cens", "trunc", "mi"),
    dpars = c("mu", "sigma", "quantile"),
    links = list(mu = mu_link(link, "quantile_elf"),
                 sigma = lk_sigma, quantile = lk_quantile),
    lpdf = function(y, dpars, aterms) {
      pq <- dpar_complement(dpars, "quantile", lk_quantile)
      elf_lpdf(y, dpars[["mu"]], dpars[["sigma"]], pq$p, pq$q, lam)
    },
    cdf_log = function(q, dpars, aterms) {
      pq <- dpar_complement(dpars, "quantile", lk_quantile)
      elf_cdf_log(q, dpars[["mu"]], dpars[["sigma"]], pq$p, pq$q, lam)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y),
      sigma = function(y, aterms) stats::sd(y) / 2,
      quantile = function(y, aterms) 0.5
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        p <- dpars[["quantile"]]
        dpars[["mu"]] + lam * dpars[["sigma"]] *
          (digamma(lam * (1 - p)) - digamma(lam * p))
      },
      var_fn = function(dpars, aterms) {
        p <- dpars[["quantile"]]
        (lam * dpars[["sigma"]])^2 *
          (trigamma(lam * (1 - p)) + trigamma(lam * p))
      }
    ),
    sim = function(dpars, aterms, n) {
      pq <- dpar_complement(dpars, "quantile", lk_quantile)
      elf_draw(n, dpars[["mu"]], dpars[["sigma"]], pq$p, pq$q, lam)
    }
  )
}

#' The ELF log density, on doubles or on the tape. `p` and `q1` are the
#' quantile and its complement, each held to full relative precision.
#'
#' Written as the asymmetric Laplace term plus a smooth correction:
#' `lambda log(1 + e^t) = max(r, 0) / sigma + lambda log(1 + e^-|t|)`
#' with `t = r / (lambda sigma)`. So `exp()` never overflows, and the
#' kink of `|u|` cancels against the kink of the correction, whose sum
#' has the derivatives of `lambda log(1 + e^t)` everywhere, `u = 0`
#' included.
#'
#' @noRd
elf_lpdf <- function(y, mu, sigma, p, q1, lam) {
  u <- (y - mu) / sigma
  au <- abs(u)
  -0.5 * (au + (p - q1) * u) - lam * log1p(exp(-au / lam)) -
    log(lam * sigma) - lbeta_ad(lam * q1, lam * p)
}

# Beyond this many kernel widths from mu, plogis() of the far tail is
# below 1e-299, the end of the range where log_ibeta_two() keeps its
# digits; each tail is continued there by its exact exponential slope.
elf_tail_t <- 690

#' The ELF log CDF and log survivor: `log I_x(a, b)` and
#' `log I_{1-x}(b, a)` at `x = plogis(t)`, `t = (q - mu) / (lambda
#' sigma)`, from `log_ibeta_two()`, which takes both tails directly.
#'
#' Past `|t| = elf_tail_t` the pair is read at the edge and continued.
#' Below, `I_x(a, b) = x^a (1 - x)^b / (a B(a, b)) (1 + O(x))` with
#' `log x = t + O(e^t)`, so `log F` falls by exactly `a` per unit of `t`
#' to the precision of a double; above, `log S` by `b`. Each other tail
#' is the complement of the continued one.
#'
#' The complement is switched in by `sign_pos()`, whose derivative is 0,
#' and not added as a difference that is 0 inside. That difference was
#' 0 in value but not in derivative: the complement's slope `S / F` is
#' huge where `F` is small, and the tape's two copies of it did not
#' cancel, so `d log F / d mu` was off by about `eps / F` and wrong once
#' `F < 1e-16` (dev/reviews/2026-11-17-squant.md, item 1). The clamp and
#' the slopes use `elf_pos0()`, so at `|t| = 690` itself the derivatives
#' are those of the pair read inside.
#'
#' @noRd
elf_cdf_log <- function(q, mu, sigma, p, q1, lam) {
  a <- lam * q1
  b <- lam * p
  t <- (q - mu) / (lam * sigma)
  tc <- t - elf_pos0(t - elf_tail_t)
  tc <- tc + elf_pos0(-elf_tail_t - tc)
  dl <- t - tc
  n <- max(length(t), length(a), length(b))
  tc <- cdf_rows(tc, n)
  dl <- cdf_rows(dl, n)
  a <- cdf_rows(a, n)
  b <- cdf_rows(b, n)
  x <- 1 / (1 + exp(-tc))
  xc <- 1 / (1 + exp(tc))
  lp <- log_ibeta_two(x, xc, a, b)
  lf_low <- lp[["lF"]] - a * elf_pos0(-dl)
  ls_high <- lp[["lS"]] - b * elf_pos0(dl)
  # 1 past the edge on that side and exactly 0 at and inside it
  w_low <- 0.5 * (1 - sign_pos(dl))
  w_high <- 0.5 * (1 - sign_pos(-dl))
  # the cap keeps the complement finite where its weight is 0, so the
  # product is 0 rather than NaN
  cap <- -1e-300
  lf <- (1 - w_high) * lf_low + w_high * log1m_exp_ad(ad_min(ls_high, cap))
  ls <- (1 - w_low) * ls_high + w_low * log1m_exp_ad(ad_min(lf_low, cap))
  list(lF = lf, lS = ls)
}

#' `max(u, 0)` with derivative 0 at `u = 0`, where `ad_pos()` has 1
#' (RTMB's `abs()` has derivative 1 at 0). It is exactly 0 for `u <= 0`
#' and `u` for `u > 0`. A term `elf_pos0(dl)` therefore passes no
#' derivative of `dl` inside the edge, where `dl` is exactly 0.
#'
#' @noRd
elf_pos0 <- function(u) u + ad_pos(-u)

#' Draws from the ELF: `mu + lambda sigma (log G_a - log G_b)` with
#' independent unit-rate gamma variables, which is `logit(Z)` for
#' `Z ~ Beta(a, b)`. A small shape's gamma variable underflows to 0, so
#' its log is drawn as `log G_{a + 1} + log(U) / a`, which has the same
#' law.
#'
#' @noRd
elf_draw <- function(n, mu, sigma, p, q1, lam) {
  a <- rep_len(lam * q1, n)
  b <- rep_len(lam * p, n)
  lga <- log(stats::rgamma(n, a + 1)) + log(stats::runif(n)) / a
  lgb <- log(stats::rgamma(n, b + 1)) + log(stats::runif(n)) / b
  mu + lam * sigma * (lga - lgb)
}

#' The ELF quantile function on doubles, used by the tests.
#'
#' `qbeta()` alone is not enough: at small shapes it is inaccurate and
#' warns (at `p = 0.01`, `smoothing = 0.05`, `pr = 0.5` it gave `t =
#' 24.6` where the root is 68.3; dev/reviews/2026-11-17-squant.md, item
#' 5). Its answer is only the start. Newton steps then solve
#' `log F(t) = log(pr)` (or `log S(t) = log(1 - pr)` above one half) on
#' `elf_cdf_log()`, in the kernel units `t`. Both logs are concave in
#' `t`, because the density is log-concave, so the steps converge from
#' any start.
#'
#' @noRd
elf_quantile <- function(pr, mu, sigma, p, lam) {
  n <- max(length(pr), length(mu), length(sigma), length(p))
  pr <- rep_len(pr, n)
  p <- rep_len(p, n)
  a <- lam * (1 - p)
  b <- lam * p
  t <- suppressWarnings(stats::qlogis(stats::qbeta(pr, a, b)))
  t[!is.finite(t)] <- 0
  up <- pr > 0.5
  goal <- ifelse(up, log1p(-pr), log(pr))
  for (i in seq_len(200L)) {
    # mu 0 and sigma 1 / lambda put y on the scale of t
    cl <- elf_cdf_log(t, 0, 1 / lam, p, 1 - p, lam)
    lv <- ifelse(up, cl[["lS"]], cl[["lF"]])
    ld <- elf_lpdf(t, 0, 1 / lam, p, 1 - p, lam)
    # d log F / dt = f / F and d log S / dt = -f / S
    slope <- ifelse(up, -1, 1) * exp(ld - lv)
    step <- (lv - goal) / slope
    t <- t - step
    if (all(abs(step) <= 1e-13 * pmax(1, abs(t)))) break
  }
  mu + lam * sigma * t
}
