#' The Royston and Parmar flexible parametric survival family
#'
#' Models the log cumulative hazard as a natural cubic spline in log
#' time. It is the flexible parametric survival model of Royston and
#' Parmar (2002), and it is parameterized exactly as
#' `flexsurv::flexsurvspline()` parameterizes it, so a `gamma` vector
#' means the same curve in both packages and the two log likelihoods are
#' the same number.
#'
#' Writing `x = log(t)`, the family fits
#'
#' \deqn{g(S(t)) = \gamma_0 + \gamma_1 x + \sum_j \gamma_{j+1} v_j(x)}
#'
#' where `v_j` are the natural cubic spline basis functions of Royston
#' and Parmar, linear beyond the boundary knots, and `g` is chosen by
#' `scale`:
#'
#' \describe{
#'   \item{`"hazard"`}{`g(S) = log(-log(S))`, the log cumulative hazard.
#'     Covariates on the first coefficient are PROPORTIONAL HAZARDS.
#'     With no interior knots this is a Weibull model.}
#'   \item{`"odds"`}{`g(S) = log(1/S - 1)`, the log cumulative odds of
#'     failure. Covariates on the first coefficient are proportional
#'     odds. With no interior knots this is a log-logistic model.}
#'   \item{`"normal"`}{`g(S) = -Phi^-1(S)`, the probit scale. With no
#'     interior knots this is a lognormal model.}
#' }
#'
#' @section The distributional parameters:
#' `mu` is `gamma_0` and every other coefficient is `gamma1`, `gamma2`,
#' and so on, all with identity links. There are `df + 1` of them.
#'
#' A formula on `mu` is the proportional-hazards (or -odds, or -probit)
#' model, because `gamma_0` shifts the whole curve. A formula on any
#' other coefficient is a TIME-VARYING effect, since that coefficient
#' multiplies a function of log time; this is what `flexsurv` spells
#' `anc =`.
#'
#' ```
#' frm(bf(t | cens(censored) ~ trt, gamma1 ~ trt),
#'     family = royston_parmar(df = 3), data = dat)
#' ```
#'
#' @section Knots:
#' The spline needs `df + 1` knots: two boundary knots and `df - 1`
#' interior ones. By default they are placed at equally spaced quantiles
#' of the log UNCENSORED times, which is Royston and Parmar's own rule
#' and flexsurv's default; the boundary knots are then the smallest and
#' largest log uncensored time.
#'
#' The quantiles cannot be taken until the response is in hand, so they
#' are taken at frame assembly through `family_finalize()`, and the
#' family object the fit carries has its knots baked in. `knots =` and
#' `bknots =` pin them instead, in which case they are taken as given
#' and no data is consulted. Both are on the LOG time scale, as
#' flexsurv's are.
#'
#' @section Monotonicity is not enforced, and the floor is not free:
#' The cumulative hazard has to increase, which means the spline's
#' derivative in `x` has to stay positive, and nothing in this
#' parameterization holds it there. flexsurv does not hold it either:
#' the model is fitted unconstrained and a fit whose spline turns over
#' inside the data range is over-parameterized, not something the
#' software should have prevented.
#'
#' What this family does differently is refuse to return `NaN` for it.
#' Where the derivative goes non-positive there is no hazard and the
#' true log density is `-Inf`, so the log density becomes a large finite
#' negative number instead: `NaN` stops the optimizer, and inside a
#' [frmtmb::mixture()] one `NaN` component poisons the log-sum-exp of
#' all of them.
#'
#' That floor is NOT inert when it is used. It keeps the fit alive and
#' it makes `logLik()` and `AIC()` a pseudo-likelihood rather than a
#' density: a 60 percent cure-fraction dataset has been measured
#' converging, without a warning, with 6 floored rows and a reported log
#' likelihood 3952 units away from the density's. Nothing in the fitted
#' object says so, because the family protocol has no hook that runs
#' when a fit finishes.
#'
#' [rp_floored()] is where that goes to be read, and it REFUSES by
#' default. Call it on every fit. `frm_curve()` and its two companions
#' call it for you, because the fitted curve of a floored fit is a floor
#' artifact too.
#'
#' @section Censoring, truncation, and the accuracy limit on log S:
#' The family declares both a density and a distribution function, so
#' `cens()` and `trunc()` both work, and right, left and interval
#' censoring all reach the likelihood through frmtmb's own machinery.
#' `cens()` takes frmtmb's coding: `0` (or `"none"`) is an observed
#' event and `1` (or `"right"`) is right censored, which is the OPPOSITE
#' of the `status` column of a `Surv()` object. Pass `1 - status`.
#'
#' Both are CONDITIONAL rather than unqualified. frmtmb forms a
#' right-censored contribution as `log(1 - F(y))` on the PROBABILITY
#' scale (`R/objective.R:100`), and core offers a family no
#' complementary log-CDF slot to hand back `log S` directly, so the
#' scored `log S` carries absolute error about
#' `.Machine$double.eps / S` whatever this family does internally.
#' Measured, by forming `F` for a given `-log S` and reading
#' `log(1 - F)` back (the three scales agree to every printed digit):
#'
#' \tabular{lll}{
#'   `-log S` \tab computed \tab absolute error \cr
#'   10   \tab -10        \tab 1.3e-13 \cr
#'   19.2 \tab -19.2      \tab 2.4e-10 \cr
#'   30   \tab -29.99983  \tab 1.7e-04 \cr
#'   36   \tab -34.94504  \tab 1.05    \cr
#'   40   \tab -35.12736  \tab 4.87
#' }
#'
#' `log S` is floored at -35.127363 whatever the model says, and past
#' `-log S` of 30 the term is FLAT: its gradient is exactly zero, so the
#' optimizer prices that row at a constant and fits the rest as if it
#' were free.
#'
#' How wrong the answer gets is a property of the DATA, not of the
#' family. A floored censored row contributes -35.127363 instead of its
#' own `-log S`, so the reported log likelihood is short by about
#' `-log S - 35` for each such row. That is thousands to tens of
#' thousands as soon as one censored time sits well past the event
#' times: two runs of the same 600-subject design, differing only in
#' seed, give 2.4e+03 and 2.166e+04. Both converged without a warning
#' and both put the treatment coefficient out by tens of percent, on
#' data flexsurv declines to fit at all.
#'
#' `-log S` is the cumulative hazard `H` on the `"hazard"` scale,
#' `log(1 + exp(eta))` on `"odds"` and `-log(Phi(-eta))` on `"normal"`.
#' [rp_floored()] checks it on every censored row and refuses past 19.2,
#' which is where `eps / S` passes 1e-8. The real fix is a core one, an
#' `lccdf` slot, and it is in `dev/spline-seam-proposal.md`.
#'
#' @param df Degrees of freedom: the number of interior knots plus one,
#'   which is Royston and Parmar's convention. `df = 1` is no interior
#'   knot at all and gives the Weibull, log-logistic or lognormal model
#'   that `scale` names. flexsurv counts the same spline with `k = df -
#'   1`.
#' @param knots Interior knots on the log time scale, given explicitly.
#'   `df` is then read off their number and any `df` argument is checked
#'   against them.
#' @param bknots The two boundary knots, on the log time scale. Defaults
#'   to the range of the log uncensored times.
#' @param scale `"hazard"`, `"odds"` or `"normal"`.
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Royston, P. and Parmar, M. K. B. (2002) Flexible parametric
#' proportional-hazards and proportional-odds models for censored
#' survival data. *Statistics in Medicine* 21, 2175-2197.
#'
#' @seealso [rp_floored()], which must be called on any fit whose data
#'   carry censoring; [frm_curve()] for reading the fitted log cumulative
#'   hazard off with a band.
#' @examples
#' set.seed(1)
#' n <- 300
#' dd <- data.frame(trt = rep(0:1, each = n / 2))
#' dd$t <- rweibull(n, shape = 1.4, scale = exp(1 - 0.5 * dd$trt))
#' dd$censored <- as.integer(dd$t > 3)
#' dd$t <- pmin(dd$t, 3)
#' fit <- frmtmb::frm(frmtmb::bf(t | cens(censored) ~ trt),
#'                    family = royston_parmar(df = 2), data = dd)
#' frmtmb::fixef(fit)$mu
#' @export
royston_parmar <- function(df = 3, knots = NULL, bknots = NULL,
                           scale = c("hazard", "odds", "normal")) {
  scale <- match.arg(scale)
  if (!is.null(knots)) {
    if (!is.numeric(knots) || anyNA(knots)) {
      stop("royston_parmar(knots = ) takes the interior knots as a ",
           "numeric vector on the LOG time scale (numeric(0) for none, ",
           "NULL to place them at quantiles), not ", class(knots)[1L],
           call. = FALSE)
    }
    knots <- sort(as.numeric(knots))
    if (!missing(df) && df != length(knots) + 1L) {
      stop("royston_parmar(): df = ", df, " and ", length(knots),
           " interior knots disagree. df is the number of interior ",
           "knots plus one, so drop df and let the knots decide",
           call. = FALSE)
    }
    df <- length(knots) + 1L
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1 ||
      df != round(df)) {
    stop("royston_parmar(df = ) must be one whole number of at least 1. ",
         "df = 1 is the spline with no interior knot, which is the ",
         "Weibull, log-logistic or lognormal model that scale names",
         call. = FALSE)
  }
  df <- as.integer(df)
  if (!is.null(bknots)) {
    if (!is.numeric(bknots) || length(bknots) != 2L || anyNA(bknots) ||
        bknots[1L] >= bknots[2L]) {
      stop("royston_parmar(bknots = ) takes the two boundary knots on ",
           "the log time scale, smallest first", call. = FALSE)
    }
    bknots <- as.numeric(bknots)
  }
  cfg <- list(df = df, knots = knots, bknots = bknots, scale = scale)
  sp_rp_family(cfg, allknots = if (!is.null(knots) && !is.null(bknots)) {
    c(bknots[1L], knots, bknots[2L])
  })
}

#' The family object, once the knots are known (or with `allknots =
#' NULL` while they are not).
#'
#' Building the same object twice, once at construction and once from
#' `family_finalize()`, is what keeps the knots on the family rather than
#' in an environment the link closures read at call time. The
#' drift-diffusion package learned that the hard way and this one starts
#' there.
#'
#' @noRd
sp_rp_family <- function(cfg, allknots) {
  ng <- cfg$df + 1L
  dpars <- c("mu", paste0("gamma", seq_len(ng - 1L)))
  links <- stats::setNames(rep(list("identity"), ng), dpars)
  gam_of <- function(dp) lapply(dpars, function(nm) dp[[nm]])
  ll <- function(y, dpars, aterms) {
    sp_rp_need_knots(allknots)
    x <- log(y)
    eta <- sp_rp_eta(sp_rp_basis(allknots, x), gam_of(dpars))
    detadx <- sp_rp_eta(sp_rp_dbasis(allknots, x), gam_of(dpars))
    # log f = log h + log S, and log h = eta + log(deta/dx) - log(t)
    # holds for every scale up to the link between eta and S
    lg <- log(sp_floor_pos(detadx)) - x
    switch(cfg$scale,
      # exp(eta) overflows past eta = 709 and an Inf in the log density
      # is a NaN gradient one line later. The cumulative hazard is
      # capped at exp(30), which is 1.1e13: a survival probability of
      # exp(-1.1e13) is not a number any fit reports, and the cap is
      # exact in double precision everywhere below it.
      hazard = eta + lg - exp(sp_cap(eta, 30)),
      odds = eta + lg - 2 * sp_log1pexp(eta),
      normal = RTMB::dnorm(eta, 0, 1, log = TRUE) + lg)
  }
  cdf <- function(q, dpars, aterms) {
    sp_rp_need_knots(allknots)
    x <- log(q)
    eta <- sp_rp_eta(sp_rp_basis(allknots, x), gam_of(dpars))
    sp_squeeze(switch(cfg$scale,
      hazard = -expm1(-exp(eta)),
      odds = 1 - 1 / (1 + exp(eta)),
      normal = RTMB::pnorm(eta, 0, 1)))
  }
  custom_family(
    "royston_parmar",
    dpars = dpars,
    links = links,
    lpdf = ll,
    lcdf = cdf,
    valid_y = function(y, aterms) sp_rp_valid_y(y),
    family_finalize = function(fam, y, aterms) {
      sp_rp_family(cfg, sp_rp_knots(cfg, y, aterms))
    },
    init_dpars = sp_rp_init(cfg, dpars),
    type = "continuous",
    post = list(
      # mu is gamma0, a spline coefficient, and reporting it as a fitted
      # response is exactly the mistake this family invites: the link is
      # the identity, so nothing downstream would notice. The mean of a
      # Royston-Parmar survival time is an integral of the fitted
      # survival function with no closed form, and the censored rows do
      # not identify its upper tail. core's cox() refuses the same
      # question for the same reason.
      mean_fn = function(dpars, aterms) {
        stop("royston_parmar: a survival time has no mean on the ",
             "response scale here. mu is gamma0, the intercept of a ",
             "spline in log time, not a fitted value, and the mean ",
             "survival time is an integral over a tail the censored ",
             "rows do not identify. predict(type = \"link\", dpar = ) ",
             "gives any spline coefficient, and frm_curve() reads the ",
             "fitted log cumulative hazard off with a band",
             call. = FALSE)
      }
    ),
    sim = function(dpars, aterms, n) {
      sp_rp_sim(cfg, allknots, dpars, n)
    })
}

#' The knots are decided at frame assembly, so a density evaluated
#' before that has none. One template, called from both densities.
#'
#' @noRd
sp_rp_need_knots <- function(allknots) {
  if (is.null(allknots)) {
    stop("royston_parmar(): this family object has no knots yet. They ",
         "are quantiles of the log uncensored times, so they are found ",
         "when frm() assembles the model frame; a density cannot be ",
         "evaluated before that. Give knots = and bknots = to pin them ",
         "at construction instead", call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
sp_rp_valid_y <- function(y) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("royston_parmar(): the response is a survival TIME, so it must ",
         "be finite and strictly positive; the spline is a function of ",
         "log(t) and log of a non-positive time is not a number",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Knot placement, at frame assembly.
#'
#' Royston and Parmar's rule and flexsurv's default: equally spaced
#' quantiles of the log times of the observed EVENTS. Censored rows carry
#' no event time, so including them would place knots where nothing was
#' seen to happen.
#'
#' @noRd
sp_rp_knots <- function(cfg, y, aterms) {
  if (!is.null(cfg$knots) && !is.null(cfg$bknots)) {
    return(c(cfg$bknots[1L], cfg$knots, cfg$bknots[2L]))
  }
  x <- log(as.numeric(y))
  cen <- aterms[["cens"]]
  if (!is.null(cen)) {
    ev <- which(cen == 0)
    if (length(ev) >= cfg$df + 1L) x <- x[ev]
  }
  if (length(unique(x)) < cfg$df + 1L) {
    stop("royston_parmar(): the spline needs df + 1 = ", cfg$df + 1L,
         " knots and the data offer only ", length(unique(x)),
         " distinct uncensored log times to place them at. Lower df",
         call. = FALSE)
  }
  bk <- cfg$bknots %||% range(x)
  ik <- cfg$knots
  if (is.null(ik)) {
    ik <- if (cfg$df > 1L) {
      unname(stats::quantile(x, seq(0, 1, length.out = cfg$df + 1L)))[
        seq_len(cfg$df - 1L) + 1L]
    } else {
      numeric(0)
    }
  }
  kn <- c(bk[1L], ik, bk[2L])
  if (any(diff(kn) <= 0)) {
    stop("royston_parmar(): the knots are not strictly increasing (",
         paste(format(kn, digits = 4), collapse = ", "),
         "). Tied log event times put two quantiles at one place; lower ",
         "df, or give knots explicitly", call. = FALSE)
  }
  kn
}

#' Starting values: the model with every interior coefficient at zero,
#' which is exactly the Weibull, log-logistic or lognormal member of the
#' family, fitted by moments on the log event times.
#'
#' @noRd
sp_rp_init <- function(cfg, dpars) {
  first <- function(y, aterms) {
    p <- sp_rp_moments(y, aterms, cfg$scale)
    p[["g0"]]
  }
  second <- function(y, aterms) {
    p <- sp_rp_moments(y, aterms, cfg$scale)
    p[["g1"]]
  }
  out <- list(mu = first)
  for (nm in dpars[-1L]) out[[nm]] <- function(y, aterms) 0
  if (length(dpars) > 1L) out[["gamma1"]] <- second
  out
}

#' @noRd
sp_rp_moments <- function(y, aterms, scale) {
  x <- log(as.numeric(y))
  cen <- aterms[["cens"]]
  if (!is.null(cen) && any(cen == 0)) x <- x[cen == 0]
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  m <- mean(x)
  # extreme-value moments for the hazard scale; the other two scales use
  # the same slope, which is the right order of magnitude for all three
  g1 <- switch(scale, hazard = 1.2825 / s, odds = 1.8138 / s, 1 / s)
  list(g0 = -g1 * (m + 0.5772 / g1), g1 = g1)
}

#' Draw survival times by inverting the fitted survival function.
#'
#' Vectorized bisection on log time rather than a per-row root find: the
#' survival function is monotone in `t` whenever the model is sane, the
#' bracket is set from the knots, and 200 halvings reach the double
#' precision floor on any bracket a survival time lives in.
#'
#' @noRd
sp_rp_sim <- function(cfg, allknots, dpars, n) {
  u <- stats::runif(n)
  target <- switch(cfg$scale,
    hazard = log(-log(u)),
    odds = log(1 / u - 1),
    normal = -stats::qnorm(u))
  gam <- lapply(seq_along(dpars), function(j) {
    v <- dpars[[j]]
    if (length(v) == 1L) rep(v, n) else v
  })
  lo <- rep(allknots[1L] - 20, n)
  hi <- rep(allknots[length(allknots)] + 20, n)
  for (i in seq_len(200L)) {
    mid <- (lo + hi) / 2
    e <- sp_rp_eta(sp_rp_basis(allknots, mid), gam)
    up <- e < target
    lo <- ifelse(up, mid, lo)
    hi <- ifelse(up, hi, mid)
  }
  exp((lo + hi) / 2)
}
