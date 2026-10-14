# Starting values for a location-scale fit from the least-squares fit
# of its location, so that the scale starts at the spread the location
# leaves and not at the spread of the response.

#' Starting values for `mu` and the scale from a least-squares fit
#'
#' The family initializers start `mu` at an intercept read off the
#' response and the scale at the response's own spread (`sd(y)` for
#' gaussian()). When a covariate explains most of the response, that
#' spread is far larger than the residual one: on `y = 1 + 1e4 z1 + b_g
#' + e` (n 101000, 50 columns, residual sd 1) sigma started at 1e4 and
#' every slope at 0, the first steps put the group sd at its boundary,
#' where the gradient on the log scale vanishes, and the fit stopped
#' 9107.66 below lme4 at code 0 with no warning (dev/optrap-p2-wide1.R,
#' dev/optstart-findings.md). Started at the coefficients of a
#' least-squares fit and the spread of its residuals, the same fit
#' reaches lme4's maximum.
#'
#' `start_ls_plan()` returns, per response it applies to, the starting
#' coefficients of `mu` (`beta`, on the link scale, `0` for an aliased
#' column), the starting value of each scale parameter on its natural
#' scale (`scale`, by dpar), the residuals (`resid`, for an initializer
#' that reads them, skew_normal()'s alpha) and the fit itself (`fit`,
#' keyed by the `mu` predictor, which the response unit reuses:
#' predictor_units()). An empty list where it applies to nothing.
#'
#' It applies to a response only where every condition below holds.
#' Anywhere else `mu` and the scale start as before this rule, bit for
#' bit. The group sds start at the response unit, which changed for an
#' `mi()` response (imputed rows are not data) and for
#' zero_inflated_asym_laplace() (its zeros are not data), so there they
#' and the units differ from 0.70.0's even where this rule does not
#' apply.
#'
#' * The family is in `start_ls_rules` and `mu` has the link the rule is
#'   written for. A family with no scale (poisson(), binomial()), an
#'   ordinal or categorical family, a mixture and a family from an
#'   extension or a user are not in it.
#' * `mu` is a linear predictor with a design of its own (not nonlinear,
#'   not a constant, not tied to another), and with something to fit
#'   besides an intercept: a covariate or an offset. An intercept-only
#'   `mu` starts where its initializer puts it already.
#' * Every scale parameter the rule sets is estimated (not held at a
#'   constant, as `se()` and `sigma = 2` do).
#' * The response carries no addition term but `weights()`, `subset()`,
#'   `index()` and `mi()`. Censored and truncated values are not the
#'   values a least-squares fit reads, and `se()` sets the residual
#'   spread itself.
#' * `mu` has no smooth (`s()`, `t2()`) and no `gp()` term. Their
#'   penalized parts enter as random effects, and least squares on the
#'   rest of the design takes the signal they share with it: on
#'   test-v11.R's scalar-on-function regression the smoothing sd then
#'   ran to its boundary and the fit stopped 0.18 below mgcv's maximum at
#'   code 0, and test-smooth-population.R's fit 7.02 below 0.70.0's
#'   (dev/optstart-sof.R, dev/optstart-findings.md).
#' * The design's condition number, columns at unit length, is at most
#'   1e5 (see below in start_ls_one()).
#'
#' A `start` or a prior does not turn the rule off. Each entry `start`
#' names is the user's and is written over the result (make_start()),
#' so a partly named `start` leaves the others at these defaults, as
#' `frm(start =)` documents; a prior on an ordinary coefficient moves no
#' starting value (test-par-template.R). A bound a prior sets can put a
#' least-squares value outside the box, and the optimizer starts from
#' its projection into the box (nlminb() does that itself: a start of 10
#' under `upper = 5` reaches the optimum 3).
#'
#' @noRd
start_ls_plan <- function(frame) {
  out <- list()
  for (lp in frame[["linpreds"]]) {
    if (!identical(lp[["dpar"]], "mu")) next
    key <- as.character(lp[["resp"]])
    if (!is.null(out[[key]])) next
    pl <- tryCatch(start_ls_one(frame, lp), error = function(e) NULL)
    if (!is.null(pl)) out[[key]] <- pl
  }
  out
}

#' The plan of start_ls_plan() for one response, given its `mu`
#' predictor, or `NULL` where the rule does not apply.
#'
#' @noRd
start_ls_one <- function(frame, lp) {
  resp <- lp[["resp"]]
  fam <- frame[["spec"]]$responses[[resp]]$family
  fname <- fam[["family"]]
  if (!is.character(fname) || length(fname) != 1L) return(NULL)
  rule <- start_ls_rules[[fname]]
  if (is.null(rule)) return(NULL)
  if (!is.null(lp[["constant"]]) || !is.null(lp[["equate"]]) ||
        !is.null(lp[["nl_body"]])) return(NULL)
  if (!identical(lp[["link"]]$name %||% "", rule$link)) return(NULL)
  X <- lp[["X"]]
  if (!design_ok(X) || length(lp[["idx"]]) != ncol(X)) return(NULL)
  if (length(lp[["smooths"]]) || length(lp[["gps"]])) return(NULL)
  cn <- colnames(X) %||% character(ncol(X))
  icpt <- match("(Intercept)", cn)
  off <- lp[["offset"]]
  has_off <- is.numeric(off) && any(off != 0)
  if (!has_off && all(cn == "(Intercept)")) return(NULL)
  # the intercept carries the shift between the mean of the link-scale
  # response and the link of its mean
  if (isTRUE(rule$needs_intercept) && is.na(icpt)) return(NULL)
  at <- frame[["aterm_values"]][[resp]]
  if (length(setdiff(names(at), c("weights", "subset", "index", "mi")))) {
    return(NULL)
  }
  y <- frame[["y"]][[resp]]
  if (!is.numeric(y) || !is.null(dim(y))) return(NULL)
  scale_lps <- list()
  for (l in frame[["linpreds"]]) {
    if (!identical(l[["resp"]], resp) ||
          !l[["dpar"]] %in% rule$dpars) next
    if (!is.null(l[["constant"]]) || !is.null(l[["equate"]]) ||
          !is.null(l[["nl_body"]])) return(NULL)
    scale_lps[[l[["dpar"]]]] <- l
  }
  if (!all(rule$dpars %in% names(scale_lps))) return(NULL)
  w <- at[["weights"]]
  if (!(is.numeric(w) && length(w) == length(y))) w <- NULL
  # an identity-link mu reads the response the unit reads, so that the
  # two share one solve; that input is NULL exactly where the unit's
  # own conditions fail (an unestimated sigma among them)
  if (identical(rule$link, "identity")) {
    inp <- predictor_unit_input(frame, lp)
    if (is.null(inp)) return(NULL)
  } else {
    z <- rule$z(y)
    ok <- is.finite(z) & response_observed(frame, resp)
    if (!is.null(w)) ok <- ok & w > 0
    inp <- list(z = z, ok = ok, w = w)
  }
  # the working weights of a rule that has them (negbinomial()'s first
  # IRLS step) multiply the likelihood weights
  lw <- inp$w
  if (is.function(rule$wt)) {
    ww <- rule$wt(y)
    lw <- if (is.null(lw)) ww else lw * ww
  }
  fit <- predictor_ls(lp, inp$z, inp$ok, lw)
  if (is.null(fit) || length(fit$coef) != ncol(X)) return(NULL)
  # An ill-conditioned design is not reached from either start, and the
  # two stops are judged differently: on test-se-flags.R's degree-8 raw
  # polynomial with (1 | g) (seed 3) the fit stopped 4.03 short from
  # the old start, where the standard-error check withholds the
  # coefficients as "illcond", and 4.10 short from this one with a
  # large-gradient warning, under which the check does not run and the
  # coefficients' SEs were 0.05 to 0.12 of the true ones
  # (dev/optstart-poly8.R). The bound and the number are newton_polish()'s
  # (lp_design_kappa()), which takes such fits on from the stop.
  if (!(lp_design_kappa(lp) <= 1e5)) return(NULL)
  r <- fit$resid
  rows <- fit$rows[fit$keep]
  wr <- if (is.null(inp$w)) rep(1, length(rows)) else inp$w[rows]
  # an exact fit leaves rounding, not a scale (predictor_resid_sd())
  if (!(sqrt(sum(wr * r^2) / sum(wr)) > 1e-10 * fit$mag)) return(NULL)
  beta <- fit$coef
  beta[is.na(beta)] <- 0
  eta <- inp$z[rows] - r
  sc <- rule$scale(r = r, w = wr, y = y[rows], eta = eta)
  if (!is.list(sc)) return(NULL)
  for (nm in names(sc)) {
    v <- sc[[nm]]
    if (nm == ".shift") next
    if (length(v) != 1L || !is.finite(v) || v <= 0) return(NULL)
  }
  if (!is.null(sc[[".shift"]])) {
    if (!is.finite(sc[[".shift"]])) return(NULL)
    beta[icpt] <- beta[icpt] + sc[[".shift"]]
    sc[[".shift"]] <- NULL
  }
  if (!all(is.finite(beta))) return(NULL)
  fits <- list()
  # the unit reads the same fit only where it solves the same problem
  if (identical(rule$link, "identity") && is.null(rule$wt)) {
    fits[[lp_key_of(frame, lp)]] <- fit
  }
  list(beta = beta, scale = sc, resid = r, fits = fits)
}

#' The key of a linear predictor in `frame$linpreds`.
#'
#' @noRd
lp_key_of <- function(frame, lp) {
  for (k in names(frame[["linpreds"]])) {
    l <- frame[["linpreds"]][[k]]
    if (identical(l[["resp"]], lp[["resp"]]) &&
          identical(l[["dpar"]], lp[["dpar"]])) return(k)
  }
  NA_character_
}

#' Weighted root mean square, the maximum-likelihood normal sd of
#' residuals whose mean the fit has taken out.
#'
#' @noRd
ls_rms <- function(r, w) sqrt(sum(w * r^2) / sum(w))

#' The shape `k` with `trigamma(k) = v`, the variance of `log y` for a
#' gamma `y` of shape `k`: the root of the two-term expansion
#' `1 / k + 1 / (2 k^2)`, then Newton steps.
#'
#' @noRd
trigamma_inverse <- function(v) {
  if (!is.finite(v) || v <= 0) return(NA_real_)
  k <- (1 + sqrt(1 + 2 * v)) / (2 * v)
  for (i in 1:20) {
    step <- (trigamma(k) - v) / psigamma(k, 2)
    k1 <- k - step
    if (!is.finite(k1)) return(NA_real_)
    if (k1 <= 0) k1 <- k / 2
    if (abs(k1 - k) <= 1e-12 * k) return(k1)
    k <- k1
  }
  k
}

#' The families start_ls_plan() reaches, each with the link of `mu` it is
#' written for, the scale parameters it sets (`dpars`) and `scale`, the
#' relation between them and the least-squares residuals `r` (weights
#' `w`, link-scale fitted values `eta`, response `y`, all on the rows
#' the fit used). `scale` returns natural-scale values by dpar, and
#' `.shift`, a correction to the intercept, where the mean of the
#' link-scale response is not the link of the mean. A family whose `mu`
#' link is not the identity has `z`, the link-scale response (`NA` where
#' a row does not inform it), and may have `wt`, working weights.
#' dev/optstart-relations.R measures each relation on simulated data.
#'
#' @noRd
start_ls_rules <- local({
  rms_sigma <- function(r, w, y, eta) list(sigma = ls_rms(r, w))
  # median absolute deviation, unscaled
  mad0 <- function(r) stats::median(abs(r - stats::median(r)))
  laplace_sigma <- function(r, w, y, eta) {
    # at quantile 0.5 the density is Laplace with scale 2 sigma, whose
    # maximum-likelihood scale is the mean absolute deviation from the
    # median
    list(sigma = sum(w * abs(r - stats::median(r))) / sum(w) / 2)
  }
  list(
    gaussian = list(link = "identity", dpars = "sigma", scale = rms_sigma),
    lognormal = list(link = "identity", dpars = "sigma", scale = rms_sigma),
    shifted_lognormal = list(link = "identity", dpars = "sigma",
                             scale = rms_sigma),
    hurdle_lognormal = list(link = "identity", dpars = "sigma",
                            scale = rms_sigma),
    student = list(link = "identity", dpars = "sigma",
                   scale = function(r, w, y, eta) {
                     # the scale of a t with the 4 degrees of freedom
                     # nu starts at, from its median absolute deviation:
                     # robust to the tails the family is there for
                     list(sigma = mad0(r) / stats::qt(0.75, 4))
                   }),
    # sigma and alpha come from the family's initializers, which read
    # the residuals
    skew_normal = list(link = "identity", dpars = character(0),
                       scale = function(r, w, y, eta) list()),
    exgaussian = list(link = "identity", dpars = c("sigma", "beta"),
                      scale = function(r, w, y, eta) {
                        # half the residual spread to each, the
                        # divisor of the family's initializer (sd(y) / 2
                        # to each). Where a group sd 30 times sigma holds
                        # the residual, the density's gradient goes NaN
                        # from larger starts: rms / sqrt(2) failed 10 of
                        # 10 such fits, rms / 2 fails 5 of 10 as 0.70.0
                        # did, and a third-moment split failed 9 of 10
                        # (dev/optstart-rev-exgfix.R,
                        # dev/optstart-exg-arms.R)
                        s <- ls_rms(r, w) / 2
                        list(sigma = s, beta = s)
                      }),
    asym_laplace = list(link = "identity", dpars = "sigma",
                        scale = laplace_sigma),
    zero_inflated_asym_laplace = list(link = "identity", dpars = "sigma",
                                      scale = laplace_sigma),
    huber = list(link = "identity", dpars = "sigma",
                 scale = function(r, w, y, eta) {
                   s <- stats::mad(r)
                   list(sigma = if (s > 0) s else ls_rms(r, w))
                 }),
    Gamma = list(link = "log", dpars = "shape", needs_intercept = TRUE,
                 z = function(y) log(y),
                 scale = function(r, w, y, eta) {
                   # var(log y) = trigamma(shape), and
                   # E log y = log mu + digamma(shape) - log(shape)
                   k <- trigamma_inverse(ls_rms(r, w)^2)
                   list(shape = k, .shift = log(k) - digamma(k))
                 }),
    weibull = list(link = "log", dpars = "shape", needs_intercept = TRUE,
                   z = function(y) log(y),
                   scale = function(r, w, y, eta) {
                     # log y is a Gumbel of scale 1 / shape:
                     # var(log y) = pi^2 / (6 shape^2), and
                     # E log y = log mu - lgamma(1 + 1 / shape)
                     #   - euler / shape
                     k <- pi / (sqrt(6) * ls_rms(r, w))
                     list(shape = k,
                          .shift = lgamma(1 + 1 / k) +
                            -digamma(1) / k)
                   }),
    beta = list(link = "logit", dpars = "phi",
                z = function(y) stats::qlogis(y),
                scale = function(r, w, y, eta) {
                  # var(y) = mu (1 - mu) / (1 + phi)
                  mu <- stats::plogis(eta)
                  phi <- sum(w * mu * (1 - mu)) / sum(w * (y - mu)^2) - 1
                  list(phi = max(phi, 0.5))
                }),
    negbinomial = list(link = "log", dpars = "shape", needs_intercept = TRUE,
                       # glm()'s first IRLS step from mustart = y + 0.1
                       z = function(y) log(y + 0.1) - 0.1 / (y + 0.1),
                       wt = function(y) y + 0.1,
                       scale = function(r, w, y, eta) {
                         # var(y) = mu + mu^2 / shape, with the
                         # intercept moved so the fitted mean is the
                         # mean
                         mu0 <- exp(eta)
                         s <- log(sum(w * y) / sum(w * mu0))
                         mu <- mu0 * exp(s)
                         ex <- sum(w * ((y - mu)^2 - mu))
                         k <- if (ex > 0) sum(w * mu^2) / ex else 10
                         list(shape = max(k, 0.1), .shift = s)
                       })
  )
})

#' Which rows of a response are data: `FALSE` where `mi()` imputes the
#' value. The frame stores an imputed value as 0, so without this the
#' zeros read as observations: on brms_missings' `bmi | mi()` the unit
#' came from a response nine of whose 25 rows were 0 (lane remeasure,
#' R1). `resp` is the response's key in `frame$y`, a name or a position.
#'
#' @noRd
response_observed <- function(frame, resp) {
  y <- frame[["y"]][[resp]]
  n <- NROW(y)
  nm <- if (is.character(resp)) resp else names(frame[["y"]])[resp]
  mm <- if (length(nm) == 1L && !is.na(nm)) frame[["mi_map"]][[nm]]
  out <- rep(TRUE, n)
  if (is.null(mm)) return(out)
  if (!is.null(mm[["obs"]])) {
    out[] <- FALSE
    out[mm[["obs"]]] <- TRUE
  } else {
    out[mm[["rows"]]] <- FALSE
  }
  out
}
