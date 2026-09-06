# A penalized coefficient block whose VALUE a nonlinear body consumes.
#
# Every other penalized smooth in frmtmb ends up multiplied by Z: the
# frame splits an mgcv basis into a null space that joins X and a wiggly
# part that becomes a random-effect block, and the objective adds `Z b`
# to a linear predictor. That is exactly right for `s(x)` and useless
# for a body that needs the spline's VALUE inside an expression -
# `exp(amp) * ps(age + shift)` - because the argument of the spline is
# itself a function of parameters, so there is no fixed Z to build.
#
# `ps()` keeps the same mixed-model split and moves the evaluation onto
# the tape. The basis is a divided difference of truncated powers, which
# is branch free: RTMB refuses `pmax()` on an advector ("Comparison is
# generally unsafe for AD types") and exports no `CondExp`, so the
# recursive Cox-de Boor form is unavailable and this one is not a
# convenience but the only construction there is.

#' A penalized spline whose value a nonlinear body consumes
#'
#' `ps()` declares a penalized coefficient block inside the body of a
#' nonlinear formula (`bf(..., nl = TRUE)` or [nlf()]) and evaluates to
#' the spline's VALUE at `expr`. It is not added to a linear predictor,
#' which is what separates it from [mgcv::s()]: the body may multiply
#' it, exponentiate it, or pass it on.
#'
#' @param expr The expression the spline is evaluated at. It may name
#'   nonlinear parameters, so the evaluation point moves with the fit -
#'   `ps(age + shift)` is a curve each subject sees shifted along its own
#'   time axis.
#' @param k Number of basis functions (default 10). At least
#'   `max(degree + 1, 4)` and at most 50; see Accuracy.
#' @param degree Degree of the B-spline pieces (default 3, cubic).
#' @param pad Fraction of the data-time range of `expr` added at each end
#'   before the knots are placed (default 0.1). The basis is exactly zero
#'   outside its knot span, so an evaluation point that the parameters
#'   push past the padded range reads as a curve value of zero rather
#'   than as an extrapolation. Widen `pad` when the transformation is
#'   large; the fit reports the coverage it achieved.
#' @param center If `TRUE` (default), the coefficients are constrained to
#'   sum to zero. B-splines are a partition of unity, so the constraint
#'   removes the curve's overall level, which is otherwise confounded
#'   with an intercept in the body.
#' @return `ps()` is a term, not a function to call. Evaluating one
#'   outside a nonlinear body is an error.
#' @section Where the pieces go:
#' The second-difference penalty `S` is eigendecomposed and the basis is
#' split the way [mgcv::smooth2random()] splits `s()`: the null space of
#' `S` joins the fixed coefficients (`beta`, or `betad` for a
#' distributional parameter), and the range space becomes one
#' random-effect block with a single variance in `theta`. The smoothing
#' parameter is that variance's inverse, so smoothness is estimated
#' jointly with every other variance component rather than chosen.
#'
#' [VarCorr()] and [ranef()] report the block under the term's own label.
#' Its standard deviation is a smoothing parameter, not a subject
#' effect.
#'
#' @section The knots are frozen:
#' The knot vector is built once, from the range of `expr` evaluated on
#' the model frame with every nonlinear parameter set to zero, padded by
#' `pad`. It is then frozen on the frame and reused on `newdata`, which
#' is [poly()]'s rule and for the same reason: a basis rebuilt from new
#' data is a different basis, and the coefficients would not mean what
#' they were fitted to mean.
#'
#' Boundary knots are DISTINCT and lie outside the padded range. The
#' divided-difference form degenerates on the repeated boundary knots
#' `splines::splineDesign()` normally uses, so `ps()` places them apart
#' rather than exposing the choice.
#'
#' This is NOT the rule of D'Alessandro, Thoresen and Sorensen (2026),
#' section 2.3.1, and `ps()` is not a reimplementation of it. That paper
#' rescales the spline argument to `[0, 1]` between bounds
#' `min(t) - 3 sigma_b` and `max(t) + 3 sigma_b`, which are smooth
#' functions of the CURRENT variance estimates, so its knots move with
#' the fit. `ps()` freezes instead, for two reasons. Moving knots make
#' the penalty matrix `S` and its eigendecomposition functions of the
#' parameters, so the split into fixed and random coefficients would
#' have to be redone inside the objective; and a basis that is rebuilt
#' between the fit and the prediction is a different basis, which is
#' the rule `poly()` follows and the reason `pad =` exists. On the
#' paper's own application the two agree closely enough that six of its
#' seven reported parameters come back at the published precision, but
#' they are related constructions rather than one.
#'
#' @section Accuracy:
#' The basis agrees with [splines::splineDesign()] to 3.1e-13 absolute at
#' `k = 12` and 1.8e-11 at `k = 40`, with an AD input giving the same
#' values bit for bit and the taped derivative matching
#' `splineDesign(derivs = 1)` to 4.5e-13 and 4.5e-11 respectively. The
#' error grows like `k^(degree - 1)` because a divided difference of
#' truncated powers cancels terms of order `(range / spacing)^degree`.
#' `k` above 50 is refused for that reason; a curve that needs more than
#' 50 basis functions needs a different basis, not a looser tolerance.
#'
#' @section What the rest of the package does with it:
#' `predict()` and [simulate()] work: the body is re-evaluated at the new
#' data or at the drawn coefficients, and the frozen basis is evaluated
#' at whatever argument comes out.
#'
#' `predict(se.fit = TRUE)` stays refused for a nonlinear predictor.
#' [frm_lp_basis()] is the route: it tapes the body and returns
#' `d eta / d coef` as a Jacobian, which is what a delta method over a
#' warped curve needs.
#'
#' `REML = TRUE`, `quadrature = TRUE`, `frmtmb_control(profile = TRUE)`
#' and multivariate ([mvbf()]) models are refused by name. The first
#' three integrate out something a `ps()` block has already put in the
#' Laplace approximation; the fourth is untested.
#'
#' The importance correction ([frm()]'s `importance =`) is refused, and
#' not by anything `ps()` declares: core refuses the correction for
#' EVERY nonlinear predictor, because a body mixes parameter values with
#' raw data columns and the corrected objective evaluates the predictor
#' once per draw.
#'
#' @references
#' D'Alessandro, M., Thoresen, M. and Sorensen, O. (2026). A
#' Semiparametric Nonlinear Mixed Effects Model with Penalized Splines
#' Using Automatic Differentiation. arXiv:2603.11728.
#'
#' Wood, S. N. (2004). Stable and efficient multiple smoothing parameter
#' estimation for generalized additive models. Journal of the American
#' Statistical Association 99, 673-686.
#' @seealso [frm_lp_basis()] for standard errors on a curve a `ps()`
#'   block builds, [frmtmb-extension-api]
#' @examples
#' set.seed(1)
#' n_id <- 40
#' d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
#' sh <- stats::rnorm(n_id, 0, 0.05)[d$id]
#' d$y <- 2 + sin(2 * pi * (d$t + sh)) + stats::rnorm(nrow(d), 0, 0.1)
#' fit <- frm(bf(y ~ lev + ps(t + shift, k = 8),
#'               lev ~ 1, shift ~ 0 + (1 | id), nl = TRUE),
#'            data = d)
#' fixef(fit)
#' @export
ps <- function(expr, k = 10, degree = 3, pad = 0.1, center = TRUE) {
  stop("ps() is a term, not a function: it declares a penalized ",
       "coefficient block and is only meaningful inside the body of a ",
       "nonlinear formula (bf(..., nl = TRUE) or nlf()). It was ",
       "evaluated as ordinary R code instead", call. = FALSE)
}

#' The truncated positive part, branch free.
#'
#' RTMB refuses `pmax()` on an advector and exports no `CondExp`, so
#' `0.5 * (e + abs(e))` is the only spelling of `(e)_+` that tapes.
#'
#' @noRd
ps_pos <- function(e) 0.5 * (e + abs(e))

#' `p^e` for a small whole `e`, by repeated multiplication.
#'
#' `^` on an advector tapes a `pow` node with a log inside it, which is
#' `NaN` at the exactly-zero values the truncated part produces on more
#' than half of every column.
#'
#' @noRd
ps_pow <- function(p, e) {
  out <- p
  for (i in seq_len(e - 1L)) out <- out * p
  out
}

#' Uniform knot vector: `nb` basis functions of order `ord` whose joint
#' support covers `[lo, hi]`, with distinct knots throughout.
#'
#' @noRd
ps_knots <- function(lo, hi, nb, ord) {
  h <- (hi - lo) / (nb + 1L - ord)
  lo + (seq_len(nb + ord) - ord) * h
}

#' Divided-difference weights of the B-spline basis:
#' `B_i(x) = sum_j W[i, j] * (t[i + j - 1] - x)_+^(ord - 1)`.
#'
#' The weights are constants, so all of the arithmetic that depends on
#' the knots happens once at frame time and the tape carries only the
#' truncated powers.
#'
#' @noRd
ps_weights <- function(knots, ord) {
  nb <- length(knots) - ord
  W <- matrix(0, nb, ord + 1L)
  for (i in seq_len(nb)) {
    tt <- knots[i + seq_len(ord + 1L) - 1L]
    for (j in seq_len(ord + 1L)) {
      W[i, j] <- (knots[i + ord] - knots[i]) / prod(tt[j] - tt[-j])
    }
  }
  W
}

#' The `n x nb` basis matrix at `x`, which may be numeric or an
#' advector.
#'
#' @noRd
ps_design <- function(x, knots, ord, W) {
  nb <- nrow(W)
  P <- lapply(seq_along(knots),
              function(l) ps_pow(ps_pos(knots[l] - x), ord - 1L))
  cols <- vector("list", nb)
  for (i in seq_len(nb)) {
    b <- W[i, 1L] * P[[i]]
    for (j in seq_len(ord + 1L)[-1L]) b <- b + W[i, j] * P[[i + j - 1L]]
    cols[[i]] <- b
  }
  do.call(cbind, cols)
}

#' The curve `sum_i cf[i] * B_i(x)`, summed per basis function rather
#' than per knot.
#'
#' The two orders are algebraically the same and numerically are not: a
#' divided difference cancels terms of order `(range / spacing)^degree`,
#' and collapsing the sum over knots first performs that cancellation
#' across the whole knot vector instead of inside each basis function's
#' `ord + 1` window. Grouping by basis keeps the measured accuracy.
#'
#' @noRd
ps_value <- function(x, knots, ord, W, cf) {
  nb <- nrow(W)
  P <- lapply(seq_along(knots),
              function(l) ps_pow(ps_pos(knots[l] - x), ord - 1L))
  val <- 0
  for (i in seq_len(nb)) {
    b <- W[i, 1L] * P[[i]]
    for (j in seq_len(ord + 1L)[-1L]) b <- b + W[i, j] * P[[i + j - 1L]]
    val <- val + cf[i] * b
  }
  val
}

#' The spline coefficients of one `ps()` term, from the null-space fixed
#' coefficients and the penalized random-effect block.
#'
#' Spelled as a loop over columns rather than as `U0 %*% b0`: the factors
#' are a numeric matrix and an advector, and a column-at-a-time
#' accumulation needs no matrix-multiply method for the mixed pair.
#'
#' @noRd
ps_coefs <- function(pt, b0, u) {
  cf <- pt[["U0"]][, 1L] * b0[1L]
  for (j in seq_along(b0)[-1L]) cf <- cf + pt[["U0"]][, j] * b0[j]
  for (j in seq_along(u)) cf <- cf + pt[["Us"]][, j] * u[j]
  cf
}

#' The per-term closures a nonlinear body finds by ordinary lookup.
#'
#' They go in `ev`, the per-call evaluation frame, and NOT in `nl_env`:
#' `nl_env` is fixed at parse time, before any coefficient exists, so a
#' closure living there would capture a stale environment and could
#' never see the current draw of the penalized block. `ev` is rebuilt
#' per objective call and already holds advectors, so a closure placed
#' there receives its argument already evaluated as an advector and
#' returns one.
#'
#' `check` is `FALSE`, `TRUE`, or a collector environment. `predict()`
#' passes `TRUE` and the closure warns where it stands, which is right
#' there because it evaluates the body exactly once. `frm_lp_basis()`
#' passes an environment, and the reason is measured rather than
#' defensive: that route enters this closure TWICE per term, once for
#' its own off-tape span pass (`R/predict.R`, before `MakeTape()`) and
#' once for the tape build. Warning in place would double-fire on every
#' call. The collector keys on the term, so the second pass overwrites
#' rather than repeats; see `ps_span_flush()`.
#'
#' RTMB itself enters once, so a reader chasing a double-fire should
#' look at the off-tape pass and not at the tape.
#'
#' @noRd
ps_env <- function(lp, pars, bvec, check = FALSE) {
  out <- list()
  for (pt in lp[["ps_terms"]] %||% list()) {
    # `eval_dpars(b = NULL)` drops the random-effect contribution, and a
    # penalized block IS a random effect: dropping it leaves the null
    # space, which is the unpenalized part of the curve
    u <- if (is.null(bvec)) rep(0, pt[["n_pen"]]) else bvec[pt[["c_idx"]]]
    cf <- ps_coefs(pt, pars[[pt[["par"]]]][pt[["beta_idx"]]], u)
    out[[pt[["fname"]]]] <- local({
      knots <- pt[["knots"]]
      ord <- pt[["ord"]]
      W <- pt[["W"]]
      cff <- cf
      term <- pt
      chk <- check
      function(x) {
        if (!isFALSE(chk)) ps_span_warning(x, term, chk)
        ps_value(x, knots, ord, W, cff)
      }
    })
  }
  out
}

#' Pull every `ps()` call out of a nonlinear body and rewrite the body to
#' call the term's own closure instead.
#'
#' The rewrite is what makes several `ps()` calls in one body possible:
#' each becomes a call to a distinctly named closure, so the objective
#' does not have to work out which block a bare `ps` meant.
#'
#' @noRd
ps_extract <- function(body, dpar) {
  terms <- list()
  walk <- function(e, inside) {
    if (!is.call(e)) return(e)
    fn <- e[[1L]]
    if (is.name(fn) && identical(as.character(fn), "ps")) {
      if (inside) {
        stop("ps() inside ps(): a penalized block cannot be the ",
             "argument of another one, because the inner block's knots ",
             "would have to be rebuilt whenever its coefficients move",
             call. = FALSE)
      }
      pt <- ps_spec(e, length(terms) + 1L, dpar)
      terms[[length(terms) + 1L]] <<- pt
      return(as.call(list(as.name(pt[["fname"]]),
                          walk(pt[["expr"]], TRUE))))
    }
    for (i in seq_along(e)) {
      # an empty argument (`m[, 1]`) is the missing-argument symbol and
      # has to be FORCED inside the handler; extracting it succeeds and
      # passing it on is what fails
      a <- tryCatch({
        ai <- e[[i]]
        force(ai)
        ai
      }, error = function(err) NULL)
      if (is.call(a)) e[[i]] <- walk(a, inside)
    }
    e
  }
  list(body = walk(body, FALSE), terms = terms)
}

#' Validate one `ps()` call and turn it into a term specification.
#'
#' `k`, `degree`, `pad` and `center` are basis settings, not model terms,
#' so they are read at parse time from constants. A value that has to be
#' computed from data would make the basis depend on the data twice.
#'
#' @noRd
ps_spec <- function(call, index, dpar) {
  mc <- tryCatch(match.call(ps, call, expand.dots = FALSE),
                 error = function(e) {
                   # `match.call()`'s own "unused argument" names the
                   # argument and nothing else, so a reader is told that
                   # `by = g` is unused without being told by what
                   stop("ps(): unknown argument in ", deparse1(call),
                        ". ps() takes expr, k, degree, pad and center. ",
                        "The penalty is always the second difference; ",
                        "`by =` and `id =`, which a factor-smooth ",
                        "surface would need, are not implemented and ",
                        "are not silently ignored", call. = FALSE)
                 })
  if (is.null(mc[["expr"]])) {
    stop("ps() needs an expression to evaluate the spline at, as its ",
         "first argument", call. = FALSE)
  }
  # A LITERAL, not an expression that happens to evaluate. Evaluating in
  # `baseenv()` looked safe and was not: `ps(x, k = length(t))` resolves
  # `t` to the transpose function there and comes back as 1, so a
  # data-dependent basis size was silently accepted as a small one.
  lit <- function(nm, default) {
    if (is.null(mc[[nm]])) return(default)
    v <- mc[[nm]]
    if (!(is.numeric(v) || is.logical(v)) || length(v) != 1L) {
      stop("ps(", nm, " = ) must be a literal constant, because the ",
           "basis is frozen before any data is read; got ",
           deparse1(mc[[nm]]), call. = FALSE)
    }
    v
  }
  degree <- as.integer(lit("degree", 3L))
  k <- as.integer(lit("k", 10L))
  pad <- as.numeric(lit("pad", 0.1))
  center <- isTRUE(lit("center", TRUE))
  if (is.na(degree) || degree < 1L || degree > 5L) {
    stop("ps(degree = ) must be a whole number between 1 and 5",
         call. = FALSE)
  }
  kmin <- max(degree + 1L, 4L)
  if (is.na(k) || k < kmin) {
    stop("ps(k = ) must be at least ", kmin, " for degree ", degree,
         ", so the penalty has a range space to put in the ",
         "random-effect block", call. = FALSE)
  }
  if (k > 50L) {
    stop("ps(k = ) above 50 is refused: the divided-difference basis ",
         "cancels terms of order (range / spacing)^degree, so its ",
         "agreement with splineDesign() degrades like k^(degree - 1) ",
         "and is already 1.8e-11 at k = 40", call. = FALSE)
  }
  if (is.na(pad) || pad < 0) {
    stop("ps(pad = ) must be a non-negative fraction of the range of ",
         "the spline argument", call. = FALSE)
  }
  list(index = index, dpar = dpar,
       fname = paste0(".frm_ps_", dpar, "_", index),
       expr = mc[["expr"]], k = k, degree = degree, ord = degree + 1L,
       pad = pad, center = center, label = deparse1(call))
}

#' Frame time: the frozen knots, the divided-difference weights and the
#' eigensplit of the second-difference penalty.
#'
#' The split follows Wood (2004) exactly, which is what
#' `mgcv::smooth2random()` does for `s()`: eigendecompose `S`, send the
#' null space to the fixed coefficients and the range space, rescaled by
#' the inverse square roots of the nonzero eigenvalues, to a
#' random-effect block whose single variance IS the inverse smoothing
#' parameter.
#'
#' @noRd
ps_build <- function(pt, mf, nlpars, env) {
  x0 <- ps_data_time(pt, mf, nlpars, env)
  rng <- range(x0)
  span <- rng[2L] - rng[1L]
  if (!(span > 0)) {
    stop("ps(): the spline argument ", deparse1(pt[["expr"]]),
         " is constant over the data, so there is no range to place ",
         "knots on", call. = FALSE)
  }
  lo <- rng[1L] - pt[["pad"]] * span
  hi <- rng[2L] + pt[["pad"]] * span
  m <- pt[["k"]]
  ord <- pt[["ord"]]
  knots <- ps_knots(lo, hi, m, ord)
  W <- ps_weights(knots, ord)

  # second differences of adjacent coefficients, the penalty the paper
  # this term was built for uses
  D2 <- matrix(0, m - 2L, m)
  for (i in seq_len(m - 2L)) D2[i, i + 0:2] <- c(1, -2, 1)
  S <- crossprod(D2)
  # the sum-to-zero constraint, absorbed before the eigensplit so that
  # the reported null space is the constrained one
  Q <- if (pt[["center"]]) {
    qr.Q(qr(matrix(1, m, 1L)), complete = TRUE)[, -1L, drop = FALSE]
  } else {
    diag(m)
  }
  Sc <- crossprod(Q, S %*% Q)
  ev <- eigen(Sc, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * sqrt(.Machine$double.eps)
  U0 <- Q %*% ev$vectors[, !keep, drop = FALSE]
  Us <- Q %*% ev$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(ev$values[keep]), sum(keep), sum(keep))
  pt[["knots"]] <- knots
  pt[["W"]] <- W
  pt[["U0"]] <- U0
  pt[["Us"]] <- Us
  pt[["n_fixed"]] <- ncol(U0)
  pt[["n_pen"]] <- ncol(Us)
  pt[["data_range"]] <- rng
  pt[["knot_range"]] <- c(knots[ord], knots[m + 1L])
  pt
}

#' The spline argument at data time: the expression evaluated on the
#' model frame with every nonlinear parameter at zero.
#'
#' A parameter has no value when the knots are placed, and the knots have
#' to be placed before the first likelihood evaluation. Zero is the
#' honest stand-in - it is where every nonlinear parameter's own linear
#' predictor starts - and `pad` is what covers the distance the fit then
#' moves.
#'
#' @noRd
ps_data_time <- function(pt, mf, nlpars, env) {
  vars <- nl_body_vars(pt[["expr"]])
  zeros <- stats::setNames(as.list(rep(0, length(vars))), vars)
  zeros <- zeros[setdiff(vars, names(mf))]
  x0 <- tryCatch(
    as.numeric(eval(pt[["expr"]], c(as.list(mf), zeros), env)),
    error = function(e) {
      stop("ps(): the spline argument ", deparse1(pt[["expr"]]),
           " could not be evaluated on the model frame: ",
           conditionMessage(e), call. = FALSE)
    })
  if (!length(x0) || anyNA(x0) || any(!is.finite(x0))) {
    stop("ps(): the spline argument ", deparse1(pt[["expr"]]),
         " is not finite everywhere on the model frame, so no knot ",
         "range exists", call. = FALSE)
  }
  x0
}

#' Whether a frame carries any `ps()` term.
#'
#' @noRd
frame_has_ps <- function(frame) {
  for (lp in frame[["linpreds"]]) {
    if (length(lp[["ps_terms"]] %||% list())) return(TRUE)
  }
  FALSE
}

#' Every `ps()` term on a frame, flattened.
#'
#' @noRd
frame_ps_terms <- function(frame) {
  out <- list()
  for (lp in frame[["linpreds"]]) {
    for (pt in lp[["ps_terms"]] %||% list()) out[[length(out) + 1L]] <- pt
  }
  out
}

#' The fitting options a `ps()` block cannot be combined with.
#'
#' Each of REML, quadrature and profiling integrates something out with
#' a Laplace approximation about a single inner mode, and a `ps()`
#' block's null-space coefficients are already inside that integral
#' through a body that is nonlinear in them. Refusing by name beats
#' returning a number nobody can interpret.
#'
#' @noRd
check_ps_fit <- function(frame, REML, quadrature, control) {
  if (!frame_has_ps(frame)) return(invisible(NULL))
  lab <- frame_ps_terms(frame)[[1L]][["label"]]
  if (isTRUE(REML)) {
    stop("REML = TRUE cannot be combined with ", lab,
         ": REML integrates the fixed coefficients out, and this term ",
         "puts its own null space among them inside a body that is ",
         "nonlinear in it", call. = FALSE)
  }
  if (isTRUE(quadrature)) {
    stop("quadrature = TRUE cannot be combined with ", lab,
         ": marginalizing this block by Gauss-Kronrod means integrating ",
         "over a whole curve, and the rule takes one scalar random ",
         "effect at a time", call. = FALSE)
  }
  if (isTRUE(control$profile)) {
    stop("frmtmb_control(profile = TRUE) cannot be combined with ", lab,
         ": profiling assumes the objective is quadratic in the fixed ",
         "coefficients, which a nonlinear body makes it not",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Say so when a prediction leaves the frozen knot span.
#'
#' The basis is a partition of unity only between the outer knots. Past
#' `knot_range` it is a PARTIAL SUM, so the curve decays smoothly, and
#' past the last knot it is exactly zero, so the body reads the curve as
#' zero and the prediction bends to whatever the rest of the body gives.
#' Both are arithmetically correct and neither is an extrapolation of
#' the fitted curve, and a growth curve that bends gently to its
#' intercept is exactly the shape a user will not question.
#'
#' `ps_coverage_warning()` covers the FITTED rows at fit end. This is
#' the same statement at `predict(newdata = )` and at
#' `frm_lp_basis(newdata = )`, the two doors that evaluate the curve
#' somewhere the fit never saw. It is placed inside the closure rather
#' than beside it because the closure is the only place the evaluated
#' argument exists: the expression may name nonlinear parameters, so the
#' value is not known until the body has been walked.
#'
#' `inherits(x, "advector")`, not `is.numeric()`, is what keeps it off
#' the tape. `is.numeric()` is TRUE for an advector, so it guards
#' nothing; it went unnoticed while `predict()`'s nonlinear branch was
#' the only armed caller, because that branch evaluates the body off the
#' tape. `frm_lp_basis()` tapes the body, and a `ps()` term whose
#' argument names a nonlinear parameter reaches this function as an
#' advector there, where `<` raises rather than answers.
#'
#' `sink` is `TRUE` to warn now, or an environment to collect into. See
#' `ps_span_flush()` for why a taped caller collects.
#'
#' @noRd
ps_span_warning <- function(x, pt, sink = TRUE) {
  if (inherits(x, "advector") || !is.numeric(x) || !length(x)) {
    return(invisible(NULL))
  }
  lo <- pt[["knot_range"]][1L]
  hi <- pt[["knot_range"]][2L]
  out <- sum(x < lo | x > hi, na.rm = TRUE)
  if (!out) return(invisible(NULL))
  if (is.environment(sink)) {
    # keyed on the term, so a tape build that re-enters the closure
    # leaves one entry rather than one per pass
    assign(pt[["label"]], list(out, length(x), pt), envir = sink)
    return(invisible(NULL))
  }
  ps_span_signal(out, length(x), pt)
}

#' One warning per collected `ps()` term, raised by the caller that
#' armed the collector rather than by the closure.
#'
#' A caller that tapes the body cannot warn from inside it: the number
#' of R-level passes a tape build makes is RTMB's business, and today it
#' is one, so warning in place would be right by accident. Collecting
#' and flushing makes the count one per term per call by construction.
#'
#' @noRd
ps_span_flush <- function(sink) {
  if (!is.environment(sink)) return(invisible(NULL))
  for (nm in ls(sink, all.names = TRUE, sorted = TRUE)) {
    do.call(ps_span_signal, get(nm, envir = sink))
  }
  invisible(NULL)
}

#' The one place the span message is written. It is assembled here
#' rather than by the caller so that the template stays a literal
#' argument of the condition constructor, which is where
#' `test-message-uniqueness.R` reads templates from.
#'
#' The warning carries a class, so a consumer can catch this one warning
#' without matching on its text. `frmtmb.spline`'s curve functions do
#' exactly that: they surface it once under their own name, and refuse a
#' feature search whose bracket leaves the span.
#'
#' @noRd
ps_span_signal <- function(out, n, pt) {
  warning(warningCondition(paste0(
    out, " of ", n, " predicted values of ",
    deparse1(pt[["expr"]]), " lie outside the frozen knot span of ",
    pt[["label"]], " [", format(pt[["knot_range"]][1L], digits = 4), ", ",
    format(pt[["knot_range"]][2L], digits = 4),
    "]. The basis is a partial sum there ",
    "and exactly zero past the outer knot, so these are not ",
    "extrapolations of the fitted curve: the curve decays and the ",
    "prediction bends to whatever the rest of the body gives. ",
    "Refit with a larger pad = to cover the range you predict on"),
    class = "frmtmb_ps_span_warning"))
  invisible(NULL)
}

#' Fit-end coverage report for `ps()` terms.
#'
#' The basis is exactly zero outside its knot span, so an evaluation
#' point the fit pushes past the padded range reads as a curve value of
#' zero. That is a cliff, not an extrapolation, and it is silent: the
#' likelihood is finite and the gradient is zero there. This says so,
#' with the number to widen.
#'
#' @noRd
ps_coverage_warning <- function(fit) {
  frame <- fit$frame
  if (!frame_has_ps(frame)) return(invisible(NULL))
  bvec <- coef_b(fit)
  vals <- NULL
  for (lp in frame[["linpreds"]]) {
    if (!length(lp[["ps_terms"]] %||% list())) next
    vals <- vals %||% eval_dpars(fit)
    ev <- c(vals[[lp[["resp"]]]][c(lp[["nl_pars"]], lp[["nl_dpar_refs"]])],
            lp[["data_list"]])
    for (pt in lp[["ps_terms"]]) {
      x <- tryCatch(as.numeric(eval(pt[["expr"]], ev, lp[["nl_env"]])),
                    error = function(e) NULL)
      if (is.null(x) || !length(x)) next
      out <- sum(x < pt[["knot_range"]][1L] | x > pt[["knot_range"]][2L])
      if (out > 0L) {
        warning(out, " of ", length(x), " fitted values of ",
                deparse1(pt[["expr"]]), " fall outside the knot span of ",
                pt[["label"]], " [", format(pt[["knot_range"]][1L],
                digits = 4), ", ",
                format(pt[["knot_range"]][2L], digits = 4),
                "], where the basis is exactly zero and its gradient ",
                "with it. Refit with a larger pad =", call. = FALSE)
      }
    }
  }
  invisible(NULL)
}
