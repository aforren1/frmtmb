#' Fit a model
#'
#' Fits a model specified with [bf()] by maximum likelihood, using the
#' Laplace approximation for random effects through RTMB.
#'
#' @param formula A `frmtmb_formula` from [bf()] (with a family attached
#'   via `+`), or a plain formula combined with the `family` argument.
#' @param data A data frame.
#' @param family A family: a `frmtmb_family`, a [stats::family] object or
#'   constructor (for example `gaussian`, `poisson`, `binomial`), or a
#'   family name as a string. Overrides a family already attached to
#'   `formula`.
#' @param REML If `TRUE`, integrate the `mu` fixed effects out of the
#'   likelihood along with the random effects (restricted maximum
#'   likelihood).
#' @param start Optional named list of starting values; components must
#'   match the parameter template (`beta`, `betad`, `theta`).
#' @param control A list from [frmtmb_control()].
#' @param se If `TRUE`, run [RTMB::sdreport()] at fit time. The default
#'   (`FALSE`) defers it until standard errors are first needed
#'   (`summary`, `vcov`, `confint`, `predict(se.fit = TRUE)`), which cuts
#'   roughly a quarter off fit time in fit-and-predict or bootstrap
#'   loops. The deferred report is cached, so nothing is computed twice.
#' @param na.action How to handle missing values, as in [stats::lm()]
#'   (default [stats::na.omit]).
#' @param lower,upper Optional named numeric vectors of hard box
#'   constraints on outer parameters (brms `lb`/`ub`), on the internal
#'   scale, e.g. `lower = c(b = 0)` for a nonlinear rate parameter.
#'   Names as in `confint()` rows.
#' @param priors Optional [set_prior()] specification. This makes the
#'   fit MAP / regularized ML (glmmTMB's `priors=` in spirit): useful
#'   for stabilizing singular variance components or separating
#'   binomials. The reported logLik/AIC then include the prior terms
#'   and are penalized quantities, and `anova()` comparisons across
#'   different priors are meaningless.
#' @param quadrature If `TRUE`, marginalize each scalar random effect by
#'   adaptive Gauss-Kronrod quadrature instead of the Laplace
#'   approximation (the `glmer(nAGQ = k)` analogue; matches it in
#'   tests). Worth it for Bernoulli responses with small clusters,
#'   where Laplace biases variance components. Scalar random-intercept
#'   models only, and not with `mi()`, `trunc()`, `REML = TRUE`, or
#'   `frmtmb_control(profile = TRUE)`. A plain Laplace fit runs first
#'   and the quadrature tape is built at its optimum: the
#'   Gauss-Kronrod rescaling is fixed when the tape is built, so the
#'   starting point decides whether the marginalized objective is
#'   finite at all. That fit also supplies the conditional modes, which
#'   the marginalized objective no longer carries, so `ranef()`,
#'   `fitted()` and `predict()` work as usual.
#' @param dry_run `"spec"` returns the parsed intermediate representation
#'   without touching `data`; `"frame"` returns the assembled design
#'   matrices and parameter template without fitting.
#' @param verbose Report fit progress; a shortcut for
#'   `control = frmtmb_control(verbose =)`, whose value wins when both
#'   are given. See [frmtmb_control()] for the levels and the output.
#' @return An object of class `frmtmb_fit`.
#' @examples
#' \dontrun{
#' data(sleepstudy, package = "lme4")
#' fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
#'               data = sleepstudy)
#' summary(fit)
#' }
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#' dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' summary(fit)
#' fixef(fit)
#' VarCorr(fit)
#'
#' # distributional regression: model sigma too
#' fit2 <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
#' anova(fit, fit2)
#' @export
frm <- function(formula, data, family = NULL, REML = FALSE, start = NULL,
                control = frmtmb_control(), se = FALSE,
                na.action = stats::na.omit, lower = NULL, upper = NULL,
                priors = NULL, quadrature = FALSE, dry_run = NULL,
                verbose = FALSE) {
  cl <- match.call()
  # frmtmb_control() leaves verbose unset (NULL), so an explicit control
  # value always wins over the frm() shortcut
  control$verbose <- control$verbose %||% verbose
  vb <- verbose_level(control)
  bform <- as_bform(formula, family)

  if (vb) t0 <- vb_now()
  spec <- parse_spec(bform)
  if (vb) vb_stage("parse", t0)
  if (identical(dry_run, "spec")) return(spec)

  if (vb) t0 <- vb_now()
  frame <- assemble_frame(spec, data, na.action = na.action,
                          sparse_x = isTRUE(control$sparse_x))
  if (vb) vb_stage("frame", t0, vb_frame_detail(frame))
  check_re_structure(spec, frame, control)
  if (identical(dry_run, "frame")) return(frame)

  fit_assembled(spec, frame, bform, cl, REML = REML, start = start,
                control = control, se = se, lower = lower, upper = upper,
                priors = priors, quadrature = quadrature)
}

# --- verbose progress reporting --------------------------------------
# Stage lines go through message(): suppressMessages() silences them and
# they can never contaminate stdout results. Every call site is guarded
# by an `if (vb)` so a default fit does no clock reads and builds no
# strings. Nothing here reaches inside the tape - per-evaluation
# printing would dominate the fit it is meant to measure.

# Resolved level: 0 = silent, 1 = stage progress, 2 = also the
# optimizer's own iteration trace.
verbose_level <- function(control) {
  v <- control$verbose
  if (is.null(v) || isFALSE(v)) return(0L)
  if (isTRUE(v)) return(1L)
  v <- suppressWarnings(as.integer(v)[1L])
  if (is.na(v) || v < 0L) 0L else v
}

vb_now <- function() proc.time()[["elapsed"]]

vb_say <- function(...) message("frmtmb: ", ...)

# "frmtmb: <stage> [1.23s]: <detail>"
vb_stage <- function(stage, t0, detail = NULL) {
  message("frmtmb: ", stage, " [", sprintf("%.2f", vb_now() - t0), "s]",
          if (is.null(detail)) "" else paste0(": ", detail))
}

vb_plural <- function(n, what) {
  paste0(n, " ", what, if (n != 1L) "s")
}

vb_frame_detail <- function(frame) {
  paste0(frame$n_obs, " obs, ",
         vb_plural(length(frame$linpreds), "linear predictor"), ", ",
         vb_plural(length(frame$re_blocks), "random-effect block"))
}

# First line of a fit: the family and every mode that changes what the
# optimizer is solving, so a slow log says which problem it is timing.
vb_fit_detail <- function(spec, REML, control, quadrature, priors) {
  fams <- vapply(spec$responses, function(r) r$family$family %||% "?", "")
  opt <- control$optimizer %||% "nlminb"
  if (is.function(opt)) opt <- "custom"
  flags <- c(if (isTRUE(REML)) "REML" else "ML",
             if (isTRUE(control$profile)) "profile",
             if (isTRUE(quadrature)) "quadrature",
             if (isTRUE(control$autoscale)) "autoscale",
             if (!is.null(priors)) "priors")
  paste0(paste(unique(fams), collapse = " + "), ", ",
         paste(flags, collapse = ", "), ", ", opt)
}

vb_opt_detail <- function(opt) {
  paste0("objective ", format(opt$objective, digits = 8),
         if (opt$convergence != 0) {
           paste0(", convergence ", opt$convergence)
         })
}

# verbose >= 2 turns on the optimizer's own iteration trace, unless the
# user already asked for one. That trace is printed by nlminb/optim
# themselves and so goes to stdout, not through message(); a custom
# optimizer receives optCtrl untouched.
vb_trace_ctrl <- function(optCtrl, optimizer) {
  optCtrl <- optCtrl %||% list()
  if (!is.null(optCtrl[["trace"]])) return(optCtrl)
  if (identical(optimizer, "nlminb")) {
    optCtrl$trace <- 1L
  } else if (identical(optimizer, "optim")) {
    optCtrl$trace <- 1L
    if (is.null(optCtrl[["REPORT"]])) optCtrl$REPORT <- 1L
  }
  optCtrl
}

# Fitting core shared by frm() and refit(): objective build through the
# convergence check. A non-NULL `template` bypasses make_start (warm
# starts when refitting to a new response).
fit_assembled <- function(spec, frame, bform, cl, REML, start, control,
                          se, lower, upper, priors, quadrature,
                          template = NULL) {
  lower_arg <- lower
  upper_arg <- upper
  vb <- verbose_level(control)
  if (vb) {
    t_fit <- vb_now()
    vb_say("fit: ", vb_fit_detail(spec, REML, control, quadrature,
                                  priors))
  }
  ascale <- if (isTRUE(control$autoscale)) autoscale_plan(frame)
  if (!is.null(ascale) && is.null(template)) {
    # two-stage warm start: fit the standardized frame, back-transform
    # the optimum, and continue below as the ordinary unscaled fit
    # (see R/autoscale.R). Doubles the (cheap) optimization. A caller
    # template (refit and friends) skips the pre-fit but keeps the
    # plan, so the optimizer and sdreport still run in natural units.
    if (vb) t0 <- vb_now()
    template <- autoscale_prefit(spec, frame, bform, cl, REML = REML,
                                 start = start, control = control,
                                 lower = lower, upper = upper,
                                 priors = priors,
                                 quadrature = quadrature, plan = ascale)
    if (vb) vb_stage("autoscale pre-fit", t0)
  }
  if (vb) t0 <- vb_now()
  nll <- build_objective(frame)
  if (!is.null(priors)) {
    # MAP / regularized ML: the optimized objective includes the prior
    # terms, so logLik/AIC are penalized quantities - documented
    ri <- resolve_prior_input(list(frame = frame, spec = spec), priors)
    if (length(ri$entries)) {
      nll0 <- nll
      nlp <- neg_log_prior_fn(ri$entries)
      nll <- function(pars) nll0(pars) + nlp(pars)
    }
    if (length(ri$lower)) {
      lower <- utils::modifyList(as.list(ri$lower),
                                 as.list(lower %||% c()))
      lower <- unlist(lower)
    }
    if (length(ri$upper)) {
      upper <- utils::modifyList(as.list(ri$upper),
                                 as.list(upper %||% c()))
      upper <- unlist(upper)
    }
  }
  if (is.null(template)) template <- make_start(frame, start)

  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- c(if (!is.null(template[["b"]])) "b",
              if (!is.null(template[["miss"]])) "miss")
  if (REML) random <- c(random, "beta")
  if (!length(random)) random <- NULL

  # A mixture likelihood is invariant to permuting its components, so
  # the mu coefficients enter a multimodal objective. Both REML and
  # profile = TRUE integrate those coefficients out with a Laplace
  # approximation about a single inner mode, which is not defined here:
  # the inner Newton solve walks between the component modes and the
  # fit either dies at "NA/NaN gradient evaluation" or reports an
  # optimum with a gradient near 1e9. Quadrature is unaffected because
  # it marginalizes the random effects, not the coefficients.
  if (has_mixture(spec)) {
    if (REML) {
      stop("REML = TRUE cannot be combined with mixture(): the ",
           "mixture likelihood is multimodal in the fixed effects ",
           "REML integrates out, so the restricted likelihood is not ",
           "defined. Use REML = FALSE", call. = FALSE)
    }
    if (isTRUE(control$profile)) {
      stop("frmtmb_control(profile = TRUE) cannot be combined with ",
           "mixture(): profiling moves the fixed effects into the ",
           "inner Laplace problem, and the mixture likelihood is ",
           "multimodal in them. Use profile = FALSE", call. = FALSE)
    }
  }

  integrate <- NULL
  if (isTRUE(quadrature)) {
    if (!is.null(template[["miss"]])) {
      stop("quadrature = TRUE cannot be combined with mi()",
           call. = FALSE)
    }
    scalar_iid <- vapply(frame$re_blocks, function(bk) {
      bk$dim == 1L && bk$covstruct %in% c("us", "diag", "homdiag")
    }, TRUE)
    if (!length(scalar_iid) || !all(scalar_iid)) {
      stop("quadrature = TRUE currently supports scalar random ",
           "intercepts only (every block must be a dim-1 us/diag term)",
           call. = FALSE)
    }
    if (REML) {
      stop("quadrature = TRUE cannot be combined with REML = TRUE",
           call. = FALSE)
    }
    # The truncation normalizer is log(F(ub) - F(lb)) over plain CDFs.
    # The Gauss-Kronrod nodes reach random-effect values where that
    # difference underflows to exactly zero while the density itself is
    # still representable, so the integrand is +Inf there and the
    # marginalized objective comes back -Inf - at the Laplace optimum
    # as well as at the starting values. Laplace never leaves the
    # neighborhood of the mode and is unaffected. Refusing beats
    # reporting logLik = +Inf as a converged fit.
    trunc_resp <- names(which(vapply(
      frame$aterm_values,
      function(a) !is.null(a$trunc_lb) || !is.null(a$trunc_ub), TRUE)))
    if (length(trunc_resp)) {
      stop("quadrature = TRUE cannot be combined with trunc() (",
           paste(trunc_resp, collapse = ", "), "): the truncation ",
           "normalizer underflows at the Gauss-Kronrod nodes and the ",
           "marginalized objective is unbounded. Use quadrature = ",
           "FALSE (Laplace), REML = TRUE, or ",
           "frmtmb_control(profile = TRUE)", call. = FALSE)
    }
    # adaptive Gauss-Kronrod marginalization per scalar random effect
    # (TMB's experimental `integrate`; the nAGQ analogue - matches
    # glmer(nAGQ = 25) in tests). Spec replicated from TMB's GK().
    integrate <- list(b = structure(
      list(dim = 1, adaptive = FALSE, debug = FALSE,
           method = "marginal_gk"),
      class = "GK"
    ))
  }
  profile_arg <- NULL
  if (isTRUE(control$profile)) {
    if (REML) {
      stop("frmtmb_control(profile = TRUE) cannot be combined with ",
           "REML = TRUE (beta is already integrated)", call. = FALSE)
    }
    if (isTRUE(quadrature)) {
      stop("frmtmb_control(profile = TRUE) cannot be combined with ",
           "quadrature = TRUE", call. = FALSE)
    }
    profile_arg <- "beta"
  }
  # Under integrate= the objective is built in two passes: the plain
  # Laplace tape first, then the marginalized one calibrated at its
  # optimum. See quad_fit() for why. The Laplace objective is kept
  # afterwards because it is the only source of the conditional modes.
  lap_obj <- if (is.null(integrate)) NULL else {
    RTMB::MakeADFun(nll, template, random = random, map = frame$map,
                    silent = TRUE)
  }
  obj <- if (is.null(integrate)) {
    RTMB::MakeADFun(nll, template, random = random, map = frame$map,
                    profile = profile_arg, silent = TRUE)
  } else lap_obj
  if (vb) {
    vb_stage("tape", t0,
             paste0(length(obj$par), " outer, ",
                    length(obj$env$random), " inner parameters"))
  }
  # control must ride along: outer_par_names drops beta under
  # profile = TRUE, and a shim without it misaligns every bound
  bounds <- resolve_bounds(list(frame = frame, REML = REML,
                                control = control), lower, upper)
  # a badly scaled coefficient has a badly scaled gradient too: judge
  # (and steer) the optimizer in per-parameter natural units
  par_units <- if (!is.null(ascale)) {
    autoscale_units(frame, ascale, names(obj$par))
  }
  # the optimizer trace rides on the optimizer's own control list, so
  # keep it out of the control stored on the fit (refit and friends
  # reuse that list and must not inherit a trace)
  ctl_opt <- control
  if (vb >= 2L) {
    ctl_opt$optCtrl <- vb_trace_ctrl(control$optCtrl, control$optimizer)
  }
  if (is.null(integrate)) {
    opt <- nl_start_context(spec, start,
                            optimize_obj(obj, ctl_opt, bounds, par_units,
                                         verbose = vb))
  } else {
    qf <- quad_fit(nll, template, random, frame$map, integrate, lap_obj,
                   ctl_opt, bounds, par_units, vb)
    obj <- qf$obj
    opt <- qf$opt
  }

  # Estimates come cheaply from the parameter list at the optimum;
  # sdreport (a quarter of typical fit time) is computed on demand
  # through sdr_of() unless se = TRUE asked for it now.
  est <- obj$env$parList(opt$par)
  for (nm in names(frame$par_template)) {
    names(est[[nm]]) <- names(frame$par_template[[nm]])
  }
  # under integrate= (quadrature), parList leaves outer components NA;
  # the optimizer vector is authoritative for them either way
  pn <- names(opt$par)
  for (cp in setdiff(unique(pn), random %||% character(0))) {
    pos <- seq_along(frame$par_template[[cp]])
    if (cp == "betad" && length(frame$betad_fixed_idx)) {
      pos <- setdiff(pos, frame$betad_fixed_idx)
    }
    est[[cp]][pos] <- unname(opt$par[pn == cp])
  }
  if (!is.null(integrate)) {
    # integrate= removes the random effects from the tape, so this
    # objective has no conditional modes to report: parList() leaves
    # them NA and slides the outer values into their slots. ranef(),
    # fitted() and predict(newdata =) all read them, so recover them
    # from the inner Newton solve of the Laplace objective at the
    # quadrature optimum.
    inner <- solved_par_list(lap_obj, opt$par)
    for (cp in random) est[[cp]][] <- inner[[cp]]
  }
  fit <- structure(
    list(spec = spec, frame = frame, obj = obj, opt = opt, sdr = NULL,
         REML = REML, estimates = est, priors = priors,
         bform = bform, call = cl,
         control = control, quadrature = isTRUE(quadrature),
         lower = lower_arg, upper = upper_arg, par_units = par_units,
         cache = new.env(parent = emptyenv())),
    class = "frmtmb_fit"
  )
  if (se) {
    if (vb) t0 <- vb_now()
    fit$cache$sdr <- autoscale_sdreport(fit)
    if (vb) vb_stage("sdreport", t0)
  }
  chk <- check_convergence(fit, control)
  if (vb) {
    vb_stage("done", t_fit,
             paste0("objective ", format(fit$opt$objective, digits = 8),
                    ", max|grad| ", format(chk$grad, digits = 3), ", ",
                    vb_plural(length(chk$warnings), "warning")))
  }
  fit
}

# Does any response carry a mixture() family?
has_mixture <- function(spec) {
  any(vapply(spec$responses,
             function(r) !is.null(r$family[["mix"]]), TRUE))
}

# parList() at an outer parameter vector with the inner problem solved
# there first. fn() runs the inner Newton solve and leaves the
# conditional modes in last.par, which is parList's default `par`.
solved_par_list <- function(obj, par) {
  obj$fn(par)
  obj$env$parList(par)
}

# Build and optimize the Gauss-Kronrod (integrate=) objective.
#
# TMBad's marginal_gk transform rescales each integrand ONCE: it finds
# the mode and curvature of the log-integrand by finite differences and
# bakes that (mu, sigma) pair into the tape as constants, at whichever
# parameter values `template` happens to hold. Everything downstream
# depends on that one calibration, so this function has to do two
# things the transform does not do for itself.
#
# 1. Tape at a sensible point. From the cold start the frozen rescaling
#    sits far from the real conditional mode, and for every family
#    whose inverse link exponentiates the linear predictor the rescaled
#    integrand then overflows: obj$fn() is NaN before the optimizer
#    takes a step (poisson, Gamma and Beta over nested scalar blocks,
#    Beta over a single one). Gaussian responses survive it only
#    because their integrand is quadratic wherever it is sampled. So
#    fit the plain Laplace objective first and tape the marginalized
#    one at that optimum: the two optima maximize the same marginal
#    likelihood, one exactly and one to O(n^-1).
#
# 2. Recalibrate when the tape expires. A frozen rescaling is only
#    trustworthy near the point it was made at, so it can run out in
#    two ways. The optimizer can walk far enough that the rescaled
#    integrand breaks (RTMB then raises "NA/NaN gradient evaluation"
#    from inside nlminb), or it can stop somewhere the tape's own
#    gradient does not vanish - a mixture whose Laplace fit collapses a
#    mixing weight to exp(-35) does that, and the reported objective is
#    then a value no neighborhood shares. Either way the answer is to
#    tape again at the best point reached and carry on, and to keep the
#    cold template as a last anchor when the Laplace optimum is the bad
#    one. Each candidate costs a tape, so they are tried in order and
#    the first stationary result wins.
quad_fit <- function(nll, template, random, map, integrate, lap_obj,
                     control, bounds, par_units, vb = 0L, rounds = 3L) {
  if (vb) t0 <- vb_now()
  lap_opt <- optimize_obj(lap_obj, control, bounds, par_units)
  if (vb) vb_stage("quadrature warm start", t0, vb_opt_detail(lap_opt))

  # A template holding the outer values `par` plus the conditional
  # modes there, which is what MakeADFun needs to calibrate the tape.
  anchor <- function(par) solved_par_list(lap_obj, par)[names(template)]

  # One build-and-optimize pass. Returns the fit, or - when the tape
  # broke mid-optimization - the best point it reached, so the caller
  # can recalibrate there.
  attempt <- function(tpl) {
    if (vb) t0 <- vb_now()
    o <- RTMB::MakeADFun(nll, tpl, random = random, map = map,
                         integrate = integrate, silent = TRUE)
    f0 <- try(o$fn(o$par), silent = TRUE)
    if (inherits(f0, "try-error") || !is.finite(f0)) return(NULL)
    if (vb) vb_stage("quadrature tape", t0)
    op <- tryCatch(optimize_obj(o, control, bounds, par_units,
                                verbose = vb),
                   error = function(e) NULL)
    if (is.null(op) || !is.finite(op$objective)) {
      pb <- o$env$last.par.best
      if (is.null(pb) || length(pb) != length(o$par) || anyNA(pb)) {
        return(NULL)
      }
      return(list(retry = stats::setNames(as.numeric(pb),
                                          names(o$par))))
    }
    g <- try(max(abs(o$gr(op$par) * (par_units %||% 1))), silent = TRUE)
    if (inherits(g, "try-error")) g <- NA_real_
    list(obj = o, opt = op, stationary = isTRUE(g < control$grad_tol))
  }

  best <- NULL
  tpl <- anchor(lap_opt$par)
  cold_left <- TRUE
  for (i in seq_len(rounds)) {
    a <- attempt(tpl)
    if (!is.null(a) && !is.null(a$obj)) {
      if (is.null(best) || a$stationary) best <- a
      if (a$stationary) break
    }
    nxt <- if (!is.null(a) && !is.null(a$retry)) {
      # the frozen rescaling expired where the optimizer walked to:
      # tape again at the best point it managed
      anchor(a$retry)
    } else if (cold_left) {
      # the anchor itself is the problem (the Laplace optimum can sit
      # on a singular variance component): the untouched template is a
      # different, sometimes better, calibration point
      cold_left <- FALSE
      template
    }
    if (is.null(nxt)) break
    tpl <- nxt
  }
  if (is.null(best)) {
    stop("quadrature = TRUE could not marginalize this model: the ",
         "Gauss-Kronrod objective broke down for this likelihood ",
         "(a non-finite objective or gradient). Refit with ",
         "quadrature = FALSE", call. = FALSE)
  }
  best[c("obj", "opt")]
}

# The joint precision is the covariance source for parameters outside
# cov.fixed: REML (beta random) and control profile = TRUE (beta inner).
needs_jp <- function(fit) {
  fit$REML || isTRUE(fit$control$profile)
}

# Memoized sdreport: the standard-error machinery (summary, vcov,
# confint, predict se.fit, diagnose) triggers it on first use.
sdr_of <- function(fit) {
  cache <- fit$cache
  if (is.null(cache$sdr)) {
    cache$sdr <- autoscale_sdreport(fit)
  }
  cache$sdr
}

#' Control parameters for frmtmb fits
#'
#' @param optimizer `"nlminb"` (default), `"optim"` (L-BFGS-B), or a
#'   function with signature `(par, fn, gr, lower, upper, control)`
#'   returning a list with elements `par`, `objective`, `convergence`
#'   (0 = success), and optionally `message` - the hook for optimx,
#'   nloptr, and friends without frmtmb depending on them. Example:
#'   ```
#'   nlopt <- function(par, fn, gr, lower, upper, control) {
#'     r <- nloptr::nloptr(par, fn, gr, lb = lower, ub = upper,
#'                         opts = list(algorithm = "NLOPT_LD_LBFGS",
#'                                     xtol_rel = 1e-10, maxeval = 2000))
#'     list(par = r$solution, objective = r$objective,
#'          convergence = as.integer(r$status < 0), message = r$message)
#'   }
#'   frm(..., control = frmtmb_control(optimizer = nlopt))
#'   ```
#' @param optCtrl Control list passed to the built-in optimizers
#'   ([stats::nlminb()] / [stats::optim()]).
#' @param restarts Number of times to restart the optimizer from the
#'   current optimum while the gradient remains above `grad_tol`.
#' @param grad_tol Warn (and restart) if the maximum absolute gradient at
#'   the optimum exceeds this value.
#' @param profile Experimental: move the primary (`beta`) coefficients
#'   into the inner (Laplace) problem, TMB's `profile` argument - the
#'   analog of `glmmTMBControl(profile = TRUE)` and `glmer(nAGQ = 0)`.
#'   Speeds up models with many fixed effects, and like those it is an
#'   approximation: estimates differ slightly from the exact fit.
#'   Coefficient covariance comes from the joint precision. Not
#'   compatible with `REML = TRUE` or `quadrature = TRUE`;
#'   profile/uniroot `confint()` and `hypothesis(method = "profile")`
#'   need a non-profiled fit.
#' @param sparse_x Build the parametric fixed-effect design matrices as
#'   sparse [Matrix::sparse.model.matrix()] objects, the analog of
#'   `glmmTMB(sparseX =)`. Worth it when a many-level fixed factor makes
#'   the dense design dominate memory; estimates are identical either
#'   way. `model.matrix()` on the fit then returns a sparse matrix.
#' @param autoscale Standardize badly scaled continuous predictors
#'   internally (the lme4 >= 1.1.37 feature): fit a copy of the model
#'   with each qualifying fixed-effect column centered and scaled
#'   (scaled only, in a linear predictor without an intercept), map
#'   that optimum back to the original parameterization exactly, and
#'   warm-start the ordinary fit there. Reported results are always on
#'   the original scale, so every downstream method works unchanged;
#'   the cost is a second (cheap) optimization. Columns qualify when
#'   they are parametric and take more than two distinct values;
#'   intercepts, factor contrasts, smooth bases, and mo()/mi() columns
#'   are never touched, and the whole step is a silent no-op when
#'   nothing qualifies. Compatible with `profile = TRUE`. Under
#'   `priors` or bounds, the first stage applies them to the scaled
#'   coefficients; the second stage is the fit that is reported.
#' @param check_nlev_1 What to do about a scalar random-effect term
#'   whose grouping factor has a single level: `"warning"` (default),
#'   `"ignore"`, or `"stop"`, following lme4's `lmerControl()` check
#'   vocabulary. Such a term has no variance to estimate - the single
#'   level is absorbed by the intercept - and its standard deviation
#'   collapses to zero. Structured blocks over several terms per level
#'   (`ar1()`, `us()`, the spatial covariance structures) are never
#'   flagged: one grouping level there is a single realization of a
#'   field, which is the normal way to write them.
#' @param check_olre What to do about an observation-level random
#'   effect - one level per row - on a gaussian, student or lognormal
#'   response: `"warning"` (default), `"ignore"`, or `"stop"`. Its
#'   variance is confounded with the residual standard deviation, so
#'   only their sum is identified and the split between them is
#'   arbitrary. The check is skipped when `sigma` is not free to absorb
#'   it - a `se()` response or a constant `sigma` - which is the
#'   random-effects meta-analysis, and for discrete families, where an
#'   observation-level term is the usual overdispersion model.
#' @param verbose Report fit progress through [message()], one terse
#'   line per stage with its elapsed seconds, so a slow fit shows where
#'   the time went. `FALSE` (default) is silent and costs nothing.
#'   `TRUE` (or `1`) reports parsing, frame assembly, the autoscale
#'   pre-fit, tape construction, each optimizer run and restart,
#'   `sdreport()` when `se = TRUE`, and a closing line with the
#'   objective, the maximum absolute gradient, and the number of
#'   convergence warnings. The fit itself opens with a line naming the
#'   family, the mode (ML or REML, plus profile, quadrature, autoscale,
#'   priors), and the optimizer. `2` adds the optimizer's own trace, by
#'   setting `optCtrl$trace` unless you set it yourself. That trace is
#'   printed by [stats::nlminb()] / [stats::optim()] to standard
#'   output, not through `message()`, and a custom optimizer function
#'   receives `optCtrl` unchanged, so `verbose` does not reach it.
#'
#'   Use `suppressMessages()` to silence a verbose fit. The refit loops
#'   in [frm_bootstrap()], [influence()], and [frm_allfit()] force
#'   `verbose` off, so a verbose fit does not make them print hundreds
#'   of lines; [refit()] and [frm_multiple()] report normally.
#' @return A list of control settings.
#' @export
frmtmb_control <- function(optimizer = "nlminb",
                           optCtrl = list(iter.max = 1000, eval.max = 1000),
                           restarts = 1, grad_tol = 1e-3,
                           profile = FALSE, sparse_x = FALSE,
                           autoscale = FALSE,
                           check_nlev_1 = c("warning", "ignore", "stop"),
                           check_olre = c("warning", "ignore", "stop"),
                           verbose = NULL) {
  # verbose stays NULL when unset, which is how frm(verbose =) knows an
  # explicit control value must win over its own shortcut
  list(optimizer = optimizer, optCtrl = optCtrl, restarts = restarts,
       grad_tol = grad_tol, profile = isTRUE(profile),
       sparse_x = isTRUE(sparse_x), autoscale = isTRUE(autoscale),
       check_nlev_1 = match.arg(check_nlev_1),
       check_olre = match.arg(check_olre),
       verbose = verbose)
}

# lme4's lmerControl runs a battery of structural checks before the fit
# and gives each one an ignore/warning/stop setting; these are the two
# that change what a frmtmb fit MEANS rather than how fast it runs.
# Both currently fit silently to an answer the user did not ask for.
# [lme4 lmerControl checks]
re_check_act <- function(what, msg) {
  switch(what %||% "warning",
         ignore = invisible(NULL),
         stop = stop(msg, call. = FALSE),
         warning(msg, call. = FALSE))
}

check_re_structure <- function(spec, frame, control) {
  gaussian_like <- c("gaussian", "student", "lognormal")
  for (bk in frame$re_blocks) {
    # smooth / gp / hsgp blocks carry a synthetic n_levels of 1 and no
    # grouping levels at all; only real grouping factors are checked
    if (is.null(bk$levels)) next
    # A structured block over several terms per level (ar1, us, the
    # spatial covstructs) is a single realization of a field, and one
    # group level is the normal way to write it; only a SCALAR term
    # loses its variance to a single level.
    if (bk$n_levels == 1L && bk$dim == 1L) {
      re_check_act(
        control$check_nlev_1,
        paste0("Grouping factor '", bk$group_name, "' in `",
               bk$term_label, "` has a single level, so its variance is ",
               "not identified and collapses to zero. Drop the term (it ",
               "is absorbed by the intercept), or set ",
               "frmtmb_control(check_nlev_1 = \"ignore\")"))
      next
    }
    lp <- frame$linpreds[[bk$components[[1L]]$lp_key]]
    resp <- spec$responses[[lp$resp]]
    fam <- resp$family
    if (is.null(fam)) next
    # se() supplies the residual sd row by row and a constant dpar pins
    # it outright; either way sigma is no longer free to absorb the
    # observation-level variance, so the two are identified. That is
    # exactly the random-effects meta-analysis, where the
    # observation-level term IS the between-study variance.
    sigma_free <- is.null(resp$aterms$se) &&
      is.null(frame$linpreds[[linpred_key(lp$resp, "sigma")]]$constant)
    if (sigma_free && bk$n_levels == frame$n_obs && bk$dim == 1L &&
        fam$family %in% gaussian_like &&
        bk$dpar %in% (fam$primary_dpars %||% "mu")) {
      re_check_act(
        control$check_olre,
        paste0("`", bk$term_label, "` gives every observation its own ",
               "random effect, and for a ", fam$family, " response that ",
               "variance is confounded with the residual sd: only their ",
               "sum is identified, so the split between them is ",
               "arbitrary. Observation-level random effects are ",
               "meaningful for discrete families (overdispersion), not ",
               "here. Set frmtmb_control(check_olre = \"ignore\") to ",
               "keep it"))
    }
  }
  invisible(NULL)
}

# One optimizer invocation, normalized to nlminb's result shape.
# par_units (autoscale) carries per-parameter magnitudes into nlminb's
# scaling hook; the custom-optimizer contract is unchanged.
run_optimizer <- function(optimizer, par, fn, gr, lower, upper, control,
                          par_units = NULL) {
  if (is.function(optimizer)) {
    res <- optimizer(par, fn, gr, lower, upper, control)
    need <- c("par", "objective", "convergence")
    if (!all(need %in% names(res))) {
      stop("A custom optimizer must return par, objective, and ",
           "convergence", call. = FALSE)
    }
    res$message <- res$message %||% ""
    return(res)
  }
  switch(optimizer,
    nlminb = stats::nlminb(par, fn, gr, control = control,
                           # PORT iterates in scale * par units
                           scale = if (is.null(par_units)) 1 else
                             1 / par_units,
                           lower = lower, upper = upper),
    optim = {
      ctl <- control[names(control) %in%
                       c("maxit", "factr", "pgtol", "trace", "REPORT")]
      if (is.null(ctl$maxit)) ctl$maxit <- 1000
      # L-BFGS-B ignores parscale (see ?optim); par_units still govern
      # the gradient-based convergence checks
      r <- stats::optim(par, fn, gr, method = "L-BFGS-B",
                        lower = lower, upper = upper, control = ctl)
      list(par = r$par, objective = r$value,
           convergence = r$convergence, message = r$message %||% "")
    },
    stop("Unknown optimizer '", optimizer,
         "' (use \"nlminb\", \"optim\", or a function)", call. = FALSE)
  )
}

optimize_obj <- function(obj, control,
                         bounds = list(lower = -Inf, upper = Inf),
                         par_units = NULL, verbose = 0L) {
  optimizer <- control$optimizer %||% "nlminb"
  # A model with no free outer parameters - every dpar pinned by a
  # constant and every design zero-column, e.g. y | trials(n) ~ 0 - is
  # already at its optimum. nlminb's PORT front end rejects the empty
  # start vector with "'d' must be a nonempty numeric (double) vector",
  # which names nothing the user wrote, so evaluate the template
  # instead and report a converged degenerate fit. [glmmTMB#1325, #1317]
  if (!length(obj$par)) {
    if (verbose) vb_say("no free parameters; evaluating the template")
    return(list(par = obj$par, objective = as.numeric(obj$fn(obj$par)),
                convergence = 0L,
                message = "no free parameters (degenerate model)"))
  }
  if (verbose) t0 <- vb_now()
  opt <- run_optimizer(optimizer, obj$par, obj$fn, obj$gr,
                       bounds$lower, bounds$upper, control$optCtrl,
                       par_units)
  if (verbose) vb_stage("optimize", t0, vb_opt_detail(opt))
  for (i in seq_len(control$restarts)) {
    g <- max(abs(obj$gr(opt$par) * (par_units %||% 1)))
    if (is.finite(g) && g < control$grad_tol) break
    if (verbose) t0 <- vb_now()
    opt2 <- run_optimizer(optimizer, opt$par, obj$fn, obj$gr,
                          bounds$lower, bounds$upper, control$optCtrl,
                          par_units)
    if (verbose) {
      vb_stage(paste0("restart ", i), t0,
               paste0(vb_opt_detail(opt2), ", from max|grad| ",
                      format(g, digits = 3)))
    }
    if (opt2$objective <= opt$objective) opt <- opt2
  }
  opt
}

# A nonlinear model's coefficients have no data-driven starting values -
# make_start() can only seed intercepts through a family's init_dpars,
# and a nonlinear mu has no design of its own - so an nl fit begins at
# zero, where most nonlinear forms are flat, singular, or undefined.
# The optimizer then dies deep inside RTMB with a message that names
# neither the model nor the remedy, so name `start=` here instead of
# leaving the user to guess. [brms#734 doctrine]
nl_start_context <- function(spec, start, expr) {
  nl <- any(vapply(spec$responses,
                   function(r) length(r$nlpars) > 0L, TRUE))
  if (!nl || !is.null(start)) return(expr)
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      stop("The nonlinear fit failed from the default zero starting ",
           "values (", conditionMessage(e), "). Nonlinear models need ",
           "starting values in the right region: supply them with the ",
           "`start` argument of frm(), e.g. start = list(beta = c(...)) ",
           "in the order of fixef(); frm(..., dry_run = \"frame\"",
           ")$par_template shows the layout", call. = FALSE)
    }),
    warning = function(w) {
      # nlminb's own "NA/NaN function evaluation" is the optimizer
      # noticing the same undefined objective the error below names;
      # letting both through would report the failure twice
      if (grepl("NA/NaN", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

make_start <- function(frame, start) {
  tpl <- frame$par_template
  for (lp in frame$linpreds) {
    if (!is.null(lp$constant)) next   # mapped; keep link(constant)
    resp <- frame$spec$responses[[lp$resp]]
    init_fn <- resp$family$init_dpars[[lp$dpar]]
    if (is.null(init_fn)) next
    icpt <- match("(Intercept)", colnames(lp$X))
    if (is.na(icpt)) next
    val <- lp$link$linkfun(init_fn(frame$y[[lp$resp]],
                                   frame$aterm_values[[lp$resp]]))
    if (is.finite(val)) tpl[[lp$par]][lp$idx[icpt]] <- val
  }
  if (!is.null(start)) {
    for (nm in names(start)) {
      if (!nm %in% names(tpl)) {
        stop("Unknown start component: '", nm, "' (template has: ",
             paste(names(tpl), collapse = ", "), ")", call. = FALSE)
      }
      if (length(start[[nm]]) != length(tpl[[nm]])) {
        stop("start$", nm, " must have length ", length(tpl[[nm]]),
             call. = FALSE)
      }
      tpl[[nm]][] <- start[[nm]]
    }
  }
  tpl
}

check_convergence <- function(fit, control) {
  # collected first, warned second, so the verbose summary line can
  # report how many diagnostics the fit raised
  msgs <- character(0)
  if (fit$opt$convergence != 0) {
    # a nonlinear model that started at zero is the likeliest cause, and
    # `start` is the only lever, so name it here [brms#734]
    nl_hint <- if (any(vapply(fit$spec$responses,
                              function(r) length(r$nlpars) > 0L, TRUE))) {
      paste0(". This is a nonlinear model; if you did not set `start`, ",
             "the fit began at zero, which is rarely in the right region")
    } else ""
    msgs <- c(msgs, paste0("Optimizer did not report convergence: ",
                           fit$opt$message, nl_hint))
  }
  # under autoscale the gradient is judged in the same natural units
  # the optimizer used (a 1e6-scale column bounds its coefficient's
  # absolute gradient near machine noise times 1e6)
  g <- if (!length(fit$opt$par)) NA_real_ else {
    try(max(abs(fit$obj$gr(fit$opt$par) * (fit$par_units %||% 1))),
        silent = TRUE)
  }
  if (inherits(g, "try-error")) g <- NA_real_
  if (is.finite(g) && g > control$grad_tol) {
    msgs <- c(msgs, paste0("Large maximum absolute gradient at the ",
                           "optimum (", format(g, digits = 3),
                           "); the fit may not have converged. ",
                           "diagnose() names the offending parameter; ",
                           "see the 'Convergence problems' section of ",
                           "vignette('diagnostics') for the remedies"))
  }
  # pdHess is only known once sdreport has run (se = TRUE); the lazy
  # path surfaces it through summary()/diagnose() instead
  sdr <- fit$cache$sdr
  if (!is.null(sdr) && !is.null(sdr$pdHess) && !isTRUE(sdr$pdHess)) {
    msgs <- c(msgs, paste0("Hessian is not positive definite; standard ",
                           "errors are unreliable. The model may be ",
                           "overparameterized"))
  }
  for (m in msgs) warning(m, call. = FALSE)
  invisible(list(grad = g, warnings = msgs))
}
