#' Fit a model
#'
#' Fits a model specified with [bf()] by maximum likelihood, using the
#' Laplace approximation for random effects through RTMB.
#'
#' @param formula A `frmtmb_formula` from [bf()] (with a family attached
#'   via `+`), or a plain formula combined with the `family` argument.
#' @param data A data frame.
#' @param data2 A named list of objects that are not columns of `data`:
#'   the adjacency matrix of `car()`, the mesh triple of `spde()`, and
#'   the matrices of `gr(prec = )`, `gr(cov = )` and `equalto()`. This
#'   is brms's `data2` argument, with one deliberate extension: brms
#'   accepts a bare name only, while frmtmb also evaluates compound
#'   expressions with `data2` in front of the data mask, so
#'   `gr(g, cov = solve(Q))` finds `Q` there.
#'   Each structural expression resolves from `data2` first, then
#'   `data`, then the formula environment; that last step is what a
#'   model written before `data2` relied on, so old code keeps working.
#'   Prefer `data2`: its objects are stored on the fit, so `saveRDS()`
#'   and a later `refit()`, `influence()` or `update()` in a fresh
#'   session do not need the calling environment to still exist.
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
#'   (default [stats::na.omit]). Rows dropped for missingness are
#'   reported in a message; wrap the call in `suppressMessages()` to
#'   silence it.
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
#'
#' @srrstats {RE1.0} Models are specified through a formula interface:
#'   [bf()] builds a `frmtmb_formula` from one or more R formulas and a
#'   family attaches with `+`. A plain formula plus `family =` is also
#'   accepted. There is no matrix-only entry point.
#' @srrstats {G2.4,G2.4a,G2.4b,G2.4c,G2.4e} Type conversion during frame
#'   assembly is explicit, never implicit. The response is converted with
#'   `as.numeric()`, or `storage.mode(y) <- "double"` for a matrix
#'   response; an ordinal factor response is converted from factor with
#'   `as.numeric()` and a two-level binomial factor with
#'   `as.numeric(y) - 1`; `mo()` category codes use `as.integer()`;
#'   grouping factors use `as.integer()` for the level index and
#'   `as.character()` for the level labels; addition terms other than
#'   `cens()` use `as.numeric()`.
#' @srrstats {G2.5} Where a factor input is expected, the expected kind is
#'   checked and documented. `mo()` requires an ordered factor and errors
#'   otherwise ("mo(): factor variables must be ordered factors"). An
#'   ordinal family takes level order as category order, and warns when
#'   the response is an unordered factor, naming the order it is about to
#'   use: that order is alphabetical unless the user set it, so the model
#'   can differ from the one intended. The fit still runs, which is what
#'   brms does. A factor response for a non-ordinal, non-binomial family
#'   is refused. The compatibility registry states the requirement in
#'   prose and it is rendered in `vignette("compatibility")`.
#' @srrstats {G2.6} One-dimensional responses are pre-processed to a plain
#'   numeric vector regardless of the class they arrive in: a factor, a
#'   one-column matrix, a `scale()`d matrix carrying attributes, or a bare
#'   vector all reach the objective as `as.numeric(as.vector(y))`.
#' @srrstats {G2.8} Pre-processing funnels every input into two internal
#'   classes before any analytic code runs: `parse_spec()` produces a
#'   data-free `frmtmb_spec`, and `assemble_frame()` produces a
#'   `frmtmb_frame` holding the design matrices and the parameter
#'   template. Every sub-function downstream of assembly sees only those
#'   two classes. Both are reachable for inspection through `dry_run`.
#' @srrstats {G2.13} Missing data is checked during pre-processing, before
#'   anything is passed to the optimizer. Rows are removed by `na.action`
#'   inside `stats::model.frame()`; assembly then errors if any missing
#'   value remains in a model variable, and errors if no complete
#'   observation is left.
#' @srrstats {G2.14,G2.14a} `na.action` lets the user choose how missing
#'   data is handled, following [stats::lm()]. `stats::na.fail` errors on
#'   missing data, `stats::na.omit` (the default) and `stats::na.exclude`
#'   drop the affected rows, and `na.exclude` pads `fitted()`,
#'   `residuals()`, `predict()`, and `simulate()` back to the input
#'   length with `NA` in the original positions.
#' @srrstats {G2.14b} Rows dropped for missingness are reported, not
#'   dropped silently: frame assembly emits one `message()` per fit
#'   giving the number of rows removed. It is a message, not a warning,
#'   so `suppressMessages()` silences it for callers that ask for
#'   `na.omit` deliberately, and `na.action()` on the fit still names the
#'   rows.
#' @srrstats {G2.15} No function assumes non-missingness by inheriting a
#'   default. The invariant is established once at the boundary: assembly
#'   errors unless every model variable is free of missing values, so
#'   downstream arithmetic operates on complete data by construction
#'   rather than by defensive `na.rm` flags that would silently change
#'   the estimand.
#' @srrstats {G2.16} Undefined values are handled separately from missing
#'   ones. The response check is explicitly written as
#'   `any(!is.finite(y) & !is.na(y))`, so `Inf`, `-Inf`, and `NaN` are
#'   rejected with their own message while `NA` is left to `na.action`.
#' @srrstats {RE2.1} The processing of missing values is controlled by an
#'   explicit parameter (`na.action`), and `NA`/`NaN` are distinguished
#'   from `Inf` as described under G2.16.
#' @srrstats {RE2.4,RE2.4a} Perfect collinearity among predictors is
#'   detected during pre-processing. Each parametric fixed-effect design
#'   is factorized with `qr()`; when the rank is short of the column
#'   count the aliased columns are named in a `message()` and dropped,
#'   and the null space of the design is stored on the fit so that
#'   `predict()` can refuse rows of `newdata` that are not estimable. A
#'   cheap sparse singular-value screen gates the dense check so that the
#'   sparse and dense backends drop the same columns.
#' @srrstats {RE3.0} Models that fail to converge raise warnings, one per
#'   diagnostic: a nonzero optimizer status (with the optimizer's own
#'   message and, for a nonlinear model, a hint that `start` was not
#'   set), a maximum absolute gradient above `grad_tol`, a Hessian that
#'   is not positive definite, and non-finite standard errors.
#' @srrstats {RE3.1} Those diagnostics are `warning()` conditions, so
#'   `suppressWarnings()` silences them, and the returned object still
#'   carries enough to identify the failure: `fit$opt$convergence` and
#'   `fit$opt$message` hold the optimizer verdict, and [diagnose()]
#'   recomputes the gradient, the Hessian verdict, separation, singular
#'   variance components, and predictor scaling on demand.
#' @srrstats {RE4.0} The return value is a model object of class
#'   `frmtmb_fit`, with the standard `stats` accessor methods defined
#'   for it.
#' @srrstats {RE4.1} `dry_run` generates the model object without fitting
#'   it: `"spec"` returns the parsed representation without touching
#'   `data`, and `"frame"` returns the assembled design matrices and the
#'   parameter template. Both are useful for batch setup and for
#'   inspecting what a formula compiled to.
#' @srrstats {RE4.4} The model specification is recoverable as a formula
#'   through `formula()`, and the original call through `fit$call`.
#' @srrstats {RE4.5} The number of observations used is returned by
#'   `nobs()`, and the rows dropped by `na.action` by `na.action()`.
#' @srrstats {RE4.8} The response and its metadata are retained on the
#'   fit: `fit$frame$y` holds the response per response name,
#'   `fit$frame$y_levels` the original factor levels for ordinal and
#'   categorical responses, and `fit$frame$aterm_values` the addition
#'   terms (`weights`, `trials`, `cens`, `trunc`, `se`, `offset`).
#'   `model.frame()` returns the stored model frame with its row names.
#' @srrstats {RE4.17} `print()` on a `frmtmb_fit` summarizes the model
#'   input (family, formula, grouping structure) and the fitted
#'   coefficients.
#' @srrstats {RE4.18} A distinct `summary()` method is implemented, and it
#'   is computationally non-trivial: with the default `se = FALSE` it
#'   triggers the deferred `RTMB::sdreport()` that produces the standard
#'   errors, then forms z and p values, variance components, group
#'   counts, information criteria, residual correlations, and smooth
#'   effective degrees of freedom. The report is cached on the fit, so a
#'   second `summary()` is cheap.
#'
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
                priors = NULL, quadrature = FALSE, data2 = list(),
                dry_run = NULL, verbose = FALSE) {
  cl <- match.call()
  data2 <- validate_data2(data2)
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
                          sparse_x = isTRUE(control$sparse_x),
                          data2 = data2)
  if (vb) vb_stage("frame", t0, vb_frame_detail(frame))
  check_re_structure(spec, frame, control)
  if (identical(dry_run, "frame")) return(frame)

  fit_assembled(spec, frame, bform, cl, REML = REML, start = start,
                control = control, se = se, lower = lower, upper = upper,
                priors = priors, quadrature = quadrature, data2 = data2)
}

# --- verbose progress reporting --------------------------------------
# Stage lines go through message(): suppressMessages() silences them and
# they can never contaminate stdout results. Every call site is guarded
# by an `if (vb)` so a default fit does no clock reads and builds no
# strings. Nothing here reaches inside the tape - per-evaluation
# printing would dominate the fit it is meant to measure.

#' Resolved level: 0 = silent, 1 = stage progress, 2 = also the
#' optimizer's own iteration trace.
#'
#' @noRd
verbose_level <- function(control) {
  v <- control$verbose
  if (is.null(v) || isFALSE(v)) return(0L)
  if (isTRUE(v)) return(1L)
  v <- suppressWarnings(as.integer(v)[1L])
  if (is.na(v) || v < 0L) 0L else v
}

#' The clock every verbose timing reads: elapsed seconds as one number.
#'
#' @noRd
vb_now <- function() proc.time()[["elapsed"]]

#' One progress line with the package prefix, through `message()`.
#'
#' @noRd
vb_say <- function(...) message("frmtmb: ", ...)

#' One timed stage line, shaped `"frmtmb: <stage> [1.23s]: <detail>"`.
#'
#' @noRd
vb_stage <- function(stage, t0, detail = NULL) {
  message("frmtmb: ", stage, " [", sprintf("%.2f", vb_now() - t0), "s]",
          if (is.null(detail)) "" else paste0(": ", detail))
}

#' A count and its noun, with the plural `s` only when the count needs
#' one.
#'
#' @noRd
vb_plural <- function(n, what) {
  paste0(n, " ", what, if (n != 1L) "s")
}

#' The size of an assembled frame in one clause: observations, linear
#' predictors, and random-effect blocks. The detail of the "frame" stage
#' line.
#'
#' @noRd
vb_frame_detail <- function(frame) {
  paste0(frame$n_obs, " obs, ",
         vb_plural(length(frame$linpreds), "linear predictor"), ", ",
         vb_plural(length(frame$re_blocks), "random-effect block"))
}

#' First line of a fit: the family and every mode that changes what the
#' optimizer is solving, so a slow log says which problem it is timing.
#'
#' @noRd
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

#' An optimizer result in one clause: the objective, plus the status code
#' when the optimizer did not report success.
#'
#' @noRd
vb_opt_detail <- function(opt) {
  paste0("objective ", format(opt$objective, digits = 8),
         if (opt$convergence != 0) {
           paste0(", convergence ", opt$convergence)
         })
}

#' verbose >= 2 turns on the optimizer's own iteration trace, unless the
#' user already asked for one. That trace is printed by nlminb/optim
#' themselves and so goes to stdout, not through message(); a custom
#' optimizer receives optCtrl untouched.
#'
#' @noRd
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

#' Fitting core shared by frm() and refit(): objective build through the
#' convergence check. A non-NULL `template` bypasses make_start (warm
#' starts when refitting to a new response).
#'
#' @noRd
fit_assembled <- function(spec, frame, bform, cl, REML, start, control,
                          se, lower, upper, priors, quadrature,
                          template = NULL, data2 = list()) {
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
    if (length(frame$autocor %||% list())) {
      # the Gauss-Kronrod rule integrates one scalar random effect
      # against a PRODUCT of per-row densities; an R-side residual is a
      # joint density over each group, so no per-row integrand exists
      stop("quadrature = TRUE cannot be combined with the residual ",
           "correlation term ", frame$autocor[[1L]]$label,
           ": the rule integrates a random effect against ",
           "per-observation densities, and this residual is a joint ",
           "density over each group. Use quadrature = FALSE (Laplace) ",
           "or REML = TRUE", call. = FALSE)
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
    opt <- fit_error_context(
      spec, start, REML, control, quadrature, priors,
      tryCatch(optimize_obj(obj, ctl_opt, bounds, par_units, verbose = vb),
               error = function(e) {
                 rs <- fit_recovery_starts(obj, nll, template, random,
                                           frame$map, frame, start,
                                           ctl_opt, bounds, par_units)
                 for (lbl in names(rs)) {
                   if (vb) {
                     vb_say("optimizer failed (", conditionMessage(e),
                            "); restarting from ", lbl)
                   }
                   op <- tryCatch(optimize_obj(obj, ctl_opt, bounds,
                                               par_units, verbose = vb,
                                               start_par = rs[[lbl]]),
                                  error = function(e2) NULL)
                   if (!is.null(op)) return(op)
                 }
                 stop(e)
               }))
  } else {
    qf <- quad_fit(nll, template, random, frame$map, integrate, lap_obj,
                   ctl_opt, bounds, par_units, frame, vb)
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
         bform = bform, call = cl, data2 = data2,
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

#' Does any response carry a mixture() family?
#'
#' @noRd
has_mixture <- function(spec) {
  any(vapply(spec$responses,
             function(r) !is.null(r$family[["mix"]]), TRUE))
}

#' parList() at an outer parameter vector with the inner problem solved
#' there first. fn() runs the inner Newton solve and leaves the
#' conditional modes in last.par, which is parList's default `par`.
#'
#' @noRd
solved_par_list <- function(obj, par) {
  obj$fn(par)
  obj$env$parList(par)
}

#' Which of two quadrature candidates quad_fit() keeps. A stationary one
#' wins outright (the loop stops there); among non-stationary ones the
#' lowest objective is the best point reached, not the first that
#' happened to tape and optimize without breaking.
#'
#' @noRd
quad_keep_best <- function(best, a) {
  if (is.null(a) || is.null(a$obj)) return(best)
  if (isTRUE(a$stationary)) return(a)
  if (is.null(best) || a$opt$objective < best$opt$objective) a else best
}

#' Build and optimize the Gauss-Kronrod (integrate=) objective.
#'
#' TMBad's marginal_gk transform rescales each integrand ONCE: it finds
#' the mode and curvature of the log-integrand by finite differences and
#' bakes that (mu, sigma) pair into the tape as constants, at whichever
#' parameter values `template` happens to hold. Everything downstream
#' depends on that one calibration, so this function has to do two
#' things the transform does not do for itself.
#'
#' 1. Tape at a sensible point. From the cold start the frozen rescaling
#'    sits far from the real conditional mode, and for every family
#'    whose inverse link exponentiates the linear predictor the rescaled
#'    integrand then overflows: obj$fn() is NaN before the optimizer
#'    takes a step (poisson, Gamma and Beta over nested scalar blocks,
#'    Beta over a single one). Gaussian responses survive it only
#'    because their integrand is quadratic wherever it is sampled. So
#'    fit the plain Laplace objective first and tape the marginalized
#'    one at that optimum: the two optima maximize the same marginal
#'    likelihood, one exactly and one to O(n^-1).
#'
#' 2. Recalibrate when the tape expires. A frozen rescaling is only
#'    trustworthy near the point it was made at, so it can run out in
#'    two ways. The optimizer can walk far enough that the rescaled
#'    integrand breaks (RTMB then raises "NA/NaN gradient evaluation"
#'    from inside nlminb), or it can stop somewhere the tape's own
#'    gradient does not vanish - a mixture whose Laplace fit collapses a
#'    mixing weight to exp(-35) does that, and the reported objective is
#'    then a value no neighborhood shares. Either way the answer is to
#'    tape again at the best point reached and carry on, and to keep the
#'    cold template as a last anchor when the Laplace optimum is the bad
#'    one. Each candidate costs a tape, so they are tried in order and
#'    the first stationary result wins; if none is stationary the
#'    candidate with the lowest objective does.
#'
#' 3. Displace the widths when neither anchor holds. What the transform
#'    freezes is a width per random effect, so a calibration that expires
#'    within a step or two is one whose widths are wrong, and the
#'    parameter that sets them is theta. Half a log-SD either way, then a
#'    whole one, is enough where displacement helps at all - it recovers
#'    single-block fits whose every un-displaced calibration dies inside
#'    the optimizer. It cannot recover a nested block, where the outer
#'    integrand is the inner rescaling's output and the objective is NaN
#'    before the optimizer takes a step.
#'
#' @noRd
quad_fit <- function(nll, template, random, map, integrate, lap_obj,
                     control, bounds, par_units, frame = NULL, vb = 0L,
                     rounds = 3L) {
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

  # Calibration points, in the order they are tried. The Laplace optimum
  # first; then the untouched template, because the anchor itself can be
  # the problem (the Laplace optimum can sit on a singular variance
  # component); then the anchor with the integrand widths displaced.
  lap_tpl <- anchor(lap_opt$par)
  seeds <- list(lap_tpl, template)
  if (length(lap_tpl$theta)) {
    seeds <- c(seeds, lapply(c(-0.5, 0.5, -1), function(s) {
      tpl <- lap_tpl
      tpl$theta <- tpl$theta + s
      tpl
    }))
  }

  best <- NULL
  for (seed in seeds) {
    tpl <- seed
    for (i in seq_len(rounds)) {
      a <- attempt(tpl)
      best <- quad_keep_best(best, a)
      if (is.null(a) || isTRUE(a$stationary) || is.null(a$retry)) break
      # the frozen rescaling expired where the optimizer walked to:
      # tape again at the best point it managed
      tpl <- anchor(a$retry)
    }
    if (!is.null(best) && isTRUE(best$stationary)) break
  }
  if (is.null(best)) {
    stop(quad_breakdown_message(frame, lap_tpl$theta), call. = FALSE)
  }
  best[c("obj", "opt")]
}

#' The refusal message when every calibration point breaks down. A bare
#' "the objective was non-finite" tells a user nothing they can act on,
#' and the two shapes this failure takes want different answers, so name
#' whichever applies: an iterated integral over nested blocks (the outer
#' integrand is the inner rescaling's output, and no calibration of the
#' outer one can repair that), and a variance component so narrow that
#' the transform's unit-step finite differences cannot measure it.
#'
#' @noRd
quad_breakdown_message <- function(frame, theta) {
  blocks <- frame$re_blocks %||% list()
  desc <- function(b) paste0("'", b$term_label, "' (", b$n_levels,
                             " levels)")
  sds <- vapply(blocks, function(b) {
    v <- tryCatch(covstruct_registry[[b$covstruct]]$vcov(
      theta[b$theta_idx], b), error = function(e) NULL)
    if (is.null(v)) NA_real_ else sqrt(v[1L, 1L])
  }, 1)
  why <- character(0)
  if (length(blocks) > 1L) {
    why <- c(why, paste0("the model asks for an iterated integral over ",
                         length(blocks), " nested blocks (",
                         paste(vapply(blocks, desc, ""), collapse = ", "),
                         ") on ", frame$n_obs, " observations, and the ",
                         "outer integrand is itself the frozen ",
                         "rescaling of the inner one"))
  }
  sing <- which(is.finite(sds) & sds < 1e-3)
  if (length(sing)) {
    why <- c(why, paste0("the Laplace optimum leaves ",
                         paste(vapply(blocks[sing], desc, ""),
                               collapse = ", "),
                         " on a singular variance component (sd ",
                         paste(format(sds[sing], digits = 2),
                               collapse = ", "),
                         "), which the transform measures by finite ",
                         "differences with a step of 1"))
  }
  if (!length(why)) {
    why <- paste0("the calibration expired on ",
                  paste(vapply(blocks, desc, ""), collapse = ", "),
                  " at every point tried")
  }
  paste0("quadrature = TRUE could not marginalize this model: the ",
         "Gauss-Kronrod objective broke down (a non-finite objective ",
         "or gradient) at every calibration point - the Laplace ",
         "optimum, the best point the optimizer reached, the cold ",
         "start, and the optimum with the integrand widths displaced. ",
         "Here ", paste(why, collapse = "; "),
         ". Refit with quadrature = FALSE: the Laplace approximation ",
         "fits this model (it is what the warm start above already did)")
}

#' The joint precision is the covariance source for parameters outside
#' cov.fixed: REML (beta random) and control profile = TRUE (beta inner).
#'
#' @noRd
needs_jp <- function(fit) {
  fit$REML || isTRUE(fit$control$profile)
}

#' Memoized sdreport: the standard-error machinery (summary, vcov,
#' confint, predict se.fit, diagnose) triggers it on first use.
#'
#' @noRd
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
#'   observation-level term is the usual overdispersion model. It is
#'   also skipped for the known-structure blocks `gr(cov = )`,
#'   `gr(prec = )` and `equalto()`: their level covariance is not
#'   proportional to the identity, so the block and the residual are
#'   separately identified. That is the animal model with one
#'   measurement per individual.
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
#'
#' @srrstats {RE2.0} Data transformations are documented and can be turned
#'   off. `autoscale` is the only transformation of predictor values, it
#'   is `FALSE` by default, and the documentation states exactly which
#'   columns qualify, which are never touched (intercepts, factor
#'   contrasts, smooth bases, `mo()`/`mi()` columns), and that results are
#'   always mapped back and reported on the original scale. `sparse_x`
#'   changes only the storage of the design and is documented as leaving
#'   estimates identical. Factor handling is delegated to
#'   `stats::model.matrix()` with the contrasts frozen at fit time and
#'   reapplied to `newdata`.
#' @srrstats {RE2.3} `autoscale = TRUE` centers and scales qualifying
#'   continuous fixed-effect columns, that is, converts them to z-scores.
#'   Centering is applied only where an intercept exists to absorb the
#'   shift, and the effect is exercised in `tests/testthat/test-autoscale.R`,
#'   which checks that the scaled and unscaled fits agree.
#' @srrstats {RE3.2} Convergence thresholds have documented defaults:
#'   `grad_tol = 1e-3` on the maximum absolute gradient at the optimum,
#'   and `optCtrl = list(iter.max = 1000, eval.max = 1000)` for the
#'   built-in optimizers. Both appear in the usage section of the manual
#'   page with an `@param` describing them.
#' @srrstats {RE3.3} Those thresholds can be set explicitly:
#'   `grad_tol` for the gradient criterion, `optCtrl` for the optimizer's
#'   own tolerances and iteration caps, `restarts` for how many times to
#'   restart from the current optimum while the gradient is still above
#'   `grad_tol`, and `optimizer` to substitute another optimizer
#'   altogether.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 20)))
#' dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'
#' # another optimizer, with its own control list
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'            control = frmtmb_control(optimizer = "optim",
#'                                     optCtrl = list(maxit = 500)))
#' fit$opt$convergence
#'
#' # a tighter gradient criterion, with restarts from the current
#' # optimum until it is met
#' frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
#'     control = frmtmb_control(grad_tol = 1e-4, restarts = 3))
#'
#' # badly scaled predictors: fit an internally standardized copy first,
#' # then warm-start the reported fit from it
#' dd$xbig <- dd$x * 1e5
#' frm(bf(y ~ xbig + (1 | g)) + gaussian(), data = dd,
#'     control = frmtmb_control(autoscale = TRUE))
#'
#' # the object is a plain list, so it can be built once and reused
#' ctrl <- frmtmb_control(check_nlev_1 = "ignore")
#' ctrl$optimizer
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

#' lme4's lmerControl runs a battery of structural checks before the fit
#' and gives each one an ignore/warning/stop setting; these are the two
#' that change what a frmtmb fit MEANS rather than how fast it runs.
#' Both currently fit silently to an answer the user did not ask for.
#' `[lme4 lmerControl checks]`
#'
#' @noRd
re_check_act <- function(what, msg) {
  switch(what %||% "warning",
         ignore = invisible(NULL),
         stop = stop(msg, call. = FALSE),
         warning(msg, call. = FALSE))
}

#' Runs those two checks over the assembled random-effect blocks, before
#' any fitting starts: a scalar term whose grouping factor has one level,
#' and an observation-level term on a response whose residual sd already
#' holds that variance. Each check reports through `re_check_act()`, so
#' the control setting decides between silence, a warning, and an error.
#'
#' @noRd
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
    # A known-structure block is not an OLRE even with one row per
    # level. Its levels are correlated through the fixed relationship
    # matrix (gr(cov = A), gr(prec = Q)) or its covariance is fixed
    # outright (equalto), so the block covariance is no longer
    # proportional to the identity and the residual sd is not a
    # reparameterization of it. That is the animal model: the additive
    # genetic variance and the residual variance are separately
    # identified precisely BECAUSE A is not the identity.
    structured <- bk$covstruct %in% c("gr_cov", "gr_prec", "equalto")
    if (sigma_free && !structured &&
        bk$n_levels == frame$n_obs && bk$dim == 1L &&
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

#' One optimizer invocation, normalized to nlminb's result shape.
#' par_units (autoscale) carries per-parameter magnitudes into nlminb's
#' scaling hook; the custom-optimizer contract is unchanged.
#'
#' @noRd
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

#' A parameter list flattened into the outer vector the optimizer
#' iterates on: the components of names(obj$par), in that order, with
#' the dpars that se() maps to constants removed. Same loop
#' autoscale_units() walks, and the inverse of the one that fills
#' `estimates` from opt$par.
#'
#' @noRd
outer_from_template <- function(tpl, obj, frame) {
  pn <- names(obj$par)
  out <- numeric(0)
  for (cp in unique(pn)) {
    v <- tpl[[cp]]
    if (is.null(v)) return(NULL)
    if (cp == "betad" && length(frame$betad_fixed_idx)) {
      v <- v[-frame$betad_fixed_idx]
    }
    out <- c(out, unname(v))
  }
  if (length(out) != length(pn)) return(NULL)
  stats::setNames(out, pn)
}

#' Starting points to try when the optimizer could not get through at
#' all. Both only ever run after a failure, so a healthy fit pays
#' nothing.
#'
#' The start itself can be unusable. The autoscale pre-fit hands back a
#' warm start, and when the standardized fit ran a correlated block to a
#' perfect correlation the mapped-back template has an infinite
#' objective and a NaN gradient - nlminb dies before its first step,
#' while the cold start make_start() would have built fits the same
#' model. So offer the cold start whenever the fit did not begin there.
#'
#' Or the path from a usable start crosses a hole. Under profile = TRUE
#' beta moves into the inner problem, and the outer objective left
#' behind can be undefined where the optimizer must pass (a Gamma shape
#' intercept on its way to exp(34) does it). The plain Laplace objective
#' over the same model - the one every other mode optimizes - has no
#' such barrier, so its optimum is a starting point on the far side of
#' the hole. Same recipe quad_fit() uses to calibrate the Gauss-Kronrod
#' tape at the Laplace optimum.
#'
#' @noRd
fit_recovery_starts <- function(obj, nll, template, random, map, frame,
                                start, control, bounds, par_units) {
  out <- list()
  differs <- function(p) {
    !is.null(p) && length(p) == length(obj$par) && all(is.finite(p)) &&
      !isTRUE(all.equal(unname(as.numeric(p)), unname(as.numeric(obj$par))))
  }
  cold <- outer_from_template(make_start(frame, start), obj, frame)
  if (differs(cold)) out[["the cold starting values"]] <- cold
  if (isTRUE(control$profile)) {
    plain <- tryCatch({
      o <- RTMB::MakeADFun(nll, template, random = random, map = map,
                           silent = TRUE)
      # bounds and units are resolved for the profiled parameterization,
      # which has no beta in it; this fit only supplies a starting point
      # and the real optimization below applies both
      op <- optimize_obj(o, control)
      outer_from_template(o$env$parList(op$par), obj, frame)
    }, error = function(e) NULL)
    if (differs(plain)) out[["the Laplace optimum"]] <- plain
  }
  # a start outside the model's bounds is rejected by the optimizer
  lapply(out, function(p) pmin(pmax(p, bounds$lower), bounds$upper))
}

#' Second chance after the optimizer aborts the whole fit.
#'
#' nlminb hands RTMB a step that lands outside the region where the
#' likelihood is defined, RTMB raises "NA/NaN gradient evaluation" from
#' inside the optimizer, and the fit dies - on models that fit perfectly
#' well from a different starting point (an ar1 block under
#' profile = TRUE, a wide smooth under autoscale). The tape is not
#' broken: it still holds the best point the line search reached, and
#' restarting there steps around the hole. quad_fit() already does
#' exactly this when a frozen Gauss-Kronrod rescaling expires the same
#' way. One retry only, and only from a point the optimizer preferred to
#' where it started: a second failure is an objective that is genuinely
#' undefined nearby, not an unlucky step.
#'
#' @noRd
optimizer_from_best <- function(obj, par, e, optimizer, bounds, control,
                                par_units, verbose = 0L) {
  pb <- obj$env$last.par.best
  # last.par.best spans the joint vector when the model has random
  # effects (profiled parameters included); lfixed() selects the outer
  # block the optimizer actually iterates on
  p0 <- if (is.null(pb)) NULL else pb[obj$env$lfixed()]
  if (is.null(p0) || length(p0) != length(par) || anyNA(p0) ||
      isTRUE(all.equal(unname(as.numeric(p0)), unname(as.numeric(par))))) {
    stop(e)
  }
  if (verbose) {
    vb_say("optimizer failed (", conditionMessage(e),
           "); restarting from the best point it reached")
  }
  run_optimizer(optimizer, stats::setNames(as.numeric(p0), names(par)),
                obj$fn, obj$gr, bounds$lower, bounds$upper,
                control$optCtrl, par_units)
}

#' The single entry point to the optimizer for every fit mode: it runs
#' the objective to an optimum, restarts from there while the gradient
#' stays above `grad_tol`, and keeps the better result. Bounds, autoscale
#' units, and the restart-from-best recovery are applied here, so no
#' caller has to repeat them.
#'
#' @noRd
optimize_obj <- function(obj, control,
                         bounds = list(lower = -Inf, upper = Inf),
                         par_units = NULL, verbose = 0L,
                         start_par = obj$par) {
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
  run <- function(par) {
    tryCatch(run_optimizer(optimizer, par, obj$fn, obj$gr,
                           bounds$lower, bounds$upper, control$optCtrl,
                           par_units),
             error = function(e) optimizer_from_best(obj, par, e, optimizer,
                                                     bounds, control,
                                                     par_units, verbose))
  }
  if (verbose) t0 <- vb_now()
  opt <- run(start_par)
  if (verbose) vb_stage("optimize", t0, vb_opt_detail(opt))
  for (i in seq_len(control$restarts)) {
    g <- max(abs(obj$gr(opt$par) * (par_units %||% 1)))
    if (is.finite(g) && g < control$grad_tol) break
    if (verbose) t0 <- vb_now()
    opt2 <- run(opt$par)
    if (verbose) {
      vb_stage(paste0("restart ", i), t0,
               paste0(vb_opt_detail(opt2), ", from max|grad| ",
                      format(g, digits = 3)))
    }
    if (opt2$objective <= opt$objective) opt <- opt2
  }
  opt
}

#' Does an optimizer failure message read like an undefined objective?
#'
#' The diagnosis is only attached to errors that actually look numerical.
#' Everything raised inside the wrapped expression lands here, including
#' a misspelled optimizer, a custom optimizer that broke its return
#' contract, and errors from user code in a custom family - none of which
#' says anything about the likelihood. Reporting those as an undefined
#' objective sends the user to the wrong remedies and buries the real
#' message under advice.
#'
#' @noRd
fit_numerical_error <- function(msg) {
  # nlminb ("NA/NaN function evaluation"), RTMB ("NA/NaN gradient
  # evaluation") and optim's L-BFGS-B finiteness checks are the
  # vocabulary of an objective that came back undefined
  grepl("NA/NaN|NaN|non-?finite|not finite|infinite|Inf\\b", msg,
        ignore.case = TRUE)
}

#' Runs the optimizer call and rewrites any failure into a message that
#' names the model and what to try next.
#'
#' Whatever the optimizer throws reaches the user raw otherwise. RTMB's
#' "NA/NaN gradient evaluation" and nlminb's PORT strings name neither
#' the model nor anything to try, and they arrive with the call stack of
#' stats::nlminb, so the user cannot even tell which of their models
#' failed. Wrap the optimizer call once and say both.
#'
#' A nonlinear model gets the sharper advice, because there the cause is
#' nearly always the start: make_start() can only seed intercepts through
#' a family's init_dpars, and a nonlinear mu has no design of its own, so
#' an nl fit begins at zero, where most nonlinear forms are flat,
#' singular, or undefined. `[brms#734 doctrine]`
#'
#' The condition carries a class of its own, `frmtmb_fit_error`, so the
#' autoscale pre-fit (an inner `fit_assembled()` that has already been
#' through here) is rethrown rather than wrapped a second time.
#'
#' @noRd
fit_error_context <- function(spec, start, REML, control, quadrature,
                              priors, expr) {
  nl <- any(vapply(spec$responses,
                   function(r) length(r$nlpars) > 0L, TRUE))
  nl_start <- nl && is.null(start)
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      if (inherits(e, "frmtmb_fit_error")) stop(e)
      emsg <- conditionMessage(e)
      numerical <- fit_numerical_error(emsg)
      msg <- if (nl_start && numerical) {
        paste0("The nonlinear fit failed from the default zero starting ",
               "values (", emsg, "). Nonlinear models ",
               "need starting values in the right region: supply them ",
               "with the `start` argument of frm(), e.g. ",
               "start = list(beta = c(...)) in the order of fixef(); ",
               "frm(..., dry_run = \"frame\")$par_template shows the ",
               "layout")
      } else if (numerical) {
        paste0("The optimizer failed on this model (",
               vb_fit_detail(spec, REML, control, quadrature, priors),
               "): ", emsg,
               ". The likelihood was undefined or unbounded somewhere ",
               "the optimizer stepped. Refit with verbose = TRUE to see ",
               "which stage broke, try another optimizer ",
               "(frmtmb_control(optimizer = \"optim\")), or start the ",
               "fit nearer the optimum with the `start` argument of ",
               "frm()")
      } else {
        # the real message first, then only what is certainly true:
        # which model it came from
        paste0(emsg, " (raised while fitting: ",
               vb_fit_detail(spec, REML, control, quadrature, priors),
               ")")
      }
      stop(structure(
        class = c("frmtmb_fit_error", "simpleError", "error", "condition"),
        list(message = msg, call = NULL)))
    }),
    warning = function(w) {
      # nlminb's own "NA/NaN function evaluation" is the optimizer
      # noticing the same undefined objective the error above names;
      # letting both through would report the failure twice
      if (nl_start && grepl("NA/NaN", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' The cold starting values: the parameter template with each linear
#' predictor's intercept seeded from the family's own initializer, then
#' whatever the user gave in `start` written over the result. Component
#' names and lengths are checked here, where the error can still name
#' the template.
#'
#' @noRd
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

#' The post-fit verdict: the optimizer status, the maximum absolute
#' gradient, and, when a report is already there, the Hessian and the
#' standard errors. Each failure is a separate warning, so
#' `suppressWarnings()` silences them and the fit is still returned.
#'
#' @noRd
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
  # Covariance verdicts are only known once sdreport has run (se =
  # TRUE); the lazy path surfaces them through vcov()/summary()/
  # diagnose() instead. pdHess does not imply usable standard errors:
  # the Cholesky it comes from succeeds on a Hessian LAPACK's solver
  # then refuses as computationally singular, and cov.fixed is NaN.
  sdr <- fit$cache$sdr
  if (!is.null(sdr) && !is.null(sdr$pdHess) && !isTRUE(sdr$pdHess)) {
    msgs <- c(msgs, paste0("Hessian is not positive definite; standard ",
                           "errors are unreliable. The model may be ",
                           "overparameterized"))
  } else if (!is.null(sdr) && length(sdr$cov.fixed) &&
             any(!is.finite(sdr$cov.fixed))) {
    msgs <- c(msgs, paste0("Some standard errors are not finite: the ",
                           "covariance could not be recovered from the ",
                           "Hessian. diagnose() names the offending ",
                           "parameters; see the 'Convergence problems' ",
                           "section of vignette('diagnostics')"))
  }
  for (m in msgs) warning(m, call. = FALSE)
  invisible(list(grad = g, warnings = msgs))
}
