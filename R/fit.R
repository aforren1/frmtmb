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
#'   models only.
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
  obj <- RTMB::MakeADFun(nll, template, random = random,
                         map = frame$map, integrate = integrate,
                         profile = profile_arg, silent = TRUE)
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
  opt <- optimize_obj(obj, ctl_opt, bounds, par_units, verbose = vb)

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
                           autoscale = FALSE, verbose = NULL) {
  # verbose stays NULL when unset, which is how frm(verbose =) knows an
  # explicit control value must win over its own shortcut
  list(optimizer = optimizer, optCtrl = optCtrl, restarts = restarts,
       grad_tol = grad_tol, profile = isTRUE(profile),
       sparse_x = isTRUE(sparse_x), autoscale = isTRUE(autoscale),
       verbose = verbose)
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
    msgs <- c(msgs, paste0("Optimizer did not report convergence: ",
                           fit$opt$message))
  }
  # under autoscale the gradient is judged in the same natural units
  # the optimizer used (a 1e6-scale column bounds its coefficient's
  # absolute gradient near machine noise times 1e6)
  g <- try(max(abs(fit$obj$gr(fit$opt$par) * (fit$par_units %||% 1))),
           silent = TRUE)
  if (inherits(g, "try-error")) g <- NA_real_
  if (is.finite(g) && g > control$grad_tol) {
    msgs <- c(msgs, paste0("Large maximum absolute gradient at the ",
                           "optimum (", format(g, digits = 3),
                           "); the fit may not have converged"))
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
