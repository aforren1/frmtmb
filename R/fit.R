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
                priors = NULL, quadrature = FALSE, dry_run = NULL) {
  cl <- match.call()
  bform <- if (inherits(formula, c("frmtmb_formula", "frmtmb_mvformula"))) {
    formula
  } else {
    bf(formula)
  }
  if (!is.null(family)) {
    fam <- as_frmtmb_family(family)
    if (inherits(bform, "frmtmb_mvformula")) {
      bform$forms <- lapply(bform$forms, function(f) {
        if (is.null(f$family)) f$family <- fam
        f
      })
    } else {
      bform$family <- fam
    }
  }

  spec <- parse_spec(bform)
  if (identical(dry_run, "spec")) return(spec)

  frame <- assemble_frame(spec, data, na.action = na.action)
  if (identical(dry_run, "frame")) return(frame)

  fit_assembled(spec, frame, bform, cl, REML = REML, start = start,
                control = control, se = se, lower = lower, upper = upper,
                priors = priors, quadrature = quadrature)
}

# Fitting core shared by frm() and refit(): objective build through the
# convergence check. A non-NULL `template` bypasses make_start (warm
# starts when refitting to a new response).
fit_assembled <- function(spec, frame, bform, cl, REML, start, control,
                          se, lower, upper, priors, quadrature,
                          template = NULL) {
  lower_arg <- lower
  upper_arg <- upper
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
  bounds <- resolve_bounds(list(frame = frame, REML = REML), lower, upper)
  opt <- optimize_obj(obj, control, bounds)

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
         lower = lower_arg, upper = upper_arg,
         cache = new.env(parent = emptyenv())),
    class = "frmtmb_fit"
  )
  if (se) {
    fit$cache$sdr <- RTMB::sdreport(obj, getJointPrecision = needs_jp(fit))
  }
  check_convergence(fit, control)
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
    cache$sdr <- RTMB::sdreport(fit$obj, getJointPrecision = needs_jp(fit))
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
#' @return A list of control settings.
#' @export
frmtmb_control <- function(optimizer = "nlminb",
                           optCtrl = list(iter.max = 1000, eval.max = 1000),
                           restarts = 1, grad_tol = 1e-3,
                           profile = FALSE) {
  list(optimizer = optimizer, optCtrl = optCtrl, restarts = restarts,
       grad_tol = grad_tol, profile = isTRUE(profile))
}

# One optimizer invocation, normalized to nlminb's result shape.
run_optimizer <- function(optimizer, par, fn, gr, lower, upper, control) {
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
                           lower = lower, upper = upper),
    optim = {
      ctl <- control[names(control) %in%
                       c("maxit", "factr", "pgtol", "trace")]
      if (is.null(ctl$maxit)) ctl$maxit <- 1000
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
                         bounds = list(lower = -Inf, upper = Inf)) {
  optimizer <- control$optimizer %||% "nlminb"
  opt <- run_optimizer(optimizer, obj$par, obj$fn, obj$gr,
                       bounds$lower, bounds$upper, control$optCtrl)
  for (i in seq_len(control$restarts)) {
    g <- max(abs(obj$gr(opt$par)))
    if (is.finite(g) && g < control$grad_tol) break
    opt2 <- run_optimizer(optimizer, opt$par, obj$fn, obj$gr,
                          bounds$lower, bounds$upper, control$optCtrl)
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
  if (fit$opt$convergence != 0) {
    warning("Optimizer did not report convergence: ", fit$opt$message,
            call. = FALSE)
  }
  g <- try(max(abs(fit$obj$gr(fit$opt$par))), silent = TRUE)
  if (!inherits(g, "try-error") && is.finite(g) && g > control$grad_tol) {
    warning("Large maximum absolute gradient at the optimum (",
            format(g, digits = 3), "); the fit may not have converged",
            call. = FALSE)
  }
  # pdHess is only known once sdreport has run (se = TRUE); the lazy
  # path surfaces it through summary()/diagnose() instead
  sdr <- fit$cache$sdr
  if (!is.null(sdr) && !is.null(sdr$pdHess) && !isTRUE(sdr$pdHess)) {
    warning("Hessian is not positive definite; standard errors are ",
            "unreliable. The model may be overparameterized",
            call. = FALSE)
  }
  invisible(fit)
}
