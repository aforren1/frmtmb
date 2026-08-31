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
#' @export
frm <- function(formula, data, family = NULL, REML = FALSE, start = NULL,
                control = frmtmb_control(), se = FALSE,
                na.action = stats::na.omit, lower = NULL, upper = NULL,
                priors = NULL, dry_run = NULL) {
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
  template <- make_start(frame, start)

  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- if (!is.null(template[["b"]])) "b" else character(0)
  if (REML) random <- c(random, "beta")
  if (!length(random)) random <- NULL

  obj <- RTMB::MakeADFun(nll, template, random = random,
                         map = frame$map, silent = TRUE)
  bounds <- resolve_bounds(list(frame = frame, REML = REML), lower, upper)
  opt <- optimize_obj(obj, control, bounds)

  # Estimates come cheaply from the parameter list at the optimum;
  # sdreport (a quarter of typical fit time) is computed on demand
  # through sdr_of() unless se = TRUE asked for it now.
  est <- obj$env$parList(opt$par)
  for (nm in names(frame$par_template)) {
    names(est[[nm]]) <- names(frame$par_template[[nm]])
  }
  fit <- structure(
    list(spec = spec, frame = frame, obj = obj, opt = opt, sdr = NULL,
         REML = REML, estimates = est, priors = priors,
         bform = bform, call = cl,
         cache = new.env(parent = emptyenv())),
    class = "frmtmb_fit"
  )
  if (se) fit$cache$sdr <- RTMB::sdreport(obj, getJointPrecision = REML)
  check_convergence(fit, control)
  fit
}

# Memoized sdreport: the standard-error machinery (summary, vcov,
# confint, predict se.fit, diagnose) triggers it on first use.
sdr_of <- function(fit) {
  cache <- fit$cache
  if (is.null(cache$sdr)) {
    cache$sdr <- RTMB::sdreport(fit$obj, getJointPrecision = fit$REML)
  }
  cache$sdr
}

#' Control parameters for frmtmb fits
#'
#' @param optCtrl Control list passed to [stats::nlminb()].
#' @param restarts Number of times to restart the optimizer from the
#'   current optimum while the gradient remains above `grad_tol`.
#' @param grad_tol Warn (and restart) if the maximum absolute gradient at
#'   the optimum exceeds this value.
#' @return A list of control settings.
#' @export
frmtmb_control <- function(optCtrl = list(iter.max = 1000, eval.max = 1000),
                           restarts = 1, grad_tol = 1e-3) {
  list(optCtrl = optCtrl, restarts = restarts, grad_tol = grad_tol)
}

optimize_obj <- function(obj, control,
                         bounds = list(lower = -Inf, upper = Inf)) {
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = control$optCtrl,
                       lower = bounds$lower, upper = bounds$upper)
  for (i in seq_len(control$restarts)) {
    g <- max(abs(obj$gr(opt$par)))
    if (is.finite(g) && g < control$grad_tol) break
    opt2 <- stats::nlminb(opt$par, obj$fn, obj$gr,
                          control = control$optCtrl,
                          lower = bounds$lower, upper = bounds$upper)
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
