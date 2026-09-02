# Custom-family checking, MCMC bridge, emmeans glue.

#' Check a custom family's log-density for AD safety
#'
#' Tapes the family's `lpdf` on test values and compares the AD gradient
#' against central finite differences. A mismatch usually means the lpdf
#' uses operations the tape cannot see (base `matrix()`/`c()` on
#' advectors, branching on parameter values, `min`/`max`, clamping).
#'
#' @param family A `frmtmb_family` (from [frmtmb_family()] /
#'   [custom_family()]).
#' @param y A response vector of test data.
#' @param dpars Named list of numeric test values, one entry per dpar
#'   (each of length 1 or `length(y)`).
#' @param aterms Named list of addition-term values (e.g. `trials`).
#' @param tol Maximum relative gradient error.
#' @return Invisibly `TRUE`; signals an error on failure.
#' @examples
#' set.seed(1)
#' y <- rpois(50, 3)
#'
#' # a hand-written poisson: check it before fitting anything with it
#' ok <- custom_family(
#'   "my_poisson", dpars = "mu", links = list(mu = "log"),
#'   lpdf = function(y, dpars, aterms) {
#'     y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
#'   },
#'   type = "discrete"
#' )
#' check_custom_family(ok, y = y, dpars = list(mu = rep(2.5, 50)))
#'
#' # base matrix() strips the advector class, so the tape sees constants
#' # and the gradient is silently wrong. The check catches it.
#' bad <- custom_family(
#'   "bad", dpars = "mu", links = list(mu = "log"),
#'   lpdf = function(y, dpars, aterms) {
#'     m <- matrix(dpars$mu, ncol = 1)
#'     y * log(m[, 1]) - m[, 1] - lgamma(y + 1)
#'   },
#'   type = "discrete"
#' )
#' try(check_custom_family(bad, y = y, dpars = list(mu = rep(2.5, 50))))
#' @export
check_custom_family <- function(family, y, dpars, aterms = list(),
                                tol = 1e-4) {
  stopifnot(inherits(family, "frmtmb_family"))
  if (!setequal(names(dpars), family$dpars)) {
    stop("`dpars` must supply test values for exactly: ",
         paste(family$dpars, collapse = ", "), call. = FALSE)
  }
  f <- function(p) -sum(family$lpdf(y, p, aterms))
  v0 <- f(lapply(dpars, as.numeric))
  if (!is.finite(v0)) {
    stop("lpdf is not finite at the test values", call. = FALSE)
  }
  obj <- tryCatch(
    RTMB::MakeADFun(f, dpars, silent = TRUE),
    error = function(e) {
      stop("Failed to tape the lpdf: ", conditionMessage(e),
           ". Typical cause: base matrix()/c() stripping the advector ",
           "class, or branching on parameter values", call. = FALSE)
    }
  )
  if (abs(obj$fn(obj$par) - v0) > 1e-8 * max(1, abs(v0))) {
    stop("Taped lpdf disagrees with its plain-numeric value (",
         format(obj$fn(obj$par)), " vs ", format(v0), "): the lpdf uses ",
         "operations that behave differently on the AD tape (base ",
         "matrix()/c() on advectors are the usual culprits)",
         call. = FALSE)
  }
  g <- as.vector(obj$gr(obj$par))
  p0 <- obj$par
  h <- 1e-6 * pmax(abs(p0), 1)
  fd <- vapply(seq_along(p0), function(i) {
    pp <- p0; pp[i] <- pp[i] + h[i]
    pm <- p0; pm[i] <- pm[i] - h[i]
    (obj$fn(pp) - obj$fn(pm)) / (2 * h[i])
  }, numeric(1))
  rel <- abs(g - fd) / pmax(abs(fd), 1)
  if (any(rel > tol)) {
    stop("AD gradient disagrees with finite differences (max relative ",
         "error ", format(max(rel), digits = 3), ")", call. = FALSE)
  }
  invisible(TRUE)
}

#' Prior specifications for frm_sample
#'
#' Priors apply on the INTERNAL parameter scale: coefficients are on
#' their link scale, and covariance parameters (`theta_*`) are the
#' unconstrained parameterization (log-SDs, scaled-Cholesky terms), so
#' `prior_normal(0, 1)` on `theta_1` is a lognormal prior on that SD.
#'
#' @param location,scale,df Prior parameters.
#' @return A `frmtmb_prior` object.
#' @examples
#' # the objects themselves are cheap descriptions
#' prior_normal(0, 2)
#' prior_t(df = 3, location = 0, scale = 1)
#'
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # names are parameter names as in the draws, or whole components.
#' # theta_1 is a log-SD, so a normal there is a lognormal on the SD.
#' ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0,
#'                  priors = list(beta = prior_normal(0, 5),
#'                                theta_1 = prior_t(3, 0, 1)))
#' summary(ds)
#' }
#' }
#' @name frmtmb-priors
NULL

#' @rdname frmtmb-priors
#' @export
prior_normal <- function(location = 0, scale = 1) {
  stopifnot(scale > 0)
  structure(list(kind = "normal", location = location, scale = scale),
            class = "frmtmb_prior")
}

#' @rdname frmtmb-priors
#' @export
prior_t <- function(df = 3, location = 0, scale = 1) {
  stopifnot(df > 0, scale > 0)
  structure(list(kind = "t", df = df, location = location, scale = scale),
            class = "frmtmb_prior")
}

#' Resolve a named prior list to per-component index/parameter vectors.
#' Names may be individual parameters (as in outer_par_names()) or whole
#' components ("beta", "betad", "theta", "thetar").
#'
#' @noRd
resolve_priors <- function(fit, priors) {
  stopifnot(is.list(priors), !is.null(names(priors)))
  tpl <- fit$frame$par_template
  comp_names <- list()
  for (cp in setdiff(names(tpl), c("b", "miss"))) {
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    if (cp == "betad" && length(fit$frame$betad_fixed_idx)) {
      v[fit$frame$betad_fixed_idx] <- NA   # mapped: no prior
    }
    comp_names[[cp]] <- v
  }
  entries <- list()
  add <- function(comp, idx, pr) {
    entries[[length(entries) + 1L]] <<- list(comp = comp, idx = idx,
                                             prior = pr)
  }
  for (nm in names(priors)) {
    pr <- priors[[nm]]
    if (!inherits(pr, "frmtmb_prior")) {
      stop("priors[['", nm, "']] must be a prior object ",
           "(prior_normal(), prior_t())", call. = FALSE)
    }
    if (nm %in% names(comp_names)) {
      idx <- which(!is.na(comp_names[[nm]]))
      add(nm, idx, pr)
      next
    }
    hit <- FALSE
    for (cp in names(comp_names)) {
      i <- which(comp_names[[cp]] == nm)
      if (length(i)) {
        add(cp, i, pr)
        hit <- TRUE
        break
      }
    }
    if (!hit) {
      stop("Unknown parameter in priors: '", nm, "'. Available: ",
           paste(unlist(comp_names)[!is.na(unlist(comp_names))],
                 collapse = ", "),
           " or component names ",
           paste(names(comp_names), collapse = ", "), call. = FALSE)
    }
  }
  entries
}

#' AD-safe negative log prior over resolved per-parameter entries
#' (each: comp, idx, dist, scale; see prior_logdens).
#'
#' @noRd
neg_log_prior_fn <- function(entries) {
  function(pars) {
    nlp <- 0
    for (e in entries) {
      nlp <- nlp - sum(prior_logdens(pars[[e$comp]][e$idx], e$dist,
                                     e$scale))
    }
    nlp
  }
}

#' Named bound specs -> full-length vectors over the outer parameters.
#'
#' @noRd
resolve_bounds <- function(fit, lower, upper) {
  nm <- outer_par_names(fit)
  mk <- function(x, fill) {
    out <- rep(fill, length(nm))
    if (is.null(x)) return(out)
    if (is.null(names(x)) || any(names(x) == "")) {
      stop("Bounds must be named numeric vectors, e.g. ",
           "lower = c(x = 0)", call. = FALSE)
    }
    bad <- setdiff(names(x), nm)
    if (length(bad)) {
      stop("Unknown parameter(s) in bounds: ",
           paste(bad, collapse = ", "), ". Available: ",
           paste(nm, collapse = ", "), call. = FALSE)
    }
    out[match(names(x), nm)] <- as.numeric(x)
    out
  }
  list(lower = mk(lower, -Inf), upper = mk(upper, Inf))
}

#' Pull a start value strictly inside the bounding box. Stan turns a
#' bound into a constrained transform, and a start AT the bound maps to
#' an infinite unconstrained value: rstan then reports "Initialization
#' failed" and names neither the parameter nor the bound. The pad is
#' relative so it survives bounds of any magnitude; a box narrower than
#' two pads collapses to its midpoint.
#'
#' @noRd
clamp_into_bounds <- function(v, lower, upper) {
  if (is.null(lower) && is.null(upper)) return(v)
  lo <- lower %||% rep(-Inf, length(v))
  hi <- upper %||% rep(Inf, length(v))
  pad_lo <- 1e-3 * pmax(1, abs(lo))
  pad_hi <- 1e-3 * pmax(1, abs(hi))
  fl <- ifelse(is.finite(lo), lo + pad_lo, -Inf)
  ce <- ifelse(is.finite(hi), hi - pad_hi, Inf)
  narrow <- is.finite(lo) & is.finite(hi) & fl > ce
  mid <- (lo + hi) / 2
  fl[narrow] <- mid[narrow]
  ce[narrow] <- mid[narrow]
  pmin(pmax(v, fl), ce)
}

#' Per-chain initial values around the ML mode: chain 1 exactly at the
#' mode (the short-warmup anchor), later chains at mode + N(0, jitter)
#' on the unconstrained scale, restoring the overdispersion Rhat needs
#' to detect chains agreeing for the wrong reason. Every chain - the
#' mode-anchored one included, because the bounds can exclude the ML
#' mode itself - is then pulled inside the bounds.
#'
#' @noRd
mode_inits <- function(mode, chains, jitter, lower = NULL, upper = NULL) {
  mode <- as.numeric(mode)
  if (!is.finite(jitter) || jitter <= 0 || chains <= 1L) {
    v <- clamp_into_bounds(mode, lower, upper)
    return(lapply(seq_len(max(chains, 1L)), function(i) v))
  }
  lapply(seq_len(chains), function(i) {
    v <- if (i == 1L) mode else mode + stats::rnorm(length(mode), 0, jitter)
    clamp_into_bounds(v, lower, upper)
  })
}

#' tmbstan widens outer-length bounds over the whole parameter vector
#' when the objective has random effects (the inner block is unbounded);
#' the inits are the full vector too, so they must be clamped against
#' the same widening.
#'
#' @noRd
mode_aligned_bounds <- function(obj, bounds, laplace, n) {
  if (is.null(bounds)) return(NULL)
  if (length(bounds$lower) == n) return(bounds)
  rnd <- obj$env$random
  if (laplace || is.null(rnd) || length(rnd) + length(bounds$lower) != n) {
    return(NULL)
  }
  lo <- rep(-Inf, n)
  hi <- rep(Inf, n)
  lo[-rnd] <- bounds$lower
  hi[-rnd] <- bounds$upper
  list(lower = lo, upper = hi)
}

#' Retape the fit's objective with priors added; parameters start at the
#' ML estimates so sampling initializes at (near) the posterior mode.
#'
#' @noRd
prior_augmented_obj <- function(fit, entries) {
  nll <- build_objective(fit$frame)
  nlp <- neg_log_prior_fn(entries)
  tpl <- fit$frame$par_template
  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- c(if (!is.null(tpl[["b"]])) "b",
              if (!is.null(tpl[["miss"]])) "miss")
  if (!length(random)) random <- NULL
  RTMB::MakeADFun(function(pars) nll(pars) + nlp(pars),
                  fit$estimates, random = random,
                  map = fit$frame$map, silent = TRUE)
}

#' Labels for the full sampled parameter vector, in template order,
#' skipping mapped entries; b kept as `b[i]`. `include_random = FALSE`
#' drops the inner components (b, miss) for laplace-marginalized draws.
#'
#' @noRd
all_par_labels <- function(fit, include_b = TRUE, include_random = TRUE) {
  tpl <- fit$frame$par_template
  out <- character(0)
  for (cp in names(tpl)) {
    v <- names(tpl[[cp]])
    if (is.null(v)) v <- paste0(cp, "_", seq_along(tpl[[cp]]))
    if (cp == "betad" && length(fit$frame$betad_fixed_idx)) {
      v <- v[-fit$frame$betad_fixed_idx]
    }
    if (cp == "miss" && !include_random) next
    if (cp == "b") {
      if (!include_b || !include_random) next
      v <- paste0("b[", seq_along(tpl[[cp]]), "]")
    }
    out <- c(out, v)
  }
  out
}

#' The core count tmbstan would actually use. rstan defaults its `cores`
#' argument to getOption("mc.cores", 1L), so a session-wide mc.cores
#' starts parallel chains just as an explicit cores = does, and the
#' Windows guard has to read the same source or it never fires.
#'
#' @noRd
stan_cores <- function(args) {
  as.numeric(args$cores %||% getOption("mc.cores", 1L))[1L]
}

#' Sample the fitted model with NUTS
#'
#' Runs [tmbstan::tmbstan()] on the fitted objective, initialized at the
#' ML estimates (which shortens warmup considerably), and returns the
#' draws with frmtmb coefficient names. Without priors this samples the
#' likelihood with flat improper priors on the outer parameters - the
#' random effects get their proper hierarchical Gaussian terms - so
#' treat the result as an ML diagnostic (see [check_laplace()]) rather
#' than a full Bayesian analysis; posteriors can be improper for
#' variance components with few groups.
#'
#' @param fit A `frmtmb_fit`.
#' @param ... Passed to [tmbstan::tmbstan()] (`chains`, `iter`,
#'   `laplace`, `cores`, ...). On Windows more than one core falls back
#'   to sequential chains with a warning: parallel chains run on socket
#'   workers, which can evaluate neither the RTMB tape nor the
#'   objective closure (the known RTMB limitation of tmbstan,
#'   tmbstan#27). The fallback also covers a core count inherited from
#'   `options(mc.cores)`, which is what rstan reads when `cores` is not
#'   given. Fork clusters on unix can, so `cores` works there.
#' @param priors Optional named list of priors (see [prior_normal()]);
#'   names are parameter names as in the draws (or whole components:
#'   `"beta"`, `"theta"`, ...). Parameters without a prior keep the flat
#'   improper default. The objective is re-taped with the prior terms
#'   added; the ML fit itself is unchanged.
#' @param lower,upper Optional named numeric vectors of hard bounds on
#'   outer parameters (brms `lb`/`ub`), applied on the internal scale
#'   through Stan's constrained transforms. Chain starting values are
#'   clamped strictly inside the bounds; a bound that excludes the ML
#'   mode itself warns, because the chains then no longer start there.
#' @param init Initialization; the default starts chain 1 exactly at
#'   the ML mode and every further chain at the mode plus a normal
#'   perturbation of sd `init_jitter` on the unconstrained scale.
#'   The mode anchor keeps warmup short; the jitter keeps the chains
#'   overdispersed enough for Rhat to retain power against
#'   multimodality (the standard objection to identical mode starts).
#'   `"random"` requests Stan's own overdispersed initialization.
#' @param init_jitter Per-chain perturbation sd for the default init;
#'   `0` starts every chain exactly at the mode. Draws from the R
#'   session's RNG, so `set.seed()` makes the inits reproducible.
#' @return An object of class `frmtmb_draws`: list with the `stanfit`,
#'   a draws matrix with named columns (`as.matrix()` method), and the
#'   originating fit.
#' @section Multimodal posteriors:
#'   For [mixture()] fits the posterior is multimodal by construction
#'   (label switching at minimum). Mode-centered inits, jittered or
#'   not, leave every chain in one symmetry branch, so Rhat cannot
#'   flag the others; use `init = "random"` there and inspect chains
#'   individually.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#' ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
#' summary(ds)
#' fixef(ds)
#' hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")
#' }
#' }
#' @export
frm_sample <- function(fit, ..., priors = NULL, lower = NULL,
                       upper = NULL, init = "last.par.best",
                       init_jitter = 0.25) {
  if (!requireNamespace("tmbstan", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    stop("frm_sample() needs the 'tmbstan' and 'rstan' packages",
         call. = FALSE)
  }
  pr_lower <- c()
  pr_upper <- c()
  obj <- fit$obj
  # a MAP fit carries its priors into sampling unless overridden
  priors <- priors %||% fit$priors
  if (!is.null(priors)) {
    ri <- resolve_prior_input(fit, priors)
    if (length(ri$entries)) obj <- prior_augmented_obj(fit, ri$entries)
    pr_lower <- ri$lower
    pr_upper <- ri$upper
  }
  # explicit lower/upper override set_prior() bounds on overlap
  lower <- utils::modifyList(as.list(pr_lower),
                             as.list(lower %||% c()))
  upper <- utils::modifyList(as.list(pr_upper),
                             as.list(upper %||% c()))
  laplace <- isTRUE(list(...)$laplace)
  bounds <- if (length(lower) || length(upper)) {
    resolve_bounds(fit, unlist(lower), unlist(upper))
  }
  if (identical(init, "last.par.best")) {
    # a singular ML mode (variance at the boundary) is exactly the
    # pathological start the mode-init criticism is about. Only the
    # log-sd components read that way: a CAR mixing proportion or an
    # AR(1) phi at its own boundary is a large theta on a converged,
    # perfectly samplable fit (see log_sd_theta_index()).
    th_all <- fit$estimates$theta %||% numeric(0)
    sd_i <- log_sd_theta_index(fit)
    ext <- sd_i[abs(th_all[sd_i]) > 8]
    if (length(ext)) {
      warning("The ML mode has an extreme covariance parameter (",
              paste(names(ext), collapse = ", "),
              "; likely a boundary/singular fit); mode initialization ",
              "starts the chains there. Consider init = \"random\", ",
              "or regularize with priors =", call. = FALSE)
    }
    lpb <- obj$env$last.par.best
    rnd <- obj$env$random
    # under laplace tmbstan samples only the outer parameters, so the
    # full-length mode has the wrong length; take the outer slice
    mode <- if (laplace && length(rnd)) lpb[-rnd] else lpb
    if (!is.null(bounds)) {
      # a bound tighter than the ML estimate makes the mode inadmissible
      # as a start; the chains then begin somewhere else, and a
      # posterior pinned against a bound is a modeling statement worth
      # hearing about
      om <- if (laplace || !length(rnd)) mode else mode[-rnd]
      viol <- which(as.numeric(om) < bounds$lower |
                      as.numeric(om) > bounds$upper)
      if (length(viol)) {
        warning("The ML mode violates the requested bound(s) on ",
                paste(outer_par_names(fit)[viol], collapse = ", "),
                "; every chain starts at the clamped value instead of ",
                "the mode. The bounded posterior is not centered on the ",
                "unconstrained ML estimate", call. = FALSE)
      }
    }
    mb <- mode_aligned_bounds(obj, bounds, laplace, length(mode))
    init <- mode_inits(mode, list(...)$chains %||% 4, init_jitter,
                       mb$lower, mb$upper)
  }
  args <- list(obj = obj, init = init, ...)
  # rstan runs parallel chains on PSOCK workers on Windows, and neither
  # the RTMB tape's external pointer nor the objective closure's
  # namespace survives the trip: every chain dies at the first internal
  # function and rstan's error names neither. Known upstream as
  # tmbstan#27; fork clusters (unix) inherit both, so only Windows
  # needs the guard.
  if (.Platform$OS.type == "windows" && stan_cores(args) > 1) {
    warning("cores > 1 is not available on Windows: parallel chains ",
            "run on socket workers, which cannot evaluate the RTMB ",
            "tape. Running the chains sequentially", call. = FALSE)
    args$cores <- 1
  }
  if (!is.null(bounds)) {
    args$lower <- bounds$lower
    args$upper <- bounds$upper
  }
  sf <- do.call(tmbstan::tmbstan, args)
  a <- rstan::extract(sf, permuted = FALSE)   # iter x chain x par
  stan_names <- dimnames(a)[[3]]
  # laplace draws skip the inner components entirely; labeling them
  # with the full template order would misattribute theta as b[i]
  labels <- all_par_labels(fit, include_random = !laplace)
  n_lab <- min(length(labels), length(stan_names))
  stan_names[seq_len(n_lab)] <- labels[seq_len(n_lab)]
  m <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
  colnames(m) <- stan_names
  structure(list(stanfit = sf, draws = m, fit = fit),
            class = "frmtmb_draws")
}

#' @export
as.matrix.frmtmb_draws <- function(x, ...) x$draws

#' @export
print.frmtmb_draws <- function(x, ...) {
  m <- x$draws
  keep <- setdiff(colnames(m),
                  c("lp__", grep("^b\\[", colnames(m), value = TRUE)))
  cat("<frmtmb_draws> ", nrow(m), " draws x ", ncol(m),
      " parameters\n\n", sep = "")
  tab <- t(vapply(keep, function(nm) {
    c(mean = mean(m[, nm]), sd = stats::sd(m[, nm]),
      `2.5%` = unname(stats::quantile(m[, nm], 0.025)),
      `97.5%` = unname(stats::quantile(m[, nm], 0.975)))
  }, numeric(4)))
  print(signif(tab, 4))
  invisible(x)
}

#' Check the Laplace/Wald approximation against NUTS
#'
#' Samples the fitted objective (see [frm_sample()]) and compares the ML
#' estimates and sdreport standard errors against posterior means and
#' SDs. Close agreement supports the Laplace approximation and Wald
#' intervals; a posterior SD much larger than the Wald SE, or a shifted
#' mean, flags parameters where they are unreliable (typically variance
#' components with few groups).
#'
#' @param fit A `frmtmb_fit`.
#' @param chains,iter Passed to [frm_sample()].
#' @param ... Passed to [frm_sample()].
#' @return A data frame (one row per outer parameter) with columns
#'   `ml`, `post_mean`, `wald_se`, `post_sd`, `z_shift`
#'   ((post_mean - ml)/post_sd) and `sd_ratio` (post_sd/wald_se).
#'
#' @srrstats {RE1.4} The assumptions the fit rests on are documented, and
#'   the consequences of violating them are both documented and testable.
#'   Random effects are integrated out by the Laplace approximation,
#'   which assumes the integrand is close to Gaussian around the
#'   conditional mode, and Wald intervals assume the log-likelihood is
#'   close to quadratic at the optimum. Both degrade in the same places:
#'   variance components estimated from few groups, and binary data in
#'   small clusters. `check_laplace()` measures the violation directly by
#'   running NUTS on the same objective and reporting the shift of the
#'   posterior mean from the maximum-likelihood estimate in posterior
#'   standard deviations (`z_shift`) and the ratio of the posterior
#'   standard deviation to the Wald standard error (`sd_ratio`), warning
#'   when either leaves its tolerance. The documented remedies are
#'   `confint(method = "profile")`, `frm_bootstrap()`, and
#'   `frm(quadrature = TRUE)`, which replaces the approximation with
#'   adaptive quadrature; the test suite checks that the quadrature fit
#'   matches `lme4::glmer(nAGQ = 25)` and GLMMadaptive in exactly the
#'   regime where the Laplace fit is biased.
#'   `vignette("diagnostics")` works through this.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#' # a binary GLMM with small clusters: the regime where the Laplace
#' # approximation and Wald intervals are least reliable
#' set.seed(4)
#' dd <- data.frame(x = rnorm(120), g = factor(rep(1:30, 4)))
#' dd$y <- rbinom(120, 1,
#'                plogis(0.3 + 0.5 * dd$x + rnorm(30, 0, 1)[dd$g]))
#' fit <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)
#'
#' cl <- check_laplace(fit, chains = 1, iter = 500, refresh = 0)
#' cl
#' # |z_shift| well above 0 or sd_ratio far from 1 marks the parameters
#' # whose Wald interval to replace with a profile or bootstrap one
#' cl[abs(cl$z_shift) > 0.3 | cl$sd_ratio > 1.3, ]
#' }
#' }
#' @export
check_laplace <- function(fit, chains = 2, iter = 1000, ...) {
  ds <- frm_sample(fit, chains = chains, iter = iter, ...)
  m <- ds$draws
  keep <- setdiff(colnames(m),
                  c("lp__", grep("^b\\[", colnames(m), value = TRUE)))
  ml <- fit$opt$par
  se <- sqrt(diag(sdr_of(fit)$cov.fixed))
  stopifnot(length(ml) == length(keep))
  post_mean <- colMeans(m[, keep, drop = FALSE])
  post_sd <- apply(m[, keep, drop = FALSE], 2, stats::sd)
  out <- data.frame(
    parameter = keep,
    ml = unname(ml),
    post_mean = unname(post_mean),
    wald_se = unname(se),
    post_sd = unname(post_sd),
    z_shift = unname((post_mean - ml) / post_sd),
    sd_ratio = unname(post_sd / se),
    row.names = NULL
  )
  flagged <- abs(out$z_shift) > 0.5 | out$sd_ratio > 1.5 |
    out$sd_ratio < 2 / 3
  if (any(flagged)) {
    message("Laplace/Wald approximation questionable for: ",
            paste(out$parameter[flagged], collapse = ", "))
  }
  out
}

#' Sample from a frmtmb fit with tmbstan (NUTS)
#'
#' Hands the fitted RTMB object to [tmbstan::tmbstan()]; from Stan's
#' point of view the model is the (Laplace-free) joint density, so all
#' parameters including random effects are sampled unless `laplace =
#' TRUE` is passed through.
#'
#' @param fit A `frmtmb_fit`.
#' @param ... Passed to [tmbstan::tmbstan()] (chains, iter, laplace, ...).
#' @return A `stanfit` object.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE)) {
#' set.seed(9)
#' dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
#' dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # The raw stanfit, for rstan and bayesplot code that wants one.
#' # Use frm_sample() instead when you want frmtmb parameter names.
#' # This run is deliberately short, so expect sampler warnings.
#' sf <- as_tmbstan(fit, chains = 1, iter = 400, refresh = 0)
#' class(sf)
#' # every parameter is sampled, random effects included, because Stan
#' # sees the joint density. Pass laplace = TRUE to integrate them out.
#' dim(as.matrix(sf))
#' }
#' }
#' @export
as_tmbstan <- function(fit, ...) {
  if (!requireNamespace("tmbstan", quietly = TRUE)) {
    stop("as_tmbstan() needs the 'tmbstan' package", call. = FALSE)
  }
  tmbstan::tmbstan(fit$obj, ...)
}

#' emmeans support: registered in .onLoad, only for the parametric fixed
#' part of the mu linear predictor of univariate, non-nl fits.
#'
#' An ordinal fit therefore lands on the LATENT scale, which is
#' emmeans's own `mode = "latent"` convention for `clm`-like models: the
#' K-1 thresholds are not coefficients of this design, so the basis has
#' no intercept and the marginal means are the latent predictor without
#' a threshold offset. Contrasts are unaffected by that offset;
#' category probabilities are a different question, answered by
#' `predict(type = "response")` and `conditional_effects()`.
#'
#' @noRd
emm_mu_linpred <- function(object) {
  if (length(object$spec$responses) > 1) {
    stop("emmeans support is univariate-only for now", call. = FALSE)
  }
  rspec <- object$spec$responses[[1]]
  lp <- object$frame$linpreds[[linpred_key(rspec$resp_name, "mu")]]
  if (is.null(lp) || !is.null(lp$nl_body)) {
    stop("emmeans support needs a linear mu predictor", call. = FALSE)
  }
  lp
}

# marginaleffects support: the four extension generics plus the
# class-whitelist option set in .onLoad. Predictions under set_coef are
# conditional on the estimated random-effect modes (the glmmTMB/lmer
# convention for delta-method slopes).

#' @exportS3Method marginaleffects::get_coef
get_coef.frmtmb_fit <- function(model, ...) {
  est <- model$estimates
  bd <- est$betad
  if (length(model$frame$betad_fixed_idx)) {
    bd <- bd[-model$frame$betad_fixed_idx]
  }
  stats::setNames(c(est$beta, bd), estimated_coef_names(model))
}

#' @exportS3Method marginaleffects::set_coef
set_coef.frmtmb_fit <- function(model, coefs, ...) {
  tpl <- model$frame$par_template
  nb <- length(tpl$beta)
  model$estimates$beta[] <- coefs[seq_len(nb)]
  if (!is.null(tpl$betad)) {
    keep <- setdiff(seq_along(tpl$betad), model$frame$betad_fixed_idx)
    model$estimates$betad[keep] <- coefs[nb + seq_along(keep)]
  }
  model
}

#' @exportS3Method marginaleffects::get_vcov
get_vcov.frmtmb_fit <- function(model, ...) {
  vcov(model)
}

#' @exportS3Method marginaleffects::get_predict
get_predict.frmtmb_fit <- function(model, newdata, type = "response",
                                   ...) {
  type <- if (identical(type, "link")) "link" else "response"
  p <- predict(model, newdata = newdata, type = type)
  if (is.matrix(p)) {
    # A categorical outcome predicts a DISTRIBUTION per row (an ordinal
    # family's K category probabilities, a multinomial's D cell means),
    # so one newdata row is several predictions. marginaleffects keys
    # those with a `group` column and repeats the rowid; without it the
    # flattened matrix was handed back as n * K unrelated rows numbered
    # 1..nK, which silently misaligns every downstream contrast.
    g <- colnames(p) %||% as.character(seq_len(ncol(p)))
    return(data.frame(
      rowid = rep(seq_len(nrow(p)), times = ncol(p)),
      group = rep(g, each = nrow(p)),
      estimate = as.numeric(p)
    ))
  }
  data.frame(rowid = seq_along(p), estimate = as.numeric(p))
}

#' @exportS3Method emmeans::recover_data
recover_data.frmtmb_fit <- function(object, ..., data = NULL) {
  lp <- emm_mu_linpred(object)
  emmeans::recover_data(object$call,
                        stats::delete.response(lp$terms),
                        na.action = NULL,
                        data = data %||% model.frame(object), ...)
}

#' @exportS3Method emmeans::emm_basis
emm_basis.frmtmb_fit <- function(object, trms, xlev, grid, ...) {
  lp <- emm_mu_linpred(object)
  m <- stats::model.frame(trms, grid, na.action = stats::na.pass,
                          xlev = xlev)
  X <- stats::model.matrix(trms, m, contrasts.arg = lp$contrasts)
  # The fitted design is not always model.matrix()'s: an ordinal family
  # drops the intercept (the K-1 thresholds take its place), and a
  # rank-deficient fit drops aliased columns. Select the fitted columns
  # by name, or the basis is not conformable with bhat and emmeans fails
  # with "Non-conformable elements in reference grid".
  pn <- lp$param_colnames[seq_len(lp$n_param_cols)]
  if (!identical(colnames(X), pn)) {
    keep <- match(pn, colnames(X))
    if (anyNA(keep)) {
      stop("emmeans support cannot rebuild the fitted design: column(s) ",
           paste(pn[is.na(keep)], collapse = ", "),
           " are missing from the reference grid", call. = FALSE)
    }
    X <- X[, keep, drop = FALSE]
  }
  idx <- lp$idx[seq_len(lp$n_param_cols)]
  bhat <- object$estimates[[lp$par]][idx]
  V <- vcov(object)[idx, idx, drop = FALSE]
  list(X = X, bhat = bhat, nbasis = matrix(NA), V = V,
       dffun = function(k, dfargs) Inf, dfargs = list(), misc = list())
}

# --- lme4::getME ------------------------------------------------------
# Only the pieces that mean the same thing here. The vocabulary stays
# small on purpose: a name that would have to be faked (Lambdat, u, the
# lme4 sparse-Cholesky machinery) is worse than a name that errors,
# because downstream code cannot tell a wrong answer from a right one.

frmtmb_getME_vocab <- c("X", "Z", "Zt", "beta", "fixef", "b", "theta",
                        "lower", "sigma", "flist", "n_rtrms",
                        "n_rfacs")

#' Blocks that carry a real grouping factor. Smooths and Gaussian
#' processes are stored as random-effect blocks too, but their "levels"
#' are basis functions, so they have no factor to report.
#'
#' @noRd
getME_group_blocks <- function(object) {
  Filter(function(bk) !bk$covstruct %in% c("smooth", "gp", "hsgp") &&
           !is.null(bk$components[[1L]]$bar),
         object$frame$re_blocks)
}

#' Labels for the random-effect coefficient vector, which is also the
#' column order of every Z: level-major within a block, one entry per
#' (level, term coefficient), the way lme4 labels Zt rows.
#'
#' @noRd
re_coef_labels <- function(frame) {
  if (!length(frame$re_blocks)) return(character(0))
  unlist(lapply(frame$re_blocks, function(bk) {
    lv <- bk$levels %||% as.character(seq_len(bk$n_levels))
    paste0(rep(lv, each = bk$dim), ".", rep(bk$cnms, bk$n_levels))
  }), use.names = FALSE)
}

#' The grouping factors as they were at fit time, rebuilt from the
#' stored model frame the same way predict() rebuilds them for newdata.
#'
#' @noRd
getME_flist <- function(object) {
  out <- list()
  for (bk in getME_group_blocks(object)) {
    comp <- bk$components[[1L]]
    lp <- object$frame$linpreds[[comp$lp_key]]
    env <- object$spec$responses[[lp$resp]]$formula_env
    gv <- tryCatch(
      as.character(eval(comp$bar[[3L]], object$frame$data_frame, env)),
      error = function(e) NULL
    )
    if (length(gv) != object$frame$n_obs) {
      stop("getME(\"flist\"): cannot rebuild the grouping factor for `",
           bk$term_label, "`", call. = FALSE)
    }
    nm <- bk$group_name
    if (is.null(out[[nm]])) out[[nm]] <- factor(gv, levels = bk$levels)
  }
  out
}

#' Extract components of a fit, lme4 style
#'
#' A small [lme4::getME()] vocabulary, for downstream code written
#' against merMod objects. Registered on lme4's generic, so call it as
#' `lme4::getME(fit, "X")` or load lme4 first.
#'
#' Supported names:
#' \describe{
#'   \item{`"X"`}{Fixed-effect design matrix of the `mu` predictor.}
#'   \item{`"Z"`, `"Zt"`}{The sparse random-effect design of the `mu`
#'     predictor and its transpose. Columns (rows of `Zt`) span the
#'     whole random-effect coefficient vector, so a block belonging to
#'     another distributional parameter contributes zero columns here.}
#'   \item{`"beta"`, `"fixef"`}{The primary (`mu`-family) fixed-effect
#'     coefficients, named. Coefficients of auxiliary distributional
#'     parameters are a separate vector; use [fixef()] for all of them.}
#'   \item{`"b"`}{Conditional modes in coefficient space, aligned with
#'     the columns of `Z`. Reduced-rank (`rr()`) blocks are expanded
#'     through their loadings, so this is not the internal parameter
#'     vector.}
#'   \item{`"theta"`}{Covariance parameters on the internal
#'     (unconstrained) scale, as in `confint()`. These are not lme4's
#'     relative-covariance-factor entries.}
#'   \item{`"lower"`}{Lower bounds on `theta`. The internal
#'     parameterization is unbounded, so this is a vector of `-Inf`, not
#'     lme4's mixture of `0` and `-Inf`. Code that tests
#'     `theta == lower` to detect a singular fit will never fire; use
#'     [diagnose()] or [VarCorr()] instead.}
#'   \item{`"sigma"`}{Residual standard deviation
#'     ([sigma.frmtmb_fit()]).}
#'   \item{`"flist"`}{The grouping factors, one per distinct grouping
#'     variable. Smooth and Gaussian-process blocks are excluded: their
#'     levels are basis functions, not groups. There is no `"assign"`
#'     attribute.}
#'   \item{`"n_rtrms"`, `"n_rfacs"`}{Number of random-effect terms and
#'     of distinct grouping factors.}
#' }
#'
#' Multivariate fits have one design per response, so `"X"`, `"Z"` and
#' `"Zt"` need `resp`; the other names answer without it.
#'
#' @param object A `frmtmb_fit`.
#' @param name One or more names from the vocabulary above. A vector
#'   returns a named list.
#' @param resp Response name, for the design extractors on a
#'   multivariate fit.
#' @param ... Unused.
#' @return The requested component, or a named list when `name` names
#'   several.
#'
#' @srrstats {RE4.13} Predictor variables and their metadata are
#'   retrievable from the fitted object. `getME()` exposes the
#'   fixed-effect design `X`, the random-effect design `Z` carrying the
#'   row names of the input data, the grouping-factor structure, and the
#'   parameter vectors, using lme4's vocabulary; `model.matrix()` returns
#'   the design and `model.frame()` the stored model frame with its row
#'   names. The frame also keeps the `terms`, `xlevels`, and `contrasts`
#'   of each linear predictor, frozen at fit time and reapplied to
#'   `newdata`, plus the frozen bases of data-dependent terms such as
#'   `poly()` and `scale()`. A name outside the vocabulary errors and
#'   lists the accepted names.
#'
#' @examples
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   set.seed(1)
#'   dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
#'   dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
#'   fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#'   # the designs, for downstream code written against merMod objects
#'   dim(lme4::getME(fit, "X"))
#'   dim(lme4::getME(fit, "Zt"))
#'
#'   # a vector of names returns a named list
#'   str(lme4::getME(fit, c("n_rtrms", "n_rfacs", "sigma")))
#'
#'   # the conditional modes in coefficient space, aligned with Z
#'   head(lme4::getME(fit, "b"))
#'   # note: "lower" is all -Inf here, because the internal covariance
#'   # parameterization is unbounded. Use diagnose() to spot a singular
#'   # fit, not theta == lower.
#'   lme4::getME(fit, "lower")
#' }
#' @exportS3Method lme4::getME
getME.frmtmb_fit <- function(object, name, resp = NULL, ...) {
  if (missing(name) || !is.character(name) || !length(name)) {
    stop("getME() needs one or more names: ",
         paste(frmtmb_getME_vocab, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(name, frmtmb_getME_vocab)
  if (length(bad)) {
    stop("getME(): unknown name(s) ", paste(bad, collapse = ", "),
         ". Supported: ", paste(frmtmb_getME_vocab, collapse = ", "),
         call. = FALSE)
  }
  if (length(name) > 1L) {
    out <- lapply(name, function(nm) getME.frmtmb_fit(object, nm, resp))
    return(stats::setNames(out, name))
  }
  # find_linpred() raises the "disambiguate with resp =" error for a
  # multivariate fit, which is exactly the guard the designs need
  mu_lp <- function() find_linpred(object, resp, "mu")
  mu_Z <- function() {
    Z <- mu_lp()$Z
    if (is.null(Z)) {
      Z <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                x = numeric(0),
                                dims = c(object$frame$n_obs,
                                         object$frame$n_c %||% 0L))
    }
    dimnames(Z) <- list(rownames(object$frame$data_frame),
                        re_coef_labels(object$frame))
    Z
  }
  theta <- object$estimates$theta %||% numeric(0)
  names(theta) <- if (length(theta)) {
    paste0("theta_", seq_along(theta))
  }
  switch(
    name,
    X = mu_lp()$X,
    Z = mu_Z(),
    Zt = Matrix::t(mu_Z()),
    beta = ,
    fixef = object$estimates$beta,
    b = {
      bv <- coef_b(object) %||% numeric(0)
      stats::setNames(as.numeric(bv), re_coef_labels(object$frame))
    },
    theta = theta,
    lower = stats::setNames(rep(-Inf, length(theta)), names(theta)),
    sigma = sigma(object),
    flist = getME_flist(object),
    n_rtrms = length(getME_group_blocks(object)),
    n_rfacs = length(getME_flist(object))
  )
}
