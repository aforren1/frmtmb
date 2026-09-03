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
  check_positive(tol, "tol")
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
  # a length-2 location built a prior that recycled against the whole
  # parameter block it was attached to, and `stopifnot(scale > 0)`
  # reported a negative scale as "scale > 0 is not TRUE", which names
  # the test rather than the argument
  check_number(location, "location")
  check_positive(scale, "scale")
  structure(list(kind = "normal", location = location, scale = scale),
            class = "frmtmb_prior")
}

#' @rdname frmtmb-priors
#' @export
prior_t <- function(df = 3, location = 0, scale = 1) {
  check_positive(df, "df")
  check_number(location, "location")
  check_positive(scale, "scale")
  structure(list(kind = "t", df = df, location = location, scale = scale),
            class = "frmtmb_prior")
}

#' Resolve a named prior list to per-component index/parameter vectors.
#' Names may be individual parameters (as in outer_par_names()) or whole
#' components ("beta", "betad", "theta", "thetar", "thetaac").
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
    if (identical(pr$kind, "lkj")) {
      # this spelling addresses parameters one at a time; the LKJ
      # density is over a block's whole correlation and needs the
      # structure's map, which only the class spelling carries
      stop("prior_lkj() addresses a block's whole correlation, so it ",
           "cannot be given by parameter name; write ",
           "set_prior(\"lkj(", format(pr$eta), ")\", class = \"cor\")",
           call. = FALSE)
    }
    if (nm %in% names(comp_names)) {
      idx <- which(!is.na(comp_names[[nm]]))
      add(nm, idx, pr)
      next
    }
    hit <- FALSE
    for (cp in names(comp_names)) {
      # both spellings: the template's own `(Intercept)` and the
      # parenthesis-free `Intercept` the draws, variables() and
      # hypothesis() all use. Priors are written against names the user
      # read off one of those surfaces
      i <- which(comp_names[[cp]] == nm |
                   par_name_bare(comp_names[[cp]]) == par_name_bare(nm))
      if (length(i)) {
        add(cp, i, pr)
        hit <- TRUE
        break
      }
    }
    if (!hit) {
      stop("Unknown parameter in priors: '", nm, "'. Available: ",
           paste(par_name_bare(unlist(comp_names))[
             !is.na(unlist(comp_names))], collapse = ", "),
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
    # the paren-tolerant addressing of confint(parm =), so a name copied
    # out of a hypothesis() expression works here too, plus the bare
    # name of an intercept-only nonlinear parameter: bounds on an ODE
    # model are written against the parameters of the dynamics (la),
    # not against their design-matrix spelling (la_(Intercept))
    pos <- apply_nlpar_alias(fit, names(x), match_par_name(names(x), nm))
    if (anyNA(pos)) {
      stop("Unknown parameter(s) in bounds: ",
           paste(names(x)[is.na(pos)], collapse = ", "), ". Available: ",
           paste(nm, collapse = ", "),
           " (parentheses may be dropped, and intercept-only nonlinear ",
           "parameters may be named bare)", call. = FALSE)
    }
    out[pos] <- as.numeric(x)
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

# --- non-centered sampling --------------------------------------------
#
# See R/covstruct.R for the factor accessors and for which structures
# have one. Everything here is SAMPLING-ONLY: the fit, its tape and every
# number it reports are untouched, because the Laplace approximation
# integrates `b` out and that integral does not care how the integrated
# variable is parameterized.

#' Why a structure has no non-centered form. Read out loud by the
#' message that names the blocks left centered, so a fallback is never
#' silent.
#'
#' @noRd
ncp_reasons <- c(
  us_t = "Student-t latent: a scale mixture, not a linear factor",
  diag_t = "Student-t latent: a scale mixture, not a linear factor",
  car = "sparse CAR precision; its factor is dense",
  spde = "sparse SPDE precision; its factor is dense",
  gr_prec = "sparse precision from gr(prec = ); its factor is dense",
  gp = "dense kernel over the observed positions",
  ou = "dense covariance over the field's positions",
  exp = "dense covariance over the field's positions",
  gau = "dense covariance over the field's positions",
  mat = "dense covariance over the field's positions",
  toep = "not positive definite everywhere in its parameterization",
  homtoep = "not positive definite everywhere in its parameterization"
)

#' Why this block stays centered.
#'
#' @noRd
ncp_reason <- function(bk) {
  cs <- bk[["covstruct"]]
  if (is_student_block(bk)) return(unname(ncp_reasons[["us_t"]]))
  # [[ ]] on a named character vector ERRORS on a missing name, so the
  # membership test comes first
  if (cs %in% names(ncp_reasons)) return(unname(ncp_reasons[[cs]]))
  reg <- covstruct_registry[[cs]]
  if (is.null(reg[["chol_sd"]]) && is.null(reg[["chol_L"]])) {
    return("no Cholesky factor is registered for this structure")
  }
  # a factor exists, but some parameter of the block is neither a
  # standard deviation nor a correlation the prior lane can name, so
  # nothing would bound the direction the transform opens
  paste("it has a parameter that is neither a standard deviation nor a",
        "correlation with a density on it, so a non-centered chain",
        "would be free in a direction no prior bounds")
}

#' The `theta` positions a resolved prior specification covers.
#'
#' @noRd
ncp_priored_theta <- function(entries) {
  unlist(lapply(entries %||% list(), function(e) {
    if (identical(e$comp, "theta")) e$idx
  })) %||% integer(0)
}

#' Which blocks a sampling run non-centers, and which it leaves alone.
#'
#' TWO conditions, and they are the same condition twice. A block is
#' non-centered only when (a) every parameter it has is a standard
#' deviation or a correlation, with a registered factor (see
#' R/covstruct.R), and (b) every one of those parameters carries a
#' PRIOR. Both are about not handing the sampler a direction it can run
#' away in.
#'
#' Non-centering gives NUTS the run of theta's whole range. Under a FLAT
#' prior on a log standard deviation that range has a flat tail at the
#' bottom: send `sd` down and `z` stays where it is, `b = sd z` goes to
#' zero, the likelihood settles on the no-random-effect model and the
#' density stops changing. The centered funnel is the only thing that
#' keeps a chain out of it. Measured on a six-group random-intercept
#' model, one chain of 2000, three seeds: sampled from the FIT, where
#' every prior is flat, the non-centered chain walks `theta` to -1e15
#' with a bulk-ESS of 1; sampled from the FORMULA, where the default
#' `student_t(3, 0, s)` on the sd makes that tail integrable, it runs at
#' a bulk-ESS of 174-284 against 3-48 centered. Adding the prior by hand
#' on the fit route does the same thing there (297-355).
#'
#' A CORRELATED block is the same argument with a different parameter.
#' Flat on frmtmb's correlation parameter is `(1 - rho^2)^-3/2`, which
#' is improper, and before 0.39 that kept every correlated block
#' centered whatever else was priored. The LKJ prior (R/priors.R) closes
#' that tail, the formula route applies `lkj(1)` by default, and the
#' rule above then admits those blocks with nothing else changed.
#'
#' `laplace = TRUE` integrates the random effects out, so there is no
#' funnel left to fix and the transform would only cost inner Newton
#' steps: that route stays centered whatever was asked for.
#'
#' @noRd
ncp_plan <- function(fit, reparameterize, laplace, entries = NULL) {
  blocks <- fit$frame$re_blocks %||% list()
  if (!isTRUE(reparameterize) || !length(blocks) || laplace) {
    return(list(idx = integer(0), centered = character(0),
                labels = character(0)))
  }
  priored <- ncp_priored_theta(entries)
  has_prior <- vapply(blocks, function(bk) {
    all(bk[["theta_idx"]] %in% priored)
  }, TRUE)
  ok <- vapply(blocks, ncp_eligible, TRUE) & has_prior
  # an rr() block's own coefficients ARE standard normal factors, so it
  # is non-centered already and has nothing to report
  done <- vapply(blocks, function(bk) {
    identical(bk[["covstruct"]], "rr")
  }, TRUE)
  lab <- vapply(blocks, function(bk) {
    paste0(bk[["term_label"]], " [", bk[["covstruct"]], "]")
  }, "")
  left <- which(!ok & !done)
  list(idx = which(ok), labels = unname(lab[ok]),
       centered = unname(vapply(left, function(i) {
         bk <- blocks[[i]]
         why <- if (!ncp_eligible(bk)) {
           ncp_reason(bk)
         } else if (block_n_cor(bk) > 0L) {
           # the correlation is the interesting half here: flat on
           # frmtmb's correlation parameter is (1 - rho^2)^-3/2, which
           # is improper, and the centered geometry was what kept a
           # chain out of that tail
           paste("its variance or correlation parameters have a flat",
                 "prior here, and flat on a correlation is",
                 "(1 - rho^2)^-3/2, improper. Give the block priors,",
                 "set_prior(class = \"sd\") and",
                 "set_prior(class = \"cor\"), which the formula",
                 "interface supplies for you")
         } else if (identical(bk[["covstruct"]], "hsgp")) {
           # the class-"sd" default cannot reach a lengthscale, so the
           # generic flat-prior advice below would never fix this block
           paste("its lengthscales share the block's theta and have a",
                 "flat prior here, which the default priors do not",
                 "cover. Prior the whole block,",
                 "set_prior(class = \"theta\"), to non-center it")
         } else {
           paste("its variance parameter has a flat prior here, and a",
                 "non-centered chain walks the flat tail that opens at",
                 "sd = 0. Give it a prior, set_prior(class = \"sd\"),",
                 "which the formula interface supplies for you")
         }
         paste0(lab[i], ": ", why)
       }, "")))
}

#' One line per block that could not be non-centered, so a partial
#' fallback is visible without being an error. `message()`, not
#' `warning()`: the run is correct either way, only slower to mix.
#'
#' @noRd
announce_ncp <- function(plan) {
  if (!length(plan$centered)) return(invisible(NULL))
  hdr <- if (!length(plan$idx)) {
    paste("frm_sample(): sampling stays centered: no random-effect",
          "block of this model has a non-centered form:")
  } else {
    paste("frm_sample(): non-centered where possible; these blocks",
          "stay centered:")
  }
  message(paste(c(hdr, paste0("  ", plan$centered)), collapse = "\n"))
  invisible(NULL)
}

#' The sampling objective with the planned blocks reparameterized, and
#' the prior terms (if any) added the way prior_augmented_obj() adds
#' them. The frame copy carrying `ncp_blocks` never leaves this
#' function, so no fitted object can acquire a non-centered tape.
#'
#' @noRd
ncp_objective <- function(fit, idx, entries) {
  frame <- fit$frame
  frame$ncp_blocks <- idx
  nll <- build_objective(frame)
  fn <- if (length(entries)) {
    nlp <- neg_log_prior_fn(entries)
    function(pars) nll(pars) + nlp(pars)
  } else {
    nll
  }
  tpl <- fit$frame$par_template
  # [[ ]] to avoid $ partial matching ("b" matching "beta" in GLMs)
  random <- c(if (!is.null(tpl[["b"]])) "b",
              if (!is.null(tpl[["miss"]])) "miss")
  if (!length(random)) random <- NULL
  RTMB::MakeADFun(fn, ncp_start_pars(fit, idx), random = random,
                  map = fit$frame$map, silent = TRUE)
}

#' The starting parameter list on the non-centered scale: `z0 = L^-1 b`
#' per reparameterized block, so that `init = "last.par.best"` still
#' starts the chains at the ML mode of the MODEL rather than at the
#' mode of one parameterization and the origin of the other.
#'
#' @noRd
ncp_start_pars <- function(fit, idx) {
  pars <- fit$estimates
  th <- pars$theta %||% numeric(0)
  bv <- pars[["b"]]
  for (i in idx) {
    bk <- fit$frame$re_blocks[[i]]
    bv[bk$b_idx] <- ncp_unscale_b(bk, bv[bk$b_idx], th[bk$theta_idx])
  }
  pars[["b"]] <- bv
  pars
}

#' Put the draws matrix back on the centered scale, draw by draw: each
#' `z` maps through the `L(theta)` of ITS OWN draw. The `b[i]` columns
#' that come out are the same quantity, in the same order and under the
#' same names, as the centered route's, which is what lets every draws
#' method downstream stay ignorant of the difference.
#'
#' @noRd
ncp_backtransform <- function(m, fit, idx) {
  pidx <- draws_par_index(fit)
  bc <- pidx[["b"]]
  tc <- pidx[["theta"]]
  blocks <- fit$frame$re_blocks
  for (r in seq_len(nrow(m))) {
    th <- unname(m[r, tc])
    for (i in idx) {
      bk <- blocks[[i]]
      cols <- bc[bk$b_idx]
      m[r, cols] <- ncp_scale_b(bk, unname(m[r, cols]),
                                th[bk$theta_idx])
    }
  }
  m
}

#' Labels for the full sampled parameter vector, in template order,
#' skipping mapped entries; b kept as `b[i]`. `include_random = FALSE`
#' drops the inner components (b, miss) for laplace-marginalized draws.
#'
#' Parentheses are dropped: `Intercept`, not `(Intercept)`. A draws
#' object is the brms-facing surface, and brms names draws
#' `b_Intercept` / `sd_g__Intercept` with no parentheses anywhere. These
#' names go straight into posterior and bayesplot, where a
#' parenthesized column has to be backquoted in every expression that
#' touches it, and they are already the vocabulary `variables()` and
#' `hypothesis()` speak. The FIT side keeps its own canonical spelling
#' (`confint()` and `vcov()` still show `(Intercept)`), and every
#' draws-side lookup accepts both (`match_par_name()`).
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
  par_name_bare(out)
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

#' Assemble the unfitted object the formula interface samples.
#'
#' @noRd
sample_assemble <- function(formula, data, family, data2, start,
                            control, na.action, REML, lower, upper) {
  if (!inherits(formula, c("formula", "frmtmb_formula",
                           "frmtmb_mvformula"))) {
    stop("frm_sample() takes a frmtmb fit or a formula; got an object ",
         "of class ", paste(class(formula), collapse = "/"),
         ". Fit with frm() first, or pass bf(y ~ x) with data =",
         call. = FALSE)
  }
  if (is.null(data)) {
    stop("frm_sample() from a formula needs data =: there is no fitted ",
         "model to take the design from", call. = FALSE)
  }
  frm(formula, data, family = family, REML = REML, start = start,
      control = control, na.action = na.action, lower = lower,
      upper = upper, data2 = data2, dry_run = "objective")
}

# --- default priors for the formula interface --------------------------
#
# WHY ONLY HERE. frm_sample(fit) is a DIAGNOSTIC: check_laplace()
# compares the Laplace/Wald approximation against the shape of the
# LIKELIHOOD's own posterior, and a default prior would change the thing
# being measured. frm_sample(formula) is not a diagnostic - the
# posterior is the whole answer - so it defaults to brms's own
# weakly-informative priors instead of to flat improper ones.
#
# THE RECIPE is brms 2.23's, read off `brms::default_prior()` on matched
# models (tests/testthat/test-sample-direct.R checks it against brms
# directly where brms is installed). Writing y* for the response
# transformed by the mu link:
#
#   b (slopes)      flat, as brms leaves them
#   Intercept       student_t(3, round(median(y*), 1), s)
#   sd              student_t(3, 0, s) on the natural sd scale
#   cor             lkj(1), uniform over correlation matrices, carried
#                   onto frmtmb's own correlation parameters with the
#                   Jacobian of that map (see R/priors.R)
#   sigma           student_t(3, 0, s) on the natural scale, when sigma
#                   is intercept-only; student_t(3, 0, 2.5) on the LOG
#                   scale when sigma carries its own predictor (brms
#                   makes the same distinction)
#   where s = max(2.5, round(mad(y*), 1))
#
# The transform is `brms:::def_scale_prior.brmsterms()`, followed line
# for line rather than approximated:
#
#   - the transformable links are identity, log, inverse, sqrt and
#     1/mu^2, and a LOG-SCALE family (lognormal and its relatives, whose
#     mu link is spelled identity but whose response is logged inside
#     the density) counts as log whatever its link says. Every other
#     link - logit, cloglog, ... - keeps location 0 and scale 2.5, which
#     is what brms reports for a bernoulli model;
#   - a response that was a FACTOR, or whose family spreads its dpars
#     over categories (ordinal, categorical, multinomial), is not
#     transformed either, so those keep 0 and 2.5. brms spells that
#     `!is_like_factor(y) && !conv_cats_dpars(x)`;
#   - the zero shift is PER ELEMENT and applies to zeros only -
#     `ifelse(y == 0, y + 0.1, y)` - and only under log/inverse/1/mu^2,
#     not under sqrt. Shifting the whole vector instead moves the
#     median: a count response with median 1 and some zeros gets
#     location log(1) = 0 from brms and log(1.1) = 0.1 from the naive
#     spelling;
#   - `round(mad(y*), 1)` raises the scale only when it is finite, and
#     `round(median(y*), 1)` sets the location only when it is finite.
#
# WHAT IS NOT MATCHED, and why, is documented in ?frm_sample and said
# out loud by the disclosure message: a correlation whose STRUCTURE has
# no LKJ density (toep(), which is not positive definite everywhere, so
# there is no correlation matrix over its parameter space to put one
# on), ORDINAL THRESHOLDS (brms
# priors them as its Intercept class; frmtmb keeps them in the `thres`
# extra-parameter vector, which set_prior() cannot address), the
# shape/phi/nu-style dispersion parameters (brms uses gamma and
# inverse-gamma densities that set_prior() does not carry), and
# multivariate models (set_prior() cannot address one response of
# several).

# Families whose response is on the log scale inside the density even
# though the mu link is spelled identity; brms's `has_logscale()`.
logscale_families <- c("lognormal", "shifted_lognormal",
                       "hurdle_lognormal")

#' brms's `def_scale_prior()` location and scale for one response.
#'
#' @noRd
default_prior_scale <- function(fit) {
  rspec <- fit$spec$responses[[1L]]
  fam <- rspec$family
  y <- fit$frame$y[[rspec$resp_name]]
  link <- if (fam$family %in% logscale_families) "log" else {
    fam$links[["mu"]]$name %||% "identity"
  }
  flat <- list(location = 0, scale = 2.5, link = link, centered = FALSE)
  if (is.null(y) || !is.numeric(y) || is.matrix(y) || !length(y)) {
    return(flat)
  }
  # a factor response (its codes are labels, not numbers) and a family
  # whose dpars run over categories both skip the transform
  if (!is.null(fit$frame$y_levels[[rspec$resp_name]])) return(flat)
  if (fam$type %in% c("ordinal", "categorical")) return(flat)
  if (!link %in% c("identity", "log", "inverse", "sqrt", "1/mu^2")) {
    return(flat)
  }
  if (link %in% c("log", "inverse", "1/mu^2")) {
    y <- ifelse(y == 0, y + 0.1, y)
  }
  yl <- switch(link, identity = y, log = log(y), inverse = 1 / y,
               sqrt = sqrt(y), `1/mu^2` = 1 / y^2)
  loc <- round(stats::median(yl), 1L)
  scl <- round(stats::mad(yl), 1L)
  list(location = if (is.finite(loc)) loc else 0,
       scale = if (is.finite(scl)) max(2.5, scl) else 2.5,
       link = link, centered = TRUE)
}

#' A `set_prior()` spec whose distribution sits on `exp(coefficient)`
#' rather than on the coefficient, with the class `"sd"` log-Jacobian.
#'
#' @noRd
natural_dpar_prior <- function(dist, dpar) {
  spec <- unclass(set_prior(dist, class = "Intercept", dpar = dpar))
  spec[[1L]]$natural <- TRUE
  structure(spec, class = "frmtmb_priorlist")
}

#' brms's default priors for the model this object holds, as a
#' `frmtmb_priorlist`, or `NULL` when the model has no slot they cover.
#'
#' @noRd
default_priors_for <- function(fit) {
  if (length(fit$spec$responses) > 1L) return(NULL)
  rspec <- fit$spec$responses[[1L]]
  ps <- default_prior_scale(fit)
  pl <- NULL
  add <- function(p) pl <<- if (is.null(pl)) p else pl + p
  st <- function(loc, scl) {
    sprintf("student_t(3, %s, %s)", format(loc), format(scl))
  }

  for (lp in fit$frame$linpreds) {
    if (!is.null(lp$constant) || !is.null(lp$nl_body)) next
    if (!"(Intercept)" %in% colnames(lp$X)) next
    if (lp$dpar %in% rspec$primary_dpars) {
      add(set_prior(st(ps$location, ps$scale), class = "Intercept"))
    } else if (identical(lp$dpar, "sigma") &&
                 identical(lp$link$name, "log")) {
      # brms scales the prior on sigma itself by the response's mad
      # only when sigma is a single number; with a predictor the
      # intercept gets the plain student_t(3, 0, 2.5) on the log scale
      if (ncol(lp$X) == 1L && is.null(lp$Z)) {
        add(natural_dpar_prior(st(0, ps$scale), "sigma"))
      } else {
        add(set_prior(st(0, 2.5), class = "Intercept", dpar = "sigma"))
      }
    }
  }

  has_sd <- any(vapply(fit$frame$re_blocks, function(bk) {
    length(covstruct_registry[[bk$covstruct]]$sd_idx(bk$dim)) > 0L
  }, TRUE))
  if (has_sd) add(set_prior(st(0, ps$scale), class = "sd"))
  # lkj(1), brms's own default: uniform over correlation matrices. Only
  # when some block's correlation HAS an LKJ density, so that a model
  # whose only correlated block is a toep() does not fail on a default
  # it cannot honor (default_prior_notes() names that block instead)
  has_cor <- any(vapply(fit$frame$re_blocks, function(bk) {
    !bk$covstruct %in% names(lkj_refusals) && !is.null(block_cor_spec(bk))
  }, TRUE))
  if (has_cor) add(set_prior("lkj(1)", class = "cor"))
  pl
}

#' What this model has that brms would prior and frmtmb deliberately
#' leaves flat. Every such slot is named out loud in the disclosure
#' message, so a family that gets no defaults never passes in silence.
#'
#' @noRd
default_prior_notes <- function(fit) {
  if (length(fit$spec$responses) > 1L) {
    return(paste("multivariate model: no defaults, because set_prior()",
                 "cannot address one response of several"))
  }
  rspec <- fit$spec$responses[[1L]]
  notes <- character(0)
  if (identical(rspec$family$type, "ordinal")) {
    notes <- c(notes, paste("no defaults for this family's thresholds",
                            "(brms priors them as its Intercept class)"))
  }
  # dispersion dpars brms gives a gamma or inverse-gamma default, which
  # set_prior() cannot express
  disp <- setdiff(intersect(names(rspec$dpars %||% list()),
                            c("shape", "phi", "nu")),
                  rspec$primary_dpars)
  if (length(disp)) {
    notes <- c(notes, paste0("no defaults for ",
                             paste(disp, collapse = ", "),
                             " (brms uses gamma / inverse-gamma)"))
  }
  # correlations get lkj(1), brms's own default, wherever a density
  # exists for the block's parameterization; what is left is named
  ungated <- unique(vapply(fit$frame$re_blocks, function(bk) {
    if (bk$covstruct %in% names(lkj_refusals)) {
      paste0(bk$term_label, " [", bk$covstruct, "]")
    } else {
      ""
    }
  }, ""))
  ungated <- ungated[nzchar(ungated)]
  if (length(ungated)) {
    notes <- c(notes, paste0("no correlation defaults for ",
                             paste(ungated, collapse = ", "),
                             " (no LKJ density fits its parameters)"))
  }
  notes
}

#' The classes a user's priorlist speaks for; a default of the same
#' class steps aside for it, which is brms's partial-override rule.
#'
#' @noRd
priorlist_classes <- function(pl) {
  unique(vapply(unclass(pl), function(s) {
    paste0(s$class, if (nzchar(s$dpar)) paste0(":", s$dpar))
  }, ""))
}

#' Drop the default specs a user specification has taken over.
#'
#' @noRd
drop_overridden <- function(defaults, user_classes) {
  keep <- Filter(function(s) {
    key <- paste0(s$class, if (nzchar(s$dpar)) paste0(":", s$dpar))
    !(key %in% user_classes)
  }, unclass(defaults))
  if (!length(keep)) NULL else structure(keep,
                                         class = "frmtmb_priorlist")
}

#' One compact line per prior class, so the call discloses what it
#' chose - and one line per slot it deliberately left flat, so a model
#' that gets few defaults, or none, says so instead of passing in
#' silence. `message()` rather than `warning()`: this is information,
#' not a fault, and `suppressMessages()` turns it off.
#'
#' @noRd
announce_default_priors <- function(pl, notes) {
  msg <- paste0("frm_sample(): default priors (brms 2.23 defaults; ",
                "priors = \"flat\" opts out)")
  for (s in unclass(pl %||% list())) {
    kind <- if (identical(s$dist$kind, "t")) "student_t" else s$dist$kind
    d <- paste0(kind, "(", paste(unlist(s$dist[-1L]), collapse = ", "),
                ")")
    lab <- if (nzchar(s$dpar)) paste0(s$class, " (", s$dpar, ")") else
      s$class
    msg <- c(msg, sprintf("  %-18s %s%s", lab, d,
                          if (isTRUE(s$natural)) "  [natural scale]"
                          else if (identical(s$class, "sd"))
                            "  [natural sd scale]"
                          else if (identical(s$class, "cor"))
                            "  [correlation matrix]" else ""))
  }
  msg <- c(msg, "  b                  (flat), as brms leaves slopes")
  for (nt in notes) {
    msg <- c(msg, paste0("  ", nt, " - see ?frm_sample"))
  }
  message(paste(msg, collapse = "\n"))
  invisible(NULL)
}

#' Merge resolved prior inputs; a later one wins per parameter. Used to
#' let a user's `priors =` take over individual parameters from the
#' defaults without discarding the rest.
#'
#' @noRd
merge_prior_inputs <- function(base, over) {
  if (is.null(base)) return(over)
  if (is.null(over)) return(base)
  # collapsed, because a class "cor" entry covers a whole correlation
  # segment and its key names every position it holds
  key <- function(e) paste0(e$comp, ".", paste(e$idx, collapse = ","))
  ent <- base$entries
  names(ent) <- vapply(ent, key, "")
  for (e in over$entries) {
    # a base entry that shares ANY position with the incoming one is
    # retired, not just the one with the same key: the default LKJ term
    # spans a whole correlation segment, and leaving it in place beside
    # a per-parameter override would apply both densities there
    hit <- vapply(ent, function(b) {
      identical(b$comp, e$comp) && length(intersect(b$idx, e$idx)) > 0L
    }, TRUE)
    if (any(hit)) ent <- ent[!hit]
    ent[[key(e)]] <- e
  }
  list(entries = unname(ent),
       lower = utils::modifyList(as.list(base$lower),
                                 as.list(over$lower)),
       upper = utils::modifyList(as.list(base$upper),
                                 as.list(over$upper)))
}

#' Resolve the priors of a formula-interface call: brms defaults, the
#' opt-out, and a user specification composed on top.
#'
#' @noRd
sample_resolve_priors <- function(fit, priors, announce = TRUE) {
  if (is.character(priors)) {
    if (!identical(priors, "flat")) {
      stop("priors = must be a set_prior() specification, a named list ",
           "of prior objects, or the string \"flat\" to sample the ",
           "likelihood with improper flat priors; got \"",
           paste(priors, collapse = "\", \""), "\"", call. = FALSE)
    }
    if (length(fit$frame$re_blocks)) {
      warning("priors = \"flat\": every variance component has a flat ",
              "prior on its log standard deviation, under which the ",
              "posterior need not be proper (it usually is not with ",
              "few groups). The chains still run and Rhat cannot see ",
              "it. Drop priors = to get the brms default priors ",
              "instead", call. = FALSE)
    }
    return(list(effective = NULL, ri = NULL))
  }
  defaults <- default_priors_for(fit)
  user_pl <- if (inherits(priors, "frmtmb_priorlist")) priors
  if (!is.null(defaults) && !is.null(user_pl)) {
    defaults <- drop_overridden(defaults, priorlist_classes(user_pl))
  }
  # announced even when there are no defaults at all: a family whose
  # slots frmtmb cannot prior must say so rather than look flat by
  # accident
  if (announce) {
    announce_default_priors(defaults, default_prior_notes(fit))
  }
  eff <- if (is.null(defaults)) user_pl else if (is.null(user_pl)) {
    defaults
  } else {
    defaults + user_pl
  }
  ri <- if (!is.null(eff)) resolve_prior_input(fit, eff)
  if (!is.null(priors) && is.null(user_pl)) {
    # the legacy named-list spelling addresses raw internal parameters;
    # it takes over exactly those and leaves the rest of the defaults
    ri <- merge_prior_inputs(ri, resolve_prior_input(fit, priors))
    eff <- eff %||% structure(list(), class = "frmtmb_priorlist")
    attr(eff, "overrides") <- priors
  }
  list(effective = eff, ri = ri)
}

#' Sample a model with NUTS
#'
#' Runs [tmbstan::tmbstan()] on the model's objective and returns the
#' draws with frmtmb parameter names. Given a [frm()] fit it samples the
#' fitted objective, initialized at the ML estimates, which shortens
#' warmup considerably. Given a formula and `data` it assembles the same
#' objective without optimizing anything first (see Sampling from a
#' formula).
#'
#' **The two routes answer different questions.** `frm_sample(fit)` is a
#' DIAGNOSTIC: it explores the LIKELIHOOD, with flat improper priors on
#' the outer parameters, so that [check_laplace()] can compare the
#' Laplace and Wald approximations against the shape of the very same
#' objective the fit maximized. `frm_sample(formula, data)` is a
#' SAMPLING tool: it samples a POSTERIOR, and defaults to brms's own
#' weakly-informative priors (see Default priors), because from a
#' formula there is no estimate for the run to be a diagnostic for.
#'
#' On the fit path, then, a parameter without a prior keeps a flat
#' improper one, and the posterior of a variance component with few
#' groups can be improper. That is the price of measuring the
#' likelihood, and it is why the fit path is a diagnostic.
#'
#' @section Sampling from a formula:
#' `frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian())`
#' parses, assembles and tapes exactly as [frm()] does, stops before the
#' optimizer, and hands the objective to Stan. Every pre-optimizer
#' refusal still applies (REML, quadrature, the mixture and [hmm()]
#' guards). There is no mode, so the default `init` is `"random"`:
#' Stan's own overdispersed initialization on the unconstrained scale,
#' inside any `lower`/`upper` bounds.
#'
#' The returned object supports the whole draws surface -
#' [summary()], [fixef()], [VarCorr()], [ranef()], [hypothesis()],
#' [posterior_epred()], [posterior_predict()], [posterior_linpred()],
#' [pp_check()], [as_draws()] - because those read the model frame and
#' one draw at a time. The methods that report a maximum-likelihood
#' quantity refuse instead of inventing one: [check_laplace()] (it
#' compares NUTS against a mode that does not exist here), and on the
#' embedded object reachable as `x$fit`, [summary()], [vcov()],
#' [confint()], [logLik()], [fixef()], [ranef()], [VarCorr()],
#' [predict()], [fitted()], [residuals()] and [simulate()].
#'
#' @section Reparameterization:
#' A random-effect block has a FUNNEL in its centered joint posterior:
#' the width of the prior on `b` is a standard deviation that is being
#' sampled at the same time, so the region NUTS must explore narrows as
#' that standard deviation shrinks, and one step size cannot fit both
#' ends of it.
#'
#' `reparameterize = TRUE` (the default) samples `z ~ N(0, I)` instead
#' and computes `b = L(theta) z` on the tape, with `L` the block's own
#' Cholesky factor, which is brms's construction. Each draw is mapped back
#' through ITS OWN `theta`, so the `b[i]` columns of the draws matrix
#' hold the same quantity in the same order under the same names as
#' `reparameterize = FALSE` gives, and every method downstream
#' ([posterior_epred()], [ranef()], [log_lik()], [loo()],
#' [conditional_effects()], [hypothesis()]) reads them without knowing
#' which route produced them. Only the `stanfit` inside the object
#' carries `z`.
#'
#' Priors are unaffected. They are declared on the OUTER parameters
#' (`beta`, `theta`, `betad`), and `theta` means the same thing on both
#' routes; nothing prior-able is reparameterized. What priors DO decide
#' is which blocks are eligible, below.
#'
#' *Which blocks, and why not all of them.* A block is non-centered when
#' two things hold, and they are the same thing twice: every parameter it
#' has is a standard deviation or a correlation with a Cholesky factor
#' registered for its structure, and every one of those parameters
#' carries a PRIOR. Both are about not handing the sampler a direction
#' it can run away in.
#' Non-centering gives NUTS the run of `theta`'s whole range, and the
#' centered funnel is what was keeping a chain out of the parts of that
#' range where the posterior is flat or improper. Removing the funnel
#' without closing those off first trades a slow chain for a wrong one.
#'
#' In practice that means the FORMULA interface, whose default priors
#' cover every standard deviation and every correlation, non-centers
#' `(1 | g)` and any block written one term at a time, `diag()` and
#' `homdiag()` blocks, [mgcv::s()] smooths, `equalto()`, `gr(cov = )`,
#' and the CORRELATED blocks `(x | g)`, `cs()`, `ar1()` and `hetar1()`.
#' A `k =` Hilbert-space [gp()] block stays
#' centered on the formula route even though its factor is diagonal:
#' its LENGTHSCALES share the block's `theta`, the default priors cover
#' only standard deviations, and the gate wants every parameter of a
#' block priored. Prior the whole block by hand
#' (`set_prior(class = "theta")`) to non-center it.
#' `rr()` is already non-centered by
#' construction, since its own coefficients are the standard normal
#' factors. A fitted model sampled with `frm_sample(fit)` has flat
#' priors by design (it is a diagnostic; see above), so it stays
#' centered unless you give its variance parameters a prior.
#'
#' The call `message()`s every block it left centered, with the reason.
#' A model made only of those samples centered throughout and says so
#' rather than failing. The reasons:
#'
#' - A FLAT PRIOR on the block's standard deviation. Send `sd` down with
#'   `z` where it is: `b = sd z` goes to zero, the likelihood settles on
#'   the model without that random effect, and the density stops
#'   changing: a flat tail with nothing to stop a chain in it. Measured
#'   on a six-group random-intercept model, one chain of 2000, three
#'   seeds: from the fit, with flat priors, the non-centered chain walks
#'   `theta` to -1e15 at a bulk-ESS of 1; from the formula, where the
#'   default `student_t(3, 0, s)` makes that tail integrable, 174 to 284
#'   against 3 to 48 centered.
#' - A FLAT PRIOR on the block's CORRELATION, which is the same
#'   argument: flat on frmtmb's unbounded correlation parameter is
#'   `(1 - rho^2)^-3/2`, improper, with all its mass at `|rho| = 1`.
#'   Before 0.39 that kept every correlated block centered, because no
#'   proper correlation prior existed to close the tail; `lkj(eta)`
#'   now does, and the formula route sets `lkj(1)` by default, so
#'   `(Days | Subject)`, `cs()`, `ar1()` and `hetar1()` non-center
#'   there like any other block. On the FIT route they still do not:
#'   its priors are flat by design.
#' - A Student-t latent (`gr(dist = "student")`): a scale mixture, not a
#'   linear factor.
#' - [car()], `spde()` and `gr(prec = )`: sparse precisions whose factor
#'   is dense.
#' - The exact [gp()] and the spatial covariances (`ou`, `exp`, `gau`,
#'   `mat`): a full factorization per gradient evaluation.
#' - `toep()`: not positive definite everywhere in its parameterization,
#'   so `b = L z` is not a bijection there.
#'
#' *What it is worth.* One chain of 2000 iterations, three seeds
#' (dev/benchmarks.md). Where each group's own data say little (the
#' regime the funnel lives in) it is decisive: 80 groups of 2 binary
#' observations run at a min-ESS of 236 against 5 centered, 55 times the
#' effective draws per second. Where the groups are informative it is a
#' wash: the `epilepsy` GLMM and an uncorrelated `sleepstudy` are within
#' noise of the centered chain either way. It is not a blanket speed-up,
#' and the blocks it declines to touch are the ones where it would have
#' done harm.
#'
#' `laplace = TRUE` integrates the random effects out, which removes the
#' funnel by itself, so that route ignores `reparameterize` entirely.
#'
#' The ML fit is untouched either way. The Laplace approximation
#' integrates `b` out, and that integral is invariant under a linear
#' change of the integrated variable, so `frm()` has nothing to change
#' and no fitted object ever carries a non-centered tape.
#'
#' @section Default priors:
#' The formula interface defaults to brms 2.23's own weakly-informative
#' priors, read off `brms::default_prior()` on matched models. Write
#' `y*` for the response transformed by the `mu` link and
#' `s = max(2.5, round(mad(y*), 1))`:
#'
#' | class | default | scale |
#' |---|---|---|
#' | `b` (slopes) | flat | - |
#' | `Intercept` | `student_t(3, round(median(y*), 1), s)` | link |
#' | `sd` | `student_t(3, 0, s)` | natural sd, log-Jacobian applied |
#' | `cor` | `lkj(1)` | correlation matrix, Jacobian applied |
#' | `sigma` (intercept only) | `student_t(3, 0, s)` | natural |
#' | `sigma` (with a predictor) | `student_t(3, 0, 2.5)` | log |
#'
#' The link is transformed only for `identity`, `log`, `inverse`,
#' `sqrt` and `1/mu^2` - brms's own list - with a log-scale family
#' ([lognormal()] and its relatives, whose `mu` link is spelled
#' `identity` but whose response is logged inside the density) counted
#' as `log`. Under any other link (`logit`, `cloglog`, ...), for a
#' response that was a factor, and for a family whose parameters run
#' over categories (ordinal, categorical, multinomial), the location is
#' 0 and the scale 2.5 - which is what brms reports for a bernoulli
#' model. Under a `log`, `inverse` or `1/mu^2` link the ZEROS of the
#' response are replaced by `0.1` before the transform, element by
#' element, exactly as brms does; the non-zero values are left alone,
#' so a count response with median 1 keeps `log(1) = 0` as its
#' location.
#'
#' The call `message()`s the priors it chose, one compact line per
#' class, AND one line per slot it deliberately left flat, so a model
#' that gets few defaults - or none - says so rather than looking flat
#' by accident. Wrap it in `suppressMessages()` to silence that.
#' `prior_summary()` on the returned draws reproduces the chosen priors
#' exactly.
#'
#' *What is deliberately NOT matched.* Each of these is named in the
#' message whenever the model has one.
#' - A `toep()` correlation stays flat. Its banded parameterization is
#'   not positive definite everywhere, so there is no correlation matrix
#'   over its whole parameter space for an LKJ density to be about.
#'   Every other correlated structure gets `lkj(1)`, brms's own default
#'   (see [set_prior()] for what the density is on frmtmb's
#'   parameters).
#' - ORDINAL THRESHOLDS stay flat. brms priors them
#'   `student_t(3, 0, 2.5)` under its `Intercept` class; frmtmb keeps
#'   them in the `thres` extra-parameter vector, which is not a design
#'   column and which [set_prior()]'s class vocabulary
#'   (`b`/`Intercept`/`sd`/`theta`) cannot address. An ordinal model
#'   still gets its `sd` defaults, and the message names the gap.
#' - The `shape`, `phi` and `nu` dispersion parameters stay flat: brms
#'   gives them gamma and inverse-gamma defaults, which [set_prior()]
#'   does not carry.
#' - MULTIVARIATE models get no defaults at all, because [set_prior()]
#'   cannot address one response of several.
#'
#' *Overriding and opting out.* A `set_prior()` specification takes over
#' the classes it names and leaves the other defaults in place, which is
#' brms's partial-override rule; a named list of prior objects takes
#' over exactly the internal parameters it names. `priors = "flat"`
#' turns the defaults off entirely and samples the likelihood, which
#' warns when the model has variance components: their flat-prior
#' posteriors need not be proper, and neither the chains nor Rhat can
#' see that.
#'
#' @param fit A `frmtmb_fit`, or a `bf()`/formula to assemble and sample
#'   directly (then `data` is required).
#' @param data A data frame of model data, when `fit` is a formula.
#' @param family Family, when `fit` is a plain formula that does not
#'   carry one (`frm_sample(bf(y ~ x), data = dd, family = poisson())`;
#'   the `+` spelling `bf(y ~ x) + poisson()` works too).
#' @param data2,start,control,na.action,REML As in [frm()]; used only on
#'   the formula path.
#' @param ... Passed to [tmbstan::tmbstan()] (`chains`, `iter`,
#'   `laplace`, `cores`, ...). On Windows more than one core falls back
#'   to sequential chains with a warning: parallel chains run on socket
#'   workers, which can evaluate neither the RTMB tape nor the
#'   objective closure (the known RTMB limitation of tmbstan,
#'   tmbstan#27). The fallback also covers a core count inherited from
#'   `options(mc.cores)`, which is what rstan reads when `cores` is not
#'   given. Fork clusters on unix can, so `cores` works there.
#' @param priors Priors: a [set_prior()] specification, or a named list
#'   of prior objects (see [prior_normal()]) whose names are parameter
#'   names as in the draws (or whole components: `"beta"`, `"theta"`,
#'   ...), or the string `"flat"`. The objective is re-taped with the
#'   prior terms added; a fitted model itself is unchanged. On the fit
#'   path a parameter without a prior keeps a flat improper one. On the
#'   formula path the brms default priors apply to whatever the
#'   specification leaves alone (see Default priors), and
#'   `priors = "flat"` opts out of them entirely.
#' @param lower,upper Optional named numeric vectors of hard bounds on
#'   outer parameters (brms `lb`/`ub`), applied on the internal scale
#'   through Stan's constrained transforms. Chain starting values are
#'   clamped strictly inside the bounds; a bound that excludes the ML
#'   mode itself warns, because the chains then no longer start there.
#'   Names as in `confint()` rows, with parentheses optional; a
#'   nonlinear parameter declared intercept-only (`la ~ 1`) may be named
#'   bare, `lower = c(la = 0)` for `la_(Intercept)`. One that carries
#'   several coefficients is refused rather than resolved to one.
#' @param init Initialization. On a fit the default
#'   (`"last.par.best"`) starts chain 1 exactly at the ML mode and every
#'   further chain at the mode plus a normal perturbation of sd
#'   `init_jitter` on the unconstrained scale. The mode anchor keeps
#'   warmup short; the jitter keeps the chains overdispersed enough for
#'   Rhat to retain power against multimodality (the standard objection
#'   to identical mode starts). `"random"` requests Stan's own
#'   overdispersed initialization, and is the default from a formula,
#'   where there is no mode to start at.
#' @param init_jitter Per-chain perturbation sd for the default init;
#'   `0` starts every chain exactly at the mode. Draws from the R
#'   session's RNG, so `set.seed()` makes the inits reproducible.
#' @param reparameterize Sample the qualifying random-effect blocks in
#'   their non-centered form (`b = L(theta) z`), which removes the
#'   funnel of the centered joint posterior. `TRUE` by default. Not
#'   every block qualifies, and the call says which did not and why;
#'   see Reparameterization, which also explains why nothing downstream
#'   can tell the two routes apart. The mode-anchored `init` still
#'   starts at the ML mode: it is mapped through
#'   `z0 = L(theta_hat)^-1 b_hat`. Seeded draws differ between the two
#'   routes, so `reparameterize = FALSE` is also the way to reproduce a
#'   run made before this default existed.
#' @return An object of class `frmtmb_draws`: list with the `stanfit`,
#'   a draws matrix with named columns (`as.matrix()` method), the
#'   originating fit, and, when any block was non-centered, a
#'   `reparam` note saying which. The `stanfit` holds the parameters as
#'   Stan sampled them, so on a non-centered run its random-effect
#'   columns are `z` while the draws matrix holds `b`.
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
#'
#' # the same model sampled straight from the formula, with no ML fit.
#' # It reports the brms default priors it chose, and prior_summary()
#' # gives them back.
#' ds2 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
#'                   chains = 1, iter = 500, refresh = 0)
#' prior_summary(ds2)
#' fixef(ds2)
#'
#' # a set_prior() specification takes over the classes it names and
#' # leaves the rest of the defaults alone
#' ds3 <- frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian(),
#'                   chains = 1, iter = 500, refresh = 0,
#'                   priors = set_prior("exponential(1)", class = "sd"))
#' prior_summary(ds3)
#' }
#' }
#' @export
frm_sample <- function(fit, data = NULL, family = NULL, ...,
                       priors = NULL, lower = NULL,
                       upper = NULL, init = NULL,
                       init_jitter = 0.25, reparameterize = TRUE,
                       data2 = list(), start = NULL,
                       control = frmtmb_control(),
                       na.action = stats::na.omit, REML = FALSE) {
  if (!requireNamespace("tmbstan", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    stop("frm_sample() needs the 'tmbstan' and 'rstan' packages",
         call. = FALSE)
  }
  if (!is.logical(reparameterize) || length(reparameterize) != 1L ||
      is.na(reparameterize)) {
    stop("reparameterize = must be TRUE or FALSE: it selects the ",
         "non-centered sampling parameterization of the random-effect ",
         "blocks that have one", call. = FALSE)
  }
  from_formula <- !inherits(fit, "frmtmb_fit")
  if (from_formula) {
    fit <- sample_assemble(fit, data, family, data2 = data2,
                           start = start, control = control,
                           na.action = na.action, REML = REML,
                           lower = lower, upper = upper)
  } else if (!is.null(data) || !is.null(family)) {
    stop("frm_sample(data =, family =) belongs to the formula ",
         "interface; the model of a fitted object is already fixed. ",
         "Drop them, or pass the formula instead of the fit",
         call. = FALSE)
  }
  # no mode to anchor on when the model was never optimized
  init <- init %||% if (from_formula) "random" else "last.par.best"
  pr_lower <- c()
  pr_upper <- c()
  obj <- fit$obj
  ri <- NULL
  if (from_formula) {
    # the formula route samples a posterior, so it defaults to brms's
    # weakly-informative priors rather than to the flat ones the
    # diagnostic route needs
    rp <- sample_resolve_priors(fit, priors)
    fit$priors <- rp$effective
    ri <- rp$ri
  } else {
    # a MAP fit carries its priors into sampling unless overridden
    priors <- priors %||% fit$priors
    if (is.character(priors)) {
      stop("priors = \"flat\" is the formula interface's opt-out from ",
           "its default priors. frm_sample() on a fit already samples ",
           "the likelihood with flat priors, so drop the argument",
           call. = FALSE)
    }
    if (!is.null(priors)) ri <- resolve_prior_input(fit, priors)
  }
  if (!is.null(ri)) {
    pr_lower <- ri$lower
    pr_upper <- ri$upper
  }
  # explicit lower/upper override set_prior() bounds on overlap
  lower <- utils::modifyList(as.list(pr_lower),
                             as.list(lower %||% c()))
  upper <- utils::modifyList(as.list(pr_upper),
                             as.list(upper %||% c()))
  laplace <- isTRUE(list(...)$laplace)
  # the priors ride on the OUTER parameters, so they compose with the
  # reparameterization rather than interacting with it: theta is theta on
  # either route, and no prior addresses b (or z) at all
  ncp <- ncp_plan(fit, reparameterize, laplace, ri$entries)
  announce_ncp(ncp)
  if (length(ncp$idx)) {
    obj <- ncp_objective(fit, ncp$idx, ri$entries %||% list())
  } else if (!is.null(ri) && length(ri$entries)) {
    obj <- prior_augmented_obj(fit, ri$entries)
  }
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
  check_tmbstan_build("frm_sample()")
  sf <- do.call(tmbstan::tmbstan, args)
  if (!length(sf@sim) || !length(sf@sim$samples)) {
    stop("frm_sample(): the sampler returned no draws (rstan printed ",
         "the cause above). A known case: a fit whose tape solves ",
         "ODEs (frm_ode) fails inside tmbstan even at the fitted ",
         "optimum, an RTMBode/tmbstan interaction under upstream ",
         "investigation - dev/feature-gaps.md lists it. A solver ",
         "failure mid-run also corrupts rstan's sampler state, so ",
         "retry in a fresh R session", call. = FALSE)
  }
  a <- rstan::extract(sf, permuted = FALSE)   # iter x chain x par
  stan_names <- dimnames(a)[[3]]
  # laplace draws skip the inner components entirely; labeling them
  # with the full template order would misattribute theta as b[i]
  labels <- all_par_labels(fit, include_random = !laplace)
  n_lab <- min(length(labels), length(stan_names))
  stan_names[seq_len(n_lab)] <- labels[seq_len(n_lab)]
  m <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
  colnames(m) <- stan_names
  if (length(ncp$idx)) m <- ncp_backtransform(m, fit, ncp$idx)
  structure(list(stanfit = sf, draws = m, fit = fit,
                 reparam = if (length(ncp$idx)) {
                   list(blocks = ncp$idx, labels = ncp$labels,
                        centered = ncp$centered)
                 }),
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
#' This is a diagnostic tool: it explores the LIKELIHOOD, with flat
#' priors, which is what makes the comparison against the ML mode and
#' its Wald standard errors meaningful. `frm_sample()` on a formula is
#' the sampling tool instead: it samples a POSTERIOR, under brms's
#' default priors. A default prior here would change the very thing
#' being measured, so `check_laplace()` never sets one.
#'
#' That is also why it samples CENTERED. `frm_sample()`'s non-centered
#' parameterization (see Reparameterization there) is offered only for
#' blocks whose variance parameters carry a prior, and here none do: the
#' flat prior that makes the comparison meaningful is exactly the one
#' that leaves a flat tail at `sd = 0` for a non-centered chain to walk
#' into. So the default costs this function nothing and changes nothing
#' about it. Give the variance parameters a prior through `priors =`
#' and the run non-centers; but then it is measuring the Laplace
#' approximation of a different posterior, which is usually not the
#' question.
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
  # the whole point is NUTS against the ML mode and its Wald standard
  # errors; without a mode there is nothing to check the sampler against
  require_fitted(fit, "check_laplace()")
  ds <- frm_sample(fit, chains = chains, iter = iter, ...)
  m <- ds$draws
  keep <- setdiff(colnames(m),
                  c("lp__", grep("^b\\[", colnames(m), value = TRUE)))
  ml <- fit$opt$par
  se <- sqrt(diag(sdr_of(fit)$cov.fixed))
  stopifnot(length(ml) == length(keep))
  post_mean <- colMeans(m[, keep, drop = FALSE])
  post_sd <- apply(m[, keep, drop = FALSE], 2, stats::sd)
  # the check is only as good as the chain: a short chain that wandered
  # inflates post_sd and reads as "Laplace questionable" when the truth
  # is "chain unusable", so the effective sample size rides along and
  # an unhealthy chain is called out as such
  ess <- if (requireNamespace("posterior", quietly = TRUE)) {
    vapply(keep, function(nm) posterior::ess_bulk(m[, nm]), numeric(1))
  } else {
    rep(NA_real_, length(keep))
  }
  out <- data.frame(
    parameter = keep,
    ml = unname(ml),
    post_mean = unname(post_mean),
    wald_se = unname(se),
    post_sd = unname(post_sd),
    z_shift = unname((post_mean - ml) / post_sd),
    sd_ratio = unname(post_sd / se),
    ess_bulk = unname(ess),
    row.names = NULL
  )
  if (any(is.finite(out$ess_bulk) & out$ess_bulk < 100)) {
    message("check_laplace(): the chain mixed too poorly to judge the ",
            "approximation (bulk ESS under 100 for ",
            paste(out$parameter[is.finite(out$ess_bulk) &
                                  out$ess_bulk < 100], collapse = ", "),
            "). Rerun with more iterations before reading z_shift or ",
            "sd_ratio")
  }
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
  check_tmbstan_build("as_tmbstan()")
  tmbstan::tmbstan(fit$obj, ...)
}

#' Refuse a tmbstan build that silently samples the wrong density.
#'
#' tmbstan generates its Stan model at INSTALL time: tools/autogen.R
#' rewrites a `std_normal_lpdf` placeholder in the stanc output into
#' the call that evaluates the TMB objective. stanc 2.39.0 (shipped by
#' StanHeaders 2.39.1, on CRAN 2026-09-02) emits TWO log_prob_impl
#' overloads where 2.32 emitted one, and autogen patches only the
#' first match. HMC reads value and gradient from the unpatched
#' reverse-mode overload, so every chain samples a standard normal
#' instead of the model, with the objective, priors and data all
#' silently absent. Proven by a one-line autogen patch flipping the
#' behavior; dev/prior-dropping-investigation.md has the full record.
#'
#' The check is static (one read of the installed model.hpp, cached
#' per session) and names the exact defect, so a user on an affected
#' build gets a refusal before any sampling rather than plausible
#' garbage after it.
#'
#' @noRd
tmbstan_build_broken <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    hpp <- system.file("model.hpp", package = "tmbstan")
    cached <<- nzchar(hpp) &&
      any(grepl("std_normal_lpdf<propto__>(y)",
                readLines(hpp, warn = FALSE), fixed = TRUE))
    cached
  }
})

check_tmbstan_build <- function(caller) {
  if (tmbstan_build_broken()) {
    stop(caller, ": this tmbstan installation was built against ",
         "StanHeaders >= 2.39, whose code generator leaves one of the ",
         "two generated log-density overloads unpatched, so EVERY ",
         "chain silently samples a standard normal instead of the ",
         "model (tmbstan tools/autogen.R replaces only the first ",
         "match). Until an upstream fix, install a binary tmbstan ",
         "build, or reinstall tmbstan with StanHeaders 2.32.10, and ",
         "distrust any draws already produced by this installation",
         call. = FALSE)
  }
  invisible(NULL)
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
    # a multi-membership row belongs to several levels at once, so there
    # is no per-observation grouping factor to report; lme4's flist has
    # no representation for one, and returning the first member would
    # be a wrong answer rather than a missing one
    if (!is.null(comp$mm)) next
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
