#' The refusal a function whose `...` would otherwise swallow the
#' retired `priors =` spelling raises instead.
#'
#' The argument is `prior`, brms's spelling and the only one. R's
#' partial matching cannot bind `priors` to it (a longer name is not a
#' prefix of a shorter one), so every direct call fails loudly on its
#' own. The exception is a function whose `...` stands ready to absorb
#' the name in silence and fit the model with no priors at all, and
#' frm_sample() is the only one in either package.
#'
#' @noRd
refuse_retired_priors <- function(dots, what) {
  if (!"priors" %in% names(dots)) return(invisible(NULL))
  stop(what, " takes `prior`, not `priors`: the argument follows brms's ",
       "spelling, and this function's `...` would otherwise pass the ",
       "old name through and fit with no priors at all. Rename it to ",
       "`prior`", call. = FALSE)
}

# --- sampler start values and bounds -----------------------------------

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
  if (!covstruct_has_chol(bk)) {
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
                 "set_prior(class = \"cor\"), which frm_sample()",
                 "supplies by default unless prior = \"flat\" turned",
                 "them off")
         } else if (identical(bk[["covstruct"]], "hsgp")) {
           # the class-"sd" default cannot reach a lengthscale, so the
           # generic flat-prior advice below would never fix this block
           paste("its lengthscales share the block's theta and have a",
                 "flat prior here, which no default covers.",
                 "Prior the whole block,",
                 "set_prior(class = \"theta\"), to non-center it")
         } else {
           paste("its variance parameter has a flat prior here, and a",
                 "non-centered chain walks the flat tail that opens at",
                 "sd = 0. Give it a prior, set_prior(class = \"sd\"),",
                 "which frm_sample() supplies by default unless",
                 "prior = \"flat\" turned it off")
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
#' Windows startup note has to read the same source or it never fires.
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

# --- default priors for frm_sample() -----------------------------------
#
# WHERE THEY APPLY. Both of frm_sample()'s routes, because both sample a
# POSTERIOR and the posterior is the whole answer on either. Until
# 0.43.0 the fit route was flat, on the argument that it exists to serve
# check_laplace(); the audit in dev/brms-vignette-audit.md measured what
# that costs a reader who ports a prior-free brm() call onto a fitted
# model, and it is not a diagnostic they get - it is the chain pinned at
# whatever boundary mode maximum likelihood found (kidney's sd(patient)
# at 3e-4), or a chain that mixes an order of magnitude worse than the
# same model sampled from the formula. check_laplace() keeps the
# unpenalized density it needs through frm_sample(.diagnostic = TRUE),
# which is now the only route that is flat by default.
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
# inverse-gamma densities that set_prior() does not carry), NONLINEAR
# parameters (class "b" in brms, which brms leaves flat), and
# multivariate models (the default scales are read off one response).

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

  nlpars <- rspec$nlpars %||% character(0)
  for (lp in fit$frame$linpreds) {
    if (!is.null(lp$constant) || !is.null(lp$nl_body)) next
    if (!"(Intercept)" %in% colnames(lp$X)) next
    # a NONLINEAR parameter's coefficients are class "b" in brms, and
    # brms leaves class "b" flat: the response's median and mad say
    # nothing about a rate or a shape sitting inside a nonlinear body,
    # so a default there would be an invented prior rather than brms's
    if (lp$dpar %in% nlpars) next
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
    length(block_sd_idx(bk)) > 0L
  }, TRUE))
  if (has_sd) add(set_prior(st(0, ps$scale), class = "sd"))
  # lkj(1), brms's own default: uniform over correlation matrices. Only
  # when some block's correlation HAS an LKJ density, so that a model
  # whose only correlated block is a toep() does not fail on a default
  # it cannot honor (default_prior_notes() names that block instead)
  has_cor <- any(vapply(fit$frame$re_blocks, function(bk) {
    identical(block_cor_prior(bk), "lkj")
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
    return(paste("multivariate model: no defaults, because the default",
                 "scales are read off ONE response; write them with",
                 "set_prior(resp = )"))
  }
  rspec <- fit$spec$responses[[1L]]
  notes <- character(0)
  nlp <- rspec$nlpars %||% character(0)
  if (length(nlp)) {
    notes <- c(notes, paste0("nonlinear parameters stay flat: ",
                             paste(nlp, collapse = ", "),
                             " (brms leaves them flat too, as class b; ",
                             "set_prior(nlpar = ) writes them)"))
  }
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
    if (identical(block_cor_prior(bk), "unsupported")) {
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
  unique(vapply(unclass(pl), prior_class_key, ""))
}

#' The class key a default and a user specification are matched on. The
#' nonlinear parameter and the response belong in it for the same
#' reason the distributional parameter does: a prior on ONE of them
#' does not speak for the same class on the others, and a default that
#' stepped aside for it would leave those slots silently flat.
#'
#' @noRd
prior_class_key <- function(s) {
  paste0(s$class, if (nzchar(s$dpar)) paste0(":", s$dpar),
         if (nzchar(s$nlpar %||% "")) paste0("/", s$nlpar),
         if (nzchar(s$resp %||% "")) paste0("@", s$resp))
}

#' Drop the specifications a later, equally explicit one supersedes:
#' same class, coef, group, response, dpar and nlpar, so the same slot.
#' `resolve_priorlist()` would have applied only the later one anyway;
#' this is what keeps `prior_summary()` from listing both, which is how
#' a MAP fit's own prior and a `prior =` on the sampling call would
#' otherwise read.
#'
#' A BOUNDS-ONLY specification supersedes nothing: it is written to
#' tighten an entry an earlier distribution created, and dropping that
#' entry would discard the density it was meant to bound.
#'
#' @noRd
drop_superseded <- function(base, over) {
  tg <- vapply(Filter(function(s) !is.null(s$dist), unclass(over)),
               spec_target, "")
  keep <- Filter(function(s) !(spec_target(s) %in% tg), unclass(base))
  if (!length(keep)) NULL else structure(keep,
                                         class = "frmtmb_priorlist")
}

#' Drop the default specs a user specification has taken over.
#'
#' @noRd
drop_overridden <- function(defaults, user_classes) {
  keep <- Filter(function(s) {
    !(prior_class_key(s) %in% user_classes)
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
                "prior = \"flat\" opts out)")
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
#' let a user's `prior =` take over individual parameters from the
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

#' The slots a prior specification names, for the message that has to
#' say which of the fit's own priors an opt-out discards. Both
#' spellings: a priorlist names classes, the legacy list names internal
#' parameters.
#'
#' @noRd
prior_target_list <- function(pl) {
  tg <- if (inherits(pl, "frmtmb_priorlist")) {
    vapply(unclass(pl), spec_target, "")
  } else {
    names(pl) %||% character(0)
  }
  paste(unique(tg), collapse = "; ")
}

#' Stack prior specifications least explicit first. `resolve_priorlist()`
#' already resolves a later specification over an earlier one, position
#' by position, so the stacking IS the precedence rule and nothing has
#' to be dropped by hand.
#'
#' The legacy named-list spelling addresses raw internal parameters, not
#' classes, so it cannot be concatenated into a priorlist; it rides
#' alongside in its own slot and is merged at the resolved-entry level.
#'
#' @noRd
prior_stack <- function(...) {
  pl <- NULL
  legacy <- list()
  for (p in list(...)) {
    p <- as_priorlist(p)
    if (is.null(p) || !length(p)) next
    if (inherits(p, "frmtmb_priorlist")) {
      if (!is.null(pl)) pl <- drop_superseded(pl, p)
      pl <- if (is.null(pl)) p else pl + p
    } else {
      legacy[[length(legacy) + 1L]] <- p
    }
  }
  list(pl = pl, legacy = legacy)
}

#' Resolve the priors of a sampling call: the brms defaults, the fit's
#' own prior, a specification given on the call, and the opt-out.
#'
#' PRECEDENCE, most explicit wins per addressed slot: `prior =` on the
#' call, then `base` (a fitted object's own prior, which a MAP fit
#' carries into sampling), then the defaults. The first two are stacked
#' in that order and resolved later-wins; the defaults step aside for a
#' whole class either of them names, which is brms's partial-override
#' rule.
#'
#' `defaults = FALSE` is the DIAGNOSTIC mode [check_laplace()] runs in:
#' no defaults, no disclosure, and no opt-out warning. That caller
#' withholds `base` as well, and for the same reason: a MAP fit's
#' penalty is already taped into `fit$obj`, so resolving it here would
#' rebuild a tape the caller already has, and adding anything else
#' would change the density it is measuring.
#'
#' @noRd
sample_resolve_priors <- function(fit, prior, base = NULL,
                                  defaults = TRUE) {
  if (is.character(prior)) {
    if (!identical(prior, "flat")) {
      stop("prior = must be a set_prior() specification, a named list ",
           "of prior objects, or the string \"flat\" to sample the ",
           "likelihood with improper flat priors; got \"",
           paste(prior, collapse = "\", \""), "\"", call. = FALSE)
    }
    if (defaults && length(fit$frame$re_blocks)) {
      warning("prior = \"flat\": every variance component has a flat ",
              "prior on its log standard deviation, under which the ",
              "posterior need not be proper (it usually is not with ",
              "few groups). The chains still run and Rhat cannot see ",
              "it. Drop prior = to get the brms default priors ",
              "instead", call. = FALSE)
    }
    # a MAP fit's penalty is taped INTO fit$obj, so "flat" is a claim
    # about that tape too: the sampled density is the bare likelihood,
    # and the fit's own prior goes with the defaults. Said out loud,
    # because the alternative reading (drop the defaults, keep the
    # fit's penalty) is equally plausible from the argument name
    if (defaults && !is.null(base)) {
      message("prior = \"flat\" drops the prior this fit was made ",
              "with (", prior_target_list(base), ") as well as the ",
              "defaults: the sampled density is the bare likelihood, ",
              "not the penalized objective frm() maximized")
    }
    return(list(effective = NULL, ri = NULL, flat = TRUE))
  }
  st <- prior_stack(base, prior)
  user_pl <- st$pl
  defs <- if (defaults) default_priors_for(fit)
  if (!is.null(defs) && !is.null(user_pl)) {
    defs <- drop_overridden(defs, priorlist_classes(user_pl))
  }
  # announced even when there are no defaults at all: a family whose
  # slots frmtmb cannot prior must say so rather than look flat by
  # accident
  if (defaults) {
    announce_default_priors(defs, default_prior_notes(fit))
  }
  eff <- if (is.null(defs)) user_pl else if (is.null(user_pl)) {
    defs
  } else {
    defs + user_pl
  }
  ri <- if (!is.null(eff)) resolve_prior_input(fit, eff)
  if (length(st$legacy)) {
    # the legacy spelling takes over exactly the internal parameters it
    # names and leaves the rest of the stack in place
    for (lg in st$legacy) {
      ri <- merge_prior_inputs(ri, resolve_prior_input(fit, lg))
    }
    eff <- eff %||% structure(list(), class = "frmtmb_priorlist")
    attr(eff, "overrides") <- Reduce(utils::modifyList, st$legacy)
  }
  list(effective = eff, ri = ri)
}

#' Sample a model with NUTS
#'
#' Runs [tmbstan::tmbstan()] on the model's objective and returns the
#' draws with frmtmb parameter names. Given a [frmtmb::frm()] fit it samples the
#' fitted objective, initialized at the ML estimates, which shortens
#' warmup considerably. Given a formula and `data` it assembles the same
#' objective without optimizing anything first (see Sampling from a
#' formula).
#'
#' **Both routes sample a posterior under the same default priors** (see
#' Default priors), brms 2.23's own weakly-informative ones. What the
#' fit adds is a starting point, not a different density: the ML mode
#' anchors the chains and shortens warmup, and the model is the one the
#' fit already carries.
#'
#' The two differ only in what else is in the stack. A fit made with
#' `frm(prior = )` is a MAP fit, and it carries its own prior into
#' sampling: the priors that apply are then this call's `prior =` first,
#' the fit's own next, and the defaults under both, most explicit
#' winning per slot it addresses. `prior = "flat"` opts out of all of
#' it and samples the bare likelihood, on either route.
#'
#' [check_laplace()] is the one caller that stays unpenalized by
#' default. It measures the Laplace and Wald approximations against the
#' shape of the very objective the fit maximized, and a default prior
#' would change the thing being measured, so it asks for the flat
#' density explicitly rather than inheriting whatever this function
#' defaults to.
#'
#' @section Sampling from a formula:
#' `frm_sample(bf(y ~ x + (1 | g)), data = dd, family = gaussian())`
#' parses, assembles and tapes exactly as [frmtmb::frm()] does, stops before the
#' optimizer, and hands the objective to Stan. Every pre-optimizer
#' refusal still applies (REML, quadrature, the mixture and [frmtmb::hmm()]
#' guards). There is no mode, so the default `init` is `"random"`:
#' Stan's own overdispersed initialization on the unconstrained scale,
#' inside any `lower`/`upper` bounds.
#'
#' The returned object supports the whole draws surface -
#' [summary()], [frmtmb::fixef()], [frmtmb::VarCorr()], [frmtmb::ranef()], [frmtmb::hypothesis()],
#' [posterior_epred()], [posterior_predict()], [posterior_linpred()],
#' [frmtmb::pp_check()], [frmtmb::as_draws()] - because those read the model frame and
#' one draw at a time. The methods that report a maximum-likelihood
#' quantity refuse instead of inventing one: [check_laplace()] (it
#' compares NUTS against a mode that does not exist here), and on the
#' embedded object reachable as `x$fit`, [summary()], [vcov()],
#' [confint()], [logLik()], [frmtmb::fixef()], [frmtmb::ranef()], [frmtmb::VarCorr()],
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
#' ([posterior_epred()], [frmtmb::ranef()], [log_lik()], [frmtmb::loo()],
#' [frmtmb::conditional_effects()], [frmtmb::hypothesis()]) reads them without knowing
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
#' In practice that means a default-prior run, on EITHER route, whose
#' priors cover every standard deviation and every correlation,
#' non-centers
#' `(1 | g)` and any block written one term at a time, `diag()` and
#' `homdiag()` blocks, [mgcv::s()] smooths, `equalto()`, `gr(cov = )`,
#' and the CORRELATED blocks `(x | g)`, `cs()`, `ar1()` and `hetar1()`.
#' A `k =` Hilbert-space `gp()` block stays
#' centered on the formula route even though its factor is diagonal:
#' its LENGTHSCALES share the block's `theta`, the default priors cover
#' only standard deviations, and the gate wants every parameter of a
#' block priored. Prior the whole block by hand
#' (`set_prior(class = "theta")`) to non-center it.
#' `rr()` is already non-centered by
#' construction, since its own coefficients are the standard normal
#' factors. What stays centered is a run whose variance parameters are
#' flat: `prior = "flat"`, and [check_laplace()], which asks for that
#' density on purpose.
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
#'   seeds: under flat priors the non-centered chain walks `theta` to
#'   -1e15 at a bulk-ESS of 1; under the default `student_t(3, 0, s)`,
#'   which makes that tail integrable, 174 to 284 against 3 to 48
#'   centered.
#' - A FLAT PRIOR on the block's CORRELATION, which is the same
#'   argument: flat on frmtmb's unbounded correlation parameter is
#'   `(1 - rho^2)^-3/2`, improper, with all its mass at `|rho| = 1`.
#'   Before 0.39 that kept every correlated block centered, because no
#'   proper correlation prior existed to close the tail; `lkj(eta)`
#'   now does, and it is a default, so
#'   `(Days | Subject)`, `cs()`, `ar1()` and `hetar1()` non-center
#'   like any other block.
#' - A Student-t latent (`gr(dist = "student")`): a scale mixture, not a
#'   linear factor.
#' - `car()`, `spde()` and `gr(prec = )`: sparse precisions whose factor
#'   is dense.
#' - The exact `gp()` and the spatial covariances (`ou`, `exp`, `gau`,
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
#' Both routes default to brms 2.23's own weakly-informative
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
#' ([frmtmb::lognormal()] and its relatives, whose `mu` link is spelled
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
#'   (see [frmtmb::set_prior()] for what the density is on frmtmb's
#'   parameters).
#' - ORDINAL THRESHOLDS stay flat. brms priors them
#'   `student_t(3, 0, 2.5)` under its `Intercept` class; frmtmb keeps
#'   them in the `thres` extra-parameter vector, which is not a design
#'   column and which [frmtmb::set_prior()]'s class vocabulary
#'   (`b`/`Intercept`/`sd`/`theta`) cannot address. An ordinal model
#'   still gets its `sd` defaults, and the message names the gap.
#' - The `shape`, `phi` and `nu` dispersion parameters stay flat: brms
#'   gives them gamma and inverse-gamma defaults, which [frmtmb::set_prior()]
#'   does not carry.
#' - NONLINEAR parameters stay flat, which is what brms does: their
#'   coefficients are class `b`, and the response's median and mad say
#'   nothing about a rate or a shape sitting inside a nonlinear body.
#'   Write them with `set_prior(nlpar = )`, as the brms nonlinear
#'   vignette does.
#' - MULTIVARIATE models get no defaults at all: the default location
#'   and scale are read off ONE response, and frmtmb does not read them
#'   per response. `set_prior(resp = )` addresses one response, so
#'   they can be written by hand.
#'
#' *Overriding and opting out.* A `set_prior()` specification takes over
#' the classes it names and leaves the other defaults in place, which is
#' brms's partial-override rule; a named list of prior objects takes
#' over exactly the internal parameters it names. `prior = "flat"`
#' turns the defaults off entirely and samples the likelihood, which
#' warns when the model has variance components: their flat-prior
#' posteriors need not be proper, and neither the chains nor Rhat can
#' see that. On a MAP fit `"flat"` reaches further than the defaults:
#' that fit's own prior is taped into its objective, and the bare
#' likelihood is what `"flat"` names, so the penalty is rebuilt out and
#' the call says which prior it dropped.
#'
#' *Three sources, most explicit first.* On a fitted model the stack is
#' this call's `prior =`, then the fit's own prior (`frm(prior = )`
#' makes a MAP fit, and that prior is part of the model the fit is),
#' then the defaults. Each addressed slot is settled by the most
#' explicit source that names it, one slot at a time: a call-level
#' `set_prior(class = "sd")` replaces the fit's `sd` prior and the `sd`
#' default while the `Intercept` default stays. `prior_summary()` on
#' the returned draws prints what the stack came to.
#'
#' @param fit A `frmtmb_fit`, or a `bf()`/formula to assemble and sample
#'   directly (then `data` is required).
#' @param data A data frame of model data, when `fit` is a formula.
#' @param family Family, when `fit` is a plain formula that does not
#'   carry one (`frm_sample(bf(y ~ x), data = dd, family = poisson())`;
#'   the `+` spelling `bf(y ~ x) + poisson()` works too).
#' @param data2,start,control,na.action,REML As in [frmtmb::frm()]; used only on
#'   the formula path.
#' @param ... Passed to [tmbstan::tmbstan()] (`chains`, `iter`,
#'   `laplace`, `cores`, ...). `cores` parallelizes over chains on
#'   every platform. On Windows the chains run on socket workers, each
#'   of which rebuilds the tape from the serialized objective closure
#'   (tmbstan retapes on the worker; the closures are self-contained),
#'   giving draws identical to a sequential run at the same seed. The
#'   per-worker startup (a new R process, package load, retape) is a
#'   fixed several seconds, so short chains gain nothing; long chains
#'   approach a per-chain speedup. A core count inherited from
#'   `options(mc.cores)`, which rstan reads when `cores` is not given,
#'   behaves the same way.
#' @param prior Priors: a [frmtmb::set_prior()] specification, or a named list
#'   of prior objects (see [frmtmb::prior_normal()]) whose names are parameter
#'   names as in the draws (or whole components: `"beta"`, `"theta"`,
#'   ...), or the string `"flat"`. The objective is re-taped with the
#'   prior terms added; a fitted model itself is unchanged. On both
#'   paths the brms default priors apply to whatever this argument, and
#'   a MAP fit's own prior, leave alone (see Default priors), and
#'   `prior = "flat"` opts out of them entirely. A `brmsprior` object
#'   built by brms's own `prior()` is translated row by row. The
#'   argument takes brms's spelling, `prior`; the `priors` of releases
#'   before 0.43 is gone rather than aliased, and because this
#'   function's `...` would otherwise swallow it, the old name is
#'   refused by name.
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
#' @param .diagnostic Internal, and named with a dot because there is no
#'   reason to set it by hand. `TRUE` turns the default priors off and
#'   silences the disclosure, and leaves the fit's own objective alone:
#'   a MAP fit's penalty is taped into it, so no prior is resolved
#'   again and no tape is rebuilt. [check_laplace()] sets it, because
#'   it measures the Laplace and Wald approximations against the
#'   density the fit maximized and a default prior would change what is
#'   being measured. Write `prior = "flat"` for the user-facing
#'   opt-out, which goes further and rebuilds a MAP fit's objective
#'   without its penalty.
#' @return An object of class `frmtmb_draws`: list with the `stanfit`,
#'   a draws matrix with named columns (`as.matrix()` method), the
#'   originating fit, and, when any block was non-centered, a
#'   `reparam` note saying which. The `stanfit` holds the parameters as
#'   Stan sampled them, so on a non-centered run its random-effect
#'   columns are `z` while the draws matrix holds `b`.
#' @section Multimodal posteriors:
#'   For [frmtmb::mixture()] fits the posterior is multimodal by construction
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
#'                   prior = set_prior("exponential(1)", class = "sd"))
#' prior_summary(ds3)
#' }
#' }
#' @export
frm_sample <- function(fit, data = NULL, family = NULL, ...,
                       prior = NULL, lower = NULL,
                       upper = NULL, init = NULL,
                       init_jitter = 0.25, reparameterize = TRUE,
                       data2 = list(), start = NULL,
                       control = frmtmb_control(),
                       na.action = stats::na.omit, REML = FALSE,
                       .diagnostic = FALSE) {
  # `...` goes to tmbstan, which would take the retired spelling as an
  # unknown sampler option and run the model with no priors at all
  refuse_retired_priors(list(...), "frm_sample()")
  prior <- as_priorlist(prior)
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
  # BOTH routes sample a POSTERIOR, so both default to brms's
  # weakly-informative priors: a literal translation of a prior-free
  # brm() call onto a fitted model otherwise pins the chain at whatever
  # boundary mode maximum likelihood found. Precedence is
  # most-explicit-wins per slot - this call's prior =, then the fit's
  # own (a MAP fit's), then the defaults - and prior = "flat" opts out
  # of all of it.
  #
  # `own` is read before it is overwritten, and `.diagnostic` leaves it
  # out of the stack on purpose: a MAP fit's penalty is taped INTO
  # fit$obj at fit time, so resolving it again here would rebuild a
  # tape that check_laplace() already has. See check_laplace().
  own <- fit[["prior"]]
  rp <- sample_resolve_priors(fit, prior,
                              base = if (!isTRUE(.diagnostic)) own,
                              defaults = !isTRUE(.diagnostic))
  fit$prior <- rp$effective
  ri <- rp$ri
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
  # every prior-carrying route rebuilds from build_objective(), the
  # BARE likelihood, and adds the resolved entries once; fit$obj is
  # sampled untouched only when nothing has to be added to it
  if (length(ncp$idx)) {
    obj <- ncp_objective(fit, ncp$idx, ri$entries %||% list())
  } else if (!is.null(ri) && length(ri$entries)) {
    obj <- prior_augmented_obj(fit, ri$entries)
  } else if (isTRUE(rp$flat) && !isTRUE(.diagnostic) && !is.null(own)) {
    # the opt-out has something to SUBTRACT here: a MAP fit's own prior
    # is part of fit$obj, and "flat" on the user-facing argument means
    # the likelihood without it. Under `.diagnostic` it means the
    # weaker thing - no prior ADDED to the density the fit maximized -
    # which is fit$obj as it stands
    obj <- prior_augmented_obj(fit, list())
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
              "starts the chains there. The default priors usually ",
              "pull a chain off that boundary; if this one stays ",
              "pinned, use init = \"random\" or a tighter prior =",
              call. = FALSE)
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
  # rstan runs parallel chains on PSOCK workers on Windows. The tape's
  # external pointer dies in serialization, but tmbstan ships the cure:
  # its sampling method calls fn() on the worker, which retapes from
  # the objective closure, and frmtmb's generated closures are
  # self-contained (data baked in, namespace references loading frmtmb
  # and RTMB on deserialization), so each worker rebuilds its own tape.
  # Verified seed-for-seed identical to sequential chains. The only
  # cost is startup, hence the note rather than a fallback.
  if (.Platform$OS.type == "windows" && stan_cores(args) > 1 &&
      (args$chains %||% 4) > 1) {
    message("parallel chains on Windows start one R process per core ",
            "and rebuild the tape in each: expect several seconds of ",
            "startup before sampling. Short chains may run faster ",
            "with cores = 1")
  }
  if (!is.null(bounds)) {
    args$lower <- bounds$lower
    args$upper <- bounds$upper
  }
  check_tmbstan_build("frm_sample()")
  sf <- do.call(tmbstan::tmbstan, args)
  if (!length(sf@sim) || !length(sf@sim$samples)) {
    # Deliberately generic about the cause. The known case is a tape
    # that calls an external solver, and the package supplying such a
    # tape documents its own failure; naming that package from here
    # would tie this message to something frmtmb.sample does not depend
    # on, for a pointer that package's own documentation gives better.
    stop("frm_sample(): the sampler returned no draws (rstan printed ",
         "the cause above). A known case: a tape that calls an ",
         "external solver can fail inside tmbstan even at the fitted ",
         "optimum, where the same tape optimizes without complaint; ",
         "the package supplying such a tape documents that. A solver ",
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
#' **It samples the density the fit maximized, and no other.** It hands
#' Stan `fit$obj` as it stands. For an ordinary fit that objective is
#' the LIKELIHOOD, so the run is flat; for a MAP fit (`frm(prior = )`)
#' the penalty is taped INTO that objective, so the run carries it,
#' which is right: the mode and the Wald standard errors this function
#' compares against come from the same penalized Hessian. "Flat" here
#' means no prior ADDED, not the fit's own penalty stripped - for the
#' bare likelihood of a MAP fit, write `frm_sample(fit, prior =
#' "flat")` instead, and read it as a different question.
#'
#' What this function never adds is [frm_sample()]'s DEFAULT priors.
#' Those defaults
#' apply on both of `frm_sample()`'s routes, and this function opts out
#' of them explicitly (`frm_sample(.diagnostic = TRUE)`) rather than
#' inheriting whatever the default is: a prior nothing in the fit ever
#' saw would change the very thing being measured.
#'
#' That is also why it samples CENTERED on an ordinary fit.
#' `frm_sample()`'s non-centered
#' parameterization (see Reparameterization there) is offered only for
#' blocks whose variance parameters carry a prior, and with the
#' defaults off none do: the
#' flat prior that makes the comparison meaningful is exactly the one
#' that leaves a flat tail at `sd = 0` for a non-centered chain to walk
#' into. So the reparameterization default costs this function nothing
#' and changes nothing about it. Give the variance parameters a prior
#' through `prior =`
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
#' @details
#' The srr standard this function used to carry (RE1.4, "the
#' assumptions the model rests on are documented and testable") is
#' claimed by `frmtmb::frm()` now that the two live in different
#' packages: the assumption is the FIT's, and core has to state it
#' whether or not this package is installed. See the Laplace
#' approximation section of `?frmtmb::frm`.
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
  # explicit, not inherited: frm_sample() defaults to the brms priors on
  # both of its routes now, and this function must keep sampling the
  # density the fit MAXIMIZED - otherwise it would compare NUTS on one
  # posterior against a mode and Wald errors from another. That density
  # is the flat likelihood for an ordinary fit and the fit's own
  # penalized one for a MAP fit, which is what the flag leaves standing
  ds <- frm_sample(fit, chains = chains, iter = iter, ...,
                   .diagnostic = TRUE)
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
#' It hands over the fit's OWN objective and nothing else: no default
#' priors, no non-centering, no named draws. [frm_sample()] is the
#' route that applies brms's default priors and returns the draws
#' surface; this one is the escape hatch to tmbstan's own arguments.
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
