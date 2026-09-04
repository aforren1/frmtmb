#' Define a model family
#'
#' Constructs a family object for use with [frm()]. The log-density
#' function must be vectorized and AD-compatible: it is evaluated on RTMB
#' 'advector' objects during taping, so it must use RTMB-overloaded
#' operations (RTMB and RTMBdist `d*` functions, plain arithmetic) and must
#' not branch on parameter values.
#'
#' @param family Character name of the family.
#' @param dpars Character vector of distributional parameter names. The
#'   first entry must be `"mu"`.
#' @param links Named list mapping each dpar to a link name (see
#'   `frmtmb:::frmtmb_links`) or a link object.
#' @param lpdf Function `(y, dpars, aterms)` returning the vectorized
#'   log-density. `dpars` is a named list of advector vectors; `aterms` is a
#'   named list of numeric addition-term values (for example `trials`).
#' @param valid_y Optional function `(y, aterms)` that signals an error for
#'   invalid responses. Called once at assembly time.
#' @param init_dpars Optional named list of functions `(y, aterms)` giving a
#'   response-scale starting value per dpar (applied to the intercept
#'   through the link).
#' @param type One of `"continuous"`, `"discrete"`, `"ordinal"`,
#'   `"categorical"`. The last two say the modelled response is a
#'   distribution over `1..K` categories rather than a number, so
#'   [fitted()] returns an `n x K` probability matrix; `"ordinal"`
#'   shares one latent predictor across the categories and
#'   `"categorical"` gives each non-reference category its own.
#' @param post Named list of numeric helper functions used by
#'   [fitted()], [predict()] and [residuals()]: `mean_fn(dpars,
#'   aterms)` (the response mean), `var_fn(dpars, aterms)` (for pearson
#'   residuals) and `dev_fn(y, dpars, aterms)` (the unit deviance
#'   `2 * (loglik of the saturated fit - loglik at the fitted value)`,
#'   for `residuals(type = "deviance")`). A family that omits one is
#'   refused by the method that needs it.
#' @param sim Optional numeric simulator `(dpars, aterms, n)` returning `n`
#'   response draws; used by [simulate()], `posterior_predict()` and
#'   [frm_simulate()]. It is stateless and rowwise: it sees the
#'   distributional parameters and nothing else. A family whose extra
#'   parameters (`extra_pars`) enter the draw declares a fourth
#'   argument `extra` instead.
#' @param sim_ctx Optional structured simulator `(ctx)` for a family
#'   whose draws are not rowwise (see Structured simulators). It takes
#'   precedence over `sim`.
#' @param sim_refusal Optional one-sentence reason why the family has no
#'   simulator, appended to the refusal each entry point raises. Use it
#'   when the omission is a decision rather than a gap.
#' @param primary_dpars Which dpars receive the main model formula
#'   (default `"mu"`). Families with several location predictors (for
#'   example multinomial's per-category `mu2`, `mu3`, ...) list them all;
#'   these live in the `beta` parameter vector and are integrated out
#'   under REML.
#' @param lcdf Optional vectorized AD log-safe CDF `(q, dpars, aterms)`
#'   returning probabilities; enables `cens()` and `trunc()` addition
#'   terms.
#' @param extra_pars Optional function `(y, aterms)` returning a named
#'   list of numeric starting vectors for family-level parameters outside
#'   the dpar system (for example ordinal thresholds). They join the
#'   parameter template under their own names and reach `lpdf` as its
#'   fourth argument.
#' @param drop_intercept If `TRUE`, the intercept column is removed from
#'   the main formula's design matrix (ordinal families: thresholds take
#'   its place).
#' @param structure Optional [frmtmb_structure()] for a family whose
#'   likelihood does not factorize over rows (a group-level [mixture()],
#'   a hidden Markov chain). It carries the non-rowwise log-likelihood,
#'   the frame block that likelihood reads, and the capability flags
#'   that say which post-fit methods the family can answer. `lpdf` stays
#'   required even then, for the rowwise contract, and may be a stub
#'   that refuses.
#' @return An object of class `frmtmb_family`.
#' @section Structured simulators:
#' Some families cannot draw a response one row at a time: a group-level
#' [mixture()] draws one class per group, an [hmm()] walks a Markov
#' chain per sequence, and a [mixture_mvn()] draw needs the class
#' covariances, which are family-level extras rather than dpars. Those
#' families supply `sim_ctx(ctx)` instead of `sim(dpars, aterms, n)`.
#'
#' `ctx` is a list with `fit` (any object carrying `spec`, `frame` and
#' `estimates` - a fitted model, one posterior draw, or the de novo
#' shim), `family`, `rspec`, `resp`, `dpars` (the evaluated numeric
#' distributional parameters), `aterms`, `n`, `extra` (the family-level
#' extra parameters) and the frame structures `autocor` and `block`
#' (the structured family's own data; see [frmtmb_structure()]).
#' Read its fields with `[[ ]]`.
#'
#' The same `sim_ctx()` serves [simulate()], `posterior_predict()` and
#' [frm_simulate()]. Because a structured draw covers whole sequences or
#' groups, `trunc()` rejection and `newdata` cannot apply to it and are
#' refused.
#' @examples
#' # a custom family is a plain R log-density over taped parameters
#' dd <- data.frame(y = rbinom(100, 5, 0.4),
#'                  size = 5, x = rnorm(100))
#' fam <- custom_family(
#'   "vbinom", dpars = "mu", links = list(mu = "logit"),
#'   lpdf = function(y, dpars, aterms) {
#'     RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
#'   },
#'   type = "discrete"
#' )
#' fit <- frm(bf(y | vint(size) ~ x) + fam, data = dd)
#' fixef(fit)
#' @srrstats {G2.0,G2.1} The family contract is asserted on both length
#'   and type before anything is built: `family` must be a length-one
#'   character vector, `dpars` a character vector of length at least one,
#'   and `lpdf` a function. A family supplied as a `stats::family()`
#'   object, a family constructor, or a name is dispatched to a single
#'   internal representation, and an unrecognized value errors naming the
#'   supported families.
#' @srrstats {RE4.12} The transform used on each linear predictor and its
#'   inverse are both available. Every distributional parameter carries a
#'   link with `linkfun`, `linkinv`, and `mu_eta` (the derivative), and
#'   `links` selects them per parameter. They are reachable through the
#'   fit with `family()`, through `insight::link_function()` and
#'   `insight::link_inverse()`, and are applied by `predict(type =)` to
#'   move between the link and response scales. Link functions are
#'   written out over plain arithmetic rather than taken from
#'   `stats::make.link()`, because the latter clamps at C level in ways
#'   the AD tape cannot see.
#'
#' @section Tape-safe scope:
#' `lpdf` and `lcdf` run with RTMB's tape-safe `c()`, `[<-` and
#' `diag<-` in scope automatically (the
#' `"c" <- RTMB::ADoverload("c")` boilerplate is spliced in unless the
#' function already binds it), so base spellings keep the
#' automatic-differentiation class. A helper the density CALLS still
#' needs its own bindings: lexical scope does not travel into other
#' functions.
#'
#' @export
frmtmb_family <- function(family, dpars, links, lpdf, valid_y = NULL,
                          init_dpars = list(), type = "continuous",
                          post = list(), sim = NULL, sim_ctx = NULL,
                          sim_refusal = NULL,
                          primary_dpars = "mu", lcdf = NULL,
                          extra_pars = NULL, drop_intercept = FALSE,
                          structure = NULL) {
  stopifnot(is.character(family), length(family) == 1,
            is.character(dpars), length(dpars) >= 1,
            is.function(lpdf))
  # `type` selects the response check, the residual scale and what
  # "response" means to predict(), so an unrecognized string used to
  # produce a family that silently behaved as none of them
  check_string_choice(type, "type",
                      c("continuous", "discrete", "ordinal", "categorical"))
  check_flag(drop_intercept, "drop_intercept")
  if (!all(primary_dpars %in% dpars)) {
    stop("`primary_dpars` must be a subset of `dpars`", call. = FALSE)
  }
  if (!setequal(names(links), dpars)) {
    stop("`links` must name every dpar exactly once", call. = FALSE)
  }
  links <- lapply(links, get_link)
  # user-written densities run with the AD overloads in scope; a
  # function that already binds them itself is left untouched
  lpdf <- ad_overload_fn(lpdf)
  if (!is.null(lcdf)) lcdf <- ad_overload_fn(lcdf)
  if (!is.null(structure) && !inherits(structure, "frmtmb_structure")) {
    stop("frmtmb_family(structure =) must come from frmtmb_structure(), ",
         "which is what declares a likelihood that does not factorize ",
         "over rows", call. = FALSE)
  }
  # base::structure(), spelled out: the `structure` ARGUMENT above masks
  # the base function inside this body
  base::structure(
    list(family = family, dpars = dpars, links = links, lpdf = lpdf,
         valid_y = valid_y, init_dpars = init_dpars, type = type,
         post = post, sim = sim, sim_ctx = sim_ctx,
         sim_refusal = sim_refusal, primary_dpars = primary_dpars,
         lcdf = lcdf, extra_pars = extra_pars,
         drop_intercept = isTRUE(drop_intercept),
         structure = structure),
    class = "frmtmb_family"
  )
}

#' @rdname frmtmb_family
#' @export
custom_family <- frmtmb_family

#' Make a response validator that rejects any value less than or equal
#' to zero. The families with support on `(0, Inf)` use it so a bad
#' response fails at assembly time with the family name in the message.
#'
#' @noRd
positive_y <- function(name) {
  force(name)
  function(y, aterms) {
    if (any(y <= 0)) {
      stop(name, ": response must be strictly positive", call. = FALSE)
    }
  }
}

#' Make a response validator that rejects a response which is not a
#' non-negative integer. The count families use it.
#'
#' @noRd
count_y <- function(name) {
  force(name)
  function(y, aterms) {
    if (any(y < 0) || any(y != round(y))) {
      stop(name, ": response must be non-negative integers", call. = FALSE)
    }
  }
}

# --- Robust (linear-predictor scale) density inputs ------------------
#
# The dpar contract hands a log-density the INVERSE-LINKED value, which
# is what every numeric post-fit path wants. On the tape it throws away
# the far tail: plogis(40) is exactly 1, so `1 - mu` is exactly 0 and a
# binomial log-density is -Inf where the truth is -40, with a gradient
# of NaN. `build_objective()` therefore stores the linear predictor
# beside each dpar under `.eta_<dpar>`, and the accessors below recover
# the quantity a robust density needs from it. Off the tape the entries
# are absent, every accessor returns NULL, and the family falls back to
# the plain form - which is correct there, because nothing off the tape
# is differentiated.

#' The log-odds behind a probability dpar, on the linear-predictor
#' scale, or NULL when the objective did not supply the linear predictor
#' or the dpar's link has no exact log-odds form (identity, inverse).
#'
#' @noRd
robust_logit <- function(dpars, link, name = "mu") {
  e <- dpars[[paste0(".eta_", name)]]
  if (is.null(e) || is.null(link$logit_eta)) return(NULL)
  link$logit_eta(e)
}

#' The log mean behind a positive-mean dpar, on the linear-predictor
#' scale, or NULL as in [robust_logit()]. Only the log link supplies it,
#' where it is the linear predictor itself.
#'
#' @noRd
robust_logmu <- function(dpars, link, name = "mu") {
  e <- dpars[[paste0(".eta_", name)]]
  if (is.null(e) || is.null(link$log_eta)) return(NULL)
  link$log_eta(e)
}

#' The log of a dpar whose link the FAMILY fixes at log (`shape`, `phi`,
#' `sigma`): the linear predictor is that log exactly. Do not call it on
#' a dpar whose link the user chooses.
#'
#' @noRd
log_dpar <- function(dpars, name) {
  e <- dpars[[paste0(".eta_", name)]]
  if (is.null(e)) log(dpars[[name]]) else e
}

#' `log(p)` and `log(1 - p)` for a mixture gate (`zi`, `hu`) or any other
#' dpar the family pins to the logit link. Both terms stay finite and
#' exactly differentiable at a gate the optimizer has pushed against 0 or
#' 1, which is where a separated zero-inflation predictor lives.
#'
#' @noRd
gate_logs <- function(dpars, name) {
  e <- dpars[[paste0(".eta_", name)]]
  if (is.null(e)) {
    p <- dpars[[name]]
    return(list(l = log(p), l1m = log(1 - p)))
  }
  list(l = log_inv_logit(e), l1m = log1m_inv_logit(e))
}

#' The pair `(mu, 1 - mu)` for a density that needs both, from the
#' log-odds when it is available. `1 - plogis(eta)` is 0 for eta beyond
#' 37 and carries one significant digit from 30 on; `exp(log(1 - p))`
#' off the log-odds is exact over the whole line.
#'
#' @noRd
mu_pair <- function(dpars, link, name = "mu") {
  lo <- robust_logit(dpars, link, name)
  if (is.null(lo)) {
    mu <- dpars[[name]]
    return(list(p = mu, q = 1 - mu))
  }
  list(p = exp(log_inv_logit(lo)), q = exp(log1m_inv_logit(lo)))
}

#' `P(a < Z < b)` for standard normal bounds. `pnorm(b) - pnorm(a)`
#' loses every significant digit when both bounds sit in the upper tail,
#' which is exactly where truncated means are computed; the mirrored
#' form keeps full precision there.
#'
#' @noRd
pnorm_diff <- function(a, b) {
  ifelse(a > 0,
         stats::pnorm(-a) - stats::pnorm(-b),
         stats::pnorm(b) - stats::pnorm(a))
}

#' Residual SD including a known `se()` component (meta-analysis):
#' `se()` alone replaces `sigma`; `se(x, sigma = TRUE)` adds them in
#' quadrature.
#'
#' @noRd
resid_sd <- function(sigma, aterms) {
  if (is.null(aterms$se)) return(sigma)
  if (isTRUE(aterms$se_sigma)) sqrt(sigma^2 + aterms$se^2) else aterms$se
}

#' `trunc()` bounds from an aterm-value list as a pair of length-n
#' numeric vectors, or `NULL` when the response is not truncated. An
#' absent bound is the family's unbounded default; the per-family
#' truncated means all handle the infinities.
#'
#' @noRd
trunc_bounds <- function(aterms, n) {
  lb <- aterms[["trunc_lb"]]
  ub <- aterms[["trunc_ub"]]
  if (is.null(lb) && is.null(ub)) return(NULL)
  list(lb = rep(lb %||% -Inf, length.out = n),
       ub = rep(ub %||% Inf, length.out = n))
}

# On a truncated response this is E[Y | lb <= Y <= ub], the quantity
# fitted(), residuals() and predict(type = "response") report; per-dpar
# predictions stay untruncated, because they describe the latent
# parameter rather than the observable.

#' @rdname frmtmb-extension-api
#' @export
response_mean <- function(fam, dpars, aterms) {
  mu <- if (!is.null(fam$post$mean_fn)) {
    fam$post$mean_fn(dpars, aterms)
  } else {
    dpars$mu
  }
  tb <- trunc_bounds(aterms, length(mu))
  if (is.null(tb)) return(mu)
  tmf <- fam$post$trunc_mean_fn
  if (is.null(tmf)) {
    stop("Family '", fam$family, "' has no truncated mean; fitted(), ",
         "residuals() and predict(type = \"response\") would report the ",
         "untruncated mean", call. = FALSE)
  }
  tmf(dpars, aterms, tb$lb, tb$ub)
}

# --- Unit deviances ---
#
# d_i = 2 * (loglik of the saturated fit - loglik at the fitted value),
# with the dispersion parameter held at its estimate, which is the
# quantity glm() calls the unit deviance. Families that are exponential
# dispersion models reproduce stats::glm()'s dev.resids exactly; the
# rest (negbinomial, beta, tweedie) use the same saturated-likelihood
# definition at a fixed shape.

#' `y log(y / mu)` under the `0 log 0 = 0` convention the saturated fit
#' needs at a zero count.
#'
#' @noRd
ylogy_mu <- function(y, mu) ifelse(y > 0, y * log(y / mu), 0)

#' Negative-binomial unit deviance at a given size (shape). nbinom1
#' feeds it the row's own size `mu / phi`.
#'
#' @noRd
nbinom_deviance <- function(y, mu, size) {
  2 * (ylogy_mu(y, mu) - (y + size) * log((y + size) / (mu + size)))
}

#' Binomial unit deviance on COUNTS out of `size` trials; both terms
#' vanish at the boundaries `y = 0` and `y = size`.
#'
#' @noRd
binomial_deviance <- function(y, mu, size) {
  2 * (ylogy_mu(y, size * mu) + ylogy_mu(size - y, size * (1 - mu)))
}

#' Gamma unit deviance; also the exponential one (shape fixed at 1).
#'
#' @noRd
gamma_deviance <- function(y, mu) 2 * ((y - mu) / mu - log(y / mu))

#' Families that define a unit deviance, for the `residuals()` refusal
#' message. Read off the registry so the list cannot drift; the
#' constructors that need arguments (multinomial) simply drop out.
#'
#' @noRd
deviance_family_names <- function() {
  nms <- names(family_registry)
  ok <- vapply(nms, function(nm) {
    fam <- tryCatch(family_registry[[nm]](), error = function(e) NULL)
    !is.null(fam) && !is.null(fam$post$dev_fn)
  }, TRUE)
  sort(unique(nms[ok]))
}

#' Deviance residuals: `sign(y - E[Y]) * sqrt(w_i d_i)`. Weights
#' multiply the unit deviance, as in `glm()`. A truncated or censored
#' response is refused: the fitted likelihood is not the family's own
#' density there, so the saturated comparison the unit deviance is built
#' on does not describe the model that was estimated.
#'
#' @noRd
deviance_residuals <- function(fam, y, dpars, aterms, n) {
  dev <- fam$post$dev_fn
  if (is.null(dev)) {
    stop("residuals(type = \"deviance\") is not available for family '",
         fam$family, "': it has no standard unit deviance. Families ",
         "with one: ", paste(deviance_family_names(), collapse = ", "),
         ". Use type = \"osa\" or dharma_residuals() instead.",
         call. = FALSE)
  }
  if (!is.null(trunc_bounds(aterms, n))) {
    stop("residuals(type = \"deviance\") is not defined for a trunc()ed ",
         "response: the unit deviance compares against the untruncated ",
         "family, not the likelihood the model was fitted with. Use ",
         "type = \"osa\", which builds its CDF on [lb, ub]", call. = FALSE)
  }
  if (!is.null(aterms$cens) && any(aterms$cens != 0)) {
    stop("residuals(type = \"deviance\") is not defined on a cens()ed ",
         "response: a censored row observes an event, not a value, so it ",
         "has no unit deviance. Use type = \"osa\"", call. = FALSE)
  }
  d <- dev(y, dpars, aterms)
  w <- aterms$weights %||% 1
  # rounding can push an exactly saturated row a few ulps below zero
  sign(y - response_mean(fam, dpars, aterms)) * sqrt(pmax(w * d, 0))
}

#' Per-observation subset of dpar / aterm vectors, for resampling the
#' rows a rejection step has not accepted yet.
#'
#' @noRd
subset_obs <- function(x, idx, n) {
  lapply(x, function(v) {
    if (is.numeric(v) && length(v) %in% c(1L, n)) {
      rep(v, length.out = n)[idx]
    } else {
      v
    }
  })
}

#' One simulated response vector, respecting `trunc()` bounds by
#' rejection: out-of-bounds draws are redrawn until every row is inside
#' its own interval. Bounds that exclude nearly all the family's mass
#' never converge, so the iteration cap reports the acceptance rate
#' instead of spinning.
#'
#' @noRd
sim_response <- function(fam, dpars, aterms, n, max_iter = 100L,
                         extra = NULL) {
  # families whose draws need parameters outside the dpar system
  # (ordinal thresholds) declare a fourth argument; the rest keep the
  # three-argument contract
  fam_sim <- if (length(formals(fam$sim)) >= 4L) {
    function(dp, av, nn) fam$sim(dp, av, nn, extra)
  } else {
    fam$sim
  }
  y <- fam_sim(dpars, aterms, n)
  tb <- trunc_bounds(aterms, n)
  if (is.null(tb)) return(y)
  drawn <- n
  bad <- which(y < tb$lb | y > tb$ub)
  it <- 0L
  while (length(bad)) {
    it <- it + 1L
    if (it > max_iter) {
      stop("trunc(): rejection sampling did not fill ", length(bad),
           " of ", n, " rows in ", max_iter, " passes (acceptance rate ",
           format((n - length(bad)) / drawn, digits = 2),
           "). The bounds exclude nearly all of the fitted ",
           "distribution's mass.", call. = FALSE)
    }
    drawn <- drawn + length(bad)
    yb <- fam_sim(subset_obs(dpars, bad, n), subset_obs(aterms, bad, n),
                  length(bad))
    y[bad] <- yb
    bad <- bad[yb < tb$lb[bad] | yb > tb$ub[bad]]
  }
  y
}

# --- the structured simulator contract --------------------------------
#
# `fam$sim(dpars, aterms, n)` is stateless and ROWWISE: it sees one
# column of numbers per distributional parameter and draws each row on
# its own. Three model classes cannot be written that way.
#
#   - family-level EXTRAS. An lca() draw needs the item-profile
#     parameters, which live outside the dpar system. That case already
#     has an answer: a four-argument `sim` receives the extra vector.
#   - GROUP or TIME structure. A group-level mixture draws one class per
#     GROUP, and an hmm() draw walks a Markov chain per SEQUENCE. Both
#     structures are resolved at frame-assembly time and live on the
#     frame, as the structured family's `blocks` entry, not on the
#     family object.
#   - CROSS-ROW dependence. An autocor residual is one multivariate draw
#     per group, so no per-row simulator exists at all.
#
# `fam$sim_ctx(ctx)` is the extension: one function per family, taking a
# CONTEXT that carries the fit-like object as well as the evaluated
# dpars. The context is built the same way from all three entry points -
# `simulate()` on a fit, `posterior_predict()` on one draw (where
# `draws_fit_at()` makes a real fit out of it), and `frm_simulate()` on
# the de novo shim - so a structured family has exactly ONE
# implementation and every entry point reaches it.
#
# `[[ ]]`, never `$`, on the context and on frame fields: `$` partial
# matching once let a `ctx$mix` read return the group structure instead.

#' The simulation context: everything a simulator may read, assembled
#' identically from a fit, from one posterior draw, and from the de novo
#' shim. `fit` is any object carrying `spec`, `frame` and `estimates`.
#'
#' @noRd
sim_context <- function(fit, rspec, dpars, aterms = NULL, n = NULL,
                        extra = NULL) {
  resp <- rspec[["resp_name"]]
  frame <- fit[["frame"]]
  list(fit = fit,
       family = rspec[["family"]],
       rspec = rspec,
       resp = resp,
       dpars = dpars,
       aterms = aterms %||% frame[["aterm_values"]][[resp]] %||% list(),
       n = n %||% frame[["n_obs"]],
       extra = extra %||% fit_extras(fit),
       autocor = frame[["autocor"]][[resp]],
       block = frame_block_of(frame, resp))
}

#' The family's structured simulator, from its structure or, for a
#' family that has one without a structure, from its own slot.
#'
#' @noRd
fam_sim_ctx <- function(fam) {
  fam_structure(fam)[["sim_ctx"]] %||% fam[["sim_ctx"]]
}

#' Whether a family can simulate at all, by either half of the contract.
#'
#' @noRd
sim_can <- function(fam) {
  !is.null(fam[["sim"]]) || !is.null(fam_sim_ctx(fam))
}

#' The family's own explanation for having no simulator, appended by
#' whichever entry point refused. Empty for a family that simply has not
#' been given one yet.
#'
#' @noRd
sim_note <- function(fam) {
  note <- fam[["sim_refusal"]]
  if (is.null(note)) "" else paste0(". ", note)
}

#' Whether the draws of this context come whole rather than row by row,
#' which is what makes `trunc()` rejection and `newdata` inapplicable.
#'
#' @noRd
sim_is_structured <- function(ctx) {
  !is.null(ctx[["autocor"]]) || !is.null(fam_sim_ctx(ctx[["family"]]))
}

#' One simulated response for a context: the single implementation every
#' entry point calls. Cross-row structure is resolved first (a residual
#' correlation belongs to the frame, not to the family), then the
#' family's structured simulator, then the rowwise contract.
#'
#' @noRd
sim_draw <- function(ctx) {
  ac <- ctx[["autocor"]]
  sf <- fam_sim_ctx(ctx[["family"]])
  if (is.null(ac) && is.null(sf)) {
    return(sim_response(ctx[["family"]], ctx[["dpars"]], ctx[["aterms"]],
                        ctx[["n"]], extra = ctx[["extra"]]))
  }
  if (!is.null(trunc_bounds(ctx[["aterms"]], ctx[["n"]]))) {
    stop("trunc() cannot be combined with a structured draw (here: '",
         ctx[["family"]][["family"]], "'): a hidden state sequence, a ",
         "group-level latent class and a correlated residual are each ",
         "drawn whole, so a row outside its bounds cannot be redrawn ",
         "on its own and the rejection step has nothing to resample",
         call. = FALSE)
  }
  if (!is.null(ac)) return(sim_autocor_rows(ctx, ac))
  sf(ctx)
}

#' The autocor branch: one multivariate residual draw per group added to
#' the mean predictor, rather than n independent family draws.
#'
#' @noRd
sim_autocor_rows <- function(ctx, ac) {
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  th <- ctx[["fit"]][["estimates"]][["thetaac"]]
  R <- autocor_cor(th[ac[["theta_idx"]]], ac)
  rep(dp[["mu"]], length.out = n) +
    autocor_draw_resid(ac, R, rep(dp[["sigma"]], length.out = n), n,
                       nu = if (isTRUE(ac[["student"]])) dp[["nu"]][1])
}

#' Gaussian family, dpars `mu` and `sigma`. A known `se()` term enters
#' the residual SD, which is what makes meta-analysis work. It supplies
#' a CDF, so `cens()` and `trunc()` apply to it.
#'
#' @noRd
fam_gaussian <- function(link = "identity") {
  frmtmb_family(
    "gaussian",
    dpars = c("mu", "sigma"),
    links = list(mu = link, sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu, resid_sd(dpars$sigma, aterms), log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::pnorm((q - dpars$mu) / resid_sd(dpars$sigma, aterms))
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) resid_sd(dpars$sigma, aterms)^2,
      # sigma is the dispersion, so it divides out of the unit deviance.
      # se() breaks that: the residual sd is row-specific, and a raw
      # squared residual would then compare rows of different precision
      # on one scale. The known variance enters exactly as a glm prior
      # weight sigma^2 / s_i^2, which is 1 without se(); se() alone maps
      # sigma out at 1, so the row weight is the usual 1 / se_i^2.
      dev_fn = function(y, dpars, aterms) {
        s <- resid_sd(dpars$sigma, aterms)
        (y - dpars$mu)^2 * (dpars$sigma / s)^2
      },
      # mu + sigma * (phi(a) - phi(b)) / (Phi(b) - Phi(a))
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        s <- resid_sd(dpars$sigma, aterms)
        a <- (lb - dpars$mu) / s
        b <- (ub - dpars$mu) / s
        dpars$mu + s * (stats::dnorm(a) - stats::dnorm(b)) /
          pnorm_diff(a, b)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnorm(n, dpars$mu, resid_sd(dpars$sigma, aterms))
    }
  )
}

#' Poisson family, single dpar `mu` (the mean). It supplies a CDF and a
#' truncated mean, so `cens()` and `trunc()` apply to it.
#'
#' @noRd
fam_poisson <- function(link = "log") {
  frmtmb_family(
    "poisson",
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars$mu, log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::ppois(q, dpars$mu)
    },
    valid_y = count_y("poisson"),
    init_dpars = list(mu = function(y, aterms) mean(y) + 0.1),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu,
      dev_fn = function(y, dpars, aterms) {
        2 * (ylogy_mu(y, dpars$mu) - (y - dpars$mu))
      },
      # sum_{lb}^{ub} y dpois(y) = mu * (F(ub-1) - F(lb-2)), over the
      # same F(ub) - F(lb-1) normalizer the likelihood uses: the
      # inclusive lower bound keeps its own mass (brms#1903)
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        mu <- dpars$mu
        mu * (stats::ppois(ub - 1, mu) - stats::ppois(lb - 2, mu)) /
          (stats::ppois(ub, mu) - stats::ppois(lb - 1, mu))
      }
    ),
    sim = function(dpars, aterms, n) stats::rpois(n, dpars$mu)
  )
}

#' Binomial family, single dpar `mu` (the success probability). The
#' number of trials comes from the `trials()` addition term and defaults
#' to one, so the response is a count out of `trials`.
#'
#' On the tape it uses `RTMB::dbinom_robust()`, TMB's log-odds
#' parameterization of the same density, which is what glmmTMB fits
#' with. `dbinom(y, n, plogis(eta))` is -Inf from `|eta| = 37` and
#' already wrong in the second decimal at 30, so a separated predictor
#' takes the optimizer into a region where it has no gradient at all.
#'
#' @noRd
fam_binomial <- function(link = "logit") {
  lk <- get_link(link)
  frmtmb_family(
    "binomial",
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      size <- aterms$trials %||% 1
      lo <- robust_logit(dpars, lk)
      if (is.null(lo)) {
        return(RTMB::dbinom(y, size, dpars$mu, log = TRUE))
      }
      RTMB::dbinom_robust(y, size, lo, log = TRUE)
    },
    valid_y = function(y, aterms) {
      size <- aterms$trials %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("binomial: response must be integer counts in [0, trials]",
             call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        size <- aterms$trials %||% 1
        p <- mean(y / size)
        min(max(p, 0.02), 0.98)
      }
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu * (aterms$trials %||% 1),
      var_fn = function(dpars, aterms) {
        (aterms$trials %||% 1) * dpars$mu * (1 - dpars$mu)
      },
      dev_fn = function(y, dpars, aterms) {
        binomial_deviance(y, dpars$mu, aterms$trials %||% 1)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, aterms$trials %||% 1, dpars$mu)
    }
  )
}

#' Gamma family in the mean parameterization, dpars `mu` and `shape`.
#' The scale is `mu / shape`, so `mu` stays the mean at any shape.
#'
#' @noRd
fam_Gamma <- function(link = "log") {
  frmtmb_family(
    "Gamma",
    dpars = c("mu", "shape"),
    links = list(mu = link, shape = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dgamma(y, shape = dpars$shape, scale = dpars$mu / dpars$shape,
                   log = TRUE)
    },
    valid_y = positive_y("Gamma"),
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      shape = function(y, aterms) {
        max(mean(y)^2 / stats::var(y), 0.1)
      }
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu^2 / dpars$shape,
      # 1 / shape is the dispersion and divides out
      dev_fn = function(y, dpars, aterms) gamma_deviance(y, dpars$mu)
    ),
    sim = function(dpars, aterms, n) {
      stats::rgamma(n, shape = dpars$shape, scale = dpars$mu / dpars$shape)
    }
  )
}

#' Lognormal family, dpars `mu` and `sigma`. As in brms, they are the
#' mean and SD on the LOG scale, which is why `mu` takes the identity
#' link. It supplies a CDF and a truncated mean.
#'
#' @noRd
fam_lognormal <- function(link = "identity") {
  frmtmb_family(
    "lognormal",
    dpars = c("mu", "sigma"),
    links = list(mu = link, sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(log(y), dpars$mu, dpars$sigma, log = TRUE) - log(y)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::pnorm((log(q) - dpars$mu) / dpars$sigma)
    },
    valid_y = positive_y("lognormal"),
    init_dpars = list(
      mu = function(y, aterms) mean(log(y)),
      sigma = function(y, aterms) stats::sd(log(y))
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        exp(dpars$mu + dpars$sigma^2 / 2)
      },
      var_fn = function(dpars, aterms) {
        (exp(dpars$sigma^2) - 1) * exp(2 * dpars$mu + dpars$sigma^2)
      },
      # E[Y] * (Phi(b - sigma) - Phi(a - sigma)) / (Phi(b) - Phi(a)),
      # a, b the log-scale standardized bounds
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        sg <- dpars$sigma
        a <- (log(pmax(lb, 0)) - dpars$mu) / sg
        b <- (log(ub) - dpars$mu) / sg
        exp(dpars$mu + sg^2 / 2) * pnorm_diff(a - sg, b - sg) /
          pnorm_diff(a, b)
      }
    ),
    sim = function(dpars, aterms, n) stats::rlnorm(n, dpars$mu, dpars$sigma)
  )
}

#' Student-t family, dpars `mu`, `sigma` and `nu`. The `logm1` link on
#' `nu` holds the degrees of freedom above one. A known `se()` term
#' enters the scale, as in the gaussian family.
#'
#' @noRd
fam_student <- function(link = "identity") {
  frmtmb_family(
    "student",
    dpars = c("mu", "sigma", "nu"),
    links = list(mu = link, sigma = "log", nu = "logm1"),
    lpdf = function(y, dpars, aterms) {
      sd_t <- resid_sd(dpars$sigma, aterms)
      RTMB::dt((y - dpars$mu) / sd_t, df = dpars$nu, log = TRUE) -
        log(sd_t)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y),
      nu = function(y, aterms) 4
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) {
        sd_t <- resid_sd(dpars$sigma, aterms)
        ifelse(dpars$nu > 2, sd_t^2 * dpars$nu / (dpars$nu - 2),
               NA_real_)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars$mu + resid_sd(dpars$sigma, aterms) * stats::rt(n, dpars$nu)
    }
  )
}

#' Negative binomial with quadratic variance (brms negbinomial, glmmTMB
#' nbinom2), dpars `mu` and `shape`. The variance is
#' `mu + mu^2 / shape`.
#'
#' @noRd
fam_negbinomial <- function(link = "log") {
  lk <- get_link(link)
  frmtmb_family(
    "negbinomial",
    dpars = c("mu", "shape"),
    links = list(mu = lk, shape = "log"),
    lpdf = function(y, dpars, aterms) {
      # dnbinom2 forms p = mu / var, which is 0 / 0 = NaN once exp(eta)
      # has underflowed; dnbinom_robust takes log(mu) and
      # log(var - mu) = 2 log(mu) - log(shape) and never divides.
      lmu <- robust_logmu(dpars, lk)
      if (is.null(lmu)) {
        return(RTMB::dnbinom2(y, dpars$mu,
                              dpars$mu + dpars$mu^2 / dpars$shape,
                              log = TRUE))
      }
      RTMB::dnbinom_robust(y, lmu, 2 * lmu - log_dpar(dpars, "shape"),
                           log = TRUE)
    },
    valid_y = count_y("negbinomial"),
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      shape = function(y, aterms) {
        m <- mean(y)
        v <- stats::var(y)
        if (v > m) max(m^2 / (v - m), 0.1) else 10
      }
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu + dpars$mu^2 / dpars$shape,
      dev_fn = function(y, dpars, aterms) {
        nbinom_deviance(y, dpars$mu, dpars$shape)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnbinom(n, size = dpars$shape, mu = dpars$mu)
    }
  )
}

#' Negative binomial with linear variance (glmmTMB nbinom1), dpars `mu`
#' and `phi`. The variance is `mu * (1 + phi)`, so `phi` is a
#' quasi-Poisson style dispersion.
#'
#' @noRd
fam_nbinom1 <- function(link = "log") {
  lk <- get_link(link)
  frmtmb_family(
    "nbinom1",
    dpars = c("mu", "phi"),
    links = list(mu = lk, phi = "log"),
    lpdf = function(y, dpars, aterms) {
      # var - mu = mu * phi, so log(var - mu) = log(mu) + log(phi)
      lmu <- robust_logmu(dpars, lk)
      if (is.null(lmu)) {
        return(RTMB::dnbinom2(y, dpars$mu, dpars$mu * (1 + dpars$phi),
                              log = TRUE))
      }
      RTMB::dnbinom_robust(y, lmu, lmu + log_dpar(dpars, "phi"),
                           log = TRUE)
    },
    valid_y = count_y("nbinom1"),
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      phi = function(y, aterms) {
        m <- mean(y)
        v <- stats::var(y)
        max(v / max(m, 0.1) - 1, 0.1)
      }
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu * (1 + dpars$phi),
      # glmmTMB's convention: the negative-binomial unit deviance with
      # the size held at the FITTED row's mu / phi. Letting the size
      # follow the saturated mean instead is not a deviance at all - the
      # nbinom1 log-likelihood in mu is not maximized at mu = y once the
      # size moves with it, and the difference goes negative.
      dev_fn = function(y, dpars, aterms) {
        nbinom_deviance(y, dpars$mu, dpars$mu / dpars$phi)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnbinom(n, size = dpars$mu / dpars$phi, mu = dpars$mu)
    }
  )
}

#' Beta family in the mean-precision parameterization, dpars `mu` and
#' `phi`. The two shapes are `mu * phi` and `(1 - mu) * phi`, and the
#' response must lie strictly inside `(0, 1)`.
#'
#' @noRd
fam_beta <- function(link = "logit") {
  lk <- get_link(link)
  frmtmb_family(
    "beta",
    dpars = c("mu", "phi"),
    links = list(mu = lk, phi = "log"),
    lpdf = function(y, dpars, aterms) {
      # the SECOND shape is (1 - mu) * phi, and 1 - plogis(eta) is
      # exactly 0 past eta = 37: dbeta at shape 0 is -Inf. Taking the
      # pair off the log-odds keeps both shapes strictly positive.
      mp <- mu_pair(dpars, lk)
      RTMB::dbeta(y, mp$p * dpars$phi, mp$q * dpars$phi, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y <= 0) || any(y >= 1)) {
        stop("beta: response must lie strictly in (0, 1)", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      phi = function(y, aterms) {
        m <- mean(y)
        max(m * (1 - m) / stats::var(y) - 1, 0.5)
      }
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) {
        dpars$mu * (1 - dpars$mu) / (1 + dpars$phi)
      },
      # 2 (ll(y; y, phi) - ll(y; mu, phi)); the terms in y alone cancel
      # (the betareg deviance-residual definition). glmmTMB returns NA
      # for beta, so there is nothing to match there.
      dev_fn = function(y, dpars, aterms) {
        mu <- dpars$mu
        ph <- dpars$phi
        2 * (lgamma(mu * ph) + lgamma((1 - mu) * ph) -
               lgamma(y * ph) - lgamma((1 - y) * ph) +
               (y - mu) * ph * log(y / (1 - y)))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbeta(n, dpars$mu * dpars$phi, (1 - dpars$mu) * dpars$phi)
    }
  )
}

#' Tweedie family for a non-negative response with a mass at zero,
#' dpars `mu`, `phi` and `power`. The `power12` link holds `power`
#' between 1 and 2, the compound Poisson-gamma range.
#'
#' @noRd
fam_tweedie <- function(link = "log") {
  frmtmb_family(
    "tweedie",
    dpars = c("mu", "phi", "power"),
    links = list(mu = link, phi = "log", power = "power12"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dtweedie(y, dpars$mu, dpars$phi, dpars$power, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y < 0)) {
        stop("tweedie: response must be non-negative", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      phi = function(y, aterms) 1,
      power = function(y, aterms) 1.5
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$phi * dpars$mu^dpars$power,
      # the standard Tweedie unit deviance (phi is the dispersion and
      # divides out); at y = 0 the first two terms vanish
      dev_fn = function(y, dpars, aterms) {
        mu <- dpars$mu
        p <- dpars$power
        2 * (y^(2 - p) / ((1 - p) * (2 - p)) -
               y * mu^(1 - p) / (1 - p) +
               mu^(2 - p) / (2 - p))
      }
    )
  )
}

#' Conway-Maxwell-Poisson family for counts that are under- or
#' over-dispersed, dpars `mu` (the mean) and `nu` (the dispersion,
#' with `nu = 1` giving the Poisson).
#'
#' @noRd
fam_compois <- function(link = "log") {
  frmtmb_family(
    "compois",
    dpars = c("mu", "nu"),
    links = list(mu = link, nu = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dcompois2(y, dpars$mu, dpars$nu, log = TRUE)
    },
    valid_y = count_y("compois"),
    init_dpars = list(
      mu = function(y, aterms) mean(y) + 0.1,
      nu = function(y, aterms) 1
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu
    )
  )
}

#' Zero-inflated Poisson, dpars `mu` (the Poisson mean) and `zi` (the
#' probability of a structural zero).
#'
#' @noRd
fam_zi_poisson <- function(link = "log") {
  frmtmb_family(
    "zero_inflated_poisson",
    dpars = c("mu", "zi"),
    links = list(mu = link, zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      # y == 0 is data, so the mixture stays branch-free in parameters.
      # The whole mixture runs in log space: log(1 - zi) is -Inf once
      # the zi predictor separates, and the Poisson's own log P(0) is
      # -mu exactly, with no exp() to underflow.
      i0 <- as.numeric(y == 0)
      g <- gate_logs(dpars, "zi")
      i0 * RTMB::logspace_add(g$l, g$l1m - dpars$mu) +
        (1 - i0) * (g$l1m + RTMB::dpois(y, dpars$mu, log = TRUE))
    },
    valid_y = count_y("zero_inflated_poisson"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars$zi) * dpars$mu,
      var_fn = function(dpars, aterms) {
        (1 - dpars$zi) * dpars$mu * (1 + dpars$zi * dpars$mu)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars$zi) * stats::rpois(n, dpars$mu)
    }
  )
}

#' Zero-inflated negative binomial (nbinom2 variance), dpars `mu`,
#' `shape` and `zi` (the probability of a structural zero).
#'
#' @noRd
fam_zi_negbinomial <- function(link = "log") {
  lk <- get_link(link)
  frmtmb_family(
    "zero_inflated_negbinomial",
    dpars = c("mu", "shape", "zi"),
    links = list(mu = lk, shape = "log", zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      g <- gate_logs(dpars, "zi")
      lmu <- robust_logmu(dpars, lk)
      if (is.null(lmu)) {
        lp0 <- dpars$shape * (log(dpars$shape) -
                                log(dpars$shape + dpars$mu))
        base <- RTMB::dnbinom2(y, dpars$mu,
                               dpars$mu + dpars$mu^2 / dpars$shape,
                               log = TRUE)
      } else {
        # log P(0) = shape * log(shape / (shape + mu)), which is
        # -shape * log(1 + mu / shape) with the ratio taken in logs
        lsh <- log_dpar(dpars, "shape")
        lp0 <- -dpars$shape *
          RTMB::logspace_add(0 * lmu, lmu - lsh)
        base <- RTMB::dnbinom_robust(y, lmu, 2 * lmu - lsh, log = TRUE)
      }
      i0 * RTMB::logspace_add(g$l, g$l1m + lp0) + (1 - i0) * (g$l1m + base)
    },
    valid_y = count_y("zero_inflated_negbinomial"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      shape = function(y, aterms) 1,
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars$zi) * dpars$mu
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars$zi) *
        stats::rnbinom(n, size = dpars$shape, mu = dpars$mu)
    }
  )
}

#' Hurdle Poisson, dpars `mu` and `hu`. `hu` is the probability of a
#' zero and the positive part is a zero-truncated Poisson, so zeros come
#' only from the hurdle.
#'
#' @noRd
fam_hurdle_poisson <- function(link = "log") {
  frmtmb_family(
    "hurdle_poisson",
    dpars = c("mu", "hu"),
    links = list(mu = link, hu = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      g <- gate_logs(dpars, "hu")
      # nonzero part is a zero-truncated poisson. Its normalizer
      # log(1 - exp(-mu)) cancels to log(0) once mu underflows below
      # the double epsilon; expm1 keeps it, and the truth there is
      # simply log(mu).
      lztrunc <- log(-expm1(-dpars$mu))
      i0 * g$l +
        (1 - i0) * (g$l1m + RTMB::dpois(y, dpars$mu, log = TRUE) -
                      lztrunc)
    },
    valid_y = count_y("hurdle_poisson"),
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars$hu) * dpars$mu / (1 - exp(-dpars$mu))
      }
    )
  )
}

#' Bernoulli family for a 0/1 response, single dpar `mu` (the success
#' probability). Robust in the linear predictor, as [fam_binomial()]
#' describes; separation is this family's normal failure mode, so the
#' saturating form is exactly the wrong one here.
#'
#' @noRd
fam_bernoulli <- function(link = "logit") {
  lk <- get_link(link)
  frmtmb_family(
    "bernoulli",
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      lo <- robust_logit(dpars, lk)
      if (is.null(lo)) return(RTMB::dbinom(y, 1, dpars$mu, log = TRUE))
      RTMB::dbinom_robust(y, 1, lo, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (!all(y %in% c(0, 1))) {
        stop("bernoulli: response must be 0/1", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) min(max(mean(y), 0.02), 0.98)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu * (1 - dpars$mu),
      dev_fn = function(y, dpars, aterms) {
        binomial_deviance(y, dpars$mu, 1)
      }
    ),
    sim = function(dpars, aterms, n) stats::rbinom(n, 1, dpars$mu)
  )
}

#' Geometric family, single dpar `mu` (the mean). It is the negative
#' binomial with the shape held at one.
#'
#' @noRd
fam_geometric <- function(link = "log") {
  lk <- get_link(link)
  frmtmb_family(
    "geometric",
    dpars = "mu",
    links = list(mu = lk),
    lpdf = function(y, dpars, aterms) {
      # var - mu = mu^2 (the negative binomial at shape 1)
      lmu <- robust_logmu(dpars, lk)
      if (is.null(lmu)) {
        return(RTMB::dnbinom2(y, dpars$mu, dpars$mu * (1 + dpars$mu),
                              log = TRUE))
      }
      RTMB::dnbinom_robust(y, lmu, 2 * lmu, log = TRUE)
    },
    valid_y = count_y("geometric"),
    init_dpars = list(mu = function(y, aterms) mean(y) + 0.1),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu * (1 + dpars$mu),
      # negbinomial with the shape fixed at 1
      dev_fn = function(y, dpars, aterms) {
        nbinom_deviance(y, dpars$mu, 1)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rnbinom(n, size = 1, mu = dpars$mu)
    }
  )
}

#' Exponential family in the mean parameterization, single dpar `mu`
#' (the mean, that is one over the rate). It supplies a CDF and a
#' truncated mean.
#'
#' @noRd
fam_exponential <- function(link = "log") {
  frmtmb_family(
    "exponential",
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms) {
      -log(dpars$mu) - y / dpars$mu
    },
    lcdf = function(q, dpars, aterms) {
      1 - exp(-q / dpars$mu)
    },
    valid_y = positive_y("exponential"),
    init_dpars = list(mu = function(y, aterms) mean(y)),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu^2,
      # Gamma with the shape fixed at 1
      dev_fn = function(y, dpars, aterms) gamma_deviance(y, dpars$mu),
      # int_lb^ub y f(y) dy = (lb + mu) e^{-lb/mu} - (ub + mu) e^{-ub/mu}
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        mu <- dpars$mu
        lo <- pmax(lb, 0)
        el <- exp(-lo / mu)
        eu <- exp(-ub / mu)
        # Inf * 0 at an absent upper bound; the term is zero there
        hi_term <- ifelse(is.finite(ub), (ub + mu) * eu, 0)
        ((lo + mu) * el - hi_term) / (el - eu)
      }
    ),
    sim = function(dpars, aterms, n) stats::rexp(n, 1 / dpars$mu)
  )
}

#' Weibull family, dpars `mu` and `shape`. As in brms, `mu` is the mean
#' and the scale is `mu / gamma(1 + 1 / shape)`. It supplies a CDF and a
#' truncated mean.
#'
#' @noRd
fam_weibull <- function(link = "log") {
  frmtmb_family(
    "weibull",
    dpars = c("mu", "shape"),
    links = list(mu = link, shape = "log"),
    lpdf = function(y, dpars, aterms) {
      # brms parameterization: mu is the mean, scale = mu/gamma(1+1/k)
      sc <- dpars$mu / exp(lgamma(1 + 1 / dpars$shape))
      RTMB::dweibull(y, shape = dpars$shape, scale = sc, log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      sc <- dpars$mu / exp(lgamma(1 + 1 / dpars$shape))
      1 - exp(-(q / sc)^dpars$shape)
    },
    valid_y = positive_y("weibull"),
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      shape = function(y, aterms) 1.2
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) {
        g1 <- exp(lgamma(1 + 1 / dpars$shape))
        g2 <- exp(lgamma(1 + 2 / dpars$shape))
        dpars$mu^2 * (g2 / g1^2 - 1)
      },
      # int_lb^ub y f(y) dy = scale * gamma(1 + 1/k) * (P(1 + 1/k, zu) -
      # P(1 + 1/k, zl)) with z = (y/scale)^k; scale * gamma(1 + 1/k) = mu
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        k <- dpars$shape
        sc <- dpars$mu / exp(lgamma(1 + 1 / k))
        zl <- (pmax(lb, 0) / sc)^k
        zu <- (ub / sc)^k
        dpars$mu * (stats::pgamma(zu, 1 + 1 / k) -
                      stats::pgamma(zl, 1 + 1 / k)) /
          (exp(-zl) - exp(-zu))
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rweibull(n, shape = dpars$shape,
                      scale = dpars$mu / gamma(1 + 1 / dpars$shape))
    }
  )
}

#' Shifted lognormal for response times, dpars `mu`, `sigma` (both on
#' the log scale) and `ndt`, the non-decision time the distribution
#' starts at.
#'
#' @noRd
fam_shifted_lognormal <- function(link = "identity") {
  frmtmb_family(
    "shifted_lognormal",
    dpars = c("mu", "sigma", "ndt"),
    links = list(mu = link, sigma = "log", ndt = "log"),
    lpdf = function(y, dpars, aterms) {
      # y <= ndt gives NaN, which the optimizer treats as a rejected
      # step; the ndt init keeps the start feasible
      RTMB::dnorm(log(y - dpars$ndt), dpars$mu, dpars$sigma, log = TRUE) -
        log(y - dpars$ndt)
    },
    valid_y = positive_y("shifted_lognormal"),
    init_dpars = list(
      mu = function(y, aterms) mean(log(y - min(y) / 2)),
      sigma = function(y, aterms) stats::sd(log(y - min(y) / 2)),
      ndt = function(y, aterms) min(y) / 2
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        dpars$ndt + exp(dpars$mu + dpars$sigma^2 / 2)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars$ndt + stats::rlnorm(n, dpars$mu, dpars$sigma)
    }
  )
}

#' Hurdle gamma for a non-negative response, dpars `mu`, `shape` and
#' `hu`. `hu` is the probability of an exact zero and the positive part
#' is a mean-parameterized gamma.
#'
#' @noRd
fam_hurdle_gamma <- function(link = "log") {
  frmtmb_family(
    "hurdle_gamma",
    dpars = c("mu", "shape", "hu"),
    links = list(mu = link, shape = "log", hu = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      yp <- y + i0   # dodge dgamma(0) = -Inf; the term carries weight 0
      g <- gate_logs(dpars, "hu")
      i0 * g$l +
        (1 - i0) * (g$l1m +
                      RTMB::dgamma(yp, shape = dpars$shape,
                                   scale = dpars$mu / dpars$shape,
                                   log = TRUE))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0)) {
        stop("hurdle_gamma: response must be non-negative", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) if (any(y > 0)) mean(y[y > 0]) else 1,
      shape = function(y, aterms) 1,
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars$hu) * dpars$mu
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars$hu) *
        stats::rgamma(n, shape = dpars$shape,
                      scale = dpars$mu / dpars$shape)
    }
  )
}

#' Hurdle lognormal for a non-negative response, dpars `mu`, `sigma`
#' (both on the log scale) and `hu`, the probability of an exact zero.
#'
#' @noRd
fam_hurdle_lognormal <- function(link = "identity") {
  frmtmb_family(
    "hurdle_lognormal",
    dpars = c("mu", "sigma", "hu"),
    links = list(mu = link, sigma = "log", hu = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      yp <- y + i0
      g <- gate_logs(dpars, "hu")
      i0 * g$l +
        (1 - i0) * (g$l1m +
                      RTMB::dnorm(log(yp), dpars$mu, dpars$sigma,
                                  log = TRUE) - log(yp))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0)) {
        stop("hurdle_lognormal: response must be non-negative",
             call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        if (any(y > 0)) mean(log(y[y > 0])) else 0
      },
      sigma = function(y, aterms) {
        if (sum(y > 0) > 1) stats::sd(log(y[y > 0])) else 1
      },
      hu = function(y, aterms) min(max(mean(y == 0), 0.05), 0.95)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars$hu) * exp(dpars$mu + dpars$sigma^2 / 2)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars$hu) *
        stats::rlnorm(n, dpars$mu, dpars$sigma)
    }
  )
}

#' Zero-inflated binomial, dpars `mu` (the success probability) and `zi`
#' (the probability of a structural zero). Trials come from the
#' `trials()` addition term.
#'
#' @noRd
fam_zi_binomial <- function(link = "logit") {
  lk <- get_link(link)
  frmtmb_family(
    "zero_inflated_binomial",
    dpars = c("mu", "zi"),
    links = list(mu = lk, zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      size <- aterms$trials %||% 1
      i0 <- as.numeric(y == 0)
      g <- gate_logs(dpars, "zi")
      lo <- robust_logit(dpars, lk)
      if (is.null(lo)) {
        lp0 <- size * log(1 - dpars$mu)
        base <- RTMB::dbinom(y, size, dpars$mu, log = TRUE)
      } else {
        # log P(0) = size * log(1 - mu), both terms off the log-odds
        lp0 <- size * log1m_inv_logit(lo)
        base <- RTMB::dbinom_robust(y, size, lo, log = TRUE)
      }
      i0 * RTMB::logspace_add(g$l, g$l1m + lp0) + (1 - i0) * (g$l1m + base)
    },
    valid_y = function(y, aterms) {
      size <- aterms$trials %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("zero_inflated_binomial: response must be integer counts ",
             "in [0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        size <- aterms$trials %||% 1
        min(max(mean(y / size), 0.05), 0.95)
      },
      zi = function(y, aterms) min(max(mean(y == 0) / 2, 0.05), 0.9)
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars$zi) * dpars$mu * (aterms$trials %||% 1)
      }
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars$zi) *
        stats::rbinom(n, aterms$trials %||% 1, dpars$mu)
    }
  )
}

#' Zero-inflated beta for a response in `[0, 1)`, dpars `mu`, `phi` and
#' `zi`, the probability of an exact zero.
#'
#' @noRd
fam_zi_beta <- function(link = "logit") {
  lk <- get_link(link)
  frmtmb_family(
    "zero_inflated_beta",
    dpars = c("mu", "phi", "zi"),
    links = list(mu = lk, phi = "log", zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      ya <- y + i0 * 0.5   # dodge dbeta(0) = -Inf; term carries weight 0
      g <- gate_logs(dpars, "zi")
      mp <- mu_pair(dpars, lk)
      i0 * g$l +
        (1 - i0) * (g$l1m +
                      RTMB::dbeta(ya, mp$p * dpars$phi,
                                  mp$q * dpars$phi, log = TRUE))
    },
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y >= 1)) {
        stop("zero_inflated_beta: response must be in [0, 1)",
             call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        if (any(y > 0)) min(max(mean(y[y > 0]), 0.05), 0.95) else 0.5
      },
      phi = function(y, aterms) 5,
      zi = function(y, aterms) min(max(mean(y == 0), 0.05), 0.9)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) (1 - dpars$zi) * dpars$mu
    ),
    sim = function(dpars, aterms, n) {
      stats::rbinom(n, 1, 1 - dpars$zi) *
        stats::rbeta(n, dpars$mu * dpars$phi, (1 - dpars$mu) * dpars$phi)
    }
  )
}

#' Asymmetric Laplace family, dpars `mu`, `sigma` and `quantile`. At a
#' fixed `quantile` the maximum-likelihood fit gives the quantile
#' regression point estimates.
#'
#' @noRd
fam_asym_laplace <- function(link = "identity") {
  frmtmb_family(
    "asym_laplace",
    dpars = c("mu", "sigma", "quantile"),
    links = list(mu = link, sigma = "log", quantile = "logit"),
    lpdf = function(y, dpars, aterms) {
      # rho_p(u) = 0.5 * (|u| + (2p - 1) u); ML at fixed quantile p
      # reproduces quantile-regression point estimates
      p <- dpars$quantile
      u <- (y - dpars$mu) / dpars$sigma
      # log(p) + log(1 - p) is the normalizer; an extreme quantile makes
      # one of them -Inf through the logit round trip
      gq <- gate_logs(dpars, "quantile")
      gq$l + gq$l1m - log(dpars$sigma) -
        0.5 * (abs(u) + (2 * p - 1) * u)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y),
      sigma = function(y, aterms) stats::sd(y) / 2,
      quantile = function(y, aterms) 0.5
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        dpars$mu + dpars$sigma * (1 - 2 * dpars$quantile) /
          (dpars$quantile * (1 - dpars$quantile))
      }
    ),
    sim = function(dpars, aterms, n) {
      p <- dpars$quantile
      dpars$mu + dpars$sigma *
        (stats::rexp(n) / p - stats::rexp(n) / (1 - p))
    }
  )
}

#' Zero-inflated asymmetric Laplace (brms spelling): a point mass at
#' zero mixed with the continuous ALD, so the mixture density is
#' well-defined without the +i0 dodge the discrete zi families need.
#'
#' @noRd
fam_zi_asym_laplace <- function(link = "identity") {
  frmtmb_family(
    "zero_inflated_asym_laplace",
    dpars = c("mu", "sigma", "quantile", "zi"),
    links = list(mu = link, sigma = "log", quantile = "logit",
                 zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      p <- dpars$quantile
      u <- (y - dpars$mu) / dpars$sigma
      gq <- gate_logs(dpars, "quantile")
      ald <- gq$l + gq$l1m - log(dpars$sigma) -
        0.5 * (abs(u) + (2 * p - 1) * u)
      g <- gate_logs(dpars, "zi")
      i0 * g$l + (1 - i0) * (g$l1m + ald)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y[y != 0]),
      sigma = function(y, aterms) {
        s <- stats::sd(y[y != 0]) / 2
        if (is.finite(s) && s > 0) s else 1
      },
      quantile = function(y, aterms) 0.5,
      zi = function(y, aterms) min(max(mean(y == 0), 0.05), 0.9)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        (1 - dpars$zi) *
          (dpars$mu + dpars$sigma * (1 - 2 * dpars$quantile) /
             (dpars$quantile * (1 - dpars$quantile)))
      }
    ),
    sim = function(dpars, aterms, n) {
      p <- dpars$quantile
      stats::rbinom(n, 1, 1 - dpars$zi) *
        (dpars$mu + dpars$sigma *
           (stats::rexp(n) / p - stats::rexp(n) / (1 - p)))
    }
  )
}

#' Huber's rho at tuning constant `k`, branch-free:
#' `rho(u) = (u^2 - max(|u| - k, 0)^2) / 2`. Inside the kink the second
#' term is zero and this is the gaussian `u^2 / 2`; outside it the
#' quadratics cancel and leave `k|u| - k^2/2`, the Laplace tail. `max(x,
#' 0)` is written `(x + |x|) / 2`, so nothing branches on a parameter
#' value - the switch between the two regimes has to happen ON the tape,
#' once per observation per gradient evaluation, and an `if` there would
#' bake one regime into the tape for good.
#'
#' @noRd
huber_rho <- function(u, k) {
  a <- abs(u)
  e <- a - k
  0.5 * (a * a - 0.25 * (e + abs(e))^2)
}

#' The normalizing constant of `exp(-rho_k(u))`:
#' `sqrt(2 pi) (2 Phi(k) - 1) + (2 / k) exp(-k^2 / 2)`, the gaussian
#' piece over `[-k, k]` plus the two exponential tails. Written through
#' `dnorm(k)` because `exp(-k^2 / 2) = sqrt(2 pi) phi(k)`; the tests
#' check it against a numeric integral at 1e-10.
#'
#' @noRd
huber_norm <- function(k) {
  sqrt(2 * pi) * (2 * stats::pnorm(k) - 1 + 2 * stats::dnorm(k) / k)
}

#' `E[U^2]` under the standard (sigma = 1) Huber density, for the
#' pearson-residual variance. Both integrals are closed form: the
#' truncated gaussian second moment plus the two exponential tails'.
#'
#' @noRd
huber_var_u <- function(k) {
  ph <- stats::dnorm(k)
  c0 <- 2 * stats::pnorm(k) - 1
  (c0 + 2 * ph * (2 / k + 2 / k^3)) / (c0 + 2 * ph / k)
}

#' Draws from the standard Huber density by inverting its CDF. This runs
#' off the tape, so the three pieces (exponential lower tail, gaussian
#' middle, exponential upper tail) branch freely.
#'
#' @noRd
rhuber_u <- function(n, k) {
  z <- huber_norm(k)
  p0 <- exp(-k^2 / 2) / (k * z)            # F(-k), and 1 - F(k)
  p <- stats::runif(n)
  u <- numeric(n)
  lo <- p < p0
  hi <- p > 1 - p0
  mid <- !lo & !hi
  u[lo] <- (log(p[lo] * k * z) - k^2 / 2) / k
  u[hi] <- (k^2 / 2 - log((1 - p[hi]) * k * z)) / k
  u[mid] <- stats::qnorm(stats::pnorm(-k) +
                           (p[mid] - p0) * z / sqrt(2 * pi))
  u
}

#' Huber's least-favorable density, dpars `mu` and `sigma`. `k` is a
#' fixed tuning constant of the FAMILY, not a dpar: it says where the
#' analyst draws the line between "residual" and "outlier", which is a
#' choice, not a quantity the data identifies. `MASS::rlm()` treats it
#' the same way. Estimating it would let the likelihood buy fit by
#' widening the gaussian core, which is the opposite of what the family
#' is for.
#'
#' @noRd
fam_huber <- function(link = "identity", k = 1.345) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("huber(k =): the tuning constant must be one finite positive ",
         "number; 1.345 (the default, and MASS::rlm()'s) gives 95% ",
         "efficiency against a gaussian", call. = FALSE)
  }
  lognorm <- log(huber_norm(k))
  varu <- huber_var_u(k)
  frmtmb_family(
    "huber",
    dpars = c("mu", "sigma"),
    links = list(mu = link, sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      s <- resid_sd(dpars$sigma, aterms)
      -log(s) - lognorm - huber_rho((y - dpars$mu) / s, k)
    },
    init_dpars = list(
      mu = function(y, aterms) stats::median(y),
      sigma = function(y, aterms) {
        s <- stats::mad(y)
        if (is.finite(s) && s > 0) s else max(stats::sd(y), 1e-3)
      }
    ),
    type = "continuous",
    post = list(
      # symmetric about mu, so the mean is mu
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) {
        resid_sd(dpars$sigma, aterms)^2 * varu
      },
      # 2 (ll at mu = y - ll at mu), the scale held fixed: 2 sigma^2
      # rho_k(u), which collapses to the gaussian (y - mu)^2 as k grows
      dev_fn = function(y, dpars, aterms) {
        s <- resid_sd(dpars$sigma, aterms)
        2 * s^2 * huber_rho((y - dpars$mu) / s, k)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars$mu + resid_sd(dpars$sigma, aterms) * rhuber_u(n, k)
    }
  )
}

# --- RTMBdist-backed families ---

#' Beta-binomial family in the mean-precision parameterization, dpars
#' `mu` and `phi`. Trials come from the `trials()` addition term.
#'
#' @noRd
fam_beta_binomial <- function(link = "logit") {
  lk <- get_link(link)
  frmtmb_family(
    "beta_binomial",
    dpars = c("mu", "phi"),
    links = list(mu = lk, phi = "log"),
    lpdf = function(y, dpars, aterms) {
      size <- aterms$trials %||% 1
      # RTMBdist has no log-odds form, so the robustness has to go into
      # the shapes: a zero second shape gives lgamma(0) = NaN
      mp <- mu_pair(dpars, lk)
      RTMBdist::dbetabinom(y, size, mp$p * dpars$phi, mp$q * dpars$phi,
                           log = TRUE)
    },
    valid_y = function(y, aterms) {
      size <- aterms$trials %||% 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("beta_binomial: response must be integer counts in ",
             "[0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        p <- mean(y / (aterms$trials %||% 1))
        min(max(p, 0.02), 0.98)
      },
      phi = function(y, aterms) 5
    ),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu * (aterms$trials %||% 1)
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rbetabinom(n, aterms$trials %||% 1,
                           dpars$mu * dpars$phi,
                           (1 - dpars$mu) * dpars$phi)
    }
  )
}

#' Skew-normal family in the mean parameterization, dpars `mu` (the
#' mean), `sigma` and `alpha` (the skewness). The `alpha` start avoids
#' zero, which is a stationary point of the likelihood.
#'
#' @noRd
fam_skew_normal <- function(link = "identity") {
  frmtmb_family(
    "skew_normal",
    dpars = c("mu", "sigma", "alpha"),
    links = list(mu = link, sigma = "log", alpha = "identity"),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dskewnorm2(y, dpars$mu, dpars$sigma, dpars$alpha,
                           log = TRUE)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y),
      alpha = function(y, aterms) {
        # alpha = 0 is a stationary point of the skew-normal likelihood
        # (singular information); start from the sample skewness side
        m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
        2 * sign(m3) + 0.5 * m3
      }
    ),
    type = "continuous",
    post = list(mean_fn = function(dpars, aterms) dpars$mu),
    sim = function(dpars, aterms, n) {
      RTMBdist::rskewnorm2(n, dpars$mu, dpars$sigma, dpars$alpha)
    }
  )
}

#' Inverse gaussian family, dpars `mu` (the mean) and `shape`. It
#' supplies a CDF and a truncated mean.
#'
#' @noRd
fam_inverse_gaussian <- function(link = "log") {
  frmtmb_family(
    "inverse.gaussian",
    dpars = c("mu", "shape"),
    links = list(mu = link, shape = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dinvgauss(y, mean = dpars$mu, shape = dpars$shape,
                          log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMBdist::pinvgauss(q, mean = dpars$mu, shape = dpars$shape)
    },
    valid_y = positive_y("inverse.gaussian"),
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      shape = function(y, aterms) mean(y)^3 / max(stats::var(y), 1e-8)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$mu^3 / dpars$shape,
      # 1 / shape is the dispersion and divides out
      dev_fn = function(y, dpars, aterms) {
        (y - dpars$mu)^2 / (y * dpars$mu^2)
      },
      # Partial moment (Jorgensen): E[Y 1{Y <= x}] =
      # mu (Phi(z1) - e^{2 lambda / mu} Phi(-z2)), with the CDF the same
      # pair added instead of subtracted; the exponential factor
      # overflows on its own, so it rides in log space
      trunc_mean_fn = function(dpars, aterms, lb, ub) {
        mu <- dpars$mu
        lam <- dpars$shape
        part <- function(x) {
          r <- sqrt(lam / x)
          mu * (stats::pnorm(r * (x / mu - 1)) -
                  exp(2 * lam / mu +
                        stats::pnorm(-r * (x / mu + 1), log.p = TRUE)))
        }
        cdf <- function(x) {
          r <- sqrt(lam / x)
          stats::pnorm(r * (x / mu - 1)) +
            exp(2 * lam / mu +
                  stats::pnorm(-r * (x / mu + 1), log.p = TRUE))
        }
        lo <- pmax(lb, 0)
        pl <- ifelse(lo > 0, part(lo), 0)
        fl <- ifelse(lo > 0, cdf(lo), 0)
        pu <- ifelse(is.finite(ub), part(ub), mu)
        fu <- ifelse(is.finite(ub), cdf(ub), 1)
        (pu - pl) / (fu - fl)
      }
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rinvgauss(n, mean = dpars$mu, shape = dpars$shape)
    }
  )
}

#' Ex-gaussian family, dpars `mu`, `sigma` and `beta`.
#' brms parameterization: `mu` is the DISTRIBUTION mean, `beta` the
#' scale of the exponential component (gaussian component sits at
#' `mu - beta`).
#'
#' @noRd
fam_exgaussian <- function(link = "identity") {
  frmtmb_family(
    "exgaussian",
    dpars = c("mu", "sigma", "beta"),
    links = list(mu = link, sigma = "log", beta = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dexgauss(y, dpars$mu - dpars$beta, dpars$sigma,
                         1 / dpars$beta, log = TRUE)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y) / 2,
      beta = function(y, aterms) stats::sd(y) / 2
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$sigma^2 + dpars$beta^2
    ),
    sim = function(dpars, aterms, n) {
      RTMBdist::rexgauss(n, dpars$mu - dpars$beta, dpars$sigma,
                         1 / dpars$beta)
    }
  )
}

# Exact 1{y == k} for k = 1..K built from arithmetic alone. RTMB
# advectors carry no comparison operators, and oneStepPredict re-tapes
# the objective with the response promoted to a parameter, which is what
# breaks the indexing forms below. At integer y this Lagrange basis is
# exact in floating point - off the diagonal one factor is exactly 0, on
# it every factor is exactly 1 - so the taped path and the data path
# agree bit for bit.
#' During oneStepPredict RTMB hands the response to the lpdf as an "osa"
#' object: the taped value in `@x` and the per-row data-term indicator
#' in `@keep`. RTMB's own densities apply the indicator through
#' dGenericOSA, so a hand-written lpdf has to do it itself or every
#' observation stays switched on and the one-step sequence collapses.
#' This returns the pair, or `NULL` when the response is plain data.
#'
#' @noRd
osa_unwrap <- function(y) {
  if (!methods::is(y, "osa")) return(NULL)
  keep <- y@keep
  if (ncol(keep) != 1L) {
    stop("osa_method = \"cdf\" is not supported for this family",
         call. = FALSE)
  }
  list(y = y@x, keep = keep[, 1])
}

#' The `1{y == k}` indicators for `k = 1..K`, one vector per category,
#' from the arithmetic Lagrange basis described above. It works with the
#' response on the AD tape, where comparison operators are not available.
#'
#' @noRd
ord_cat_sel <- function(y, K) {
  lapply(seq_len(K), function(k) {
    s <- 1
    for (j in seq_len(K)) {
      if (j != k) s <- s * ((y - j) / (k - j))
    }
    s
  })
}

#' `1{j < y}` for `j = 1..K-1`, from the same basis.
#'
#' @noRd
ord_cat_below <- function(sel, K) {
  lapply(seq_len(K - 1L), function(j) {
    b <- sel[[j + 1L]]
    if (j + 1L < K) for (k in (j + 2L):K) b <- b + sel[[k]]
    b
  })
}

#' Cumulative ordinal: response in `1..K` (or an ordered factor). The
#' linear predictor has no intercept; K-1 ordered thresholds take its
#' place, parameterized as (tau_1, log increments) in `extra_pars`.
#'
#' @noRd
fam_cumulative <- function(link = "logit") {
  Fcdf <- switch(link,
    logit = function(x) 1 / (1 + exp(-x)),
    probit = function(x) RTMB::pnorm(x),
    stop("cumulative() supports links 'logit' and 'probit'", call. = FALSE)
  )
  # the logistic CDF has an exact log-space difference (see below); the
  # normal one does not, because RTMB::pnorm carries no log.p
  robust <- identical(link, "logit")
  frmtmb_family(
    "cumulative",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      "[<-" <- RTMB::ADoverload("[<-")
      raw <- extra$tau_raw
      K1 <- length(raw)
      tau <- rep(raw[1], K1)
      if (K1 > 1) {
        for (k in 2:K1) tau[k] <- tau[k - 1] + exp(raw[k])
      }
      eta <- dpars$mu
      K <- K1 + 1L
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        # OSA re-tape: pick the category probability arithmetically
        sel <- ord_cat_sel(ov$y, K)
        Fk <- lapply(seq_len(K1), function(k) Fcdf(tau[k] - eta))
        dens <- sel[[1]] * Fk[[1]]
        if (K1 > 1) {
          for (k in 2:K1) dens <- dens + sel[[k]] * (Fk[[k]] - Fk[[k - 1]])
        }
        return(log(dens + sel[[K]] * (1 - Fk[[K1]])) * ov$keep)
      }
      iK <- as.numeric(y == K)
      i1 <- as.numeric(y == 1)
      if (!robust) {
        up <- Fcdf(tau[pmin(y, K1)] - eta) * (1 - iK) + iK
        lo <- Fcdf(tau[pmax(y - 1, 1)] - eta) * (1 - i1)
        return(log(up - lo))
      }
      # log(F(a) - F(b)) for the logistic F, entirely in log space:
      # F(a) - F(b) = (e^-b - e^-a) / ((1 + e^-a)(1 + e^-b)), and the
      # numerator's logspace_sub is well conditioned because its two
      # arguments differ by tau_k - tau_{k-1}, which does not move with
      # eta. The difference of two saturated CDFs loses every digit
      # from |eta| = 20 and is exactly 0 by 40.
      out <- i1 * log_inv_logit(tau[1] - eta) +
        iK * log1m_inv_logit(tau[K1] - eta)
      if (K1 >= 2L) {
        # data-only clamp into the interior categories, so the masked
        # rows still evaluate a legal (strictly ordered) threshold pair
        ym <- pmin(pmax(y, 2L), K1)
        a <- tau[ym] - eta
        b <- tau[ym - 1L] - eta
        out <- out + (1 - i1 - iK) *
          (RTMB::logspace_sub(-b, -a) + log_inv_logit(a) +
             log_inv_logit(b))
      }
      out
    },
    valid_y = function(y, aterms) {
      if (any(y < 1) || any(y != round(y)) || length(unique(y)) < 2) {
        stop("cumulative: response must be an ordered factor or integers ",
             "1..K with at least 2 observed categories", call. = FALSE)
      }
    },
    type = "ordinal",
    extra_pars = function(y, aterms) {
      K <- max(y)
      p <- cumsum(tabulate(y, K) / length(y))[-K]
      p <- pmin(pmax(p, 0.01), 0.99)
      tau0 <- stats::qlogis(p)
      incr <- pmax(diff(tau0), 0.05)
      list(tau_raw = c(tau0[1], log(incr)))
    },
    sim = ord_sim("cumulative", ordered = TRUE, link = link),
    drop_intercept = TRUE
  )
}

#' Shared scaffolding for the sequential ordinal families: an
#' n x (K-1) matrix of `(tau_j - eta_i)` or `(eta_i - tau_j)`, and
#' data-only indicator matrices selecting the observed category
#' (branch-free).
#'
#' @noRd
ord_indicators <- function(y, K1) {
  n <- length(y)
  jj <- rep(seq_len(K1), each = n)
  yy <- rep(y, K1)
  list(
    sel = matrix(as.numeric(yy == jj), n, K1),   # j == y (y <= K-1)
    below = matrix(as.numeric(jj < yy), n, K1)   # j < y
  )
}

#' The n x (K-1) matrix of `tau_j - eta_i` the ordinal lpdfs work over.
#' It broadcasts by matrix multiplication because `rep()` would strip the
#' advector class.
#'
#' @noRd
ord_eta_mat <- function(eta, tau, n, K1) {
  # broadcast via matmul: rep() strips the advector class
  TM <- RTMB::matrix(1, n, 1) %*% RTMB::matrix(tau, 1, K1)
  TM - eta   # column-wise recycling: row i is tau_j - eta_i
}

#' Sequential (sratio/cratio) log-density with the response on the tape.
#' Column-at-a-time so nothing indexes an advector and nothing needs the
#' n x (K-1) matrices, which the data path builds for speed.
#'
#' @noRd
ord_seq_lpdf_ad <- function(y, eta, tau, K1, cs, Fcdf, stopping) {
  sel <- ord_cat_sel(y, K1 + 1L)
  below <- ord_cat_below(sel, K1 + 1L)
  out <- 0
  for (j in seq_len(K1)) {
    Mj <- tau[j] - eta
    if (!is.null(cs)) Mj <- Mj - cs[, j]
    Pj <- if (stopping) Fcdf(Mj) else 1 - Fcdf(-Mj)
    out <- out + sel[[j]] * log(Pj) + below[[j]] * log(1 - Pj)
  }
  out
}

#' Make the response validator the ordinal families share: the response
#' must be integers `1..K`, and at least two categories must be
#' observed.
#'
#' @noRd
ord_valid_y <- function(name) {
  function(y, aterms) {
    if (any(y < 1) || any(y != round(y)) || length(unique(y)) < 2) {
      stop(name, ": response must be an ordered factor or integers ",
           "1..K with at least 2 observed categories", call. = FALSE)
    }
  }
}

#' Starting values for the K-1 ordinal thresholds, from the observed
#' cumulative category frequencies. With `ordered = TRUE` they come back
#' in the (first threshold, log increments) parameterization.
#'
#' @noRd
ord_tau_init <- function(y, ordered = TRUE) {
  K <- max(y)
  p <- cumsum(tabulate(y, K) / length(y))[-K]
  p <- pmin(pmax(p, 0.01), 0.99)
  tau0 <- stats::qlogis(p)
  if (!ordered) return(list(tau_raw = tau0))
  incr <- pmax(diff(tau0), 0.05)
  list(tau_raw = c(tau0[1], log(incr)))
}

# --- ordinal simulators ----------------------------------------------
# The lpdfs are taped, so they answer "how likely is this category"; a
# simulator needs the whole category distribution instead. These build
# it in plain doubles, one branch per family, matching each lpdf term
# for term.

#' cumulative and sratio store ordered thresholds as (tau_1, log
#' increments); cratio and acat store them raw. This returns the
#' thresholds themselves.
#'
#' @noRd
ord_tau_from_raw <- function(raw, ordered) {
  if (!ordered || length(raw) < 2L) return(raw)
  c(raw[1L], raw[1L] + cumsum(exp(raw[-1L])))
}

#' The plain numeric CDF for an ordinal link. The simulators run off the
#' tape, so they use this instead of the AD-safe version.
#'
#' @noRd
ord_num_cdf <- function(link) {
  switch(link, logit = stats::plogis, probit = stats::pnorm,
         stop("no numeric CDF for link '", link, "'", call. = FALSE))
}

#' n x K matrix of category probabilities, one branch per ordinal
#' family, in plain doubles.
#'
#' @noRd
ord_cat_probs <- function(family, eta, tau, cs, link) {
  n <- length(eta)
  K1 <- length(tau)
  K <- K1 + 1L
  if (identical(family, "acat")) {
    # P(y=r) proportional to exp((r-1) eta - cumsum(tau)[r]); the row
    # maximum comes out before exp() so a wide eta cannot overflow
    ct0 <- c(0, cumsum(tau))
    E <- outer(eta, seq_len(K) - 1) - matrix(ct0, n, K, byrow = TRUE)
    if (!is.null(cs)) {
      acc <- rep(0, n)
      for (r in seq.int(2L, K)) {
        acc <- acc + cs[, r - 1L]
        E[, r] <- E[, r] + acc
      }
    }
    ex <- exp(E - apply(E, 1L, max))
    return(ex / rowSums(ex))
  }
  Fcdf <- ord_num_cdf(link)
  M <- matrix(tau, n, K1, byrow = TRUE) - eta   # tau_j - eta_i
  if (!is.null(cs)) M <- M - cs
  if (identical(family, "cumulative")) {
    # P(y=k) = F(tau_k - eta) - F(tau_{k-1} - eta), with the two
    # boundary values pinned at 0 and 1: a column-wise difference of the
    # K+1 cumulative probabilities
    Fm <- cbind(0, Fcdf(M), 1)
    return(Fm[, -1L, drop = FALSE] - Fm[, -ncol(Fm), drop = FALSE])
  }
  # sequential families: h_j is the chance of stopping at category j
  # given survival past j-1, on each family's own scale
  h <- if (identical(family, "sratio")) Fcdf(M) else 1 - Fcdf(-M)
  P <- matrix(0, n, K)
  surv <- rep(1, n)
  for (j in seq_len(K1)) {
    P[, j] <- surv * h[, j]
    surv <- surv * (1 - h[, j])
  }
  P[, K] <- surv
  P
}

#' Make the simulator for an ordinal family. It draws one category per
#' row by inverse-CDF sampling of the full category distribution, which
#' the taped log-density does not give.
#'
#' @noRd
ord_sim <- function(family, ordered, link) {
  function(dpars, aterms, n, extra) {
    tau <- ord_tau_from_raw(extra$tau_raw, ordered)
    P <- ord_cat_probs(family, rep(dpars$mu, length.out = n), tau,
                       dpars$.cs, link)
    K <- ncol(P)
    # inverse-CDF sampling, one uniform per row
    cp <- t(apply(P, 1L, cumsum))
    if (n == 1L) cp <- matrix(cp, 1L, K)
    pmin(1L + rowSums(cp < stats::runif(n)), K)
  }
}

#' The AD-safe CDF for an ordinal link, for use inside a taped
#' log-density. An unsupported link errors with the family name.
#'
#' @noRd
ord_link_cdf <- function(name, link) {
  switch(link,
    logit = function(x) 1 / (1 + exp(-x)),
    probit = function(x) RTMB::pnorm(x),
    stop(name, "() supports links 'logit' and 'probit'", call. = FALSE)
  )
}

#' The sequential families' log-density from log-space hazards. `lstop`
#' gives `log P(stop at j)` and `lgo` gives `log P(continue past j)`,
#' both as functions of the `tau_j - eta` matrix column; sratio and
#' cratio differ only in which way round they are. One column at a time
#' rather than one matrix operation, because `logspace_add()` works on
#' the flat vector and would drop the dimensions the indicator matrices
#' need.
#'
#' @noRd
ord_log_hazard_sum <- function(M, ind, K1, lstop, lgo) {
  out <- 0
  for (j in seq_len(K1)) {
    Mj <- M[, j]
    out <- out + ind$sel[, j] * lstop(Mj) + ind$below[, j] * lgo(Mj)
  }
  out
}

#' Stopping ratio (brms sratio): `P(y=k) = F(tau_k - eta) *
#' prod_{j<k} (1 - F(tau_j - eta))`; ordered thresholds like cumulative.
#'
#' @noRd
fam_sratio <- function(link = "logit") {
  Fcdf <- ord_link_cdf("sratio", link)
  robust <- identical(link, "logit")
  frmtmb_family(
    "sratio",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      "[<-" <- RTMB::ADoverload("[<-")
      raw <- extra$tau_raw
      K1 <- length(raw)
      tau <- rep(raw[1], K1)
      if (K1 > 1) for (k in 2:K1) tau[k] <- tau[k - 1] + exp(raw[k])
      n <- length(y)
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        return(ord_seq_lpdf_ad(ov$y, dpars$mu, tau, K1, dpars$.cs, Fcdf,
                               stopping = TRUE) * ov$keep)
      }
      M <- ord_eta_mat(dpars$mu, tau, n, K1)
      if (!is.null(dpars$.cs)) M <- M - dpars$.cs
      ind <- ord_indicators(y, K1)
      if (robust) {
        # log F and log(1 - F) straight out of logspace_add: the naive
        # pair is -Inf on whichever side the logistic CDF saturated,
        # which for sratio is the negative eta tail
        return(ord_log_hazard_sum(M, ind, K1, log_inv_logit,
                                  log1m_inv_logit))
      }
      P <- Fcdf(M)
      ones <- rep(1, K1)   # rowSums strips the advector class
      as.vector((log(P) * ind$sel) %*% ones) +
        as.vector((log(1 - P) * ind$below) %*% ones)
    },
    valid_y = ord_valid_y("sratio"),
    type = "ordinal",
    extra_pars = function(y, aterms) ord_tau_init(y, ordered = TRUE),
    sim = ord_sim("sratio", ordered = TRUE, link = link),
    drop_intercept = TRUE
  )
}

#' Continuation ratio (brms cratio): `P(y=k) = (1 - F(eta - tau_k)) *
#' prod_{j<k} F(eta - tau_j)`; unordered thresholds.
#'
#' @noRd
fam_cratio <- function(link = "logit") {
  Fcdf <- ord_link_cdf("cratio", link)
  robust <- identical(link, "logit")
  frmtmb_family(
    "cratio",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      tau <- extra$tau_raw
      K1 <- length(tau)
      n <- length(y)
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        return(ord_seq_lpdf_ad(ov$y, dpars$mu, tau, K1, dpars$.cs, Fcdf,
                               stopping = FALSE) * ov$keep)
      }
      M <- ord_eta_mat(dpars$mu, tau, n, K1)   # tau_j - eta
      if (!is.null(dpars$.cs)) M <- M - dpars$.cs
      ind <- ord_indicators(y, K1)
      if (robust) {
        # P = F(-M), so log(1 - P) = log F(M) and log P = log(1 - F(M));
        # cratio saturates on the positive eta tail, the mirror of
        # sratio
        return(ord_log_hazard_sum(M, ind, K1, log_inv_logit,
                                  log1m_inv_logit))
      }
      P <- Fcdf(-M)                            # F(eta + cs_j - tau_j)
      ones <- rep(1, K1)   # rowSums strips the advector class
      as.vector((log(1 - P) * ind$sel) %*% ones) +
        as.vector((log(P) * ind$below) %*% ones)
    },
    valid_y = ord_valid_y("cratio"),
    type = "ordinal",
    extra_pars = function(y, aterms) ord_tau_init(y, ordered = FALSE),
    sim = ord_sim("cratio", ordered = FALSE, link = link),
    drop_intercept = TRUE
  )
}

#' Adjacent category (brms acat, logit link): `P(y=k)` proportional to
#' `exp(sum_{j<k} (eta - tau_j))`; unordered thresholds.
#'
#' @noRd
fam_acat <- function(link = "logit") {
  if (!identical(link, "logit")) {
    stop("acat() supports the 'logit' link only", call. = FALSE)
  }
  frmtmb_family(
    "acat",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms, extra) {
      tau <- extra$tau_raw
      K <- length(tau) + 1L
      n <- length(y)
      eta <- dpars$mu
      "c" <- RTMB::ADoverload("c")
      ct0 <- c(0, cumsum(tau))                 # length K
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        sel <- ord_cat_sel(ov$y, K)
        acc <- 0 * eta
        num <- 0
        den <- NULL
        for (r in seq_len(K)) {
          Er <- (r - 1) * eta - ct0[r]
          if (!is.null(dpars$.cs) && r >= 2L) {
            acc <- acc + dpars$.cs[, r - 1L]
            Er <- Er + acc
          }
          num <- num + sel[[r]] * Er
          # logsumexp fold: the top category's exponent is (K-1) * eta,
          # so a plain sum of exp() overflows at eta = 709 / (K - 1)
          den <- if (is.null(den)) Er else RTMB::logspace_add(den, Er)
        }
        return((num - den) * ov$keep)
      }
      # E[i, r] = (r-1) * eta_i - cumsum tau, r = 1..K; broadcast by
      # matmul (rep() strips the advector class)
      Rm <- matrix(rep(seq_len(K) - 1L, each = n), n, K)
      E <- Rm * eta -
        RTMB::matrix(1, n, 1) %*% RTMB::matrix(ct0, 1, K)
      if (!is.null(dpars$.cs)) {
        "[<-" <- RTMB::ADoverload("[<-")
        CS <- dpars$.cs
        acc <- 0 * eta
        for (r in seq.int(2L, K)) {
          acc <- acc + CS[, r - 1L]
          E[, r] <- E[, r] + acc
        }
      }
      jj <- rep(seq_len(K), each = n)
      S <- matrix(as.numeric(rep(y, K) == jj), n, K)
      ones <- rep(1, K)   # rowSums strips the advector class
      den <- E[, 1L]
      for (r in seq_len(K)[-1L]) den <- RTMB::logspace_add(den, E[, r])
      as.vector((E * S) %*% ones) - den
    },
    valid_y = ord_valid_y("acat"),
    type = "ordinal",
    extra_pars = function(y, aterms) ord_tau_init(y, ordered = FALSE),
    sim = ord_sim("acat", ordered = FALSE, link = link),
    drop_intercept = TRUE
  )
}

#' Finite mixture families
#'
#' `mixture(fam1, fam2, ...)` builds a K-component mixture: each
#' component keeps its own distributional parameters, suffixed by the
#' component index (`mu1`, `sigma1`, `mu2`, ...), and the mixing
#' proportions come from `theta1 ... theta{K-1}` (multinomial-logit
#' against the last component, each with its own linear predictor - so
#' mixing weights may depend on covariates). The main model formula
#' applies to every component mean; override per component with
#' `bf(y ~ x, mu2 ~ 1)`.
#'
#' The likelihood is a parameter-branch-free logsumexp, so Laplace
#' machinery is untouched; the usual finite-mixture ML caveats apply
#' instead: the likelihood is invariant to component relabeling, so the
#' component means are initialized on spread-out response quantiles,
#' and multimodality is real (compare starts, or order the intercepts
#' through `lower`/`upper`). Component families with extra parameters
#' (ordinal) are not supported.
#'
#' With `groups = ~g` the mixture moves to the group level (latent
#' classes): every observation of a group shares one class draw, and
#' the marginal likelihood sums the class assignment per group.
#' Continuous random effects, smooths, and gp() terms are allowed in
#' the component formulas - the class sum happens conditional on the
#' latent effects, so one Laplace approximation integrates them
#' (growth-mixture models). Random effects written in a component
#' formula are class-specific by construction; the Laplace
#' approximation of the class-mixture integrand is not exact even for
#' gaussian responses (a fraction of a log-likelihood unit in typical
#' well-separated problems). `quadrature = TRUE` makes the integral
#' numerically exact when the per-group integrand is univariate (one
#' scalar random intercept, in one class); with class-specific
#' intercepts in several classes the coordinates couple and quadrature
#' remains approximate - use `frmtmb.sample::check_laplace()` to judge.
#' Mixing-weight predictors are evaluated at each group's first row
#' (use group-constant covariates). [mixture_probs()] returns the
#' posterior class probabilities per group (or per observation for
#' ordinary mixtures), conditional on the random-effect modes.
#'
#' @param ... Two or more component families.
#' @param groups Optional one-sided formula naming the latent-class
#'   grouping factor.
#' @return A `frmtmb_family`.
#' @examples
#' # two well-separated gaussian components
#' set.seed(3)
#' dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)),
#'                  x = rnorm(160))
#' fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
#' # one mu and sigma per component, plus the mixing weight theta1
#' fixef(fit)
#' # posterior class probability per observation
#' head(mixture_probs(fit))
#'
#' # the mixing weight can take its own predictor
#' frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()), data = dd)
#'
#' \donttest{
#' # latent classes: every observation of a group shares one class
#' set.seed(4)
#' n_g <- 40
#' cls <- rep(c(1, 2), each = n_g / 2)
#' dg <- data.frame(g = factor(rep(seq_len(n_g), each = 5)))
#' dg$y <- rnorm(nrow(dg), c(0, 4)[cls[as.integer(dg$g)]], 1)
#' fg <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
#'           data = dg)
#' head(mixture_probs(fg))   # one row per group, not per observation
#' }
#' @export
mixture <- function(..., groups = NULL) {
  comps <- lapply(list(...), as_frmtmb_family)
  K <- length(comps)
  if (K < 2L) {
    stop("mixture() needs at least two component families", call. = FALSE)
  }
  for (cp in comps) {
    if (!is.null(cp$extra_pars) || isTRUE(cp$drop_intercept)) {
      stop("mixture() does not support component family '", cp$family,
           "'", call. = FALSE)
    }
    if (!"mu" %in% cp$dpars) {
      stop("mixture() components need a 'mu' parameter", call. = FALSE)
    }
  }

  dpars <- character(0)
  links <- list()
  for (k in seq_len(K)) {
    for (dp in comps[[k]]$dpars) {
      nm <- paste0(dp, k)
      dpars <- c(dpars, nm)
      links[[nm]] <- comps[[k]]$links[[dp]]
    }
  }
  for (k in seq_len(K - 1L)) {
    nm <- paste0("theta", k)
    dpars <- c(dpars, nm)
    links[[nm]] <- "identity"
  }

  comp_dpars <- function(dpars_all, k) {
    stats::setNames(
      lapply(comps[[k]]$dpars, function(dp) dpars_all[[paste0(dp, k)]]),
      comps[[k]]$dpars
    )
  }
  # log mixing weights: multinomial logit, last component reference
  log_pi <- function(dpars_all) {
    Ts <- lapply(seq_len(K - 1L), function(k) {
      dpars_all[[paste0("theta", k)]]
    })
    Ts[[K]] <- 0 * dpars_all$mu1
    lse <- Ts[[1L]]
    for (k in seq.int(2L, K)) lse <- RTMB::logspace_add(lse, Ts[[k]])
    lapply(Ts, function(t_) t_ - lse)
  }

  init <- list()
  for (k in seq_len(K)) {
    kk <- k
    init[[paste0("mu", k)]] <- local({
      k_ <- kk
      function(y, aterms) stats::quantile(y, k_ / (K + 1), names = FALSE)
    })
    for (dp in setdiff(comps[[k]]$dpars, "mu")) {
      fn <- comps[[k]]$init_dpars[[dp]]
      if (!is.null(fn)) init[[paste0(dp, k)]] <- fn
    }
  }

  types <- unique(vapply(comps, `[[`, "", "type"))
  fam <- frmtmb_family(
    paste0("mixture(", paste(vapply(comps, `[[`, "", "family"),
                             collapse = ", "), ")"),
    dpars = dpars,
    links = links,
    lpdf = function(y, dpars, aterms) {
      lp <- log_pi(dpars)
      ll <- NULL
      for (k in seq_len(K)) {
        llk <- comps[[k]]$lpdf(y, comp_dpars(dpars, k), aterms) + lp[[k]]
        ll <- if (is.null(ll)) llk else RTMB::logspace_add(ll, llk)
      }
      ll
    },
    valid_y = function(y, aterms) {
      for (cp in comps) {
        if (!is.null(cp$valid_y)) cp$valid_y(y, aterms)
      }
    },
    init_dpars = init,
    type = if (length(types) == 1L) types else "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        lp <- log_pi(dpars)
        out <- 0
        for (k in seq_len(K)) {
          mk <- comps[[k]]$post$mean_fn
          if (is.null(mk)) return(NULL)
          out <- out + exp(lp[[k]]) * mk(comp_dpars(dpars, k), aterms)
        }
        out
      }
    ),
    sim = function(dpars, aterms, n) {
      lp <- log_pi(dpars)
      P <- vapply(lp, function(l) rep(exp(l), length.out = n),
                  numeric(n))
      ks <- vapply(seq_len(n), function(i) {
        sample.int(K, 1L, prob = P[i, ])
      }, integer(1))
      out <- numeric(n)
      for (k in seq_len(K)) {
        sk <- comps[[k]]$sim
        if (is.null(sk)) stop("Mixture component ", k, " ('",
                              comps[[k]]$family,
                              "') has no simulator", call. = FALSE)
        idx <- which(ks == k)
        if (length(idx)) {
          dk <- lapply(comp_dpars(dpars, k), function(v) {
            rep(v, length.out = n)[idx]
          })
          out[idx] <- sk(dk, aterms, length(idx))
        }
      }
      out
    },
    primary_dpars = paste0("mu", seq_len(K))
  )
  # internals for the objective's group-level branch and for
  # mixture_probs(): per-component log-densities and log mixing weights
  fam$mix <- list(
    K = K,
    comp_lpdf = function(y, dpars, aterms, k) {
      comps[[k]]$lpdf(y, comp_dpars(dpars, k), aterms)
    },
    comp_dpars = comp_dpars,
    comp_sim = function(dpars_k, aterms, n, k) {
      sk <- comps[[k]]$sim
      if (is.null(sk)) {
        stop("Mixture component ", k, " ('", comps[[k]]$family,
             "') has no simulator for class-wise draws", call. = FALSE)
      }
      sk(dpars_k, aterms, n)
    },
    log_pi = log_pi
  )
  if (is.null(groups)) {
    # A rowwise mixture still refuses the two fitting options that
    # expand about a single inner mode, and that refusal is a property
    # of the family, so it travels with it rather than sitting in a
    # named gate the core has to run.
    fam[["structure"]] <- frmtmb_structure(
      latent_probs = function(fit, block) mixture_posterior(fit),
      supports = structure_supports_all(reml = FALSE, profile = FALSE),
      refusals = mixture_multimodal_refusals("a mixture() family")
    )
  } else {
    if (!inherits(groups, "formula") || length(groups) != 2L) {
      stop("groups must be a one-sided formula: groups = ~g",
           call. = FALSE)
    }
    fam$mix_groups <- groups
    # A class belongs to the GROUP, so neither the likelihood nor the
    # draw is rowwise: the group's per-observation log-densities sum
    # BEFORE the logsumexp over classes, and the rowwise `sim` above
    # would give one group two classes. That is what a structure is.
    fam[["structure"]] <- mixture_structure(fam$mix)
  }
  fam
}

#' The two fitting options every mixture-type family refuses, in its own
#' name.
#'
#' A mixture likelihood is invariant to permuting its components, so the
#' location coefficients enter a multimodal objective. Both `REML` and
#' `profile = TRUE` integrate those coefficients out with a Laplace
#' approximation about a single inner mode, which is not defined here:
#' the inner Newton solve walks between the component modes and the fit
#' either dies at "NA/NaN gradient evaluation" or reports an optimum
#' with a gradient near 1e9. Quadrature is unaffected, because it
#' marginalizes the random effects, not the coefficients.
#'
#' `what` is the family as the sentence should name it. These used to be
#' one message listing every mixture-type family, raised from one gate
#' in fit.R that knew all three by name; each family states its own now.
#'
#' @noRd
mixture_multimodal_refusals <- function(what) {
  list(
    reml = paste0(
      "REML = TRUE cannot be combined with ", what, ": the mixture ",
      "likelihood is multimodal in the fixed effects REML integrates ",
      "out, so the restricted likelihood is not defined. Use ",
      "REML = FALSE"),
    profile = paste0(
      "frmtmb_control(profile = TRUE) cannot be combined with ", what,
      ": profiling moves the fixed effects into the inner Laplace ",
      "problem, and the mixture likelihood is multimodal in them. Use ",
      "profile = FALSE")
  )
}

#' The structured-family protocol for a group-level (latent-class)
#' mixture: `mixture(groups = ~g)`.
#'
#' This is the protocol's reference implementation, and it is short on
#' purpose. Only the likelihood, the grouping data it reads and the draw
#' are not rowwise. Everything else about a mixture IS rowwise - its
#' per-row mean is the mixing-weighted mean of the component means, one
#' row at a time - so the four capability flags that depend on that mean
#' are the ones it opts into, `fitted_mean` stays `NULL`, and the core
#' keeps using the family's own `post$mean_fn`.
#'
#' @noRd
mixture_structure <- function(mx) {
  force(mx)
  frmtmb_structure(
    frame_vars = function(fam) list(fam[["mix_groups"]][[2L]]),
    check_spec = mixture_check_spec,
    frame_block = mixture_frame_block,
    unit = "a group-level mixture (mixture(groups = ))",
    loglik = function(y, dpars, aterms, weights, block, extra) {
      lps_pi <- mx[["log_pi"]](dpars)
      total <- NULL
      for (k in seq_len(mx[["K"]])) {
        ll_k <- mx[["comp_lpdf"]](y, dpars, aterms, k)
        g_k <- as.vector(block[["Gt"]] %*% (weights * ll_k)) +
          lps_pi[[k]][block[["first"]]]
        total <- if (is.null(total)) g_k else {
          RTMB::logspace_add(total, g_k)
        }
      }
      sum(total)
    },
    latent_probs = function(fit, block) mixture_posterior(fit),
    sim_ctx = mixture_sim_groups,
    # The four the doc names are the four that follow from the per-row
    # mean being rowwise. Most of the rest are TRUE not because a
    # group-level mixture can do them but because something ELSE
    # already refuses them in words of its own, and a flag here would
    # only shadow the better message: the missing CDF states cens(),
    # trunc() and mi(). `osa` is the one refusal with no better home:
    # the group branch of the objective registers no observation
    # vector, so there is nothing to step through. `cluster_robust` is
    # allowed and then checked, because the sandwich needs every
    # mixture group to sit inside one cluster rather than needing no
    # groups at all.
    # multivariate = FALSE states what mixture_check_spec() enforces
    # upstream (a latent class belongs to the group; a second response
    # would need a joint class process): the flag must not promise what
    # a future consumer of it would then wrongly allow
    supports = structure_supports_all(reml = FALSE, profile = FALSE,
                                      osa = FALSE, multivariate = FALSE),
    refusals = mixture_multimodal_refusals("a mixture() family")
  )
}

#' The model shape a group-level mixture cannot carry: a latent class
#' belongs to the GROUP, and a second response would need a joint class
#' process across responses.
#'
#' `cens()`, `trunc()` and `mi()` are refused too, but not here. A
#' mixture density carries no CDF, so the generic addition-term guard
#' refuses the first two by naming the missing CDF, and `mi()` by naming
#' the families it needs - and those are the messages a user actually
#' meets, because a `check_spec` runs BEFORE those guards and a
#' duplicate refusal here would only shadow the more specific one.
#'
#' @noRd
mixture_check_spec <- function(resp, spec, av) {
  if (length(spec$responses) > 1L || spec$rescor) {
    stop("Group-level mixtures support univariate models",
         call. = FALSE)
  }
  invisible(NULL)
}

#' The group incidence matrix and the row indices a group-level mixture
#' folds its per-observation densities over. All data, resolved once.
#'
#' @noRd
mixture_frame_block <- function(resp, spec, av, mf, y, n) {
  gv <- factor(eval(resp$family[["mix_groups"]][[2L]], mf,
                    resp$formula_env))
  G <- Matrix::sparseMatrix(i = seq_len(n), j = as.integer(gv),
                            x = 1, dims = c(n, nlevels(gv)))
  list(G = G, Gt = Matrix::t(G), first = match(levels(gv), gv),
       gindex = as.integer(gv), levels = levels(gv))
}

#' A group-level mixture draw: one class per group from the group's own
#' mixing weights, then each row from its group's component. The rowwise
#' `sim` cannot express this - it would resample the class per row and
#' wash the grouping out - which is why `groups =` installs this
#' instead.
#'
#' @noRd
mixture_sim_groups <- function(ctx) {
  mg <- ctx[["block"]]
  fam <- ctx[["family"]]
  mx <- fam[["mix"]]
  if (is.null(mg)) {
    stop("mixture(groups =) simulated without the group structure: the ",
         "frame carries no latent-class grouping for response '",
         ctx[["resp"]], "'. Rebuild the model frame from the same ",
         "formula and data", call. = FALSE)
  }
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  av <- ctx[["aterms"]]
  K <- mx[["K"]]
  lps <- mx[["log_pi"]](dp)
  Pg <- vapply(lps, function(l) {
    exp(rep(l, length.out = n)[mg[["first"]]])
  }, numeric(length(mg[["first"]])))
  kg <- vapply(seq_len(nrow(Pg)), function(g_) {
    sample.int(K, 1L, prob = Pg[g_, ])
  }, integer(1))
  kk <- kg[mg[["gindex"]]]
  ys <- numeric(n)
  for (k in seq_len(K)) {
    idx <- which(kk == k)
    if (!length(idx)) next
    dk <- lapply(mx[["comp_dpars"]](dp, k), function(v) {
      rep(v, length.out = n)[idx]
    })
    ys[idx] <- mx[["comp_sim"]](dk, av, length(idx), k)
  }
  ys
}

#' Posterior class probabilities of a mixture fit
#'
#' For an ordinary [mixture()] or [mixture_mvn()] fit, one row per
#' observation; for a group-level mixture (`groups = ~g`), one row per
#' group.
#'
#' @param fit A `frmtmb_fit` with a mixture family.
#' @return A matrix of class probabilities (rows sum to one).
#' @examples
#' set.seed(3)
#' dd <- data.frame(y = c(rnorm(80, 0, 1), rnorm(80, 5, 1)))
#' fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
#'
#' p <- mixture_probs(fit)
#' head(p)
#' rowSums(p)[1:3]                    # rows sum to one
#'
#' # the hard assignment, and how well it recovers the truth
#' cl <- max.col(p)
#' table(cl, truth = rep(1:2, each = 80))
#' @export
mixture_probs <- function(fit) {
  rspec <- single_response(fit, "mixture_probs()")
  if (is.null(rspec$family$mix)) {
    stop("mixture_probs() needs a mixture() family fit", call. = FALSE)
  }
  latent_probs(fit)
}

#' The posterior class probabilities of any family implementing the
#' `fam$mix` component interface, which is `mixture()`, `mixture_mvn()`
#' and `lca()`. This is the structures' `latent_probs` slot for all
#' three, and `mixture_probs()` and `lca_probs()` reach it through the
#' `latent_probs()` generic.
#'
#' @noRd
mixture_posterior <- function(fit) {
  rspec <- single_response(fit, "mixture_probs()")
  fam <- rspec$family
  dp <- eval_dpars(fit)[[rspec$resp_name]]
  av <- fit$frame$aterm_values[[rspec$resp_name]]
  yv <- fit$frame$y[[rspec$resp_name]]
  lps_pi <- fam$mix$log_pi(dp)
  K <- fam$mix$K
  # the group incidence of a group-level mixture, or NULL for a rowwise
  # one, where the posterior is per observation
  mg <- frame_block_of(fit$frame, rspec$resp_name)
  w <- av$weights %||% 1
  # matrix responses (mixture_mvn) have one density per ROW, so NROW,
  # not length; their class covariances live in the extra parameters
  ex <- if (!is.null(fam$extra_pars)) {
    fit$estimates[fit$frame$extra_names]
  }
  M <- vapply(seq_len(K), function(k) {
    ll_k <- if (is.null(ex)) {
      fam$mix$comp_lpdf(yv, dp, av, k)
    } else {
      fam$mix$comp_lpdf(yv, dp, av, k, ex)
    }
    if (!is.null(mg)) {
      as.vector(Matrix::t(mg$G) %*% (w * ll_k)) + lps_pi[[k]][mg$first]
    } else {
      ll_k + rep(lps_pi[[k]], length.out = length(ll_k))
    }
  }, numeric(if (!is.null(mg)) length(mg$first) else NROW(yv)))
  P <- exp(M - apply(M, 1, max))
  P <- P / rowSums(P)
  rownames(P) <- if (!is.null(mg)) mg$levels
  colnames(P) <- paste0("class", seq_len(K))
  P
}

# mclust's covariance taxonomy for mixture_mvn(). mclust writes
# Sigma_k = lambda_k * D_k * A_k * D_k' (volume * orientation * shape);
# each letter of the model name says whether volume, shape and
# orientation are Equal across classes or Vary. The subset here is the
# one whose orientation is either the identity (spherical and diagonal
# models) or completely shared/free, which is exactly the subset a
# log-SD / scaled-Cholesky parameterization expresses without
# constrained eigenvector machinery.
mvn_cov_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
                    "EEE", "VVV")

#' The extras a covariance model needs and how they assemble class k's
#' covariance. `pars` is the (name, length) template used for the error
#' messages and the documentation; `init` turns the response into the
#' starting values; `sigma` runs on the tape, so it must stay AD-safe.
#'
#' @noRd
mvn_cov_spec <- function(model, K, D) {
  if (!is.character(model) || length(model) != 1L ||
        !model %in% mvn_cov_models) {
    stop("mixture_mvn(): unknown covariance model '",
         paste(model, collapse = ", "), "'. Supported models: ",
         paste(mvn_cov_models, collapse = ", "), call. = FALSE)
  }
  us_len <- as.integer(D + D * (D - 1L) / 2L)
  # per-column response SDs; the spherical and volume-shape models
  # collapse them to their log-scale mean. unname(): response column
  # names would otherwise leak into the parameter template and give
  # confint()/frm_sample() ragged parameter labels.
  base_ls <- function(y) unname(log(pmax(apply(y, 2, stats::sd), 1e-3)))
  # diagonal covariance from log-SDs (length D, or length 1 recycled)
  diag_S <- function(ls) {
    "[<-" <- RTMB::ADoverload("[<-")
    S <- RTMB::matrix(0, D, D)
    one <- length(ls) == 1L
    for (j in seq_len(D)) S[j, j] <- exp(2 * ls[if (one) 1L else j])
    S
  }
  # volume-shape diagonal: the shape's last log entry is minus the sum
  # of the free ones, which is what pins det(A) = 1 and keeps the
  # volume identified separately from the shape
  volshape_S <- function(vol, sh) {
    "[<-" <- RTMB::ADoverload("[<-")
    S <- RTMB::matrix(0, D, D)
    a_last <- -sum(sh)
    for (j in seq_len(D)) {
      S[j, j] <- exp(2 * (vol[1] + (if (j < D) sh[j] else a_last)))
    }
    S
  }
  cls <- function(prefix) paste0(prefix, seq_len(K))
  spread <- function(nms, v) {
    stats::setNames(rep(list(v), length(nms)), nms)
  }
  free_shape <- function(y) {
    ls <- base_ls(y)
    (ls - mean(ls))[seq_len(D - 1L)]
  }
  switch(
    model,
    EII = list(
      init = function(y) list(sigmaraw = mean(base_ls(y))),
      sigma = function(extra, k) diag_S(extra[["sigmaraw"]])
    ),
    VII = list(
      init = function(y) spread(cls("sigmaraw"), mean(base_ls(y))),
      sigma = function(extra, k) diag_S(extra[[paste0("sigmaraw", k)]])
    ),
    EEI = list(
      init = function(y) list(sigmaraw = base_ls(y)),
      sigma = function(extra, k) diag_S(extra[["sigmaraw"]])
    ),
    VEI = list(
      init = function(y) {
        c(spread(cls("sigmavol"), mean(base_ls(y))),
          list(sigmashape = free_shape(y)))
      },
      sigma = function(extra, k) {
        volshape_S(extra[[paste0("sigmavol", k)]], extra[["sigmashape"]])
      }
    ),
    EVI = list(
      init = function(y) {
        c(list(sigmavol = mean(base_ls(y))),
          spread(cls("sigmashape"), free_shape(y)))
      },
      sigma = function(extra, k) {
        volshape_S(extra[["sigmavol"]], extra[[paste0("sigmashape", k)]])
      }
    ),
    VVI = list(
      init = function(y) spread(cls("sigmaraw"), base_ls(y)),
      sigma = function(extra, k) diag_S(extra[[paste0("sigmaraw", k)]])
    ),
    EEE = list(
      init = function(y) {
        list(sigmaraw = c(base_ls(y), numeric(us_len - D)))
      },
      sigma = function(extra, k) us_sigma(extra[["sigmaraw"]], D)
    ),
    VVV = list(
      init = function(y) {
        spread(cls("sigmaraw"), c(base_ls(y), numeric(us_len - D)))
      },
      sigma = function(extra, k) us_sigma(extra[[paste0("sigmaraw", k)]], D)
    )
  )
}

#' Multivariate gaussian mixture family
#'
#' `mixture_mvn(K, D)` does model-based clustering of an n x D matrix
#' response (mclust-style): K classes, each with its own D-dimensional
#' mean and a D x D covariance from mclust's model taxonomy. Every
#' class mean is a full
#' linear predictor - the main model formula applies to all of them -
#' so cluster means may depend on covariates, which mclust cannot do.
#' The location dpars are named `mu<k>d<j>` (class k, response column
#' j) and are individually overridable, e.g. `bf(Y ~ x, mu2d1 ~ 1)`
#' (all except the first, `mu1d1`). Mixing weights are `theta1 ...
#' theta{K-1}`, multinomial logit against class K, each with its own
#' linear predictor - so gating on covariates works like [mixture()].
#'
#' Class covariances are family-level extra parameters, covariate-free,
#' and their structure follows `model`, mclust's volume-shape-orientation
#' taxonomy for `Sigma_k = lambda_k * D_k * A_k * D_k'`:
#'
#' \tabular{lll}{
#'   `EII` \tab spherical, equal volume
#'     \tab `sigmaraw`, one log-SD \cr
#'   `VII` \tab spherical, varying volume
#'     \tab `sigmaraw<k>`, one log-SD each \cr
#'   `EEI` \tab diagonal, equal volume and shape
#'     \tab `sigmaraw`, D log-SDs \cr
#'   `VEI` \tab diagonal, varying volume, equal shape
#'     \tab `sigmavol<k>` plus `sigmashape` (D - 1) \cr
#'   `EVI` \tab diagonal, equal volume, varying shape
#'     \tab `sigmavol` plus `sigmashape<k>` (D - 1 each) \cr
#'   `VVI` \tab diagonal, free
#'     \tab `sigmaraw<k>`, D log-SDs each \cr
#'   `EEE` \tab one shared full covariance
#'     \tab `sigmaraw`, a `us` block \cr
#'   `VVV` \tab free full covariance per class (default)
#'     \tab `sigmaraw<k>`, one `us` block each
#' }
#'
#' A `us` block is D log-SDs then the scaled-Cholesky correlation
#' entries, as in [frm()]'s `us()` covariance structure. The
#' `sigmashape` vectors hold the first D - 1 log-shape entries; the last
#' is minus their sum, which is what fixes `det(A_k) = 1`. Log-SDs start
#' at the per-column response SDs and correlations at zero; class means
#' start on spread-out per-column response quantiles to break the label
#' symmetry. The usual finite-mixture ML caveats apply: the
#' likelihood is invariant to
#' relabeling and can be multimodal (compare starts via
#' [frm_allfit()]). [mixture_probs()] returns posterior class
#' probabilities per row; [fitted()] returns the n x D mixture-mean
#' matrix. Covariances take no linear predictor (no covariance
#' regression), and the models with a class-varying eigenvector basis
#' (`EEV`, `VEV`, `EVE`, `VEE`, `VVE`, `EVV`) are not available.
#' `cens()`/`trunc()`, [mvbf()], and `simulate()` are not supported.
#'
#' @param K Number of mixture classes (at least 2).
#' @param D Number of response columns (at least 2; for `D = 1` use
#'   `mixture(gaussian(), ...)`).
#' @param model Covariance model name from mclust's vocabulary: one of
#'   `"EII"`, `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`,
#'   `"VVV"` (the default, a free covariance per class).
#' @return A `frmtmb_family`.
#' @examples
#' set.seed(1)
#' Y <- rbind(matrix(rnorm(60, 0), ncol = 2),
#'            matrix(rnorm(60, 4), ncol = 2))
#' dd <- data.frame(row = seq_len(nrow(Y)))
#' dd$Y <- Y
#' fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
#' fixef(fit)
#' head(mixture_probs(fit))
#' # a shared spherical covariance (mclust's EII, k-means-like)
#' frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = "EII"), data = dd)
#' @export
mixture_mvn <- function(K, D, model = "VVV") {
  if (missing(K) || missing(D) || K < 2 || D < 2) {
    stop("mixture_mvn() needs K >= 2 classes and D >= 2 response ",
         "columns (for D = 1 use mixture(gaussian(), ...))",
         call. = FALSE)
  }
  K <- as.integer(K)
  D <- as.integer(D)
  cspec <- mvn_cov_spec(model, K, D)
  mu_names <- paste0("mu", rep(seq_len(K), each = D),
                     "d", rep(seq_len(D), K))
  dpars <- c(mu_names, paste0("theta", seq_len(K - 1L)))
  links <- stats::setNames(rep(list("identity"), length(dpars)), dpars)

  # n x D class-mean matrix from the class's D location dpars
  class_mean <- function(dpars_all, k, n) {
    "[<-" <- RTMB::ADoverload("[<-")
    M <- RTMB::matrix(0, n, D)
    for (j in seq_len(D)) {
      M[, j] <- dpars_all[[paste0("mu", k, "d", j)]]
    }
    M
  }
  # log mixing weights: multinomial logit, last class reference
  log_pi <- function(dpars_all) {
    Ts <- lapply(seq_len(K - 1L), function(k) {
      dpars_all[[paste0("theta", k)]]
    })
    Ts[[K]] <- 0 * dpars_all[[mu_names[1]]]
    lse <- Ts[[1L]]
    for (k in seq.int(2L, K)) lse <- RTMB::logspace_add(lse, Ts[[k]])
    lapply(Ts, function(t_) t_ - lse)
  }
  # per-row class log-density; extra carries the raw covariance
  # parameters, whose layout the covariance model decides
  comp_lpdf <- function(y, dpars_all, aterms, k, extra) {
    Sk <- cspec$sigma(extra, k)
    R <- y - class_mean(dpars_all, k, nrow(y))
    RTMB::dmvnorm(R, 0, Sk, log = TRUE)
  }

  init <- list()
  for (k in seq_len(K)) {
    for (j in seq_len(D)) {
      init[[paste0("mu", k, "d", j)]] <- local({
        k_ <- k
        j_ <- j
        function(y, aterms) {
          stats::quantile(y[, j_], k_ / (K + 1), names = FALSE)
        }
      })
    }
  }

  fam <- frmtmb_family(
    paste0("mixture_mvn(K = ", K, ", D = ", D, ", model = \"",
           model, "\")"),
    dpars = dpars,
    links = links,
    lpdf = function(y, dpars, aterms, extra) {
      lp <- log_pi(dpars)
      ll <- NULL
      for (k in seq_len(K)) {
        llk <- comp_lpdf(y, dpars, aterms, k, extra) + lp[[k]]
        ll <- if (is.null(ll)) llk else RTMB::logspace_add(ll, llk)
      }
      ll
    },
    valid_y = function(y, aterms) {
      if (!is.matrix(y) || ncol(y) != D) {
        stop("mixture_mvn(K = ", K, ", D = ", D, "): response must be ",
             "an n x ", D, " numeric matrix", call. = FALSE)
      }
    },
    init_dpars = init,
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) {
        lp <- log_pi(dpars)
        n <- length(dpars[[mu_names[1]]])
        out <- 0
        for (k in seq_len(K)) {
          Mk <- vapply(seq_len(D), function(j) {
            dpars[[paste0("mu", k, "d", j)]]
          }, numeric(n))
          out <- out + exp(lp[[k]]) * Mk
        }
        out
      }
    ),
    extra_pars = function(y, aterms) {
      # identical covariance starts across classes; the quantile-spread
      # mean inits break the label symmetry
      cspec$init(y)
    },
    primary_dpars = mu_names
  )
  # internals for mixture_probs(): same shape as mixture()'s, with the
  # extra (covariance) parameters as a fifth comp_lpdf argument
  fam$mix <- list(
    K = K,
    D = D,
    cov_model = model,
    comp_lpdf = comp_lpdf,
    comp_dpars = function(dpars_all, k) {
      stats::setNames(
        lapply(seq_len(D), function(j) {
          dpars_all[[paste0("mu", k, "d", j)]]
        }),
        paste0("mud", seq_len(D))
      )
    },
    # class k's covariance from the fit's extra parameters; the same
    # assembler the objective uses, so it also runs off the tape
    sigma = cspec$sigma,
    log_pi = log_pi
  )
  # the class covariances are family-level EXTRAS, not dpars, so the
  # rowwise contract cannot see them; the structured one can
  fam$sim_ctx <- mvn_sim_rows
  # the likelihood is rowwise (one D-variate density per row), so the
  # structure carries nothing but the two multimodality refusals
  fam[["structure"]] <- frmtmb_structure(
    latent_probs = function(fit, block) mixture_posterior(fit),
    supports = structure_supports_all(reml = FALSE, profile = FALSE),
    refusals = mixture_multimodal_refusals("a mixture_mvn() family")
  )
  fam
}

#' One `mixture_mvn()` draw: a class per row from its mixing weights,
#' then a D-variate normal about that class's mean with that class's
#' covariance, assembled from the fit's extra parameters by the same
#' `sigma()` the likelihood uses.
#'
#' @noRd
mvn_sim_rows <- function(ctx) {
  mx <- ctx[["family"]][["mix"]]
  K <- mx[["K"]]
  D <- mx[["D"]]
  n <- ctx[["n"]]
  dp <- ctx[["dpars"]]
  ex <- ctx[["extra"]]
  P <- vapply(mx[["log_pi"]](dp), function(l) {
    rep(exp(l), length.out = n)
  }, numeric(n))
  ks <- vapply(seq_len(n), function(i) {
    sample.int(K, 1L, prob = P[i, ])
  }, integer(1))
  out <- matrix(NA_real_, n, D)
  for (k in seq_len(K)) {
    idx <- which(ks == k)
    if (!length(idx)) next
    # upper Cholesky: Z %*% R has covariance R'R = Sigma_k
    R <- chol(as.matrix(mx[["sigma"]](ex, k)))
    M <- matrix(NA_real_, length(idx), D)
    for (j in seq_len(D)) {
      M[, j] <- rep(dp[[paste0("mu", k, "d", j)]], length.out = n)[idx]
    }
    Z <- matrix(stats::rnorm(length(idx) * D), length(idx), D)
    out[idx, ] <- M + Z %*% R
  }
  out
}

#' von Mises family for a circular response in `(-pi, pi]`, dpars `mu`
#' (the mean direction, through the tan-half link) and `kappa` (the
#' concentration, log link). brms's parameterization exactly.
#'
#' The density needs `log I0(kappa)`, the modified Bessel function of the
#' first kind of order zero. `base::besselI()` is not AD-safe, but RTMB
#' carries its own `besselI` method for advectors, and `RTMBdist::dvm()`
#' already uses it in the exponentially scaled form
#' (`log I0(k) = log besselI(k, 0, expon.scaled = TRUE) + k`), which is
#' what keeps a large concentration from overflowing. So the whole
#' density, and its exact derivative, come off the tape with no series
#' approximation of our own.
#'
#' @noRd
fam_von_mises <- function(link = "tan_half") {
  frmtmb_family(
    "von_mises",
    dpars = c("mu", "kappa"),
    links = list(mu = link, kappa = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMBdist::dvm(y, dpars$mu, dpars$kappa, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y < -pi) || any(y > pi)) {
        stop("von_mises: response must be angles in radians on ",
             "(-pi, pi]; wrap the response first, for example ",
             "atan2(sin(y), cos(y))", call. = FALSE)
      }
    },
    init_dpars = list(
      # the resultant vector's direction and length: the circular
      # analogues of the mean and (through Fisher's approximation of
      # A(kappa) = I1/I0) the precision
      mu = function(y, aterms) atan2(mean(sin(y)), mean(cos(y))),
      kappa = function(y, aterms) {
        rbar <- sqrt(mean(sin(y))^2 + mean(cos(y))^2)
        rbar <- min(max(rbar, 0.02), 0.99)
        if (rbar < 0.53) {
          2 * rbar + rbar^3 + 5 * rbar^5 / 6
        } else if (rbar < 0.85) {
          -0.4 + 1.39 * rbar + 0.43 / (1 - rbar)
        } else {
          1 / (rbar^3 - 4 * rbar^2 + 3 * rbar)
        }
      }
    ),
    type = "continuous",
    post = list(
      # the mean DIRECTION, which is what brms's posterior_epred()
      # reports for this family; a circular response has no mean in the
      # arithmetic sense
      mean_fn = function(dpars, aterms) dpars$mu,
      # the circular variance 1 - A(kappa), on (0, 1)
      var_fn = function(dpars, aterms) {
        k <- dpars$kappa
        1 - besselI(k, 1, expon.scaled = TRUE) /
          besselI(k, 0, expon.scaled = TRUE)
      }
    ),
    sim = function(dpars, aterms, n) {
      rvon_mises(n, rep(dpars$mu, length.out = n),
                 rep(dpars$kappa, length.out = n))
    }
  )
}

#' von Mises draws with a per-row mean direction and concentration, by
#' Best and Fisher's (1979) wrapped-Cauchy rejection scheme. `RTMBdist`
#' ships `rvm()`, but it takes scalar parameters only, and a
#' distributional fit (or `dharma_residuals()`, which draws hundreds of
#' replicates) needs one draw per row of a varying `kappa`. Rejection
#' runs vectorized over the rows still outstanding.
#'
#' @noRd
rvon_mises <- function(n, mu, kappa) {
  mu <- rep(mu, length.out = n)
  kappa <- rep(kappa, length.out = n)
  out <- numeric(n)
  # kappa near zero is the uniform distribution on the circle, and the
  # rejection constants below divide by it
  flat <- kappa < 1e-8
  if (any(flat)) out[flat] <- stats::runif(sum(flat), -pi, pi)
  todo <- which(!flat)
  if (!length(todo)) return(out)
  k <- kappa[todo]
  a <- 1 + sqrt(1 + 4 * k^2)
  b <- (a - sqrt(2 * a)) / (2 * k)
  r <- (1 + b^2) / (2 * b)
  left <- seq_along(todo)
  it <- 0L
  while (length(left)) {
    it <- it + 1L
    if (it > 1000L) {
      stop("von Mises simulation did not converge for ", length(left),
           " of ", n, " rows", call. = FALSE)
    }
    m <- length(left)
    z <- cos(pi * stats::runif(m))
    f <- (1 + r[left] * z) / (r[left] + z)
    c_ <- k[left] * (r[left] - f)
    u2 <- stats::runif(m)
    ok <- c_ * (2 - c_) - u2 > 0 | log(c_ / u2) + 1 - c_ >= 0
    if (any(ok)) {
      idx <- left[ok]
      sgn <- sign(stats::runif(sum(ok)) - 0.5)
      out[todo[idx]] <- mu[todo[idx]] + sgn * acos(pmin(pmax(f[ok], -1), 1))
    }
    left <- left[!ok]
  }
  # back onto (-pi, pi], the support the family declares
  (out + pi) %% (2 * pi) - pi
}

# --- categorical (nominal) responses ---------------------------------

#' Multinomial-logit log-density over a CATEGORY INDEX response `1..K`.
#' `eta` is the list of K-1 non-reference linear predictors in category
#' order, category 1 being the reference with `eta_1 = 0`; the shared
#' code is factored out so the taped path (`ord_cat_sel`, which needs no
#' comparison operators) and the data path agree term for term.
#'
#' @noRd
cat_logit_lpdf <- function(y, eta, K, sel = NULL) {
  if (is.null(sel)) sel <- lapply(seq_len(K), function(k) {
    as.numeric(y == k)
  })
  # the denominator accumulates in LOG space from the reference
  # category's implicit zero, so a wide predictor cannot overflow it;
  # summing exp(eta_j) first loses the model at eta = 709
  lden <- 0 * eta[[1L]]
  num <- 0
  for (j in seq_len(K - 1L)) {
    lden <- RTMB::logspace_add(lden, eta[[j]])
    num <- num + sel[[j + 1L]] * eta[[j]]
  }
  num - lden
}

#' Categorical (nominal) family over a `1..K` category index, category 1
#' the reference. One linear predictor per non-reference category, all
#' receiving the main model formula unless overridden.
#'
#' `dpar_names` is the brms spelling `mu<Level>` when the levels are
#' known and the positional `mu2 ... muK` when only `K` is.
#'
#' @noRd
fam_categorical_impl <- function(dpar_names, levels = NULL,
                                 link = "logit") {
  if (!identical(link, "logit")) {
    stop("categorical: the multinomial logit is the only link, so a ",
         "family object carrying link '", link, "' cannot be used here",
         call. = FALSE)
  }
  K <- length(dpar_names) + 1L
  fam <- frmtmb_family(
    "categorical",
    dpars = dpar_names,
    links = stats::setNames(rep(list("identity"), K - 1L), dpar_names),
    lpdf = function(y, dpars, aterms) {
      eta <- lapply(dpar_names, function(nm) dpars[[nm]])
      ov <- osa_unwrap(y)
      if (!is.null(ov)) {
        # the response is on the tape: pick the category arithmetically
        return(cat_logit_lpdf(NULL, eta, K,
                              sel = ord_cat_sel(ov$y, K)) * ov$keep)
      }
      cat_logit_lpdf(y, eta, K)
    },
    valid_y = function(y, aterms) {
      if (any(y < 1) || any(y > K) || any(y != round(y))) {
        stop("categorical: response must be a factor with ", K,
             " levels, or integer category codes 1..", K, call. = FALSE)
      }
      if (length(unique(y)) < 2) {
        stop("categorical: the response takes only one value; there is ",
             "nothing to model", call. = FALSE)
      }
    },
    type = "categorical",
    sim = function(dpars, aterms, n) {
      eta <- lapply(dpar_names, function(nm) {
        rep(dpars[[nm]], length.out = n)
      })
      E <- cbind(0, do.call(cbind, eta))
      P <- exp(E - apply(E, 1L, max))
      P <- P / rowSums(P)
      cp <- t(apply(P, 1L, cumsum))
      if (n == 1L) cp <- matrix(cp, 1L, K)
      pmin(1L + rowSums(cp < stats::runif(n)), K)
    },
    primary_dpars = dpar_names
  )
  fam$cat_levels <- levels
  fam$cat_K <- K
  fam
}

#' Placeholder returned by a bare `categorical()`: the category count is
#' a property of the data, and the formula grammar is resolved before any
#' data arrives, so `frm()` swaps this for the real family through
#' `resolve_deferred_families()` once it can see the response. It carries
#' one `mu` so a spec built from it is still well formed; anything that
#' reaches the likelihood with it in place says so.
#'
#' @noRd
fam_categorical_deferred <- function(link = "logit") {
  fam <- frmtmb_family(
    "categorical",
    dpars = "mu",
    links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms) {
      stop("categorical(): the response categories were never resolved. ",
           "They are read from the data by frm(); on another entry point ",
           "name them, categorical(levels = c(\"a\", \"b\", \"c\")) or ",
           "categorical(K = 3)", call. = FALSE)
    },
    type = "categorical"
  )
  fam$defer <- function(formula, data) {
    lv <- categorical_levels(formula, data)
    fam_categorical_impl(paste0("mu", lv[-1L]), levels = lv, link = link)
  }
  fam
}

#' The response's category levels, for the deferred `categorical()`. The
#' response expression is evaluated against the data with the addition
#' terms stripped; a character response becomes a factor here, and says
#' which order it took.
#'
#' @noRd
categorical_levels <- function(formula, data) {
  ri <- parse_response(formula)
  y <- tryCatch(eval(ri$resp, data, environment(formula) %||% globalenv()),
                error = function(e) NULL)
  if (is.null(y)) {
    stop("categorical(): the response '", deparse1(ri$resp),
         "' could not be evaluated on the data to find its categories. ",
         "Name them instead: categorical(levels = c(\"a\", \"b\"))",
         call. = FALSE)
  }
  lv <- categorical_y_levels(y, deparse1(ri$resp))
  if (length(lv) < 2L) {
    stop("categorical(): the response '", deparse1(ri$resp),
         "' has fewer than two categories", call. = FALSE)
  }
  lv
}

#' The category labels of a categorical response, in the order that fixes
#' the reference category and the dpar names. A character vector is
#' coerced with a message naming the order it took, because that order is
#' the model; a factor keeps its own levels.
#'
#' @noRd
categorical_y_levels <- function(y, label) {
  if (is.character(y) || is.logical(y)) {
    lv <- levels(factor(y))
    message("Categorical response '", label, "' is a ",
            if (is.logical(y)) "logical" else "character",
            " vector; it is read as a factor with levels ",
            paste(lv, collapse = ", "), " and '", lv[1L],
            "' as the reference category. Set the order with factor() ",
            "if that is not what you want.")
    return(lv)
  }
  if (is.factor(y)) return(levels(y))
  NULL
}

#' Swap every deferred family of a bform for the concrete one its data
#' implies. Called by `frm()` before parsing, so a per-category dpar
#' formula (`bf(y ~ x, mub ~ z)`) reaches the parser with the dpar
#' already in the family's vocabulary.
#'
#' @noRd
resolve_deferred_families <- function(bform, data) {
  if (is.null(data)) return(bform)
  one <- function(f) {
    if (is.null(f$family) || is.null(f$family$defer)) return(f)
    f$family <- f$family$defer(f$formula, data)
    f
  }
  if (inherits(bform, "frmtmb_mvformula")) {
    bform$forms <- lapply(bform$forms, one)
  } else {
    bform <- one(bform)
  }
  bform
}

# --- Cox proportional hazards ----------------------------------------

#' B-spline basis of a given order (degree + 1) over a full knot vector,
#' by the Cox-de Boor recursion. Written out rather than taken from
#' `splines::splineDesign()` so the baseline hazard needs no dependency
#' outside the ones the package already carries; it agrees with
#' `splines2` to machine precision, which the tests assert.
#'
#' @noRd
bspline_basis <- function(x, knots, order) {
  n <- length(x)
  nb1 <- length(knots) - 1L
  # order 1: the indicator of [t_j, t_{j+1}); the last non-degenerate
  # interval closes on the right so the upper boundary is covered
  B <- matrix(0, n, nb1)
  last <- max(which(knots < knots[length(knots)]))
  for (j in seq_len(nb1)) {
    B[, j] <- if (j == last) {
      as.numeric(x >= knots[j] & x <= knots[j + 1L])
    } else {
      as.numeric(x >= knots[j] & x < knots[j + 1L])
    }
  }
  if (order == 1L) return(B[, seq_len(length(knots) - order), drop = FALSE])
  for (k in 2:order) {
    Bn <- matrix(0, n, length(knots) - k)
    for (j in seq_len(ncol(Bn))) {
      d1 <- knots[j + k - 1L] - knots[j]
      d2 <- knots[j + k] - knots[j + 1L]
      # a repeated knot gives a zero-width span; its term drops out
      t1 <- if (d1 > 0) (x - knots[j]) / d1 * B[, j] else 0
      t2 <- if (d2 > 0) (knots[j + k] - x) / d2 * B[, j + 1L] else 0
      Bn[, j] <- t1 + t2
    }
    B <- Bn
  }
  B
}

#' The knot placement brms uses for the Cox baseline: internal knots on
#' equally spaced response quantiles, boundary knots just outside the
#' observed range so every observation sits strictly inside the basis
#' (`splines2`'s default through brms's `bhaz_basis_matrix()`).
#'
#' @noRd
bhaz_spec <- function(y, df, degree, intercept) {
  n_int <- df - degree - as.integer(intercept)
  if (n_int < 0L) {
    stop("cox(): df must be at least degree + 1 (", degree + 1L,
         ") with an intercept in the baseline basis", call. = FALSE)
  }
  rng <- range(y)
  d <- rng[2L] - rng[1L]
  boundary <- c(max(rng[1L] - d / 50, 0), rng[2L] + d / 50)
  internal <- if (n_int > 0L) {
    p <- seq(0, 1, length.out = n_int + 2L)[-c(1L, n_int + 2L)]
    unname(stats::quantile(y, probs = p))
  } else {
    numeric(0)
  }
  list(internal = internal, boundary = boundary, degree = degree,
       intercept = intercept, df = df)
}

#' The full knot vector of a baseline-hazard spec: boundary knots
#' repeated to the spline's order at each end.
#'
#' @noRd
bhaz_knots <- function(sp, extra = 0L) {
  k <- sp$degree + 1L + extra
  c(rep(sp$boundary[1L], k), sp$internal, rep(sp$boundary[2L], k))
}

#' The M-spline basis of the baseline hazard: the B-spline basis scaled
#' so each function integrates to one over its support, which is what
#' makes a simplex of coefficients a density-like hazard.
#'
#' @noRd
mspline_design <- function(x, sp) {
  k <- sp$degree + 1L
  tt <- bhaz_knots(sp)
  N <- bspline_basis(x, tt, k)
  M <- N
  for (j in seq_len(ncol(N))) M[, j] <- N[, j] * k / (tt[j + k] - tt[j])
  if (!sp$intercept) M <- M[, -1L, drop = FALSE]
  M
}

#' The I-spline basis: `I_j(x)` is the integral of `M_j` from the lower
#' boundary knot, so it is the CUMULATIVE baseline hazard's basis and it
#' is monotone by construction. It equals the reverse cumulative sum of
#' the order-(degree + 2) B-spline basis on the once-more-repeated
#' boundary knots (Ramsay 1988); the tests check it against both
#' `splines2::iSpline()` and a numeric integral of `mspline_design()`.
#'
#' @noRd
ispline_design <- function(x, sp) {
  N2 <- bspline_basis(x, bhaz_knots(sp, extra = 1L), sp$degree + 2L)
  R <- t(apply(N2, 1L, function(v) rev(cumsum(rev(v)))))
  if (length(x) == 1L) R <- matrix(R, 1L, ncol(N2))
  I <- R[, -1L, drop = FALSE]
  if (!sp$intercept) I <- I[, -1L, drop = FALSE]
  I
}

#' Cox proportional hazards with a flexible (M-spline) baseline hazard,
#' brms's construction. The baseline hazard is `Zbhaz %*% s` and the
#' cumulative baseline hazard `Zcbhaz %*% s` over the SAME simplex `s`,
#' with `Zcbhaz` the I-spline integral of the M-spline `Zbhaz`. The
#' simplex is what identifies the baseline against the intercept, and it
#' rides in the parameter vector as `K - 1` unconstrained values through
#' a softmax with the first component pinned at zero.
#'
#' @noRd
fam_cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE) {
  df <- as.integer(df)
  degree <- as.integer(degree)
  intercept <- isTRUE(intercept)
  # the simplex from its K-1 free coordinates; softmax with the first
  # component pinned at zero covers the whole open simplex
  sbhaz <- function(raw) {
    "c" <- RTMB::ADoverload("c")   # base c() strips the advector class
    e <- exp(c(0, raw))
    e / sum(e)
  }
  fam <- frmtmb_family(
    "cox",
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms, extra) {
      s <- sbhaz(extra$sbhaz_raw)
      bhaz <- as.vector(aterms$Zbhaz %*% s)
      cbhaz <- as.vector(aterms$Zcbhaz %*% s)
      # log h(t) + log S(t), with h(t) = h0(t) * mu and mu = exp(eta)
      log(bhaz) + log(dpars$mu) - cbhaz * dpars$mu
    },
    lcdf = function(q, dpars, aterms, extra) {
      # F(t) = 1 - exp(-H0(t) * mu). Every bound a censored or truncated
      # row asks about is plain data, so its I-spline row block was
      # built once at frame assembly and is looked up by value here.
      cbhaz <- as.vector(cox_cbhaz_design(q, aterms) %*%
                           sbhaz(extra$sbhaz_raw))
      1 - exp(-cbhaz * dpars$mu)
    },
    valid_y = function(y, aterms) {
      if (any(y <= 0)) {
        stop("cox: the response is a survival time and must be ",
             "strictly positive", call. = FALSE)
      }
      if (length(unique(y)) < df) {
        stop("cox: fewer distinct event times (", length(unique(y)),
             ") than baseline basis functions (", df,
             "); lower df in cox(df = )", call. = FALSE)
      }
    },
    init_dpars = list(
      # the baseline integrates to one over the observed window, so the
      # hazard ratio starts at the crude event rate
      mu = function(y, aterms) {
        ev <- if (is.null(aterms$cens)) 1 else mean(aterms$cens == 0)
        max(ev, 0.05) / mean(y)
      }
    ),
    type = "continuous",
    post = list(
      # mu is the hazard RATIO, not a mean, and reporting it as one is
      # exactly the mistake this family invites. brms refuses the same
      # question (posterior_epred has no cox method).
      mean_fn = function(dpars, aterms) {
        stop("cox: a survival time has no mean on the response scale ",
             "here - the model fits a hazard, and the mean survival ",
             "time would be an integral over the baseline that the ",
             "censored rows do not identify. predict(type = \"link\") ",
             "gives the log hazard ratio and cox_baseline() the fitted ",
             "baseline weights", call. = FALSE)
      }
    ),
    extra_pars = function(y, aterms) {
      list(sbhaz_raw = rep(0, ncol(aterms$Zbhaz) - 1L))
    },
    # a deliberate omission, not a gap: every entry point repeats this
    # reason after its own refusal
    sim_refusal = paste0(
      "Drawing a survival time means inverting the cumulative baseline ",
      "hazard, and cox() carries no quantile function for it: the ",
      "I-spline baseline is identified only on the observed time ",
      "window, so a draw beyond the last event time has no defined ",
      "distribution"
    )
  )
  # the baseline bases are data, not parameters: they are built once
  # from the observed times and ride with the response's addition terms
  fam$aterm_data <- function(y, aterms) {
    sp <- bhaz_spec(y, df, degree, intercept)
    out <- list(Zbhaz = mspline_design(y, sp), bhaz_spec = sp,
                cox_cb = list(cox_cb_entry(y, sp)))
    for (nm in c("cens_y2", "trunc_lb", "trunc_ub")) {
      if (!is.null(aterms[[nm]])) {
        out$cox_cb[[length(out$cox_cb) + 1L]] <-
          cox_cb_entry(aterms[[nm]], sp)
      }
    }
    out$Zcbhaz <- out$cox_cb[[1L]]$Z
    out
  }
  fam$cox_sbhaz <- sbhaz
  fam
}

#' One cached (bound, I-spline design) pair. A bound outside the
#' boundary knots is clamped first: below the lower knot the cumulative
#' baseline hazard is zero and above the upper one it is its total, and
#' the raw B-spline recursion returns an all-zero row for both, which
#' would read as "no risk accumulated" at an upper bound.
#'
#' @noRd
cox_cb_entry <- function(q, sp) {
  qq <- pmin(pmax(as.numeric(q), sp$boundary[1L]), sp$boundary[2L])
  list(q = as.numeric(q), Z = ispline_design(qq, sp))
}

#' The cumulative-baseline-hazard design for a bound the Cox likelihood
#' asks about. Frame assembly cached one per bound; anything else (a
#' post-fit query) is built on the spot.
#'
#' @noRd
cox_cbhaz_design <- function(q, aterms) {
  qn <- as.numeric(q)
  for (e in aterms$cox_cb %||% list()) {
    if (identical(e$q, qn)) return(e$Z)
  }
  cox_cb_entry(qn, aterms$bhaz_spec)$Z
}

#' Call a family's log-CDF, passing the family-level extra parameters to
#' the families that declare them. The Cox baseline lives there, so its
#' survivor function needs them; every other CDF keeps the three-argument
#' contract.
#'
#' @noRd
fam_lcdf <- function(fam, q, dpars, aterms, extra) {
  if (length(formals(fam$lcdf)) >= 4L) {
    fam$lcdf(q, dpars, aterms, extra)
  } else {
    fam$lcdf(q, dpars, aterms)
  }
}

#' The fitted baseline-hazard simplex of a `cox()` fit.
#'
#' @param fit A `frmtmb_fit` with a [cox()] family.
#' @return The `Kbhaz` M-spline weights, summing to one.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200))
#' dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
#' fit <- frm(bf(time ~ x), family = cox(), data = dd)
#' cox_baseline(fit)
#' @export
cox_baseline <- function(fit) {
  rspec <- single_response(fit, "cox_baseline()")
  fam <- rspec$family
  if (is.null(fam$cox_sbhaz)) {
    stop("cox_baseline() needs a cox() family fit", call. = FALSE)
  }
  s <- fam$cox_sbhaz(fit$estimates[["sbhaz_raw"]])
  stats::setNames(s, paste0("s", seq_along(s)))
}

#' Matrix-response multinomial: y is an n x K count matrix, category 1
#' is the reference. One linear predictor per non-reference category
#' (`mu2`, ..., `muK`), all receiving the main model formula unless
#' overridden.
#'
#' @noRd
fam_multinomial <- function(K) {
  if (missing(K) || K < 2) {
    stop("multinomial() needs the number of categories, e.g. ",
         "multinomial(K = 3)", call. = FALSE)
  }
  dpn <- paste0("mu", seq_len(K)[-1])
  frmtmb_family(
    "multinomial",
    dpars = dpn,
    links = stats::setNames(rep(list("identity"), K - 1L), dpn),
    lpdf = function(y, dpars, aterms) {
      # log-space denominator, as in cat_logit_lpdf(): a count response
      # multiplies it by the row total, so an overflow there is fatal
      lden <- 0 * dpars[[dpn[1L]]]
      for (k in dpn) lden <- RTMB::logspace_add(lden, dpars[[k]])
      nvec <- rowSums(y)
      ll <- -nvec * lden
      for (j in seq_along(dpn)) {
        ll <- ll + y[, j + 1L] * dpars[[dpn[j]]]
      }
      # multinomial coefficient: data-only
      ll + lgamma(nvec + 1) - rowSums(lgamma(y + 1))
    },
    valid_y = function(y, aterms) {
      if (!is.matrix(y) || ncol(y) != K) {
        stop("multinomial(K = ", K, "): response must be an n x ", K,
             " count matrix", call. = FALSE)
      }
      if (any(y < 0) || any(y != round(y))) {
        stop("multinomial: response must be non-negative integer counts",
             call. = FALSE)
      }
      tr <- aterms$trials
      if (!is.null(tr) && any(rowSums(y) != tr)) {
        stop("multinomial: rowSums(response) must equal trials()",
             call. = FALSE)
      }
    },
    type = "discrete",
    sim = function(dpars, aterms, n) {
      size <- aterms$trials
      if (is.null(size)) {
        stop("simulate(): a multinomial fit needs trials() to know how ",
             "many draws each row gets", call. = FALSE)
      }
      denom <- 1
      for (k in dpn) denom <- denom + exp(dpars[[k]])
      P <- matrix(0, n, K)
      P[, 1L] <- rep(1 / denom, length.out = n)
      for (j in seq_along(dpn)) {
        P[, j + 1L] <- rep(exp(dpars[[dpn[j]]]) / denom, length.out = n)
      }
      size <- rep(size, length.out = n)
      out <- matrix(0L, n, K)
      for (i in seq_len(n)) {
        out[i, ] <- stats::rmultinom(1L, size[i], P[i, ])
      }
      out
    },
    primary_dpars = dpn
  )
}

family_registry <- list(
  gaussian    = fam_gaussian,
  poisson     = fam_poisson,
  binomial    = fam_binomial,
  Gamma       = fam_Gamma,
  lognormal   = fam_lognormal,
  student     = fam_student,
  negbinomial = fam_negbinomial,
  nbinom2     = fam_negbinomial,
  nbinom1     = fam_nbinom1,
  beta        = fam_beta,
  Beta        = fam_beta,
  tweedie     = fam_tweedie,
  compois     = fam_compois,
  zero_inflated_poisson     = fam_zi_poisson,
  zero_inflated_negbinomial = fam_zi_negbinomial,
  hurdle_poisson            = fam_hurdle_poisson,
  multinomial               = fam_multinomial,
  cumulative                = fam_cumulative,
  beta_binomial             = fam_beta_binomial,
  skew_normal               = fam_skew_normal,
  inverse.gaussian          = fam_inverse_gaussian,
  exgaussian                = fam_exgaussian,
  bernoulli                 = fam_bernoulli,
  geometric                 = fam_geometric,
  exponential               = fam_exponential,
  weibull                   = fam_weibull,
  shifted_lognormal         = fam_shifted_lognormal,
  hurdle_gamma              = fam_hurdle_gamma,
  hurdle_lognormal          = fam_hurdle_lognormal,
  zero_inflated_binomial    = fam_zi_binomial,
  zero_inflated_beta        = fam_zi_beta,
  asym_laplace              = fam_asym_laplace,
  zero_inflated_asym_laplace = fam_zi_asym_laplace,
  huber                     = fam_huber,
  sratio                    = fam_sratio,
  cratio                    = fam_cratio,
  acat                      = fam_acat,
  von_mises                 = fam_von_mises,
  cox                       = fam_cox,
  categorical               = fam_categorical_deferred
)

# An unknown value errors and names the supported families.

#' @rdname frmtmb-extension-api
#' @export
as_frmtmb_family <- function(x) {
  if (inherits(x, "frmtmb_family")) return(x)
  # A bare constructor (`family = cumulative`, `family = gaussian`) is
  # legal in brms and in stats::glm, so call it with its defaults. The
  # frmtmb constructors return a `frmtmb_family` rather than a
  # `stats::family`, so the result is re-checked here.
  if (is.function(x)) {
    x <- x()
    if (inherits(x, "frmtmb_family")) return(x)
  }
  if (inherits(x, "family")) {
    ctor <- family_registry[[x$family]]
    if (is.null(ctor)) {
      stop("Unsupported family object: '", x$family,
           "'. Currently supported: ",
           paste(unique(names(family_registry)), collapse = ", "),
           call. = FALSE)
    }
    return(ctor(link = x$link))
  }
  if (is.character(x) && length(x) == 1) {
    ctor <- family_registry[[x]]
    if (is.null(ctor)) {
      stop("Unsupported family name: '", x, "'. Currently supported: ",
           paste(unique(names(family_registry)), collapse = ", "),
           call. = FALSE)
    }
    return(ctor())
  }
  stop("Cannot interpret `family` of class ",
       paste(class(x), collapse = "/"),
       ": pass a family constructor or its name as the `family` ",
       "argument of frm(), or attach it with `bf(...) + gaussian()`",
       call. = FALSE)
}

#' Additional response families
#'
#' Family constructors without a [stats::family] equivalent, following
#' brms naming. `gaussian()`, `poisson()`, `binomial()`, and `Gamma()`
#' from 'stats' are accepted directly by [frm()] and `+`.
#'
#' An ordinal family (`cumulative()`, `sratio()`, `cratio()`, `acat()`)
#' takes the response's level order as the category order. Supply an
#' ordered factor, or integer codes `1..K`: an unordered factor is
#' accepted, as brms accepts it, but warns and names the order it is
#' about to use, which is alphabetical unless the levels were set.
#'
#' @section Categorical (nominal) responses:
#' `categorical()` fits a multinomial logit to an unordered factor. The
#' FIRST level is the reference category, its linear predictor is held
#' at zero, and each remaining level gets its own predictor named
#' `mu<Level>`, as in brms. The main model formula applies to every one
#' of them, and any single category is overridden by naming it:
#' `bf(y ~ x, mub ~ z)` gives category `"b"` its own predictor and
#' leaves the rest on `~ x`.
#'
#' The categories come from the data, so `frm()` reads them off the
#' response before it parses the formula. Constructing the family away
#' from a data set (to inspect it, or to reach the parser through
#' another entry point) needs them stated: `categorical(levels = c("a",
#' "b", "c"))`, or `categorical(K = 3)` for a response already coded
#' `1..K`, which names the dpars `mu2 ... muK`. A character or logical
#' response is coerced to a factor with a message naming the level
#' order, because that order is the model.
#'
#' `fitted()` and `predict(type = "response")` return the `n x K` matrix
#' of category probabilities, columns named by the response's own
#' levels and rows summing to one - the same convention the ordinal
#' families follow. `predict(type = "link")` and `predict(dpar =)` give
#' the per-category latent predictors, which is where `se.fit` lives.
#' `simulate()` draws factor levels. The same likelihood is available on
#' a count-matrix response as `multinomial(K)`, and a one-hot matrix
#' gives an identical log-likelihood.
#'
#' @section Circular responses:
#' `von_mises()` models an angle in radians on `(-pi, pi]`. Its `mu` is
#' the mean direction and takes the `tan_half` link, which maps the
#' whole line onto that interval; `kappa` is the concentration and takes
#' a log link, with `kappa = 0` the uniform distribution on the circle.
#' Both are brms's choices. `fitted()` and `predict(type = "response")`
#' report the mean direction. The normalizing constant needs
#' `log I0(kappa)`, which RTMB differentiates exactly through its own
#' `besselI` method, so nothing here is a series approximation.
#' Residuals are differences of angles and are NOT wrapped, so read
#' `residuals()` on a von Mises fit with that in mind.
#'
#' @section Cox proportional hazards:
#' `cox()` is the flexible-parametric proportional hazards model brms
#' fits: the baseline hazard is an M-spline in time,
#' `h0(t) = sum_j s_j M_j(t)`, and the cumulative baseline hazard is the
#' I-spline integral of the same basis over the same weights, which
#' makes it monotone by construction. The weights `s` form a simplex -
#' that is what identifies the baseline against the intercept - and are
#' estimated as `sbhaz_raw`, their `Kbhaz - 1` free softmax
#' coordinates; [cox_baseline()] returns the simplex itself. The
#' default basis is brms's: `df = 5` cubic M-splines with an intercept,
#' internal knots on response quantiles and boundary knots just outside
#' the observed range.
#'
#' The hazard is `h0(t) exp(eta)`, so a coefficient is a log hazard
#' ratio, exactly as in `survival::coxph()`. Right, left, and interval
#' censoring come through the ordinary `cens()` addition term: an event
#' contributes the density and a censored observation the survivor
#' function, which is what this family's log-density and log-CDF are.
#' Random effects are the point: `time | cens(c) ~ x + (1 | g)` is a
#' frailty model, and the Laplace approximation integrates the frailties
#' out.
#'
#' The baseline is semiparametric only in spirit - it has `df`
#' parameters, not one per event time - so coefficients agree with
#' `coxph()` closely rather than exactly. A survival response has no
#' mean, so `fitted()` and `predict(type = "response")` are refused;
#' `predict(type = "link")` gives the log hazard ratio. `simulate()` is
#' not available.
#'
#' Maximum likelihood often puts one or more baseline weights ON the
#' simplex boundary, at exactly zero. Their softmax coordinates then run
#' off to minus infinity along a flat ridge, the Hessian is singular in
#' those directions, and the optimizer reports singular convergence even
#' though the gradient is zero and the regression coefficients are at
#' their optimum - [diagnose()] names `sbhaz_raw` as the culprit. This
#' is what an unpenalized flexible baseline does; brms does not meet it
#' because its Dirichlet prior keeps the weights interior. Lower `df`
#' until the baseline is one the data supports, and read
#' [cox_baseline()] to see which weights collapsed.
#'
#' @section Quantile regression inference:
#' `asym_laplace()` and `zero_inflated_asym_laplace()` fit quantile
#' regression through a WORKING likelihood: at a fixed `quantile` the
#' point estimates are consistent quantile estimates (they match
#' `quantreg::rq()`), but the asymmetric Laplace is not the data's
#' true density, so Wald standard errors and `confint()` intervals
#' computed from it are not calibrated. This is a property of the
#' asymmetric-Laplace approach, shared by brms. Use
#' [frm_bootstrap()] for intervals you can defend. The check
#' function's kink can also produce a benign false-convergence
#' warning near the optimum; `frm_allfit()` confirms the fit when in
#' doubt.
#'
#' @section Robust regression:
#' `huber()` fits Huber's least-favorable distribution: gaussian within
#' `k` residual standard deviations of `mu` and Laplace outside, so a
#' far-out point pulls on the fit with a bounded influence instead of
#' its squared distance. It is a proper normalized density, not a
#' penalty, so this is ordinary maximum likelihood and `logLik()`,
#' `AIC()` and the likelihood-ratio machinery all mean what they say.
#'
#' `k` is a fixed constant of the family, `huber(k = 1.345)`, not a
#' distributional parameter. It states where the analyst draws the line
#' between a residual and an outlier, which is a modelling choice rather
#' than something the data identifies; `MASS::rlm()` treats it the same
#' way. Estimating it would let the likelihood buy fit by widening the
#' gaussian core, which is the opposite of the point. As `k` grows the
#' family collapses to `gaussian()`.
#'
#' Point estimates track `MASS::rlm(psi = psi.huber)` closely but not
#' exactly, and the difference is the scale: `rlm()` fixes the scale at
#' a MAD-type estimate and iterates the location, while `huber()`
#' estimates `sigma` by maximum likelihood jointly with `mu`. The
#' coefficients agree to about `1e-2` on well-behaved data; the
#' `sigma` estimates need not.
#'
#' The working-likelihood caveat that applies to `asym_laplace()`
#' applies here for the same reason. If the data are not actually
#' Huber-distributed - and the family is chosen precisely because the
#' error distribution is unknown - then the model is misspecified, the
#' information matrix is not the variance of the score, and Wald
#' standard errors and `confint()` intervals from it are not
#' calibrated. The point estimates stay consistent for the location.
#' Use [frm_bootstrap()] for intervals you can defend.
#' `cens()` and `trunc()` are unavailable: the CDF is a three-piece
#' function of a parameter-dependent residual and has no branch-free
#' form for the tape.
#'
#' `rho` has a kink at `|u| = k`, so the objective is only piecewise
#' smooth and the optimizer often stops with a maximum absolute
#' gradient around `1e-4` and the accompanying false-convergence
#' warning. That is the kink, not a bad fit: the same thing happens to
#' `asym_laplace()`. The estimates satisfy Huber's own estimating
#' equations, `X' psi(u) = 0` with `psi(u) = min(max(u, -k), k)`, to
#' the same order. [frm_allfit()] confirms the fit when in doubt.
#'
#' @param link Link for `mu`.
#' @return A `frmtmb_family` object.
#' @examples
#' set.seed(4)
#' n <- 120
#' dd <- data.frame(x = rnorm(n))
#'
#' # heavier tails than gaussian(), with an estimated df
#' dd$y <- 1 + 0.8 * dd$x + rt(n, df = 4)
#' fixef(frm(bf(y ~ x) + student(), data = dd))
#'
#' # counts with more spread than poisson() allows
#' dd$cnt <- rnbinom(n, mu = exp(0.5 + 0.4 * dd$x), size = 2)
#' fit <- frm(bf(cnt ~ x) + negbinomial(), data = dd)
#' fixef(fit)$mu
#'
#' # a zero-inflated count: the zi dpar gets its own predictor
#' dd$zi <- ifelse(runif(n) < 0.3, 0, dd$cnt)
#' frm(bf(zi ~ x, zi ~ 1) + zero_inflated_poisson(), data = dd)
#'
#' # an ordered response: level order is the category order
#' dd$grade <- cut(1 + 0.8 * dd$x + rlogis(n), 3,
#'                 labels = c("low", "mid", "high"), ordered_result = TRUE)
#' frm(bf(grade ~ x) + cumulative(), data = dd)
#'
#' # a proportion in (0, 1)
#' dd$p <- plogis(0.2 + 0.6 * dd$x + rnorm(n, 0, 0.3))
#' frm(bf(p ~ x) + Beta(), data = dd)
#'
#' # bounded influence: a few wild points barely move the slope
#' dd$rob <- 1 + 0.8 * dd$x + rnorm(n)
#' dd$rob[1:5] <- dd$rob[1:5] + 30
#' fixef(frm(bf(rob ~ x), family = huber(), data = dd))$mu
#' fixef(frm(bf(rob ~ x), family = gaussian(), data = dd))$mu
#'
#' # an unordered factor: one predictor per non-reference category,
#' # named after the level it belongs to
#' dd$pick <- factor(sample(c("ale", "stout", "lager"), n, TRUE))
#' cat_fit <- frm(bf(pick ~ x), family = categorical(), data = dd)
#' fixef(cat_fit)                     # mulager and mustout; ale is the
#'                                    # reference
#' head(fitted(cat_fit))              # n x K category probabilities
#'
#' # one category may take its own predictor
#' dd$w <- rnorm(n)
#' frm(bf(pick ~ x, mustout ~ w), family = categorical(), data = dd)
#'
#' # an angle: mu is the mean direction, kappa the concentration
#' dd$angle <- atan2(sin(0.5 + dd$x), cos(0.5 + dd$x))
#' vm_fit <- frm(bf(angle ~ x), family = von_mises(), data = dd)
#' head(fitted(vm_fit))               # the mean direction, in radians
#'
#' # proportional hazards with a spline baseline; (1 | g) is a frailty
#' dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
#' dd$out <- rbinom(n, 1, 0.3)        # 1 = right censored
#' cox_fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
#' fixef(cox_fit)$mu                  # log hazard ratios
#' cox_baseline(cox_fit)              # the baseline hazard weights
#' @name frmtmb-families
NULL

#' @rdname frmtmb-families
#' @export
student <- function(link = "identity") fam_student(link)

#' @rdname frmtmb-families
#' @export
lognormal <- function(link = "identity") fam_lognormal(link)

#' @rdname frmtmb-families
#' @export
negbinomial <- function(link = "log") fam_negbinomial(link)

#' @rdname frmtmb-families
#' @export
nbinom1 <- function(link = "log") fam_nbinom1(link)

#' @rdname frmtmb-families
#' @export
Beta <- function(link = "logit") fam_beta(link)

#' @rdname frmtmb-families
#' @export
tweedie <- function(link = "log") fam_tweedie(link)

#' @rdname frmtmb-families
#' @export
compois <- function(link = "log") fam_compois(link)

#' @rdname frmtmb-families
#' @export
zero_inflated_poisson <- function(link = "log") fam_zi_poisson(link)

#' @rdname frmtmb-families
#' @export
zero_inflated_negbinomial <- function(link = "log") fam_zi_negbinomial(link)

#' @rdname frmtmb-families
#' @export
hurdle_poisson <- function(link = "log") fam_hurdle_poisson(link)

#' @rdname frmtmb-families
#' @param K For `multinomial()`: number of response categories (columns
#'   of the count-matrix response); category 1 is the reference.
#' @export
multinomial <- function(K) fam_multinomial(K)

#' @rdname frmtmb-families
#' @export
cumulative <- function(link = "logit") fam_cumulative(link)

#' @rdname frmtmb-families
#' @export
beta_binomial <- function(link = "logit") fam_beta_binomial(link)

#' @rdname frmtmb-families
#' @export
skew_normal <- function(link = "identity") fam_skew_normal(link)

#' @rdname frmtmb-families
#' @export
exgaussian <- function(link = "identity") fam_exgaussian(link)

#' @rdname frmtmb-families
#' @export
bernoulli <- function(link = "logit") fam_bernoulli(link)

#' @rdname frmtmb-families
#' @export
geometric <- function(link = "log") fam_geometric(link)

#' @rdname frmtmb-families
#' @export
exponential <- function(link = "log") fam_exponential(link)

#' @rdname frmtmb-families
#' @export
weibull <- function(link = "log") fam_weibull(link)

#' @rdname frmtmb-families
#' @export
shifted_lognormal <- function(link = "identity") fam_shifted_lognormal(link)

#' @rdname frmtmb-families
#' @export
hurdle_gamma <- function(link = "log") fam_hurdle_gamma(link)

#' @rdname frmtmb-families
#' @export
hurdle_lognormal <- function(link = "identity") fam_hurdle_lognormal(link)

#' @rdname frmtmb-families
#' @export
zero_inflated_binomial <- function(link = "logit") fam_zi_binomial(link)

#' @rdname frmtmb-families
#' @export
zero_inflated_beta <- function(link = "logit") fam_zi_beta(link)

#' @rdname frmtmb-families
#' @export
asym_laplace <- function(link = "identity") fam_asym_laplace(link)

#' @rdname frmtmb-families
#' @export
zero_inflated_asym_laplace <- function(link = "identity") {
  fam_zi_asym_laplace(link)
}

#' @rdname frmtmb-families
#' @param k For `huber()`: Huber's tuning constant, the residual size in
#'   units of `sigma` where the density stops being gaussian and becomes
#'   Laplace. It is FIXED, not estimated - the default 1.345 is
#'   `MASS::rlm()`'s, which gives 95% efficiency against a gaussian.
#' @export
huber <- function(link = "identity", k = 1.345) fam_huber(link, k)

#' @rdname frmtmb-families
#' @export
sratio <- function(link = "logit") fam_sratio(link)

#' @rdname frmtmb-families
#' @export
cratio <- function(link = "logit") fam_cratio(link)

#' @rdname frmtmb-families
#' @export
acat <- function(link = "logit") fam_acat(link)

#' @rdname frmtmb-families
#' @export
von_mises <- function(link = "tan_half") fam_von_mises(link)

#' @rdname frmtmb-families
#' @param levels For `categorical()`: the response's category labels, in
#'   the order that fixes the reference category (the first) and the
#'   dpar names. Only needed when the family is built away from the
#'   data; [frm()] reads them off the response.
#' @export
categorical <- function(link = "logit", levels = NULL, K = NULL) {
  if (!identical(link, "logit")) {
    stop("categorical() supports the 'logit' link only", call. = FALSE)
  }
  if (!is.null(levels)) {
    levels <- as.character(levels)
    if (length(levels) < 2L || anyDuplicated(levels)) {
      stop("categorical(levels =): needs at least two distinct category ",
           "labels", call. = FALSE)
    }
    return(fam_categorical_impl(paste0("mu", levels[-1L]), levels, link))
  }
  if (!is.null(K)) {
    if (length(K) != 1L || is.na(K) || K < 2) {
      stop("categorical(K =): needs at least two categories",
           call. = FALSE)
    }
    K <- as.integer(K)
    return(fam_categorical_impl(paste0("mu", seq_len(K)[-1L]),
                                levels = NULL, link = link))
  }
  fam_categorical_deferred(link)
}

#' @rdname frmtmb-families
#' @param df For `cox()`: the number of M-spline basis functions in the
#'   baseline hazard (brms's `bhaz(df = 5)` default).
#' @param degree For `cox()`: the spline degree of that basis (cubic by
#'   default).
#' @param intercept For `cox()`: keep the basis function that is
#'   non-zero at the lower boundary knot.
#' @export
cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE) {
  # df and degree size the I-spline baseline basis, so a length-2 value
  # used to build a basis of one shape and record another
  check_count(df, "df", min = 1L)
  check_count(degree, "degree", min = 1L)
  check_flag(intercept, "intercept")
  fam_cox(link, df = df, degree = degree, intercept = intercept)
}

#' @export
print.frmtmb_family <- function(x, ...) {
  links <- vapply(x$links, function(l) l$name, "")
  cat("<frmtmb family> ", x$family, "\n", sep = "")
  cat("  dpars: ", paste0(x$dpars, " (", links[x$dpars], ")",
                          collapse = ", "), "\n", sep = "")
  invisible(x)
}
