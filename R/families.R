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
#' @param type One of `"continuous"`, `"discrete"`, `"ordinal"`.
#' @param post Named list of numeric helper functions used by
#'   [fitted()], [predict()] and [residuals()]: `mean_fn(dpars,
#'   aterms)` (the response mean), `var_fn(dpars, aterms)` (for pearson
#'   residuals) and `dev_fn(y, dpars, aterms)` (the unit deviance
#'   `2 * (loglik of the saturated fit - loglik at the fitted value)`,
#'   for `residuals(type = "deviance")`). A family that omits one is
#'   refused by the method that needs it.
#' @param sim Optional numeric simulator `(dpars, aterms, n)` returning `n`
#'   response draws; used by [simulate()].
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
#' @return An object of class `frmtmb_family`.
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
#' @export
frmtmb_family <- function(family, dpars, links, lpdf, valid_y = NULL,
                          init_dpars = list(), type = "continuous",
                          post = list(), sim = NULL,
                          primary_dpars = "mu", lcdf = NULL,
                          extra_pars = NULL, drop_intercept = FALSE) {
  stopifnot(is.character(family), length(family) == 1,
            is.character(dpars), length(dpars) >= 1,
            is.function(lpdf))
  if (!all(primary_dpars %in% dpars)) {
    stop("`primary_dpars` must be a subset of `dpars`", call. = FALSE)
  }
  if (!setequal(names(links), dpars)) {
    stop("`links` must name every dpar exactly once", call. = FALSE)
  }
  links <- lapply(links, get_link)
  structure(
    list(family = family, dpars = dpars, links = links, lpdf = lpdf,
         valid_y = valid_y, init_dpars = init_dpars, type = type,
         post = post, sim = sim, primary_dpars = primary_dpars,
         lcdf = lcdf, extra_pars = extra_pars,
         drop_intercept = isTRUE(drop_intercept)),
    class = "frmtmb_family"
  )
}

#' @rdname frmtmb_family
#' @export
custom_family <- frmtmb_family

positive_y <- function(name) {
  force(name)
  function(y, aterms) {
    if (any(y <= 0)) {
      stop(name, ": response must be strictly positive", call. = FALSE)
    }
  }
}

count_y <- function(name) {
  force(name)
  function(y, aterms) {
    if (any(y < 0) || any(y != round(y))) {
      stop(name, ": response must be non-negative integers", call. = FALSE)
    }
  }
}

# P(a < Z < b) for standard normal bounds. pnorm(b) - pnorm(a) loses
# every significant digit when both bounds sit in the upper tail, which
# is exactly where truncated means are computed; the mirrored form keeps
# full precision there.
pnorm_diff <- function(a, b) {
  ifelse(a > 0,
         stats::pnorm(-a) - stats::pnorm(-b),
         stats::pnorm(b) - stats::pnorm(a))
}

# Residual SD including a known se() component (meta-analysis): se()
# alone replaces sigma; se(x, sigma = TRUE) adds them in quadrature.
resid_sd <- function(sigma, aterms) {
  if (is.null(aterms$se)) return(sigma)
  if (isTRUE(aterms$se_sigma)) sqrt(sigma^2 + aterms$se^2) else aterms$se
}

# trunc() bounds from an aterm-value list as a pair of length-n numeric
# vectors, or NULL when the response is not truncated. An absent bound
# is the family's unbounded default; the per-family truncated means all
# handle the infinities.
trunc_bounds <- function(aterms, n) {
  lb <- aterms[["trunc_lb"]]
  ub <- aterms[["trunc_ub"]]
  if (is.null(lb) && is.null(ub)) return(NULL)
  list(lb = rep(lb %||% -Inf, length.out = n),
       ub = rep(ub %||% Inf, length.out = n))
}

# Expected response at the given dpar values. On a truncated response
# this is E[Y | lb <= Y <= ub], the quantity fitted(), residuals() and
# predict(type = "response") report; per-dpar predictions stay
# untruncated because they describe the latent parameter.
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

# y log(y / mu) under the 0 log 0 = 0 convention the saturated fit needs
# at a zero count.
ylogy_mu <- function(y, mu) ifelse(y > 0, y * log(y / mu), 0)

# Negative-binomial unit deviance at a given size (shape). nbinom1 feeds
# it the row's own size mu / phi.
nbinom_deviance <- function(y, mu, size) {
  2 * (ylogy_mu(y, mu) - (y + size) * log((y + size) / (mu + size)))
}

# Binomial unit deviance on COUNTS out of `size` trials; both terms
# vanish at the boundaries y = 0 and y = size.
binomial_deviance <- function(y, mu, size) {
  2 * (ylogy_mu(y, size * mu) + ylogy_mu(size - y, size * (1 - mu)))
}

# Gamma unit deviance; also the exponential one (shape fixed at 1).
gamma_deviance <- function(y, mu) 2 * ((y - mu) / mu - log(y / mu))

# Families that define a unit deviance, for the residuals() refusal
# message. Read off the registry so the list cannot drift; the
# constructors that need arguments (multinomial) simply drop out.
deviance_family_names <- function() {
  nms <- names(family_registry)
  ok <- vapply(nms, function(nm) {
    fam <- tryCatch(family_registry[[nm]](), error = function(e) NULL)
    !is.null(fam) && !is.null(fam$post$dev_fn)
  }, TRUE)
  sort(unique(nms[ok]))
}

# Deviance residuals: sign(y - E[Y]) * sqrt(w_i d_i). Weights multiply
# the unit deviance, as in glm(). A truncated or censored response is
# refused: the fitted likelihood is not the family's own density there,
# so the saturated comparison the unit deviance is built on does not
# describe the model that was estimated.
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

# Per-observation subset of dpar / aterm vectors, for resampling the
# rows a rejection step has not accepted yet.
subset_obs <- function(x, idx, n) {
  lapply(x, function(v) {
    if (is.numeric(v) && length(v) %in% c(1L, n)) {
      rep(v, length.out = n)[idx]
    } else {
      v
    }
  })
}

# One simulated response vector, respecting trunc() bounds by rejection:
# out-of-bounds draws are redrawn until every row is inside its own
# interval. Bounds that exclude nearly all the family's mass never
# converge, so the iteration cap reports the acceptance rate instead of
# spinning.
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

fam_binomial <- function(link = "logit") {
  frmtmb_family(
    "binomial",
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms) {
      size <- aterms$trials %||% 1
      RTMB::dbinom(y, size, dpars$mu, log = TRUE)
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

fam_negbinomial <- function(link = "log") {
  frmtmb_family(
    "negbinomial",
    dpars = c("mu", "shape"),
    links = list(mu = link, shape = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnbinom2(y, dpars$mu, dpars$mu + dpars$mu^2 / dpars$shape,
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

fam_nbinom1 <- function(link = "log") {
  frmtmb_family(
    "nbinom1",
    dpars = c("mu", "phi"),
    links = list(mu = link, phi = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnbinom2(y, dpars$mu, dpars$mu * (1 + dpars$phi), log = TRUE)
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

fam_beta <- function(link = "logit") {
  frmtmb_family(
    "beta",
    dpars = c("mu", "phi"),
    links = list(mu = link, phi = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dbeta(y, dpars$mu * dpars$phi, (1 - dpars$mu) * dpars$phi,
                  log = TRUE)
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

fam_zi_poisson <- function(link = "log") {
  frmtmb_family(
    "zero_inflated_poisson",
    dpars = c("mu", "zi"),
    links = list(mu = link, zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      # y == 0 is data, so the mixture stays branch-free in parameters
      i0 <- as.numeric(y == 0)
      p0 <- exp(-dpars$mu)
      i0 * log(dpars$zi + (1 - dpars$zi) * p0) +
        (1 - i0) * (log(1 - dpars$zi) + RTMB::dpois(y, dpars$mu, log = TRUE))
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

fam_zi_negbinomial <- function(link = "log") {
  frmtmb_family(
    "zero_inflated_negbinomial",
    dpars = c("mu", "shape", "zi"),
    links = list(mu = link, shape = "log", zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      p0 <- exp(dpars$shape * (log(dpars$shape) -
                                 log(dpars$shape + dpars$mu)))
      base <- RTMB::dnbinom2(y, dpars$mu,
                             dpars$mu + dpars$mu^2 / dpars$shape,
                             log = TRUE)
      i0 * log(dpars$zi + (1 - dpars$zi) * p0) +
        (1 - i0) * (log(1 - dpars$zi) + base)
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

fam_hurdle_poisson <- function(link = "log") {
  frmtmb_family(
    "hurdle_poisson",
    dpars = c("mu", "hu"),
    links = list(mu = link, hu = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      # nonzero part is a zero-truncated poisson
      i0 * log(dpars$hu) +
        (1 - i0) * (log(1 - dpars$hu) +
                      RTMB::dpois(y, dpars$mu, log = TRUE) -
                      log(1 - exp(-dpars$mu)))
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

fam_bernoulli <- function(link = "logit") {
  frmtmb_family(
    "bernoulli",
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms) {
      RTMB::dbinom(y, 1, dpars$mu, log = TRUE)
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

fam_geometric <- function(link = "log") {
  frmtmb_family(
    "geometric",
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnbinom2(y, dpars$mu, dpars$mu * (1 + dpars$mu), log = TRUE)
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

fam_hurdle_gamma <- function(link = "log") {
  frmtmb_family(
    "hurdle_gamma",
    dpars = c("mu", "shape", "hu"),
    links = list(mu = link, shape = "log", hu = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      yp <- y + i0   # dodge dgamma(0) = -Inf; the term carries weight 0
      i0 * log(dpars$hu) +
        (1 - i0) * (log(1 - dpars$hu) +
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

fam_hurdle_lognormal <- function(link = "identity") {
  frmtmb_family(
    "hurdle_lognormal",
    dpars = c("mu", "sigma", "hu"),
    links = list(mu = link, sigma = "log", hu = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      yp <- y + i0
      i0 * log(dpars$hu) +
        (1 - i0) * (log(1 - dpars$hu) +
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

fam_zi_binomial <- function(link = "logit") {
  frmtmb_family(
    "zero_inflated_binomial",
    dpars = c("mu", "zi"),
    links = list(mu = link, zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      size <- aterms$trials %||% 1
      i0 <- as.numeric(y == 0)
      p0 <- (1 - dpars$mu)^size
      i0 * log(dpars$zi + (1 - dpars$zi) * p0) +
        (1 - i0) * (log(1 - dpars$zi) +
                      RTMB::dbinom(y, size, dpars$mu, log = TRUE))
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

fam_zi_beta <- function(link = "logit") {
  frmtmb_family(
    "zero_inflated_beta",
    dpars = c("mu", "phi", "zi"),
    links = list(mu = link, phi = "log", zi = "logit"),
    lpdf = function(y, dpars, aterms) {
      i0 <- as.numeric(y == 0)
      ya <- y + i0 * 0.5   # dodge dbeta(0) = -Inf; term carries weight 0
      i0 * log(dpars$zi) +
        (1 - i0) * (log(1 - dpars$zi) +
                      RTMB::dbeta(ya, dpars$mu * dpars$phi,
                                  (1 - dpars$mu) * dpars$phi, log = TRUE))
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
      log(p) + log(1 - p) - log(dpars$sigma) -
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

# --- RTMBdist-backed families ---

fam_beta_binomial <- function(link = "logit") {
  frmtmb_family(
    "beta_binomial",
    dpars = c("mu", "phi"),
    links = list(mu = link, phi = "log"),
    lpdf = function(y, dpars, aterms) {
      size <- aterms$trials %||% 1
      RTMBdist::dbetabinom(y, size, dpars$mu * dpars$phi,
                           (1 - dpars$mu) * dpars$phi, log = TRUE)
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

# brms parameterization: mu is the DISTRIBUTION mean, beta the scale of
# the exponential component (gaussian component sits at mu - beta).
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
# During oneStepPredict RTMB hands the response to the lpdf as an "osa"
# object: the taped value in @x and the per-row data-term indicator in
# @keep. RTMB's own densities apply the indicator through dGenericOSA,
# so a hand-written lpdf has to do it itself or every observation stays
# switched on and the one-step sequence collapses.
osa_unwrap <- function(y) {
  if (!methods::is(y, "osa")) return(NULL)
  keep <- y@keep
  if (ncol(keep) != 1L) {
    stop("osa_method = \"cdf\" is not supported for this family",
         call. = FALSE)
  }
  list(y = y@x, keep = keep[, 1])
}

ord_cat_sel <- function(y, K) {
  lapply(seq_len(K), function(k) {
    s <- 1
    for (j in seq_len(K)) {
      if (j != k) s <- s * ((y - j) / (k - j))
    }
    s
  })
}

# 1{j < y} for j = 1..K-1, from the same basis.
ord_cat_below <- function(sel, K) {
  lapply(seq_len(K - 1L), function(j) {
    b <- sel[[j + 1L]]
    if (j + 1L < K) for (k in (j + 2L):K) b <- b + sel[[k]]
    b
  })
}

# Cumulative ordinal: response in 1..K (or an ordered factor). The linear
# predictor has no intercept; K-1 ordered thresholds take its place,
# parameterized as (tau_1, log increments) in `extra_pars`.
fam_cumulative <- function(link = "logit") {
  Fcdf <- switch(link,
    logit = function(x) 1 / (1 + exp(-x)),
    probit = function(x) RTMB::pnorm(x),
    stop("cumulative() supports links 'logit' and 'probit'", call. = FALSE)
  )
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
      up <- Fcdf(tau[pmin(y, K1)] - eta) * (1 - iK) + iK
      lo <- Fcdf(tau[pmax(y - 1, 1)] - eta) * (1 - i1)
      log(up - lo)
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

# Shared scaffolding for the sequential ordinal families: an
# n x (K-1) matrix of (tau_j - eta_i) or (eta_i - tau_j), and data-only
# indicator matrices selecting the observed category (branch-free).
ord_indicators <- function(y, K1) {
  n <- length(y)
  jj <- rep(seq_len(K1), each = n)
  yy <- rep(y, K1)
  list(
    sel = matrix(as.numeric(yy == jj), n, K1),   # j == y (y <= K-1)
    below = matrix(as.numeric(jj < yy), n, K1)   # j < y
  )
}

ord_eta_mat <- function(eta, tau, n, K1) {
  # broadcast via matmul: rep() strips the advector class
  TM <- RTMB::matrix(1, n, 1) %*% RTMB::matrix(tau, 1, K1)
  TM - eta   # column-wise recycling: row i is tau_j - eta_i
}

# Sequential (sratio/cratio) log-density with the response on the tape.
# Column-at-a-time so nothing indexes an advector and nothing needs the
# n x (K-1) matrices, which the data path builds for speed.
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

ord_valid_y <- function(name) {
  function(y, aterms) {
    if (any(y < 1) || any(y != round(y)) || length(unique(y)) < 2) {
      stop(name, ": response must be an ordered factor or integers ",
           "1..K with at least 2 observed categories", call. = FALSE)
    }
  }
}

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

# cumulative and sratio store ordered thresholds as (tau_1, log
# increments); cratio and acat store them raw.
ord_tau_from_raw <- function(raw, ordered) {
  if (!ordered || length(raw) < 2L) return(raw)
  c(raw[1L], raw[1L] + cumsum(exp(raw[-1L])))
}

ord_num_cdf <- function(link) {
  switch(link, logit = stats::plogis, probit = stats::pnorm,
         stop("no numeric CDF for link '", link, "'", call. = FALSE))
}

# n x K matrix of category probabilities.
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

ord_link_cdf <- function(name, link) {
  switch(link,
    logit = function(x) 1 / (1 + exp(-x)),
    probit = function(x) RTMB::pnorm(x),
    stop(name, "() supports links 'logit' and 'probit'", call. = FALSE)
  )
}

# Stopping ratio (brms sratio): P(y=k) = F(tau_k - eta) *
# prod_{j<k} (1 - F(tau_j - eta)); ordered thresholds like cumulative.
fam_sratio <- function(link = "logit") {
  Fcdf <- ord_link_cdf("sratio", link)
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

# Continuation ratio (brms cratio): P(y=k) = (1 - F(eta - tau_k)) *
# prod_{j<k} F(eta - tau_j); unordered thresholds.
fam_cratio <- function(link = "logit") {
  Fcdf <- ord_link_cdf("cratio", link)
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

# Adjacent category (brms acat, logit link): P(y=k) proportional to
# exp(sum_{j<k} (eta - tau_j)); unordered thresholds.
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
        den <- 0
        for (r in seq_len(K)) {
          Er <- (r - 1) * eta - ct0[r]
          if (!is.null(dpars$.cs) && r >= 2L) {
            acc <- acc + dpars$.cs[, r - 1L]
            Er <- Er + acc
          }
          num <- num + sel[[r]] * Er
          den <- den + exp(Er)
        }
        return((num - log(den)) * ov$keep)
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
      as.vector((E * S) %*% ones) - log(as.vector(exp(E) %*% ones))
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
#' remains approximate - use [check_laplace()] to judge.
#' Mixing-weight predictors are evaluated at each group's first row
#' (use group-constant covariates). [mixture_probs()] returns the
#' posterior class probabilities per group (or per observation for
#' ordinary mixtures), conditional on the random-effect modes.
#'
#' @param ... Two or more component families.
#' @param groups Optional one-sided formula naming the latent-class
#'   grouping factor.
#' @return A `frmtmb_family`.
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
        if (is.null(sk)) stop("Component '", comps[[k]]$family,
                              "' has no simulator", call. = FALSE)
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
        stop("Component '", comps[[k]]$family, "' has no simulator",
             call. = FALSE)
      }
      sk(dpars_k, aterms, n)
    },
    log_pi = log_pi
  )
  if (!is.null(groups)) {
    if (!inherits(groups, "formula") || length(groups) != 2L) {
      stop("groups must be a one-sided formula: groups = ~g",
           call. = FALSE)
    }
    fam$mix_groups <- groups
  }
  fam
}

#' Posterior class probabilities of a mixture fit
#'
#' For an ordinary [mixture()] or [mixture_mvn()] fit, one row per
#' observation; for a group-level mixture (`groups = ~g`), one row per
#' group.
#'
#' @param fit A `frmtmb_fit` with a mixture family.
#' @return A matrix of class probabilities (rows sum to one).
#' @export
mixture_probs <- function(fit) {
  rspec <- uni_resp(fit, "mixture_probs()")
  fam <- rspec$family
  if (is.null(fam$mix)) {
    stop("mixture_probs() needs a mixture() family fit", call. = FALSE)
  }
  dp <- eval_dpars(fit)[[rspec$resp_name]]
  av <- fit$frame$aterm_values[[rspec$resp_name]]
  yv <- fit$frame$y[[rspec$resp_name]]
  lps_pi <- fam$mix$log_pi(dp)
  K <- fam$mix$K
  mg <- fit$frame$mix_g[[rspec$resp_name]]
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

# The extras a covariance model needs and how they assemble class k's
# covariance. `pars` is the (name, length) template used for the error
# messages and the documentation; `init` turns the response into the
# starting values; `sigma` runs on the tape, so it must stay AD-safe.
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
  fam
}

# Matrix-response multinomial: y is an n x K count matrix, category 1 is
# the reference. One linear predictor per non-reference category (mu2,
# ..., muK), all receiving the main model formula unless overridden.
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
      denom <- 1
      for (k in dpn) denom <- denom + exp(dpars[[k]])
      nvec <- rowSums(y)
      ll <- -nvec * log(denom)
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
  sratio                    = fam_sratio,
  cratio                    = fam_cratio,
  acat                      = fam_acat
)

as_frmtmb_family <- function(x) {
  if (inherits(x, "frmtmb_family")) return(x)
  if (is.function(x)) x <- x()
  if (inherits(x, "family")) {
    ctor <- family_registry[[x$family]]
    if (is.null(ctor)) {
      stop("Unsupported family: '", x$family, "'. Currently supported: ",
           paste(unique(names(family_registry)), collapse = ", "),
           call. = FALSE)
    }
    return(ctor(link = x$link))
  }
  if (is.character(x) && length(x) == 1) {
    ctor <- family_registry[[x]]
    if (is.null(ctor)) {
      stop("Unsupported family: '", x, "'. Currently supported: ",
           paste(unique(names(family_registry)), collapse = ", "),
           call. = FALSE)
    }
    return(ctor())
  }
  stop("Cannot interpret `family` of class ",
       paste(class(x), collapse = "/"), call. = FALSE)
}

#' Additional response families
#'
#' Family constructors without a [stats::family] equivalent, following
#' brms naming. `gaussian()`, `poisson()`, `binomial()`, and `Gamma()`
#' from 'stats' are accepted directly by [frm()] and `+`.
#'
#' @param link Link for `mu`.
#' @return A `frmtmb_family` object.
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
sratio <- function(link = "logit") fam_sratio(link)

#' @rdname frmtmb-families
#' @export
cratio <- function(link = "logit") fam_cratio(link)

#' @rdname frmtmb-families
#' @export
acat <- function(link = "logit") fam_acat(link)

#' @export
print.frmtmb_family <- function(x, ...) {
  links <- vapply(x$links, function(l) l$name, "")
  cat("<frmtmb family> ", x$family, "\n", sep = "")
  cat("  dpars: ", paste0(x$dpars, " (", links[x$dpars], ")",
                          collapse = ", "), "\n", sep = "")
  invisible(x)
}
