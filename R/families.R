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
#' @param post Named list of numeric helper functions `mean_fn(dpars,
#'   aterms)` and `var_fn(dpars, aterms)`, used by [fitted()] and
#'   [residuals()].
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

# Residual SD including a known se() component (meta-analysis): se()
# alone replaces sigma; se(x, sigma = TRUE) adds them in quadrature.
resid_sd <- function(sigma, aterms) {
  if (is.null(aterms$se)) return(sigma)
  if (isTRUE(aterms$se_sigma)) sqrt(sigma^2 + aterms$se^2) else aterms$se
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
      var_fn = function(dpars, aterms) resid_sd(dpars$sigma, aterms)^2
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
      var_fn = function(dpars, aterms) dpars$mu
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
      var_fn = function(dpars, aterms) dpars$mu^2 / dpars$shape
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
      var_fn = function(dpars, aterms) dpars$mu + dpars$mu^2 / dpars$shape
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
      var_fn = function(dpars, aterms) dpars$mu * (1 + dpars$phi)
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
      var_fn = function(dpars, aterms) dpars$phi * dpars$mu^dpars$power
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
      var_fn = function(dpars, aterms) dpars$mu * (1 - dpars$mu)
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
      var_fn = function(dpars, aterms) dpars$mu * (1 + dpars$mu)
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
      var_fn = function(dpars, aterms) dpars$mu^2
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
      var_fn = function(dpars, aterms) dpars$mu^3 / dpars$shape
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
#' @param ... Two or more component families.
#' @return A `frmtmb_family`.
#' @export
mixture <- function(...) {
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
  frmtmb_family(
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
