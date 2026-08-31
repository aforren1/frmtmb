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

fam_gaussian <- function(link = "identity") {
  frmtmb_family(
    "gaussian",
    dpars = c("mu", "sigma"),
    links = list(mu = link, sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      RTMB::pnorm((q - dpars$mu) / dpars$sigma)
    },
    init_dpars = list(
      mu = function(y, aterms) mean(y),
      sigma = function(y, aterms) stats::sd(y)
    ),
    type = "continuous",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu,
      var_fn = function(dpars, aterms) dpars$sigma^2
    ),
    sim = function(dpars, aterms, n) stats::rnorm(n, dpars$mu, dpars$sigma)
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
      RTMB::dt((y - dpars$mu) / dpars$sigma, df = dpars$nu, log = TRUE) -
        log(dpars$sigma)
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
        ifelse(dpars$nu > 2, dpars$sigma^2 * dpars$nu / (dpars$nu - 2),
               NA_real_)
      }
    ),
    sim = function(dpars, aterms, n) {
      dpars$mu + dpars$sigma * stats::rt(n, dpars$nu)
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
  cumulative                = fam_cumulative
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

#' @export
print.frmtmb_family <- function(x, ...) {
  links <- vapply(x$links, function(l) l$name, "")
  cat("<frmtmb family> ", x$family, "\n", sep = "")
  cat("  dpars: ", paste0(x$dpars, " (", links[x$dpars], ")",
                          collapse = ", "), "\n", sep = "")
  invisible(x)
}
