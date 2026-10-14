# The bin-2P adapters: brms 2.23.0's internal per-family helpers
# log_lik_<fun>(i, prep) and posterior_predict_<fun>(i, prep), and its
# ordinal densities d<ordinal>(x, eta, thres, link), answered by
# frmtmb's own family objects. dev/brmsport-bin2-run.R stages this file
# as helper-brms-bin2-adapters.R beside the package's test helpers;
# dev/brmsport-bin2-adapters-spec.R says which blocks reach it.
#
# The adapter is THIN on purpose. A pass must be frmtmb's answer, so the
# adapter only moves numbers: it names the family as brms names it, reads
# column i of brms's prep, and hands frmtmb's row composition the
# result. The density is frmtmb's row_lpdf(), which the objective and
# log_lik() run (it folds cens() and trunc() into the family's lpdf),
# and the draw is frmtmb's sim_draw(), which every simulation entry
# point runs (it adds trunc() rejection to the family's sim). An ordinal
# family's density and draw are its exported ord_lpdf_make() and
# ord_sim_make(). Anything brms's prep carries that this mapping does
# not read is refused rather than dropped, so an assertion can never
# hold on a model the adapter silently simplified.
#
# One thing the adapter computes itself: brms's log_lik_weight(), the
# case weight times the row's log-likelihood. frmtmb applies weights
# outside the family too (in the objective and in log_lik() of a fit),
# so a row that asserts the weighted value tests this line and not
# frmtmb; its verdict says so (log_lik:32).

# A frmtmb dpar whose brms name differs, per frmtmb family name. brms
# calls com_poisson's dispersion `shape`; frmtmb's compois calls it `nu`.
brms_adapt_dpar_names <- list(compois = c(nu = "shape"))

# The prep fields a univariate, independent-row family reads, plus an
# ordinal family's thresholds (thres). The rest (ac: autocorrelation;
# mvpars: multivariate; lncor: logistic-normal) change what a row's
# density is.
brms_adapt_prep_fields <- c("dpars", "data", "family", "ndraws", "nobs",
                            "refcat", "thres")

# The ordinal families. Their exported family object builds its row
# density and its draw from given thresholds, ord_lpdf_make(tmap) and
# ord_sim_make(tmap), where tmap maps the free parameters to the
# thresholds; identity takes brms's thresholds as they are.
brms_adapt_ordinal <- c("cumulative", "sratio", "cratio", "acat",
                        "hurdle_cumulative")

# brms's addition-term data, under frmtmb's aterm names
brms_adapt_aterm_names <- c(trials = "trials", cens = "cens",
                            cens_y2 = "rcens", trunc_lb = "lb",
                            trunc_ub = "ub")
brms_adapt_data_fields <- c("Y", "weights", "ncat", brms_adapt_aterm_names)

# frmtmb's own refusal names a family it does not have, so the message
# a not-holding row records is frmtmb's. An ordinal family takes brms's
# link (prep$family$link: logit, probit, ...), the CDF of its latent
# variable.
brms_adapt_family <- function(fun, prep = NULL) {
  link <- prep[["family"]][["link"]]
  if (fun %in% brms_adapt_ordinal && !is.null(link)) {
    return(frmtmb::brmsfamily(fun, link = link))
  }
  frmtmb::brmsfamily(fun)
}

brms_adapt_check_prep <- function(prep, fun) {
  if (!is.null(prep[["thres"]]) && !fun %in% brms_adapt_ordinal) {
    stop("brms_adapt: prep$thres on a family that is not ordinal",
         call. = FALSE)
  }
  extra <- setdiff(names(prep), brms_adapt_prep_fields)
  if (length(extra)) {
    stop("brms_adapt: the prep carries ", paste(extra, collapse = ", "),
         ", which this adapter does not map", call. = FALSE)
  }
  extra <- setdiff(names(prep[["data"]]), brms_adapt_data_fields)
  if (length(extra)) {
    stop("brms_adapt: prep$data carries ", paste(extra, collapse = ", "),
         ", which this adapter does not map", call. = FALSE)
  }
}

# brms's get_dpar(): column i of a draws x observations matrix, else the
# per-draw vector or scalar as it is; recycled to the number of draws,
# which is prep$ndraws when brms set it and the longest dpar otherwise
brms_adapt_dpars <- function(fam, i, prep) {
  rn <- brms_adapt_dpar_names[[fam[["family"]]]]
  out <- list()
  for (d in fam[["dpars"]]) {
    b <- if (d %in% names(rn)) rn[[d]] else d
    x <- prep[["dpars"]][[b]]
    if (is.null(x)) {
      stop("brms_adapt: prep$dpars has no '", b, "', which frmtmb's ",
           fam[["family"]], " dpar '", d, "' reads", call. = FALSE)
    }
    out[[d]] <- if (is.null(dim(x))) x else x[, i]
  }
  n <- prep[["ndraws"]] %||% max(lengths(out))
  list(dpars = lapply(out, rep_len, length.out = n), n = n)
}

# brms reads every addition term at [i], as log_lik_censor() and
# log_lik_truncate() do, so a scalar bound is NA past the first row in
# brms and here alike.
#
# A term the frmtmb family does not declare in accepts_aterms: brms's
# helper for a family without trials never reads prep$data$trials (the
# zero-inflated block hands one prep to every family), so trials is
# dropped; brms's helpers DO apply cens and trunc, so those are refused
# by the family's own declaration rather than handed to a density that
# has no CDF to apply them with. An uncensored row (cens 0) is brms's
# plain density, which needs no CDF, so it is not refused.
brms_adapt_aterms <- function(fam, i, prep, n) {
  acc <- fam[["accepts_aterms"]]
  term <- c(trials = "trials", cens = "cens", cens_y2 = "cens",
            trunc_lb = "trunc", trunc_ub = "trunc")
  av <- list()
  for (k in names(brms_adapt_aterm_names)) {
    v <- prep[["data"]][[brms_adapt_aterm_names[[k]]]]
    if (is.null(v)) next
    if (!is.null(acc) && !term[[k]] %in% acc) {
      if (identical(k, "trials")) next
      if (identical(k, "cens") && isTRUE(v[i] == 0)) next
      stop("brms_adapt: frmtmb's ", fam[["family"]], " does not accept ",
           term[[k]], "() (accepts_aterms: ", paste(acc, collapse = ", "),
           ")", call. = FALSE)
    }
    av[[k]] <- rep_len(v[i], n)
  }
  av
}

# The ordinal path: the thresholds differ per draw, so the family's row
# function runs once per draw with that draw's dpars and thresholds
brms_adapt_ordinal_draws <- function(dp, prep, make, call) {
  th <- matrix(prep[["thres"]][["thres"]], nrow = dp$n)
  f <- make(identity)
  vapply(seq_len(dp$n), function(s) {
    as.numeric(call(f, s, lapply(dp$dpars, `[`, s), th[s, ]))
  }, numeric(1))
}

# One log-likelihood value per draw at observation i, weighted as brms's
# log_lik_weight() weights it (and as frmtmb's log_lik() does, outside
# the row density)
brms_adapt_log_lik <- function(fun, i, prep) {
  brms_adapt_check_prep(prep, fun)
  fam <- brms_adapt_family(fun, prep)
  dp <- brms_adapt_dpars(fam, i, prep)
  y <- rep_len(prep[["data"]][["Y"]][i], dp$n)
  if (fun %in% brms_adapt_ordinal) {
    ll <- brms_adapt_ordinal_draws(dp, prep, fam[["ord_lpdf_make"]],
                                   function(lp, s, dps, th) {
      lp(y[s], dps, list(), list(tau_raw = th))
    })
  } else {
    av <- brms_adapt_aterms(fam, i, prep, dp$n)
    ll <- frmtmb::row_lpdf(fam, y, y, dp$dpars, av, NULL)
  }
  w <- prep[["data"]][["weights"]]
  if (!is.null(w)) ll <- ll * w[i]
  as.numeric(ll)
}

# One predictive draw per posterior draw at observation i. brms ignores
# weights here, and so does frmtmb's simulator.
brms_adapt_posterior_predict <- function(fun, i, prep) {
  brms_adapt_check_prep(prep, fun)
  fam <- brms_adapt_family(fun, prep)
  dp <- brms_adapt_dpars(fam, i, prep)
  if (fun %in% brms_adapt_ordinal) {
    return(brms_adapt_ordinal_draws(dp, prep, fam[["ord_sim_make"]],
                                    function(sm, s, dps, th) {
      sm(dps, list(), 1L, list(tau_raw = th))
    }))
  }
  ctx <- list(family = fam, dpars = dp$dpars,
              aterms = brms_adapt_aterms(fam, i, prep, dp$n), n = dp$n,
              autocor = NULL, extra = NULL, max_iter = NULL)
  as.vector(frmtmb::sim_draw(ctx))
}

# brms's d<ordinal>(x, eta, thres, link): the probability of each
# category in x, one row per draw. eta is a vector (one linear predictor
# per draw) or a draws x thresholds matrix (category-specific effects).
# These four families read eta only through thres - eta per threshold,
# so a matrix eta is passed as mu = 0 with the thresholds thres - eta:
# the same numbers in frmtmb's arguments, not a different model.
brms_adapt_dord <- function(fun, x, eta, thres, link = "logit") {
  fam <- frmtmb::brmsfamily(fun, link = link)
  lp <- fam[["ord_lpdf_make"]](identity)
  thres <- as.matrix(thres)
  out <- matrix(NA_real_, nrow(thres), length(x))
  for (s in seq_len(nrow(thres))) {
    if (is.matrix(eta)) {
      mu <- 0
      th <- thres[s, ] - eta[s, ]
    } else {
      mu <- eta[s]
      th <- thres[s, ]
    }
    for (k in seq_along(x)) {
      out[s, k] <- exp(as.numeric(lp(x[k], list(mu = mu, disc = 1),
                                     list(), list(tau_raw = th))))
    }
  }
  out
}

# brms's tests source reference helpers from tests/testthat/helpers,
# which dev/brms-suite/ holds; a staged run's directory is
# dev/brmsport-bin2-run/<package>
brms_adapt_suite_dir <- local({
  root <- Sys.getenv("FRMTMB_PORT_ROOT", unset = "")
  if (nzchar(root)) {
    file.path(root, "dev/brms-suite/brms/tests/testthat")
  } else {
    normalizePath("../../brms-suite/brms/tests/testthat", winslash = "/",
                  mustWork = FALSE)
  }
})

# Stands in for brms's skip_if_not_installed() on a package that only
# brms's own implementation of the family needs (extraDistr behind
# brms's beta-binomial density, betareg behind its xbeta). frmtmb's
# families need neither, so the port runs; a REFERENCE density that
# needs the package still fails on its own line, and that row records
# the failure.
brms_adapt_brms_needs <- function(pkg) invisible(pkg)

# G1, fix round 2: brms's exported d and r functions of the families
# frmtmb has, answered by the same thin read as the ordinal densities.
# Each shim keeps brms's signature, maps its arguments onto the frmtmb
# family's dpars (and trials), and calls row_lpdf() or sim_draw(); the
# density and the draw are frmtmb's. Only pinv_gaussian() has a p shim:
# inverse.gaussian is the one family here with a CDF, reached as
# row_lpdf() of a left-censored row (log F). The other p, and every q,
# function is left unbound, so its row records "could not find
# function" and takes its verdict from dev/brmsport-bin2-verdicts-2p.R.
brms_g1_d <- function(fun, x, dp, av = list(), log = FALSE, cens = NULL) {
  fam <- frmtmb::brmsfamily(fun)
  n <- max(lengths(c(list(x), dp, av)))
  x <- rep_len(x, n)
  dp <- lapply(dp, rep_len, length.out = n)
  av <- lapply(av, rep_len, length.out = n)
  if (!is.null(cens)) av$cens <- rep_len(cens, n)
  ll <- as.numeric(frmtmb::row_lpdf(fam, x, x, dp, av, NULL))
  if (log) ll else exp(ll)
}
brms_g1_r <- function(fun, n, dp, av = list()) {
  fam <- frmtmb::brmsfamily(fun)
  dp <- lapply(dp, rep_len, length.out = n)
  av <- lapply(av, rep_len, length.out = n)
  as.vector(frmtmb::sim_draw(list(family = fam, dpars = dp, aterms = av,
                                  n = n, autocor = NULL, extra = NULL,
                                  max_iter = NULL)))
}
brms_g1_shims <- list(
  dstudent_t = function(x, df, mu = 0, sigma = 1, log = FALSE) {
    brms_g1_d("student", x, list(mu = mu, sigma = sigma, nu = df),
              log = log)
  },
  rstudent_t = function(n, df, mu = 0, sigma = 1) {
    brms_g1_r("student", n, list(mu = mu, sigma = sigma, nu = df))
  },
  dvon_mises = function(x, mu, kappa, log = FALSE) {
    brms_g1_d("von_mises", x, list(mu = mu, kappa = kappa), log = log)
  },
  rvon_mises = function(n, mu, kappa) {
    brms_g1_r("von_mises", n, list(mu = mu, kappa = kappa))
  },
  dskew_normal = function(x, mu = 0, sigma = 1, alpha = 0, log = FALSE) {
    brms_g1_d("skew_normal", x, list(mu = mu, sigma = sigma, alpha = alpha),
              log = log)
  },
  rskew_normal = function(n, mu = 0, sigma = 1, alpha = 0) {
    brms_g1_r("skew_normal", n, list(mu = mu, sigma = sigma, alpha = alpha))
  },
  dexgaussian = function(x, mu, sigma, beta, log = FALSE) {
    brms_g1_d("exgaussian", x, list(mu = mu, sigma = sigma, beta = beta),
              log = log)
  },
  rexgaussian = function(n, mu, sigma, beta) {
    brms_g1_r("exgaussian", n, list(mu = mu, sigma = sigma, beta = beta))
  },
  dinv_gaussian = function(x, mu = 1, shape = 1, log = FALSE) {
    brms_g1_d("inverse.gaussian", x, list(mu = mu, shape = shape),
              log = log)
  },
  pinv_gaussian = function(q, mu = 1, shape = 1, lower.tail = TRUE,
                           log.p = FALSE) {
    lp <- brms_g1_d("inverse.gaussian", q, list(mu = mu, shape = shape),
                    log = TRUE, cens = if (lower.tail) -1 else 1)
    if (log.p) lp else exp(lp)
  },
  rinv_gaussian = function(n, mu = 1, shape = 1) {
    brms_g1_r("inverse.gaussian", n, list(mu = mu, shape = shape))
  },
  dbeta_binomial = function(x, size, mu, phi, log = FALSE) {
    brms_g1_d("beta_binomial", x, list(mu = mu, phi = phi),
              list(trials = size), log = log)
  },
  rbeta_binomial = function(n, size, mu, phi) {
    brms_g1_r("beta_binomial", n, list(mu = mu, phi = phi),
              list(trials = size))
  },
  dasym_laplace = function(x, mu = 0, sigma = 1, quantile = 0.5,
                           log = FALSE) {
    brms_g1_d("asym_laplace", x, list(mu = mu, sigma = sigma,
                                      quantile = quantile), log = log)
  },
  rasym_laplace = function(n, mu = 0, sigma = 1, quantile = 0.5) {
    brms_g1_r("asym_laplace", n, list(mu = mu, sigma = sigma,
                                      quantile = quantile))
  },
  dzero_inflated_poisson = function(x, lambda, zi, log = FALSE) {
    brms_g1_d("zero_inflated_poisson", x, list(mu = lambda, zi = zi),
              log = log)
  },
  dzero_inflated_negbinomial = function(x, mu, shape, zi, log = FALSE) {
    brms_g1_d("zero_inflated_negbinomial", x,
              list(mu = mu, shape = shape, zi = zi), log = log)
  },
  dzero_inflated_binomial = function(x, size, prob, zi, log = FALSE) {
    brms_g1_d("zero_inflated_binomial", x, list(mu = prob, zi = zi),
              list(trials = size), log = log)
  },
  dzero_inflated_beta_binomial = function(x, size, mu, phi, zi,
                                          log = FALSE) {
    brms_g1_d("zero_inflated_beta_binomial", x,
              list(mu = mu, phi = phi, zi = zi), list(trials = size),
              log = log)
  },
  # brms's (shape1, shape2) is frmtmb's mean-precision (mu, phi)
  dzero_inflated_beta = function(x, shape1, shape2, zi, log = FALSE) {
    brms_g1_d("zero_inflated_beta", x,
              list(mu = shape1 / (shape1 + shape2), phi = shape1 + shape2,
                   zi = zi), log = log)
  },
  dhurdle_poisson = function(x, lambda, hu, log = FALSE) {
    brms_g1_d("hurdle_poisson", x, list(mu = lambda, hu = hu), log = log)
  },
  dhurdle_negbinomial = function(x, mu, shape, hu, log = FALSE) {
    brms_g1_d("hurdle_negbinomial", x, list(mu = mu, shape = shape,
                                            hu = hu), log = log)
  },
  # brms's (shape, scale) is frmtmb's mean-shape (mu = shape * scale)
  dhurdle_gamma = function(x, shape, scale, hu, log = FALSE) {
    brms_g1_d("hurdle_gamma", x, list(mu = shape * scale, shape = shape,
                                      hu = hu), log = log)
  },
  dhurdle_lognormal = function(x, mu, sigma, hu, log = FALSE) {
    brms_g1_d("hurdle_lognormal", x, list(mu = mu, sigma = sigma,
                                          hu = hu), log = log)
  })

# G2, fix 2b: brms's inv_link_<ordinal>(x, link) on a draws x
# observations x thresholds array, the category probabilities of each
# draw and observation. brms's x is thres - eta for cumulative and
# sratio and eta - thres for cratio and acat (as its own d<ordinal>()
# calls them), so each slice goes through brms_adapt_dord() with mu = 0
# and the thresholds x, or -x. The same ordinal family objects answer.
brms_adapt_ilord <- function(fun, x, link) {
  d <- dim(x)
  ncat <- d[3] + 1L
  sgn <- if (fun %in% c("cratio", "acat")) -1 else 1
  out <- array(NA_real_, c(d[1], d[2], ncat))
  for (j in seq_len(d[2])) {
    th <- matrix(sgn * x[, j, ], nrow = d[1])
    out[, j, ] <- brms_adapt_dord(fun, seq_len(ncat), rep(0, d[1]), th,
                                  link = link)
  }
  out
}
brms_g2_shims <- lapply(
  stats::setNames(nm = c("cumulative", "sratio", "cratio", "acat")),
  function(f) {
    force(f)
    function(x, link = "logit") brms_adapt_ilord(f, x, link)
  })
