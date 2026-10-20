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
# mvpars: multivariate) change what a row's density is.
#
# lncor is allowed at the top level because brms 2.23.0 does not read it
# there: get_Sigma(prep, cor_name = "lncor") reads prep$dpars$lncor,
# which is where the adapter reads it too (brms_adapt_cat_lncor()).
# brms's logistic-normal test blocks set prep$lncor, so brms evaluates
# them with uncorrelated log ratios (dev/bin2port-probe1.R shows brms's
# value does not move with prep$lncor and does with prep$dpars$lncor);
# reading prep$lncor here would evaluate a different model from brms's.
brms_adapt_prep_fields <- c("dpars", "data", "family", "ndraws", "nobs",
                            "refcat", "thres", "lncor")

# The ordinal families. Their exported family object builds its row
# density and its draw from given thresholds, ord_lpdf_make(tmap) and
# ord_sim_make(tmap), where tmap maps the free parameters to the
# thresholds; identity takes brms's thresholds as they are.
brms_adapt_ordinal <- c("cumulative", "sratio", "cratio", "acat",
                        "hurdle_cumulative")

# brms's addition-term data, under frmtmb's aterm names; dec is
# frmtmb.eam's dec() for the wiener family, brms's dec() too
brms_adapt_aterm_names <- c(trials = "trials", cens = "cens",
                            cens_y2 = "rcens", trunc_lb = "lb",
                            trunc_ub = "ub", dec = "dec")
brms_adapt_data_fields <- c("Y", "weights", "ncat", brms_adapt_aterm_names)

# frmtmb's own refusal names a family it does not have, so the message
# a not-holding row records is frmtmb's. An ordinal family takes brms's
# link (prep$family$link: logit, probit, ...), the CDF of its latent
# variable. brms's wiener family is frmtmb.eam's, whose namespace
# registers it, and its dec() term, with core when it loads.
brms_adapt_family <- function(fun, prep = NULL) {
  if (identical(fun, "wiener")) loadNamespace("frmtmb.eam")
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
            trunc_lb = "trunc", trunc_ub = "trunc", dec = "dec")
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

# The categorical-like families: one predictor per category. brms's
# prep holds them as mu1, mu2, ... (get_Mu() takes every dpar whose
# name starts with "mu", in the prep's order) and inserts the reference
# category's zero at prep$refcat (insert_refcat()); frmtmb's family
# names its predictors after the categories (mu2, mu3 with category 1
# the reference). The map is by position, as brms's own read is.
brms_adapt_catlike <- c("categorical", "multinomial", "dirichlet",
                        "dirichlet2", "logistic_normal")

# frmtmb's family as frm() builds it: the bare constructor brms's block
# names, resolved through its own `defer` slot on a response of
# prep$data$ncat categories, which is brms's count. The response only
# names the categories (1..K, as brms's seq_len(ncat)); its values are
# placeholders. A multinomial() without K is refused by frmtmb 0.71.0
# ("needs the number of categories"), and the row records that refusal.
brms_adapt_cat_family <- function(fun, prep) {
  K <- prep[["data"]][["ncat"]]
  if (is.null(K)) {
    stop("brms_adapt: prep$data has no ncat, the category count brms ",
         "reads", call. = FALSE)
  }
  Y <- prep[["data"]][["Y"]]
  if (is.matrix(Y) && ncol(Y) != K) {
    stop("brms_adapt: prep$data$Y has ", ncol(Y), " columns and ncat is ",
         K, call. = FALSE)
  }
  ref <- prep[["refcat"]]
  args <- list(fun)
  if (fun %in% c("dirichlet", "logistic_normal")) {
    # category k of 1..K is labeled "k", so the index names it
    args$refcat <- if (isTRUE(is.na(ref))) NA else as.character(ref)
  }
  if (identical(fun, "categorical") && !identical(as.integer(ref), 1L)) {
    stop("brms_adapt: frmtmb's categorical() takes the first category as ",
         "the reference; the prep has refcat = ", format(ref),
         call. = FALSE)
  }
  fam <- do.call(frmtmb::brmsfamily, args)
  if (!is.null(fam[["defer"]])) {
    d <- data.frame(row = seq_len(K))
    d$Y <- switch(fun,
                  categorical = seq_len(K),
                  multinomial = matrix(1L, K, K),
                  matrix(1 / K, K, K))
    fam <- fam[["defer"]](Y ~ 1, d)
  }
  fam
}

# The dpars of a categorical-like family at observation i: brms's mu*
# (and logistic_normal's sigma*) by position onto frmtmb's per-category
# names, every other dpar of the frmtmb family by its own name
brms_adapt_cat_dpars <- function(fam, i, prep) {
  sx <- fam[["simplex"]]
  pos <- function(pat, f_names) {
    b <- grep(pat, names(prep[["dpars"]]), value = TRUE)
    if (length(b) != length(f_names)) {
      stop("brms_adapt: the prep has ", length(b), " '", pat, "' dpars (",
           paste(b, collapse = ", "), ") and frmtmb's ", fam[["family"]],
           " has ", length(f_names), " (", paste(f_names, collapse = ", "),
           ")", call. = FALSE)
    }
    stats::setNames(b, f_names)
  }
  mu_f <- if (is.null(sx)) fam[["dpars"]] else sx[["mu"]]
  map <- pos("^mu", mu_f)
  if (!is.null(sx[["sigma"]])) map <- c(map, pos("^sigma", sx[["sigma"]]))
  out <- list()
  for (d in fam[["dpars"]]) {
    b <- if (d %in% names(map)) map[[d]] else d
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

# The logistic-normal correlations, per draw, as frmtmb's chart
# (lncor_raw, the row-normalized Cholesky chart of its correlation
# matrix). brms reads prep$dpars$lncor, a draws x D(D-1)/2 matrix filled
# as get_cov_matrix() fills it (for i in 2..D, j in 1..i-1), and none
# means uncorrelated. us_theta_cor() is frmtmb's own inverse of the
# chart its density reads, so this moves numbers into frmtmb's
# coordinates and nothing else.
brms_adapt_cat_lncor <- function(fam, prep, n) {
  D <- length(fam[["simplex"]][["pred"]])
  cors <- prep[["dpars"]][["lncor"]]
  chart <- function(C) frmtmb:::us_theta_cor(C)
  if (D < 2L) return(NULL)
  if (is.null(cors)) return(list(chart(diag(D))))
  cors <- matrix(cors, nrow = n)
  lapply(seq_len(n), function(s) {
    C <- diag(D)
    k <- 0L
    for (a in 2:D) for (b in seq_len(a - 1L)) {
      k <- k + 1L
      C[a, b] <- C[b, a] <- cors[s, k]
    }
    chart(C)
  })
}

# Run a categorical-like family's row function over the draws: once for
# all of them, or once per draw where the logistic-normal correlation
# differs per draw, with that draw's dpars and chart
brms_adapt_cat_run <- function(fam, prep, dp, f) {
  th <- if (identical(fam[["family"]], "logistic_normal")) {
    brms_adapt_cat_lncor(fam, prep, dp$n)
  }
  if (length(th) <= 1L) {
    return(f(dp$dpars, seq_len(dp$n), list(lncor_raw = th[[1L]])))
  }
  out <- lapply(seq_len(dp$n), function(s) {
    f(lapply(dp$dpars, `[`, s), s, list(lncor_raw = th[[s]]))
  })
  if (is.matrix(out[[1L]])) do.call(rbind, out) else unlist(out)
}

brms_adapt_cat_log_lik <- function(fun, i, prep) {
  fam <- brms_adapt_cat_family(fun, prep)
  dp <- brms_adapt_cat_dpars(fam, i, prep)
  Y <- prep[["data"]][["Y"]]
  av <- brms_adapt_aterms(fam, i, prep, dp$n)
  brms_adapt_cat_run(fam, prep, dp, function(dps, rows, extra) {
    y <- if (is.matrix(Y)) {
      matrix(Y[i, ], length(rows), ncol(Y), byrow = TRUE)
    } else {
      rep_len(Y[i], length(rows))
    }
    avs <- lapply(av, `[`, rows)
    as.numeric(frmtmb::row_lpdf(fam, y, y, dps, avs,
                                if (length(extra[[1L]])) extra))
  })
}

brms_adapt_cat_posterior_predict <- function(fun, i, prep) {
  fam <- brms_adapt_cat_family(fun, prep)
  dp <- brms_adapt_cat_dpars(fam, i, prep)
  av <- brms_adapt_aterms(fam, i, prep, dp$n)
  out <- brms_adapt_cat_run(fam, prep, dp, function(dps, rows, extra) {
    frmtmb::sim_draw(list(family = fam, dpars = dps,
                          aterms = lapply(av, `[`, rows), n = length(rows),
                          autocor = NULL,
                          extra = if (length(extra[[1L]])) extra,
                          max_iter = NULL))
  })
  # brms's draws carry no category names
  if (is.matrix(out)) unname(out) else as.vector(out)
}

# One log-likelihood value per draw at observation i, weighted as brms's
# log_lik_weight() weights it (and as frmtmb's log_lik() does, outside
# the row density)
brms_adapt_log_lik <- function(fun, i, prep) {
  brms_adapt_check_prep(prep, fun)
  if (fun %in% brms_adapt_catlike) {
    ll <- brms_adapt_cat_log_lik(fun, i, prep)
    w <- prep[["data"]][["weights"]]
    if (!is.null(w)) ll <- ll * w[i]
    return(ll)
  }
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
  if (fun %in% brms_adapt_catlike) {
    return(brms_adapt_cat_posterior_predict(fun, i, prep))
  }
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
  out <- as.vector(frmtmb::sim_draw(ctx))
  # the adapter answers brms's helper, and of brms's helpers
  # posterior_predict_wiener() alone returns its draws as an ndraws x 1
  # matrix (rcontinuous(types = "q")); the values are frmtmb.eam's
  if (identical(fun, "wiener")) out <- matrix(out, ncol = 1L)
  out
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

# brms 2.23.0 skips its wiener posterior_predict block outright ('skip
# as long as RWiener fails on R-devel for 3.6.0'): a fault of brms's
# RWiener backend, which frmtmb.eam's simulator does not use, so the
# port runs the block and the row records frmtmb.eam's answer
brms_adapt_brms_skip <- function(why) invisible(why)

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
  },
  # bin2port: brms's frechet functions take (loc, scale, shape); frmtmb's
  # frechet(), like brms's frechet family, takes the mean mu and the
  # shape nu (brms_g1_frechet_dp())
  dfrechet = function(x, loc = 0, scale = 1, shape = 1, log = FALSE) {
    brms_g1_d("frechet", x - loc, brms_g1_frechet_dp(scale, shape),
              log = log)
  },
  pfrechet = function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE,
                      log.p = FALSE) {
    brms_g1_p("frechet", q - loc, brms_g1_frechet_dp(scale, shape),
              lower.tail, log.p)
  },
  rfrechet = function(n, loc = 0, scale = 1, shape = 1) {
    loc + brms_g1_r("frechet", n, brms_g1_frechet_dp(scale, shape))
  },
  dgen_extreme_value = function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {
    brms_g1_d("gen_extreme_value", x, list(mu = mu, sigma = sigma, xi = xi),
              log = log)
  },
  pgen_extreme_value = function(q, mu = 0, sigma = 1, xi = 0,
                                lower.tail = TRUE, log.p = FALSE) {
    brms_g1_p("gen_extreme_value", q, list(mu = mu, sigma = sigma, xi = xi),
              lower.tail, log.p)
  },
  rgen_extreme_value = function(n, mu = 0, sigma = 1, xi = 0) {
    brms_g1_r("gen_extreme_value", n, list(mu = mu, sigma = sigma, xi = xi))
  })

# A p function through the family's CDF, as pinv_gaussian() reads it: a
# left-censored row is log F, a right-censored one log S
brms_g1_p <- function(fun, q, dp, lower.tail, log.p) {
  lp <- brms_g1_d(fun, q, dp, log = TRUE, cens = if (lower.tail) -1 else 1)
  if (log.p) lp else exp(lp)
}

# frmtmb's frechet() is parameterized by its mean, mu = scale *
# gamma(1 - 1 / nu), which is finite only for shape > 1. A shape at or
# below 1 has no mean, so there is no frmtmb parameter to hand the
# density: it is refused here rather than passed on as mu = Inf, where
# a length-only assertion would hold on NaN values.
brms_g1_frechet_dp <- function(scale, shape) {
  if (any(shape <= 1)) {
    stop("brms_g1: frmtmb's frechet() takes the mean mu = scale * ",
         "gamma(1 - 1 / nu), which is infinite at shape <= 1 (here ",
         "shape = ", paste(utils::head(shape[shape <= 1], 3), collapse = ", "),
         "); brms's d/p/r functions take the scale and answer there",
         call. = FALSE)
  }
  list(mu = scale * gamma(1 - 1 / shape), nu = shape)
}

# bin2port: brms's dcategorical(x, eta) and inv_link_categorical(x,
# refcat = 1), the category probabilities of each draw (and observation),
# answered by frmtmb's categorical family of K categories. frmtmb's
# reference category is the first, with its predictor at 0; the softmax
# is unchanged by a common shift, so brms's full eta (dcategorical() is
# called with refcat = NULL, every column given) becomes the predictors
# eta_k - eta_1. The family's row density at y = k is log P(Y = k).
brms_g2_cat_probs <- function(eta) {
  K <- ncol(eta)
  fam <- frmtmb::categorical(K = K)
  dp <- stats::setNames(lapply(seq_len(K)[-1L], function(k) {
    eta[, k] - eta[, 1L]
  }), fam[["dpars"]])
  n <- nrow(eta)
  vapply(seq_len(K), function(k) {
    y <- rep(k, n)
    exp(as.numeric(frmtmb::row_lpdf(fam, y, y, dp, list(), NULL)))
  }, numeric(n))
}
brms_g2_cat_shims <- list(
  dcategorical = function(x, eta, log = FALSE) {
    if (is.null(dim(eta))) eta <- matrix(eta, nrow = 1)
    p <- matrix(brms_g2_cat_probs(eta), nrow = nrow(eta))[, x, drop = FALSE]
    if (log) log(p) else p
  },
  # an array of draws x observations x (K - 1), category 1 the reference,
  # which is brms's default refcat = 1
  inv_link_categorical = function(x, log = FALSE, refcat = 1) {
    if (!identical(refcat, 1)) {
      stop("brms_g2: the shim takes brms's default refcat = 1 only",
           call. = FALSE)
    }
    d <- dim(x)
    out <- array(NA_real_, c(d[1], d[2], d[3] + 1L))
    for (j in seq_len(d[2])) {
      eta <- cbind(0, matrix(x[, j, ], nrow = d[1]))
      out[, j, ] <- brms_g2_cat_probs(eta)
    }
    if (log) log(out) else out
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

# Value checks (lane bin2port). brms wrote most of the blocks above as
# "runs without errors" and asserts a length or a dim; the value a line
# computed through the adapter is then never compared with anything.
# Where the line before a shape assertion computed a log density or a
# CDF, dev/brmsport-bin2-gen.R adds a check after the assertion: the
# same line, as brms wrote it, evaluated against brms's OWN function
# (brms:::log_lik_<family>(), or brms's d/p function under the name the
# header binds to a shim), on the same arguments, compared with what
# the adapter returned. The check is recorded as kind "value" under the
# brms row's id and is never part of that row's verdict.
#
# brms_value_mark() saves the random number state before the line, so a
# line that draws its arguments (dvon_mises(runif(n, -pi, pi), ...))
# draws the same ones for brms; brms_value_check() restores the state
# the block was in afterwards, so the rows after it see the stream they
# see without the check.
brms_value_state <- new.env(parent = emptyenv())
brms_value_mark <- function() {
  brms_value_state$seed <- get0(".Random.seed", envir = globalenv(),
                                inherits = FALSE)
  invisible(NULL)
}

# The reference environment: brms's own functions under the names the
# header binds to shims, and a prep whose family is brms's family of the
# same name where the block's family call reached frmtmb's constructor
# (brms's block calls brms's categorical(); the port calls frmtmb's)
brms_value_env <- function(env) {
  ref <- new.env(parent = env)
  ns <- asNamespace("brms")
  for (nm in c(names(brms_g1_shims), names(brms_g2_cat_shims))) {
    if (exists(nm, envir = env) && exists(nm, envir = ns)) {
      assign(nm, get(nm, envir = ns), envir = ref)
    }
  }
  prep <- get0("prep", envir = env)
  if (is.list(prep) && inherits(prep[["family"]], "frmtmb_family")) {
    prep$family <- brms::brmsfamily(prep[["family"]][["family"]])
    assign("prep", prep, envir = ref)
  }
  ref
}

brms_value_compare <- function(cur, ref) {
  cur <- as.numeric(unlist(cur))
  ref <- as.numeric(unlist(ref))
  if (length(cur) != length(ref)) {
    return(list(held = FALSE, msg = sprintf(
      "DIFFERS: %d values from the adapter, %d from brms", length(cur),
      length(ref))))
  }
  if (!identical(is.na(cur), is.na(ref))) {
    return(list(held = FALSE, msg = sprintf(
      "DIFFERS: NA in %d values from the adapter, %d from brms",
      sum(is.na(cur)), sum(is.na(ref)))))
  }
  ok <- !is.na(cur) & is.finite(ref)
  d <- abs(cur[ok] - ref[ok])
  same_inf <- identical(cur[!is.finite(ref) & !is.na(ref)],
                        ref[!is.finite(ref) & !is.na(ref)])
  rel <- if (length(d)) max(d / pmax(abs(ref[ok]), .Machine$double.xmin))
  # all.equal() is a MEAN relative difference, so one small entry far
  # off hides among large ones (review 2026-10-14, finding 2: 0.01 on an
  # entry of 0.263 among values to 5.8e6 held); every entry must also
  # be within testthat's tolerance of its own size, or of 1 below 1
  bound <- sqrt(.Machine$double.eps) * pmax(abs(ref[ok]), 1)
  held <- same_inf && isTRUE(all.equal(ref, cur)) && all(d <= bound)
  list(held = held, msg = sprintf(
    "%s: %d values (%d NA, %d infinite), max abs diff %.3g, max rel %.3g",
    if (held) "HOLDS" else "DIFFERS", length(cur), sum(is.na(cur)),
    sum(!is.finite(ref) & !is.na(ref)), if (length(d)) max(d) else 0,
    if (length(d)) rel else 0))
}

# `expected` is the outcome the ledger recorded ("holds", "differs" or
# "error"; "" before the first ledger run); a tier run asserts it. A
# check that does not hold also carries `msg_has`, the text its manual
# verdict was written for (dev/brmsport-bin2-valuechecks-manual.tsv),
# and the tier asserts that too: a defect fixed while the label stays
# "differs" must not keep its verdict.
#
# `adapted` is given where brms's call sits inside the assertion itself
# (expect_equal(length(brms:::log_lik_wiener(i, prep)), ns)): the
# adapter's value is then that call as the port runs it, evaluated here
# from the same random number state as brms's
brms_value_check <- function(id, var, expected, code, adapted = NULL,
                             msg_has = "") {
  call <- substitute(code)
  acall <- substitute(adapted)
  env <- parent.frame()
  now <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  res <- list(held = FALSE, vacuous = FALSE, msg = "", caught = "",
              raw_held = FALSE)
  if (is.null(acall)) {
    cur <- get0(var, envir = env, inherits = FALSE)
    stale <- brms_stale_reads(as.name(var), env)
  } else {
    if (!is.null(brms_value_state$seed)) {
      assign(".Random.seed", brms_value_state$seed, envir = globalenv())
    }
    a_env <- new.env(parent = env)
    cur <- tryCatch(suppressWarnings(suppressMessages({
      eval(acall, a_env)
      get(var, envir = a_env, inherits = FALSE)
    })), error = function(e) NULL)
    stale <- brms_stale_reads(acall, env)
  }
  if (length(stale) || is.null(cur)) {
    res$msg <- "NOT CHECKED: the adapter's value does not exist"
    res$vacuous <- TRUE
  } else {
    ref_env <- brms_value_env(env)
    if (!is.null(brms_value_state$seed)) {
      assign(".Random.seed", brms_value_state$seed, envir = globalenv())
    }
    out <- tryCatch(
      suppressWarnings(suppressMessages({
        eval(call, ref_env)
        list(ok = TRUE, value = get(var, envir = ref_env, inherits = FALSE))
      })),
      error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    if (out$ok) {
      cmp <- brms_value_compare(cur, out$value)
      res$held <- cmp$held
      res$raw_held <- cmp$held
      res$msg <- cmp$msg
    } else {
      res$msg <- brms_port_clean(paste("REFERENCE ERROR:", out$msg))
    }
  }
  if (!is.null(now)) assign(".Random.seed", now, envir = globalenv())
  brms_value_state$seed <- NULL
  outcome <- if (res$held) "holds" else if (grepl("^DIFFERS", res$msg)) {
    "differs"
  } else {
    "error"
  }
  if (nzchar(Sys.getenv("FRMTMB_BRMSPORT_RECORD"))) {
    brms_port_record("value", id, expected, res)
    testthat::succeed()
  } else if (!identical(outcome, expected)) {
    testthat::fail(sprintf(
      "brms %s value check is recorded as '%s' and is now '%s': %s", id,
      expected, outcome, res$msg))
  } else if (nzchar(msg_has) && !grepl(msg_has, res$msg, fixed = TRUE)) {
    # the same label for a different reason: the manual verdict was
    # written for the message it names, so it is stale
    testthat::fail(sprintf(paste(
      "brms %s value check is still '%s', but its message no longer",
      "contains '%s', for which its verdict was written: %s"), id,
      outcome, msg_has, res$msg))
  } else {
    testthat::succeed()
  }
  invisible(res)
}
