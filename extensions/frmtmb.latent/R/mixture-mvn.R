# Multivariate normal mixtures: model-based clustering of a matrix
# response, in the manner of mclust.
#
# This family lived in frmtmb until it moved here, beside lca() and
# hmm(): all three are models whose latent variable is a discrete class
# or state. It was already written against the family protocol that
# frmtmb exports, with one exception, the unstructured covariance
# assembler, which frmtmb now exports as us_sigma() for any family
# with a free covariance.

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
    frm_stop("mixture_mvn(): unknown covariance model '",
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
#' response, in the manner of mclust: K classes, each with its own
#' D-dimensional mean and a D x D covariance from mclust's model
#' taxonomy. Every class mean is a full linear predictor, and the main
#' model formula applies to all of them, so cluster means may depend on
#' covariates, which mclust cannot do. The location dpars are named
#' `mu<k>d<j>` (class k, response column j). Each one except the first,
#' `mu1d1`, takes a formula of its own, for example
#' `bf(Y ~ x, mu2d1 ~ 1)`. The mixing weights are `theta1 ...
#' theta{K-1}`, a multinomial logit against class K, each with its own
#' linear predictor. So gating on covariates works as in
#' [frmtmb::mixture()].
#'
#' brms has no multivariate normal mixture. Its nearest spelling,
#' `bf(mvbind(y1, y2) ~ 1)` with `family = mixture(gaussian, gaussian)`,
#' gives each response a mixture with a class variable of its own, which
#' is a different model; frmtmb fits that one in core.
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
#' entries, as in the `us()` covariance structure of [frmtmb::frm()];
#' [frmtmb::us_sigma()] assembles it. The `sigmashape` vectors hold the
#' first D - 1 log-shape entries; the last is minus their sum, which is
#' what fixes `det(A_k) = 1`. Log-SDs start at the per-column response
#' SDs and correlations at zero; class means start on spread-out
#' per-column response quantiles to break the label symmetry. The usual
#' finite-mixture ML caveats apply: the likelihood is invariant to
#' relabeling and can be multimodal (compare starts with
#' [frmtmb::frm_allfit()]).
#'
#' [frmtmb::mixture_probs()] returns the posterior class probabilities
#' per row, and `fitted()` returns the n x D mixture-mean matrix.
#' `simulate()`, `posterior_predict()` and [frmtmb::frm_simulate()]
#' draw a class per row from its mixing weights, then a D-variate
#' normal about that class's mean; a draw is an n x D matrix. A draw at
#' `newdata` is refused. Covariances take no linear predictor (no
#' covariance regression), and the models with a class-varying
#' eigenvector basis (`EEV`, `VEV`, `EVE`, `VEE`, `VVE`, `EVV`) are not
#' available. `cens()`, `trunc()`, [frmtmb::mvbf()], REML and
#' `profile = TRUE` are refused.
#'
#' @param K Number of mixture classes (at least 2).
#' @param D Number of response columns (at least 2; for `D = 1` use
#'   `mixture(gaussian(), ...)`).
#' @param model Covariance model name from mclust's vocabulary: one of
#'   `"EII"`, `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`,
#'   `"VVV"` (the default, a free covariance per class).
#' @return A `frmtmb_family`.
#' @seealso [lca()] for latent classes over categorical items,
#'   [frmtmb::mixture()] for a mixture of univariate families.
#' @examples
#' set.seed(1)
#' Y <- rbind(matrix(rnorm(60, 0), ncol = 2),
#'            matrix(rnorm(60, 4), ncol = 2))
#' dd <- data.frame(row = seq_len(nrow(Y)))
#' dd$Y <- Y
#' fit <- frm(bf(Y ~ 1), family = mixture_mvn(K = 2, D = 2), data = dd)
#' fixef(fit)
#' head(mixture_probs(fit))
#' # a shared spherical covariance (mclust's EII, like k-means)
#' frm(bf(Y ~ 1), family = mixture_mvn(K = 2, D = 2, model = "EII"),
#'     data = dd)
#'
#' # the alternative spelling adds the family to the formula
#' frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
#' @export
mixture_mvn <- function(K, D, model = "VVV") {
  if (missing(K) || missing(D) || K < 2 || D < 2) {
    frm_stop("mixture_mvn() needs K >= 2 classes and D >= 2 response ",
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
    accepts_aterms = "weights",
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
        frm_stop("mixture_mvn(K = ", K, ", D = ", D, "): response must be ",
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
      },
      # the same softmax reporting scale mixture() gives its mixing
      # weights, for the same reason
      dpar_response = list(
        dpars = paste0("theta", seq_len(K - 1L)),
        value = function(dpars, dnm) {
          exp(log_pi(dpars)[[as.integer(sub("^theta", "", dnm))]])
        },
        deriv = function(dpars, dnm) {
          p <- exp(log_pi(dpars)[[as.integer(sub("^theta", "", dnm))]])
          p * (1 - p)
        }
      )
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
  fam[["mix"]] <- list(
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
  fam[["sim_ctx"]] <- mvn_sim_rows
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

## ---- the compatibility matrix's mixture_mvn() rows -------------------

# These rows came from frmtmb's R/compat.R with the family, unchanged,
# for the reason hmm()'s and lca()'s live beside theirs: a rule naming
# mixture_mvn() is only true when this package is loaded.

#' @noRd
mixture_mvn_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  r("REML", "mixture_mvn", "refused",
    "Refused for the same reason as mixture(): the classes are exchangeable, so the restricted likelihood is not defined.")
  r("quadrature", "mixture_mvn", "refused",
    "Refused in practice: mixture_mvn() clusters a matrix response and carries no scalar random-intercept block, so the scalar-intercept guard rejects it.")
  r("profile", "mixture_mvn", "refused",
    "Refused for the same reason as mixture(): the inner problem would be multimodal.")
  r("cens()", "mixture_mvn", "refused",
    "Refused: mixture_mvn() has no CDF.")
  r("trunc()", "mixture_mvn", "refused",
    "Refused: mixture_mvn() has no CDF.")
  r("mixture_mvn", "simulate", "works",
    "Since v0.36 a structured simulator (fam$sim_ctx) draws a class per row from the gating weights and then a D-variate normal about that class's mean with that class's covariance, assembled from the family-level extras by the same sigma() the likelihood uses. A draw is an n x D matrix. simulate(), posterior_predict() and frm_simulate() all reach the same implementation.")
  r("mixture_mvn", "mvbf", "refused",
    "Refused: the family already takes a matrix response.")
  r("mixture_mvn", "kind:aterm", "untested", "")
  r("group:autocor", "mixture_mvn", "refused",
    "Refused for the same reason as mixture().")
  b$rules()
}
