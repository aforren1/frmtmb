# Latent class analysis (poLCA's measurement model) as a frmtmb family.
#
# The model is the classic one: J polytomous indicator items per
# subject, one latent class per subject, and conditional independence
# of the items within a class,
#
#   P(y_i | x_i) = sum_c w_c(x_i) * prod_j pi_{j,c,y_ij}
#
# with class-membership weights w_c from a multinomial logit on
# covariates (poLCA's "latent class regression", Linzer & Lewis 2011,
# JSS 42(10)) and class-conditional item profiles pi_{j,c,.} as free
# simplexes.
#
# The contained shape is deliberate, the same choice mixture_mvn()
# made: one matrix response, one class variable, no random effects.
# The general form (a mixture over an mvbf with rescor) stays
# unscheduled; see dev/feature-gaps.md.
#
# Two things carry the likelihood. The class-membership weights are
# ORDINARY theta dpars, so the gating linear predictors, their
# priors and their standard errors all come from the existing
# machinery for free - that IS latent class regression. The item
# profiles are family extra parameters (covariate-free by
# construction, as in poLCA), stored as reference-category logits and
# turned into log-probabilities by a logspace_add fold, which is
# stable without ever branching on a parameter.

#' The response as an n x J item matrix.
#'
#' A one-item model reaches the family as a plain vector, because the
#' frame drops an n x 1 matrix response (scale() and friends return
#' one, glmmTMB#937).
#'
#' @noRd
lca_matrix <- function(y) if (is.matrix(y)) y else matrix(y, ncol = 1L)

#' Resolve the per-item category counts.
#'
#' `ncat = NULL` infers them as the largest observed code per item,
#' which is what poLCA does (`apply(y, 2, max)`). A category that is
#' never observed at the top of an item's range is therefore invisible;
#' pass `ncat` to declare it.
#'
#' @noRd
lca_resolve_ncat <- function(y, ncat) {
  J <- ncol(y)
  if (is.null(ncat)) {
    return(as.integer(apply(y, 2L, max, na.rm = TRUE)))
  }
  ncat <- as.integer(ncat)
  if (length(ncat) == 1L) ncat <- rep(ncat, J)
  if (length(ncat) != J) {
    stop("lca(ncat =): ", length(ncat), " category count(s) for ", J,
         " item column(s); give one value per item, or a single value ",
         "for equal-length items", call. = FALSE)
  }
  ncat
}

#' Item codes and observed-item mask, in the shape the log-density
#' gathers with. A missing item response becomes code 1 with mask 0, so
#' the term is switched off by DATA rather than by a parameter branch:
#' the tape has no idea a value was ever missing.
#'
#' @noRd
lca_codes <- function(y) {
  J <- ncol(y)
  miss <- is.na(y)
  codes <- vector("list", J)
  mask <- vector("list", J)
  for (j in seq_len(J)) {
    cj <- as.integer(y[, j])
    mj <- miss[, j]
    cj[mj] <- 1L
    codes[[j]] <- cj
    # [[<- NULL would DELETE the slot and renumber the rest of the list
    if (any(mj)) mask[[j]] <- as.numeric(!mj)
  }
  list(codes = codes, mask = mask)
}

#' The extra-parameter name of item `j`'s logits.
#'
#' One vector per ITEM, rather than one flat vector for the whole
#' model, is what makes the family stateless. The item structure is
#' then readable off the parameters themselves - `J` is how many `pi<j>`
#' entries there are, and item `j` has `length(pi<j>) / K + 1`
#' categories - so nothing has to remember the shape of the data a
#' family object was last used on. The alternative (caching the
#' resolved counts in the family's environment at fit time) silently
#' mis-reports the profiles of the first fit when one saved `lca()`
#' object is reused on a second, differently shaped data set.
#'
#' @noRd
lca_par_name <- function(j) paste0("pi", j)

#' The per-item category counts implied by the extra parameters.
#' Addressed by name, never by position: `mo()` and `cs()` terms add
#' extras of their own (`zeta<n>`, `bcs<n>`) alongside these.
#'
#' @noRd
lca_ncat_from_extra <- function(extra, K) {
  nms <- grep("^pi[0-9]+$", names(extra), value = TRUE)
  J <- length(nms)
  if (!J) {
    stop("lca(): the fit carries no item-profile parameters; this is ",
         "not an lca() fit's parameter list", call. = FALSE)
  }
  vapply(seq_len(J), function(j) {
    as.integer(length(extra[[lca_par_name(j)]]) / K + 1L)
  }, integer(1))
}

#' Log class-conditional category probabilities for one item and one
#' class, from the reference-category logits.
#'
#' The normalizer is a `logspace_add` fold rather than
#' `log(sum(exp()))`: a logit of 40 overflows the naive form while the
#' fold stays exact, and neither version may look at the parameter's
#' value to decide.
#'
#' @noRd
lca_log_probs <- function(a) {
  "c" <- RTMB::ADoverload("c")
  lse <- 0                      # the reference category's logit
  for (m in seq_along(a)) lse <- RTMB::logspace_add(lse, a[m])
  c(0, a) - lse
}

#' Per-item, per-class log-probability tables from the item parameter
#' vectors: `out[[j]][[c]]` is a length-`ncat[j]` (advector) vector.
#' Within an item the layout is class-major, then the non-reference
#' categories 2..C_j.
#'
#' @noRd
lca_profile_tables <- function(extra, K) {
  ncat <- lca_ncat_from_extra(extra, K)
  lapply(seq_along(ncat), function(j) {
    pj <- extra[[lca_par_name(j)]]
    fj <- ncat[j] - 1L
    lapply(seq_len(K), function(k) {
      lca_log_probs(pj[(k - 1L) * fj + seq_len(fj)])
    })
  })
}

#' Deterministic starting profiles.
#'
#' A categorical likelihood has no analog of `mixture()`'s
#' response-quantile means, and poLCA's answer is random restarts. The
#' analog used here is a deterministic one: score each subject by the
#' mean of its item codes rescaled to `[0, 1]`, cut the scores at `K`
#' equal quantiles, and take each slice's smoothed empirical category
#' proportions as that class's starting profile. The labels are
#' therefore ordered by that score - class 1 is the low-score end - so
#' a given data set always starts, and usually ends, with the same
#' labeling. Laplace smoothing keeps a slice that never saw a category
#' off the -Inf boundary.
#'
#' @noRd
lca_init_extras <- function(y, ncat, K) {
  n <- nrow(y)
  J <- ncol(y)
  span <- pmax(ncat - 1L, 1L)
  sc <- rowMeans(sweep(y - 1, 2L, span, "/"), na.rm = TRUE)
  fin <- is.finite(sc)
  sc[!fin] <- if (any(fin)) stats::median(sc[fin]) else 0
  # rank-based cut: ties in the score cannot collapse a slice to zero
  # rows the way quantile() breakpoints can on coarse binary items
  slice <- pmin(1L + ((rank(sc, ties.method = "first") - 1L) * K) %/% n, K)
  out <- vector("list", J)
  for (j in seq_len(J)) {
    v <- numeric(0)
    for (k in seq_len(K)) {
      cj <- y[slice == k, j]
      cj <- cj[!is.na(cj)]
      cnt <- tabulate(cj, nbins = ncat[j]) + 1   # Laplace smoothing
      p <- cnt / sum(cnt)
      v <- c(v, log(p[-1L]) - log(p[1L]))
    }
    out[[j]] <- v
  }
  stats::setNames(out, vapply(seq_len(J), lca_par_name, ""))
}

#' Log mixing weights from the theta dpars: multinomial logit against
#' the last class, exactly `mixture()`'s gating.
#'
#' @noRd
lca_log_pi <- function(dpars_all, K) {
  Ts <- lapply(seq_len(K - 1L), function(k) dpars_all[[paste0("theta", k)]])
  Ts[[K]] <- 0 * dpars_all[["theta1"]]
  lse <- Ts[[1L]]
  if (K > 1L) for (k in seq.int(2L, K)) lse <- RTMB::logspace_add(lse, Ts[[k]])
  lapply(Ts, function(t_) t_ - lse)
}

#' Per-subject log-density of the items GIVEN class k, summed over
#' items. This is the mixture component density, so `mixture_probs()`
#' reaches it through `fam$mix$comp_lpdf` and gets posterior class
#' probabilities with no LCA-specific code.
#'
#' @noRd
lca_comp_lpdf <- function(y, K, extra, k) {
  y <- lca_matrix(y)
  tabs <- lca_profile_tables(extra, K)
  if (length(tabs) != ncol(y)) {
    # a saved lca() family object reused across two differently shaped
    # data sets used to reach here with the first fit's item structure
    stop("lca(): the fit carries ", length(tabs), " item profile(s) ",
         "but the response has ", ncol(y), " item column(s)",
         call. = FALSE)
  }
  cc <- lca_codes(y)
  S <- 0
  for (j in seq_along(tabs)) {
    term <- tabs[[j]][[k]][cc$codes[[j]]]
    if (!is.null(cc$mask[[j]])) term <- term * cc$mask[[j]]
    S <- S + term
  }
  S
}

#' Latent class analysis
#'
#' `lca(K)` fits the classic latent class measurement model, the one
#' 'poLCA' fits (Linzer and Lewis 2011): `J` polytomous indicator items
#' per subject, one latent class `c` in `1..K` per subject, and the
#' items conditionally independent given the class,
#'
#' \deqn{P(y_i) = \sum_c w_c \prod_j \pi_{j,c,y_{ij}}.}
#'
#' The response is a MATRIX, one row per subject and one column per
#' item, holding integer category codes `1..C_j`. Write it as
#' `cbind(item1, item2, ...)` on the left of the formula, or attach a
#' matrix column to the data (`dd$Y <- data.matrix(dd[items])`) and
#' name it. Items may have different numbers of categories.
#'
#' The class-membership weights are the `theta1 ... theta{K-1}` dpars,
#' multinomial logit against class `K`, and the main model formula
#' applies to every one of them. So `bf(cbind(a, b, c) ~ 1)` is the
#' plain measurement model and `bf(cbind(a, b, c) ~ age + educ)` is
#' poLCA's latent class regression: covariates on class membership come
#' from the ordinary linear-predictor machinery, with the usual
#' `fixef()`, `confint()` and `hypothesis()` on top. Individual gating
#' predictors are overridable as `bf(Y ~ x, theta2 ~ 1)` (all but
#' `theta1`, which the main formula owns).
#'
#' The item profiles are family extra parameters, held as
#' reference-category logits, one vector `pi<j>` per item (item `j`'s
#' `K * (C_j - 1)` free logits, class-major), and reported
#' on the probability scale by [lca_profiles()]. They take no linear
#' predictor - a covariate acts on class membership, never on an item's
#' conditional response probability, which is what makes the classes
#' interpretable as measurement.
#'
#' @section Labeling and starting values:
#' A latent class likelihood is invariant to relabeling the classes and
#' is genuinely multimodal, so the answer depends on where the
#' optimizer starts. poLCA uses random restarts (`nrep`). The starting
#' values here are deterministic instead: subjects are scored by the
#' mean of their item codes rescaled to `[0, 1]`, cut into `K`
#' equal-count slices by that score, and each slice's smoothed
#' empirical category proportions become one class's starting profile.
#' Class 1 is therefore the low-score end and class `K` the high-score
#' end, and re-running the same data gives the same labeling.
#'
#' That fixes reproducibility, not multimodality. Do what
#' `poLCA(..., nrep = 10)` does and compare starts before reading a
#' solution: perturb the item parameters and refit,
#'
#' ```
#' p0 <- fit$frame$par_template[fit$frame$extra_names]
#' refits <- replicate(10, simplify = FALSE,
#'   frm(bf(Y ~ 1), family = lca(K = 3), data = dd,
#'       start = lapply(p0, function(v) v + rnorm(length(v)))))
#' sapply(refits, logLik)
#' ```
#'
#' and keep the best. [frm_allfit()] is a different check: it re-runs
#' the four optimizers from the SAME start, so it tests the optimizer,
#' not the surface. Ordering the gating intercepts through `lower` and
#' `upper` is the way to pin the labeling itself.
#'
#' Asking for more classes than the data hold drives item
#' probabilities to 0 and 1. The optimizer then reports singular
#' convergence and the standard errors come back `NaN`, which is the
#' boundary showing through rather than a fault; poLCA lands on the
#' same solutions. Compare `AIC()` and `BIC()` across `K` and read
#' [lca_probs()]'s entropy before committing to a `K`.
#'
#' @section Missing item responses:
#' With `na.rm = TRUE` (the default, poLCA's default) a subject with
#' any missing item is dropped by the usual `na.action`, and the
#' message naming the dropped rows is the ordinary one. With
#' `na.rm = FALSE` the subject is kept and only the missing item's
#' factor leaves that subject's likelihood, which is poLCA's
#' `na.rm = FALSE` behavior; the mask is data, so the tape never sees a
#' branch. A subject missing EVERY item contributes a constant and is
#' better dropped.
#'
#' @section What is and is not supported:
#' Post-fit, [lca_probs()] gives posterior class-membership
#' probabilities per subject (with the relative-entropy classification
#' diagnostic attached), [lca_profiles()] gives the item profile table,
#' and `simulate()` draws a class per subject and then its items.
#' `fitted()`, `predict(type = "response")` and `residuals()` are
#' refused: the response is a matrix of nominal codes, so an "expected
#' item code" would be an average of arbitrary labels. Read
#' [lca_probs()] and [lca_profiles()] instead. `predict()` itself
#' returns the gating linear predictor (`theta1` by default, any
#' `theta` with `dpar =`), on `newdata` as well.
#'
#' The gating coefficients are ordinary fixed effects, so `fixef()`,
#' `confint()` (Wald, profile and uniroot), `hypothesis()`,
#' `set_prior()`, `lower`/`upper` bounds and
#' `frmtmb.sample::frm_sample()` all work on
#' them; `anova()` compares nested gating formulas at one `K`.
#'
#' Refused in this version: random effects and smooths anywhere in the
#' model (latent classes plus continuous random effects is the
#' growth-mixture shape - use `mixture(..., groups = ~g)` for that),
#' `REML` and `frmtmb_control(profile = TRUE)` (the classes are
#' exchangeable, so there is no single inner mode), `quadrature`, every
#' addition term (`weights()`, `cens()`, `trunc()`, `se()`, `mi()`,
#' `trials()`), [mvbf()], and `residuals(type = "osa")`.
#'
#' @param K Number of latent classes (at least 2).
#' @param ncat Number of categories per item: a vector of length `J`, a
#'   single value for equally-sized items, or `NULL` (the default) to
#'   infer each item's count as its largest observed code, which is
#'   what poLCA does.
#' @param na.rm If `TRUE` (default) subjects with any missing item are
#'   dropped by `na.action`; if `FALSE` they are kept and each missing
#'   item's term is masked out of that subject's likelihood.
#' @return A `frmtmb_family`.
#' @references Linzer, D. A. and Lewis, J. B. (2011). poLCA: An R
#'   Package for Polytomous Variable Latent Class Analysis. *Journal of
#'   Statistical Software*, 42(10), 1-29.
#' @seealso [lca_probs()], [lca_profiles()], [mixture()] for a mixture
#'   of continuous responses, [mixture_mvn()] for model-based
#'   clustering of a numeric matrix.
#' @examples
#' set.seed(1)
#' # four binary items measuring two well-separated classes
#' n <- 300
#' cl <- rbinom(n, 1, 0.4) + 1
#' pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
#' Y <- matrix(0L, n, 4)
#' for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
#' dd <- data.frame(x = rnorm(n))
#' dd$Y <- Y
#'
#' fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
#' lca_profiles(fit)
#' head(lca_probs(fit))
#'
#' # latent class regression: covariates gate class membership
#' frm(bf(Y ~ x), family = lca(K = 2), data = dd)
#' @export
lca <- function(K, ncat = NULL, na.rm = TRUE) {
  if (missing(K) || length(K) != 1L || !is.finite(K) || K != round(K) ||
        K < 2) {
    stop("lca() needs a single whole number of latent classes, at ",
         "least 2, e.g. lca(K = 3)", call. = FALSE)
  }
  K <- as.integer(K)
  if (!is.null(ncat)) {
    if (!is.numeric(ncat) || !length(ncat) || any(!is.finite(ncat)) ||
          any(ncat != round(ncat)) || any(ncat < 2)) {
      stop("lca(ncat =) must be whole numbers of at least 2, one per ",
           "item", call. = FALSE)
    }
  }
  # na.rm decides whether a missing item is masked out of that
  # subject's likelihood or the subject is dropped, so isTRUE() reading
  # a mistake as FALSE changes the estimand
  check_flag(na.rm, "na.rm")
  na_rm <- isTRUE(na.rm)
  # Nothing below writes to this closure. A family object is a value
  # that may be built once and handed to several fits, so the item
  # structure is carried by the PARAMETERS (one pi<j> vector per item)
  # and re-read from them wherever it is needed, rather than cached
  # here at fit time.
  dpn <- paste0("theta", seq_len(K - 1L))
  fam <- frmtmb_family(
    paste0("lca(K = ", K, ")"),
    dpars = dpn,
    links = stats::setNames(rep(list("identity"), K - 1L), dpn),
    lpdf = function(y, dpars, aterms, extra) {
      lw <- lca_log_pi(dpars, K)
      ll <- NULL
      for (k in seq_len(K)) {
        llk <- lca_comp_lpdf(y, K, extra, k) + lw[[k]]
        ll <- if (is.null(ll)) llk else RTMB::logspace_add(ll, llk)
      }
      ll
    },
    valid_y = function(y, aterms) {
      if (is.list(y) || is.character(y) ||
            (is.matrix(y) && !ncol(y))) {
        stop("lca(K = ", K, "): the response must be a matrix of item ",
             "codes, one row per subject and one column per item; write ",
             "cbind(item1, item2, ...) ~ ... or attach a matrix column ",
             "with data.matrix()", call. = FALSE)
      }
      y <- lca_matrix(y)
      bad <- !is.na(y) & (y != round(y) | y < 1)
      if (any(bad)) {
        stop("lca(): item responses must be whole-number category codes ",
             "1..C_j (a factor column is converted with as.integer() or ",
             "data.matrix(), which uses its level order)", call. = FALSE)
      }
      if (!na_rm && any(rowSums(!is.na(y)) == 0L)) {
        stop("lca(na.rm = FALSE): ", sum(rowSums(!is.na(y)) == 0L),
             " subject(s) have no observed item at all and carry no ",
             "information; remove them", call. = FALSE)
      }
      nc <- lca_resolve_ncat(y, ncat)
      obs <- vapply(seq_len(ncol(y)), function(j) {
        length(unique(y[!is.na(y[, j]), j]))
      }, integer(1))
      if (any(obs < 2L)) {
        bad_j <- which(obs < 2L)
        stop("lca(): item column(s) ",
             paste(bad_j, collapse = ", "),
             " take fewer than two distinct values, so no class can be ",
             "told apart by them; drop the item(s)", call. = FALSE)
      }
      over <- vapply(seq_len(ncol(y)), function(j) {
        max(y[, j], na.rm = TRUE) > nc[j]
      }, TRUE)
      if (any(over)) {
        stop("lca(ncat =): item column(s) ",
             paste(which(over), collapse = ", "),
             " hold codes above the declared category count",
             call. = FALSE)
      }
      for (at in c("weights", "cens", "trunc_lb", "trunc_ub", "se",
                   "trials")) {
        if (!is.null(aterms[[at]])) {
          stop("lca() does not support addition terms on the response; ",
               "the item matrix carries no per-row weight, window or ",
               "known standard error to attach one to", call. = FALSE)
        }
      }
      invisible(NULL)
    },
    init_dpars = list(),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) {
        stop("An lca() fit has no fitted mean: the response is a matrix ",
             "of nominal item codes, so averaging them would average ",
             "arbitrary labels. Use lca_probs() for posterior class ",
             "membership and lca_profiles() for the item profiles",
             call. = FALSE)
      }
    ),
    extra_pars = function(y, aterms) {
      y <- lca_matrix(y)
      lca_init_extras(y, lca_resolve_ncat(y, ncat), K)
    },
    sim = function(dpars, aterms, n, extra) {
      lw <- lca_log_pi(dpars, K)
      P <- vapply(lw, function(l) exp(rep(l, length.out = n)), numeric(n))
      ks <- vapply(seq_len(n), function(i) {
        sample.int(K, 1L, prob = P[i, ])
      }, integer(1))
      # the fit's own item parameters say how many items there are and
      # how many categories each has
      tabs <- lca_profile_tables(extra, K)
      out <- matrix(0L, n, length(tabs))
      for (j in seq_along(tabs)) {
        for (k in seq_len(K)) {
          idx <- which(ks == k)
          if (!length(idx)) next
          pk <- exp(tabs[[j]][[k]])
          out[idx, j] <- sample.int(length(pk), length(idx),
                                    replace = TRUE, prob = pk)
        }
      }
      out
    },
    primary_dpars = dpn
  )
  # mixture()'s component interface, deliberately: one implementation of
  # the posterior class probabilities then serves mixture(),
  # mixture_mvn() and lca() alike, which is why mixture_probs() works on
  # an lca fit and lca_probs() is documented as the same matrix.
  fam$mix <- list(
    K = K,
    comp_lpdf = function(y, dpars, aterms, k, extra) {
      lca_comp_lpdf(y, K, extra, k)
    },
    comp_dpars = function(dpars_all, k) list(),
    log_pi = function(dpars_all) lca_log_pi(dpars_all, K)
  )
  fam$lca <- list(K = K, na_rm = na_rm)
  fam[["structure"]] <- lca_structure(na_rm)
  fam
}

#' Everything the core needs to know about an `lca()` response, as the
#' one object the structured-family protocol reads
#' (`frmtmb_structure()`).
#'
#' An LCA likelihood DOES factorize: one subject is one row, and the
#' family's own `lpdf` computes it. So there is no `loglik` here and the
#' structure is a capability declaration, which is the whole reason a
#' rowwise family would carry one: three of the refusals below used to
#' be branches the core took by name (`fam$na_response` at frame
#' assembly, the `has_mixture()` gate in fit.R, an `is_lca_family()`
#' test in predict.R), and none of them was ever a property of the
#' method that raised it.
#'
#' @noRd
lca_structure <- function(na_rm) {
  frmtmb_structure(
    # with na.rm = FALSE a missing item is masked out of that subject's
    # likelihood rather than costing the subject its row, so the NA is
    # data the family reads
    keep_na = !na_rm,
    check_frame = function(spec, frame) {
      check_lca_structure(spec, frame[["linpreds"]])
    },
    # an LCA implements mixture()'s component interface, so the shared
    # posterior serves it with no LCA-specific code
    latent_probs = function(fit, block) mixture_posterior(fit),
    supports = structure_supports_all(reml = FALSE, profile = FALSE,
                                      osa = FALSE),
    refusals = c(
      mixture_multimodal_refusals("an lca() family"),
      list(osa = paste0(
        "residuals(type = \"osa\") is not available for an lca() ",
        "fit: one observation is a subject's whole item response ",
        "pattern, not a single value with a conditional CDF to step ",
        "through. Use lca_probs() to see how sharply each subject ",
        "is classified"))
    )
  )
}

#' Whether a response carries an `lca()` family.
#'
#' @noRd
is_lca_family <- function(fam) !is.null(fam[["lca"]])

#' Structural refusals an `lca()` fit cannot state from `valid_y()`,
#' because they live in the linear predictors rather than the response.
#' This is the structure's `check_frame`, so it runs once the assembled
#' predictors exist.
#'
#' @noRd
check_lca_structure <- function(spec, linpreds) {
  lca_resps <- vapply(spec$responses, function(r) is_lca_family(r$family),
                      TRUE)
  if (!any(lca_resps)) return(invisible(NULL))
  # mvbf() is already refused upstream, by the extra-parameter guard in
  # the response loop, which runs before the predictors exist
  for (lp in linpreds) {
    if (!is.null(lp$Z) || length(lp$smooth %||% list())) {
      stop("lca() does not support random effects, smooths or gp() ",
           "terms in the class-membership predictor ('", lp$dpar,
           "'). A latent class with continuous random effects is the ",
           "growth-mixture model; mixture(..., groups = ~g) fits that ",
           "shape", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' The fitted item profiles of an `lca()` fit
#'
#' The class-conditional item-response probabilities `pi[j, c, k]`,
#' which are the main output of a latent class analysis: for each item,
#' a `K x C_j` table of the probability of each category given each
#' class. This is poLCA's `$probs`.
#'
#' @param fit A `frmtmb_fit` with an [lca()] family.
#' @return A named list with one `K x C_j` probability matrix per item,
#'   of class `frmtmb_lca_profiles`, carrying the estimated class sizes
#'   (the mean prior class probabilities, poLCA's `$P`) in the
#'   `class_sizes` attribute.
#' @seealso [lca()], [lca_probs()]
#' @examples
#' set.seed(2)
#' n <- 200
#' cl <- rbinom(n, 1, 0.5) + 1
#' pr <- rbind(c(0.9, 0.85, 0.8), c(0.1, 0.15, 0.2))
#' Y <- matrix(0L, n, 3)
#' for (j in 1:3) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
#' dd <- data.frame(row = seq_len(n))
#' dd$Y <- Y
#' fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
#'
#' pf <- lca_profiles(fit)
#' pf
#' pf[[1]]                       # item 1's K x C table
#' attr(pf, "class_sizes")
#' @export
lca_profiles <- function(fit) {
  if (!inherits(fit, "frmtmb_fit")) {
    stop("lca_profiles() takes a fitted model from frm()", call. = FALSE)
  }
  rspec <- single_response(fit, "lca_profiles()")
  fam <- rspec$family
  if (!is_lca_family(fam)) {
    stop("lca_profiles() needs a fit with an lca() family", call. = FALSE)
  }
  K <- fam$lca$K
  # the item structure comes from THIS fit's parameters, so one saved
  # lca() object handed to several differently shaped data sets reports
  # each fit correctly
  ex <- fit$estimates[fit$frame$extra_names %||% character(0)]
  tabs <- lca_profile_tables(ex, K)
  yv <- lca_matrix(fit$frame$y[[rspec$resp_name]])
  item_names <- colnames(yv) %||% paste0("item", seq_along(tabs))
  out <- lapply(seq_along(tabs), function(j) {
    m <- t(vapply(tabs[[j]], exp, numeric(length(tabs[[j]][[1L]]))))
    dimnames(m) <- list(paste0("class", seq_len(K)),
                        paste0("cat", seq_len(ncol(m))))
    m
  })
  names(out) <- item_names
  dp <- eval_dpars(fit)[[rspec$resp_name]]
  lw <- fam$mix$log_pi(dp)
  sizes <- vapply(lw, function(l) mean(exp(rep(l, length.out = nrow(yv)))),
                  numeric(1))
  structure(out, class_sizes = stats::setNames(
    sizes, paste0("class", seq_len(K))),
    class = "frmtmb_lca_profiles")
}

#' @export
print.frmtmb_lca_profiles <- function(x, digits = 4, ...) {
  cs <- attr(x, "class_sizes")
  cat("<lca profiles> ", length(cs), " classes, ", length(x),
      " items\n", sep = "")
  cat("\nEstimated class sizes (mean prior probability):\n")
  print(round(cs, digits))
  y <- unclass(x)
  attr(y, "class_sizes") <- NULL
  for (nm in names(y)) {
    cat("\n", nm, ":\n", sep = "")
    print(round(y[[nm]], digits))
  }
  invisible(x)
}

#' Posterior class membership of an `lca()` fit
#'
#' One row per subject, one column per latent class, rows summing to
#' one: the probability that a subject belongs to each class given its
#' observed item responses and its gating covariates. This is poLCA's
#' `$posterior`.
#'
#' The relative entropy of the classification is attached as the
#' `entropy` attribute: `1 - sum(-p log p) / (n log K)`, which is 1
#' for a partition with no ambiguity and 0 when every subject is
#' equally likely to be in any class. Values above about 0.8 are the
#' usual rule of thumb for classes worth naming.
#'
#' This is [mixture_probs()] under an LCA-specific name and check; the
#' two return the same matrix for an [lca()] fit.
#'
#' @param fit A `frmtmb_fit` with an [lca()] family.
#' @return A `n x K` matrix of posterior class probabilities with an
#'   `entropy` attribute.
#' @seealso [lca()], [lca_profiles()], [mixture_probs()]
#' @examples
#' set.seed(3)
#' n <- 200
#' cl <- rbinom(n, 1, 0.5) + 1
#' pr <- rbind(c(0.9, 0.85, 0.8, 0.9), c(0.1, 0.15, 0.2, 0.1))
#' Y <- matrix(0L, n, 4)
#' for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
#' dd <- data.frame(row = seq_len(n))
#' dd$Y <- Y
#' fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
#'
#' p <- lca_probs(fit)
#' head(p)
#' attr(p, "entropy")            # classification quality, 0 to 1
#' table(max.col(p), truth = cl) # the modal assignment, up to relabeling
#' @export
lca_probs <- function(fit) {
  if (!inherits(fit, "frmtmb_fit")) {
    stop("lca_probs() takes a fitted model from frm()", call. = FALSE)
  }
  rspec <- single_response(fit, "lca_probs()")
  if (!is_lca_family(rspec$family)) {
    stop("lca_probs() needs a fit with an lca() family; mixture_probs() ",
         "covers mixture() and mixture_mvn()", call. = FALSE)
  }
  P <- latent_probs(fit)
  lp <- ifelse(P > 0, log(P), 0)
  ent <- 1 + sum(P * lp) / (nrow(P) * log(ncol(P)))
  attr(P, "entropy") <- ent
  P
}

## ---- the compatibility matrix's lca() rows ---------------------------

# These rows live beside the family rather than in R/compat.R, for the
# same reason hmm()'s do: the family and everything the package claims
# about it travel together.

#' @noRd
lca_compat_rules <- function() {
  b <- compat_rule_builder()
  r <- b$r
  # the three fitting modes lca() refuses. They sat in compat.R's REML,
  # quadrature and profile sections while one gate in fit.R stated all
  # three mixture-type refusals in one message; lca() states its own now.
  r("REML", "lca", "refused",
    "Refused for the same reason as mixture(): the latent classes are exchangeable, so the likelihood is multimodal in the fixed effects REML integrates out and the restricted likelihood is not defined.")
  r("quadrature", "lca", "refused",
    "Refused in practice: lca() refuses random effects outright, so there is no block for the scalar-intercept guard to accept and it rejects the fit.")
  r("profile", "lca", "refused",
    "Refused for the same reason as mixture(): profiling moves the gating coefficients into an inner Laplace problem that has one mode per class labeling.")
  r("lca", "kind:aterm", "refused",
    "Refused: an lca() response is a matrix of item codes with no per-row weight, censoring window, trial count or known standard error to attach an addition term to. One message covers the whole set.")
  r("lca", "kind:covstruct", "refused",
    "Refused: v1 has no random effects anywhere in the model. Latent classes plus continuous random effects is the growth-mixture shape, which mixture(..., groups = ~g) already fits.")
  r("lca", "kind:special", "untested",
    "s(), t2() and gp() are refused with the random-effect message, because each builds a random-effect block; mo() is verified. mi() and cs() in the gating predictor are untested.")
  r("lca", "mo()", "works",
    "Verified by a tiny fit: a monotonic predictor enters the gating linear predictor like any other term.")
  r("lca", "s()", "refused",
    "Refused: a smooth is a random-effect block, and lca() refuses random effects.")
  r("lca", "t2()", "refused",
    "Refused: a tensor smooth is a random-effect block.")
  r("lca", "gp_pred()", "refused",
    "Refused: gp() builds a random-effect block.")
  r("lca", "mvbf", "refused",
    "Refused by the extra-parameter guard: the item profiles are family extra parameters, which multivariate fits do not carry. A second response would need its own class variable, which is the general mixture-over-mvbf model.")
  r("lca", "rescor", "refused",
    "Refused: rescor = TRUE requires all responses to be gaussian, and the error names lca() among the ones that are not.")
  r("lca", "mixture", "refused",
    "Refused: the two families cannot both own the response. lca() IS a mixture, over conditionally independent categorical items.")
  r("lca", "mixture_mvn", "refused",
    "Refused: one response carries one family.")
  r("lca", "nl", "refused",
    "Refused: nl = TRUE requires a family with a single mu location parameter, and lca()'s location parameters are the gating thetas.")
  r("lca", "fitted", "refused",
    "Refused: the response is a matrix of nominal item codes, so there is no mean to fit. lca_probs() and lca_profiles() are the post-fit surface.")
  r("lca", "residuals", "refused",
    "Refused for the same reason as fitted(): no fitted mean, so no residual.")
  r("lca", "residuals_osa", "refused",
    "Refused: one observation is a subject's whole item response pattern, not a value with a univariate conditional CDF to step through.")
  r("lca", "predict", "conditional",
    "predict() returns the gating linear predictor (theta1 by default, any theta with dpar =), including on newdata. type = \"response\" is refused with the fitted() message.")
  r("lca", "simulate", "works",
    "Verified: simulate() draws a class per subject from its gating weights and then its items from that class's profile, and returns an n x J matrix per draw. Refitting a 4000-subject draw recovers the profiles it came from.")
  r("lca", "confint_profile", "works",
    "Verified: profile intervals on a gating coefficient run.")
  r("lca", "emmeans", "untested", "")
  r("lca", "prior", "works",
    "Verified by a tiny fit: a prior on the gating coefficients is an ordinary penalty on the outer problem.")
  r("lca", "bounds", "works",
    "Verified by a tiny fit. Bounds on the gating coefficients are also the one way to order the classes and pin the labeling.")
  r("lca", "sparse_x", "works",
    "Verified by a tiny fit: the gating design is an ordinary fixed-effect design.")
  r("lca", "kind:mode", "untested",
    "See the REML, profile, quadrature, prior, bounds and sparse_x rules; autoscale and verbose are untested.")
  b$rules()
}

frmtmb_register_compat(features = c(lca = "structure"),
                       rules = lca_compat_rules)
