# frmtmb_frame -> negative log-likelihood closure for RTMB::MakeADFun.
#
# Discipline (see SPEC.md section 2): everything observation-length is
# vectorized; loops run only over model structure (blocks, linear
# predictors, responses), which is resolved at tape time; nothing branches
# on parameter values. Data referenced by the closure is baked into the
# tape as constants.

#' A truncation bound evaluated at the censored rows. An absent bound is
#' the scalar 0 / 1 the untruncated case wants, and a present one is a
#' full-length CDF vector, so the two need different subsetting.
#'
#' @noRd
bound_rows <- function(v, i) if (length(v) == 1L) v else v[i]

#' Turns an assembled `frmtmb_frame` into the negative log-likelihood
#' closure that `RTMB::MakeADFun()` tapes. The returned function takes the
#' parameter list and gives back one value, with the design matrices, the
#' response, and the addition terms baked in as constants. It is a closure
#' and not compiled code because RTMB tapes plain R: the model structure
#' is resolved once here, outside the function the tape records, so the
#' taped code loops over data only through vectorized operations.
#'
#' @noRd
build_objective <- function(frame) {
  lps <- frame$linpreds
  blocks <- frame$re_blocks
  block_fns <- lapply(blocks, function(bk) covstruct_registry[[bk$covstruct]]$nll)
  spec <- frame$spec
  resps <- spec$responses
  rescor <- isTRUE(spec$rescor)
  y <- frame$y
  atv <- frame$aterm_values
  n <- frame$n_obs
  acs <- frame$autocor %||% list()

  extra_names <- frame$extra_names %||% character(0)

  # Cluster-robust scores (R/sandwich.R) need the per-cluster pieces of
  # the objective as a function of one extra parameter each, so that
  # `obj$gr()` at a cluster indicator IS that cluster's score. The hook
  # is a per-response row -> cluster index; it is NULL for every
  # ordinary fit and is only ever set on a private copy of the frame,
  # never on the one stored in the fit. vcov_cluster() refuses every
  # likelihood whose contribution is not the per-row `sum(w * ll)`
  # below (autocor, hmm, mixtures, rescor), so this is the only place a
  # cluster weight has to enter.
  clw_idx <- frame[["cluster_w"]]

  function(pars) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    nll <- 0

    for (i in seq_along(blocks)) {
      bk <- blocks[[i]]
      nll <- nll - block_fns[[i]](pars[["b"]][bk$b_idx],
                                  pars[["theta"]][bk$theta_idx], bk)
    }

    extra <- NULL
    if (length(extra_names)) {
      extra <- lapply(stats::setNames(extra_names, extra_names),
                      function(nm) pars[[nm]])
    }

    # coefficient-space vector for the Z products (rr blocks expand
    # their factors through the loadings)
    bvec <- if (isTRUE(frame$has_rr)) {
      expand_b(frame, pars[["b"]], pars[["theta"]])
    } else {
      pars[["b"]]
    }

    # mi(): observed-or-latent value vectors per imputation response
    mivals <- list()
    for (vn in names(frame$mi_map %||% list())) {
      mm_ <- frame$mi_map[[vn]]
      xv <- y[[vn]]
      xv[mm_$rows] <- pars[["miss"]][mm_$idx]
      mivals[[vn]] <- xv
      if (!is.null(mm_$se)) {
        # measurement model: observed values scatter around the latent
        # truth with known SD (brms me())
        nll <- nll - sum(RTMB::dnorm(y[[vn]][mm_$obs], xv[mm_$obs],
                                     mm_$se[mm_$obs], log = TRUE))
      }
    }

    dparv <- list()
    for (lp in lps) {
      if (!is.null(lp$nl_body)) {
        # nonlinear predictor: arbitrary R code over the nlpar values and
        # raw data columns, evaluated straight onto the tape
        ev <- c(dparv[[lp$resp]], lp$data_list)
        # the body is taped once, so the handler costs nothing per
        # gradient; it exists because a body name that resolved to an
        # environment object instead of a column fails here, far from
        # the cause
        eta <- tryCatch(eval(lp$nl_body, ev, lp$nl_env),
                        error = function(e) nl_body_error(e, lp))
        dparv[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
        dparv[[lp$resp]][[paste0(".eta_", lp$dpar)]] <- eta
        next
      }
      # as.vector, not drop: it collapses both Matrix and advector results
      eta <- if (ncol(lp$X)) {
        as.vector(lp$X %*% pars[[lp$par]][lp$idx])
      } else {
        rep(0, n)   # threshold-only ordinal model
      }
      if (!is.null(lp$Z)) {
        eta <- eta + as.vector(lp$Z %*% bvec)
      }
      if (!is.null(lp$offset)) {
        eta <- eta + lp$offset
      }
      # monotonic terms: scale coefficient (in beta, zero X column)
      # times D * cumulative simplex at the observed category
      for (mi in lp$mo %||% list()) {
        zeta <- exp(c(0, pars[[mi$zeta]]))
        zeta <- zeta / sum(zeta)
        cz0 <- c(0, cumsum(zeta))
        term <- mi$D * cz0[mi$codes + 1L]
        if (!is.null(mi$mult)) term <- term * mi$mult
        eta <- eta + pars[[lp$par]][lp$idx[mi$col]] * term
      }
      # mi(x) terms: coefficient times the observed-or-latent values
      for (mt in lp$mi %||% list()) {
        xv <- mivals[[mt$var]] %||% y[[mt$var]]
        if (!is.null(mt$mult)) xv <- xv * mt$mult
        eta <- eta + pars[[lp$par]][lp$idx[mt$col]] * xv
      }
      # cs(x) terms: n x (K-1) threshold-specific offsets, consumed by
      # the sequential ordinal lpdfs through dpars$.cs
      if (length(lp$cs %||% list())) {
        CS <- 0
        for (ct in lp$cs) {
          bcs <- pars[[ct$par]]
          CS <- CS + RTMB::matrix(ct$vals, n, 1) %*%
            RTMB::matrix(bcs, 1, length(bcs))
        }
        dparv[[lp$resp]][[".cs"]] <- CS
      }
      dparv[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
      # The linear predictor rides along beside the inverse-linked value
      # under a reserved `.eta_` name. An inverse link saturates -
      # plogis(40) is exactly 1, exp(-800) is exactly 0 - so a density
      # that recomputes 1 - mu or log(mu) from the dpar reads -Inf or
      # NaN with a useless gradient in a region the linear predictor
      # itself describes perfectly well. The families that can use the
      # eta scale (see `robust_logit()`, `robust_logmu()`,
      # `gate_logs()`) read it from here, and undoing the link to
      # recover it would defeat the purpose. Only the taped objective
      # supplies it: the numeric post-fit paths (fitted(), simulate(),
      # the CDFs) pass the dpar values alone, and there nothing is
      # differentiated, so the saturation is harmless.
      dparv[[lp$resp]][[paste0(".eta_", lp$dpar)]] <- eta
    }

    if (rescor) {
      # gaussian joint likelihood: standardized residuals against a
      # constant correlation matrix keeps this vectorized even with
      # observation-varying sigmas
      K <- length(resps)
      zvec <- NULL
      logsig <- 0
      for (r in names(resps)) {
        zr <- (y[[r]] - dparv[[r]]$mu) / dparv[[r]]$sigma
        zvec <- if (is.null(zvec)) zr else c(zvec, zr)
        logsig <- logsig + sum(log(dparv[[r]]$sigma))
      }
      Zstd <- RTMB::matrix(zvec, n, K)
      C <- us_chol_cor(pars[["thetar"]], K)
      nll <- nll - sum(RTMB::dmvnorm(Zstd, 0, C, log = TRUE)) + logsig
    } else {
      for (r in names(resps)) {
        fam <- resps[[r]]$family
        w <- atv[[r]]$weights %||% 1
        if (!is.null(clw_idx)) w <- w * pars[["clw"]][clw_idx[[r]]]
        if (!is.null(acs[[r]])) {
          # R-side residual correlation: the response's density is a
          # joint (multivariate normal / t) one per group, not a
          # product over rows, so it replaces fam$lpdf entirely. Every
          # aterm that would reshape a per-row contribution was refused
          # at frame assembly, so nothing below this point applies.
          ac <- acs[[r]]
          R <- autocor_cor(pars[["thetaac"]][ac$theta_idx], ac)
          sg <- dparv[[r]]$sigma
          z <- (y[[r]] - dparv[[r]]$mu) / sg
          lsig <- if (length(sg) == 1L) n * log(sg) else sum(log(sg))
          nu <- if (isTRUE(ac$student)) dparv[[r]]$nu[1] else NULL
          nll <- nll - autocor_loglik(z, R, ac, lsig, nu)
          next
        }
        if (!is.null(fam[["hmm"]]) && !is.null(frame[["hmm_g"]][[r]])) {
          # hidden Markov: the discrete state sequence is summed out
          # EXACTLY by the forward recursion, per sequence, inside the
          # Laplace approximation that integrates any random effects in
          # the state predictors. Nothing here is per row, so none of
          # the aterm machinery below applies - the frame refused every
          # term that would have needed it.
          nll <- nll - hmm_loglik_ad(fam, dparv[[r]],
                                     frame[["hmm_g"]][[r]],
                                     atv[[r]], extra)
          next
        }
        if (!is.null(fam$mix) && !is.null(frame$mix_g[[r]])) {
          # latent-class (group-level) mixture: one class draw per
          # group, so the group's per-observation log-densities sum
          # BEFORE the logsumexp over classes
          mg <- frame$mix_g[[r]]
          lps_pi <- fam$mix$log_pi(dparv[[r]])
          total <- NULL
          for (k in seq_len(fam$mix$K)) {
            ll_k <- fam$mix$comp_lpdf(y[[r]], dparv[[r]], atv[[r]], k)
            g_k <- as.vector(mg$Gt %*% (w * ll_k)) +
              lps_pi[[k]][mg$first]
            total <- if (is.null(total)) g_k else {
              RTMB::logspace_add(total, g_k)
            }
          }
          nll <- nll - sum(total)
          next
        }
        # OBS() drives simulation/OSA machinery, but registers data under
        # the DEPARSED ARGUMENT EXPRESSION: in this loop every response
        # would collide on "y[[r]]" and silently swap data. Univariate
        # fits only; it also has no matrix method.
        yobs <- if (!is.null(mivals[[r]])) {
          # mi() response: the density is evaluated at the
          # observed-or-latent values (the latent entries' Gaussian
          # contribution IS this term)
          mivals[[r]]
        } else if (length(resps) == 1L && !is.matrix(y[[r]])) {
          # a stable symbol gives OSA machinery its observation.name
          .frm_obs <- y[[r]]
          RTMB::OBS(.frm_obs)
        } else {
          y[[r]]
        }
        ll <- if (is.null(fam$extra_pars)) {
          fam$lpdf(yobs, dparv[[r]], atv[[r]])
        } else {
          fam$lpdf(yobs, dparv[[r]], atv[[r]], extra)
        }
        # Truncation bounds are resolved BEFORE the censoring block: a
        # censored row under trunc() observes its event INSIDE the
        # window, so the censored numerators need the same F(lb), F(ub)
        # the normalizer uses. Composing the two the other way round
        # (censor against the whole line, then divide by the window)
        # is a different, wrong likelihood.
        trunc_on <- !is.null(atv[[r]]$trunc_lb) ||
          !is.null(atv[[r]]$trunc_ub)
        Fub <- 1
        Flb <- 0
        if (trunc_on) {
          lb <- atv[[r]]$trunc_lb
          if (!is.null(lb) && identical(fam$type, "discrete")) {
            # inclusive lower bound: P(lb <= Y <= ub) needs F(lb - 1)
            # (brms#1903 off-by-one)
            if (any(lb < 1)) {
              stop("Discrete truncation needs lb >= 1 (lb = 0 is no ",
                   "truncation)", call. = FALSE)
            }
            lb <- lb - 1
          }
          if (!is.null(atv[[r]]$trunc_ub)) {
            Fub <- fam_lcdf(fam, atv[[r]]$trunc_ub, dparv[[r]], atv[[r]],
                            extra)
          }
          if (!is.null(lb)) {
            Flb <- fam_lcdf(fam, lb, dparv[[r]], atv[[r]], extra)
          }
        }
        if (!is.null(atv[[r]]$cens)) {
          # censoring codes are data, so grouped sub-assignment (one
          # vectorized [<- per group) replaces the density on censored
          # rows without parameter branching or 0 * -Inf hazards.
          # Right censoring is P(y < Y <= ub), left censoring is
          # P(lb <= Y <= y); without trunc() those collapse to the
          # familiar 1 - F(y) and F(y), because Fub = 1 and Flb = 0.
          # The discrete F(lb - 1) convention carries over unchanged:
          # F(y) already includes the point y, which is what a
          # left-censored count observes.
          cen <- atv[[r]]$cens
          Fv <- fam_lcdf(fam, y[[r]], dparv[[r]], atv[[r]], extra)
          i_r <- which(cen == 1)
          i_l <- which(cen == -1)
          i_i <- which(cen == 2)
          if (length(i_r)) ll[i_r] <- log(bound_rows(Fub, i_r) - Fv[i_r])
          if (length(i_l)) ll[i_l] <- log(Fv[i_l] - bound_rows(Flb, i_l))
          if (length(i_i)) {
            F2 <- fam_lcdf(fam, atv[[r]]$cens_y2, dparv[[r]],
                             atv[[r]], extra)
            # an interval is already a windowed difference of CDFs, so
            # only the division by the window mass below is missing; an
            # interval reaching outside [lb, ub] is a data contradiction
            # (the response could not have been observed there at all)
            ll[i_i] <- log(F2[i_i] - Fv[i_i])
          }
        }
        if (trunc_on) {
          ll <- ll - log(Fub - Flb)
        }
        nll <- nll - sum(w * ll)
      }
    }
    nll
  }
}
