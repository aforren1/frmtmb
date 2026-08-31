# frmtmb_frame -> negative log-likelihood closure for RTMB::MakeADFun.
#
# Discipline (see SPEC.md section 2): everything observation-length is
# vectorized; loops run only over model structure (blocks, linear
# predictors, responses), which is resolved at tape time; nothing branches
# on parameter values. Data referenced by the closure is baked into the
# tape as constants.
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

  extra_names <- frame$extra_names %||% character(0)

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

    dparv <- list()
    for (lp in lps) {
      if (!is.null(lp$nl_body)) {
        # nonlinear predictor: arbitrary R code over the nlpar values and
        # raw data columns, evaluated straight onto the tape
        ev <- c(dparv[[lp$resp]], lp$data_list)
        eta <- eval(lp$nl_body, ev, lp$nl_env)
        dparv[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
        next
      }
      # as.vector, not drop: it collapses both Matrix and advector results
      eta <- if (ncol(lp$X)) {
        as.vector(lp$X %*% pars[[lp$par]][lp$idx])
      } else {
        rep(0, n)   # threshold-only ordinal model
      }
      if (!is.null(lp$Z)) {
        eta <- eta + as.vector(lp$Z %*% pars[["b"]])
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
        eta <- eta + pars[[lp$par]][lp$idx[mi$col]] * mi$D *
          cz0[mi$codes + 1L]
      }
      dparv[[lp$resp]][[lp$dpar]] <- lp$link$linkinv(eta)
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
        # OBS() drives simulation/OSA machinery, but registers data under
        # the DEPARSED ARGUMENT EXPRESSION: in this loop every response
        # would collide on "y[[r]]" and silently swap data. Univariate
        # fits only; it also has no matrix method.
        yobs <- if (length(resps) == 1L && !is.matrix(y[[r]])) {
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
        if (!is.null(atv[[r]]$cens)) {
          # censoring codes are data, so grouped sub-assignment (one
          # vectorized [<- per group) replaces the density on censored
          # rows without parameter branching or 0 * -Inf hazards
          cen <- atv[[r]]$cens
          Fv <- fam$lcdf(y[[r]], dparv[[r]], atv[[r]])
          i_r <- which(cen == 1)
          i_l <- which(cen == -1)
          i_i <- which(cen == 2)
          if (length(i_r)) ll[i_r] <- log(1 - Fv[i_r])
          if (length(i_l)) ll[i_l] <- log(Fv[i_l])
          if (length(i_i)) {
            F2 <- fam$lcdf(atv[[r]]$cens_y2, dparv[[r]], atv[[r]])
            ll[i_i] <- log(F2[i_i] - Fv[i_i])
          }
        }
        if (!is.null(atv[[r]]$trunc_lb) || !is.null(atv[[r]]$trunc_ub)) {
          disc <- identical(fam$type, "discrete")
          lb <- atv[[r]]$trunc_lb
          if (!is.null(lb) && disc) {
            # inclusive lower bound: P(lb <= Y <= ub) needs F(lb - 1)
            # (brms#1903 off-by-one)
            if (any(lb < 1)) {
              stop("Discrete truncation needs lb >= 1 (lb = 0 is no ",
                   "truncation)", call. = FALSE)
            }
            lb <- lb - 1
          }
          Fub <- if (!is.null(atv[[r]]$trunc_ub)) {
            fam$lcdf(atv[[r]]$trunc_ub, dparv[[r]], atv[[r]])
          } else 1
          Flb <- if (!is.null(lb)) {
            fam$lcdf(lb, dparv[[r]], atv[[r]])
          } else 0
          ll <- ll - log(Fub - Flb)
        }
        nll <- nll - sum(w * ll)
      }
    }
    nll
  }
}
