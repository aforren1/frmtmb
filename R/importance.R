# Importance-sampling correction of the Laplace approximation, after
# Skaug and Fournier (2006, Computational Statistics & Data Analysis
# 51:699-709) and ADMB's `-is` option.
#
# WHY this and not more quadrature. `quadrature = TRUE` hands the
# integral to TMBad's marginal_gk transform, whose C++ carries a scalar
# mu and sigma per integrand: a block of dimension two returns NaN, and
# the `dim` in the spec is never read. It also freezes its rescaling
# once, so a stale calibration BREAKS the objective, which is what
# quad_fit() spends five seeds and three recalibration rounds fighting.
# A frozen Gaussian importance proposal cannot break that way. It stays
# a valid proposal wherever theta walks and only loses efficiency, and
# the effective sample size measures how much. Any block dimension, and
# the remaining error is a number the fit reports rather than an
# assumption it makes.
#
# THE ESTIMATOR. For group g with random-effect block u_g of dimension
# q, the marginal likelihood is the integral of exp(-f_g(theta, u_g))
# over u_g, where f_g is the group's joint negative log-density: its
# rows' likelihood plus its block's prior. Take the Laplace Gaussian as
# the proposal, centered at the conditional mode u_hat_g with
# covariance H_g^-1, and L_g L_g' = H_g^-1. With standard normal draws
# z_{g,i} and u_{g,i} = u_hat_g + L_g z_{g,i},
#
#   log m_g(theta) = log mean_i exp(-f_g(theta, u_{g,i}) + |z_{g,i}|^2 / 2)
#                    + (q / 2) log(2 pi) + log|L_g|
#
# and the objective is -sum_g log m_g(theta), which is an unbiased
# estimate of the exact marginal log-likelihood at EVERY theta, not
# only near the anchor. The draws are numeric constants baked into the
# closure, so the tape's parameters are the fixed effects, the
# dispersion and the covariance parameters alone: `random = NULL`, no
# inner problem, and sdreport works the ordinary way.
#
# STACKED, NEVER LOOPED. For draw i every group's block is set to its
# i-th draw at once. `Z U` is one sparse product forming the whole
# n x N random-effect contribution as a CONSTANT, the family's density
# is called ONCE on a length n*N vector with the response and the
# addition terms recycled, and the per-group aggregation is one sparse
# indicator product. A loop over draws would put N copies of the
# likelihood on the tape as N separate calls.

#' The covariance structures whose block density factorizes over
#' grouping levels AND is Gaussian in the level's coefficients. Both
#' properties are load-bearing: factorization is what makes the
#' integral a product of per-group integrals, and Gaussianity is what
#' lets `imp_prior_terms()` recover the per-level density from the
#' registry's summed one. Everything else is refused by name in
#' `check_importance_scope()`, and `imp_plan()` measures the
#' factorization again on the fitted Hessian rather than trusting this
#' list.
#'
#' @noRd
imp_covstructs <- c("us", "diag", "homdiag", "cs", "homcs", "toep",
                    "homtoep", "ar1", "hetar1", "ou", "exp", "gau", "mat")

#' Structures that carry a grouping factor but correlate its LEVELS,
#' so the integral does not factorize over groups at all. Named
#' separately from the plain unsupported set because the remedy
#' differs: there is no draw count that fixes these.
#'
#' @noRd
imp_crosslevel <- c("gr_cov", "gr_prec", "equalto", "car", "spde")

#' Every refusal `importance =` owes the user, raised before any tape
#' is built so the message names the conflict rather than a downstream
#' symptom. The style follows the `quadrature = TRUE` guards in
#' R/fit.R: say what cannot be combined, and say what to do instead.
#'
#' @noRd
check_importance_scope <- function(spec, frame, template, REML, quadrature,
                                   control) {
  if (isTRUE(quadrature)) {
    stop("`importance` and `quadrature = TRUE` are two different ",
         "corrections of the same approximation and cannot both be ",
         "asked for. Adaptive Gauss-Kronrod marginalizes one scalar ",
         "random intercept at a time; the importance correction ",
         "reweights the Laplace Gaussian and takes any block ",
         "dimension. Keep one", call. = FALSE)
  }
  if (REML) {
    stop("`importance` cannot be combined with REML = TRUE. The ",
         "correction reweights the integral over the random effects ",
         "only, and REML puts the fixed effects into that same ",
         "integral, where the reweighting has no proposal for them. ",
         "Use REML = FALSE", call. = FALSE)
  }
  if (isTRUE(control$profile)) {
    stop("`importance` cannot be combined with ",
         "frmtmb_control(profile = TRUE). Profiling moves the fixed ",
         "effects into the inner Laplace problem, and the corrected ",
         "objective has no inner problem left to move them into. Drop ",
         "profile = TRUE", call. = FALSE)
  }
  if (!is.null(template[["miss"]])) {
    stop("`importance` cannot be combined with mi(). The imputed ",
         "values are latent variables of their own, with no grouping ",
         "factor to give them a per-group proposal. Fit with ",
         "importance = 0", call. = FALSE)
  }
  if (length(frame[["autocor"]] %||% list())) {
    stop("`importance` cannot be combined with the residual ",
         "correlation term ", frame[["autocor"]][[1L]]$label,
         ": the correction resamples a random effect against a ",
         "PRODUCT of per-row densities, and this residual is one joint ",
         "density over each group, so no per-row integrand exists. Use ",
         "importance = 0, or REML = TRUE", call. = FALSE)
  }
  if (isTRUE(spec$rescor)) {
    stop("`importance` cannot be combined with rescor: the responses ",
         "share one multivariate normal density per row, which is not ",
         "the rowwise density the correction resamples. Fit the ",
         "responses without rescor, or use importance = 0",
         call. = FALSE)
  }
  if (length(spec$responses) != 1L) {
    stop("`importance` supports one response, and this model has ",
         length(spec$responses), " (",
         paste(names(spec$responses), collapse = ", "), "). A ",
         "multivariate model spreads its groups over several ",
         "likelihood terms, which the first version does not gather. ",
         "Fit the responses one at a time, or use importance = 0",
         call. = FALSE)
  }
  for (rn in names(spec$responses)) {
    if (!is.null(fam_structure(spec$responses[[rn]]$family)[["loglik"]])) {
      stop("`importance` cannot correct the '",
           spec$responses[[rn]]$family[["family"]], "' family: it ",
           "supplies its own log-likelihood, which does not factorize ",
           "over rows, so a group's rows have no separable integrand ",
           "to resample. This is the same restriction quadrature has. ",
           "Use importance = 0", call. = FALSE)
    }
  }
  for (lp in frame[["linpreds"]]) {
    if (!is.null(lp[["nl_body"]])) {
      stop("`importance` cannot correct a nonlinear predictor (",
           lp[["resp"]], ".", lp[["dpar"]], "). A nonlinear body mixes ",
           "parameter values with raw data columns, and the corrected ",
           "objective evaluates the predictor once per draw, where a ",
           "column would recycle silently against a longer vector. Use ",
           "importance = 0", call. = FALSE)
    }
    if (length(lp[["cs"]] %||% list())) {
      stop("`importance` cannot correct a cs() term (", lp[["resp"]],
           ".", lp[["dpar"]], "): its threshold-specific offsets are a ",
           "matrix per observation, which the first version does not ",
           "stack over draws. Use importance = 0", call. = FALSE)
    }
  }
  if (!is.null(frame[["map"]][["b"]])) {
    stop("`importance` cannot correct a model whose random-effect ",
         "coefficients are mapped: the proposal is built per grouping ",
         "level from the conditional Hessian, and a map merges or ",
         "fixes coefficients across levels. Use importance = 0",
         call. = FALSE)
  }
  blocks <- frame[["re_blocks"]] %||% list()
  if (!length(blocks)) {
    stop("`importance` needs a random-effect block to correct, and this ",
         "model has none. The correction reweights an integral over ",
         "random effects, and a model without any has no such integral ",
         "and no approximation to improve on. Use importance = 0",
         call. = FALSE)
  }
  labs <- vapply(blocks, function(bk) bk[["term_label"]], "")
  for (bk in blocks) {
    if (is.null(bk[["levels"]])) {
      stop("`importance` needs a grouping factor, and `",
           bk[["term_label"]], "` has none: a smooth, a Gaussian ",
           "process or an HSGP basis is one field over all ",
           "observations, so there are no independent groups to give ",
           "separate proposals to. Use importance = 0", call. = FALSE)
    }
  }
  # SEVERAL BLOCKS, ONE FACTOR. Distributional regression writes them
  # by construction: `(1 | g)` in mu and `(1 | g)` in sigma are two
  # blocks, never one, because only an |ID| key merges terms. They
  # still factorize over the levels of g, so a level's coefficients
  # from every block are drawn together from one joint proposal and
  # the integral stays a product over groups. Blocks over DIFFERENT
  # factors do not: the integral is nested, and no per-group proposal
  # exists for it.
  grps <- vapply(blocks, function(bk) bk[["group_name"]] %||% NA_character_, "")
  if (length(unique(grps)) != 1L) {
    stop("`importance` takes several random-effect blocks only when ",
         "they share ONE grouping factor, and this model spreads ",
         length(blocks), " blocks over ", length(unique(grps)),
         " factors (",
         paste0("`", labs, "` over ", grps, collapse = ", "),
         "). Crossed or nested factors make the marginal likelihood ",
         "one integral over every factor at once, which does not ",
         "split into the per-group integrals the correction resamples. ",
         "Fit with importance = 0", call. = FALSE)
  }
  lv <- blocks[[1L]][["levels"]]
  for (i in seq_along(blocks)[-1L]) {
    if (identical(blocks[[i]][["levels"]], lv)) next
    odd <- setdiff(union(lv, blocks[[i]][["levels"]]),
                   intersect(lv, blocks[[i]][["levels"]]))
    stop("`importance` needs every block over `", grps[[1L]],
         "` to carry the same grouping levels, and `", labs[[1L]],
         "` has ", length(lv), " where `", labs[[i]], "` has ",
         length(blocks[[i]][["levels"]]),
         if (length(odd)) paste0(" (first difference: '", odd[[1L]], "')"),
         ". The proposal draws one level's coefficients from every ",
         "block at once, so a level carried by only some of them has ",
         "no joint Gaussian to be drawn from. Fit with importance = 0",
         call. = FALSE)
  }
  for (bk in blocks) {
    if (bk[["covstruct"]] %in% imp_crosslevel) {
      stop("`importance` cannot correct the '", bk[["covstruct"]],
           "' structure in `", bk[["term_label"]], "`: it correlates the ",
           "grouping LEVELS with each other through a supplied ",
           "relationship or neighbor matrix, so the marginal likelihood ",
           "is one integral over every level at once and does not split ",
           "into per-group ones. Use importance = 0", call. = FALSE)
    }
    if (!bk[["covstruct"]] %in% imp_covstructs) {
      stop("`importance` cannot correct the '", bk[["covstruct"]],
           "' structure in `", bk[["term_label"]], "`. The correction ",
           "recovers each level's prior density from the block density, ",
           "which needs that density to be Gaussian in the level's ",
           "coefficients and independent between levels; the supported ",
           "structures are ",
           paste(imp_covstructs, collapse = ", "),
           ". Use importance = 0", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' WHERE A GROUP'S COEFFICIENTS LIVE. The one object the whole
#' several-blocks generalization turns on.
#'
#' `b` is level-major WITHIN a block, so a level's `dim` coefficients
#' are consecutive there; but the blocks are CONCATENATED, so with more
#' than one block over the factor a group owns one such run per block
#' and the runs sit `n_levels * dim` apart. A group's coefficients are
#' then scattered, and the contiguous spelling `(k - 1) * q +
#' seq_len(q)` is wrong in every place it appeared.
#'
#' `idx` is that index: `q` rows by `n_group` columns, column k holding
#' the positions group k owns, blocks in order. Every consumer reads a
#' column of it (the Hessian slice, the draw placement, the
#' verification swap) or a row of it (the half-norms, the
#' per-coefficient draw slices), so the sites cannot drift apart the
#' way five separate spellings could. `pos_group` is its transpose as a
#' lookup, one group label per position, which is what the Hessian
#' block-diagonality gate needs.
#'
#' Positions are `c_idx`, the Z column space. It coincides with the `b`
#' parameter space for every structure in `imp_covstructs`; `rr`, the
#' one structure where the two differ, is refused by the whitelist, and
#' `imp_plan()`'s count check would catch it anyway because its two
#' spaces have different lengths.
#'
#' @noRd
imp_layout <- function(blocks) {
  q <- vapply(blocks, function(bk) as.integer(bk[["dim"]]), 1L)
  ng <- blocks[[1L]][["n_levels"]]
  offset <- cumsum(c(0L, q))[seq_along(q)]
  qt <- sum(q)
  idx <- matrix(0L, qt, ng)
  for (m in seq_along(blocks)) {
    # level-major within the block is exactly what makes a block's
    # contribution one reshape: column k of the reshape is level k's
    # run, and it lands on the block's own rows of the group's vector
    idx[offset[[m]] + seq_len(q[[m]]), ] <-
      matrix(blocks[[m]][["c_idx"]], q[[m]], ng)
  }
  pos_group <- integer(qt * ng)
  pos_group[as.vector(idx)] <- rep(seq_len(ng), each = qt)
  list(blocks = blocks, q = q, offset = offset, qt = qt, n_group = ng,
       nb = qt * ng, idx = idx, pos_group = pos_group,
       levels = blocks[[1L]][["levels"]],
       group_name = blocks[[1L]][["group_name"]],
       label = paste(vapply(blocks, function(bk) bk[["term_label"]], ""),
                     collapse = " + "))
}

#' Which grouping level each observation's random-effect contribution
#' belongs to, read off the sparsity of `Z` rather than off the model
#' frame: `Z` is what the objective actually multiplies, so a row that
#' reaches two levels is visible here and nowhere else.
#'
#' Gathered across EVERY block, not one, so the same pass that catches
#' a multi-membership term also catches two blocks that disagree about
#' which level a row belongs to. That second reading is what makes
#' "one factor with one level set" a measured property of `Z` rather
#' than an inference from two matching `group_name` strings.
#'
#' A row reaching NO level (an all-zero `Z` row) is assigned to the
#' first group on purpose. Its density does not depend on `u` at all,
#' so it enters that group's `f_g` as a constant, and a constant
#' factors straight back out of a log-mean-exp: the total is unchanged
#' and the group's importance weights are shifted by the same amount,
#' so its effective sample size is unchanged too.
#'
#' @noRd
imp_group_map <- function(frame, lay) {
  n <- frame[["n_obs"]]
  ng <- lay[["n_group"]]
  rows <- list()
  levs <- list()
  labs <- list()
  for (bk in lay[["blocks"]]) {
    ci <- bk[["c_idx"]]
    # level-major: a level's `dim` coefficients are consecutive in c_idx
    col_level <- rep(seq_len(bk[["n_levels"]]), each = bk[["dim"]])
    for (lp in frame[["linpreds"]]) {
      if (is.null(lp[["Z"]])) next
      tz <- methods::as(lp[["Z"]][, ci, drop = FALSE], "TsparseMatrix")
      if (!length(tz@i)) next
      rows[[length(rows) + 1L]] <- tz@i + 1L
      levs[[length(levs) + 1L]] <- col_level[tz@j + 1L]
      labs[[length(labs) + 1L]] <- rep(bk[["term_label"]], length(tz@i))
    }
  }
  rr <- unlist(rows, use.names = FALSE) %||% integer(0)
  ll <- unlist(levs, use.names = FALSE) %||% integer(0)
  tl <- unlist(labs, use.names = FALSE) %||% character(0)
  # one entry per (row, level) pair; a row that keeps two of them
  # touches two groups and the integral does not factorize. The key is
  # formed in DOUBLE arithmetic: n * (ng + 1) overflows an integer on a
  # design of a few million rows, and an NA key would silently drop the
  # very clash this is looking for.
  keep <- !duplicated(as.numeric(rr) * (ng + 1) + ll)
  rr <- rr[keep]
  ll <- ll[keep]
  tl <- tl[keep]
  clash <- unique(rr[duplicated(rr)])
  if (length(clash)) {
    stop("`importance` needs every observation to belong to one ",
         "grouping level, and ", length(clash), " of them reach ",
         "several through `",
         paste(unique(tl[rr %in% clash]), collapse = "`, `"),
         "` (first at row ", min(clash), "). A multi-membership term ",
         "shares a row between groups, so the marginal likelihood does ",
         "not split into per-group integrals. Use importance = 0",
         call. = FALSE)
  }
  row_level <- integer(n)
  row_level[rr] <- ll
  row_level[row_level == 0L] <- 1L
  list(row_level = row_level,
       S = Matrix::sparseMatrix(i = row_level, j = seq_len(n),
                                x = rep(1, n),
                                dims = c(ng, n)))
}

#' Antithetic standard normal draws: `n_draw / 2` columns and their
#' negatives. Measured on the package's own probe design, the pairing
#' cuts the variance of the corrected objective by two to three times
#' at every draw count tried, which is worth more than the draws it
#' costs, so it is not optional.
#'
#' The stream is private. The draws must neither depend on nor disturb
#' the session's random state, because "same seed, same answer" is a
#' documented property of the fit and a user who sets a seed for a
#' simulation must get the same simulation whether or not a fit ran.
#'
#' @noRd
imp_draws <- function(nb, n_draw, seed) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())))
  }
  set.seed(seed)
  z <- matrix(stats::rnorm(nb * (n_draw %/% 2L)), nb, n_draw %/% 2L)
  cbind(z, -z)
}

#' The number of draws actually taken. Antithetic pairing needs an even
#' count, so an odd request is rounded UP: rounding down could silently
#' turn `importance = 1` into no draws at all.
#'
#' @noRd
imp_n_draw <- function(n) 2L * ((as.integer(n) + 1L) %/% 2L)

#' The proposal, frozen at one outer parameter vector: conditional
#' modes, per-level Hessian blocks, their Cholesky factors, and the
#' draws they map the standard normals through.
#'
#' The modes and the Hessian come from the Laplace objective's own
#' inner Newton solve, so the proposal IS the Laplace Gaussian and
#' nothing has to be re-derived. `spHess(random = TRUE)` is the
#' conditional Hessian of the joint negative log-density; under the
#' scope above it is block diagonal by grouping level, and that is
#' checked here rather than assumed, because it is the one property
#' that makes the whole per-group construction valid.
#'
#' @noRd
imp_plan <- function(lap_obj, frame, lay, par, n_draw, seed) {
  lap_obj$fn(par)                       # inner solve leaves the modes
  full <- lap_obj$env$last.par
  uhat <- as.numeric(full[lap_obj$env$random])
  hess <- lap_obj$env$spHess(full, random = TRUE)
  # q is the group's TOTAL dimension, summed over the blocks: one
  # proposal per group, over every coefficient that group owns. The
  # blocks are independent given theta, but the LIKELIHOOD couples them
  # within a group (a group's rows read its mu block and its sigma
  # block at once), so the joint Hessian is the right proposal and its
  # cross-block entries are load-bearing, not noise.
  q <- lay[["qt"]]
  ng <- lay[["n_group"]]
  nb <- lay[["nb"]]
  idx <- lay[["idx"]]
  if (length(uhat) != nb) {
    stop("`importance` expected ", nb, " conditional modes for `",
         lay[["label"]], "` and the Laplace fit reports ",
         length(uhat), ". The proposal is built per grouping level, so ",
         "the two counts have to agree", call. = FALSE)
  }
  th <- methods::as(hess, "TsparseMatrix")
  lev <- lay[["pos_group"]]
  cross <- lev[th@i + 1L] != lev[th@j + 1L]
  if (any(cross)) {
    worst <- max(abs(th@x[cross]))
    scale <- max(abs(Matrix::diag(hess)))
    if (worst > 1e-8 * scale) {
      stop("`importance` needs the conditional Hessian of `",
           lay[["label"]], "` to be block diagonal by grouping ",
           "level, and it is not: the largest entry linking two levels ",
           "is ", format(worst, digits = 3), " against a diagonal of ",
           format(scale, digits = 3), ". The groups are then not ",
           "independent and the marginal likelihood does not split ",
           "into per-group integrals. Use importance = 0",
           call. = FALSE)
    }
  }
  z <- imp_draws(nb, n_draw, seed)
  nd <- ncol(z)
  u <- matrix(0, nb, nd)
  logdet_l <- numeric(ng)
  if (q == 1L) {
    # one scalar block and nothing else: H is a diagonal, so L is one
    # reciprocal square root per level and no factorization is needed.
    # Still addressed through `idx`, so the fast path reads the same
    # index as the general one rather than re-deriving the layout.
    p1 <- idx[1L, ]
    hd <- Matrix::diag(hess)[p1]
    if (any(!is.finite(hd)) || any(hd <= 0)) {
      stop("`importance` found a non-positive conditional variance for ",
           sum(!is.finite(hd) | hd <= 0), " level(s) of `",
           lay[["label"]], "`, so no Gaussian proposal exists ",
           "there. The Laplace fit is at a singular variance ",
           "component; refit it first", call. = FALSE)
    }
    sdv <- 1 / sqrt(hd)
    u[p1, ] <- uhat[p1] + sdv * z[p1, , drop = FALSE]
    logdet_l <- log(sdv)
  } else {
    for (k in seq_len(ng)) {
      ii <- idx[, k]
      hk <- as.matrix(hess[ii, ii, drop = FALSE])
      pk <- tryCatch(chol(hk), error = function(e) NULL)
      if (is.null(pk)) {
        stop("`importance` could not factor the conditional Hessian of ",
             "level '", lay[["levels"]][k], "' in `", lay[["label"]],
             "`: it is not positive definite, so that group has no ",
             "Gaussian proposal. The Laplace fit sits on a singular ",
             "variance component; refit it first", call. = FALSE)
      }
      # L = P^-1 with P'P = H gives L L' = H^-1, upper triangular and
      # backsolved rather than inverted
      lk <- backsolve(pk, diag(q))
      u[ii, ] <- uhat[ii] + lk %*% z[ii, , drop = FALSE]
      logdet_l[k] <- -sum(log(diag(pk)))
    }
  }
  # 0.5 |z_{g,i}|^2, one value per group and draw. Row j of `idx` is
  # the j-th coefficient of every group at once, which is what makes
  # this a q-term loop over G x N matrices and not a loop over groups.
  zsq <- matrix(0, ng, nd)
  for (j in seq_len(q)) {
    zsq <- zsq + z[idx[j, ], , drop = FALSE]^2
  }
  zsq <- 0.5 * zsq
  list(u = u, zsq = zsq, logdet_l = logdet_l, n_draw = nd, par = par,
       uhat = uhat, q = q, n_group = ng, levels = lay[["levels"]],
       idx = idx)
}

#' Repeat an observation-length quantity once per draw. Anything that
#' is not observation-length is a scalar the density recycles by
#' itself, so it passes through untouched.
#'
#' @noRd
imp_expand <- function(v, ridx, n) {
  if (is.null(v)) return(NULL)
  if (is.list(v)) return(lapply(v, imp_expand, ridx = ridx, n = n))
  if (is.matrix(v)) {
    if (nrow(v) != n) return(v)
    return(v[ridx, , drop = FALSE])
  }
  if (length(v) == n) return(v[ridx])
  v
}

#' One level's log prior density, per group and draw, recovered FROM
#' the registry density rather than reimplemented beside it.
#'
#' `covstruct_registry[[k]]$nll` returns the block's log density summed
#' over every level, and the per-group integrand needs one level's
#' worth. Writing each structure's density a second time here would be
#' two copies of the same mathematics, free to drift. Instead, evaluate
#' the registry's own density on a ONE-LEVEL block at the canonical
#' patterns: it is Gaussian in the level's coefficients (the scope
#' guarantees it), so the value at zero is the normalizer and the
#' differences at the unit vectors and their pairs are the entries of
#' the precision matrix. q + q(q+1)/2 + 1 small evaluations, taped
#' once, and the gaussian identity test pins the result exactly.
#'
#' With several blocks over the factor the group's log prior is their
#' SUM: the blocks are independent given theta, so the group's
#' precision is block diagonal and no cross terms appear, and each
#' block reads only its own slice of the group's coefficients. Each
#' term is extracted from its own block's registry density at its own
#' theta.
#'
#' @noRd
imp_prior_terms <- function(lay, u_slices) {
  per <- lapply(seq_along(lay[["blocks"]]), function(m) {
    imp_block_prior(lay[["blocks"]][[m]],
                    u_slices[lay[["offset"]][[m]] +
                               seq_len(lay[["q"]][[m]])])
  })
  tix <- lapply(lay[["blocks"]], function(bk) bk[["theta_idx"]])
  function(theta) {
    out <- per[[1L]](theta[tix[[1L]]])
    for (m in seq_along(per)[-1L]) {
      out <- out + per[[m]](theta[tix[[m]]])
    }
    out
  }
}

#' One block's contribution to that sum, on the coefficient slices it
#' owns. Factored out of `imp_prior_terms()` so the extraction is
#' written once whatever the block count.
#'
#' It needs no `ADoverload()` of its own: the arithmetic here is `*`
#' and `+` on advectors, and `covstruct_registry[[k]]$nll` installs the
#' overloads it needs itself (R/covstruct.R). That is worth stating
#' rather than assuming, because an overload lost in a factoring-out
#' fails far from its cause.
#'
#' @noRd
imp_block_prior <- function(bk, u_slices) {
  bfn <- covstruct_registry[[bk[["covstruct"]]]]$nll
  q <- bk[["dim"]]
  bk1 <- bk
  bk1[["n_levels"]] <- 1L
  zero <- numeric(q)
  unit <- lapply(seq_len(q), function(j) replace(zero, j, 1))
  pair <- list()
  if (q > 1L) {
    for (j in seq_len(q - 1L)) {
      for (k in (j + 1L):q) {
        pair[[paste0(j, ".", k)]] <- replace(zero, c(j, k), 1)
      }
    }
  }
  function(theta) {
    c0 <- bfn(zero, theta, bk1)
    # dg[[j]] = -P[j, j] / 2
    dg <- lapply(unit, function(e) bfn(e, theta, bk1) - c0)
    quad <- dg[[1L]] * (u_slices[[1L]] * u_slices[[1L]])
    if (q > 1L) {
      for (j in 2:q) {
        quad <- quad + dg[[j]] * (u_slices[[j]] * u_slices[[j]])
      }
      for (j in seq_len(q - 1L)) {
        for (k in (j + 1L):q) {
          # -P[j, k], from the value at e_j + e_k less both diagonals
          cjk <- bfn(pair[[paste0(j, ".", k)]], theta, bk1) - c0 -
            dg[[j]] - dg[[k]]
          quad <- quad + cjk * (u_slices[[j]] * u_slices[[k]])
        }
      }
    }
    c0 + quad
  }
}

#' The theta-only negative log-likelihood closure: `-sum_g log m_g`,
#' with the draws, their random-effect contribution to every linear
#' predictor, and the per-group indicator all baked in as constants.
#'
#' Returns the closure and, beside it, an `amat` that hands back the
#' G x N matrix of log importance weights at a numeric parameter list.
#' The diagnostics read that matrix, so the effective sample sizes and
#' the Monte Carlo standard error describe the objective that was
#' optimized and not a second computation of it.
#'
#' @noRd
build_importance_objective <- function(frame, lay, gmap, plan) {
  lps <- frame[["linpreds"]]
  spec <- frame[["spec"]]
  rn <- names(spec$responses)[[1L]]
  fam <- spec$responses[[rn]]$family
  n <- frame[["n_obs"]]
  nd <- plan[["n_draw"]]
  ng <- plan[["n_group"]]
  q <- plan[["q"]]
  ridx <- rep.int(seq_len(n), nd)
  yraw <- imp_expand(frame[["y"]][[rn]], ridx, n)
  av <- imp_expand(frame[["aterm_values"]][[rn]], ridx, n)
  wts <- av[["weights"]] %||% 1
  extra_names <- frame[["extra_names"]] %||% character(0)
  smat <- gmap[["S"]]

  # The draws are constants, so Z U is a constant too: the whole
  # n x N random-effect contribution to each linear predictor is one
  # sparse product, formed here, once, off the tape.
  zu <- lapply(lps, function(lp) {
    if (is.null(lp[["Z"]])) return(NULL)
    as.vector(as.matrix(lp[["Z"]] %*% plan[["u"]]))
  })
  # one G x N slice per within-GROUP coefficient, gathered by row j of
  # the scattered index rather than by a stride, because a group's
  # coefficients are only contiguous when there is one block. The
  # quadratic form multiplies pairs of them, and the products are
  # formed inside the closure so only q slices are held on the R side
  u_slices <- lapply(seq_len(q), function(j) {
    plan[["u"]][plan[["idx"]][j, ], , drop = FALSE]
  })
  prior_fn <- imp_prior_terms(lay, u_slices)
  zsq <- plan[["zsq"]]
  # A per-group constant subtracted before the exponential, so that the
  # weights sit near 1 at the anchor. It cancels exactly from the
  # result, and being a CONSTANT it puts no branch on parameter values,
  # which a running maximum would.
  offs <- matrix(0, ng, nd)
  const <- sum(q / 2 * log(2 * pi) + plan[["logdet_l"]])

  amat <- function(pars) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    extra <- NULL
    if (length(extra_names)) {
      extra <- lapply(stats::setNames(extra_names, extra_names),
                      function(nm) pars[[nm]])
    }
    dpv <- list()
    for (k in seq_along(lps)) {
      lp <- lps[[k]]
      # the fixed part is the same for every draw, so it is built once
      # at length n and expanded; only Z u varies from draw to draw
      eta <- lp_eta_fixed(lp, pars, n, list(), frame[["y"]])[ridx]
      if (!is.null(zu[[k]])) eta <- eta + zu[[k]]
      dpv[[lp[["dpar"]]]] <- lp[["link"]]$linkinv(eta)
      dpv[[paste0(".eta_", lp[["dpar"]])]] <- eta
    }
    ll <- wts * row_lpdf(fam, yraw, yraw, dpv, av, extra)
    agg <- smat %*% RTMB::matrix(ll, n, nd)
    # The sparse product gives an advector matrix on the tape and a
    # Matrix object in the numeric path, and a Matrix loses its
    # dimensions through exp(). The test is on CLASS, so it resolves
    # while the tape is being built and leaves no branch on it.
    if (!inherits(agg, "advector")) agg <- as.matrix(agg)
    agg + prior_fn(pars[["theta"]]) + zsq
  }

  list(amat = amat,
       fn = function(pars) {
         a <- amat(pars)
         # RTMB::rowSums, not the base one: base's drops the advector
         # class and the tape then ends in a bare complex vector
         -sum(offs[, 1L] + log(RTMB::rowSums(exp(a - offs)) / nd)) - const
       },
       set_offset = function(value) {
         offs <<- matrix(value, ng, nd)
         invisible(NULL)
       })
}

#' Per-group effective sample size and the Monte Carlo standard error
#' of the corrected objective, both from the same log-weight matrix the
#' objective sums.
#'
#' `ESS_g = (sum w)^2 / sum w^2` on the group's normalized weights. The
#' standard error follows from the delta method: `log m_g` has variance
#' `1/ESS_g - 1/N`, and the groups are independent, so the total is the
#' square root of their sum.
#'
#' @noRd
imp_ess <- function(a, n_draw) {
  a <- as.matrix(a)
  w <- exp(a - apply(a, 1L, max))
  ess <- rowSums(w)^2 / rowSums(w^2)
  list(ess = ess / n_draw,
       mcse = sqrt(sum(pmax(1 / ess - 1 / n_draw, 0))))
}

#' The pin that keeps the stacked reimplementation honest.
#'
#' The closure above rebuilds the linear predictors and the group sums
#' itself, so it could in principle disagree with the objective it is
#' correcting for some family, link or addition term. It cannot
#' disagree silently, and the check is PER GROUP: a pair of errors that
#' cancelled in the total would pass a check on the total, and the
#' per-group values are exactly what the effective sample sizes and the
#' log-mean-exp rest on.
#'
#' Per group without needing a reference value: take two draw columns
#' and compare the DIFFERENCE of the group's joint log-density between
#' them. Setting one group to its column-j draw while every other group
#' stays at its column-k draw leaves every other group's term unchanged,
#' so it cancels from the difference and what is left is that one
#' group's. The total across groups is checked as well, because a
#' constant added to every group cancels from the differences but not
#' from the sum.
#'
#' Costs one evaluation of the plain objective per group checked, of a
#' cheap R closure and not a tape. Every group is checked up to 64 of
#' them, and an evenly spaced 64 beyond that, which bounds the cost on
#' a model with thousands of levels while still catching the kind of
#' error that is a code fault rather than a data accident.
#'
#' What it does NOT catch, stated so nobody trusts it further than it
#' goes: an offset that is identical in every draw AND sums to zero
#' across groups cancels from the differences and from the total alike.
#' The plain objective only ever returns sums, so no per-group ABSOLUTE
#' reference can be recovered from it, and that residue is the price.
#' It is not a shape any plausible bug takes: an error that varies with
#' the draw, or that fails to cancel across groups, is caught.
#'
#' @noRd
imp_verify <- function(io, nll, plan, template, par) {
  pars <- imp_par_list(template, par)
  a <- as.matrix(io$amat(pars))
  nd <- plan[["n_draw"]]
  ng <- plan[["n_group"]]
  refuse <- function(what, ours, theirs) {
    stop("`importance` could not reproduce this model's joint ",
         "log-density from its per-group pieces (", what, ": ",
         format(ours, digits = 10), " against ",
         format(theirs, digits = 10), "). The correction integrates a ",
         "group at a time, so the two have to agree exactly. Report ",
         "this model, and fit it with importance = 0 meanwhile",
         call. = FALSE)
  }
  agrees <- function(ours, theirs) {
    is.finite(ours) && is.finite(theirs) &&
      abs(ours - theirs) <= 1e-6 * max(1, abs(theirs))
  }
  # the total, at each of two columns: catches an offset shared by
  # every group, which the per-group differences below cannot see
  for (j in unique(c(1L, nd))) {
    joint <- pars
    joint[["b"]] <- plan[["u"]][, j]
    ours <- sum(a[, j] - plan[["zsq"]][, j])
    theirs <- -nll(joint)
    if (!agrees(ours, theirs)) refuse(paste("draw", j), ours, theirs)
  }
  if (nd < 2L) return(invisible(NULL))
  # and each group on its own
  base <- pars
  base[["b"]] <- plan[["u"]][, nd]
  nll_base <- nll(base)
  mine <- (a[, 1L] - plan[["zsq"]][, 1L]) -
    (a[, nd] - plan[["zsq"]][, nd])
  gs <- if (ng <= 64L) {
    seq_len(ng)
  } else {
    unique(round(seq(1, ng, length.out = 64L)))
  }
  for (g in gs) {
    swapped <- base
    # the group's positions in `b`, scattered across the blocks: swap
    # the WHOLE group and nothing else, or the difference below stops
    # isolating one group's term
    idx <- plan[["idx"]][, g]
    swapped[["b"]][idx] <- plan[["u"]][idx, 1L]
    theirs <- nll_base - nll(swapped)
    if (!agrees(mine[g], theirs)) {
      refuse(paste0("group '", plan[["levels"]][g], "'"), mine[g], theirs)
    }
  }
  invisible(NULL)
}

#' The outer parameter template the corrected tape is built on: every
#' component of the Laplace template except the random effects, which
#' are now constants inside the closure.
#'
#' @noRd
imp_template <- function(template, random) {
  template[setdiff(names(template), random)]
}

#' A full parameter list at an outer vector, for the numeric paths
#' (verification and diagnostics) that call the closure directly.
#'
#' @noRd
imp_par_list <- function(template, par) {
  out <- template
  pn <- names(par)
  for (cp in unique(pn)) {
    if (is.null(out[[cp]])) next
    out[[cp]][] <- unname(par[pn == cp])
  }
  out
}

#' Fit by iterated importance correction: freeze the proposal at the
#' Laplace optimum, optimize theta, refreeze at the new theta, optimize
#' again.
#'
#' Refreezing is what keeps the weights from degenerating. The estimate
#' is unbiased at every theta whatever the proposal, so a round that
#' moves far is not wrong, only noisy; refreezing buys the accuracy
#' back.
#'
#' The loop stops at a fixed point: a theta at which the objective
#' anchored AT theta is stationary. The freshly anchored gradient says
#' that directly, and the parameter move is the cruder backstop for a
#' round that simply went nowhere; either ends the loop, and the round
#' cap bounds the cost. What the fit reports comes from the LAST
#' freeze, which is anchored at the estimate.
#'
#' @noRd
importance_fit <- function(nll, template, random, map, lap_obj, control,
                           bounds, par_units, frame, n_draw, vb) {
  if (vb) t0 <- vb_now()
  lap_opt <- optimize_obj(lap_obj, control, bounds, par_units)
  if (vb) vb_stage("importance warm start", t0, vb_opt_detail(lap_opt))
  lay <- imp_layout(frame[["re_blocks"]])
  gmap <- imp_group_map(frame, lay)
  otpl <- imp_template(template, random)
  omap <- map[intersect(names(map), names(otpl))]
  if (!length(omap)) omap <- NULL
  seed <- control$importance_seed %||% 1L
  rounds <- max(1L, as.integer(control$importance_rounds %||% 5L))

  # One freeze: build the proposal at `at`, centre the weights on it,
  # and tape. The offset is the anchor's own row maxima, so the
  # exponentials start at 1 and stay in range as the optimizer walks.
  # It cancels from the value exactly, and being a constant it puts no
  # branch on parameter values the way a running maximum would.
  freeze <- function(at, verify) {
    if (vb) tf <- vb_now()
    plan <- imp_plan(lap_obj, frame, lay, at, n_draw, seed)
    io <- build_importance_objective(frame, lay, gmap, plan)
    amat <- as.matrix(io$amat(imp_par_list(template, at)))
    io$set_offset(apply(amat, 1L, max))
    if (verify) imp_verify(io, nll, plan, template, at)
    obj <- RTMB::MakeADFun(io$fn, otpl, map = omap, silent = TRUE)
    if (vb) {
      vb_stage("importance tape", tf,
               paste0(plan[["n_draw"]], " draws over ", plan[["n_group"]],
                      " groups"))
    }
    list(plan = plan, io = io, obj = obj, amat = amat)
  }

  par <- lap_opt$par
  # The first proposal is the Laplace Gaussian at the Laplace optimum,
  # which is where it is closest to the conditional it approximates.
  fz <- freeze(par, TRUE)
  # The corrected objective AT the Laplace estimates. Every later round
  # is judged against it, because the whole point of the correction is
  # to do better than the fit it corrects.
  start_ess <- imp_ess(fz$amat, fz$plan[["n_draw"]])
  start_value <- fz$obj$fn(par)
  opt <- NULL
  moves <- numeric(0)
  grads <- numeric(0)
  for (k in seq_len(rounds)) {
    if (vb) t0 <- vb_now()
    # start_par, not the tape's own par: MakeADFun sets that from the
    # START TEMPLATE, so without it every round would re-optimize from
    # a cold start and the Laplace warm start would be used for the
    # anchor and thrown away for the optimizer. On a ridge that also
    # lets consecutive rounds land in different basins, which is the
    # very failure the divergence guard exists to catch.
    opt <- optimize_obj(fz$obj, control, bounds, par_units, verbose = vb,
                        start_par = par)
    moved <- max(abs((opt$par - par) * (par_units %||% 1)))
    moves <- c(moves, moved)
    par <- opt$par
    # Every round ends by refreezing AT its own estimate. That serves
    # twice: it is the next round's proposal, and if the loop stops
    # here it is the tape the fit reports from. Reporting from the
    # anchor is not a detail. The optimizer minimizes the exact
    # objective plus a Monte Carlo error that is EXACTLY zero at the
    # anchor and grows away from it (measured: zero variance over seeds
    # at the anchor for a gaussian response, at every draw count), so
    # the value it stops at is optimistic by construction and by an
    # amount no draw count removes. At the anchor the weights are as
    # near equal as this proposal can make them, which is also what
    # makes the effective sample sizes below describe the estimate that
    # is actually reported.
    fz <- freeze(par, FALSE)
    # The fixed point this loop is looking for is a theta at which the
    # objective anchored AT theta is stationary, so the gradient of the
    # fresh tape is the criterion that says so directly; the
    # parameter move is the cruder one that catches a round which
    # simply did not go anywhere. Either ends the loop.
    grad <- try(max(abs(fz$obj$gr(par) * (par_units %||% 1))),
                silent = TRUE)
    if (inherits(grad, "try-error")) grad <- NA_real_
    grads <- c(grads, grad)
    if (vb) {
      # the round's own effective sample sizes, at its own anchor: a
      # proposal that has stopped covering the integrand shows up here,
      # round by round, rather than only in the closing line
      re <- imp_ess(fz$amat, fz$plan[["n_draw"]])
      vb_stage("importance round", t0,
               paste0(vb_opt_detail(opt), ", moved ",
                      format(moved, digits = 3), ", re-anchored grad ",
                      format(grad, digits = 3), ", ESS min ",
                      format(min(re$ess), digits = 3), " median ",
                      format(stats::median(re$ess), digits = 3)))
    }
    if (moved < imp_move_tol) break
    if (is.finite(grad) && grad < control$grad_tol) break
  }

  opt$objective <- fz$obj$fn(par)
  ess <- imp_ess(fz$amat, fz$plan[["n_draw"]])
  grad <- grads[length(grads)]
  # The iteration must not have made the fit WORSE than the correction
  # it started from. Both values estimate the same exact marginal
  # likelihood, each at its own theta and each from a proposal frozen
  # there, so they compare directly. A rise means the fixed-point
  # iteration diverged, which is what an unidentified covariance does
  # to it: the surface has a ridge, the optimizer wanders along it, and
  # the rounds chase each other out to a collapsed variance component
  # where every weight is trivially equal and the effective sample
  # sizes look perfect. Refusing beats reporting a corrected fit that
  # is worse than the Laplace fit it corrected.
  worse <- opt$objective - start_value
  if (is.finite(worse) && worse > imp_worse_tol(ess$mcse, start_ess$mcse)) {
    stop("`importance` did not converge on this model: the corrected ",
         "negative log-likelihood ROSE from ",
         format(start_value, digits = 8), " at the Laplace estimates ",
         "to ", format(opt$objective, digits = 8), " after ",
         length(moves), " rounds, so the iteration moved away from the ",
         "answer instead of toward it. That is what an unidentified ",
         "covariance does to it: check whether the Laplace fit itself ",
         "converged and whether `", lay[["label"]], "` has enough ",
         "rows per group to identify every variance component",
         call. = FALSE)
  }
  if (vb) {
    vb_stage("importance report", vb_now(),
             paste0(vb_opt_detail(opt), ", min ESS ",
                    format(min(ess$ess), digits = 3), ", MCSE ",
                    format(ess$mcse, digits = 3)))
  }
  # the layout goes out with the fit because the reporting sites name a
  # level and a grouping factor, and it is the object that owns the one
  # level set every block over the factor shares
  list(obj = fz$obj, opt = opt, plan = fz$plan, ess = ess, io = fz$io,
       lay = lay,
       rounds = length(moves), moved = moves[length(moves)], moves = moves,
       seed = seed, grad = grad, grads = grads,
       start_value = start_value, start_ess_min = min(start_ess$ess),
       capped = length(moves) >= rounds &&
         moves[length(moves)] >= imp_move_tol &&
         !(is.finite(grad) && grad < control$grad_tol))
}

#' How far the parameters may move between rounds and still count as a
#' fixed point. Fixed rather than a control knob because there is
#' nothing to tune: 1e-3 in natural parameter units is far below the
#' standard error of any parameter in a model small enough to need this
#' correction, and far below the Monte Carlo standard error of the
#' objective itself, so a smaller tolerance would only buy rounds that
#' chase noise.
#'
#' @noRd
imp_move_tol <- 1e-3

#' The effective sample size, as a fraction of the draw count, below
#' which the proposal is reported as having stopped covering its
#' integrand.
#'
#' Placed by measurement, and stated against the anchor the fit
#' actually reports from (the corrected estimate, not the Laplace
#' optimum; the two differ, and on a hard design by a lot).
#'
#' Ordinary designs sit far above it: the probe (60 groups of 8
#' Bernoulli rows, correlated slope) reports 0.93, the scalar-intercept
#' design 0.97, a gaussian response exactly 1.00, and the design in
#' `vignette("diagnostics")` 0.64.
#'
#' A proposal that has genuinely stopped working sits far below it:
#' displacing the probe's two log standard deviations by half a unit,
#' already worth about half a unit of log-likelihood, drops the worst
#' group to 0.009, some thirty times under.
#'
#' Between those two regimes the threshold is meant to FIRE, and does.
#' Forty groups of three binary rows with a standard deviation of two
#' is genuinely hard, and over seeds 1, 2, 3, 7, 11, 101 and 601 at a
#' thousand draws it reports 0.49, 0.35, 0.08, 0.16, 0.60, 0.62 and
#' 0.32: the warning fires on the hard half and stays quiet on the
#' rest. That is the diagnostic working, not a false alarm, and it is
#' why the number is not calibrated to clear every design anyone can
#' build.
#'
#' @noRd
imp_ess_floor <- 0.25


#' What the fit remembers about its correction: enough to say how it
#' was computed (draws, seed, rounds), how far it can be trusted (the
#' effective sample sizes and the Monte Carlo standard error), and
#' whether the rounds settled.
#'
#' @noRd
imp_record <- function(imp) {
  if (is.null(imp)) return(NULL)
  list(draws = imp$plan[["n_draw"]], seed = imp$seed, rounds = imp$rounds,
       moved = imp$moved, moves = imp$moves, capped = imp$capped,
       grad = imp$grad, grads = imp$grads,
       start_value = imp$start_value,
       ess_min = min(imp$ess$ess), ess_median = stats::median(imp$ess$ess),
       ess = stats::setNames(imp$ess$ess, imp$lay[["levels"]]),
       mcse = imp$ess$mcse)
}

#' Report the groups whose proposal has degenerated. A warning and not
#' an error: the estimate is still unbiased, it is only noisy, and the
#' remedy (more draws) is the user's to weigh.
#'
#' Named by the GROUPING FACTOR and not by a term label, because one
#' proposal covers a level's coefficients from every block over that
#' factor at once: there is one effective sample size per level, not
#' one per block.
#'
#' @noRd
imp_ess_warning <- function(ess, lay, n_draw, floor_at) {
  bad <- which(ess < floor_at)
  if (!length(bad)) return(invisible(NULL))
  shown <- utils::head(bad, 5L)
  warning("The importance proposal covers ", length(bad), " of ",
          length(ess), " groups of `", lay[["group_name"]],
          "` poorly: ", paste0("'", lay[["levels"]][shown], "' (",
                               format(ess[shown], digits = 2), ")",
                               collapse = ", "),
          if (length(bad) > length(shown)) ", ...",
          ". These are effective sample sizes as a fraction of the ",
          n_draw, " draws, and below ", floor_at,
          " the group's estimate rests on a handful of them. Raise ",
          "`importance`, or fit the Laplace model first and check that ",
          "it converged", call. = FALSE)
  invisible(NULL)
}

#' The one line `print()` and `summary()` give the correction. It says
#' what was computed (how many draws, how many rounds) and how far it
#' can be trusted (the Monte Carlo standard error of the reported
#' log-likelihood, and the worst group's effective sample size), so a
#' reader never has to take the corrected number on faith.
#'
#' @noRd
imp_report_line <- function(imp) {
  paste0("Marginal likelihood: importance-corrected, ", imp[["draws"]],
         " draws per group in ", imp[["rounds"]],
         if (imp[["rounds"]] == 1L) " round" else " rounds",
         " (MCSE ", format(imp[["mcse"]], digits = 2),
         ", min ESS ", format(imp[["ess_min"]], digits = 2),
         " of 1)")
}

#' How much the corrected objective may RISE across the rounds and
#' still count as noise rather than divergence.
#'
#' Measured on nine designs. Every fit whose iteration settled came out
#' the same or better (deltas -0.128 to -0.009), and the only rises
#' among them were the gaussian designs' optimizer drift, at most
#' +0.010. The two diverging fits rose by +0.352 and +3.227. The floor
#' of 0.1 sits ten times above the worst legitimate rise and three
#' times below the smallest divergence, and it is also a tenth of the
#' log-likelihood unit that one AIC unit is worth, so nothing that
#' could change a model comparison passes it. The Monte Carlo term
#' covers a low draw count, where the two estimates are noisy in their
#' own right.
#'
#' @noRd
imp_worse_tol <- function(mcse_final, mcse_start) {
  max(3 * mcse_final, 3 * mcse_start, 0.1)
}

#' A corrected fit's own frozen proposal, rebuilt from the fit.
#'
#' The proposal is a deterministic function of the frame, the reported
#' estimates, the draw count and the seed, all of which the fit
#' records, so rebuilding it reproduces the very proposal the fit
#' reported from rather than a fresh one. Checked: the effective sample
#' sizes it gives back at the estimate equal `fit$importance$ess_min`
#' exactly.
#'
#' Rebuilt rather than stored because the plan holds the whole draw
#' matrix, which is the largest object in the fit by far and would have
#' to be serialized with every corrected model.
#'
#' @noRd
imp_frozen_proposal <- function(fit) {
  fr <- fit[["frame"]]
  lay <- imp_layout(fr[["re_blocks"]])
  tpl <- fit[["estimates"]]
  lap <- RTMB::MakeADFun(build_objective(fr), tpl, random = "b",
                         map = fr[["map"]], silent = TRUE)
  plan <- imp_plan(lap, fr, lay, fit$opt$par, fit$importance$draws,
                   fit$importance$seed)
  list(io = build_importance_objective(fr, lay, imp_group_map(fr, lay), plan),
       plan = plan, template = tpl)
}

#' The effective sample sizes of that frozen proposal, read at an
#' arbitrary outer parameter vector. One numeric evaluation of the
#' log-weight matrix, no tape.
#'
#' @noRd
imp_ess_at <- function(prop, at) {
  imp_ess(as.matrix(prop$io$amat(imp_par_list(prop$template, at))),
          prop$plan[["n_draw"]])
}

#' Warn when a profile bound was computed where the frozen proposal no
#' longer covers the integrand.
#'
#' `confint(method = "profile")` walks one parameter away from the
#' estimate and reoptimizes the rest ON THE FIT'S OWN TAPE. That tape
#' carries a proposal frozen at the estimate, so the further the walk
#' goes the worse the weights get, and a covariance parameter can walk
#' far: measured, `theta_1` moved by 0.43 and the worst group's
#' effective sample size fell to 0.09 while the fit's reported 0.93
#' described the anchor and said nothing about the path.
#'
#' The proposal is NOT refrozen along the path, deliberately. Doing so
#' would make the profiled objective a different random function at
#' every grid point, and the interpolation the profile uses to find its
#' crossing assumes one smooth deterministic curve, which is exactly
#' what freezing the draws buys. Measured on a 40-group design, it
#' would also cost at least three times the profile it replaces. So the
#' curve stays frozen and the bound is reported with a warning when the
#' weights say it should not be trusted.
#'
#' @noRd
imp_profile_ess_warn <- function(fit, prop, parname, i, bounds) {
  floor_at <- fit$control$importance_ess %||% imp_ess_floor
  worst <- Inf
  where <- NA_real_
  for (b in bounds) {
    if (!is.finite(b)) next
    at <- fit$opt$par
    at[i] <- b
    e <- min(imp_ess_at(prop, at)$ess)
    if (is.finite(e) && e < worst) {
      worst <- e
      where <- b
    }
  }
  if (!is.finite(worst) || worst >= floor_at) return(invisible(NULL))
  warning("The profile bound ", format(where, digits = 6), " for '",
          parname, "' lies where this fit's importance proposal has ",
          "stopped covering the integrand: the worst group holds ",
          format(worst, digits = 2), " of its draws there, against a ",
          "threshold of ", floor_at, ". The proposal is frozen at the ",
          "estimate, so a profile that walks far from it is reweighting ",
          "with the wrong Gaussian and this bound is not reliable. Use ",
          "confint(method = \"wald\"), or refit with more draws",
          call. = FALSE)
  invisible(NULL)
}
