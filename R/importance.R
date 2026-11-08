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
# A Gaussian importance proposal cannot break that way. It stays a valid
# proposal wherever theta walks and only loses efficiency, and the
# effective sample size measures how much. Any block dimension, and
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
# only near the anchor. The standard normals are numeric constants baked
# into the closure, so the tape's parameters are the fixed effects, the
# dispersion and the covariance parameters alone: `random = NULL`, no
# inner problem, and sdreport works the ordinary way. u_hat_g and L_g
# are functions of theta on the tape, to first order around the anchor
# (build_importance_objective() says why). The Jacobian of the map from
# z to u is L_g(theta), so log|L_g| is on the tape with them.
#
# STACKED, NEVER LOOPED. For draw i every group's block is set to its
# i-th draw at once. The G x N draws of each coefficient reach the n x N
# predictor through one index, the family's density is called ONCE on a
# length n*N vector with the response and the addition terms recycled,
# and the per-group aggregation is one sparse indicator product. A loop
# over draws would put N copies of the likelihood on the tape as N
# separate calls.

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
imp_crosslevel <- function() {
  ents <- cs_entries()[cs_names()]
  c("equalto", names(ents)[vapply(ents, function(e) {
    isTRUE(e[["levels_coupled"]])
  }, TRUE)])
}

#' Every refusal `importance =` owes the user, raised before any tape
#' is built so the message names the conflict rather than a downstream
#' symptom. The style follows the `quadrature = TRUE` guards in
#' R/fit.R: say what cannot be combined, and say what to do instead.
#'
#' @noRd
check_importance_scope <- function(spec, frame, template, REML, quadrature,
                                   control) {
  if (isTRUE(quadrature)) {
    frm_stop("`importance` and `quadrature = TRUE` are two different ",
             "corrections of the same approximation and cannot both be ",
             "asked for. Adaptive Gauss-Kronrod marginalizes one scalar ",
             "random intercept at a time; the importance correction ",
             "reweights the Laplace Gaussian and takes any block ",
             "dimension. Keep one", call. = FALSE)
  }
  if (length(frame_latres(frame))) {
    frm_stop("`importance` cannot be combined with a latent residual: the ",
             "correction draws one proposal per group of a grouping factor, ",
             "and these latents are a chain with no such factor. Fit with ",
             "importance = 0", call. = FALSE)
  }
  if (frame_has_mo_re(frame)) {
    frm_stop("`importance` cannot be combined with a group-level mo() ",
             "term: the correction forms Z u once, off the tape, and that ",
             "term's column moves with its simplex. Fit with importance = 0",
             call. = FALSE)
  }
  if (REML) {
    frm_stop("`importance` cannot be combined with REML = TRUE. The ",
             "correction reweights the integral over the random effects ",
             "only, and REML puts the fixed effects into that same ",
             "integral, where the reweighting has no proposal for them. ",
             "Use REML = FALSE", call. = FALSE)
  }
  if (isTRUE(control$profile)) {
    frm_stop("`importance` cannot be combined with ",
             "frmtmb_control(profile = TRUE). Profiling moves the fixed ",
             "effects into the inner Laplace problem, and the corrected ",
             "objective has no inner problem left to move them into. Drop ",
             "profile = TRUE", call. = FALSE)
  }
  if (!is.null(template[["miss"]])) {
    frm_stop("`importance` cannot be combined with ",
             if (is.null(frame[["me"]])) "mi(). The imputed" else
               "me() or mi(). The latent",
             " values are latent variables of their own, with no grouping ",
             "factor to give them a per-group proposal. Fit with ",
             "importance = 0", call. = FALSE)
  }
  if (length(frame[["autocor"]] %||% list())) {
    frm_stop("`importance` cannot be combined with the residual ",
             "correlation term ", frame[["autocor"]][[1L]]$label,
             ": the correction resamples a random effect against a ",
             "PRODUCT of per-row densities, and this residual is one joint ",
             "density over each group (cov = TRUE) or gives each row a ",
             "mean that reads the residuals of earlier rows (cov = FALSE), ",
             "so no per-row integrand exists. Use importance = 0, or ",
             "REML = TRUE", call. = FALSE)
  }
  if (isTRUE(spec$rescor)) {
    frm_stop("`importance` cannot be combined with rescor: the responses ",
             "share one multivariate normal density per row, which is not ",
             "the rowwise density the correction resamples. Fit the ",
             "responses without rescor, or use importance = 0",
             call. = FALSE)
  }
  if (length(spec$responses) != 1L) {
    frm_stop("`importance` supports one response, and this model has ",
             length(spec$responses), " (",
             paste(names(spec$responses), collapse = ", "), "). A ",
             "multivariate model spreads its groups over several ",
             "likelihood terms, which the first version does not gather. ",
             "Fit the responses one at a time, or use importance = 0",
             call. = FALSE)
  }
  for (rn in names(spec$responses)) {
    st <- fam_structure(spec$responses[[rn]]$family)
    if (is.null(st[["loglik"]])) next
    # The correction needs the PIECES of the likelihood, one per
    # grouping level; a whole-response loglik is their sum with the
    # pieces thrown away. A family that HAS them declares one of the
    # factorization slots, and this refusal is for one that does not.
    if (is.null(st[["loglik_row"]]) && is.null(st[["loglik_group"]])) {
      frm_stop("`importance` cannot correct the '",
               spec$responses[[rn]]$family[["family"]], "' family: it ",
               "supplies its own log-likelihood, which returns one number ",
               "for the whole response, so a group's rows have no ",
               "separable integrand to resample. A family whose likelihood ",
               "does factorize over its groups says so with ",
               "frmtmb_structure(loglik_group = ) or (loglik_row = ), and ",
               "this one declares neither. Use importance = 0",
               call. = FALSE,
               package = frm_family_package(spec$responses[[rn]]$family))
    }
  }
  for (lp in frame[["linpreds"]]) {
    if (!is.null(lp[["nl_body"]])) {
      frm_stop("`importance` cannot correct a nonlinear predictor (",
               lp[["resp"]], ".", lp[["dpar"]], "). A nonlinear body mixes ",
               "parameter values with raw data columns, and the corrected ",
               "objective evaluates the predictor once per draw, where a ",
               "column would recycle silently against a longer vector. Use ",
               "importance = 0", call. = FALSE)
    }
    if (length(lp[["cs"]] %||% list())) {
      frm_stop("`importance` cannot correct a cs() term (", lp[["resp"]],
               ".", lp[["dpar"]], "): its threshold-specific offsets are a ",
               "matrix per observation, which the first version does not ",
               "stack over draws. Use importance = 0", call. = FALSE)
    }
  }
  if (!is.null(frame[["map"]][["b"]])) {
    frm_stop("`importance` cannot correct a model whose random-effect ",
             "coefficients are mapped: the proposal is built per grouping ",
             "level from the conditional Hessian, and a map merges or ",
             "fixes coefficients across levels. Use importance = 0",
             call. = FALSE)
  }
  blocks <- frame[["re_blocks"]] %||% list()
  if (!length(blocks)) {
    frm_stop("`importance` needs a random-effect block to correct, and this ",
             "model has none. The correction reweights an integral over ",
             "random effects, and a model without any has no such integral ",
             "and no approximation to improve on. Use importance = 0",
             call. = FALSE)
  }
  labs <- vapply(blocks, function(bk) bk[["term_label"]], "")
  for (bk in blocks) {
    if (is.null(bk[["levels"]])) {
      frm_stop("`importance` needs a grouping factor, and `",
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
    frm_stop("`importance` takes several random-effect blocks only when ",
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
    frm_stop("`importance` needs every block over `", grps[[1L]],
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
    if (bk[["covstruct"]] %in% imp_crosslevel()) {
      frm_stop("`importance` cannot correct the '", bk[["covstruct"]],
               "' structure in `", bk[["term_label"]], "`: it correlates the ",
               "grouping LEVELS with each other through a supplied ",
               "relationship or neighbor matrix, so the marginal likelihood ",
               "is one integral over every level at once and does not split ",
               "into per-group ones. Use importance = 0", call. = FALSE)
    }
    if (!bk[["covstruct"]] %in% imp_covstructs) {
      frm_stop("`importance` cannot correct the '", bk[["covstruct"]],
               "' structure in `", bk[["term_label"]], "`. The correction ",
               "recovers each level's prior density from the block density, ",
               "which needs that density to be Gaussian in the level's ",
               "coefficients and independent between levels; the supported ",
               "structures are ",
               paste(imp_covstructs, collapse = ", "),
               ". Use importance = 0", call. = FALSE)
    }
  }
  # LAST, because it is the only check that needs the group map, and
  # building the map has refusals of its own that belong after every
  # cheaper one. A structured family reaching here declares a
  # factorization slot; what is left to check is that the units it
  # factorizes into are the ones the correction integrates over.
  for (rn in names(spec$responses)) {
    fam <- spec$responses[[rn]]$family
    st <- fam_structure(fam)
    if (is.null(st[["loglik"]])) next
    imp_check_family_grouping(frame, fam, rn)
  }
  invisible(NULL)
}

#' The family's independent units and the correction's groups must be
#' the SAME partition of the rows.
#'
#' They are different objects that coincide in the case the seam was
#' built for: `rw_delta(subject = id)` under `(1 | id)` gives one family
#' unit per grouping level. The same family under `(1 | item)` gives a
#' per-SUBJECT likelihood against a per-ITEM proposal, and summing the
#' first into the second is not defined. Getting it wrong is silent -
#' the arithmetic goes through and the numbers are meaningless - so the
#' check is required rather than advisory.
#'
#' Equality, not refinement, in either direction. A family unit spanning
#' two groups has no separable integrand; a group holding two family
#' units would have the family run one recursion over both, which
#' concatenates two sequences that never met.
#'
#' @noRd
imp_check_family_grouping <- function(frame, fam, rn) {
  blk <- frame_block_of(frame, rn)
  codes <- structure_group_codes(blk)
  lay <- imp_layout(frame[["re_blocks"]])
  labs <- vapply(lay[["blocks"]], function(bk) bk[["term_label"]], "")
  if (is.null(codes)) {
    frm_stop("`importance` cannot correct the '", fam[["family"]],
             "' family here: it declares a factorization of its likelihood ",
             "but its frame block carries no `group`, so there is nothing ",
             "to check the grouping of `", labs[[1L]], "` against. The ",
             "correction would be summing the family's pieces into groups ",
             "it has no reason to belong to. Use importance = 0",
             call. = FALSE,
             package = frm_family_package(fam))
  }
  rl <- imp_group_map(frame, lay)[["row_level"]]
  pairs <- length(unique(paste(codes, rl, sep = ":")))
  if (pairs == length(unique(codes)) && pairs == length(unique(rl))) {
    return(invisible(NULL))
  }
  frm_stop("`importance` needs the '", fam[["family"]], "' family's own ",
           "units and the grouping of `", labs[[1L]],
           "` to be the same partition of the rows, and they are not: the ",
           "family has ", length(unique(codes)), " unit(s) against ",
           length(unique(rl)), " grouping level(s), covering ", pairs,
           " combination(s). The correction resamples one grouping level ",
           "at a time and adds up the family's pieces inside it, which ",
           "means something only when a piece belongs to exactly one ",
           "level. Group the model on what the family groups on, or use ",
           "importance = 0", call. = FALSE,
           package = frm_family_package(fam))
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
    frm_stop("`importance` needs every observation to belong to one ",
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

#' The proposal, anchored at one outer parameter vector: conditional
#' modes, per-level Hessian blocks, their Cholesky factors, the draws
#' they map the standard normals through, and the rates at which the
#' tape moves the proposal away from the anchor.
#'
#' The modes and the Hessian come from the Laplace objective's own
#' inner Newton solve, so the proposal AT THE ANCHOR is the Laplace
#' Gaussian and nothing has to be re-derived. `spHess(random = TRUE)` is the
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
    frm_stop("`importance` expected ", nb, " conditional modes for `",
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
      frm_stop("`importance` needs the conditional Hessian of `",
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
  hess_b <- imp_block_filler(hess, lay)(hess)
  # The Hessian splits in two: the prior's precision, which is a known
  # function of theta and goes on the tape, and the likelihood's
  # curvature at the mode. A likelihood that is not log-concave can have a
  # negative direction here, and the proposal then has its negative part
  # cut off so that any theta the optimizer tries still has a positive
  # definite Gaussian.
  theta_a <- imp_par_list(frame[["par_template"]], par,
                          frame[["map"]])[["theta"]]
  prior_t <- imp_prior_scales(lay, theta_a)
  prior_fn <- imp_prior_terms(lay, prior_t)
  prec <- imp_prior_precision(prior_fn(theta_a), q)
  hlik <- array(0, c(q, q, ng))
  if (q == 1L) {
    # one scalar block and nothing else: H is a diagonal, so L is one
    # reciprocal square root per level and no factorization is needed.
    # Still addressed through `idx`, so the fast path reads the same
    # index as the general one rather than re-deriving the layout.
    p1 <- idx[1L, ]
    hd <- hess_b[1L, 1L, ]
    if (any(!is.finite(hd)) || any(hd <= 0)) {
      frm_stop("`importance` found a non-positive conditional variance for ",
               sum(!is.finite(hd) | hd <= 0), " level(s) of `",
               lay[["label"]], "`, so no Gaussian proposal exists ",
               "there. The Laplace fit is at a singular variance ",
               "component; refit it first", call. = FALSE)
    }
    hlik[1L, 1L, ] <- pmax(hd - prec[1L, 1L], 0)
    sdv <- 1 / sqrt(hlik[1L, 1L, ] + prec[1L, 1L])
    u[p1, ] <- uhat[p1] + sdv * z[p1, , drop = FALSE]
    logdet_l <- log(sdv)
  } else {
    for (k in seq_len(ng)) {
      ii <- idx[, k]
      hk <- hess_b[, , k]
      pk <- tryCatch(chol(hk), error = function(e) NULL)
      if (is.null(pk)) {
        frm_stop("`importance` could not factor the conditional Hessian of ",
                 "level '", lay[["levels"]][k], "' in `", lay[["label"]],
                 "`: it is not positive definite, so that group has no ",
                 "Gaussian proposal. The Laplace fit sits on a singular ",
                 "variance component; refit it first", call. = FALSE)
      }
      hl <- imp_psd_part(hk - prec)
      hlik[, , k] <- hl
      pk <- chol(hl + prec)
      # L = P^-1 with P'P = H gives L L' = H^-1, upper triangular and
      # backsolved rather than inverted
      lk <- backsolve(pk, diag(q))
      u[ii, ] <- uhat[ii] + lk %*% z[ii, , drop = FALSE]
      logdet_l[k] <- -sum(log(diag(pk)))
    }
  }
  sl <- imp_anchor_slopes(lap_obj, full, lay, hlik, prec, frame, prior_fn)
  # 0.5 |z_{g,i}|^2, one value per group and draw. Row j of `idx` is
  # the j-th coefficient of every group at once, which is what makes
  # this a q-term loop over G x N matrices and not a loop over groups.
  zsq <- matrix(0, ng, nd)
  for (j in seq_len(q)) {
    zsq <- zsq + z[idx[j, ], , drop = FALSE]^2
  }
  zsq <- 0.5 * zsq
  list(u = u, z = z, zsq = zsq, logdet_l = logdet_l, n_draw = nd, par = par,
       uhat = uhat, q = q, n_group = ng, levels = lay[["levels"]],
       idx = idx, hlik = hlik, prec = prec, prior_t = prior_t,
       cross_d = sl[["cross"]],
       hlik_rate = sl[["rate"]])
}

#' The per-group blocks of a block-diagonal sparse matrix, as a
#' `q x q x G` array.
#'
#' `imp_anchor_slopes()` reads this matrix several times per coordinate, and
#' its sparsity pattern does not change, so the map from its entries to the
#' array is built once here and each read fills it from `@x`: the triplet
#' conversion is the part worth not repeating.
#'
#' @noRd
imp_block_filler <- function(m, lay) {
  q <- lay[["qt"]]
  ng <- lay[["n_group"]]
  # where a position sits inside its group, which `idx` holds by column
  within <- integer(lay[["nb"]])
  within[as.vector(lay[["idx"]])] <- rep(seq_len(q), ng)
  group <- lay[["pos_group"]]
  # the entries of a compressed column matrix in the order @x holds them
  i <- m@i + 1L
  j <- rep(seq_len(ncol(m)), diff(m@p))
  same <- group[i] == group[j]
  dims <- c(q, q, ng)
  lin <- (group[i[same]] - 1L) * q * q + (within[j[same]] - 1L) * q +
    within[i[same]]
  mirror <- (group[i[same]] - 1L) * q * q + (within[i[same]] - 1L) * q +
    within[j[same]]
  sym <- methods::is(m, "symmetricMatrix")
  function(mm) {
    out <- array(0, dims)
    out[lin] <- mm@x[same]
    if (sym) out[mirror] <- mm@x[same]
    out
  }
}

#' The symmetric part of a small matrix with its negative eigenvalues
#' set to zero. Returned untouched when it has none, so a log-concave
#' likelihood keeps the Laplace Hessian to the last bit.
#'
#' @noRd
imp_psd_part <- function(m) {
  m <- 0.5 * (m + t(m))
  e <- eigen(m, symmetric = TRUE)
  if (all(e$values >= 0)) return(m)
  e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
}

#' The prior's precision matrix of one group, from `imp_prior_terms()`.
#'
#' @noRd
imp_prior_precision <- function(terms, q) {
  out <- matrix(0, q, q)
  for (j in seq_len(q)) out[j, j] <- -2 * terms[["diag"]][[j]]
  for (e in terms[["off"]]) {
    out[e[["j"]], e[["k"]]] <- out[e[["k"]], e[["j"]]] <- -e[["value"]]
  }
  out
}

#' How the mode's score and the likelihood's curvature move with the
#' outer parameters, at the anchor.
#'
#' `cross` is `d^2 f / du dpar`, `nb x p`, the likelihood's slope for the
#' one Newton step the proposal takes. Its covariance-parameter columns
#' stay zero: those parameters reach the joint density through the prior
#' alone, whose precision is on the tape at its exact value.
#'
#' `rate` is `G x p`, the rate at which each group's likelihood curvature
#' scales as the parameter moves AND THE MODE MOVES WITH IT. The total
#' derivative is the one that matters, because the proposal sits at the
#' moved mode: a step in a log-link intercept moves the mode by minus that
#' step, the predictor at the mode stays where it is, and the curvature
#' there does not change, where its partial derivative equals itself
#' (measured: taking the partial halved the tape's curvature along the
#' intercept, 4.1 against 8.1, and left two fits of 648 rising). It is
#' the least-squares scalar that takes the group's `Hlik` to that
#' derivative, so that the curvature stays positive semi-definite at any
#' theta. On a gaussian response a dispersion scales it exactly, at -2 per
#' log sigma, and nothing else moves it.
#'
#' A difference of the AD gradient and the sparse Hessian, because the
#' tape offers those cheaply and no mixed derivative at all, and because
#' nothing here has to be exact: a proposal that does not depend on the
#' weights leaves the estimator unbiased, so an error here costs efficiency
#' and never validity. For the same reason a coordinate whose difference is
#' not finite stays at zero, which leaves the proposal where it is along
#' that coordinate, rather than a refusal. A fixed effect's difference is
#' one-sided, since at the step below it reaches about one part in 1e4 and
#' costs half the Hessians, which set the cost of a freeze at many outer
#' parameters. A dispersion's stays central: a gaussian fit's fixed point
#' is exact only as far as its rate of -2 is.
#'
#' @noRd
imp_anchor_slopes <- function(lap_obj, full, lay, hlik, prec, frame, prior_fn) {
  env <- lap_obj$env
  rnd <- env$random
  fix <- setdiff(seq_along(full), rnd)
  q <- lay[["qt"]]
  ng <- lay[["n_group"]]
  idx <- lay[["idx"]]
  is_theta <- names(full)[fix] == "theta"
  is_beta <- names(full)[fix] == "beta"
  theta_of <- function(at) unname(at[fix][is_theta])
  # f() and spHess() record the point they read, and a later reader of
  # last.par wants the anchor the inner solve left there
  on.exit(env$last.par <- full, add = TRUE)
  grad_u <- function(at) as.numeric(env$f(at, order = 1L))[rnd]
  fill <- imp_block_filler(env$spHess(full, random = TRUE), lay)
  prec_v <- rep(as.vector(prec), ng)
  # the likelihood's part of the group blocks at a point, the same split
  # imp_plan() makes at the anchor. Only theta moves the prior's part.
  hlik_at <- function(at, theta_moved) {
    pv <- if (theta_moved) {
      rep(as.vector(imp_prior_precision(prior_fn(theta_of(at)), q)), ng)
    } else {
      prec_v
    }
    matrix(fill(env$spHess(at, random = TRUE)) - pv, q * q, ng)
  }
  # H^-1 of each group, once: every coordinate solves against it
  hinv <- array(0, c(q, q, ng))
  for (k in seq_len(ng)) hinv[, , k] <- chol2inv(chol(hlik[, , k] + prec))
  hsolve <- function(d) {
    out <- numeric(length(d))
    for (j in seq_len(q)) {
      for (k in seq_len(q)) {
        out[idx[j, ]] <- out[idx[j, ]] + hinv[j, k, ] * d[idx[k, ]]
      }
    }
    out
  }
  cross <- matrix(0, length(rnd), length(fix))
  rate <- matrix(0, ng, length(fix))
  hl <- matrix(hlik, q * q, ng)
  den <- colSums(hl^2)
  steps <- imp_outer_steps(frame, full[fix])
  g0 <- if (any(is_beta)) grad_u(full)
  h0 <- if (any(is_beta)) hlik_at(full, FALSE)
  for (j in seq_along(fix)) {
    h <- steps[[j]]
    up <- full
    up[[fix[[j]]]] <- up[[fix[[j]]]] + h
    if (is_beta[[j]]) {
      d <- (grad_u(up) - g0) / h
      if (!all(is.finite(d))) next
      cross[, j] <- d
      up[rnd] <- up[rnd] - hsolve(d) * h
      r <- colSums((hlik_at(up, FALSE) - h0) * hl) / (h * den)
    } else {
      dn <- full
      dn[[fix[[j]]]] <- dn[[fix[[j]]]] - h
      d <- (grad_u(up) - grad_u(dn)) / (2 * h)
      if (!all(is.finite(d))) next
      if (!is_theta[[j]]) cross[, j] <- d
      du <- -hsolve(d) * h
      up[rnd] <- up[rnd] + du
      dn[rnd] <- dn[rnd] - du
      r <- colSums((hlik_at(up, is_theta[[j]]) -
                      hlik_at(dn, is_theta[[j]])) * hl) / (2 * h * den)
    }
    r[!is.finite(r)] <- 0
    rate[, j] <- r
  }
  list(cross = cross, rate = rate)
}

#' Each outer parameter's step for those differences, in units of the
#' predictor it moves.
#'
#' A coefficient of column `x` moves the predictor by `h |x|`, so the step
#' that keeps the predictor's move at `imp_cross_step` is that over the
#' largest `|x|`, whatever the column's center and scale. On `year =
#' round(2000 + 5 z)` a step of `1e-4` moves the predictor by 0.2, the slope's
#' cross derivative comes out wrong, and with it the proposal off the anchor:
#' slope standard errors up to 1.47 times the centered fit's, or capped
#' rounds, or a refusal (measured over six data seeds).
#'
#' A predictor in a response's units moves in those units. Without that, a
#' Student fit on `y x 1e-6` steps by 100 of its own scale, its slopes come out
#' wrong, and it stops on 2 of 4 data seeds.
#'
#' No step reads a parameter's own size. A parameter with no column moves
#' the predictor or a density's scale one for one or not at all, so its size
#' says nothing about how far it may step. An ordinal's first threshold takes
#' a column's center: at about 1600, a step of its own size moved the predictor
#' by 0.16, and the slope's SE read 4.56 times the centered fit's, silently.
#'
#' @noRd
imp_outer_steps <- function(frame, at) {
  h <- rep(imp_cross_step, length(at))
  pos <- imp_outer_pos(frame[["par_template"]], at, frame[["map"]])
  units <- predictor_units(frame)
  scale <- lapply(frame[["par_template"]], function(v) rep(NA_real_, length(v)))
  for (key in names(frame[["linpreds"]])) {
    lp <- frame[["linpreds"]][[key]]
    cp <- lp[["par"]]
    X <- lp[["X"]]
    if (is.null(scale[[cp]]) || is.null(X) || !ncol(X)) next
    if (length(lp[["idx"]]) != ncol(X)) next
    biggest <- imp_col_max(X)
    biggest[!is.finite(biggest) | biggest <= 0] <- 1
    u <- if (key %in% names(units)) units[[key]] else 1
    scale[[cp]][lp[["idx"]]] <- pmin(scale[[cp]][lp[["idx"]]], u / biggest,
                                     na.rm = TRUE)
  }
  pn <- names(at)
  for (cp in names(pos)) {
    s <- scale[[cp]][pos[[cp]]]
    known <- !is.na(s)
    h[pn == cp][known] <- imp_cross_step * s[known]
  }
  h
}

#' The largest absolute entry of each column of a design, dense or sparse.
#'
#' @noRd
imp_col_max <- function(X) {
  if (!methods::is(X, "sparseMatrix")) {
    return(apply(abs(as.matrix(X)), 2L, max))
  }
  tx <- methods::as(X, "TsparseMatrix")
  out <- numeric(ncol(X))
  o <- order(abs(tx@x))
  # the last assignment to a column is its largest entry
  out[tx@j[o] + 1L] <- abs(tx@x[o])
  out
}

#' How far those differences move the predictor. Well below the scale a
#' density's curvature changes on and well above the rounding of its
#' gradient.
#'
#' @noRd
imp_cross_step <- 1e-4

#' Where each outer parameter sits in its template component. A `map`
#' fixes or shares entries, so a component's outer values are not its
#' entries one for one: the position is the first entry of each level
#' the map keeps, in level order, which is the order the tape holds them.
#'
#' @noRd
imp_outer_pos <- function(template, par, map = NULL) {
  pn <- names(par)
  out <- list()
  for (cp in unique(pn)) {
    mp <- map[[cp]]
    out[[cp]] <- if (is.null(mp)) {
      seq_along(template[[cp]])
    } else {
      fm <- factor(mp)
      match(seq_len(nlevels(fm)), as.integer(fm))
    }
    if (length(out[[cp]]) != sum(pn == cp)) {
      frm_stop("`importance` read ", sum(pn == cp), " outer values for `",
               cp, "`, whose template and map keep ", length(out[[cp]]),
               ". The two have to agree for the proposal to follow the ",
               "parameters", call. = FALSE)
    }
  }
  out
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
#' `cs_entry(k)$nll` returns the block's log density summed
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
imp_prior_terms <- function(lay, scales) {
  per <- lapply(seq_along(lay[["blocks"]]), function(m) {
    imp_block_prior(lay[["blocks"]][[m]], scales[[m]])
  })
  tix <- lapply(lay[["blocks"]], function(bk) bk[["theta_idx"]])
  q <- lay[["q"]]
  offset <- lay[["offset"]]
  function(theta) {
    c0 <- 0
    dg <- list()
    off <- list()
    for (m in seq_along(per)) {
      tm <- per[[m]](theta[tix[[m]]])
      c0 <- c0 + tm[["c0"]]
      dg <- c(dg, tm[["diag"]])
      for (e in tm[["off"]]) {
        e[["j"]] <- e[["j"]] + offset[[m]]
        e[["k"]] <- e[["k"]] + offset[[m]]
        off[[length(off) + 1L]] <- e
      }
    }
    list(c0 = c0, diag = dg, off = off)
  }
}

#' Where to read each coefficient's precision: per block, a `t` per
#' coefficient at which the block's density at `t e_j` lies about one nat
#' below its value at zero, at the anchor's theta.
#'
#' The precision is a difference from the normalizer, and at `t = 1` that
#' difference rounds once `1 / sd^2` nears the normalizer's `eps`: relative
#' error 1.9e-3 at sd 1e6, 0.29 at 1e7, and it reads 0 at 1e8. A fit on
#' `y x 1e6` then stopped with a request to report the model. The block is
#' Gaussian in its coefficients, so one density at any `t` reads the same
#' number. A search on doubles, off the tape, once per anchor.
#'
#' @noRd
imp_prior_scales <- function(lay, theta) {
  lapply(lay[["blocks"]], function(bk) {
    bfn <- cs_entry(bk[["covstruct"]])$nll
    q <- bk[["dim"]]
    bk1 <- bk
    bk1[["n_levels"]] <- 1L
    th <- theta[bk[["theta_idx"]]]
    zero <- numeric(q)
    c0 <- bfn(zero, th, bk1)
    vapply(seq_len(q), function(j) {
      t0 <- 1
      for (i in seq_len(20L)) {
        v <- c0 - bfn(replace(zero, j, t0), th, bk1)
        if (!is.finite(v) || v > 1e-2) break
        t0 <- t0 * 1e3
      }
      if (is.finite(v) && v > 0) t0 / sqrt(v) else 1
    }, 0)
  })
}

#' One block's contribution to that sum, on the coefficient slices it
#' owns. Factored out of `imp_prior_terms()` so the extraction is
#' written once whatever the block count.
#'
#' It needs no `ADoverload()` of its own: the arithmetic here is `*`
#' and `+` on advectors, and `cs_entry(k)$nll` installs the
#' overloads it needs itself (R/covstruct.R). That is worth stating
#' rather than assuming, because an overload lost in a factoring-out
#' fails far from its cause.
#'
#' @noRd
imp_block_prior <- function(bk, t0) {
  bfn <- cs_entry(bk[["covstruct"]])$nll
  q <- bk[["dim"]]
  bk1 <- bk
  bk1[["n_levels"]] <- 1L
  zero <- numeric(q)
  unit <- lapply(seq_len(q), function(j) replace(zero, j, t0[[j]]))
  pair <- list()
  if (q > 1L) {
    for (j in seq_len(q - 1L)) {
      for (k in (j + 1L):q) {
        pair[[paste0(j, ".", k)]] <- replace(zero, c(j, k), t0[c(j, k)])
      }
    }
  }
  function(theta) {
    c0 <- bfn(zero, theta, bk1)
    # at[[j]] = -P[j, j] t_j^2 / 2
    at <- lapply(unit, function(e) bfn(e, theta, bk1) - c0)
    dg <- lapply(seq_len(q), function(j) at[[j]] / t0[[j]]^2)
    off <- list()
    if (q > 1L) {
      for (j in seq_len(q - 1L)) {
        for (k in (j + 1L):q) {
          # -P[j, k] t_j t_k, from the value at both less both diagonals
          off[[length(off) + 1L]] <- list(
            j = j, k = k,
            value = (bfn(pair[[paste0(j, ".", k)]], theta, bk1) - c0 -
                       at[[j]] - at[[k]]) / (t0[[j]] * t0[[k]]))
        }
      }
    }
    list(c0 = c0, diag = dg, off = off)
  }
}

#' The theta-only negative log-likelihood closure: `-sum_g log m_g`,
#' with the standard normals, the design values that carry them into
#' every linear predictor, and the per-group indicator all baked in as
#' constants.
#'
#' Returns the closure and, beside it, an `amat` that hands back the
#' G x N matrix of log importance weights at a numeric parameter list.
#' The diagnostics read that matrix, so the effective sample sizes and
#' the Monte Carlo standard error describe the objective that was
#' optimized and not a second computation of it.
#'
#' @noRd
build_importance_objective <- function(frame, lay, gmap, plan,
                                       neg_log_prior = NULL) {
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
  # the rows mi_ml_frame() took out, as the objective's weight has them
  gone <- frame[["mi_map"]][[rn]][["unread"]]
  if (length(gone)) {
    keep <- rep(1, n)
    keep[gone] <- 0
    wts <- wts * imp_expand(list(k = keep), ridx, n)$k
  }
  extra_names <- frame[["extra_names"]] %||% character(0)
  smat <- gmap[["S"]]
  # A structured family's own factorization, resolved HERE, off the
  # tape and once: the branch is on the family object and the frame
  # block, both of them data, so nothing below puts a condition on a
  # parameter value.
  st <- fam_structure(fam)
  blk <- frame_block_of(frame, rn)
  st_row <- st[["loglik_row"]]
  st_grp <- st[["loglik_group"]]
  # Group values come back in the family's own level order and the
  # correction counts groups in its own. The two orders describe the
  # same partition (check_importance_scope() refuses otherwise) and
  # need not agree, so the permutation between them is a sparse product
  # like the aggregation it replaces, formed once here.
  pmat <- NULL
  if (!is.null(st_grp)) {
    codes <- structure_group_codes(blk)
    rl <- gmap[["row_level"]]
    to <- rl[match(seq_len(max(codes)), codes)]
    pmat <- Matrix::sparseMatrix(i = to, j = seq_along(to),
                                 x = rep(1, length(to)),
                                 dims = c(ng, length(to)))
  }

  # The draws reach each linear predictor through Z, one coefficient at
  # a time. A row belongs to one group (imp_group_map() refuses a row that
  # does not), so coefficient j of the row's group has one design value on
  # that row, and the n x N predictor reads the G x N draws of j through
  # one index. A coefficient that is 1 on every row needs no product.
  rl <- gmap[["row_level"]]
  rlidx <- rep(rl, nd) + rep(ng * (seq_len(nd) - 1L), each = n)
  z_design <- lapply(lps, function(lp) {
    if (is.null(lp[["Z"]])) return(NULL)
    lapply(seq_len(q), function(j) {
      zv <- as.vector(lp[["Z"]][, plan[["idx"]][j, ], drop = FALSE] %*%
                        rep(1, ng))
      if (all(zv == 0)) NULL else if (all(zv == 1)) 1 else zv[ridx]
    })
  })
  # one G x N slice of standard normals per within-GROUP coefficient,
  # gathered by row j of the scattered index rather than by a stride,
  # because a group's coefficients are only contiguous when there is one
  # block
  z_slices <- lapply(seq_len(q), function(j) {
    plan[["z"]][plan[["idx"]][j, ], , drop = FALSE]
  })
  uhat_g <- lapply(seq_len(q), function(j) {
    plan[["uhat"]][plan[["idx"]][j, ]]
  })
  hlik <- plan[["hlik"]]
  prec_a <- plan[["prec"]]
  prior_fn <- imp_prior_terms(lay, plan[["prior_t"]])
  zsq <- plan[["zsq"]]
  # A per-group constant subtracted before the exponential, so that the
  # weights sit near 1 at the anchor. It cancels exactly from the
  # result, and being a CONSTANT it puts no branch on parameter values,
  # which a running maximum would.
  shift <- new.env(parent = emptyenv())
  shift$offs <- matrix(0, ng, nd)
  # Outer parameters that reach the likelihood, against the anchor.
  tmpl <- frame[["par_template"]]
  opos <- imp_outer_pos(tmpl, plan[["par"]], frame[["map"]])
  anchor <- unname(plan[["par"]])
  cross_d <- plan[["cross_d"]]
  hlik_rate <- plan[["hlik_rate"]]
  gidx <- rep.int(seq_len(ng), nd)

  # THE PROPOSAL FOLLOWS THE PARAMETERS. With the draws frozen in place, a
  # step d in an intercept moves every group's integrand by d under a
  # proposal that stays put. The curvature the optimizer then sees is the
  # group's conditional information times (1 - mean z^2), a sample
  # statistic far larger than the marginal curvature it has to resolve
  # (measured: exact 34, and 5 to 70 over seeds at 200 draws, negative at
  # one), and the re-anchored iteration diverges wherever it lands too
  # small, however good the weights are at the anchor. So the optimizer's
  # theta gets its own Gaussian: one Newton step from the anchor's mode
  # toward theta's, with the prior's precision exact and the likelihood's
  # curvature and slope carried to first order. That removes the (z^2 - 1)
  # terms from the gradient at the anchor: on a gaussian response the
  # weights stay equal to first order for every draw, and the fixed point
  # is the exact optimum at any draw count. At a given anchor the proposal
  # is a function of theta alone, so the estimate is unbiased at each
  # theta. From round 2 the anchor is the optimum over the same standard
  # normals, which measurably does not move the error. They stay fixed, so
  # the objective stays one smooth function.
  proposal <- function(pars) {
    "c" <- RTMB::ADoverload("c")
    tm <- prior_fn(pars[["theta"]])
    th <- pars[[names(opos)[[1L]]]][opos[[1L]]]
    for (cp in names(opos)[-1L]) th <- c(th, pars[[cp]][opos[[cp]]])
    dth <- th - anchor
    hscale <- exp(as.vector(hlik_rate %*% dth))
    prec <- function(j, k) {
      if (j == k) return(-2 * tm[["diag"]][[j]])
      for (e in tm[["off"]]) {
        if (e[["j"]] == j && e[["k"]] == k) return(-e[["value"]])
      }
      NULL
    }
    # the upper Cholesky factor R of H(theta), one G-vector per entry, and
    # its inverse: a q(q+1)/2 loop, never a loop over groups
    rr <- rinv <- vector("list", q * q)
    at <- function(j, k) (k - 1L) * q + j
    for (j in seq_len(q)) {
      for (k in j:q) {
        pjk <- prec(j, k)
        s <- hlik[j, k, ] * hscale
        if (!is.null(pjk)) s <- s + pjk
        for (l in seq_len(j - 1L)) s <- s - rr[[at(l, j)]] * rr[[at(l, k)]]
        rr[[at(j, k)]] <- if (k == j) sqrt(s) else s / rr[[at(j, j)]]
      }
    }
    for (j in rev(seq_len(q))) {
      rinv[[at(j, j)]] <- 1 / rr[[at(j, j)]]
      for (k in seq_len(q)[seq_len(q) > j]) {
        s <- rr[[at(j, j + 1L)]] * rinv[[at(j + 1L, k)]]
        for (l in seq_len(q)[seq_len(q) > j + 1L & seq_len(q) <= k]) {
          s <- s + rr[[at(j, l)]] * rinv[[at(l, k)]]
        }
        rinv[[at(j, k)]] <- -s * rinv[[at(j, j)]]
      }
    }
    logdet <- -log(rr[[at(1L, 1L)]])
    for (j in seq_len(q)[-1L]) logdet <- logdet - log(rr[[at(j, j)]])
    # the score of the joint density at the anchor's mode under theta: the
    # prior's moved precision on that mode, and the likelihood's slope
    lik <- as.vector(cross_d %*% dth)
    rhs <- lapply(seq_len(q), function(j) {
      s <- lik[plan[["idx"]][j, ]]
      for (k in seq_len(q)) {
        pjk <- prec(min(j, k), max(j, k))
        if (is.null(pjk)) next
        s <- s + (pjk - prec_a[j, k]) * uhat_g[[k]]
      }
      s
    })
    # H^-1 rhs through R^-1 R^-T, and then the draws. Only upper entries
    # exist, so each sum runs over them alone.
    half <- lapply(seq_len(q), function(k) {
      s <- rinv[[at(1L, k)]] * rhs[[1L]]
      for (j in seq_len(k)[-1L]) s <- s + rinv[[at(j, k)]] * rhs[[j]]
      s
    })
    u <- lapply(seq_len(q), function(j) {
      s <- rinv[[at(j, j)]] * half[[j]]
      for (k in seq_len(q)[seq_len(q) > j]) {
        s <- s + rinv[[at(j, k)]] * half[[k]]
      }
      center <- uhat_g[[j]] - s
      draw <- center[gidx] + rinv[[at(j, j)]][gidx] * z_slices[[j]]
      for (k in seq_len(q)[seq_len(q) > j]) {
        draw <- draw + rinv[[at(j, k)]][gidx] * z_slices[[k]]
      }
      draw
    })
    list(u = u, logdet = logdet, prior = tm)
  }

  weights <- function(pars) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    extra <- NULL
    if (length(extra_names)) {
      extra <- lapply(stats::setNames(extra_names, extra_names),
                      function(nm) pars[[nm]])
    }
    pr <- proposal(pars)
    u <- pr[["u"]]
    dpv <- list()
    for (k in seq_along(lps)) {
      lp <- lps[[k]]
      # the fixed part is the same for every draw, so it is built once
      # at length n and expanded; only Z u varies from draw to draw
      eta <- lp_eta_fixed(lp, pars, n, list(), frame[["y"]])[ridx]
      for (j in seq_along(z_design[[k]])) {
        zv <- z_design[[k]][[j]]
        if (is.null(zv)) next
        eta <- eta + zv * u[[j]][rlidx]
      }
      dpv[[lp[["dpar"]]]] <- lp[["link"]]$linkinv(eta)
      dpv[[paste0(".eta_", lp[["dpar"]])]] <- eta
    }
    agg <- if (!is.null(st_grp)) {
      # one value per family unit per draw, already the group sums the
      # rowwise path has to form; the family decides what a weight
      # means to a likelihood that is not rowwise, so it is passed and
      # not applied here
      pmat %*% RTMB::matrix(st_grp(yraw, dpv, av, wts, blk, extra),
                            ncol(pmat), nd)
    } else if (!is.null(st_row)) {
      smat %*% RTMB::matrix(st_row(yraw, dpv, av, wts, blk, extra), n, nd)
    } else {
      smat %*% RTMB::matrix(wts * row_lpdf(fam, yraw, yraw, dpv, av, extra),
                            n, nd)
    }
    # The sparse product gives an advector matrix on the tape and a
    # Matrix object in the numeric path, and a Matrix loses its
    # dimensions through exp(). The test is on CLASS, so it resolves
    # while the tape is being built and leaves no branch on it.
    if (!inherits(agg, "advector")) agg <- as.matrix(agg)
    tm <- pr[["prior"]]
    logp <- tm[["c0"]] + tm[["diag"]][[1L]] * (u[[1L]] * u[[1L]])
    for (j in seq_len(q)[-1L]) {
      logp <- logp + tm[["diag"]][[j]] * (u[[j]] * u[[j]])
    }
    for (e in tm[["off"]]) {
      logp <- logp + e[["value"]] * (u[[e[["j"]]]] * u[[e[["k"]]]])
    }
    list(a = agg + logp + zsq, logdet = pr[["logdet"]])
  }

  list(amat = function(pars) weights(pars)[["a"]],
       fn = function(pars) {
         offs <- shift$offs
         w <- weights(pars)
         # RTMB::rowSums, not the base one: base's drops the advector
         # class and the tape then ends in a bare complex vector
         # A prior reads the outer parameters only, so it factors out of
         # every group's integral.
         val <- -sum(offs[, 1L] +
                       log(RTMB::rowSums(exp(w[["a"]] - offs)) / nd)) -
           sum(q / 2 * log(2 * pi) + w[["logdet"]])
         if (is.null(neg_log_prior)) val else val + neg_log_prior(pars)
       },
       neg_log_prior = neg_log_prior,
       set_offset = function(value) {
         shift$offs <- matrix(value, ng, nd)
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
imp_verify <- function(io, nll, plan, template, par, map = NULL) {
  pars <- imp_par_list(template, par, map)
  a <- as.matrix(io$amat(pars))
  nd <- plan[["n_draw"]]
  ng <- plan[["n_group"]]
  refuse <- function(what, ours, theirs) {
    frm_stop("`importance` could not reproduce this model's joint ",
             "log-density from its per-group pieces (", what, ": ",
             format(ours, digits = 10), " against ",
             format(theirs, digits = 10), "). The correction integrates a ",
             "group at a time, so the two have to agree exactly. Report ",
             "this model, and fit it with importance = 0 meanwhile",
             call. = FALSE)
  }
  agrees <- function(ours, theirs, size = 0) {
    is.finite(ours) && is.finite(theirs) &&
      abs(ours - theirs) <= 1e-6 * max(1, abs(theirs), size)
  }
  # the total, at each of two columns: catches an offset shared by
  # every group, which the per-group differences below cannot see
  for (j in unique(c(1L, nd))) {
    joint <- pars
    joint[["b"]] <- plan[["u"]][, j]
    ours <- sum(a[, j] - plan[["zsq"]][, j])
    if (!is.null(io$neg_log_prior)) ours <- ours - io$neg_log_prior(pars)
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
    # the difference of two sums rounds at the sums' size, not its own:
    # measured 5e-6 on an uncentered slope, at a block condition of 4e12
    if (!agrees(mine[g], theirs, abs(nll_base))) {
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
imp_par_list <- function(template, par, map = NULL) {
  out <- template
  pn <- names(par)
  for (cp in unique(pn)) {
    if (is.null(out[[cp]])) next
    mp <- map[[cp]]
    if (is.null(mp)) {
      out[[cp]][] <- unname(par[pn == cp])
    } else {
      # the entries a map fixes keep their template value; recycling the
      # free ones over them would evaluate another model
      fm <- as.integer(factor(mp))
      keep <- !is.na(fm)
      out[[cp]][keep] <- unname(par[pn == cp])[fm[keep]]
    }
  }
  out
}

#' Fit by iterated importance correction: freeze the proposal at the
#' Laplace optimum, optimize theta, refreeze at the new theta, optimize
#' again.
#'
#' Between refreezes the proposal follows the optimizer's theta to first
#' order (`build_importance_objective()`), which is what makes a round's
#' optimum close to the exact one. Refreezing resets the point that first
#' order is taken around, and the estimate is unbiased at every theta
#' whatever the proposal, so a round that moves far is not wrong, only
#' noisy.
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
  # The rounds start where the Laplace route ends. A false collapse of a
  # variance component is a stop the Laplace route leaves through
  # sd_up_escape(), and nothing later leaves it: at log sd -130 every weight
  # is equal. A column's center makes that stop: on an uncentered ordinal, 15
  # of 40 data sets stayed there, 0.15 to 3.33 nats below the Laplace fit.
  lap_esc <- sd_up_escape(
    structure(
      list(spec = frame[["spec"]], frame = frame, obj = lap_obj, opt = lap_opt,
           REML = FALSE, quadrature = FALSE, control = control,
           estimates = solved_par_list(lap_obj, lap_opt$par), prior = NULL,
           cache = new.env(parent = emptyenv())),
      class = "frmtmb_fit"),
    control, bounds, par_units)
  if (!is.null(lap_esc)) lap_opt <- lap_esc
  if (vb) vb_stage("importance warm start", t0, vb_opt_detail(lap_opt))
  lay <- imp_layout(frame[["re_blocks"]])
  gmap <- imp_group_map(frame, lay)
  otpl <- imp_template(template, random)
  omap <- map[intersect(names(map), names(otpl))]
  if (!length(omap)) omap <- NULL
  seed <- control$importance_seed %||% 1L
  rounds <- max(1L, as.integer(control$importance_rounds %||% 5L))

  # One freeze: build the proposal at `at`, center the weights on it,
  # and tape. The offset is the anchor's own row maxima, so the
  # exponentials start at 1 and stay in range as the optimizer walks.
  # It cancels from the value exactly, and being a constant it puts no
  # branch on parameter values the way a running maximum would.
  freeze <- function(at, verify, n_draw) {
    if (vb) tf <- vb_now()
    plan <- imp_plan(lap_obj, frame, lay, at, n_draw, seed)
    io <- build_importance_objective(frame, lay, gmap, plan,
                                     attr(nll, "neg_log_prior"))
    amat <- as.matrix(io$amat(imp_par_list(template, at, map)))
    io$set_offset(apply(amat, 1L, max))
    if (verify) imp_verify(io, nll, plan, template, at, map)
    obj <- RTMB::MakeADFun(io$fn, otpl, map = omap, silent = TRUE)
    if (vb) {
      vb_stage("importance tape", tf,
               paste0(plan[["n_draw"]], " draws over ", plan[["n_group"]],
                      " groups"))
    }
    list(plan = plan, io = io, obj = obj, amat = amat)
  }

  # one complete iteration at one draw count, so that a refusal can
  # measure the count that works by running the same thing
  iterate <- function(n_draw) {
    par <- lap_opt$par
    # The first proposal is the Laplace Gaussian at the Laplace optimum,
    # which is where it is closest to the conditional it approximates.
    fz <- freeze(par, TRUE, n_draw)
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
      # anchor and thrown away for the optimizer.
      opt <- optimize_obj(fz$obj, control, bounds, par_units, verbose = vb,
                          start_par = par)
      moved <- max(abs((opt$par - par) * (par_units %||% 1)))
      moves <- c(moves, moved)
      par <- opt$par
      # Every round ends by refreezing AT its own estimate. That serves
      # twice: it is the next round's proposal, and if the loop stops
      # here it is the tape the fit reports from. Reporting from the
      # anchor is not a detail. A round's optimum is the minimum of the
      # exact objective plus a Monte Carlo error, so the value it stops
      # at is optimistic, and a fresh proposal at that point gives the
      # weights the fit can stand behind. It is also what makes the
      # effective sample sizes below describe the estimate that is
      # actually reported.
      fz <- freeze(par, FALSE, n_draw)
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
    # The iteration must not have made the fit WORSE than the correction
    # it started from. Both values estimate the same exact marginal
    # likelihood, each at its own theta and each from a proposal anchored
    # there, so they compare directly.
    worse <- opt$objective - start_value
    last_grad <- grads[length(grads)]
    capped <- length(moves) >= rounds &&
      moves[length(moves)] >= imp_move_tol &&
      !(is.finite(last_grad) && last_grad < control$grad_tol)
    rose <- is.finite(worse) &&
      worse > imp_worse_tol(ess$mcse, start_ess$mcse)
    list(fz = fz, opt = opt, par = par, moves = moves, grads = grads,
         ess = ess, start_value = start_value, start_ess = start_ess,
         worse = worse, tol = imp_worse_tol(ess$mcse, start_ess$mcse),
         rose = rose, capped = capped, settled = !rose && !capped)
  }

  it <- iterate(n_draw)
  if (it$rose) {
    imp_refuse_draws(it, n_draw, iterate, lay, frame[["n_obs"]],
                     control$importance_ess %||% imp_ess_floor)
  }
  fz <- it$fz
  opt <- it$opt
  moves <- it$moves
  grads <- it$grads
  ess <- it$ess
  start_value <- it$start_value
  start_ess <- it$start_ess
  grad <- grads[length(grads)]
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
       capped = it$capped)
}

#' The refusal for a fit that came out worse than the correction it
#' started from. It runs the same iteration at larger draw counts, so a
#' count it names is one it saw settle: not worse than where it started, not
#' capped, and with every group above the effective sample size floor. It
#' says what it tried and no more. A larger count that it did
#' not try, or that ran out of memory, is not evidence about that count.
#'
#' @noRd
imp_refuse_draws <- function(it, n_draw, iterate, lay, n_obs, ess_floor) {
  top <- min(imp_refuse_max, floor(imp_refuse_cells / n_obs))
  worst <- min(it$ess$ess, it$start_ess$ess)
  ask <- paste0(
    "`importance = ", n_draw, "` did not settle on this model. Optimizing ",
    "over its draws reached estimates where a new estimate of the negative ",
    "log-likelihood is ", format(it$worse, digits = 3), " worse than at ",
    "the Laplace estimates, and the Monte Carlo error explains ",
    format(it$tol, digits = 3), ". The worst group's effective sample size ",
    "is ", format(worst, digits = 2), " of its draws, which is ",
    format(worst * n_draw, digits = 2), " of ", n_draw, ". ")
  tried <- integer(0)
  why_not <- list()
  worst_tried <- NA_real_
  m <- 2L * n_draw
  while (length(tried) < imp_refuse_doublings && m <= top) {
    next_try <- tryCatch(iterate(m), error = function(e) e)
    if (inherits(next_try, "error")) {
      # bad_alloc on a smaller machine: the message is what this is for
      why_not <- list(m = m, msg = conditionMessage(next_try))
      break
    }
    tried <- c(tried, m)
    # A count whose worst group rests on a handful of draws is not one to
    # name, though it passes: on a slope near separation, 128 draws passed at
    # a worst group of 0.012 with an error of 1.87 nats, and 256 to 1024 capped
    if (next_try$settled && min(next_try$ess$ess) >= ess_floor) {
      frm_stop(ask, "The same fit settled at importance = ", m,
               ", with a worst group of ",
               format(min(next_try$ess$ess), digits = 2), " and a Monte Carlo ",
               "standard error of ", format(next_try$ess$mcse, digits = 2),
               ". Use importance = ", m, call. = FALSE)
    }
    worst_tried <- min(next_try$ess$ess, next_try$start_ess$ess)
    m <- 2L * m
  }
  frm_stop(ask,
           if (length(tried)) {
             paste0("It did not settle with every group above an effective ",
                    "sample size of ", ess_floor, " at importance = ",
                    paste(tried, collapse = ", "),
                    " either (worst group ", format(worst_tried, digits = 2),
                    " at the last). ")
           },
           if (length(why_not)) {
             paste0("importance = ", why_not[["m"]], " could not run: ",
                    why_not[["msg"]], ". ")
           } else if (m > imp_refuse_max) {
             paste0("No count above ", m %/% 2L, " was tried: this check ",
                    "stops at ", imp_refuse_max, " draws. ")
           } else if (m > top) {
             paste0("No count above ", m %/% 2L, " was tried: on ", n_obs,
                    " rows it would tape more than ",
                    format(imp_refuse_cells, scientific = TRUE),
                    " rows times draws. ")
           } else {
             paste0("No count above ", m %/% 2L, " was tried. ")
           },
           "Check whether the Laplace fit converged and whether `",
           lay[["label"]], "` has enough rows per group to identify every ",
           "variance component", call. = FALSE)
}

#' How many times that refusal doubles the draw count, the largest count
#' it tries, and the most rows times draws it tapes. A refused fit pays
#' for them, so they stay bounded. Measured on 2400 poisson rows at 8192
#' draws: a tape holds about 1.1 KB per row and draw, and the ladder from
#' 1024 to 8192 took 22.1 GB and 33.8 times the fit that failed, where the
#' fit at 512 took 2.36 GB. Under 2e6 cells the worst ladder, 244 rows from
#' 512 draws, took 3.3 GB and 15.3 times the fit.
#'
#' @noRd
imp_refuse_doublings <- 4L
imp_refuse_max <- 8192L
imp_refuse_cells <- 2e6

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
#' `vignette("diagnostics")` 0.63. Over the 648 fits of six families
#' behind `imp_worse_tol()`, at 64 to 2000 draws, the worst group runs
#' from 0.27 up, and its median is 0.97. Forty groups of three binary
#' rows with a standard deviation of two is genuinely hard, and over
#' seeds 1, 2, 3, 7, 11, 101 and 601 at a thousand draws it reports 0.90,
#' 0.75, 0.41, 0.50, 0.28, 0.83 and 0.85.
#'
#' A proposal that has stopped working sits below it, and draws move it
#' back or do not: a correlated slope in `mu` plus a `sigma` block over
#' 12 groups of 10 rows holds 0.14 at 500 draws and 0.35 at 2000. On 30
#' groups of two null binary rows with a correlated slope, the worst group
#' at 500 draws runs from 0.019 to 1 over ten data sets, and the one at 0.019
#' holds 0.043 at 2000, because its integrand's tail is the prior's and the
#' Gaussian's is narrower. There the estimate rests on a handful of draws
#' whatever the count, and the warning is the point.
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
  frm_warning("The importance proposal covers ", length(bad), " of ",
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
#' Measured over 648 fits (six families, three block shapes and three
#' group sizes at 64 to 2000 draws): the rise runs from -0.365 to +0.078
#' and never passes 0.15 of this tolerance. The floor of 0.1 is a tenth of
#' the log-likelihood unit that one AIC unit is worth, so nothing that
#' could change a model comparison passes it. The Monte Carlo term covers
#' a low draw count, where the two estimates are noisy in their own
#' right.
#'
#' @noRd
imp_worse_tol <- function(mcse_final, mcse_start) {
  max(3 * mcse_final, 3 * mcse_start, 0.1)
}

#' A corrected fit's own proposal, rebuilt from the fit.
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
imp_fit_proposal <- function(fit) {
  fr <- fit[["frame"]]
  lay <- imp_layout(fr[["re_blocks"]])
  tpl <- fit[["estimates"]]
  lap <- RTMB::MakeADFun(build_objective(fr), tpl, random = "b",
                         map = fr[["map"]], silent = TRUE)
  plan <- imp_plan(lap, fr, lay, fit$opt$par, fit$importance$draws,
                   fit$importance$seed)
  list(io = build_importance_objective(fr, lay, imp_group_map(fr, lay), plan),
       plan = plan, template = tpl, map = fr[["map"]])
}

#' The effective sample sizes of that proposal, read at an
#' arbitrary outer parameter vector. One numeric evaluation of the
#' log-weight matrix, no tape.
#'
#' @noRd
imp_ess_at <- function(prop, at) {
  imp_ess(as.matrix(prop$io$amat(imp_par_list(prop$template, at,
                                              prop$map))),
          prop$plan[["n_draw"]])
}

#' Warn when a profile bound was computed where the fit's proposal no
#' longer covers the integrand.
#'
#' `confint(method = "profile")` walks one parameter away from the
#' estimate and reoptimizes the rest ON THE FIT'S OWN TAPE. That tape's
#' proposal follows theta to first order around the estimate, so the
#' further the walk goes the worse the weights can get, while the fit's
#' reported effective sample sizes describe the anchor and say nothing
#' about the path.
#'
#' The anchor is NOT moved along the path, deliberately. Doing so would
#' make the profiled objective a different random function at every grid
#' point, and the interpolation the profile uses to find its crossing
#' assumes one smooth deterministic curve, which is exactly what fixing
#' the draws buys. Measured on a 40-group design, it would also cost at
#' least three times the profile it replaces. So the anchor stays and the
#' bound is reported with a warning when the weights say it should not
#' be trusted.
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
  frm_warning("The profile bound ", format(where, digits = 6), " for '",
              parname, "' lies where this fit's importance proposal has ",
              "stopped covering the integrand: the worst group holds ",
              format(worst, digits = 2), " of its draws there, against a ",
              "threshold of ", floor_at, ". The proposal is anchored at the ",
              "estimate, so a profile that walks far from it is reweighting ",
              "with the wrong Gaussian and this bound is not reliable. Use ",
              "confint(method = \"wald\"), or refit with more draws",
              call. = FALSE)
  invisible(NULL)
}
