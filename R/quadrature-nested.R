# Exact quadrature over two nested random-effect blocks.
#
# WHY NOT TMBad's marginal_gk. One marginal_gk pass over every effect
# integrates only the effects whose integrand splits into
# one-dimensional factors and leaves the rest to the Laplace
# approximation. One pass per block, inner block first, is exact, but
# the transform copies a row's subgraph for every pair of nodes it
# nests: 12 to 551 MB per row (dev/fuzzfix-findings.md, "To file").
# So the iterated integral is written on the tape here instead, as the
# importance correction writes its per-group integral (R/importance.R).
#
# THE RULE. With outer level a (effect u_a) and the inner levels c in
# it (effects v_c), the marginal likelihood is a product over a of
#
#   m_a = int N(u_a) prod_c [int N(v_c) prod_{i in c} p_i dv_c]
#                    prod_{i in a, in no c} p_i du_a.
#
# Both integrals take Gauss-Hermite rules centered and scaled by the
# joint Laplace Gaussian at an anchor: u_a = uhat_a + s_a x_k and, given
# u_a, v_c = vhat_c - r_c (u_a - uhat_a) + t_c y_j, where t_c^-2 = H_cc,
# r_c = H_ca / H_cc and s_a^-2 = H_aa - sum_c H_ca^2 / H_cc are the
# conditional and marginal precisions of that Gaussian. The nodes are
# constants on the tape, so a row's predictor at a node is its fixed
# part plus a constant, and the rule stays a rule for the same
# integral wherever the optimizer walks; it only loses adaptation,
# which re-anchoring restores (quad_nest_fit()).
#
# STACKED, NEVER LOOPED. Every row is repeated once per node pair and
# the family's density is called ONCE on that long vector; the sums over
# a cluster's rows and the log-sum-exp over nodes are sparse indicator
# products and row sums.

#' Gauss-Hermite rule for the standard normal weight (Golub-Welsch).
#'
#' `sum(w * g(x))` approximates `E g(X)`, `X ~ N(0, 1)`, exactly for a
#' polynomial `g` of degree below `2 k`.
#'
#' @noRd
quad_gh_rule <- function(k) {
  if (k == 1L) return(list(x = 0, w = 1))
  i <- seq_len(k - 1L)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  x <- e$values[o]
  w <- e$vectors[1L, o]^2
  # the rule is symmetric; averaging the halves removes the rounding of
  # the eigensolver from the odd moments
  x <- (x - rev(x)) / 2
  w <- (w + rev(w)) / 2
  list(x = x, w = w / sum(w))
}

#' The nodes per level that `frmtmb_control(quad_nodes = )` leaves
#' unset: outer, then inner (dev/nestquad-findings.md, "Nodes").
#'
#' @noRd
quad_nest_nodes_default <- c(11L, 15L)

#' The level of each row in one block, read off the sparsity of every
#' predictor's `Z`, 0 for a row that reads none. A row that reads two
#' levels of one block (a multi-membership term) is refused: its effect
#' is a sum over levels, and the integral then does not factor by level.
#'
#' @noRd
quad_block_levels <- function(frame, bk) {
  n <- frame[["n_obs"]]
  col_level <- rep(seq_len(bk[["n_levels"]]), each = bk[["dim"]])
  rr <- integer(0)
  ll <- integer(0)
  for (lp in frame[["linpreds"]]) {
    Z <- lp[["Z"]]
    if (is.null(Z)) next
    cols <- intersect(bk[["c_idx"]], seq_len(ncol(Z)))
    if (!length(cols)) next
    tz <- methods::as(Matrix::drop0(Z[, cols, drop = FALSE]),
                      "TsparseMatrix")
    rr <- c(rr, tz@i + 1L)
    ll <- c(ll, col_level[match(cols, bk[["c_idx"]])][tz@j + 1L])
  }
  keep <- !duplicated(as.numeric(rr) * (bk[["n_levels"]] + 1) + ll)
  rr <- rr[keep]
  ll <- ll[keep]
  if (anyDuplicated(rr)) {
    frm_stop("quadrature = TRUE needs every row to read at most one level ",
             "of `", block_term_as_written(bk), "`, and row ",
             min(rr[duplicated(rr)]), " reads several. A multi-membership ",
             "effect is a sum over levels, so the integral does not ",
             "factor by level. Fit with quadrature = FALSE, the Laplace ",
             "approximation", call. = FALSE)
  }
  out <- integer(n)
  out[rr] <- ll
  out
}

#' Decide whether two blocks that share rows can be integrated as nested
#' blocks, and lay out the rows and levels if so. Everything this route
#' cannot do is refused here by name, before any tape.
#'
#' Nested means every level of the inner block sits inside one level of
#' the outer block: then, given the outer effect, the inner effects of
#' different levels touch disjoint rows, and the integral over them is a
#' product of one-dimensional integrals. When each block's levels nest in
#' the other's (one grouping factor read by two predictors), the first
#' block is the outer one. Crossed blocks have neither property, and
#' their integral does not factor at all.
#'
#' @noRd
quad_nest_plan <- function(frame, depth) {
  blocks <- frame[["re_blocks"]]
  labs <- vapply(blocks, block_term_as_written, "")
  quote_labs <- function(i) paste0("`", labs[i], "`", collapse = ", ")
  spec <- frame[["spec"]]
  if (length(spec$responses) != 1L) {
    frm_stop("quadrature = TRUE integrates nested blocks for one response, ",
             "and this model has ", length(spec$responses), " (",
             paste(names(spec$responses), collapse = ", "), "). Fit with ",
             "quadrature = FALSE, the Laplace approximation", call. = FALSE)
  }
  # the predictors each block enters, read off Z
  reads <- lapply(blocks, function(bk) {
    sort(names(Filter(function(lp) {
      Z <- lp[["Z"]]
      if (is.null(Z)) return(FALSE)
      cols <- intersect(bk[["c_idx"]], seq_len(ncol(Z)))
      length(cols) > 0L && any(Z[, cols, drop = FALSE] != 0)
    }, frame[["linpreds"]])))
  })
  two_preds <- length(unique(reads)) > 1L
  if (length(blocks) != 2L) {
    by_split <- Filter(function(bk) !is.null(bk[["by"]]), blocks)
    why <- if (length(by_split)) {
      paste0("`", block_term_as_written(by_split[[1L]]), "` is a by-split ",
             "term, one block per level of its by-variable, and the rule ",
             "takes one outer and one inner block")
    } else if (two_preds) {
      paste0("they enter different predictors (",
             paste(unique(unlist(reads)), collapse = ", "), "), and the ",
             "rule takes one outer and one inner block on the same ",
             "predictors")
    } else if (depth > 2L) {
      paste0("a row reads ", depth, " of them, and each further nested ",
             "level multiplies the nodes every row is evaluated at")
    } else {
      "the rule takes one outer and one inner block"
    }
    frm_stop("quadrature = TRUE integrates two random-effect blocks on a ",
             "row, one nested in the other, and this model has ",
             length(blocks), " blocks (", quote_labs(seq_along(blocks)),
             "): ", why, ". Fit with quadrature = FALSE, the Laplace ",
             "approximation", call. = FALSE)
  }
  lev <- lapply(blocks, quad_block_levels, frame = frame)
  # inner levels per outer level, and the reverse, over the rows that
  # read both
  both <- lev[[1L]] > 0L & lev[[2L]] > 0L
  pr <- unique(cbind(lev[[1L]][both], lev[[2L]][both]))
  in_one <- function(inner) {
    all(tabulate(pr[, inner], blocks[[inner]][["n_levels"]]) <= 1L)
  }
  inner <- if (in_one(2L)) 2L else if (in_one(1L)) 1L else 0L
  if (inner == 0L) {
    # crossing needs both directions: a level of each block that reaches
    # several levels of the other, so neither nests in the other
    lv_count <- function(k) tabulate(pr[, k], blocks[[k]][["n_levels"]])
    c1 <- lv_count(1L)
    c2 <- lv_count(2L)
    ex1 <- which.max(c1)
    ex2 <- which.max(c2)
    g1 <- blocks[[1L]][["group_name"]]
    g2 <- blocks[[2L]][["group_name"]]
    frm_stop("quadrature = TRUE integrates two blocks that share rows only ",
             "when one is nested in the other, and ", quote_labs(1:2),
             " are crossed: level '", blocks[[2L]][["levels"]][ex2], "' of ",
             g2, " reaches ", c2[ex2], " levels of ", g1, ", so ", g2,
             " does not nest in ", g1, ", and level '",
             blocks[[1L]][["levels"]][ex1], "' of ", g1, " reaches ",
             c1[ex1], " levels of ", g2, ", so ", g1, " does not nest in ",
             g2, ". The integral over crossed effects does not factor into ",
             "one-dimensional integrals, one level inside another. Fit ",
             "with quadrature = FALSE, the Laplace approximation",
             call. = FALSE)
  }
  # A block on another predictor than its partner's, such as (1 | g) in
  # mu and in sigma, makes the width of one integrand change with the
  # other effect. The joint Gaussian gives the inner rule one width for
  # every outer node, and the rule then converged slowly: on 72 gaussian
  # rows the log-likelihood error at 11 x 15 nodes was 7e-3, against
  # 7e-6 for (1 | ga/gb) on Bernoulli rows (dev/nestquad-findings.md,
  # "Blocks on two predictors").
  if (two_preds) {
    frm_stop("quadrature = TRUE integrates two nested blocks when they ",
             "enter the same predictors, and ", quote_labs(1:2), " enter ",
             "different ones. An effect on another predictor, such as a ",
             "scale, changes the width of the other effect's integrand ",
             "from node to node, and the rule then converges slowly. Fit ",
             "with quadrature = FALSE, the Laplace approximation",
             call. = FALSE)
  }
  outer <- 3L - inner
  bo <- blocks[[outer]]
  bi <- blocks[[inner]]
  for (bk in blocks) {
    if (!bk[["covstruct"]] %in% c("us", "diag", "homdiag")) {
      frm_stop("quadrature = TRUE integrates nested blocks with Gaussian ",
               "effects, and `", block_term_as_written(bk), "` is a ",
               "Student-t block. Its effects are integrated exactly by ",
               "quadrature only as the one block of a row. Fit with ",
               "quadrature = FALSE, the Laplace approximation", call. = FALSE)
    }
  }
  if (!is.null(frame[["map"]][["b"]])) {
    frm_stop("quadrature = TRUE cannot integrate nested blocks whose ",
             "random-effect coefficients are mapped: the rule is built per ",
             "level, and a map merges or fixes coefficients across levels. ",
             "Fit with quadrature = FALSE, the Laplace approximation",
             call. = FALSE)
  }
  rn <- names(spec$responses)[[1L]]
  fam <- spec$responses[[rn]]$family
  if (!is.null(fam_structure(fam)[["loglik"]])) {
    frm_stop("quadrature = TRUE integrates nested blocks against per-row ",
             "densities, and the '", fam[["family"]], "' family supplies ",
             "one log-likelihood for the whole response. Fit with ",
             "quadrature = FALSE, the Laplace approximation", call. = FALSE,
             package = frm_family_package(fam))
  }
  n <- frame[["n_obs"]]
  for (lp in frame[["linpreds"]]) {
    what <- paste0(lp[["resp"]], ".", lp[["dpar"]])
    if (!is.null(lp[["nl_body"]])) {
      frm_stop("quadrature = TRUE cannot integrate nested blocks under a ",
               "nonlinear predictor (", what, "): the rule evaluates every ",
               "row once per node pair, and a body that reads a raw data ",
               "column would recycle it against the longer vector. Fit ",
               "with quadrature = FALSE, the Laplace approximation",
               call. = FALSE)
    }
    if (length(lp[["cs"]] %||% list())) {
      frm_stop("quadrature = TRUE cannot integrate nested blocks with a ",
               "cs() term (", what, "): its threshold offsets are a matrix ",
               "per row, which the rule does not repeat per node. Fit with ",
               "quadrature = FALSE, the Laplace approximation", call. = FALSE)
    }
    if (NROW(lp[["X"]]) != n) {
      frm_stop("quadrature = TRUE needs every predictor of nested blocks to ",
               "read all ", n, " rows, and ", what, " reads ",
               NROW(lp[["X"]]), ". Fit with quadrature = FALSE, the Laplace ",
               "approximation", call. = FALSE)
    }
  }
  lo <- lev[[outer]]
  li <- lev[[inner]]
  # the outer level of each inner level; an inner level whose rows read
  # no outer level has no outer effect in its integrand, so any outer
  # level holds it exactly: its integral is the same at every outer node
  # and factors out of that level's log-sum-exp
  inner_outer <- rep(1L, bi[["n_levels"]])
  w <- li > 0L & lo > 0L
  inner_outer[li[w]] <- lo[w]
  row_outer <- ifelse(li > 0L, inner_outer[pmax(li, 1L)], pmax(lo, 1L))
  list(outer = bo, inner = bi, row_outer = row_outer, row_inner = li,
       inner_outer = inner_outer, labels = labs[c(outer, inner)])
}

#' The joint Laplace Gaussian at the outer parameters `par`, in the
#' coordinates the rule needs: modes, the inner levels' conditional
#' scale and regression on their outer effect, and the outer levels'
#' marginal scale.
#'
#' The conditional Hessian must have no entry outside the pairs the
#' nesting allows (a level with itself, an inner level with its outer
#' level). That is what makes the integral factor, so it is measured
#' here rather than assumed.
#'
#' @noRd
quad_nest_anchor <- function(lap_obj, plan, par) {
  lap_obj$fn(par)
  full <- lap_obj$env$last.par
  b <- as.numeric(full[lap_obj$env$random])
  hs <- lap_obj$env$spHess(full, random = TRUE)
  h <- methods::as(Matrix::forceSymmetric(hs, uplo = "U"), "TsparseMatrix")
  ia <- plan$outer[["b_idx"]]
  ic <- plan$inner[["b_idx"]]
  na <- length(ia)
  nc <- length(ic)
  role <- integer(length(b))
  lvl <- integer(length(b))
  role[ia] <- 1L
  lvl[ia] <- seq_len(na)
  role[ic] <- 2L
  lvl[ic] <- seq_len(nc)
  i <- h@i + 1L
  j <- h@j + 1L
  x <- h@x
  dg <- i == j
  haa <- numeric(na)
  hcc <- numeric(nc)
  haa[lvl[i[dg & role[i] == 1L]]] <- x[dg & role[i] == 1L]
  hcc[lvl[i[dg & role[i] == 2L]]] <- x[dg & role[i] == 2L]
  off <- !dg
  # an off-diagonal entry is allowed only between an inner level and
  # its own outer level
  ci <- ifelse(role[i] == 2L, i, j)
  ai <- ifelse(role[i] == 2L, j, i)
  ok <- off & role[i] != role[j] & role[ci] == 2L &
    plan$inner_outer[lvl[ci]] == lvl[ai]
  hca <- numeric(nc)
  hca[lvl[ci[ok]]] <- x[ok]
  bad <- off & !ok & x != 0
  scale <- max(abs(x[dg]), 1)
  if (any(bad) && max(abs(x[bad])) > 1e-8 * scale) {
    frm_stop("quadrature = TRUE needs the effects of ", plan$labels[[1L]],
             " and ", plan$labels[[2L]], " to be independent across ",
             "levels given the data, and the conditional Hessian links two ",
             "levels the nesting keeps apart (largest entry ",
             format(max(abs(x[bad])), digits = 3), " against a diagonal of ",
             format(scale, digits = 3), "). The integral then does not ",
             "factor. Fit with quadrature = FALSE, the Laplace approximation",
             call. = FALSE)
  }
  # the outer level's marginal precision: its own, less what each of its
  # inner levels explains
  pa <- haa
  shed <- rowsum(hca^2 / hcc, plan$inner_outer)
  at <- as.integer(rownames(shed))
  pa[at] <- pa[at] - shed[, 1L]
  if (any(!is.finite(hcc) | hcc <= 0) || any(!is.finite(pa) | pa <= 0)) {
    frm_stop("quadrature = TRUE found no Gaussian to center the rule on: ",
             "the Laplace fit leaves a level of ", plan$labels[[1L]], " or ",
             plan$labels[[2L]], " with a non-positive conditional ",
             "precision. The fit sits on a singular variance component; ",
             "fit with quadrature = FALSE, the Laplace approximation",
             call. = FALSE)
  }
  list(par = par, full = full, uhat = b[ia], vhat = b[ic],
       s = 1 / sqrt(pa), t = 1 / sqrt(hcc), r = hca / hcc)
}

#' The nested-quadrature negative log-likelihood closure on the outer
#' parameters, with the nodes and every row-to-level map baked in.
#'
#' The nodes are frozen on the standardized scale, z = u / sd: the tape
#' reads u = sd(theta) z, and the prior is the weight of the rule, so
#' theta reaches the integrand through the likelihood alone. Frozen on
#' the effects' own scale instead, the rule kept the prior density at
#' fixed points, which grows like 1 / sd as an sd shrinks toward a node
#' near 0: the optimizer found that, and a Bernoulli fit reported a
#' log-likelihood of +29.6 (dev/nestquad-findings.md, "Fix round 1").
#' At the anchor the two forms are the same rule.
#'
#' `nodes` is `c(outer, inner)`. Returns the closure `fn`, `prime()`,
#' which fixes the log-sum-exp offsets at a numeric parameter list (they
#' cancel exactly, and as constants they put no branch on the tape), and
#' `joint()`, the joint log-density at one node pair, which
#' `quad_nest_verify()` compares with the model's own objective.
#'
#' @noRd
quad_nest_objective <- function(frame, plan, anc, nodes, pars_a,
                                neg_log_prior = NULL) {
  lps <- frame[["linpreds"]]
  spec <- frame[["spec"]]
  rn <- names(spec$responses)[[1L]]
  fam <- spec$responses[[rn]]$family
  n <- frame[["n_obs"]]
  ko <- as.integer(nodes[[1L]])
  ki <- as.integer(nodes[[2L]])
  gx <- quad_gh_rule(ko)
  gy <- quad_gh_rule(ki)
  # the log of each node's weight against the integrand itself:
  # int e^h dz = tau int e^(h + x^2 / 2) sqrt(2 pi) phi(x) dx
  lwx <- log(gx$w) + gx$x^2 / 2 + 0.5 * log(2 * pi)
  lwy <- log(gy$w) + gy$x^2 / 2 + 0.5 * log(2 * pi)
  bo <- plan$outer
  bi <- plan$inner
  na <- bo[["n_levels"]]
  nc <- bi[["n_levels"]]
  li <- plan$row_inner
  lo <- plan$row_outer
  io <- plan$inner_outer
  r1 <- which(li > 0L)
  r0 <- which(li == 0L)
  n1 <- length(r1)
  n0 <- length(r0)

  # Gaussian per-level priors, read from the registry's own density; the
  # coefficient of u^2 is -1 / (2 sd^2), so it gives the sd on the tape
  theta_a <- pars_a[["theta"]]
  pri_o <- imp_block_prior(bo, imp_prior_scales(list(blocks = list(bo)),
                                                theta_a)[[1L]])
  pri_i <- imp_block_prior(bi, imp_prior_scales(list(blocks = list(bi)),
                                                theta_a)[[1L]])
  sd_of <- function(pr) sqrt(-0.5 / pr[["diag"]][[1L]])
  sds <- function(pars) {
    th <- pars[["theta"]]
    list(o = sd_of(pri_o(th[bo[["theta_idx"]]])),
      i = sd_of(pri_i(th[bi[["theta_idx"]]])))
  }
  sd0 <- sds(pars_a)

  # node values on the effects' scale at the anchor, then standardized:
  # u[a, k], and v[c, k, j] given the outer node k (the joint Gaussian's
  # conditional mean and scale). Centering each (c, k) on its own
  # conditional mode instead, found by Newton steps, was measured no
  # better on Poisson rows and worse on Bernoulli rows
  # (dev/nestquad-findings.md, "Nodes").
  U <- anc$uhat + outer(anc$s, gx$x)
  vck <- anc$vhat - outer(anc$r * anc$s[io], gx$x)
  V <- array(rep(as.vector(vck), ki), c(nc, ko, ki)) +
    array(rep(anc$t, ko * ki) * rep(gy$x, each = nc * ko), c(nc, ko, ki))
  zU <- U / sd0[["o"]]
  zV <- V / sd0[["i"]]

  # rows in stacked order: row fastest, then the outer node. The inner
  # node is the loop in rows() below, one block per node: the short
  # dimension, so each block's temporaries are garbage before the next
  # is built (dev/nestquad-findings.md, "Memory")
  idx1 <- rep.int(r1, ko)
  k1 <- rep(seq_len(ko), each = n1)
  idx0 <- rep.int(r0, ko)
  zcol <- function(Z, bk) {
    cols <- intersect(bk[["c_idx"]], seq_len(ncol(Z)))
    if (!length(cols)) return(NULL)
    zv <- as.vector(Z[, cols, drop = FALSE] %*% rep(1, length(cols)))
    if (all(zv == 0)) NULL else zv
  }
  # what the standardized nodes add to each predictor per unit sd, a
  # constant per stacked row: the outer one shared by every inner node,
  # the inner one a column per inner node
  add_o1 <- add_i1 <- add_o0 <- list()
  zu1 <- zU[cbind(lo[idx1], k1)]
  zu0 <- zU[cbind(lo[idx0], rep(seq_len(ko), each = n0))]
  for (key in names(lps)) {
    Z <- lps[[key]][["Z"]]
    if (is.null(Z)) next
    za <- zcol(Z, bo)
    zi <- zcol(Z, bi)
    if (!is.null(za)) {
      add_o1[[key]] <- za[idx1] * zu1
      if (n0) add_o0[[key]] <- za[idx0] * zu0
    }
    if (!is.null(zi)) {
      a1 <- matrix(0, n1 * ko, ki)
      for (j in seq_len(ki)) {
        a1[, j] <- zi[idx1] * zV[cbind(li[idx1], k1, j)]
      }
      add_i1[[key]] <- a1
    }
  }
  rm(zu1, zu0, k1)
  yv <- frame[["y"]][[rn]]
  av <- frame[["aterm_values"]][[rn]]
  wts <- av[["weights"]] %||% rep(1, n)
  if (length(wts) == 1L) wts <- rep(wts, n)
  gone <- frame[["mi_map"]][[rn]][["unread"]]
  if (length(gone)) wts[gone] <- 0
  y1 <- imp_expand(yv, idx1, n)
  av1 <- imp_expand(av, idx1, n)
  w1 <- wts[idx1]
  y0 <- imp_expand(yv, idx0, n)
  av0 <- imp_expand(av, idx0, n)
  w0 <- wts[idx0]
  extra_names <- frame[["extra_names"]] %||% character(0)

  m_in <- Matrix::sparseMatrix(i = li[r1], j = seq_len(n1), x = rep(1, n1),
                               dims = c(nc, n1))
  m_out <- Matrix::sparseMatrix(i = io, j = seq_len(nc), x = rep(1, nc),
                                dims = c(na, nc))
  m_0 <- Matrix::sparseMatrix(i = lo[r0], j = seq_len(n0), x = rep(1, n0),
                              dims = c(na, n0))
  extra_of <- function(pars) {
    extra <- NULL
    if (length(extra_names)) {
      extra <- lapply(stats::setNames(extra_names, extra_names),
                      function(nm) pars[[nm]])
    }
    resp_extras(frame, extra, rn)
  }

  # every weight of the rule is a constant: the node weight, the
  # standard normal density at the standardized node, and the rule's
  # scale on the standardized axis
  cy <- rep(lwy, each = nc * ko) + stats::dnorm(as.vector(zV), log = TRUE)
  cx <- rep(lwx, each = na) + stats::dnorm(as.vector(zU), log = TRUE)
  logt <- rep(log(anc$t / sd0[["i"]]), ko)
  logs <- log(anc$s / sd0[["o"]])
  st <- new.env(parent = emptyenv())
  st$off_i <- rep(0, nc * ko * ki)
  st$off_o <- rep(0, na * ko)

  # densities of every stacked row: the row log-likelihood sums per
  # inner level and node pair (nc x ko ki), and the inner-free rows' per
  # outer level and node (na x ko)
  rows <- function(pars) {
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    extra <- extra_of(pars)
    sd <- sds(pars)
    eta1 <- list()
    dp0 <- list()
    for (key in names(lps)) {
      lp <- lps[[key]]
      eta <- lp_eta_fixed(lp, pars, n, list(), frame[["y"]])
      e1 <- eta[idx1]
      if (!is.null(add_o1[[key]])) e1 <- e1 + sd[["o"]] * add_o1[[key]]
      eta1[[key]] <- e1
      if (n0) {
        e0 <- eta[idx0]
        if (!is.null(add_o0[[key]])) e0 <- e0 + sd[["o"]] * add_o0[[key]]
        dp0[[lp[["dpar"]]]] <- lp[["link"]]$linkinv(e0)
        dp0[[paste0(".eta_", lp[["dpar"]])]] <- e0
      }
    }
    # an advector matrix on the tape, a Matrix object in the numeric
    # path; the class test resolves while the tape is built
    dense <- function(m) {
      as.vector(if (inherits(m, "advector")) m else as.matrix(m))
    }
    s1 <- vector("list", ki)
    for (j in seq_len(ki)) {
      dp1 <- list()
      for (key in names(lps)) {
        lp <- lps[[key]]
        e1 <- eta1[[key]]
        if (!is.null(add_i1[[key]])) {
          e1 <- e1 + sd[["i"]] * add_i1[[key]][, j]
        }
        dp1[[lp[["dpar"]]]] <- lp[["link"]]$linkinv(e1)
        dp1[[paste0(".eta_", lp[["dpar"]])]] <- e1
      }
      s1[[j]] <- dense(m_in %*% RTMB::matrix(
        w1 * row_lpdf(fam, y1, y1, dp1, av1, extra), n1, ko))
    }
    s1 <- do.call(c, s1)
    s0 <- if (n0) {
      dense(m_0 %*% RTMB::matrix(
        w0 * row_lpdf(fam, y0, y0, dp0, av0, extra), n0, ko))
    }
    list(s1 = s1, s0 = s0)
  }
  # log-sum-exp over the columns of an m x k layout held as a vector;
  # `off` is that layout's constant offset
  lse_rows <- function(v, off, m, k) {
    log(RTMB::rowSums(exp(RTMB::matrix(v - off, m, k)))) +
      off[seq_len(m)]
  }
  outer_of <- function(rw, off_i) {
    log_i <- lse_rows(rw$s1 + cy, off_i, nc * ko, ki) + logt
    tt <- m_out %*% RTMB::matrix(log_i, nc, ko)
    if (!inherits(tt, "advector")) tt <- as.matrix(tt)
    tt <- as.vector(tt)
    if (n0) tt <- tt + rw$s0
    tt + cx
  }
  fn <- function(pars) {
    rw <- rows(pars)
    tt <- outer_of(rw, st$off_i)
    val <- -sum(lse_rows(tt, st$off_o, na, ko) + logs)
    if (is.null(neg_log_prior)) val else val + neg_log_prior(pars)
  }
  prime <- function(pars) {
    rw <- rows(pars)
    m <- matrix(rw$s1 + cy, nc * ko, ki)
    st$off_i <- rep(apply(m, 1L, max), ki)
    tt <- matrix(outer_of(rw, st$off_i), na, ko)
    st$off_o <- rep(apply(tt, 1L, max), ko)
    invisible(NULL)
  }
  joint <- function(pars, k, j) {
    rw <- rows(pars)
    sd <- sds(pars)
    s1 <- matrix(rw$s1, nc * ko, ki)[(k - 1L) * nc + seq_len(nc), j]
    vv <- sd[["i"]] * zV[, k, j]
    uu <- sd[["o"]] * zU[, k]
    th <- pars[["theta"]]
    po <- pri_o(th[bo[["theta_idx"]]])
    pi_ <- pri_i(th[bi[["theta_idx"]]])
    val <- sum(s1) +
      sum(pi_[["c0"]] + pi_[["diag"]][[1L]] * vv^2) +
      sum(po[["c0"]] + po[["diag"]][[1L]] * uu^2)
    if (n0) val <- val + sum(matrix(rw$s0, na, ko)[, k])
    b <- numeric(length(anc$full[names(anc$full) == "b"]))
    b[bo[["b_idx"]]] <- uu
    b[bi[["b_idx"]]] <- vv
    list(value = val, b = b)
  }
  list(fn = fn, prime = prime, joint = joint, nodes = c(ko, ki),
       n_stack = n1 * ko * ki + n0 * ko)
}

#' The pin that keeps the stacked reimplementation honest: the joint
#' log-density at three node pairs, against the model's own objective
#' with the random effects set to those nodes. The pairs are the centre
#' and two opposite corners, so a block read on the wrong rows or with
#' the wrong design value shows at one of them.
#'
#' @noRd
quad_nest_verify <- function(qo, nll, template, par, map = NULL) {
  pars <- imp_par_list(template, par, map)
  ko <- qo$nodes[[1L]]
  ki <- qo$nodes[[2L]]
  nlp <- attr(nll, "neg_log_prior")
  at <- unique(list(c((ko + 1L) %/% 2L, (ki + 1L) %/% 2L), c(1L, ki),
                    c(ko, 1L)))
  for (kj in at) {
    jt <- qo$joint(pars, kj[[1L]], kj[[2L]])
    full <- pars
    full[["b"]] <- jt$b
    ours <- jt$value
    if (!is.null(nlp)) ours <- ours - nlp(pars)
    theirs <- -nll(full)
    if (!(is.finite(ours) && is.finite(theirs) &&
            abs(ours - theirs) <= 1e-6 * max(1, abs(theirs)))) {
      frm_stop("quadrature = TRUE could not reproduce this model's joint ",
               "log-density from its nested pieces (node pair ",
               kj[[1L]], ", ", kj[[2L]], ": ", format(ours, digits = 10),
               " against ", format(theirs, digits = 10), "). Report this ",
               "model, and fit it with quadrature = FALSE meanwhile",
               call. = FALSE)
    }
  }
  invisible(NULL)
}

#' The size of a taped objective: operations, values and input indices.
#'
#' TMB exports no accessor for this; `TMB:::info()` reads it from the
#' tape, and RTMB's tapes carry the same entry point. NULL when it cannot
#' be read, and the guard then falls back to its count of stacked rows.
#'
#' @noRd
quad_tape_size <- function(obj) {
  info <- tryCatch(utils::getFromNamespace("info", "TMB")(obj$env$ADFun),
                   error = function(e) NULL)
  if (is.null(info)) return(NULL)
  v <- c(ops = info[["opstack_size"]], values = info[["values_size"]],
         inputs = info[["inputs_size"]])
  if (length(v) != 3L || anyNA(v)) return(NULL)
  as.numeric(v)
}

#' Predict the memory of the nested tape before building it, and refuse
#' above `frmtmb_control(quad_max_gb = )`.
#'
#' The tape is linear in the node counts: every stacked row and every
#' (inner level, node pair) cell costs the same operations at any node
#' count. So three small probe tapes, at nodes (1, 1), (1, 2) and (2, 1),
#' give the cost of one more inner node layer and of one more outer node,
#' and the full tape is extrapolated from them. The sizes do not depend
#' on the values the nodes take, so the probes run on a placeholder
#' anchor, before the Laplace warm start. The peak the fit adds is not
#' the tape: taping holds several times its bytes for a moment, by a
#' factor that depends on the density's operations and is measured per
#' family (`quad_nest_mem_cost`).
#'
#' @noRd
quad_nest_guard <- function(frame, plan, nodes, template, random, max_gb,
                            omap, neg_log_prior = NULL) {
  otpl <- imp_template(template, random)
  anc <- quad_nest_probe_anchor(plan, template)
  probe <- function(k) {
    qo <- quad_nest_objective(frame, plan, anc, k, template, neg_log_prior)
    obj <- RTMB::MakeADFun(qo$fn, otpl, map = omap, silent = TRUE)
    quad_tape_size(obj)
  }
  ko <- nodes[[1L]]
  ki <- nodes[[2L]]
  key <- quad_nest_family_key(frame)
  s11 <- probe(c(1L, 1L))
  tape <- NULL
  if (!is.null(s11)) {
    s12 <- probe(c(1L, 2L))
    s21 <- probe(c(2L, 1L))
    inner <- s12 - s11
    outer <- s21 - s11 - inner
    tape <- sum((s11 + inner * (ko * ki - 1) + outer * (ko - 1)) *
                  quad_nest_entry_bytes)
  }
  bytes <- quad_nest_bytes(tape, plan, nodes, key)
  out <- list(bytes = bytes, tape = tape, key = key)
  if (bytes <= max_gb * 2^30) return(invisible(out))
  frm_stop("quadrature = TRUE over ", plan$labels[[1L]], " and ",
           plan$labels[[2L]], " at ", ko, " x ", ki, " nodes would need ",
           "about ", format(bytes / 2^30, digits = 3), " GB of memory, ",
           "above the limit of ", format(max_gb), " GB ",
           "(frmtmb_control(quad_max_gb = )). Its tape repeats every row ",
           "once per node pair, so the memory grows with the rows times ",
           ko, " x ", ki, ". Use fewer nodes (frmtmb_control(quad_nodes = ",
           ")), raise quad_max_gb if the machine has the memory, or fit ",
           "with quadrature = FALSE, the Laplace approximation",
           call. = FALSE)
}

#' A placeholder anchor for the guard's probe tapes: the tape's size does
#' not read the node values, only their count. The values are nonzero,
#' because RTMB folds an operation with a constant zero operand, which
#' would shrink a probe below the tape it stands for.
#'
#' @noRd
quad_nest_probe_anchor <- function(plan, template) {
  na <- plan$outer[["n_levels"]]
  nc <- plan$inner[["n_levels"]]
  nb <- length(template[["b"]])
  list(par = NULL, full = stats::setNames(numeric(nb), rep("b", nb)),
       uhat = rep(0.5, na), vhat = rep(0.5, nc), s = rep(1, na),
       t = rep(1, nc), r = rep(0.25, nc))
}

#' The key the measured memory factors are read under: the family, and
#' whether any row is censored, which changes the density's operations.
#'
#' @noRd
quad_nest_family_key <- function(frame) {
  rn <- names(frame[["spec"]]$responses)[[1L]]
  fam <- frame[["spec"]]$responses[[rn]]$family[["family"]] %||% ""
  cen <- frame[["aterm_values"]][[rn]][["cens"]]
  if (!is.null(cen) && any(cen != 0)) paste0(fam, "+cens") else fam
}

#' The memory a nested fit adds at its peak: the larger of the predicted
#' tape's bytes times the family's measured factor and the stacked rows
#' times the family's measured cost per stacked row, plus a fixed part.
#' Both terms are needed. Taping holds several times the final tape
#' beside it, and the factor depends on the density's operations, which
#' the tape's size measures within a family. A censored Beta row is taped
#' through `pbeta()` and then mostly eliminated, so its final tape is
#' small while taping it took 540 kB per stacked row whatever the share
#' of censored rows: only the count of stacked rows predicts it.
#' Without the tape's size only the second term is used.
#'
#' @noRd
quad_nest_bytes <- function(tape, plan, nodes, key) {
  li <- plan$row_inner
  stacked <- sum(li > 0L) * nodes[[1L]] * nodes[[2L]] +
    sum(li == 0L) * nodes[[1L]]
  cost <- quad_nest_mem_cost[[key]]
  cens <- grepl("+cens", key, fixed = TRUE)
  known <- !is.null(cost)
  if (!known) {
    # a family not measured: the largest factor measured, and for a
    # censored one the largest cost per row measured
    cost <- c(factor = max(vapply(quad_nest_mem_cost, `[[`, 1, "factor")),
              row_kb = if (cens) quad_nest_mem_cens_kb else 0)
  }
  per_row <- stacked * cost[["row_kb"]] * 1024
  if (is.null(tape)) {
    if (!known && !cens) per_row <- stacked * quad_nest_mem_plain_kb * 1024
    return(quad_nest_peak_fixed + per_row)
  }
  quad_nest_peak_fixed + max(cost[["factor"]] * tape, per_row)
}

#' Bytes per operation, value and input index of a taped objective; the
#' fixed part of a nested fit's peak; and, per family, the factor from
#' the tape's bytes to the peak and the peak per stacked row in kB.
#' Each is the largest measured on 3000 rows (20 to 1000 for the costly
#' ones) at 11 x 15 nodes, times 1.3 (dev/nestquad-findings.md, "Memory,
#' fix round 1"). A family not listed gets the largest factor, and when
#' censored the largest censored cost per row; without the tape's size,
#' an unlisted uncensored family gets the largest uncensored cost per
#' row.
#'
#' @noRd
quad_nest_entry_bytes <- c(ops = 8, values = 16, inputs = 4)
quad_nest_peak_fixed <- 128 * 2^20
quad_nest_mem_cost <- list(
  bernoulli = c(factor = 9.5, row_kb = 1.1),
  binomial = c(factor = 8.5, row_kb = 1.2),
  lognormal = c(factor = 8.6, row_kb = 2.3),
  von_mises = c(factor = 10.7, row_kb = 2.8),
  beta_binomial = c(factor = 9.6, row_kb = 11.4),
  poisson = c(factor = 8.3, row_kb = 1.7),
  gaussian = c(factor = 9.3, row_kb = 2.2),
  cumulative = c(factor = 9.3, row_kb = 2.6),
  Gamma = c(factor = 10.3, row_kb = 2.9),
  beta = c(factor = 10.9, row_kb = 7.8),
  negbinomial = c(factor = 9.9, row_kb = 5.1),
  hurdle_poisson = c(factor = 8.5, row_kb = 2.4),
  zero_inflated_negbinomial = c(factor = 10.3, row_kb = 4.6),
  skew_normal = c(factor = 11.1, row_kb = 4.5),
  weibull = c(factor = 11.3, row_kb = 4.4),
  student = c(factor = 42.9, row_kb = 14.0),
  exgaussian = c(factor = 10.8, row_kb = 67.6),
  `lognormal+cens` = c(factor = 10.3, row_kb = 2.6),
  `poisson+cens` = c(factor = 20.4, row_kb = 28.6),
  `beta+cens` = c(factor = 30.0, row_kb = 711)
)
quad_nest_mem_cens_kb <- 711
quad_nest_mem_plain_kb <- 67.6

#' Fit by nested quadrature: anchor the rule at the Laplace optimum,
#' optimize the frozen tape inside a trust region, re-anchor at that
#' optimum, and repeat until the freshly anchored objective is
#' stationary.
#'
#' The frozen tape is a rule for the integral only near its anchor, and
#' its optimum can lie far outside that region: on a near-separable
#' Bernoulli the first round stepped by 7114 in parameter units. So each
#' round optimizes inside a box of the trust radius around the anchor,
#' and is accepted only when the objective re-anchored at its end is no
#' worse than the one re-anchored at its start. Both values are the rule
#' anchored at its own point, the best value of the integral there, so
#' the accepted values fall and the loop cannot cycle; the unguarded loop
#' moved by 2 to 17 every round for 20 rounds on that design
#' (dev/nestquad-findings.md, "Fix round 1"). An accepted round that
#' reaches the edge of the box doubles the radius; a rejected one cuts it
#' to a quarter of its step, and three rejections in a row end the loop:
#' the rule's own error then hides any smaller improvement. A loop that
#' ends without a stationary anchor warns with class
#' `frmtmb_quad_unsettled`.
#'
#' @noRd
quad_nest_fit <- function(nll, template, random, map, lap_obj, control,
                          bounds, par_units, frame, plan, vb = 0L,
                          rounds = 25L) {
  nodes <- control$quad_nodes %||% quad_nest_nodes_default
  if (length(nodes) == 1L) nodes <- rep(nodes, 2L)
  otpl <- imp_template(template, random)
  omap <- map[intersect(names(map), names(otpl))]
  if (!length(omap)) omap <- NULL
  nlp <- attr(nll, "neg_log_prior")
  if (vb) t0 <- vb_now()
  mem <- quad_nest_guard(frame, plan, nodes, template, random,
                         control$quad_max_gb %||% 4, omap, nlp)
  if (vb) {
    vb_stage("nested quadrature memory", t0,
             paste0("about ", format(mem$bytes / 2^20, digits = 3),
                    " MB predicted (tape ",
                    format((mem$tape %||% NA) / 2^20, digits = 3),
                    " MB, family ", mem$key, ")"))
    t0 <- vb_now()
  }
  lap_opt <- optimize_obj(lap_obj, control, bounds, par_units)
  if (vb) vb_stage("quadrature warm start", t0, vb_opt_detail(lap_opt))
  freeze <- function(at, first) {
    if (vb) tf <- vb_now()
    anc <- quad_nest_anchor(lap_obj, plan, at)
    pl <- imp_par_list(template, at, map)
    qo <- quad_nest_objective(frame, plan, anc, nodes, pl, nlp)
    qo$prime(pl)
    if (first) quad_nest_verify(qo, nll, template, at, map)
    # the numeric passes above leave stacked-length garbage behind; the
    # tape below is the peak, so it starts from a collected heap
    gc(verbose = FALSE)
    obj <- RTMB::MakeADFun(qo$fn, otpl, map = omap, silent = TRUE)
    if (vb) {
      vb_stage("nested quadrature tape", tf,
               paste0(nodes[[1L]], " x ", nodes[[2L]], " nodes, ",
                      qo$n_stack, " stacked rows"))
    }
    list(obj = obj, value = obj$fn(at))
  }
  units <- par_units %||% 1
  grad_at <- function(fz, at) {
    g <- try(max(abs(fz$obj$gr(at) * units)), silent = TRUE)
    if (inherits(g, "try-error")) NA_real_ else g
  }
  par <- lap_opt$par
  np <- length(par)
  lo <- rep_len(bounds$lower %||% -Inf, np)
  hi <- rep_len(bounds$upper %||% Inf, np)
  fz <- freeze(par, TRUE)
  grad <- grad_at(fz, par)
  settled <- is.finite(grad) && grad < control$grad_tol
  opt <- lap_opt
  moved <- 0
  radius <- 1
  rejected <- 0L
  k <- 0L
  while (!settled && k < rounds && rejected < 3L) {
    k <- k + 1L
    if (vb) t0 <- vb_now()
    box <- list(lower = pmax(lo, par - radius / units),
                upper = pmin(hi, par + radius / units))
    opt <- optimize_obj(fz$obj, control, box, par_units, verbose = vb,
                        start_par = par)
    cand <- opt$par
    step <- max(abs((cand - par) * units))
    f_old <- fz$value
    # the old tape is released before the new one is built: two alive
    # at once doubled the peak memory (dev/nestquad-findings.md)
    fz <- NULL
    gc(verbose = FALSE)
    fz <- freeze(cand, FALSE)
    ok <- is.finite(fz$value) && fz$value <= f_old + 1e-8 * max(1, abs(f_old))
    if (ok) {
      moved <- step
      par <- cand
      grad <- grad_at(fz, par)
      settled <- is.finite(grad) && grad < control$grad_tol
      if (step >= 0.9 * radius) radius <- 2 * radius
      rejected <- 0L
    } else {
      radius <- step / 4
      rejected <- rejected + 1L
      fz <- NULL
      gc(verbose = FALSE)
      fz <- freeze(par, FALSE)
      grad <- grad_at(fz, par)
    }
    if (vb) {
      vb_stage("nested quadrature round", t0,
               paste0(if (ok) "accepted" else "rejected", " a step of ",
                      format(step, digits = 3), ", objective ",
                      format(fz$value, digits = 8), ", re-anchored grad ",
                      format(grad, digits = 3), ", radius ",
                      format(radius, digits = 3)))
    }
    if (ok && step < 1e-8) break
  }
  opt$par <- par
  opt$objective <- fz$value
  quad_nest_check_bound(frame, opt$objective, nlp, template, par, map)
  rule_err <- NA_real_
  big <- NULL
  if (!settled) {
    # the size of the rule's own error at the estimate: the same integral
    # by a larger rule, on doubles (no tape). On the review's two
    # unsettled designs the fit was 1.7e-2 and 6.2e-2 from a reference
    big <- as.integer(ifelse(nodes < 21L, 21L, nodes + 10L))
    rule_err <- tryCatch({
      anc <- quad_nest_anchor(lap_obj, plan, par)
      pl <- imp_par_list(template, par, map)
      qb <- quad_nest_objective(frame, plan, anc, big, pl, nlp)
      qb$prime(pl)
      abs(qb$fn(pl) - fz$value)
    }, error = function(e) NA_real_)
    frm_warning("quadrature = TRUE did not settle the nested rule: after ",
                k, " round(s) the rule anchored at the estimate still has ",
                "a gradient of ", format(grad, digits = 3), " there, and ",
                "the last round moved the estimates by ",
                format(moved, digits = 3), ". The log-likelihood is the ",
                "rule anchored at the reported estimate, but that estimate ",
                "is not its maximum, and the rule's own error there is of ",
                "the order of ", format(rule_err, digits = 2), " (the ",
                "log-likelihood at ", nodes[[1L]], " x ", nodes[[2L]],
                " nodes against ", big[[1L]], " x ", big[[2L]], "). More ",
                "nodes (frmtmb_control(quad_nodes = )) are not a sure ",
                "remedy on such a design: the error need not fall as they ",
                "rise", call. = FALSE, class = "frmtmb_quad_unsettled")
  }
  list(obj = fz$obj, opt = opt,
       info = list(nodes = nodes, rounds = k, settled = settled,
                   moved = moved, grad = grad, rule_error = rule_err,
                   rule_error_nodes = big, memory = mem))
}

#' A log-likelihood above its family's bound is not reported. For a
#' discrete response every row's probability is at most 1, so the
#' integral's logarithm is at most 0; a rule that gives more has been
#' evaluated where it does not integrate the model.
#'
#' @noRd
quad_nest_check_bound <- function(frame, objective, nlp, template, par,
                                  map) {
  rn <- names(frame[["spec"]]$responses)[[1L]]
  fam <- frame[["spec"]]$responses[[rn]]$family
  if (!identical(fam[["type"]], "discrete")) return(invisible())
  pr <- if (is.null(nlp)) 0 else nlp(imp_par_list(template, par, map))
  ll <- pr - objective
  if (!is.finite(ll) || ll > 1e-8 * frame[["n_obs"]]) {
    frm_stop("quadrature = TRUE reached a log-likelihood of ",
             format(ll, digits = 6), " on a discrete response, where every ",
             "row's probability is at most 1 and the log-likelihood at ",
             "most 0. The nested rule was evaluated where it does not ",
             "integrate the model. Use more nodes ",
             "(frmtmb_control(quad_nodes = )), or fit with quadrature = ",
             "FALSE, the Laplace approximation", call. = FALSE,
             class = "frmtmb_quad_bound")
  }
  invisible()
}

#' Refuse a method that reads a nested quadrature fit's tape away from
#' its estimate: TMB's profile and root search on it. The tape is a rule
#' for the integral only near the point it is anchored at; one profile
#' interval came out (0.520, 0.953) around an estimate of 1.171
#' (dev/nestquad-findings.md, "Fix round 1").
#'
#' @noRd
quad_nested_refuse_profile <- function(fit, what) {
  if (is.null(fit[["quad_nested"]])) return(invisible())
  frm_stop(what, " is not available on a fit whose quadrature = TRUE ",
           "integrates nested blocks: it would read the fit's tape away ",
           "from the estimate, where the rule is not anchored and does not ",
           "integrate the model, and give a wrong interval. Use the Wald ",
           "method, whose standard errors come from the rule at the ",
           "estimate, or the bootstrap, which refits by the same rule",
           call. = FALSE, class = "frmtmb_quad_nested_unsupported")
}
