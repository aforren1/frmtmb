# Two similarity-based models whose likelihood does NOT factorize over
# rows: the probability of one response reads every other item in the
# same set, because a similarity has to be normalized over the set
# before it is a decision probability.
#
# That is what frmtmb_structure() is for. See
# dev/structured-family-protocol.md for the slot contracts and
# inst/rl/rw-delta.R for the worked example the protocol documents.
# vignette("bayesian-cognitive-modeling") shows these in place;
# tests/testthat/test-bcm-gcm.R and tests/testthat/test-bcm-simple.R
# fit this code. Load it with:
#
#   source(system.file("bcm", "binomial-extras.R", package = "frmtmb"))
#   source(system.file("bcm", "similarity.R", package = "frmtmb"))
#
# binomial-extras.R comes first because SIMPLE caps its recall
# probability with bcm_cap1() from there.

## ---- bcm-set-block ----
# Rows grouped into sets, each set's rows in item order, padded to the
# longest set. This is inst/rl/rw-delta.R's block with trials replaced
# by items, and the same two rules apply: `mask` is 1 on a real item and
# 0 on a pad, and the matrix() before each transpose is load-bearing,
# because vapply() drops to a plain vector when there is one item and
# t() of a vector would read every set as one more item of a single set.
bcm_set_block <- function(group, order_by, n, what) {
  gv <- factor(group)
  if (anyNA(gv)) {
    stop(what, ": every row needs a set, and ", sum(is.na(gv)),
         " row(s) have none", call. = FALSE)
  }
  rows <- lapply(split(seq_len(n), gv), function(r) r[order(order_by[r])])
  key <- paste(as.integer(gv), order_by, sep = "|")
  if (anyDuplicated(key)) {
    stop(what, ": item positions must be unique within a set; two rows ",
         "share one, and a normalized similarity has no way to tell ",
         "them apart", call. = FALSE)
  }
  len <- lengths(rows)
  ni <- max(len)
  idx <- t(matrix(vapply(rows, function(r) c(r, rep(r[1L], ni - length(r))),
                         integer(ni)), nrow = ni))
  mask <- t(matrix(vapply(len, function(l) as.numeric(seq_len(ni) <= l),
                          numeric(ni)), nrow = ni))
  list(idx = idx, mask = mask, len = len, n_set = length(rows),
       n_item = ni, levels = levels(gv), n = n)
}

# Broadcast a dpar the core handed back as one value. rep(v, each = )
# strips the advector class; multiplying by a vector of ones does not.
bcm_bcast <- function(v, n) if (length(v) == 1L) v * rep(1, n) else v

## ---- bcm-gcm ----
# The generalized context model of Nosofsky, as Lee and Wagenmakers
# chapter 17 states it. Stimulus i is called category A with
# probability
#
#   r_i = b S1_i / (b S1_i + (1 - b) S2_i)
#
# where S1_i and S2_i sum the similarity of i to the stimuli of each
# category, and the similarity of i to j decays exponentially in a
# weighted distance:  s_ij = exp(-c (w d1_ij + (1 - w) d2_ij)).
#
# `c` is the generalization gradient and `w` the attention weight. The
# two distance matrices and the category assignment are the design of
# the experiment, not columns of the data frame, so they ride on the
# family the way hmm()'s state count does.
#
# The book bounds c on (0, 5). That bound is a prior, not part of the
# likelihood, so the family uses a log link and the fit is free to leave
# the interval; the test asserts that it does not.
bcm_gcm <- function(stimulus, d1, d2, a, b = 0.5) {
  stim_expr <- substitute(stimulus)
  d1 <- as.matrix(d1)
  d2 <- as.matrix(d2)
  if (!identical(dim(d1), dim(d2)) || nrow(d1) != ncol(d1)) {
    stop("bcm_gcm(): d1 and d2 must be the same square matrix of ",
         "pairwise distances", call. = FALSE)
  }
  if (length(a) != nrow(d1) || !all(a %in% c(1, 2))) {
    stop("bcm_gcm(): a gives each stimulus's category as 1 or 2, one ",
         "per row of d1", call. = FALSE)
  }
  fam <- frmtmb_family(
    "bcm_gcm",
    dpars = c("c", "w"),
    links = list(c = "log", w = "logit"),
    primary_dpars = "c",
    type = "discrete",
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A bcm_gcm() stimulus's decision probability is a similarity ",
           "normalized over the whole stimulus set, so the family has no ",
           "row-wise log density. Use logLik() for the total, or ",
           "fitted() for the per-stimulus category probabilities",
           call. = FALSE)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("bcm_gcm(): the response is the number of category A ",
             "decisions, an integer count in [0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(c = function(y, aterms) 1,
                      w = function(y, aterms) 0.5),
    structure = bcm_gcm_structure())
  fam[["gcm"]] <- list(stim_expr = stim_expr,
                       d1 = d1, d2 = d2, a = as.integer(a), b = b)
  fam
}

# The category probability of every row. One function, called from the
# loglik slot with AD parameters and from fitted_mean() with plain ones,
# which is what keeps the two from drifting apart.
bcm_gcm_probs <- function(block, cc, ww) {
  A1 <- as.numeric(block[["a"]] == 1L)
  A2 <- as.numeric(block[["a"]] == 2L)
  D1 <- block[["D1"]]
  D2 <- block[["D2"]]
  S1 <- 0
  S2 <- 0
  # one iteration per stimulus in the set, each a whole-column vector
  # operation over the rows: the shape RTMB rewards, and the reason
  # there is no element-wise assignment into a taped matrix here
  for (j in seq_len(ncol(D1))) {
    sj <- exp(-cc * (ww * D1[, j] + (1 - ww) * D2[, j]))
    S1 <- S1 + A1[j] * sj
    S2 <- S2 + A2[j] * sj
  }
  b <- block[["b"]]
  (b * S1) / (b * S1 + (1 - b) * S2)
}

bcm_gcm_structure <- function() {
  frmtmb_structure(
    frame_vars = function(fam) list(fam[["gcm"]][["stim_expr"]]),
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("bcm_gcm() supports univariate models only", call. = FALSE)
      }
      if (is.null(av[["trials"]])) {
        stop("bcm_gcm() needs trials(): the response counts category A ",
             "decisions out of a known number of presentations",
             call. = FALSE)
      }
    },
    frame_block = function(resp, spec, av, mf, y, n) {
      g <- resp$family[["gcm"]]
      sv <- eval(g[["stim_expr"]], mf, resp$formula_env)
      stim <- if (is.factor(sv)) as.integer(sv) else as.integer(sv)
      if (anyNA(stim) || any(stim < 1) || any(stim > nrow(g[["d1"]]))) {
        stop("bcm_gcm(): every row's stimulus must index a row of d1, ",
             "so 1 to ", nrow(g[["d1"]]), call. = FALSE)
      }
      # the design of the experiment travels in the block, because the
      # loglik slot sees the block and not the family
      list(stim = stim, n = n, a = g[["a"]], b = g[["b"]],
           D1 = g[["d1"]][stim, , drop = FALSE],
           D2 = g[["d2"]][stim, , drop = FALSE])
    },
    loglik = function(y, dpars, aterms, weights, block, extra) {
      n <- block[["n"]]
      r <- bcm_gcm_probs(block,
                         bcm_bcast(dpars[["c"]], n),
                         bcm_bcast(dpars[["w"]], n))
      sum(weights * RTMB::dbinom(y, aterms[["trials"]], r, log = TRUE))
    },
    unit = "one stimulus's decisions",
    fitted_mean = function(fit, block) {
      rs <- single_response(fit, "a bcm_gcm() fit")
      dp <- eval_dpars(fit)[[rs$resp_name]]
      av <- fit$frame[["aterm_values"]][[rs$resp_name]]
      n <- block[["n"]]
      r <- bcm_gcm_probs(block,
                         rep(as.numeric(dp[["c"]]), length.out = n),
                         rep(as.numeric(dp[["w"]]), length.out = n))
      as.numeric(r) * av[["trials"]]
    },
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a bcm_gcm() fit is not available ",
        "for newdata: a stimulus's decision probability is a similarity ",
        "normalized over the whole stimulus set, and newdata carries no ",
        "set. Use fitted() on the training data"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a bcm_gcm() fit: ",
        "the expected response depends on the whole stimulus set, which ",
        "the synthetic grid this function builds does not have"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a bcm_gcm() fit: ",
        "the tape holds the whole set with no registered observation ",
        "vector. Use type = 'pearson'")))
}

## ---- bcm-simple ----
# The SIMPLE model of Brown, Neath and Chater, as Lee and Wagenmakers
# chapter 15 states it. Item i of a list is recalled with probability
#
#   theta_i = sum_j inv_logit(s (disc_ij - t)),
#   disc_ij = sim_ij / sum_k sim_ik,
#   sim_ij  = exp(-c |log m_i - log m_j|)
#
# where m_i is how long ago item i was presented. So a serial position
# curve, not a per-item rate: the discriminability of an item is its
# similarity relative to every other item in the SAME list, which is why
# this is a structured family and not a nonlinear formula.
#
# `m` is a column of the data frame, because it is a property of the
# row, and the list is the grouping factor.
# `link_t` is an argument because SIMPLE_2 needs it to be the identity:
# there the threshold is a1 * w + a2, a linear function of list length,
# and a logit link would not be that function. SIMPLE_1 keeps the logit,
# which is what holds one free threshold per list inside (0, 1).
bcm_simple <- function(list_id, position, m, link_t = "logit") {
  list_expr <- substitute(list_id)
  pos_expr <- substitute(position)
  m_expr <- substitute(m)
  fam <- frmtmb_family(
    "bcm_simple",
    dpars = c("c", "s", "t"),
    links = list(c = "log", s = "log", t = link_t),
    primary_dpars = "c",
    type = "discrete",
    lpdf = function(y, dpars, aterms, extra = NULL) {
      stop("A bcm_simple() item's recall probability is a similarity ",
           "normalized over the whole list, so the family has no ",
           "row-wise log density. Use logLik() for the total, or ",
           "fitted() for the per-item recall probabilities",
           call. = FALSE)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("bcm_simple(): the response is the number of correct ",
             "recalls, an integer count in [0, trials]", call. = FALSE)
      }
    },
    # The cap makes an infeasible START fatal rather than merely slow:
    # where theta reaches one, a row with fewer recalls than attempts
    # has a log density of -Infinity and a gradient of NaN. A threshold
    # of 0.3 is inside that region for a ten-item list, so the default
    # is high enough to keep every serial position under one.
    init_dpars = list(c = function(y, aterms) 20,
                      s = function(y, aterms) 10,
                      t = function(y, aterms) 0.6),
    structure = bcm_simple_structure())
  fam[["simple"]] <- list(list_expr = list_expr, pos_expr = pos_expr,
                          m_expr = m_expr)
  fam
}

# The recall probability of every row, capped at one the way the book
# caps it. bcm_cap1() is in inst/bcm/binomial-extras.R and says why the
# cap has to be carried rather than assumed inactive: without it a
# serial position where nearly everyone recalled drives the likelihood
# up without bound.
bcm_simple_probs <- function(block, cc, ss, tt) {
  idx <- block[["idx"]]
  msk <- block[["mask"]]
  lm <- block[["logm"]]
  theta <- 0 * cc
  for (g in seq_len(nrow(idx))) {
    r <- idx[g, msk[g, ] == 1]
    lmg <- lm[r]
    cg <- cc[r[1L]]
    sg <- ss[r[1L]]
    tg <- tt[r[1L]]
    L <- length(r)
    # the L-by-L similarity matrix, built one column at a time so that
    # nothing is assigned element-wise into a taped object
    sim <- vector("list", L)
    tot <- 0
    for (j in seq_len(L)) {
      sim[[j]] <- exp(-cg * abs(lmg - lmg[j]))
      tot <- tot + sim[[j]]
    }
    th <- 0
    for (j in seq_len(L)) {
      th <- th + RTMB::plogis(sg * (sim[[j]] / tot - tg))
    }
    theta[r] <- bcm_cap1(th)
  }
  theta
}

bcm_simple_structure <- function() {
  frmtmb_structure(
    frame_vars = function(fam) {
      list(fam[["simple"]][["list_expr"]], fam[["simple"]][["pos_expr"]],
           fam[["simple"]][["m_expr"]])
    },
    check_spec = function(resp, spec, av) {
      if (length(spec$responses) > 1L || isTRUE(spec$rescor)) {
        stop("bcm_simple() supports univariate models only", call. = FALSE)
      }
      if (is.null(av[["trials"]])) {
        stop("bcm_simple() needs trials(): the response counts correct ",
             "recalls out of a known number of attempts", call. = FALSE)
      }
    },
    frame_block = function(resp, spec, av, mf, y, n) {
      sp <- resp$family[["simple"]]
      gv <- eval(sp[["list_expr"]], mf, resp$formula_env)
      pv <- as.numeric(eval(sp[["pos_expr"]], mf, resp$formula_env))
      mv <- as.numeric(eval(sp[["m_expr"]], mf, resp$formula_env))
      if (any(mv <= 0)) {
        stop("bcm_simple(): m is a time since presentation and enters ",
             "through its logarithm, so it must be positive",
             call. = FALSE)
      }
      blk <- bcm_set_block(gv, pv, n, "bcm_simple()")
      blk[["logm"]] <- log(mv)
      blk
    },
    loglik = function(y, dpars, aterms, weights, block, extra) {
      n <- block[["n"]]
      theta <- bcm_simple_probs(block,
                                bcm_bcast(dpars[["c"]], n),
                                bcm_bcast(dpars[["s"]], n),
                                bcm_bcast(dpars[["t"]], n))
      sum(weights * RTMB::dbinom(y, aterms[["trials"]], theta, log = TRUE))
    },
    unit = "one list",
    fitted_mean = function(fit, block) {
      rs <- single_response(fit, "a bcm_simple() fit")
      dp <- eval_dpars(fit)[[rs$resp_name]]
      av <- fit$frame[["aterm_values"]][[rs$resp_name]]
      n <- block[["n"]]
      th <- bcm_simple_probs(block,
                             rep(as.numeric(dp[["c"]]), length.out = n),
                             rep(as.numeric(dp[["s"]]), length.out = n),
                             rep(as.numeric(dp[["t"]]), length.out = n))
      as.numeric(th) * av[["trials"]]
    },
    refusals = list(
      newdata_response = paste0(
        "predict(type = 'response') on a bcm_simple() fit is not ",
        "available for newdata: an item's recall probability is its ",
        "similarity relative to every other item of the same list, and ",
        "newdata carries no list. Use fitted() on the training data"),
      conditional_effects = paste0(
        "conditional_effects() is not available for a bcm_simple() fit: ",
        "the expected response depends on the whole list, which the ",
        "synthetic grid this function builds does not have"),
      osa = paste0(
        "residuals(type = 'osa') is not available for a bcm_simple() ",
        "fit: the tape holds the whole list with no registered ",
        "observation vector. Use type = 'pearson'")))
}
