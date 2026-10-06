# The draws surface in brms's shapes.
#
# frmtmb.sample is held to brms's OUTPUT, not only to its argument names
# (dev/brmsmatch-findings.md, dev/brmsnames-findings.md). The helpers
# here are the pieces every brms-shaped accessor shares: one route from
# the stored draws matrix to a posterior draws object, brms's `pars` and
# `subset` aliases with brms's own deprecation warnings, and the layout
# of the group-level coefficients that `ranef()`, `coef()` and
# `hypothesis(scope = )` read.

#' The stored draws as an iterations x chains x variables array.
#'
#' `frm_sample()` binds the chains in order, so the matrix is already
#' chain-major and reshapes without a permutation. This is the raw
#' array bayesplot's `mcmc_*` functions read; the public `as.array()`
#' is brms's instead.
#'
#' @noRd
draws_raw_array <- function(x) {
  m <- x$draws
  nc <- draws_nchains(x)
  if (nc <= 1L || nrow(m) %% nc != 0L) nc <- 1L
  array(as.vector(m), c(nrow(m) %/% nc, nc, ncol(m)),
        dimnames = list(NULL, paste0("chain:", seq_len(nc)), colnames(m)))
}

#' The draws as a `posterior::draws_array`, restricted the way brms's
#' `as_draws_*()` restrict: `variable` by exact name unless `regex`.
#'
#' Through the chain-separated array, so that a convergence diagnostic
#' computed on the result sees the chains. `inc_warmup = TRUE` is
#' refused, because the matrix keeps the post-warmup draws alone.
#'
#' @noRd
draws_as_array <- function(x, variable = NULL, regex = FALSE,
                           inc_warmup = FALSE, what = "this function") {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    frm_stop(what, " needs the 'posterior' package", call. = FALSE)
  }
  check_flag(regex, "regex")
  check_flag(inc_warmup, "inc_warmup")
  if (inc_warmup) {
    frm_stop(what, " cannot honor inc_warmup = TRUE: the draws matrix holds ",
             "the post-warmup draws only, which is what frm_sample() keeps. ",
             "The warmup, if the sampler saved it, is in `x$stanfit`",
             call. = FALSE)
  }
  a <- posterior::as_draws_array(draws_raw_array(x))
  if (!is.null(variable)) {
    a <- draws_subset_variable(a, variable, regex, what)
  }
  a
}

#' `posterior::subset_draws()` by variable, with an exact name the draws
#' do not carry refused here, as a classed refusal, rather than by
#' posterior in its own unclassed error. A regular expression that
#' matches nothing selects nothing in posterior, and still does.
#'
#' @noRd
draws_subset_variable <- function(a, variable, regex, what) {
  if (!isTRUE(regex)) {
    miss <- setdiff(variable, posterior::variables(a))
    if (length(miss)) {
      frm_stop(what, ": the following variables are missing in the draws ",
               "object: ", paste0("'", miss, "'", collapse = ", "),
               ". variables() lists every name; a regular expression ",
               "needs regex = TRUE", call. = FALSE)
    }
  }
  posterior::subset_draws(a, variable = variable, regex = regex)
}

#' The draws columns brms reports on the natural scale: a distributional
#' parameter nobody wrote a formula for is brms's `sigma`, `shape` or
#' `sigma_ya`, the parameter itself, so its column holds the parameter
#' and not its link. The sampler works on the link scale, and every
#' reader that hands a draw back to the model inverts these columns
#' first (`draws_internal_matrix()`). Memoized on the fit's cache
#' environment, because `draws_fit_at()` asks once per draw.
#'
#' `ordinal` holds the population-level parameters outside the
#' coefficient table (see `draws_ordinal_cols()`).
#'
#' @noRd
draws_natural_cols <- function(fit) {
  cache <- fit$cache
  if (is.environment(cache) && !is.null(cache$brms_natural_cols_v2)) {
    return(cache$brms_natural_cols_v2)
  }
  tab <- brms_coef_table(fit)
  smp <- attr(tab, "simplex")
  i <- setdiff(which(tab$natural),
               unlist(lapply(smp, `[[`, "pos")))
  # an equated dpar, bf(sigma1 = "sigma2"), has no parameter of its own
  # and is a copy of its target's column, as brms's transformed
  # parameter is
  equated <- lapply(attr(tab, "equated"), function(e) {
    list(name = e$name, from = tab$brms[e$row])
  })
  # a gp() term's sd and length scales: brms's sdgp_ and lscale_ draws
  # are on brms's scales, which the internal log theta is not
  gpn <- gp_brms_natural(fit)
  out <- list(names = c(tab$brms[i], gpn$names),
              linkinv = c(attr(tab, "linkinv")[i], gpn$linkinv),
              linkfun = c(attr(tab, "linkfun")[i], gpn$linkfun),
              simplex = smp, equated = equated,
              ordinal = draws_ordinal_cols(fit))
  # the columns added after the sampled ones: a mixture's last weight,
  # an equated dpar, and the thresholds and delta an ordinal block
  # reports beyond its internal parameters
  out$extra <- c(unlist(lapply(smp, function(s) {
    s$names[length(s$names)]
  })), vapply(equated, `[[`, "", "name"),
  unlist(lapply(out$ordinal, function(o) {
    o$names[-seq_along(o$internal)]
  })))
  if (is.environment(cache)) cache$brms_natural_cols_v2 <- out
  out
}

#' An ordinal fit's thresholds and its `cs()` coefficients, under the
#' names `variables(fit)` gives them.
#'
#' The sampler holds them in their template components, `tau_raw` (for
#' `cumulative()` the first threshold and log increments) and one
#' `bcs<j>` per `cs()` design column, and `brms_par_labels()` names them
#' so: `tau_raw_1`, `bcs2_1`. brms and `variables(fit)` call them
#' `b_Intercept[1]` and `bcs_fcb[1]`, and brms's `b_Intercept[k]` draws
#' are the thresholds themselves. Each block therefore carries the map
#' to the reported values and its inverse, both declared by the family
#' through core's `brms_fixef_rows()`. A block whose family declares no
#' inverse keeps the internal names, because a column must not carry a
#' name that its values do not have.
#'
#' A threshold structure other than `"flexible"` holds fewer internal
#' parameters than it reports: `"equidistant"` holds two per vector and
#' reports every threshold and brms's `delta`, and `"sum_to_zero"`
#' holds one fewer than its thresholds. The first `length(internal)`
#' names replace the internal columns in place and the rest are added
#' after every sampled column, as a mixture's last weight is, so no
#' sampled column moves; `inv` reads them all back.
#'
#' @noRd
draws_ordinal_cols <- function(fit) {
  # the formula route's fit has no estimates of its own; the blocks are
  # read off the layout, and any vector of the right length serves
  if (!length(fit$estimates)) fit$estimates <- fit$frame[["par_template"]]
  dl <- ord_delta_info(fit)
  out <- list()
  for (e in brms_fixef_rows(fit)$extra) {
    if (!is.function(e$inv)) next
    dd <- Filter(function(d) identical(d$comp, e$comp), dl)
    out[[length(out) + 1L]] <- draws_ordinal_block(e, dd)
  }
  out
}

#' One block of `draws_ordinal_cols()`. A function of its own so that
#' each block's maps close over that block's thresholds and `delta`s: a
#' closure built in the loop body reads the loop's variables when it
#' is called, after the loop has moved on to the last block, which gave
#' every block of a multivariate model the last response's count and
#' `delta` maps.
#'
#' @noRd
draws_ordinal_block <- function(e, dd) {
  K <- length(e$names)
  map <- e$map
  inv <- e$inv
  list(internal = paste0(e$comp, "_", seq_along(e$raw),
                         recycle0 = TRUE),
       names = c(paste0(e$cls, "_", e$names),
                 vapply(dd, `[[`, "", "name")),
       map = function(r) {
         c(map(r), vapply(dd, function(d) d$value(r[d$idx]), 0))
       },
       inv = function(v) inv(v[seq_len(K)]))
}

#' `f` applied to every row of `M`, one vector of the same length each.
#'
#' @noRd
draws_rowmap <- function(M, f) {
  k <- if (nrow(M)) length(f(M[1L, ])) else ncol(M)
  v <- vapply(seq_len(nrow(M)), function(r) as.numeric(f(M[r, ])),
              numeric(k))
  matrix(v, nrow(M), k, byrow = TRUE,
         dimnames = if (k == ncol(M)) dimnames(M))
}

#' Map the natural-scale columns of a draws matrix from the link scale
#' (`inverse = FALSE`, done once when the draws are stored) or back to
#' it (`inverse = TRUE`).
#'
#' A mixture's weights map jointly: the `K - 1` log-ratio columns become
#' brms's `theta1 ... theta<K-1>`, and brms's `thetaK`, which the sampler
#' has no parameter for, is added after every sampled column and before
#' `lp__`, so each sampled column keeps its position. The inverse reads
#' all `K` shares and drops the added column again.
#'
#' An ordinal block is renamed as well as mapped, in place. The inverse
#' maps the values back and leaves the names, because the model reads a
#' draw by position. A draws object stored before the renaming carries
#' the internal names on the internal scale, and the inverse leaves it
#' as it is.
#'
#' @noRd
draws_to_natural <- function(m, fit, inverse = FALSE) {
  nc <- draws_natural_cols(fit)
  if (inverse && length(nc$equated)) {
    eq <- match(vapply(nc$equated, `[[`, "", "name"), colnames(m))
    if (length(eq <- eq[!is.na(eq)])) m <- m[, -eq, drop = FALSE]
  }
  for (o in nc$ordinal) {
    j <- match(if (inverse) o$names else o$internal, colnames(m))
    if (anyNA(j)) next
    R <- length(o$internal)
    v <- draws_rowmap(m[, j, drop = FALSE], if (inverse) o$inv else o$map)
    if (inverse) {
      m[, j[seq_len(R)]] <- v
      if (length(j) > R) m <- m[, -j[-seq_len(R)], drop = FALSE]
      next
    }
    m[, j] <- v[, seq_len(R), drop = FALSE]
    colnames(m)[j] <- o$names[seq_len(R)]
    if (length(o$names) > R) {
      # before lp__, which brms and the sampler both keep last
      at <- match("lp__", colnames(m), nomatch = ncol(m) + 1L) - 1L
      nm_all <- append(colnames(m), o$names[-seq_len(R)], after = at)
      m <- cbind(m[, seq_len(at), drop = FALSE],
                 v[, -seq_len(R), drop = FALSE],
                 m[, seq_len(ncol(m) - at) + at, drop = FALSE])
      colnames(m) <- nm_all
    }
  }
  for (s in nc$simplex) {
    K <- length(s$names)
    j <- match(s$names[-K], colnames(m))
    if (anyNA(j)) next
    if (inverse) {
      jk <- match(s$names[K], colnames(m))
      if (is.na(jk)) {
        frm_stop("Internal error: the draws lack the mixture weight '",
                 s$names[K], "'. Please report it with the model formula",
                 call. = FALSE)
      }
      m[, j] <- s$to_logits(m[, c(j, jk), drop = FALSE])
      m <- m[, -jk, drop = FALSE]
    } else {
      P <- s$to_simplex(m[, j, drop = FALSE])
      m[, j] <- P[, -K]
      # before lp__, which brms and the sampler both keep last
      at <- match("lp__", colnames(m), nomatch = ncol(m) + 1L) - 1L
      nm_all <- append(colnames(m), s$names[K], after = at)
      m <- cbind(m[, seq_len(at), drop = FALSE], P[, K],
                 m[, seq_len(ncol(m) - at) + at, drop = FALSE])
      colnames(m) <- nm_all
    }
  }
  for (k in seq_along(nc$names)) {
    j <- which(colnames(m) == nc$names[k])
    if (!length(j)) next
    if (length(j) > 1L) {
      frm_stop("Internal error: the draws carry the column '", nc$names[k],
               "' twice. Please report it with the model formula",
               call. = FALSE)
    }
    f <- if (inverse) nc$linkfun[[k]] else nc$linkinv[[k]]
    m[, j] <- f(m[, j])
  }
  if (!inverse) {
    # an equated dpar is added the way thetaK is: after every sampled
    # column, before lp__, so no sampled column moves
    for (e in nc$equated) {
      j <- match(e$from, colnames(m))
      if (is.na(j) || e$name %in% colnames(m)) next
      at <- match("lp__", colnames(m), nomatch = ncol(m) + 1L) - 1L
      nm_all <- append(colnames(m), e$name, after = at)
      m <- cbind(m[, seq_len(at), drop = FALSE], m[, j],
                 m[, seq_len(ncol(m) - at) + at, drop = FALSE])
      colnames(m) <- nm_all
    }
  }
  m
}

#' The stored draws (or the rows `rows` of them) on the scale the model
#' reads: the natural-scale columns mapped back to their link.
#'
#' @noRd
draws_internal_matrix <- function(x, rows = NULL) {
  m <- x$draws
  if (!is.null(rows)) m <- m[rows, , drop = FALSE]
  draws_to_natural(m, x$fit, inverse = TRUE)
}

#' brms's `as.matrix()`, `as.array()` and `as.data.frame()` front end:
#' `pars` is the deprecated alias of `variable` and `subset` of `draw`,
#' each accepted with brms's own warning, and `regex`, `fixed` and
#' `inc_warmup` pass through as brms's `...` passes them. So do
#' `iteration` and `chain`, which brms hands to
#' `posterior::subset_draws()` with `draw`: `as.array(fit, chain = 1)`
#' is the first chain's iterations.
#'
#' `format` converts the draws array before the draw subset, because
#' brms subsets in the format the method returns.
#'
#' @noRd
draws_accessor_args <- function(x, pars, variable, draw, subset, what,
                                format = identity, ...) {
  frm_check_dots(..., .allow = c("regex", "fixed", "inc_warmup",
                                 "iteration", "chain"))
  dots <- list(...)
  if (!anyNA(pars)) {
    frm_warning("Argument 'pars' is deprecated. Please use 'variable' ",
                "instead.", call. = FALSE)
    variable <- draws_extract_pars(pars, colnames(x$draws),
                                   fixed = dots$fixed %||% FALSE)
  }
  if (!is.null(subset)) {
    frm_warning("Argument 'subset' is deprecated. Please use argument ",
                "'draw' instead.", call. = FALSE)
    draw <- subset
  }
  a <- format(draws_as_array(x, variable, dots$regex %||% FALSE,
                             dots$inc_warmup %||% FALSE, what))
  if (!is.null(draw) || !is.null(dots$iteration) || !is.null(dots$chain)) {
    a <- suppressMessages(posterior::subset_draws(
      a, draw = draw, iteration = dots$iteration, chain = dots$chain))
  }
  a
}

#' brms's `fixef_pars()`: the population-level coefficient prefixes.
#'
#' @noRd
draws_fixef_regex <- "^b((|s|cs|sp|mo|me|mi|m))_"

#' A derived draws matrix in the layout brms's `as.matrix()` returns
#' (named `draw` and `variable` margins, the chain count as an
#' attribute), so that `VarCorr(summary = FALSE)` hands back the same
#' object brms builds from its stored `sd_` draws.
#'
#' @noRd
draws_derived_matrix <- function(x, M, names) {
  nc <- nchains(x)
  if (nrow(M) %% nc != 0L) nc <- 1L
  arr <- array(as.vector(M), c(nrow(M) %/% nc, nc, ncol(M)),
               dimnames = list(NULL, NULL, names))
  out <- unclass(posterior::as_draws_matrix(posterior::as_draws_array(arr)))
  colnames(out) <- names
  out
}

#' Column positions of the parameters `hyp_vals_only()` reads, and their
#' components, for draws with or without the random effects.
#'
#' `draws_par_index()` assumes every template component is present,
#' which draws from `frm_sample(laplace = TRUE)` are not: they skip `b`
#' and `miss`, and every later position moves.
#'
#' @noRd
draws_outer_index <- function(x) {
  fit <- x$fit
  tpl <- fit$frame[["par_template"]]
  skip <- if (draws_is_laplace(x)) c("b", "miss") else character(0)
  pos <- 0L
  idx <- list()
  for (cp in names(tpl)) {
    if (cp %in% skip) next
    len <- length(tpl[[cp]])
    if (cp == "betad" && length(fx <- fit$frame[["betad_fixed_idx"]])) {
      len <- len - length(fx)
    }
    idx[[cp]] <- pos + seq_len(len)
    pos <- pos + len
  }
  comps <- c("beta", "betad", "theta", "thetaac", "thetar")
  cols <- unlist(idx[comps], use.names = FALSE)
  comp <- rep(comps, vapply(comps, function(k) length(idx[[k]]), 1L))
  list(cols = cols, comp = comp)
}

#' Positions of brms's correlation parameters inside a K x K matrix, in
#' brms's order (`brms:::get_cornames()`: row i from 2, column j below
#' i), which is not R's column-major lower triangle once K exceeds 3.
#'
#' @noRd
draws_cor_index <- function(K) {
  out <- integer(0)
  for (i in seq_len(K)[-1]) {
    for (j in seq_len(i - 1)) out <- c(out, (j - 1L) * K + i)
  }
  out
}

#' brms's `get_cor_matrix()`: draws x K x K correlation matrices from
#' draws of the correlations in brms's order.
#'
#' @noRd
draws_cor_array <- function(cor, size) {
  out <- aperm(array(diag(1, size), dim = c(size, size, nrow(cor))),
               perm = c(3, 1, 2))
  k <- 0
  for (i in seq_len(size)[-1]) {
    for (j in seq_len(i - 1)) {
      k <- k + 1
      out[, j, i] <- out[, i, j] <- cor[, k]
    }
  }
  out
}

#' brms's `get_cov_matrix()`: draws x K x K covariance matrices from sd
#' draws and lower-triangle correlation draws.
#'
#' @noRd
draws_cov_array <- function(sd, cor) {
  sd <- as.matrix(sd)
  size <- ncol(sd)
  out <- aperm(array(diag(1, size), dim = c(size, size, nrow(sd))),
               perm = c(3, 1, 2))
  for (i in seq_len(size)) out[, i, i] <- sd[, i]^2
  k <- 0
  for (i in seq_len(size)[-1]) {
    for (j in seq_len(i - 1)) {
      k <- k + 1
      out[, j, i] <- out[, i, j] <- cor[, k] * sd[, i] * sd[, j]
    }
  }
  out
}

#' `varcorr_values()` at every draw.
#'
#' @noRd
draws_varcorr_values <- function(x, lay) {
  fit <- draws_base_fit(x)
  oi <- draws_outer_index(x)
  m <- draws_internal_matrix(x)
  lapply(seq_len(nrow(m)), function(i) {
    varcorr_values(fit, unname(m[i, oi$cols]), oi$comp, lay)
  })
}

#' Where each group-level coefficient draw lives, per grouping factor in
#' brms's `ranef()` layout: the levels (brms's, an interaction group's
#' parts joined with `_`), brms's coefficient names (with a
#' distributional parameter's prefix, `sigma_Intercept`), the renamed
#' coefficient `pars` matches, and a levels x coefficients matrix of
#' column positions in the draws matrix.
#'
#' A block whose sampled values are not the coefficients themselves (a
#' reduced-rank block samples factor scores) has `NA` positions, and
#' `block` and `pos` say which block and which coefficient to compute
#' them from; see `draws_ranef_fill()`.
#'
#' A `gr(g, by = f)` term is one block per by-level, each over the levels
#' of `g` in that by-level, and `VarCorr()` names its coefficients per
#' by-level (`Intercept:fa`). brms's `ranef()` has one column per
#' coefficient over every level of `g`, so the by-levels of one
#' coefficient fill one column here, each at its own levels.
#'
#' @noRd
draws_ranef_layout <- function(x) {
  fit <- draws_base_fit(x)
  lay <- varcorr_layout(fit)$groups
  if (!length(lay)) return(list())
  blocks <- fit$frame[["re_blocks"]]
  laplace <- draws_is_laplace(x)
  b_pos <- if (laplace) NULL else draws_par_index(fit)$b
  out <- list()
  for (g in names(lay)) {
    L <- lay[[g]]
    bk1 <- blocks[[L$block[1L]]]
    by1 <- NULL
    for (bk in blocks[L$block]) by1 <- by1 %||% bk[["by"]]
    levels <- by1[["group_levels"]] %||% bk1[["levels"]]
    parts <- lapply(seq_along(L$rnames), function(k) {
      pr <- brms_re_parts(fit, blocks[[L$block[k]]])
      list(coef = pr$coef[L$pos[k]], rcoef = pr$rcoef[L$pos[k]])
    })
    # one column per coefficient, the by-levels of a by-split term merged
    key <- vapply(seq_along(L$rnames), function(k) {
      by <- blocks[[L$block[k]]][["by"]]
      if (is.null(by)) paste0("\r", k) else
        paste0(by$key, "\r", parts[[k]]$rcoef)
    }, "")
    col_of <- match(key, unique(key))
    first <- !duplicated(key)
    bare <- vapply(parts[first], `[[`, "", "coef")
    coefs <- ifelse(vapply(blocks[L$block[first]], function(bk) {
      is.null(bk[["by"]])
    }, NA), L$rnames[first], vapply(parts[first], `[[`, "", "rcoef"))
    cols <- matrix(NA_integer_, length(levels), length(bare))
    for (k in seq_along(L$rnames)) {
      bk <- blocks[[L$block[k]]]
      p <- L$pos[k]
      if (is.null(b_pos) || !brms_block_has_r(bk)) next
      li <- match(levels, bk[["levels"]])
      ok <- !is.na(li)
      bi <- bk[["b_idx"]][(li[ok] - 1L) * bk[["dim"]] + p]
      cols[which(ok), col_of[k]] <- b_pos[bi]
    }
    bkl <- bk1
    bkl[["levels"]] <- levels
    lv <- brms_levels(bkl)
    # brms keeps each level's by-level on its level names
    # (frame_re(): attr(levels, "by")), and ranef() carries them over
    if (!is.null(by1)) attr(lv, "by") <- gsub("[ \t\r\n]+", "", by1$level_by)
    out[[g]] <- list(levels = lv, raw_levels = levels,
                     coefs = coefs, coef = bare, cols = cols,
                     block = L$block[first], pos = L$pos[first])
  }
  out
}

#' Fill the coefficients of one `ranef()` group that no draws column
#' holds, per draw, from the coefficient vector the model itself builds
#' (`expand_b()`: a reduced-rank block's loadings times its factor
#' scores). Draws taken with the random effects integrated out have no
#' group-level draws at all, and those are refused by name rather than
#' returned as `NA`.
#'
#' @noRd
draws_ranef_fill <- function(x, L, sel, A, g) {
  miss <- sel[colSums(is.na(L$cols[, sel, drop = FALSE])) > 0L]
  if (!length(miss)) return(A)
  if (draws_is_laplace(x)) {
    frm_stop(draws_caller("ranef()"), " has no draws of the group-level ",
             "coefficients of '", g, "': these draws come from ",
             "frm_sample(laplace = TRUE), which integrates the random ",
             "effects out. Sample without laplace = TRUE for them",
             call. = FALSE)
  }
  fit <- draws_base_fit(x)
  frame <- fit$frame
  idx <- draws_index(x)
  m <- draws_internal_matrix(x)
  blocks <- frame[["re_blocks"]]
  for (i in seq_len(nrow(m))) {
    cvec <- expand_b(frame, unname(m[i, idx$b]), unname(m[i, idx$theta]))
    for (k in miss) {
      bk <- blocks[[L$block[k]]]
      li <- match(L$raw_levels, bk[["levels"]])
      ok <- !is.na(li)
      j <- match(k, sel)
      A[i, which(ok), j] <-
        cvec[bk[["c_idx"]][(li[ok] - 1L) * bk[["dim"]] + L$pos[k]]]
    }
  }
  A
}

#' One hypothesis summarized over a draws vector, brms's way
#' (`brms:::eval_hypothesis()`): mean and SD or median and MAD, brms's
#' interval, the evidence ratio of a directional claim, and its
#' posterior probability.
#'
#' @noRd
draws_hyp_summary <- function(v, dir, alpha, robust) {
  lo <- if (dir == "two.sided") alpha / 2 else alpha
  er <- switch(dir,
               less = sum(v < 0) / (length(v) - sum(v < 0)),
               greater = sum(v > 0) / (length(v) - sum(v > 0)),
               NA_real_)
  list(est = if (robust) stats::median(v) else mean(v),
       err = if (robust) stats::mad(v) else stats::sd(v),
       lwr = unname(stats::quantile(v, lo)),
       upr = unname(stats::quantile(v, 1 - lo)),
       er = er, pp = if (is.infinite(er)) 1 else er / (1 + er))
}

#' `hypothesis(scope = "ranef")` and `scope = "coef"`: brms's
#' `hypothesis_coef()`, which evaluates each hypothesis on every level
#' of one grouping factor, over that level's coefficient draws, with
#' the coefficient names as the variables.
#'
#' @noRd
draws_hypothesis_coef <- function(x, hypothesis, group, scope, alpha,
                                  robust) {
  co <- if (scope == "ranef") ranef(x, summary = FALSE) else
    coef(x, summary = FALSE)
  if (!length(group) || !group %in% names(co)) {
    frm_stop("'group' should be one of ", paste(names(co), collapse = ", "),
             call. = FALSE)
  }
  A <- co[[group]]
  levels <- dimnames(A)[[2L]]
  coefs <- dimnames(A)[[3L]]
  labels <- hyp_labels(hypothesis)
  # brms's hypothesis_coef(): class "" over the coefficient names, with
  # brms's renaming, so a coefficient carrying `:` reads as one name
  hp <- hyp_parse_all(hypothesis, coefs, class = NULL)
  rows <- list()
  samp <- list()
  for (h in seq_along(hypothesis)) {
    for (l in seq_along(levels)) {
      vals <- stats::setNames(lapply(seq_along(coefs), function(k) {
        A[, l, k]
      }), coefs)
      v <- hyp_eval_in(hp$exprs[[h]], vals)
      s <- draws_hyp_summary(as.numeric(v), hp$dir[h], alpha, robust)
      rows[[length(rows) + 1L]] <- c(s, list(dir = hp$dir[h],
                                              group = levels[l],
                                              label = labels[h]))
      samp[[length(samp) + 1L]] <- as.numeric(v)
    }
  }
  get <- function(nm) vapply(rows, function(r) r[[nm]], rows[[1L]][[nm]])
  dir <- get("dir")
  out <- hyp_brms_result(
    get("label"), get("est"), get("err"), get("lwr"), get("upr"),
    get("er"), get("pp"), dir, get("pp") > 1 - alpha,
    hyp_samples_frame(do.call(cbind, samp), length(samp)),
    hyp_samples_frame(matrix(NA, nrow(A), length(samp)), length(samp)),
    "", alpha, list(method = "posterior"))
  hs <- out$hypothesis
  hs$Group <- factor(get("group"), levels)
  out$hypothesis <- hs[c("Group", setdiff(names(hs), "Group"))]
  out
}
