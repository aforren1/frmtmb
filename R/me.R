# brms me(): noise-free predictors measured with known error.
#
# `me(x, sdx)` says that the column `x` is a noisy reading of a latent
# value Xn with known measurement SD `sdx`, and the model uses Xn, not
# x. brms 2.23.0 (stan_Xme(), frame_me(), data_Xme()) writes
#
#   Xn_k ~ normal(meanme_k, sdme_k)       one latent value per row
#   x_k  ~ normal(Xn_k, sdx_k)            the measurement model
#
# with `meanme` and `sdme` estimated, one pair per me() variable. With
# `gr = g` there is one latent value per LEVEL of g instead, and x and
# sdx must be constant within a level. Several me() terms over the same
# grouping have correlated latent values unless set_mecor(FALSE) says
# otherwise: the rows (Xn_1, ..., Xn_M) are multivariate normal with the
# `corme` correlation matrix.
#
# Here the latent values are inner parameters, integrated out by the
# Laplace approximation beside the random effects. They share the `miss`
# component of the parameter template with mi(), so every post-fit path
# that already treats the mi() latent values as inner treats these the
# same way. The hyperparameters are three outer components of their
# own: `meanme` (natural scale), `logsdme` (log of sdme) and `thetame`
# (the unconstrained correlation parameters of us_chol_cor()).

#' The formals of brms's `me()`, used to name a call's arguments the way
#' brms's own `match.call()` does.
#'
#' @noRd
me_formals <- function(x, sdx, gr = NULL) NULL

#' One `me()` call, read: the noisy variable, its measurement SD, the
#' optional grouping variable, and the names brms derives from them.
#'
#' `key` is the call as written, which is what brms's `get_uni_me()`
#' compares, so `me(x, sx)` and `me(x, sdx = sx)` are two different
#' calls that use one variable and are refused, as brms refuses them.
#' `coef` is brms's `rename(paste0("me", xname))`, the suffix of
#' `meanme_<coef>`, `sdme_<coef>` and `Xme_<coef>`.
#'
#' @noRd
parse_me_call <- function(e) {
  mc <- tryCatch(match.call(me_formals, e), error = function(err) {
    frm_stop("Cannot read ", deparse1(e), ": me() takes the noisy ",
             "variable, its measurement SD and an optional grouping ",
             "variable, me(x, sdx) or me(x, sdx, gr = g)", call. = FALSE)
  })
  if (is.null(mc[["x"]])) {
    frm_stop("Argument 'x' is missing in function 'me': ", deparse1(e),
             call. = FALSE)
  }
  if (is.null(mc[["sdx"]])) {
    # brms's own words (get_me_noise())
    frm_stop("Argument 'sdx' is missing in function 'me': ", deparse1(e),
             ". Give the known measurement SD, me(x, sdx)", call. = FALSE)
  }
  gr <- mc[["gr"]]
  if (!is.null(gr)) {
    if (!is.name(gr)) {
      frm_stop("me(gr = ) takes the name of one grouping variable, not ",
               "`", deparse1(gr), "`. Build the column first, then name ",
               "it", call. = FALSE)
    }
    gr <- as.character(gr)
  }
  xname <- deparse1(mc[["x"]])
  list(key = deparse1(e), x = mc[["x"]], sdx = mc[["sdx"]], gr = gr,
       xname = xname, coef = brms_rename(paste0("me", xname)))
}

#' The factors of one `:` interaction, in written order.
#'
#' @noRd
me_split_colon <- function(e) {
  if (is.call(e) && identical(e[[1L]], as.name(":"))) {
    c(me_split_colon(e[[2L]]), me_split_colon(e[[3L]]))
  } else {
    list(e)
  }
}

#' Is this expression a bare `me()` call?
#'
#' @noRd
is_me_call <- function(e) {
  is.call(e) && identical(e[[1L]], as.name("me"))
}

#' One top-level formula term that calls `me()`, split into its special
#' parts and the plain terms it implies.
#'
#' The term is expanded the way `terms()` expands it, which is how brms
#' reads its special terms, so `me(x, sx) * z` is the three terms
#' `me(x, sx)`, `z` and `me(x, sx):z`. Each label that uses `me()` must
#' be a product of `me()` calls and plain numeric factors. A product of
#' several `me()` calls multiplies their latent values, as brms's
#' `me(x, sx):me(z, sz)` does. The label keeps `terms()` spelling, and
#' brms's `rename()` of it is the coefficient name, `bsp_mexsx:z`.
#'
#' @noRd
parse_me_term <- function(tm) {
  if (any(c("|", "||") %in% all.names(tm))) {
    frm_stop("me() is not supported inside a group-level term: ",
             deparse1(tm), ". A varying slope on a noise-free variable ",
             "has no implementation here; write the population-level ",
             "term me(x, sdx) on its own", call. = FALSE)
  }
  labs <- cs_term_labels(tm)
  if (is.null(labs)) {
    frm_stop("Cannot expand the term ", deparse1(tm), " that uses me()",
             call. = FALSE)
  }
  entries <- list()
  rest <- list()
  for (lab in labs) {
    e <- str2lang(lab)
    if (!calls_function(e, "me")) {
      rest[[length(rest) + 1L]] <- e
      next
    }
    fac <- me_split_colon(e)
    is_me <- vapply(fac, is_me_call, TRUE)
    for (f in fac[!is_me]) {
      if (calls_function(f, "me")) {
        frm_stop("me() must be a term of its own or a factor of a ':' ",
                 "interaction, not part of another expression: ",
                 deparse1(f), " in the term '", lab, "'. Transform the ",
                 "noisy variable inside the call instead, me(log(x), sdx), ",
                 "where the measurement SD then belongs to log(x)",
                 call. = FALSE)
      }
      if (any(vapply(c("mo", "mi"), function(s) calls_function(f, s), TRUE))) {
        frm_stop("me() cannot share an interaction with mo() or mi(): '",
                 lab, "'. Each of those may appear as its own term in the ",
                 "same formula", call. = FALSE)
      }
    }
    # the plain factors are kept apart: each is checked as one numeric
    # column, and their product is the multiplier. Joined with `:` they
    # would evaluate as R's sequence operator
    entries[[length(entries) + 1L]] <- list(
      calls = lapply(fac[is_me], parse_me_call), mult = fac[!is_me],
      label = lab, order = length(fac))
  }
  list(entries = entries, rest = rest)
}

#' Every distinct `me()` call of a model, in brms's order.
#'
#' brms keeps ONE set of latent values per distinct call for the whole
#' model, shared by every distributional and nonlinear parameter and by
#' every response that uses the call (`get_uni_me()` over all special
#' terms), and it refuses a variable written in two different calls.
#' Both rules are reproduced here, with brms's own message.
#'
#' @noRd
me_registry <- function(spec) {
  calls <- list()
  for (resp in spec$responses) {
    for (dp in resp[["dpars"]]) {
      for (ent in me_entries_in_brms_order(dp[["meterms"]])) {
        for (cl in ent$calls) {
          if (!cl$key %in% names(calls)) calls[[cl$key]] <- cl
        }
      }
    }
  }
  if (!length(calls)) return(list())
  xn <- vapply(calls, `[[`, "", "xname")
  if (anyDuplicated(xn)) {
    x1 <- xn[duplicated(xn)][1L]
    frm_stop("Variable '", x1, "' is used in different calls to 'me'.\n",
             "Associated calls are: ",
             paste0("'", names(calls)[xn == x1], "'", collapse = ", "),
             call. = FALSE)
  }
  unname(calls)
}

#' The `me()` entries of one linear predictor in brms's `terms()` order:
#' every main effect before every interaction, written order otherwise,
#' as for `mo()` and `mi()`. A label written twice is one term.
#'
#' @noRd
me_entries_in_brms_order <- function(entries) {
  if (!length(entries)) return(list())
  labs <- vapply(entries, `[[`, "", "label")
  entries <- entries[!duplicated(labs)]
  ord <- vapply(entries, `[[`, 1L, "order")
  entries[order(ord)]
}

#' The variables one dpar's `me()` terms need in the combined model
#' frame, as names, so that a row with a missing value in any of them is
#' dropped with the rest of the row (brms does the same).
#'
#' @noRd
me_frame_vars <- function(dp) {
  v <- character(0)
  for (ent in dp[["meterms"]] %||% list()) {
    for (cl in ent$calls) {
      v <- c(v, all.vars(cl$x), all.vars(cl$sdx), cl$gr)
    }
    for (f in ent$mult) v <- c(v, all.vars(f))
  }
  unique(v)
}

#' The measurement data and the latent layout of every `me()` call.
#'
#' Evaluates each call's variable and SD against the model frame with
#' brms's checks and messages (`get_me_values()`, `get_me_noise()`,
#' `data_Xme()`), groups the calls by their `gr` variable in order of
#' first appearance, and lays the latent values out after the `n_miss`
#' values the `mi()` responses already hold. A group of several calls is
#' correlated when `mecor` is TRUE, and then owns
#' `M * (M - 1) / 2` entries of `thetame`.
#'
#' @noRd
me_build_frame <- function(spec, mf, env, n, n_miss) {
  calls <- me_registry(spec)
  if (!length(calls)) return(NULL)
  mecor <- !isFALSE(spec[["mecor"]])
  terms <- list()
  for (k in seq_along(calls)) {
    cl <- calls[[k]]
    xv <- eval(cl$x, mf, env)
    if (!is.numeric(xv) || is.factor(xv)) {
      frm_stop("Noisy variables should be numeric: '", cl$xname, "' in ",
               cl$key, " is ", class(xv)[1L], call. = FALSE)
    }
    xv <- as.numeric(xv)
    if (length(xv) != n) {
      frm_stop("The noisy variable of ", cl$key, " has ", length(xv),
               " values for ", n, " rows", call. = FALSE)
    }
    sv <- eval(cl$sdx, mf, env)
    if (!is.numeric(sv) || is.factor(sv)) {
      frm_stop("Measurement error should be numeric: ",
               deparse1(cl$sdx), " in ", cl$key, " is ", class(sv)[1L],
               call. = FALSE)
    }
    sv <- as.numeric(sv)
    if (length(sv) == 1L) sv <- rep(sv, n)
    if (length(sv) != n) {
      frm_stop("The measurement SD of ", cl$key, " has ", length(sv),
               " values for ", n, " rows", call. = FALSE)
    }
    if (any(sv <= 0)) {
      frm_stop("Measurement error should be positive: ", deparse1(cl$sdx),
               " in ", cl$key, " has ", sum(sv <= 0), " value(s) <= 0. ",
               "A variable measured without error is an ordinary ",
               "predictor", call. = FALSE)
    }
    J <- NULL
    levels <- NULL
    if (!is.null(cl$gr)) {
      gv <- eval(as.name(cl$gr), mf, env)
      levels <- levels(factor(gv))
      J <- match(as.character(gv), levels)
      # brms keeps the first value of each level, after checking that
      # every row of the level carries the same one
      first <- match(seq_along(levels), J)
      for (l in seq_along(levels)) {
        take <- J == l
        if (length(unique(xv[take])) > 1L || length(unique(sv[take])) > 1L) {
          frm_stop("Measured values and measurement error should be ",
                   "unique for each group. Occured for level '", levels[l],
                   "' of group '", cl$gr, "' in ", cl$key, call. = FALSE)
        }
      }
      xv <- xv[first]
      sv <- sv[first]
    }
    terms[[k]] <- c(cl, list(Xn = xv, noise = sv, J = J, levels = levels))
  }
  gkey <- vapply(terms, function(t) t$gr %||% "", "")
  groups <- list()
  off <- n_miss
  n_theta <- 0L
  for (g in unique(gkey)) {
    K <- which(gkey == g)
    N <- length(terms[[K[1L]]]$Xn)
    M <- length(K)
    idx <- matrix(off + seq_len(N * M), N, M)
    off <- off + N * M
    cor <- mecor && M > 1L
    th_idx <- integer(0)
    if (cor) {
      th_idx <- n_theta + seq_len(M * (M - 1L) / 2L)
      n_theta <- n_theta + length(th_idx)
    }
    for (j in seq_along(K)) terms[[K[j]]]$idx <- idx[, j]
    groups[[length(groups) + 1L]] <- list(
      gr = g, K = K, N = N, idx = idx, cor = cor, th_idx = th_idx,
      J = terms[[K[1L]]]$J, levels = terms[[K[1L]]]$levels,
      Xn = unlist(lapply(terms[K], `[[`, "Xn"), use.names = FALSE),
      noise = unlist(lapply(terms[K], `[[`, "noise"), use.names = FALSE))
  }
  coefs <- vapply(terms, `[[`, "", "coef")
  start_sd <- vapply(terms, function(t) {
    s <- stats::sd(t$Xn)
    if (is.finite(s) && s > 0) log(s) else 0
  }, 0)
  list(terms = terms, groups = groups, mecor = mecor,
       n_latent = off - n_miss,
       latent_init = unlist(lapply(groups, function(g) g$Xn),
                            use.names = FALSE),
       meanme = stats::setNames(vapply(terms, function(t) mean(t$Xn), 0),
                                paste0("meanme_", coefs)),
       logsdme = stats::setNames(start_sd, paste0("logsdme_", coefs)),
       n_thetame = n_theta)
}

#' The `me()` columns of one linear predictor: a zero placeholder column
#' per term, as for `mo()` and `mi()`, whose value the objective and the
#' post-fit paths supply from the latent values.
#'
#' @noRd
me_lp_columns <- function(dp, me, X, mf, env) {
  info <- list()
  if (is.null(me)) return(list(X = X, info = info))
  keys <- vapply(me$terms, `[[`, "", "key")
  for (ent in me_entries_in_brms_order(dp[["meterms"]])) {
    mult <- NULL
    for (f in ent$mult) {
      v <- check_special_mult(eval(f, mf, env), f, "me")
      mult <- if (is.null(mult)) v else mult * v
    }
    lab <- brms_rename(ent$label)
    X <- cbind(X, matrix(0, nrow(X), 1, dimnames = list(NULL, lab)))
    info[[length(info) + 1L]] <- list(
      terms = match(vapply(ent$calls, `[[`, "", "key"), keys),
      col = ncol(X), label = lab, written = ent$label,
      xvars = unique(unlist(lapply(ent$calls, function(cl) all.vars(cl$x)))),
      mult = mult, mult_factors = ent$mult,
      mult_expr = if (length(ent$mult)) {
        Reduce(function(a, b) call("*", a, b), ent$mult)
      })
  }
  list(X = X, info = info)
}

#' The latent values of every `me()` call and their log-density.
#'
#' Runs on the tape (inside the objective) and off it (post-fit, with
#' the estimates), which is why it only uses operations RTMB overloads.
#' `values[[k]]` is call k's latent value per ROW of the data, expanded
#' through the level index when the call has `gr =`.
#'
#' The density is brms's: the measurement model `normal(Xn | Xme,
#' noise)` for every latent value, and the latent model, which is
#' `normal(Xme | meanme, sdme)` per call or, for a correlated group,
#' multivariate normal per row with correlation `us_chol_cor(thetame)`.
#' brms samples the standardized `zme` and pays no Jacobian; the latent
#' values here ARE the inner parameters, so the density is on them
#' directly, and the two joint densities differ by the constant
#' `N * sum(log(sdme))` of that change of variables.
#'
#' @noRd
me_latent <- function(me, pars) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  miss <- pars[["miss"]]
  mu <- pars[["meanme"]]
  lsd <- pars[["logsdme"]]
  ll <- 0
  values <- vector("list", length(me$terms))
  for (g in me$groups) {
    K <- g$K
    N <- g$N
    M <- length(K)
    xl <- miss[as.vector(g$idx)]
    ll <- ll + sum(RTMB::dnorm(g$Xn, xl, g$noise, log = TRUE))
    mrep <- rep(mu[K], each = N)
    srep <- exp(rep(lsd[K], each = N))
    if (g$cor) {
      Z <- RTMB::matrix((xl - mrep) / srep, N, M)
      C <- us_chol_cor(pars[["thetame"]][g$th_idx], M)
      ll <- ll + sum(RTMB::dmvnorm(Z, 0, C, log = TRUE)) - N * sum(lsd[K])
    } else {
      ll <- ll + sum(RTMB::dnorm(xl, mrep, srep, log = TRUE))
    }
    for (j in seq_len(M)) {
      v <- miss[g$idx[, j]]
      values[[K[j]]] <- if (is.null(g$J)) v else v[g$J]
    }
  }
  list(ll = ll, values = values)
}

#' One `me()` column's value per row: the product of its calls' latent
#' values and of its plain multiplier.
#'
#' @noRd
me_col_value <- function(mt, values) {
  v <- values[[mt$terms[1L]]]
  for (k in mt$terms[-1L]) v <- v * values[[k]]
  if (!is.null(mt$mult)) v <- v * mt$mult
  v
}

#' `me()` column values in new data. The observed value stands in for
#' the latent one, as the mean of brms's own draw does: brms predicts new
#' rows from `normal(x, sdx)` draws, whose expectation is `x`, and it
#' has no latent value for a row it did not fit. So `sdx` is not needed
#' here, and a missing `x` is refused rather than predicted at the
#' population mean.
#'
#' @noRd
me_newdata_value <- function(fit, mt, newdata, env, nd_mult) {
  terms <- fit$frame[["me"]]$terms
  v <- 1
  for (k in mt$terms) {
    t <- terms[[k]]
    xv <- tryCatch(eval(t$x, newdata, env), error = function(e) NULL)
    if (is.null(xv) || !is.numeric(xv) || anyNA(xv)) {
      frm_stop(t$key, ": newdata must supply complete numeric values of ",
               "the noisy variable '", t$xname, "'. A new row has no ",
               "latent value, so its observed value is used", call. = FALSE)
    }
    v <- v * as.numeric(xv)
  }
  for (f in mt$mult_factors) v <- v * nd_mult(f)
  v
}

#' The natural-scale hyperparameters of a fit's `me()` terms, brms's
#' `meanme_<coef>`, `sdme_<coef>` and `corme__<coef1>__<coef2>`, with
#' delta-method standard errors from the joint covariance of the
#' internal components and intervals at level `prob`. An SD's interval
#' is the Wald interval of its log, mapped back, so it stays positive as
#' a random-effect SD's does. `NULL` without `me()` terms.
#'
#' @noRd
me_hyper_table <- function(fit, prob = 0.95) {
  me <- fit$frame[["me"]]
  if (is.null(me)) return(NULL)
  comps <- me_comps(fit)
  est <- fit$estimates[comps]
  lens <- lengths(est)
  at <- function(p, cp) {
    if (!cp %in% comps) return(numeric(0))
    k <- match(cp, comps)
    p[sum(lens[seq_len(k - 1L)]) + seq_len(lens[[k]])]
  }
  f <- function(p) {
    me_hyper_values(me, at(p, "meanme"), at(p, "logsdme"),
                    at(p, "thetame"))
  }
  p0 <- unlist(est, use.names = FALSE)
  val <- f(p0)
  map <- outer_par_map(fit)
  V <- tryCatch(suppressWarnings(stats::vcov(fit, full = TRUE)),
                error = function(e) NULL)
  se <- rep(NA_real_, length(val))
  pos <- which(map$comp %in% comps)
  if (!is.null(V) && length(pos) == length(p0)) {
    h <- 1e-6
    J <- vapply(seq_along(p0), function(j) {
      up <- p0
      dn <- p0
      up[j] <- up[j] + h
      dn[j] <- dn[j] - h
      (f(up) - f(dn)) / (2 * h)
    }, numeric(length(val)))
    J <- matrix(J, nrow = length(val))
    Vh <- V[pos, pos, drop = FALSE]
    se <- sqrt(pmax(0, diag(J %*% Vh %*% t(J))))
  }
  z <- stats::qnorm(1 - (1 - prob) / 2)
  lo <- val - z * se
  hi <- val + z * se
  sd_row <- startsWith(names(val), "sdme_")
  lo[sd_row] <- val[sd_row] * exp(-z * se[sd_row] / val[sd_row])
  hi[sd_row] <- val[sd_row] * exp(z * se[sd_row] / val[sd_row])
  m <- cbind(Estimate = val, Est.Error = se, lo, hi)
  colnames(m)[3:4] <- brms_ci_cols(prob)
  rownames(m) <- names(val)
  m
}

#' The parameter-template components that hold the `me()`
#' hyperparameters, in template order.
#'
#' @noRd
me_comps <- function(fit) {
  intersect(c("meanme", "logsdme", "thetame"),
            names(fit$frame[["par_template"]]))
}

#' brms's names and values of the `me()` hyperparameters at one
#' parameter vector: `meanme_<coef>` and `sdme_<coef>` per call, then
#' `corme__<coef1>__<coef2>` (`corme_<gr>__...` for a grouped term) per
#' correlated pair, which is brms's `rename_Xme()` spelling.
#'
#' @noRd
me_hyper_values <- function(me, meanme, logsdme, thetame) {
  coefs <- vapply(me$terms, `[[`, "", "coef")
  out <- c(stats::setNames(meanme, paste0("meanme_", coefs)),
           stats::setNames(exp(logsdme), paste0("sdme_", coefs)))
  for (g in me$groups) {
    if (!g$cor) next
    M <- length(g$K)
    C <- us_chol_cor(thetame[g$th_idx], M)
    for (a in seq_len(M - 1L)) {
      for (b in seq(a + 1L, M)) {
        out[[paste0("corme", if (nzchar(g$gr)) paste0("_", g$gr), "__",
                    coefs[g$K[a]], "__", coefs[g$K[b]])]] <- C[a, b]
      }
    }
  }
  out
}

#' Put the `me()` hyperparameters into a hypothesis environment, see
#' `hyp_env_vals()`.
#'
#' @noRd
me_put_hyper <- function(fit, vals, comp, put) {
  me <- fit$frame[["me"]]
  if (is.null(me)) return(invisible(NULL))
  v <- me_hyper_values(me, vals[comp == "meanme"], vals[comp == "logsdme"],
                       vals[comp == "thetame"])
  for (nm in names(v)) put(nm, v[[nm]])
  invisible(NULL)
}

#' Control the correlation of noise-free latent variables
#'
#' Several [me()][frmtmb-me] terms model their latent (noise-free)
#' values as correlated by default, which is brms's default too. Add
#' `set_mecor(FALSE)` to a formula to model them as independent.
#'
#' @param mecor `TRUE` (the default) to estimate the correlation between
#'   the latent values of the `me()` terms that share a grouping, `FALSE`
#'   to fix it at zero.
#' @return An object of class `frmtmb_mecor`, which `+` adds to a
#'   formula built with [bf()] or [mvbf()].
#' @seealso [frmtmb-me] for the model that `me()` specifies.
#' @examples
#' set.seed(3)
#' n <- 120
#' tx <- rnorm(n)
#' tz <- 0.6 * tx + rnorm(n, 0, 0.8)
#' d <- data.frame(x = tx + rnorm(n, 0, 0.3), z = tz + rnorm(n, 0, 0.3),
#'                 sx = 0.3, sz = 0.3)
#' d$y <- 1 + 0.5 * tx - 0.4 * tz + rnorm(n, 0, 0.5)
#'
#' # independent latent values
#' f0 <- bf(y ~ me(x, sx) + me(z, sz)) + set_mecor(FALSE)
#' fit0 <- frm(f0, data = d)
#' @export
set_mecor <- function(mecor = TRUE) {
  check_flag(mecor, "mecor")
  structure(list(mecor = mecor), class = "frmtmb_mecor")
}

#' Noise-free predictors with known measurement error
#'
#' `me(x, sdx)` in a formula says that the column `x` is a noisy
#' measurement of a latent, noise-free value, with known measurement
#' standard deviation `sdx`, and that the model uses the latent value.
#' It is brms's `me()` term with brms's meaning.
#'
#' The model is
#' \deqn{x_i \sim N(\tilde{x}_i, sdx_i), \qquad
#'   \tilde{x}_i \sim N(meanme, sdme)}
#'
#' where `meanme` and `sdme` are estimated. The coefficient of the term
#' is `bsp_me<x><sdx>` in brms's names, for example `bsp_mexsx` for
#' `me(x, sx)`.
#'
#' `me()` is not a function. It is a formula special that the parser
#' reads, like `mo()` and `mi()`.
#'
#' @section Arguments of `me()`:
#' * `x`: the noisy variable, a numeric column or an expression of
#'   columns, such as `log(x)`.
#' * `sdx`: the known measurement SD, a column, an expression or a
#'   positive constant.
#' * `gr`: optional. The name of a grouping column. There is then one
#'   latent value per level of `gr`, and `x` and `sdx` must be constant
#'   within each level.
#'
#' @section Estimation:
#' The latent values are integrated out by the Laplace approximation,
#' together with the random effects. For a gaussian response with a
#' linear `me()` term the marginal likelihood is multivariate normal,
#' and the approximation is exact. An interaction of two `me()` terms,
#' or a non-gaussian response, makes it approximate.
#'
#' Several `me()` terms that share a grouping have correlated latent
#' values by default, as in brms. Use [set_mecor()] to make them
#' independent.
#'
#' `me()` can appear on its own, in a `:` or `*` interaction with
#' numeric variables or with other `me()` terms, and in the formula of
#' any distributional or nonlinear parameter. It cannot appear in a
#' group-level term, in a nonlinear formula body, inside another
#' function call, or in the same interaction as `mo()` or `mi()`.
#'
#' @section Parameters and methods:
#' The hyperparameters are the template components `meanme`, `logsdme`
#' (the log of `sdme`) and `thetame` (the correlation parameters). The
#' Noise-free Terms section of [summary()] reports them on brms's scale
#' as `meanme_<coef>`, `sdme_<coef>` and `corme__<coef1>__<coef2>`.
#' Priors on brms's classes `meanme`, `sdme` and `corme` are refused;
#' class `"b"` addresses the `me()` coefficients by their names, for
#' example `coef = "mexsx"`.
#'
#' [anova()] compares only fits with the same `me()` calls, because each
#' call puts its noisy variable's measurements into the likelihood.
#' [AIC()] does not check this, so compare by AIC only across fits that
#' share their `me()` calls.
#'
#' [fitted()], [predict()], [residuals()] and [simulate()] on the fitted
#' data use the estimated latent values (the conditional modes). On new
#' data the observed value of the noisy variable stands in for the latent
#' value. brms draws that value from `N(x, sdx)`, so the two agree on
#' the mean of a linear term, and frmtmb does not add the measurement
#' noise to a prediction interval.
#'
#' @name frmtmb-me
#' @aliases me
#' @seealso [set_mecor()]; `mi()` for a measurement model that has
#'   covariates of its own, `bf(x | mi(sdx) ~ z)`.
#' @examples
#' set.seed(1)
#' n <- 200
#' tx <- rnorm(n, 1, 0.8)
#' d <- data.frame(x = tx + rnorm(n, 0, 0.4), sx = 0.4)
#' d$y <- 2 + 0.7 * tx + rnorm(n, 0, 0.5)
#'
#' fit <- frm(bf(y ~ me(x, sx)) + gaussian(), data = d)
#' fixef(fit)
#'
#' # the naive slope is attenuated by the measurement error
#' coef(lm(y ~ x, data = d))[["x"]]
NULL
