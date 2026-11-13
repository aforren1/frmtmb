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

#' brms's check that every row of a `me(gr = )` level carries one measured
#' value and one SD, at the fit and on new data.
#'
#' @noRd
me_check_level_values <- function(xv, sv, J, levels, cl) {
  for (l in seq_along(levels)) {
    take <- J == l
    # a missing value is that row's own NA, and says nothing of the level
    if (length(unique(stats::na.omit(xv[take]))) > 1L ||
          length(unique(stats::na.omit(sv[take]))) > 1L) {
      frm_stop("Measured values and measurement error should be ",
               "unique for each group. Occured for level '", levels[l],
               "' of group '", cl$gr, "' in ", cl$key, call. = FALSE)
    }
  }
  invisible(NULL)
}

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
  # a missing sdx is refused where the data are read (me_build_frame()),
  # as brms refuses it in standata() and answers default_prior()
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
me_build_frame <- function(spec, mf, env, n, n_miss, check_noise = TRUE) {
  calls <- me_registry(spec)
  if (!length(calls)) return(NULL)
  mecor <- !isFALSE(spec[["mecor"]])
  terms <- list()
  for (k in seq_along(calls)) {
    cl <- calls[[k]]
    if (is.null(cl$sdx)) {
      if (check_noise) {
        # brms's own words (get_me_noise())
        frm_stop("Argument 'sdx' is missing in function 'me': ", cl$key,
                 ". Give the known measurement SD, me(x, sdx)",
                 call. = FALSE)
      }
      # a prior table reads no noise, as brms's default_prior() does
      # not; the placeholder keeps the latent layout
      cl$sdx <- 1
    }
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
      me_check_level_values(xv, sv, J, levels, cl)
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

#' `me()` column values in new data.
#'
#' A row the fit did not see has no latent value. `fit$me_draw`, an
#' environment one replicate or posterior draw holds, asks for brms's answer
#' (prepare_predictions(), "sample new values of latent variables"): a new
#' latent value per draw from `normal(x, sdx)`, one per row, or one per level
#' under `gr =`, independent between calls whatever `corme` says. Without it
#' the observed value stands in for the latent one, which is the mean of that
#' draw and what a summary at the estimate reads.
#'
#' A latent value is one number however many columns of the predictor read
#' it (`me(x, sdx) + me(x, sdx):z`), however many times one draw evaluates the
#' predictor, and whichever response or dpar reads it, hence the store, keyed
#' by the call. brms draws anew for each predictor, against its own Stan
#' program, which holds one `Xme_1`: at the training rows the observed
#' correlation of two responses' residuals is 0.734, one shared latent predicts
#' 0.739 and brms's replicates give -0.001 (dev/review-predread-mejoint.R, seed
#' 7701). So this is a stated divergence.
#'
#' @noRd
me_newdata_value <- function(fit, lp, mt, newdata, env, nd_mult) {
  terms <- fit$frame[["me"]]$terms
  v <- 1
  for (k in mt$terms) {
    v <- v * me_newdata_latent(fit, lp, terms[[k]], k, newdata, env)
  }
  for (f in mt$mult_factors) v <- v * nd_mult(f)
  v
}

#' The observed values of one `me()` call on the rows of `newdata`.
#'
#' @noRd
me_newdata_x <- function(t, newdata, env) {
  xv <- tryCatch(eval(t$x, newdata, env), error = function(e) NULL)
  if (is.null(xv) || !(is.numeric(xv) || (is.logical(xv) && all(is.na(xv))))) {
    frm_stop(t$key, ": newdata must supply numeric values of ",
             "the noisy variable '", t$xname, "'. A new row has no ",
             "latent value, so it is read from the observed one",
             call. = FALSE)
  }
  if (anyNA(xv)) me_warn_missing(t, t$xname, sum(is.na(xv)))
  as.numeric(xv)
}

#' Say once per reader's call that rows of `newdata` read `NA`.
#'
#' A row with no `x` or `sdx` reads `NA`, as a plain predictor's does, and the
#' fit's own data can hold such rows: refusing stopped `fitted(fit, newdata = )`
#' on the 198 rows that have values (dev/review-predread-r1-check5.R). brms
#' reads `NaN` and warns.
#'
#' One reader builds the design many times, once per finite-difference step and
#' per draw, so the mark sits in the frame of the outermost function of this
#' family of packages, which ends with the reader's call and needs no state in
#' the session.
#'
#' @noRd
me_warn_missing <- function(t, what, n) {
  top <- NULL
  for (i in seq_len(sys.nframe())) {
    ns <- topenv(environment(sys.function(i)))
    if (isNamespace(ns) && startsWith(getNamespaceName(ns), "frmtmb")) {
      top <- sys.frame(i)
      break
    }
  }
  key <- paste(t$key, what)
  seen <- if (!is.null(top)) get0(".me_missing_said", top, inherits = FALSE)
  if (key %in% seen) return(invisible(NULL))
  if (!is.null(top)) assign(".me_missing_said", c(seen, key), envir = top)
  frm_warning(t$key, ": ", n, " row(s) of newdata have no value of '", what,
              "', so they read NA", call. = FALSE,
              class = "frmtmb_me_missing")
  invisible(NULL)
}

#' One `me()` call's value on the rows of `newdata`: the observed value,
#' checked, or under `fit$me_draw` the drawn latent value.
#'
#' @noRd
me_newdata_latent <- function(fit, lp, t, k, newdata, env) {
  xv <- me_newdata_x(t, newdata, env)
  store <- fit[["me_draw"]]
  if (is.null(store)) {
    # a row with no SD has no draw, so no reader reads a value for it
    xv[is.na(me_newdata_observed(t, newdata, env)$sd)] <- NA_real_
    return(xv)
  }
  key <- as.character(k)
  old <- store[[key]]
  if (!is.null(old) && identical(old$x, xv)) return(old$value)
  o <- me_newdata_observed(t, newdata, env)
  val <- rep(NA_real_, length(xv))
  ok <- which(!is.na(o$x) & !is.na(o$sd))
  if (length(ok)) {
    lev <- unique(o$J[ok])
    first <- ok[match(lev, o$J[ok])]
    drawn <- stats::rnorm(length(first), o$x[first], o$sd[first])
    val[ok] <- drawn[match(o$J[ok], lev)]
  }
  store[[key]] <- list(x = xv, value = val)
  val
}

#' One `me()` call on the rows of `newdata`, with brms's checks on them
#' (`data_Xme()`, `get_me_noise()`): the observed value, its SD, and the
#' index of the row's latent value, the row's own or its level's under
#' `gr =`.
#'
#' @noRd
me_newdata_observed <- function(t, newdata, env) {
  xv <- me_newdata_x(t, newdata, env)
  sv <- tryCatch(eval(t$sdx, newdata, env), error = function(e) NULL)
  if (is.null(sv) || !(is.numeric(sv) || (is.logical(sv) && all(is.na(sv))))) {
    frm_stop(t$key, ": newdata must supply numeric values of the ",
             "measurement SD '", deparse1(t$sdx), "'. A new row's latent ",
             "value reads it", call. = FALSE)
  }
  sv <- rep_len(as.numeric(sv), length(xv))
  if (anyNA(sv)) me_warn_missing(t, deparse1(t$sdx), sum(is.na(sv)))
  if (any(sv <= 0, na.rm = TRUE)) {
    frm_stop("Measurement error should be positive: ", deparse1(t$sdx),
             " in ", t$key, " has ", sum(sv <= 0, na.rm = TRUE),
             " value(s) <= 0 in newdata", call. = FALSE)
  }
  J <- seq_along(xv)
  if (!is.null(t$gr)) {
    gv <- newdata[[t$gr]]
    if (is.null(gv)) stop_newdata_missing(t$gr)
    gc <- as.character(gv)
    # a missing level is a new level of its own, as brms makes it (#706), and
  # a display grid holds the group at NA, so each of its rows is one
    gc[is.na(gc)] <- paste0("new_", seq_len(sum(is.na(gc))), "__")
    lev <- unique(gc)
    J <- match(gc, lev)
    me_check_level_values(xv, sv, J, lev, t)
  }
  list(x = xv, sd = sv, J = J)
}

#' A `me(gr = )` level the fit did not see, refused unless allowed, as brms's
#' validate_newdata() refuses it. A missing level is a new one.
#'
#' `fit$me_new_levels` allows it alone: frmtmb.sample's draws take a new
#' level's latent value from `normal(x, sdx)` per draw, as brms does, and
#' refuse a group-level effect's new level, which they do not draw.
#'
#' @noRd
me_check_new_levels <- function(fit, lp, newdata, allow_new_levels) {
  if (isTRUE(allow_new_levels) || isTRUE(fit[["me_new_levels"]])) {
    return(invisible(NULL))
  }
  terms <- fit$frame[["me"]]$terms
  for (k in unique(unlist(lapply(lp[["me"]], `[[`, "terms")))) {
    t <- terms[[k]]
    if (is.null(t$gr) || is.null(newdata[[t$gr]])) next
    gc <- as.character(newdata[[t$gr]])
    new <- unique(gc[!gc %in% t$levels])
    if (length(new)) {
      stop_new_levels(paste0("Levels ", paste0("'", new, "'", collapse = ", "),
                             " of grouping factor '", t$gr, "' in ", t$key,
                             " cannot be found in the fitted model"),
                      "Consider setting argument 'allow_new_levels' to TRUE")
    }
  }
  invisible(NULL)
}

#' What a new row's `me()` draw adds to one predictor's variance.
#'
#' `frm_linpred(type = "link")` reads the observed value, which is exactly the
#' mean of the predictor under that draw, and leaves the spread of the draw to a
#' variance that is not coefficient uncertainty, as a new grouping level's. One
#' entry per call: `k`, `g`, `d eta / d Xme` at the observed values, the SD and
#' the row's level. A term linear in the call has exactly `(b sdx)^2`. A product
#' of calls has one more entry, `k = NA`, for the part no one call carries, so
#' the variance is exact: `(b sx sz)^2` on `me(x, sx):me(z, sz)`.
#'
#' Read through `sdx`, which `newdata` must therefore hold, as brms's
#' standata() requires it for every prediction.
#'
#' `fit$me_one_draw_by_value`, set by a difference curve, reads two rows with
#' one value and one SD, to a tolerance (`me_same_value()`), as one draw, which
#' is what a contrast at one position means; every other reader takes a draw per
#' row, as brms's draws do.
#'
#' @noRd
me_newdata_extra <- function(fit, lp, newdata, env, nd_mult, n) {
  if (!length(lp[["me"]]) || isTRUE(fit[["me_no_extra"]])) return(list())
  terms <- fit$frame[["me"]]$terms
  est <- fit$estimates
  obs <- new.env(parent = emptyenv())
  x_of <- function(k) {
    key <- as.character(k)
    if (is.null(obs[[key]])) {
      obs[[key]] <- me_newdata_observed(terms[[k]], newdata, env)
    }
    obs[[key]]
  }
  g <- list()
  rest <- numeric(n)
  rest_J <- NULL
  for (mt in lp[["me"]]) {
    b <- est[[lp[["par"]]]][lp[["idx"]][mt$col]]
    m <- 1
    for (fct in mt$mult_factors) m <- m * nd_mult(fct)
    first <- numeric(n)
    for (k in mt$terms) {
      others <- 1
      for (o in setdiff(mt$terms, k)) others <- others * x_of(o)$x
      key <- as.character(k)
      g[[key]] <- (g[[key]] %||% numeric(n)) + b * m * others
      first <- first + (others * x_of(k)$sd)^2
    }
    if (length(mt$terms) > 1L) {
      # Var(prod X_k) is prod(x^2 + s^2) - prod(x^2) for independent calls;
      # what the calls' first-order parts leave of it
      ob <- lapply(mt$terms, x_of)
      all_v <- Reduce(`*`, lapply(ob, function(o) o$x^2 + o$sd^2))
      all_x <- Reduce(`*`, lapply(ob, function(o) o$x^2))
      rest <- rest + (b * m)^2 * pmax(all_v - all_x - first, 0)
      rest_J <- ob[[1L]]$J
      if (!all(vapply(ob, function(o) identical(o$J, rest_J), NA))) {
        rest_J <- seq_len(n)
      }
    }
  }
  by_value <- isTRUE(fit[["me_one_draw_by_value"]])
  out <- lapply(names(g), function(key) {
    x <- x_of(as.integer(key))
    J <- x$J
    if (by_value && is.null(terms[[as.integer(key)]]$gr)) {
      id <- paste(me_same_value(x$x), me_same_value(x$sd))
      J <- match(id, id)
    }
    list(k = as.integer(key), g = g[[key]], sd = x$sd, J = J)
  })
  if (any(rest > 0, na.rm = TRUE)) {
    out[[length(out) + 1L]] <- list(k = NA_integer_, g = sqrt(rest),
                                    sd = rep(1, n), J = rest_J)
  }
  out
}

#' Which values of `v` are one value, as an id per position. Exact equality made
#' the standard error of a difference jump by 20 times: 0.0485 at `dx = 0`,
#' 0.9979 at `dx = 1e-12`, on a grid built by arithmetic
#' (dev/review-predread-r1-check2.R). The tolerance is `all.equal()`'s, relative
#' to the pair and to the column's scale, so a column near 0 has a scale too.
#'
#' @noRd
me_same_value <- function(v) {
  id <- rep(NA_integer_, length(v))
  ok <- which(!is.na(v))
  if (!length(ok)) return(id)
  o <- ok[order(v[ok])]
  s <- v[o]
  tol <- sqrt(.Machine$double.eps) *
    pmax(abs(s[-1L]), abs(s[-length(s)]), mean(abs(s)))
  id[o] <- cumsum(c(TRUE, diff(s) > tol))
  id
}

#' How one response's predictors load the draws of a new row, or `NULL` when
#' no predictor reads one.
#'
#' A new row reads three kinds of draw that are not coefficient uncertainty:
#' a `me()` value, a latent residual on `mu`, and the effects of a grouping
#' level the fit did not see (`allow_new_levels`), as brms draws each of them.
#' All three move the predictors linearly and are normal, so per row they are
#' one normal vector over the predictors they reach. Its covariance is
#' factored per row with the lead predictor first (`chol_rows()`), so that
#' `dirs[[1]]` alone moves the lead and each direction is an independent
#' standard normal.
#'
#' `key` is the lead predictor, `s2` its variance per row, `dirs` each
#' predictor's movement per unit of each direction. A `me()` call is one draw
#' whichever predictors read it; two calls, a latent and a level's effects
#' are independent draws. A block's effects reach every predictor that holds
#' one of its components through the same draw. Through 0.79.1 a second
#' predictor moved by its regression on the lead alone, which is exact only
#' when one call carries both: on a `hurdle_poisson()` row with a latent the
#' draw's sd read 0.733 against 0.925 (dev/readers-probe2.R).
#'
#' @noRd
me_row_loads <- function(fit, rspec, newdata, use_re, anl, dpar = NULL) {
  if (is.null(newdata)) return(NULL)
  rn <- rspec$resp_name
  fit[["me_no_extra"]] <- NULL
  sl <- latres_row_sd(fit, rn, dpar)
  dps <- if (is.null(dpar)) names(rspec$dpars) else dpar
  # a body reads its nonlinear parameters' draws: conditional_effects() asks
  # for mu alone, and its band missed a me() value in `a ~ me(x, sx)`
  dps <- me_body_parts(fit, rn, dps)
  n <- nrow(newdata)
  # one entry per independent standard normal source: its load on each key
  src <- list()
  for (d in dps) {
    key <- linpred_key(rn, d)
    lp <- fit$frame[["linpreds"]][[key]]
    if (is.null(lp) || !is.null(lp[["nl_body"]])) next
    need_nl <- isTRUE(anl) && use_re && me_reads_groups(fit, key)
    if (length(lp[["me"]]) || need_nl) {
      ed <- lp_eta_design(fit, lp, newdata, use_re, anl)
      for (e in ed[["me_parts"]]) {
        id <- if (is.na(e$k)) paste0("rest\r", key) else paste0("me\r", e$k)
        src <- me_src_add(src, id, key, e$g * e$sd, n)
      }
      if (need_nl) {
        B <- extra_var_blocks(lp_extra_var(fit, ed, use_re)$new_levels,
                              ed[["n"]])
        for (bn in names(B)) {
          L <- B[[bn]]$M %*% psd_root(B[[bn]]$S)
          for (cc in seq_len(ncol(L))) {
            if (any(L[, cc] != 0)) {
              src <- me_src_add(src, paste0("re\r", bn, "\r", cc), key,
                                L[, cc], n)
            }
          }
        }
      }
    }
    if (identical(d, "mu") && !is.null(sl)) {
      src <- me_src_add(src, "latent", key, rep_len(sl, n), n)
    }
  }
  if (!length(src)) return(NULL)
  used <- unique(unlist(lapply(src, names), use.names = FALSE))
  keys <- linpred_key(rn, dps)
  keys <- keys[keys %in% used]
  D <- length(keys)
  # the per-row covariance over the keys, its lower triangle
  C <- lapply(seq_len(D), function(i) lapply(seq_len(i), function(j) {
    numeric(n)
  }))
  for (s in src) {
    for (i in seq_len(D)) {
      li <- s[[keys[i]]]
      if (is.null(li)) next
      for (j in seq_len(i)) {
        lj <- s[[keys[j]]]
        if (!is.null(lj)) C[[i]][[j]] <- C[[i]][[j]] + li * lj
      }
    }
  }
  L <- chol_rows(C, D)
  dirs <- list()
  for (j in seq_len(D)) {
    col <- lapply(seq_len(D), function(i) {
      if (i >= j) L[[i]][[j]] else numeric(n)
    })
    if (!any(vapply(col, function(v) any(v != 0), NA))) next
    dirs[[length(dirs) + 1L]] <- stats::setNames(col, keys)
  }
  if (!length(dirs)) return(NULL)
  list(key = keys[1L], s2 = C[[1L]][[1L]], dirs = dirs)
}

#' The dpars `dps` and every parameter a body among them reads, through
#' nested bodies, in that order.
#'
#' @noRd
me_body_parts <- function(fit, rn, dps) {
  out <- character(0)
  todo <- dps
  while (length(todo)) {
    d <- todo[1L]
    todo <- todo[-1L]
    if (d %in% out) next
    out <- c(out, d)
    lp <- fit$frame[["linpreds"]][[linpred_key(rn, d)]]
    if (!is.null(lp[["nl_body"]])) {
      todo <- c(todo, lp[["nl_pars"]], lp[["nl_dpar_refs"]])
    }
  }
  out
}

#' One source's load `l` on predictor `key`, added to the list of sources.
#'
#' @noRd
me_src_add <- function(src, id, key, l, n) {
  if (is.null(src[[id]])) src[[id]] <- list()
  src[[id]][[key]] <- (src[[id]][[key]] %||% numeric(n)) + l
  src
}

#' Whether a predictor holds a group-level component whose new level has a
#' marginal law: a smooth, `gp()`, `hsgp()` or latent block has none.
#'
#' @noRd
me_reads_groups <- function(fit, key) {
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    if (bk[["covstruct"]] %in% c("smooth", "gp", "hsgp", "latres")) next
    for (cp in bk[["components"]]) {
      if (identical(cp$lp_key, key)) return(TRUE)
    }
  }
  FALSE
}

#' A root `R` of a covariance `S`, `R R' = S`, that holds where `S` is
#' singular, as a correlation of 1 or an sd of 0 makes it.
#'
#' @noRd
psd_root <- function(S) {
  e <- eigen(as.matrix(S), symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
}

#' The lower Cholesky factor of a covariance given per row, `C[[i]][[j]]` for
#' `j <= i`. A pivot whose remaining variance is a rounding of what was there
#' is 0, so a predictor that earlier directions carry whole opens none.
#'
#' @noRd
chol_rows <- function(C, D) {
  L <- lapply(seq_len(D), function(i) vector("list", i))
  for (j in seq_len(D)) {
    s <- C[[j]][[j]]
    for (k in seq_len(j - 1L)) s <- s - L[[j]][[k]]^2
    # one draw read by two predictors leaves the second a share of about
    # 1e-16 of its variance here, which is no direction
    ok <- C[[j]][[j]] > 0 & s > 1e-10 * C[[j]][[j]]
    piv <- ifelse(ok, sqrt(pmax(s, 0)), 0)
    L[[j]][[j]] <- piv
    for (i in seq_len(D - j) + j) {
      v <- C[[i]][[j]]
      for (k in seq_len(j - 1L)) v <- v - L[[i]][[k]] * L[[j]][[k]]
      L[[i]][[j]] <- ifelse(ok, v / ifelse(ok, piv, 1), 0)
    }
  }
  L
}

#' The marginal sd a latent residual adds to `mu` of a new row, or `NULL` where
#' the read is not of `mu` or the fit has none.
#'
#' @noRd
latres_row_sd <- function(fit, rn, dpar) {
  if (!is.null(dpar) && !identical(dpar, "mu")) return(NULL)
  bk <- frame_latres(fit$frame)[[rn]]
  if (is.null(bk) || isTRUE(fit$frame[["latres_off"]])) return(NULL)
  latres_marginal_sd(fit$estimates[["theta"]][bk[["theta_idx"]]],
                     bk[["aux"]])
}

#' The value of `fit` at a point of the draws: `what(fit)` with every
#' predictor of `loads` moved. `z` is a list with one entry per direction,
#' each one number or one per row; a bare `z` is the first direction's.
#' `d` moves the lead predictor on top.
#'
#' @noRd
me_node_value <- function(fit, loads, what) {
  function(z, d = 0) {
    if (!is.list(z)) z <- list(z)
    sh <- fit
    moved <- list()
    for (j in seq_along(loads$dirs)) {
      zj <- if (j <= length(z)) z[[j]] else 0
      for (key in names(loads$dirs[[j]])) {
        moved[[key]] <- (moved[[key]] %||% 0) + loads$dirs[[j]][[key]] * zj
      }
    }
    moved[[loads$key]] <- moved[[loads$key]] + d
    sh[["eta_shift"]] <- list(by_key = moved)
    what(sh)
  }
}

# Gauss-Hermite nodes per direction of a new row's draws, by their number: a
# node costs one prediction, about 1 ms, so a product grid keeps to a few
# hundred. The slope's integral scales the coefficients' share alone and
# takes fewer.
me_gh_size <- c(41L, 21L, 11L, 7L)
me_gh_slope_size <- c(41L, 9L, 5L, 3L)

#' The mean and second moment of `val` over `r` independent standard normals,
#' on a product Gauss-Hermite grid of `K` nodes each. A node whose value is
#' not a finite number says nothing, and a cell that lost real weight reads
#' `NaN`, as `response_mean_gh()` reads it.
#'
#' @noRd
me_gh_moments <- function(val, r, K, second = TRUE) {
  gh <- gh_std_normal(K)
  idx <- as.matrix(expand.grid(rep(list(seq_len(K)), r)))
  w <- apply(idx, 1L, function(i) prod(gh$w[i]))
  # a product weight below 1e-16 of the largest moves no digit of the mean
  keep <- w > 1e-16 * max(w)
  idx <- idx[keep, , drop = FALSE]
  w <- w[keep]
  t1 <- 0
  t2 <- 0
  lost <- 0
  for (k in seq_along(w)) {
    v <- val(as.list(gh$x[idx[k, ]]))
    bad <- !is.finite(v)
    if (any(bad)) {
      v[bad] <- 0
      lost <- lost + w[k] * bad
    }
    t1 <- t1 + w[k] * v
    if (second) t2 <- t2 + w[k] * v^2
  }
  m1 <- t1 / (1 - lost)
  m2 <- t2 / (1 - lost)
  m1[lost >= latres_gh_lost] <- NaN
  m2[lost >= latres_gh_lost] <- NaN
  list(m1 = m1, m2 = m2)
}

#' Whether the exact mean over a new row's draws is infinite, decided from
#' the model and not from the numbers: a lognormal mean is
#' `exp(mu + sigma^2 / 2)`, and a normal draw on `sigma`'s log link gives
#' `E[exp(c exp(2 u))]`, which diverges at any sd of `u`. That holds for a
#' lognormal component's `sigma<k>` of a mixture, and for a draw on a
#' nonlinear parameter that a body for `sigma` reads: the first form of this
#' rule missed both, and read 2.5e207 and 4.7e82 with no warning. Every other
#' mean this reads is finite: a log-linear one is `exp(eta + v / 2)`.
#'
#' @noRd
me_mean_diverges <- function(fit, rspec, loads, dpar) {
  if (!is.null(dpar)) return(FALSE)
  sig <- me_lognormal_sigmas(rspec$family)
  if (!length(sig)) return(FALSE)
  reach <- me_body_parts(fit, rspec$resp_name, sig)
  moved <- unique(unlist(lapply(loads$dirs, names), use.names = FALSE))
  any(vapply(moved, function(k) {
    isTRUE(fit$frame[["linpreds"]][[k]][["dpar"]] %in% reach)
  }, NA))
}

#' The `sigma` dpars of a family's lognormal parts: `sigma` of a lognormal,
#' shifted or hurdle lognormal, and `sigma<k>` of a mixture whose component
#' `k` is one, read off the family's name, `mixture(lognormal, gaussian)`.
#'
#' @noRd
me_lognormal_sigmas <- function(fam) {
  nm <- fam[["family"]] %||% ""
  if (is.null(fam[["mix"]])) {
    return(if (grepl("lognormal", nm, fixed = TRUE)) "sigma" else character())
  }
  inner <- sub("^[^(]*[(](.*)[)]$", "\\1", nm)
  comps <- trimws(strsplit(inner, ",", fixed = TRUE)[[1L]])
  k <- which(grepl("lognormal", comps, fixed = TRUE))
  if (!length(k)) character() else paste0("sigma", k)
}

#' Whether the value read is a constant times `exp(eta)` of the one predictor
#' the draws move, so that its mean over them is `exp(v / 2)` times the value
#' at the observed values, `v` the draws' variance there, and its second
#' moment `exp(2 v)` times the square: a log link whose inverse is the mean
#' (poisson, negbinomial, a log-link Gamma), a lognormal's `mu`, or a dpar
#' read on a log link. A truncated mean, a lag SAR term and a latent that
#' joins after the link are not.
#'
#' @noRd
me_log_linear <- function(fit, rspec, loads, dpar) {
  moved <- unique(unlist(lapply(loads$dirs, names), use.names = FALSE))
  if (length(moved) != 1L) return(FALSE)
  rn <- rspec$resp_name
  read <- dpar %||% "mu"
  if (!identical(moved, linpred_key(rn, read))) return(FALSE)
  lp <- fit$frame[["linpreds"]][[moved]]
  if (is.null(lp) || !is.null(lp[["nl_body"]])) return(FALSE)
  fam <- rspec$family
  if (!is.null(dpar)) {
    return(is.null(dpar_report_hook(fam, dpar, rspec)) &&
             identical(dpar_response_link(rspec, lp, dpar)[["name"]], "log"))
  }
  if (has_trunc(rspec) || autocor_is_seam(fit$frame[["autocor"]][[rn]])) {
    return(FALSE)
  }
  if (isTRUE(frame_latres(fit$frame)[[rn]][["aux"]][["on_response"]])) {
    return(FALSE)
  }
  link <- lp[["link"]][["name"]] %||% ""
  (mean_is_mu(fam) && identical(link, "log")) ||
    (identical(fam[["family"]], "lognormal") && identical(link, "identity"))
}

# The steps me_keep_dirs() moves a direction by, in sds, irrational and
# irrationally related, so no period of a bounded body meets all of them.
me_keep_steps <- c(sqrt(2) / 2, (1 + sqrt(5)) / 2 / 1.3, sqrt(3))

#' The draws that move what is read: a direction whose unit step changes no
#' cell is left out, so it costs no nodes, and past `length(me_gh_size)`
#' directions the ones that move it most are kept. Through the first round
#' of this rule a mixture kept its first four by dpar order, two of which
#' moved nothing: the mean read 1.3395 against 1.0736.
#'
#' @noRd
me_keep_dirs <- function(loads, val, v0, rn) {
  r <- length(loads$dirs)
  if (r < 2L) return(loads)
  # steps of sqrt(2) / 2, the golden ratio / 1.3 and sqrt(3) sds each way, no
  # two of them rational multiples of each other: a rational step can land on
  # a period, as sin(a)^2 at a = 0 with an sd of pi or 10 pi, where steps of
  # 0.5, 1 and 1.7 all read sin(0)^2 and the draw, which doubles the mean,
  # was dropped
  infl <- vapply(seq_len(r), function(j) {
    d <- 0
    for (st in c(-1, 1) %o% me_keep_steps) {
      z <- as.list(numeric(r))
      z[[j]] <- st
      e <- abs(as.matrix(val(z)) - v0) / pmax(abs(v0), .Machine$double.xmin)
      e[!is.finite(e)] <- 0
      d <- max(d, e)
    }
    d
  }, 0)
  keep <- which(infl > 1e-12)
  if (!length(keep)) keep <- 1L
  cap <- length(me_gh_size)
  if (length(keep) > cap) {
    frm_warning("A new row of '", rn, "' reads ", length(keep),
                " independent draws that move its value, and the mean is ",
                "integrated over the ", cap, " that move it most. The rest ",
                "stay at 0 in the mean and in Est.Error",
                call. = FALSE, class = "frmtmb_draw_dims")
    keep <- sort(keep[order(-infl[keep])][seq_len(cap)])
  }
  loads$dirs <- loads$dirs[keep]
  loads
}

#' The value over several directions of draws, read along one: per row, the
#' direction in which the value moves most at the observed values, and the
#' variance along it. A band is the draw's quantiles along it, widened to the
#' whole variance; with one direction it is the draw itself.
#'
#' @noRd
me_band_dir <- function(val, r, K = 41L) {
  hz <- 1e-3
  g <- lapply(seq_len(r), function(j) {
    up <- as.list(numeric(r))
    dn <- up
    up[[j]] <- hz
    dn[[j]] <- -hz
    as.matrix((val(up) - val(dn)) / (2 * hz))
  })
  # a category display reads the direction of the category that moves most
  nrm <- Reduce(`+`, lapply(g, function(x) x^2))
  col <- max.col(nrm, ties.method = "first")
  rows <- seq_len(nrow(nrm))
  gc <- lapply(g, function(x) x[cbind(rows, col)])
  len <- sqrt(Reduce(`+`, lapply(gc, function(x) x^2)))
  a <- lapply(seq_len(r), function(j) {
    ifelse(len > 0, gc[[j]] / ifelse(len > 0, len, 1), as.numeric(j == 1L))
  })
  val1 <- function(t, d = 0) {
    if (is.list(t)) t <- t[[1L]]
    val(lapply(a, function(aj) aj * t), d)
  }
  mo <- me_gh_moments(val1, 1L, K)
  list(val = val1, var = pmax(mo$m2 - mo$m1^2, 0))
}

#' The warning that a new row's exact mean over its draws is infinite, said
#' once per reader's call with the largest change two fewer nodes make.
#'
#' @noRd
me_warn_not_finite <- function(rn, moved) {
  if (me_said(paste("not_finite", rn))) return(invisible(NULL))
  frm_warning("The exact expected response of '", rn, "' at a new row is ",
              "infinite: a lognormal's mean exp(mu + sigma^2 / 2) over a ",
              "normal draw on sigma (a me() value, a latent or a new level) ",
              "diverges in the draw's far tail. Estimate and Est.Error are ",
              "the ", me_gh_div, "-node Gauss-Hermite mean and sd, which ",
              "reach ", format(max(gh_std_normal(me_gh_div)$x), digits = 3),
              " sds of the draw, as a sampler's draws do. On the rows of ",
              "this call ", me_gh_div - 2L, " nodes move that mean by at ",
              "most ", format(moved, digits = 2), " of itself: where that is ",
              "not small, no mean of draws is stable either, and it grows ",
              "with their number. The Q columns are the draw's quantiles",
              call. = FALSE, class = "frmtmb_mean_not_finite")
}

#' The warning of the net under the rule: a mean whose rule did not settle.
#'
#' @noRd
me_warn_unsettled <- function(rn, before, after) {
  if (me_said(paste("unsettled", rn))) return(invisible(NULL))
  frm_warning("The expected response of '", rn, "' at a new row does not ",
              "settle: two fewer Gauss-Hermite nodes move its mean by ",
              format(before, digits = 2), " of itself, or it is not finite, ",
              "so its exact mean over the row's draws may be infinite. ",
              "Estimate and Est.Error are the ", me_gh_div, "-node mean and ",
              "sd, which reach ",
              format(max(gh_std_normal(me_gh_div)$x), digits = 3), " sds of ",
              "the draw, as a sampler's draws do; ", me_gh_div - 2L,
              " nodes move that by ", format(after, digits = 2), ". The Q ",
              "columns are the draw's quantiles",
              call. = FALSE,
              class = c("frmtmb_mean_unsettled", "frmtmb_mean_not_finite"))
}

#' The largest relative change of a mean between `K` and `K - 2` nodes over
#' the rows whose value at the draws' mean is finite; `Inf` where the mean
#' there is not finite, or the smaller rule's is not.
#'
#' @noRd
me_rule_change <- function(val, r, K, m1, v0) {
  if (any(is.finite(as.matrix(v0)) & !is.finite(as.matrix(m1)))) return(Inf)
  me_change(me_gh_moments(val, r, K - 2L, second = FALSE)$m1, m1, v0)
}

#' The largest relative change from `m1` to `m2` over the rows whose value at
#' the draws' mean is finite; `Inf` where either is not finite there.
#'
#' @noRd
me_change <- function(m2, m1, v0) {
  ok <- is.finite(as.matrix(v0))
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)
  if (any(ok & !(is.finite(m1) & is.finite(m2)))) return(Inf)
  ch <- abs(m2 / m1 - 1)[ok & m1 != 0]
  max(c(0, ch))
}

#' The relative change from `m1` to `m2` per row, over the cells of the row:
#' `Inf` where either is not finite while the value at the draws' mean is,
#' 0 on a row whose value there is not finite.
#'
#' @noRd
me_row_change <- function(m2, m1, v0) {
  ok <- is.finite(as.matrix(v0))
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)
  ch <- abs(m2 / m1 - 1)
  ch[m1 == 0 & m2 == 0] <- 0
  ch[!(is.finite(m1) & is.finite(m2))] <- Inf
  ch[!ok] <- 0
  apply(ch, 1L, max)
}

#' The moments `mo` with the rows `rows` taken from `mb`.
#'
#' @noRd
me_rows_take <- function(mo, mb, rows) {
  for (nm in c("m1", "m2")) {
    a <- as.matrix(mo[[nm]])
    b <- as.matrix(mb[[nm]])
    if (identical(dim(a), dim(b))) {
      a[rows, ] <- b[rows, ]
      mo[[nm]] <- a
    }
  }
  mo
}

# The net's bound on me_rule_change(): a hurdle mean on 21 nodes at a link
# sd of 3 settles to 1e-4, an exploding one by orders of magnitude.
me_net_tol <- 1e-3

#' Whether `key` was said already in the reader's call this runs under, and
#' mark it said: the mark sits in the frame of the outermost function of
#' this family of packages, as `me_warn_missing()` places it, so one call of
#' `conditional_effects()` over several grids says it once.
#'
#' @noRd
me_said <- function(key) {
  top <- NULL
  for (i in seq_len(sys.nframe())) {
    ns <- topenv(environment(sys.function(i)))
    if (isNamespace(ns) && startsWith(getNamespaceName(ns), "frmtmb")) {
      top <- sys.frame(i)
      break
    }
  }
  if (is.null(top)) return(FALSE)
  seen <- get0(".me_said", top, inherits = FALSE)
  if (key %in% seen) return(TRUE)
  assign(".me_said", c(seen, key), envir = top)
  FALSE
}

# The rule for a mean whose exact value is infinite: 11 nodes reach 5.19
# sds, about what 4000 to a million draws reach. The 41-node rule's outer
# node, at 11.6 sds, adds about 1e16 on a lognormal with sigma ~ (1 | g) at
# sd 0.27, where a mean of 4e6 draws reads 2.09.
me_gh_div <- 11L

#' The relative change of the mean and of the sd between the rule and the
#' rule of two fewer nodes, the error estimate a wide draw needs past two
#' directions: at three directions and a link sd of 3 the 11-node rule's sd
#' is a third low.
#'
#' @noRd
me_rule_check <- function(val, r, K, mo, rn) {
  lo <- me_gh_moments(val, r, K - 2L)
  sd_k <- sqrt(pmax(mo$m2 - mo$m1^2, 0))
  sd_l <- sqrt(pmax(lo$m2 - lo$m1^2, 0))
  dm <- abs(lo$m1 / mo$m1 - 1)
  ds <- abs(sd_l / sd_k - 1)
  dm[!is.finite(dm)] <- 0
  ds[!is.finite(ds)] <- 0
  if (max(dm) > me_rule_tol[1L] || max(ds) > me_rule_tol[2L]) {
    frm_warning("The mean over a new row of '", rn, "' is integrated over ",
                r, " directions of draws on ", K, " nodes each, and the ",
                "rule does not settle: ", K - 2L, " nodes move the mean by ",
                format(max(dm), digits = 2), " and the sd by ",
                format(max(ds), digits = 2), " of themselves. Est.Error and ",
                "the Q columns are less exact than that",
                call. = FALSE, class = "frmtmb_draw_rule")
  }
  invisible(NULL)
}

# What me_rule_check() lets pass: a relative change of the mean, of the sd.
me_rule_tol <- c(1e-3, 0.02)

#' The value of a new row's mean at probability `p` of its `me()` draw and its
#' coefficients: the mean at that quantile of the predictor. It is the draw's
#' own
#' quantile, whatever the family's mean is, so it holds on a mean that no link
#' inverse gives (`hurdle`, `lognormal`), and it stays in the support. `se_coef`
#' is the coefficient part of the standard error on the response scale at the
#' observed value.
#'
#' @noRd
me_band_end <- function(mm, se_coef, p) {
  # the coefficients' share against the draw's own variance on the response
  # scale: the first predictor's gradient alone understates the draw when others
  # read the term, up to 3.3 times too wide on me() in mu and zi
  # (dev/review-predread-r1-band.R)
  extra <- ifelse(mm$var > 0, as.numeric(se_coef)^2 / mm$var, 0)
  if (!is.null(mm$var1)) {
    # several draws read along one direction: the others' variance widens it
    # as the coefficients' does
    extra <- ifelse(mm$var1 > 0,
                    (as.numeric(se_coef)^2 + mm$var - mm$var1) / mm$var1, 0)
  }
  # a mean the draws do not hold has no share to take: the draw's own quantile
  extra[!is.finite(extra)] <- 0
  # the direction of the draw as a whole: mm$grad moves the first predictor
  # alone, and a zi that rises with the same term can turn the mean the other
  # way, which swapped the two ends
  up <- drop(mm$val(1)) >= drop(mm$val(-1))
  z <- stats::qnorm(p) * sqrt(1 + extra) * ifelse(up, 1, -1)
  v <- drop(mm$val(z))
  v[!is.finite(v)] <- NA_real_
  v
}

#' Whether a new row's `me()` variance is read at all: a `fit` with this mark
#' reports the coefficient part alone, so a mean can add the integral's.
#'
#' @noRd
me_coef_only <- function(fit) {
  fit[["me_no_extra"]] <- TRUE
  fit
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
  se_flush_deferred(fit)
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
                error = keep_missing(NULL))
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
#' correlated pair, which is brms's `rename_Xme()` spelling. The calls
#' are listed group by group, brms's `meanme_1`, `meanme_2` order, and
#' not in the order the formula writes them: on `me(xg, gr = g) + me(z)
#' + me(wg, gr = g)` brms has `meanme_mexg meanme_mewg meanme_mez`.
#'
#' @noRd
me_hyper_values <- function(me, meanme, logsdme, thetame) {
  coefs <- vapply(me$terms, `[[`, "", "coef")
  ord <- unlist(lapply(me$groups, `[[`, "K"), use.names = FALSE)
  if (length(ord) != length(coefs)) ord <- seq_along(coefs)
  out <- c(stats::setNames(meanme[ord], paste0("meanme_", coefs[ord])),
           stats::setNames(exp(logsdme[ord]), paste0("sdme_", coefs[ord])))
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
#'
#' # correlated latent values, with an lkj() prior on the correlation
#' fit1 <- frm(bf(y ~ me(x, sx) + me(z, sz)), data = d,
#'             prior = set_prior("lkj(2)", class = "corme"))
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
#' Class `"b"` addresses the `me()` coefficients by their names, for
#' example `coef = "mexsx"`.
#'
#' @section Priors:
#' The hyperparameters have brms's prior classes. [default_prior()]
#' lists brms's rows for them:
#' * `"meanme"`: the latent mean, one row per `me()` variable under a
#'   class row (`coef = "mex"` for `me(x, sdx)`).
#' * `"sdme"`: the latent standard deviation itself, bounded below by 0,
#'   with the same rows. The density carries the log-Jacobian of the
#'   log map, as class `"sd"` does.
#' * `"corme"`: `lkj(eta)` on the correlation matrix of the terms that
#'   share a `gr =`, one class row and, for a grouped set, a row with
#'   `group =` its grouping variable.
#'
#' On [frm()] a prior on them makes the fit a MAP estimate, and the
#' defaults apply nothing. `frmtmb.sample::frm_sample()` applies brms's
#' defaults: flat on `meanme` and on `sdme`, and `lkj(1)` on `corme`.
#' The three classes take no `resp`, `dpar` or `nlpar`, because one set
#' of latent values serves the whole model, as in brms.
#'
#' [anova()] compares only fits with the same `me()` calls, because each
#' call puts its noisy variable's measurements into the likelihood.
#' [AIC()] does not check this, so compare by AIC only across fits that
#' share their `me()` calls.
#'
#' [fitted()], [predict()], [residuals()] and [simulate()] on the fitted data
#' use the estimated latent values (the conditional modes). A new row has no
#' latent value, and brms draws one from `N(x, sdx)` per draw. So do [predict()]
#' and [simulate()], per replicate. [fitted()] reports the mean of the expected
#' response over the draw, with its variance in `Est.Error` and its quantiles in
#' the `Q` columns. [frm_linpred()] stays at the observed value and adds the
#' variance to `se.fit`. New data must hold `sdx`, as in brms. Two responses or
#' two dpars read one latent value, a stated divergence from brms, which draws
#' one per predictor. Thus a category probability moves, as it is a joint
#' function of both predictors: 0.3336 against brms's 0.3106. A row with no `x`
#' or `sdx` reads `NA`, with a warning of class `frmtmb_me_missing`. brms reads
#' `NaN` and warns.
#'
#' @return `me()` is a formula term, not a free-standing function:
#'   `bf()` reads it at parse time, and the value it contributes is the
#'   noise-free predictor's coefficient together with the latent values
#'   and the measurement-model hyperparameters, reachable through
#'   [fixef()] and the Noise-free Terms section of [summary()]. This
#'   page itself documents the term grammar and returns nothing.
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
#'
#' # brms's prior rows for the latent mean and sd
#' default_prior(bf(y ~ me(x, sx)) + gaussian(), data = d)
#' fitp <- frm(bf(y ~ me(x, sx)) + gaussian(), data = d,
#'             prior = set_prior("normal(1, 1)", class = "meanme") +
#'               set_prior("exponential(1)", class = "sdme"))
#' prior_summary(fitp)
NULL
