# brms's `sample_prior`: draws from the priors beside the posterior
# draws ("yes"), and a run with the likelihood left out ("only").
#
# HOW BRMS DOES IT, read from brms 2.23.0's `stan_rngprior()` and
# `rename_pars()` and checked against fitted models
# (dev/priordraws-findings.md, section 1). Each prior STATEMENT of the
# Stan program gets one draw per iteration in generated quantities,
# independent of the likelihood: a class-wide prior on a vector of
# coefficients is one statement and one column (`prior_b`), a vector
# with a coefficient-specific prior anywhere in it is one statement per
# coefficient (`prior_b_x`), a flat prior is no statement and no column.
# A parameter with a bound is drawn from the truncated prior (brms
# rejects and redraws). A correlation matrix gets one off-diagonal
# element of an LKJ draw, a simplex a Dirichlet draw. The intercept
# drawn is the centered one, `prior_Intercept`.
#
# WHY DIRECT DRAWS HERE. brms's draws are independent draws from each
# prior, made by its random number generators. The same draws are made
# here in R, from the same densities, which gives independent draws at
# no cost. Running the sampler on the prior instead would give
# autocorrelated draws of the same thing for the price of a second run
# (dev/priordraws-findings.md, section 3, has the measurement).
#
# "only" is different: brms samples every parameter of the model from
# its prior with the likelihood left out, the group-level coefficients
# included, so its draws are a run of the sampler. frm_sample() does
# the same with core's `prior_only` objective.

#' brms's `validate_sample_prior()`: `NULL` is `"no"`, a flag is
#' `"yes"` or `"no"`, and a string is matched against the three
#' options.
#'
#' @noRd
check_sample_prior <- function(sample_prior) {
  opts <- c("no", "yes", "only")
  if (is.null(sample_prior)) return("no")
  if (!is.character(sample_prior)) {
    if (!is.logical(sample_prior) || length(sample_prior) != 1L ||
          is.na(sample_prior)) {
      frm_stop("sample_prior = must be \"no\", \"yes\" or \"only\" (or ",
               "TRUE or FALSE), not ", deparse1(sample_prior),
               call. = FALSE)
    }
    return(if (sample_prior) "yes" else "no")
  }
  if (length(sample_prior) != 1L || is.na(sample_prior) ||
        is.na(i <- pmatch(sample_prior, opts))) {
    frm_stop("sample_prior = must be \"no\", \"yes\" or \"only\", not ",
             deparse1(sample_prior), call. = FALSE)
  }
  opts[i]
}

#' The bounds brms declares on a distributional parameter without a
#' formula (`dpar_bounds()`), by its class. A prior draw of it is
#' truncated there.
#'
#' @noRd
prior_dpar_bounds <- function(cls, fit, rspec) {
  switch(cls,
    sigma = , shape = , phi = , kappa = , beta = , bs = , disc = c(0, NA),
    nu = c(1, NA),
    zi = , hu = , zoi = , coi = , bias = , quantile = c(0, 1),
    ndt = {
      y <- fit$frame[["y"]][[rspec$resp_name]]
      if (is.matrix(y)) y <- y[, 1L]
      y <- y[is.finite(y)]
      c(0, if (length(y)) min(y) else NA)
    },
    c(NA, NA))
}

#' The rows of a `validate_prior()` table that one predictor owns.
#'
#' @noRd
prior_px_rows <- function(tab, rr, px, nl) {
  tab$resp == rr &
    (if (nl) tab$nlpar == px & !nzchar(tab$dpar) else
      tab$dpar == (if (identical(px, "mu")) "" else px) &
        !nzchar(tab$nlpar))
}

#' One prior statement of the plan. `prior` is the density string or
#' `""` for a flat prior; `stan` is brms's Stan name of the parameter,
#' which brms's refusal under `sample_prior = "only"` quotes.
#'
#' @noRd
prior_stmt <- function(name, stan, prior, lb = NA, ub = NA,
                       type = "real", dim = 1L) {
  if (is.na(prior) || identical(prior, "(flat)")) prior <- ""
  list(name = name, stan = stan, prior = prior,
       lb = suppressWarnings(as.numeric(lb)),
       ub = suppressWarnings(as.numeric(ub)), type = type,
       dim = as.integer(dim))
}

#' The first row of `tab[sel, ]`, or `NULL`.
#'
#' @noRd
prior_tab_row <- function(tab, sel) {
  i <- which(sel)
  if (length(i)) tab[i[1L], , drop = FALSE]
}

#' The prior statements of a model, in the order brms's Stan program
#' writes them (`stan_predictor()`: per response, the nonlinear
#' parameters, then each distributional parameter in the family's order,
#' its population-level coefficients, intercept, monotonic simplexes
#' and coefficients, category-specific and smooth coefficients, smooth
#' standard deviations and residual correlation terms; then a mixture's
#' weights, the residual correlations and the group-level terms). Flat
#' slots are kept, with `prior = ""`, because the `"only"` refusal reads
#' them.
#'
#' `tab` is `validate_prior()`'s table of the effective prior and `cn`
#' the posterior draws' column names, which give the simplex lengths and
#' the smooth and correlation names under brms's spelling.
#'
#' @noRd
prior_plan <- function(fit, tab, cn, mix_theta = FALSE) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  for (col in c("class", "coef", "group", "resp", "dpar", "nlpar",
                "prior", "source")) {
    tab[[col]] <- as.character(tab[[col]] %||% rep("", nrow(tab)))
    tab[[col]][is.na(tab[[col]])] <- ""
  }
  ct <- brms_coef_table(fit)
  multi <- length(fit$spec$responses) > 1L
  blocks <- fit$frame[["re_blocks"]]
  out <- list()
  add <- function(s) out[[length(out) + 1L]] <<- s
  intercepts <- character(0)
  # the smooth blocks, each with brms's sds_ base name: the draws carry
  # them in block order, one base per smooth term
  sm_blocks <- which(vapply(blocks, function(bk) {
    identical(bk[["covstruct"]], "smooth")
  }, NA))
  sds_bases <- unique(sub("_[0-9]+$", "", grep("^sds_", cn, value = TRUE)))
  sm_label <- vapply(blocks[sm_blocks], function(bk) {
    bk[["group_name"]] %||% bk[["term_label"]] %||% ""
  }, "")
  sm_lp <- lapply(blocks[sm_blocks], function(bk) {
    for (cp in bk[["components"]]) {
      lp <- fit$frame[["linpreds"]][[cp[["lp_key"]]]]
      if (!is.null(lp)) return(lp)
    }
    NULL
  })
  # one smooth term is one base, keyed by its predictor as well as its
  # label: s(x) in the location and s(x) in sigma are two terms, while
  # the penalties of one term are blocks under one label
  sm_key <- paste(sm_label, vapply(sm_lp, function(lp) {
    paste(lp[["resp"]] %||% "", lp[["dpar"]] %||% "")
  }, ""))
  sm_base <- sds_bases[match(sm_key, unique(sm_key))]
  for (r in names(fit$spec$responses)) {
    rspec <- fit$spec$responses[[r]]
    rr <- if (multi) r else ""
    fam <- rspec$family
    nls <- rspec$nlpars %||% character(0)
    pxs <- c(nls, setdiff(fam[["dpars"]] %||% character(0), nls))
    for (px in pxs) {
      nl <- px %in% nls
      rows <- which(ct$resp %in% r & ct$dpar %in% px)
      if (!length(rows)) next
      if (any(ct$natural[rows])) {
        nm <- ct$brms[rows[1L]]
        cls <- sub("[0-9]*(_.*)?$", "", nm)
        # a mixture's weights are one simplex, drawn after the dpars;
        # a GEV xi is brms's transformed tmp_xi, which brms never draws
        if (!is.null(fam[["mix"]]) && identical(cls, "theta")) next
        if (identical(fam[["family"]], "gen_extreme_value") &&
              identical(cls, "xi")) next
        row <- prior_tab_row(tab, tab$class == px & tab$resp == rr &
                               !nzchar(tab$coef) & !nzchar(tab$dpar) &
                               !nzchar(tab$nlpar) & !nzchar(tab$group))
        if (is.null(row)) next
        bd <- prior_dpar_bounds(cls, fit, rspec)
        add(prior_stmt(nm, nm, row$prior,
                       if (is.na(row$lb)) bd[1L] else row$lb,
                       if (is.na(row$ub)) bd[2L] else row$ub))
        next
      }
      own <- prior_px_rows(tab, rr, px, nl)
      has_icpt <- any(own & tab$class == "Intercept" & !nzchar(tab$coef))
      bcoef <- tab$coef[own & tab$class == "b" & nzchar(tab$coef)]
      # each coefficient's vector (b, bsp, bcs, bs with the predictor's
      # prefix) is its brms name less "_<coef>"
      vec <- character(0)
      coef <- character(0)
      stem <- ""
      for (i in rows) {
        bn <- ct$brms[i]
        if (has_icpt && grepl("_Intercept$", bn) &&
              grepl("(^|_)[(]Intercept[)]$", ct$internal[i])) {
          stem <- sub("_Intercept$", "", bn)
          next
        }
        hit <- bcoef[endsWith(bn, paste0("_", bcoef))]
        if (!length(hit)) next
        hit <- hit[which.max(nchar(hit))]
        vec <- c(vec, substr(bn, 1L, nchar(bn) - nchar(hit) - 1L))
        coef <- c(coef, hit)
      }
      p <- if (length(vec)) sub("^b(s|sp|cs)?", "", vec[1L]) else
        sub("^b", "", stem)
      vstem <- sub("_.*$", "", vec)
      cls_row <- prior_tab_row(tab, own & tab$class == "b" &
                                 !nzchar(tab$coef) & !nzchar(tab$group))
      b_stmts <- function(st) {
        for (v in unique(vec[vstem == st])) {
          cs <- coef[vec == v]
          crow <- lapply(cs, function(cc) {
            prior_tab_row(tab, own & tab$class == "b" & tab$coef == cc)
          })
          per <- any(vapply(crow, function(x) {
            !is.null(x) && identical(x$source, "user") &&
              !x$prior %in% c("", "(flat)")
          }, NA))
          lb <- if (!is.null(cls_row)) cls_row$lb else NA
          ub <- if (!is.null(cls_row)) cls_row$ub else NA
          if (per) {
            for (j in seq_along(cs)) {
              pr <- if (is.null(crow[[j]])) "" else crow[[j]]$prior
              add(prior_stmt(paste0(v, "_", cs[j]), v, pr, lb, ub))
            }
          } else {
            pr <- if (!is.null(cls_row)) cls_row$prior else ""
            # a vectorized coefficient row carries the class's density
            if (pr %in% c("", "(flat)") && length(crow) &&
                  !is.null(crow[[1L]]) &&
                  identical(crow[[1L]]$source, "(vectorized)")) {
              pr <- crow[[1L]]$prior
            }
            add(prior_stmt(v, v, pr, lb, ub))
          }
        }
      }
      b_stmts("b")
      if (has_icpt) {
        irow <- prior_tab_row(tab, own & tab$class == "Intercept" &
                                !nzchar(tab$coef))
        icoef <- own & tab$class == "Intercept" & nzchar(tab$coef) &
          tab$source == "user" & !tab$prior %in% c("", "(flat)")
        if (!any(icoef)) {
          add(prior_stmt(paste0("Intercept", p), paste0("Intercept", p),
                         irow$prior, irow$lb, irow$ub))
        } else {
          # an ordinal model with a threshold prior of its own: brms
          # writes one statement per threshold, its own row's density or
          # the class's, drawn as prior_Intercept__k (measured,
          # dev/priordraws-brms-me.R)
          K <- sum(startsWith(cn, paste0("b", p, "_Intercept[")))
          for (k in seq_len(K)) {
            krow <- prior_tab_row(tab, own & tab$class == "Intercept" &
                                    tab$coef == as.character(k))
            pr <- if (!is.null(krow) &&
                        !krow$prior %in% c("", "(flat)")) {
              krow$prior
            } else irow$prior
            add(prior_stmt(paste0("Intercept", p, "__", k),
                           paste0("Intercept", p), pr, irow$lb, irow$ub))
          }
        }
        intercepts <- c(intercepts, p)
      }
      for (i in which(own & tab$class == "simo")) {
        nm <- paste0("simo", p, "_", tab$coef[i])
        D <- sum(startsWith(cn, paste0(nm, "[")))
        # a simplex without a row of its own is sampled under brms's
        # default dirichlet(1) (mo_simplex_nlp())
        pr <- tab$prior[i]
        if (pr %in% c("", "(flat)")) pr <- "dirichlet(1)"
        if (D > 0L) {
          add(prior_stmt(nm, nm, pr, type = "simplex", dim = D))
        }
      }
      b_stmts("bsp")
      b_stmts("bcs")
      b_stmts("bs")
      k <- 0L
      for (j in seq_along(sm_blocks)) {
        lp <- sm_lp[[j]]
        if (is.null(lp) || !identical(lp[["resp"]] %||% r, r) ||
              !identical(lp[["dpar"]], px)) next
        if (j > 1L && sm_key[j] %in% sm_key[seq_len(j - 1L)]) next
        k <- k + 1L
        g <- sm_label[j]
        sg <- tab$class == "sd" & tab$group == g & !nzchar(tab$coef)
        srow <- prior_tab_row(tab, sg & own) %||% prior_tab_row(tab, sg) %||%
          prior_tab_row(tab, tab$class == "sd" & !nzchar(tab$group) &
                          !nzchar(tab$coef) & tab$resp %in% c(rr, ""))
        if (is.null(srow) || is.na(sm_base[j])) next
        add(prior_stmt(sm_base[j], paste0("sds", p, "_", k), srow$prior,
                       if (is.na(srow$lb)) 0 else srow$lb, srow$ub))
      }
      # brms's residual correlation terms are part of the location's
      # predictor
      if (px %in% c("mu", (rspec$primary_dpars %||% "mu")[1L]) && !nl) {
        ac_cls <- c("ar", "ma", "sderr", "cosy", "cortime", "lagsar",
                    "errorsar", "car", "sdcar", "rhocar")
        ac_rows <- which(tab$class %in% ac_cls & tab$resp == rr &
                           !nzchar(tab$coef) & !nzchar(tab$group))
        for (i in ac_rows[!duplicated(tab$class[ac_rows])]) {
          nm <- paste0(tab$class[i], if (multi) paste0("_", r))
          add(prior_stmt(nm, nm, tab$prior[i], tab$lb[i], tab$ub[i]))
        }
      }
    }
    if (!is.null(fam[["mix"]]) && isTRUE(mix_theta)) {
      dps <- fam[["dpars"]] %||% character(0)
      K <- sum(grepl("^theta[0-9]+$", dps)) + 1L
      if (K > 1L) {
        add(prior_stmt(paste0("theta", if (multi) paste0("_", r)),
                       "theta", "dirichlet(1)", type = "simplex", dim = K))
      }
    }
  }
  if (isTRUE(fit$spec$rescor)) {
    row <- prior_tab_row(tab, tab$class == "rescor")
    if (!is.null(row)) {
      add(prior_stmt("rescor", "Lrescor", row$prior, type = "lkj",
                     dim = length(fit$spec$responses)))
    }
  }
  # group-level terms: one standard-deviation statement per grouping
  # factor, per coefficient where a coefficient has a prior of its own,
  # and one correlation statement
  seen <- character(0)
  id <- 0L
  for (bi in seq_along(blocks)) {
    bk <- blocks[[bi]]
    if (bi %in% sm_blocks || !length(bk[["theta_idx"]])) next
    if (isTRUE(frmtmb_covstruct_info(bk[["covstruct"]])[["field"]])) next
    g <- brms_group_name(bk)
    id <- id + 1L
    if (g %in% seen) next
    seen <- c(seen, g)
    rn <- brms_re_rnames(fit, bk)
    grow <- prior_tab_row(tab, tab$class == "sd" & tab$group == g &
                            !nzchar(tab$coef))
    if (is.null(grow)) next
    crow <- lapply(rn, function(cc) {
      prior_tab_row(tab, tab$class == "sd" & tab$group == g &
                      tab$coef == cc)
    })
    per <- any(vapply(crow, function(x) {
      !is.null(x) && identical(x$source, "user") &&
        !x$prior %in% c("", "(flat)")
    }, NA))
    lb <- if (is.na(grow$lb)) 0 else grow$lb
    if (per) {
      for (j in seq_along(rn)) {
        pr <- if (is.null(crow[[j]]) ||
                    crow[[j]]$prior %in% c("", "(flat)")) {
          grow$prior
        } else crow[[j]]$prior
        add(prior_stmt(paste0("sd_", g, "__", rn[j]), paste0("sd_", id),
                       pr, lb, grow$ub))
      }
    } else {
      add(prior_stmt(paste0("sd_", g), paste0("sd_", id), grow$prior, lb,
                     grow$ub))
    }
    if (length(rn) > 1L && any(startsWith(cn, paste0("cor_", g, "__")))) {
      crw <- prior_tab_row(tab, tab$class == "cor" & tab$group == g) %||%
        prior_tab_row(tab, tab$class == "cor" & !nzchar(tab$group))
      if (!is.null(crw)) {
        add(prior_stmt(paste0("cor_", g), paste0("L_", id), crw$prior,
                       type = "lkj", dim = length(rn)))
      }
    }
  }
  for (s in me_prior_stmts(fit, tab)) add(s)
  attr(out, "intercepts") <- unique(intercepts)
  out
}

#' The prior statements of a model's `me()` terms, which brms's Stan
#' program writes after every other one (`stan_Xme()`,
#' dev/meme-out/brms-order.txt). Per group of terms that share a `gr =`,
#' numbered `i` in order of first appearance: `meanme_<i>` and
#' `sdme_<i>`, each a vector over the group's terms, then `Lme_<i>` where
#' the group is correlated. A vector with a coefficient's own prior is
#' one statement per term, which brms draws as `prior_meanme_<i>__<k>`
#' and renames to `prior_meanme_<i>_meanme_<coef>`: its `rename_prior()`
#' replaces the trailing `_<k>` with the term's whole name
#' (dev/meme-out/brms-names.txt). A class-wide one keeps `meanme_<i>`,
#' and the correlation keeps brms's Stan name `Lme_<i>`, which brms's
#' renaming does not reach either.
#'
#' @noRd
me_prior_stmts <- function(fit, tab) {
  me <- fit$frame[["me"]]
  if (is.null(me)) return(list())
  out <- list()
  coefs <- vapply(me$terms, `[[`, "", "coef")
  flat <- function(x) is.null(x) || x$prior %in% c("", "(flat)")
  for (i in seq_along(me$groups)) {
    g <- me$groups[[i]]
    cs <- coefs[g$K]
    for (cls in c("meanme", "sdme")) {
      stan <- paste0(cls, "_", i)
      crow <- prior_tab_row(tab, tab$class == cls & !nzchar(tab$coef))
      krow <- lapply(cs, function(cc) {
        prior_tab_row(tab, tab$class == cls & tab$coef == cc)
      })
      per <- any(vapply(krow, function(x) {
        !flat(x) && identical(x$source, "user")
      }, NA))
      lb <- if (!is.null(crow)) crow$lb else NA
      ub <- if (!is.null(crow)) crow$ub else NA
      # sdme is declared vector<lower=0>, whatever its row says
      if (identical(cls, "sdme") && is.na(lb)) lb <- 0
      cpr <- if (is.null(crow)) "" else crow$prior
      if (per) {
        for (j in seq_along(cs)) {
          pr <- if (flat(krow[[j]])) cpr else krow[[j]]$prior
          out[[length(out) + 1L]] <- prior_stmt(
            paste0(stan, "_", cls, "_", cs[j]), stan, pr, lb, ub)
        }
      } else {
        out[[length(out) + 1L]] <- prior_stmt(stan, stan, cpr, lb, ub)
      }
    }
    if (isTRUE(g$cor)) {
      crw <- prior_tab_row(tab, tab$class == "corme" & tab$group == g$gr) %||%
        prior_tab_row(tab, tab$class == "corme" & !nzchar(tab$group))
      if (!is.null(crw)) {
        out[[length(out) + 1L]] <- prior_stmt(
          paste0("Lme_", i), paste0("Lme_", i), crw$prior, type = "lkj",
          dim = length(g$K))
      }
    }
  }
  out
}

#' The density object of a prior string, read from the effective prior
#' it came from, or parsed by `set_prior()` when it is not there.
#'
#' @noRd
prior_dist_of <- function(s, eff) {
  for (e in unclass(eff %||% list())) {
    if (identical(e$prior, s)) return(e$dist)
  }
  unclass(set_prior(s))[[1L]]$dist
}

#' `n` draws from one density, truncated to `[lb, ub]` where either is
#' finite. Truncation is by the inverse CDF, which gives brms's
#' rejection sampler's distribution without its loop.
#'
#' @noRd
rprior_real <- function(dist, n, lb = NA, ub = NA, what = "") {
  k <- dist$kind
  f <- switch(k,
    normal = list(
      r = function(n) stats::rnorm(n, dist$location, dist$scale),
      p = function(x) stats::pnorm(x, dist$location, dist$scale),
      q = function(u) stats::qnorm(u, dist$location, dist$scale)),
    t = list(
      r = function(n) dist$location + dist$scale * stats::rt(n, dist$df),
      p = function(x) stats::pt((x - dist$location) / dist$scale, dist$df),
      q = function(u) dist$location + dist$scale * stats::qt(u, dist$df)),
    exponential = list(
      r = function(n) stats::rexp(n, dist$rate),
      p = function(x) stats::pexp(x, dist$rate),
      q = function(u) stats::qexp(u, dist$rate)),
    logistic = list(
      r = function(n) stats::rlogis(n, dist$location, dist$scale),
      p = function(x) stats::plogis(x, dist$location, dist$scale),
      q = function(u) stats::qlogis(u, dist$location, dist$scale)),
    gamma = list(
      r = function(n) stats::rgamma(n, dist$shape, rate = dist$rate),
      p = function(x) stats::pgamma(x, dist$shape, rate = dist$rate),
      q = function(u) stats::qgamma(u, dist$shape, rate = dist$rate)),
    # Stan's inv_gamma(alpha, beta) is 1 / gamma(alpha, rate = beta)
    inv_gamma = list(
      r = function(n) 1 / stats::rgamma(n, dist$shape, rate = dist$scale),
      p = function(x) {
        ifelse(x <= 0, 0, stats::pgamma(1 / pmax(x, 0), dist$shape,
                                        rate = dist$scale,
                                        lower.tail = FALSE))
      },
      q = function(u) 1 / stats::qgamma(u, dist$shape, rate = dist$scale,
                                        lower.tail = FALSE)),
    beta = list(
      r = function(n) stats::rbeta(n, dist$shape1, dist$shape2),
      p = function(x) stats::pbeta(x, dist$shape1, dist$shape2),
      q = function(u) stats::qbeta(u, dist$shape1, dist$shape2)),
    uniform = list(
      r = function(n) stats::runif(n, dist$lower, dist$upper),
      p = function(x) stats::punif(x, dist$lower, dist$upper),
      q = function(u) stats::qunif(u, dist$lower, dist$upper)),
    lognormal = list(
      r = function(n) stats::rlnorm(n, dist$location, dist$scale),
      p = function(x) stats::plnorm(x, dist$location, dist$scale),
      q = function(u) stats::qlnorm(u, dist$location, dist$scale)),
    NULL)
  if (is.null(f)) {
    frm_stop("sample_prior = \"yes\": no draws from the ", k, " prior of ",
             what, call. = FALSE)
  }
  lo <- is.finite(lb)
  hi <- is.finite(ub)
  if (!lo && !hi) return(f$r(n))
  pl <- if (lo) f$p(lb) else 0
  pu <- if (hi) f$p(ub) else 1
  if (!(pu > pl)) {
    frm_stop("sample_prior = \"yes\": the prior of ", what, " has no mass ",
             "between its bounds [", lb, ", ", ub, "], so there is nothing ",
             "to draw", call. = FALSE)
  }
  x <- f$q(stats::runif(n, pl, pu))
  # the inverse CDF can round onto or past a bound
  if (lo) x <- pmax(x, lb)
  if (hi) x <- pmin(x, ub)
  x
}

#' The prior draws of a plan as a matrix with brms's column names
#' (`prior_<name>`, `prior_<name>[k]` for a simplex), in brms's column
#' order: `reorder_pars()` puts a parameter whose STAN name ends in
#' `_Intercept` first and keeps the rest in the order of the Stan
#' program. It reads the names before brms renames them, so only the
#' location's own `prior_Intercept` moves: `prior_sd_g__Intercept` was
#' `prior_sd_1__1` there (dev/priordraws-findings.md, section 1).
#'
#' @noRd
prior_draws_make <- function(plan, n, eff) {
  plan <- Filter(function(s) nzchar(s$prior), plan)
  if (!length(plan)) return(NULL)
  # a statement per coefficient is prior_<par>__<k> in brms's generated
  # quantities, whatever the parameter, so it never moves
  sn <- vapply(plan, function(s) {
    if (identical(s$name, s$stan)) paste0("prior_", s$stan) else ""
  }, "")
  plan <- plan[order(!grepl("_Intercept(_[[:digit:]]+)?$", sn))]
  cols <- lapply(plan, function(s) {
    d <- prior_dist_of(s$prior, eff)
    nm <- paste0("prior_", s$name)
    if (identical(s$type, "simplex")) {
      # a simplex of one weight (a two-value mo()) is 1 in every draw,
      # as brms's dirichlet_rng() of length one is
      if (s$dim == 1L) {
        return(matrix(1, n, 1L, dimnames = list(NULL, paste0(nm, "[1]"))))
      }
      a <- rep_len(d$alpha %||% 1, s$dim)
      g <- vapply(a, function(ak) stats::rgamma(n, ak), numeric(n))
      g <- matrix(g, n, s$dim)
      w <- g / rowSums(g)
      colnames(w) <- paste0(nm, "[", seq_len(s$dim), "]")
      return(w)
    }
    if (identical(s$type, "lkj")) {
      # the marginal of one off-diagonal element of an LKJ(eta)
      # correlation matrix of dimension M is 2 Beta(a, a) - 1,
      # a = eta - 1 + M / 2; brms keeps element [1, 2] of a draw
      a <- d$eta - 1 + s$dim / 2
      v <- 2 * stats::rbeta(n, a, a) - 1
      return(matrix(v, n, 1L, dimnames = list(NULL, nm)))
    }
    v <- rprior_real(d, n, s$lb, s$ub, what = s$name)
    matrix(v, n, 1L, dimnames = list(NULL, nm))
  })
  do.call(cbind, cols)
}

#' brms's refusal of `sample_prior = "only"` on a model with a flat
#' prior, naming the parameter of the first such statement as brms's
#' Stan program names it. brms refuses only a parameter with no prior
#' statement at all; a vector of coefficients some of which have a prior
#' and some not is refused here too, because sampling the flat ones from
#' nothing has no proper result.
#'
#' @noRd
prior_only_check <- function(plan) {
  flat <- Filter(function(s) !nzchar(s$prior), plan)
  if (!length(flat)) return(invisible(NULL))
  prior_only_stop(flat[[1L]]$stan)
}

#' brms's words for that refusal, naming parameter `par`, with `why`
#' after them where frmtmb's reason is not brms's.
#'
#' @noRd
prior_only_stop <- function(par, why = NULL) {
  frm_stop("Sampling from priors is not possible as some parameters have ",
           "no proper priors. Error occurred for parameter '", par, "'.",
           if (!is.null(why)) paste0(" ", why), call. = FALSE)
}

#' The outer parameters that carry no density under the effective prior
#' and the model's own extra densities, by name, for the `"only"`
#' refusal: a backstop for parameters brms has no class for (a
#' structured covariance's own parameters, say), which the plan cannot
#' see.
#'
#' @noRd
prior_only_uncovered <- function(fit, entries, mix_theta = FALSE) {
  tpl <- fit$frame[["par_template"]]
  map <- fit$frame[["map"]] %||% list()
  # the outer parameters' labels, in the layout draws_par_index() gives
  # without the random effects
  nms <- brms_par_labels(fit, include_random = FALSE)
  pidx <- draws_par_index(fit, laplace = TRUE)
  covered <- character(0)
  for (e in entries) {
    # brms's (flat) held as an entry (the map's Jacobian alone, or an
    # LKJ made flat on the Cholesky factor) is no proper density
    if (identical(e$dist$kind, "flat") || isTRUE(e$dist$flat)) next
    j <- pidx[[e$comp]]
    if (is.null(j)) next
    idx <- e$idx
    if (identical(e$comp, "betad")) {
      fx <- fit$frame[["betad_fixed_idx"]] %||% integer(0)
      keep <- setdiff(seq_along(tpl[["betad"]]), fx)
      idx <- match(idx, keep)
    }
    covered <- c(covered, nms[j[idx[!is.na(idx)]]])
  }
  # a mo() simplex and a mixture's weights carry brms's dirichlet()
  # outside the resolved entries (sample_extra_nlp())
  for (tm in simplex_frame_terms(fit)) {
    covered <- c(covered, nms[pidx[[tm$zeta]] %||% integer(0)])
  }
  if (isTRUE(mix_theta)) {
    ct <- brms_coef_table(fit)
    th <- which(ct$natural & grepl("^theta[0-9]+", ct$brms))
    covered <- c(covered, nms[th])
  }
  fixed <- unlist(lapply(names(map), function(cp) {
    j <- pidx[[cp]]
    m <- map[[cp]]
    if (is.null(j) || length(m) != length(j)) return(character(0))
    nms[j[is.na(m)]]
  }))
  out <- setdiff(nms, c(covered, fixed))
  out[!is.na(out)]
}

#' Prior draws
#'
#' brms's `prior_draws()`: the draws [frm_sample()] made from the
#' priors with `sample_prior = "yes"`, under brms's names without the
#' `prior_` prefix.
#'
#' Each column is one prior statement of brms's Stan program for the
#' model. A prior on a whole class is one column (`b`, `sd_g`, `sigma`);
#' where a coefficient has a prior of its own, each coefficient of that
#' vector is a column (`b_x`, `sd_g__Intercept`). A flat prior has no
#' column, so with brms's defaults there is no `b`. A bounded parameter
#' is drawn from its prior truncated to the bounds (`sigma` and `sd` at
#' 0, `nu` at 1). A correlation is one off-diagonal element of an
#' LKJ draw (`cor_g`, `rescor`), a `mo()` simplex or a mixture's weights
#' a Dirichlet draw (`simo_mox1[1]`, ...). `Intercept` is drawn from the
#' prior of the intercept at the centered predictors, which is where
#' brms puts that prior.
#'
#' With `variable`, each name gets the column whose name is the longest
#' match at its start, and the draws in random order, as brms does, so
#' that two names sharing a column are not paired draw by draw: `b_x`
#' and `b_z` both read `b`. A name nothing matches is left out. The
#' population-level intercept of a model with a centered intercept
#' (`b_Intercept`) has no prior draws: brms warns and leaves it out, and
#' so does this method.
#'
#' Draws made with `sample_prior = "no"` (the default) or `"only"` have
#' no prior draws, and the result is an empty data frame. Under
#' `"only"` the draws themselves are draws from the prior.
#'
#' The draws are made directly from each prior with R's random number
#' generators, as brms makes them in Stan's generated quantities, so they
#' are independent and do not depend on the likelihood. With `seed`
#' given to [frm_sample()] they are reproducible.
#'
#' @param x A `frmtmb_draws` from [frm_sample()].
#' @param variable Names of the variables whose prior draws to return,
#'   or `NULL` for every prior draw.
#' @param pars brms's deprecated alias of `variable`; it warns.
#' @param ... Unused.
#' @return A data frame with one column per prior (or per name in
#'   `variable`) and one row per post-warmup draw.
#' @seealso [hypothesis()][frmtmb::hypothesis], which divides the
#'   posterior density of a point hypothesis by its prior density read
#'   from these draws.
#' @examples
#' \donttest{
#' if (requireNamespace("tmbstan", quietly = TRUE) &&
#'     requireNamespace("rstan", quietly = TRUE) &&
#'     !frmtmb.sample:::tmbstan_build_broken()) {
#'   set.seed(9)
#'   dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#'   dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#'   ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
#'                    prior = set_prior("normal(0, 1)", class = "b"),
#'                    sample_prior = "yes", chains = 1, iter = 500,
#'                    seed = 1, refresh = 0)
#'   head(prior_draws(ds))
#'   # the prior draws of one coefficient are the class's
#'   summary(prior_draws(ds, "b_x")$b_x)
#'   # the Savage-Dickey evidence ratio of a point hypothesis reads them
#'   hypothesis(ds, "x = 0")
#' }
#' }
#' @export
prior_draws <- function(x, ...) UseMethod("prior_draws")

#' @rdname prior_draws
#' @exportS3Method brms::prior_draws
#' @export
prior_draws.frmtmb_draws <- function(x, variable = NULL, pars = NULL,
                                     ...) {
  frm_check_dots(...)
  if (!is.null(pars)) {
    frm_warning("Argument 'pars' is deprecated. Please use argument ",
                "'variable' instead.", call. = FALSE)
    variable <- pars
  }
  if (!is.null(variable)) variable <- as.character(variable)
  cn <- colnames(x$draws)
  pn <- unique(cn[grepl("^prior_", cn)])
  if (!length(pn)) return(data.frame(NULL))
  draws <- as.data.frame(x$draws[, pn, drop = FALSE])
  rownames(draws) <- NULL
  names(draws) <- sub("^prior_", "", pn)
  if (is.null(variable)) return(draws)
  icpt <- x$prior_info$intercepts %||% character(0)
  one <- function(v) {
    ml <- vapply(names(draws), function(nm) {
      attr(regexpr(paste0("^", prior_escape(nm)), v), "match.length")
    }, 1L)
    if (max(ml) == -1L || prior_ignore(v, icpt)) return(NULL)
    take <- match(max(ml), ml)
    # in random order, so that two names reading one column are not
    # paired draw by draw
    structure(list(draws[sample(nrow(draws)), take]), names = v)
  }
  out <- Filter(Negate(is.null), lapply(variable, one))
  data.frame(out, check.names = FALSE)
}

#' brms's `escape_all()`: every regular-expression metacharacter of a
#' name escaped.
#'
#' @noRd
prior_escape <- function(x) {
  for (s in c(".", "*", "+", "?", "^", "$", "(", ")", "[", "]", "|")) {
    x <- gsub(s, paste0("\\", s), x, fixed = TRUE)
  }
  x
}

#' brms's `ignore_prior()`: a population-level intercept of a predictor
#' with a centered intercept has no prior draws, and brms warns and
#' leaves it out. `icpt` holds those predictors' prefixes (`""`,
#' `"_sigma"`). The test is brms's string arithmetic, which leaves an
#' ordinal threshold `b_Intercept[1]` to the prefix match.
#'
#' @noRd
prior_ignore <- function(v, icpt) {
  if (!grepl("^b_.*Intercept($|[[])", v) || !length(icpt)) return(FALSE)
  p <- sub("_Intercept($|[[])", "", sub("^b", "", v))
  out <- p %in% icpt
  if (out) {
    frm_warning("Sampling from the prior of an overall intercept is not ",
                "possible by default. See the documentation of the ",
                "'sample_prior' argument in help('frm_sample').",
                call. = FALSE)
  }
  out
}
