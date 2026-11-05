# brms's latent residuals: `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()`
# on a family with no residual `y - mu`.
#
# WHAT brms FITS. brms 2.23.0 marks gaussian and student alone with the
# "residuals" special. On every other family stan_ac() adds a vector of
# latent residuals to the linear predictor of mu:
#
#   eta_i = ... + err_i,   zerr ~ std_normal(),   sderr > 0
#
# under cov = TRUE, cosy() and unstr(), with one correlated draw per group,
#
#   err_g = sderr * L_g * zerr_g,   L_g L_g' = C,
#
# C at unit innovation variance for ar(), ma(), arma() (so sderr is the
# innovation sd, as sigma is for the natural residuals), and a unit-diagonal
# correlation for cosy() and unstr(). The groups read leading blocks of one
# Cholesky factor, so the lag is counted in ROWS of a group, as brms counts
# it; only unstr() reads the time levels (Jtime_tg). Under cov = FALSE only
# the AR part exists, and the program is different:
#
#   err = sderr * zerr                      (iid)
#   mu_i += sum_j ar_j err_{i-j}            (a moving average of err)
#
# with `ar` boxed on (-1, 1) entrywise. That sum enters on the link scale
# where the family has a built-in density on that scale, and after the
# inverse link everywhere else (stan_eta_transform(), `latres_on_response()`).
#
# WHY NOT brms's PARAMETERIZATION. One latent per row is a random effect the
# inner problem has to solve at every outer step, and brms's err = L z makes
# its Hessian dense within a group: I + L' W L with L dense. That is k^3 per
# group per Newton step, and one series of a thousand rows would be a million
# stored entries. The chain structure is sparse in other coordinates:
#
#   ar1    s_t | s_{t-1} ~ N(phi s_{t-1}, 1), s_1 ~ N(0, 1 / (1 - phi^2));
#          err_t = sderr s_t. Hessian: tridiagonal.
#   ma1    s_0, ..., s_k iid; err_t = sderr (s_t + theta s_{t-1}). One extra
#          latent per group. Hessian: tridiagonal.
#   arma1  s_0, ..., s_k the ar1 chain; err_t = sderr (s_t + theta s_{t-1}),
#          because (1 - phi B)(1 + theta B) s = (1 + theta B) e. Tridiagonal.
#   condar z iid; err_t = sderr (z_t + sum_j ar_j z_{t-j}). Band of 2p + 1.
#   cosy   z_t, u_g iid; err_t = sderr (sqrt(1 - rho) z_t + sqrt(rho) u_g).
#          One extra latent per group. Hessian: an arrow per group.
#   unstr  u_g of length d iid; err_t = sderr L[J_t, ] u_g. A dense d x d
#          block per group, brms's cost, and no Cholesky factor of a pattern's
#          submatrix on the tape.
#
# THE LAPLACE VALUE DOES NOT MOVE. The likelihood reads the latents only
# through err = A s, and s is gaussian. The inner Hessian is Q_s + A'WA, and
# log|Q_s + A'WA| - log|Q_s| = log|I + A Q_s^-1 A' W| by Sylvester's identity,
# which reads A only through Cov(err). So every coordinate system with the
# same Cov(err) gives one Laplace approximation, whatever its dimension.
# Check C of the identity with brms's log_prob measures the joint density;
# dev/latres-findings.md measures this value too.
#
# THE ROUTE. The block is a field of the covariance-structure registry, so
# everything that reads a field (se.fit, the outer covariance, frm_sample())
# gets the latents with no second path. All its parameters live in its theta
# (log sderr first), because expand_b() and the delta method read theta.
# frame$autocor stays what it is, the terms that replace the response density:
# a latent term does not replace it. `frame_ac_terms()` is the one reader that
# sees both, for the prior classes.

#' Whether a parsed autocorrelation term on family `fam` is brms's latent
#' residual: anything but `fcor()` and a registered term, on anything but the
#' two families brms gives a residual. Exported nowhere: the frame asks once.
#'
#' @noRd
autocor_is_latent_on <- function(ac, fam) {
  !is.null(ac) && !autocor_is_seam(ac) &&
    !fam[["family"]] %in% c("gaussian", "student")
}

#' The refusals of a latent residual, in brms's words.
#'
#' `stan_ac()` and `.arma()` say these; each was read from a run of brms
#' 2.23.0 (dev/latres-grid.R). What the mixture and multi-category families
#' say is in `pull_autocor()`, which refuses before the frame sees a family.
#'
#' @noRd
check_autocor_latent <- function(resp, ac) {
  fn <- ac[["fn"]]
  fam <- resp$family
  if (!is.null(fam[["mix"]])) {
    frm_stop(fn, "(): Autocorrelation terms cannot be applied in mixture ",
             "models. brms refuses the same call", call. = FALSE)
  }
  # brms's guard is is.btnl() on mu: a nonlinear shape or zi beside a linear
  # mu has a linear predictor to add to, and brms builds it
  if ("mu" %in% resp$nl_dpars) {
    frm_stop(fn, "(): Latent residuals are not implemented for non-linear ",
             "models. brms refuses the same call: '", fam[["family"]],
             "' has no residual y - mu, so the term would be a latent ",
             "vector added to the linear predictor of mu, and a nonlinear ",
             "body has none", call. = FALSE)
  }
  # a response matrix with one mu by the parser's reading, which brms's
  # program for the same call does not compile
  if (fam[["family"]] %in% c("categorical", multinomial_families) ||
        fam_is_simplex(fam)) {
    frm_stop(fn, "(): Explicit covariance terms can only be specified on ",
             "'mu', and '", fam[["family"]], "' has one predictor per ",
             "category. brms refuses the same model", call. = FALSE)
  }
  if (ac[["struct"]] %in% c("ar", "ma", "arma")) {
    if (isTRUE(ac[["cov"]]) && (ac[["p"]] > 1L || ac[["q"]] > 1L)) {
      frm_stop(fn, "(): Covariance formulation of ARMA structures is only ",
               "possible for effects of maximal order one. brms refuses the ",
               "same call", call. = FALSE)
    }
    if (isFALSE(ac[["cov"]]) && ac[["q"]] > 0L) {
      frm_stop(fn, "(): Please set cov = TRUE when modeling MA structures ",
               "for this family. brms refuses the same call: AR and MA ",
               "structures cannot be distinguished when using a single ",
               "vector of latent residuals", call. = FALSE)
    }
  }
  if (!is.null(fam_structure(fam)[["loglik"]])) {
    frm_stop(fn, "(): the '", fam[["family"]], "' family has a likelihood of ",
             "its own that does not read a linear predictor row by row, so ",
             "no latent residual can be added to it", call. = FALSE)
  }
  ac[["latent"]] <- TRUE
  ac
}

# brms's `sbi_<link>` specials on the families frmtmb has, in the shape
# stan_eta_transform() reads them: the family's density is written on the link
# scale, so no inverse link stands between the predictor and the ARMA loop.
# A censored or truncated row has no such density on any of these, so
# `latres_on_response()` reads this table only without them. test-latres-brms
# compares it with brms's own stan_eta_transform().
latres_builtin_links <- list(
  log = c("poisson", "negbinomial", "negbinomial2", "geometric", "compois",
          "cox", "zero_inflated_poisson", "zero_inflated_negbinomial",
          "hurdle_poisson", "hurdle_negbinomial"),
  logit = c("bernoulli", "binomial", "zero_inflated_binomial")
)

#' Whether the AR part of a `cov = FALSE` latent residual adds to `mu` after
#' the inverse link, brms's `stan_eta_transform()`: the program writes
#' `mu = inv_link(mu)` before the ARMA loop on a non-identity link, unless the
#' family's density reads the link scale or the family applies its link itself
#' (the ordinal and multi-category families' "joint_link").
#'
#' @noRd
latres_on_response <- function(fam, av) {
  link <- fam[["links"]][["mu"]][["name"]] %||% "identity"
  if (identical(link, "identity")) return(FALSE)
  if (identical(fam[["type"]], "ordinal")) return(FALSE)
  cdf <- !is.null(av[["cens"]]) || !is.null(av[["trunc_lb"]]) ||
    !is.null(av[["trunc_ub"]])
  # cox() alone has its density on the log scale under a cdf too
  builtin <- fam[["family"]] %in% latres_builtin_links[[link]] &&
    (!cdf || identical(fam[["family"]], "cox"))
  !builtin
}

#' The range of `mu` where a latent residual can carry it out of the family's
#' support: `"pos"`, `"unit"`, or `NULL` where it cannot. It can on a link that
#' does not map onto the support (`identity`, `inverse`, `1/mu^2`, which a
#' latent per row can push across zero) and under a post-link AR part
#' (`on_response`), whatever the link. Stan rejects such a draw; see
#' `latres_guard()`.
#'
#' @noRd
latres_mu_range <- function(fam, on_response) {
  link <- fam[["links"]][["mu"]][["name"]] %||% "identity"
  if (!on_response && !link %in% c("identity", "inverse", "1/mu^2")) {
    return(NULL)
  }
  nm <- fam[["family"]]
  if (nm %in% c("poisson", "negbinomial", "negbinomial2", "geometric",
                "compois", "Gamma", "exponential", "weibull", "frechet",
                "inverse.gaussian", "zero_inflated_poisson",
                "zero_inflated_negbinomial", "hurdle_poisson",
                "hurdle_negbinomial", "hurdle_gamma")) {
    return("pos")
  }
  if (nm %in% c("bernoulli", "binomial", "beta", "beta_binomial",
                "zero_inflated_binomial", "zero_inflated_beta",
                "zero_inflated_beta_binomial", "zero_one_inflated_beta")) {
    return("unit")
  }
  NULL
}

#' `mu` unchanged where it is in the range of its link, and NaN elsewhere,
#' on the tape too. Stan's density rejects a mean outside its range, and brms
#' can put one there, since the AR part is added after the inverse link. A
#' density on doubles says NaN as well, but on the tape RTMB folds `0 * log(mu)`
#' at a zero count (dev/rtmb-pitfalls.md, item 22), so a row with `y = 0` would
#' keep a finite density at a negative rate.
#'
#' @noRd
latres_guard <- function(mu, a) {
  rg <- a[["mu_range"]]
  if (is.null(rg)) return(mu)
  out <- mu + nan_unless_pos(mu)
  if (identical(rg, "unit")) out <- out + nan_unless_pos(1 - mu)
  out
}

#' The fixed data of a latent residual block: which kind of chain, and the
#' gathers its tape reads. Everything here is index data, built once.
#'
#' `tidx` is the row's time level, `gidx` its group. `lag_idx[[i]]` is the row
#' `i` steps earlier in its group, by position, or `n + 1` where there is none.
#'
#' @noRd
latres_aux <- function(ac, gidx, tidx, n, on_response) {
  st <- ac[["struct"]]
  ps <- autocor_cond_positions(gidx, tidx, n, max(ac[["p"]], 1L))
  G <- length(ps$len)
  kind <- if (st == "cosy") {
    "cosy"
  } else if (st == "unstr") {
    "unstr"
  } else if (isFALSE(ac[["cov"]])) {
    "condar"
  } else if (ac[["q"]] == 0L) {
    "ar1"
  } else if (ac[["p"]] == 0L) {
    "ma1"
  } else {
    "arma1"
  }
  # `first` holds every group's first row, in the order pos_rows lists the
  # groups, so position t in that order is group t of the `n + g` latents
  first <- ps$pos_rows[[1L]]
  prev <- ps$lag_idx[[1L]]
  aux <- list(kind = kind, ac = ac, n = n, G = G, m = n, first = first,
              group = gidx, tidx = tidx, mu_range = ac[["mu_range"]],
              pos_rows = ps$pos_rows)
  if (kind %in% c("ar1", "arma1", "ma1")) {
    if (kind != "ar1") {
      prev[first] <- n + seq_len(G)
      aux$m <- n + G
    }
    aux$prev <- prev
    aux$start <- if (kind == "ar1") first else n + seq_len(G)
    aux$stepped <- setdiff(seq_len(n), if (kind == "ar1") first)
  } else if (kind == "condar") {
    aux$lag_idx <- ps$lag_idx
    aux$on_response <- isTRUE(on_response)
  } else if (kind == "cosy") {
    g_of <- integer(n)
    for (t in seq_along(ps$pos_rows)) {
      g_of[ps$pos_rows[[t]]] <- seq_along(ps$pos_rows[[t]])
    }
    aux$u <- n + g_of
    aux$m <- n + G
  } else {
    d <- ac[["d"]]
    g_of <- integer(n)
    for (t in seq_along(ps$pos_rows)) {
      g_of[ps$pos_rows[[t]]] <- seq_along(ps$pos_rows[[t]])
    }
    aux$J <- tidx
    aux$g_of <- g_of
    aux$d <- d
    aux$m <- G * d
  }
  aux
}

#' Build a latent residual term once its group and time index are known:
#' `autocor_block()`'s bookkeeping, with `autocor_cond_positions()` in place of
#' the patterns, since nothing here factorizes a covariance per pattern.
#'
#' @noRd
latres_block <- function(ac, resp, gidx, ti, n, on_response) {
  if (isFALSE(ac[["cov"]])) {
    autocor_check_lag(ac[["fn"]], max(tabulate(gidx)), ac[["p"]])
  } else if (max(tabulate(gidx)) < 2L) {
    # the term has its parameter, but one row per group shows only its product
    # with sderr: brms builds this and its chains walk the ridge
    frm_stop(ac[["fn"]], "(): every group has one row, so the correlation of ",
             "two rows of a group is never seen and cannot be told from ",
             "sderr. Check the grouping variable; (1 | row) is the model of ",
             "one latent per row with no correlation", call. = FALSE)
  }
  ac[["npar"]] <- autocor_npar(ac)
  ac[["resp"]] <- resp
  ac[["n_groups"]] <- max(gidx)
  ac[["aux"]] <- latres_aux(ac, gidx, ti$idx, n, on_response)
  ac
}

#' The coefficients of a latent term from the part of its block's theta after
#' log sderr. Each chart is the bound brms declares: (-1, 1) on an AR or MA
#' coefficient, entrywise under cov = FALSE too, and (0, 1) on cosy, whose
#' negative range a shared latent cannot reach. Tape-safe.
#'
#' @noRd
latres_coefs <- function(th, a) {
  switch(a[["kind"]],
    ar1 = list(ar = autocor_pacf(th[1L])),
    ma1 = list(ma = autocor_pacf(th[1L])),
    arma1 = list(ar = autocor_pacf(th[1L]), ma = autocor_pacf(th[2L])),
    condar = list(ar = 2 / (1 + exp(-th)) - 1),
    cosy = list(rho = 1 / (1 + exp(-th[1L]))),
    list())
}

#' Log density of a latent block's `b`, with the normalizing constant, so the
#' joint at the modes differs from Stan's only by the map (see the top of this
#' file). The AR chain's stationary start reads `1 - phi^2 = 1 / (1 + t^2)`
#' exactly, where `1 - phi^2` loses its digits as `phi` nears 1.
#'
#' @noRd
latres_lpdf <- function(b, theta, blk) {
  a <- blk[["aux"]]
  kind <- a[["kind"]]
  ll <- -0.5 * a[["m"]] * log(2 * pi)
  if (!kind %in% c("ar1", "arma1")) return(ll - 0.5 * sum(b * b))
  t1 <- theta[2L]
  phi <- autocor_pacf(t1)
  st <- a[["stepped"]]
  e <- b[st] - phi * b[a[["prev"]][st]]
  s0 <- b[a[["start"]]]
  ll - 0.5 * sum(e * e) - 0.5 * sum(s0 * s0) / (1 + t1 * t1) -
    0.5 * length(s0) * log(1 + t1 * t1)
}

#' The latent residual of every row on the link scale, in row order, from the
#' block's `b` and its theta: `expand_b()`'s map. Tape-safe and numeric-safe.
#' An `on_response` cov = FALSE term keeps its AR part out of this, in
#' `latres_shift()`.
#'
#' @noRd
latres_expand <- function(b, blk, theta) {
  "c" <- RTMB::ADoverload("c")
  a <- blk[["aux"]]
  n <- a[["n"]]
  s <- exp(theta[1L])
  th <- theta[-1L]
  cf <- latres_coefs(th, a)
  switch(a[["kind"]],
    ar1 = s * b,
    ma1 = s * (b[seq_len(n)] + cf$ma * b[a[["prev"]]]),
    arma1 = s * (b[seq_len(n)] + cf$ma * b[a[["prev"]]]),
    condar = if (isTRUE(a[["on_response"]])) {
      s * b
    } else {
      s * (b + latres_lag_sum(b, cf$ar, a))
    },
    cosy = s * (sqrt(1 - cf$rho) * b[seq_len(n)] + sqrt(cf$rho) * b[a[["u"]]]),
    unstr = {
      d <- a[["d"]]
      Lv <- as.vector(us_chol_L(th, d))
      out <- 0
      for (k in seq_len(d)) {
        state <- (a[["g_of"]] - 1L) * d + k
        out <- out + Lv[(k - 1L) * d + a[["J"]]] * b[state]
      }
      s * out
    })
}

#' `sum_j ar_j z_{t-j}` over the rows, each lag one gather of `z` padded with
#' a zero for a row with no such lag.
#'
#' @noRd
latres_lag_sum <- function(z, ar, a) {
  "c" <- RTMB::ADoverload("c")
  zp <- c(z, 0 * z[1L])
  out <- 0 * z
  for (j in seq_along(ar)) out <- out + ar[j] * zp[a[["lag_idx"]][[j]]]
  out
}

#' The part of a cov = FALSE latent residual that brms adds to `mu` after the
#' inverse link, or `NULL` for a term with none. The objective and the post-fit
#' readers of `mu` take it from here, so they cannot disagree.
#'
#' @noRd
latres_shift <- function(b, blk, theta) {
  a <- blk[["aux"]]
  if (!isTRUE(a[["on_response"]])) return(NULL)
  cf <- latres_coefs(theta[-1L], a)
  exp(theta[1L]) * latres_lag_sum(b, cf$ar, a)
}

#' d err / d b as triplets at the block's theta, for `rr_jacobians()`. The
#' arithmetic of `latres_expand()`, with numbers.
#'
#' @noRd
latres_jacobian <- function(blk, theta) {
  a <- blk[["aux"]]
  n <- a[["n"]]
  s <- exp(theta[1L])
  cf <- latres_coefs(theta[-1L], a)
  rows <- seq_len(n)
  switch(a[["kind"]],
    ar1 = list(i = rows, j = rows, x = rep(s, n)),
    ma1 = ,
    arma1 = list(i = c(rows, rows), j = c(rows, a[["prev"]]),
                 x = c(rep(s, n), rep(s * cf$ma, n))),
    condar = {
      ii <- rows
      jj <- rows
      xx <- rep(s, n)
      if (!isTRUE(a[["on_response"]])) {
        for (k in seq_along(cf$ar)) {
          has <- which(a[["lag_idx"]][[k]] <= n)
          ii <- c(ii, has)
          jj <- c(jj, a[["lag_idx"]][[k]][has])
          xx <- c(xx, rep(s * cf$ar[k], length(has)))
        }
      }
      list(i = ii, j = jj, x = xx)
    },
    cosy = list(i = c(rows, rows), j = c(rows, a[["u"]]),
                x = c(rep(s * sqrt(1 - cf$rho), n),
                      rep(s * sqrt(cf$rho), n))),
    unstr = {
      d <- a[["d"]]
      L <- us_chol_L(theta[-1L], d)
      ii <- integer(0)
      jj <- integer(0)
      xx <- numeric(0)
      for (k in seq_len(d)) {
        hit <- which(a[["J"]] >= k)
        ii <- c(ii, hit)
        jj <- c(jj, (a[["g_of"]][hit] - 1L) * d + k)
        xx <- c(xx, s * L[cbind(a[["J"]][hit], k)])
      }
      list(i = ii, j = jj, x = xx)
    })
}

#' One draw of the block's whole `b` from its law at `theta`.
#'
#' @noRd
latres_draw <- function(theta, blk) {
  a <- blk[["aux"]]
  z <- stats::rnorm(a[["m"]])
  if (!a[["kind"]] %in% c("ar1", "arma1")) return(z)
  t1 <- theta[2L]
  phi <- autocor_pacf(t1)
  b <- numeric(a[["m"]])
  b[a[["start"]]] <- z[a[["start"]]] * sqrt(1 + t1 * t1)
  # a row's predecessor is in an earlier generation
  for (r in a[["steps"]] %||% latres_steps(a)) {
    b[r] <- phi * b[a[["prev"]][r]] + z[r]
  }
  b
}

#' The rows of a chain by generation: a row's generation is its position in
#' its group, and only `ar1` has first rows that the start fills.
#'
#' @noRd
latres_steps <- function(a) {
  pr <- a[["pos_rows"]]
  if (identical(a[["kind"]], "ar1")) pr[-1L] else pr
}

#' Natural-scale coefficients of a latent term under brms's names. `th` is the
#' block's theta after log sderr.
#'
#' @noRd
latres_natural_coefs <- function(th, a) {
  ac <- a[["ac"]]
  cf <- latres_coefs(th, a)
  switch(a[["kind"]],
    ar1 = c("ar[1]" = cf$ar),
    ma1 = c("ma[1]" = cf$ma),
    arma1 = c("ar[1]" = cf$ar, "ma[1]" = cf$ma),
    condar = stats::setNames(cf$ar, sprintf("ar[%d]", seq_along(cf$ar))),
    cosy = c(cosy = cf$rho),
    unstr = {
      C <- us_chol_cor(th, a[["d"]])
      pr <- which(lower.tri(C), arr.ind = TRUE)
      stats::setNames(C[lower.tri(C)],
                      paste0("cortime__", ac[["time_levels"]][pr[, 2L]], "__",
                             ac[["time_levels"]][pr[, 1L]]))
    })
}

#' A latent block's rows for `confint_varcorr()` and the summary: `sderr`
#' first, then the term's coefficients, each on the scale its interval is built
#' on, with the Jacobian by central differences. The chart of each is a smooth
#' function of its own entries, and this runs per report.
#'
#' @noRd
latres_natural <- function(theta, blk) {
  a <- blk[["aux"]]
  cor_like <- a[["kind"]] != "cosy"
  th <- theta[-1L]
  val <- function(tt) {
    v <- unname(latres_natural_coefs(tt, a))
    if (a[["kind"]] == "cosy") return(stats::qlogis(v))
    atanh(pmin(pmax(v, -varcorr_cor_clamp), varcorr_cor_clamp))
  }
  nm <- names(latres_natural_coefs(th, a))
  est <- val(th)
  J <- vapply(seq_along(th), function(i) {
    h <- 1e-5 * max(abs(th[i]), 1)
    tp <- th
    tp[i] <- tp[i] + h
    tn <- th
    tn[i] <- tn[i] - h
    (val(tp) - val(tn)) / (2 * h)
  }, numeric(length(est)))
  J <- matrix(J, nrow = length(est))
  list(term = c("sderr", nm),
       type = c("sd", rep(if (cor_like) "cor" else "prop", length(nm))),
       est = c(theta[1L], est),
       grad = rbind(c(1, numeric(length(th))), cbind(0, J)))
}

#' The latent residual block of each response, by response name.
#'
#' @noRd
frame_latres <- function(frame) {
  out <- list()
  for (bk in frame[["re_blocks"]] %||% list()) {
    if (identical(bk[["covstruct"]], "latres")) {
      out[[bk[["aux"]][["ac"]][["resp"]]]] <- bk
    }
  }
  out
}

#' Whether a term builds no dense matrix over its time levels, so the cap on
#' their number does not apply: every latent chain but `unstr()`.
#'
#' @noRd
autocor_block_sparse <- function(ac) {
  isTRUE(ac[["latent"]]) && !identical(ac[["struct"]], "unstr")
}

#' Every autocorrelation term of a frame that a prior class reaches, by
#' response: `frame$autocor`, and each latent residual's term with `theta_idx`
#' pointing past its log sderr. The prior table and the resolver read this, so
#' the two kinds of term share one set of class rules.
#'
#' @noRd
frame_ac_terms <- function(frame) {
  out <- frame[["autocor"]] %||% list()
  for (bk in frame_latres(frame)) {
    ac <- bk[["aux"]][["ac"]]
    ac[["theta_idx"]] <- bk[["theta_idx"]][-1L]
    ac[["theta_comp"]] <- "theta"
    out[[ac[["resp"]]]] <- ac
  }
  out
}

#' The template component a term's parameters live in.
#'
#' @noRd
ac_comp <- function(ac) ac[["theta_comp"]] %||% "thetaac"

#' The post-link part of a latent term at a fit's estimates.
#'
#' @noRd
latres_shift_of <- function(fit, bk) {
  latres_shift(fit$estimates[["b"]][bk[["b_idx"]]], bk,
               fit$estimates[["theta"]][bk[["theta_idx"]]])
}

#' brms's refusal of new rows under a cov = FALSE latent residual, whose iid
#' latents have no law that places a new row after an old one. A grid path has
#' dropped the term by then (`autocor_cond_strip()`), as brms's
#' `incl_autocor = FALSE` does.
#'
#' @noRd
latres_refuse_newdata <- function(fit, newdata) {
  if (is.null(newdata) || isTRUE(fit$frame[["latres_off"]])) {
    return(invisible(NULL))
  }
  for (bk in frame_latres(fit$frame)) {
    if (identical(bk[["aux"]][["kind"]], "condar")) {
      frm_stop("Cannot predict new latent residuals when using cov = FALSE ",
               "in autocor terms. brms refuses the same call: ",
               bk[["aux"]][["ac"]][["label"]], " has one independent latent ",
               "per fitted row, with no law for a new row. Fit cov = TRUE, ",
               "or predict the rows the model was fitted on", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' What a replicate needs to draw the latent residuals of new rows: the term's
#' chain over `newdata`'s own groups and row order, as brms reads them from the
#' new data. Only `unstr()` reads the time levels, and a level the fit never
#' saw has no row or column in its matrix.
#'
#' @noRd
latres_new_spec <- function(fit, rspec, newdata) {
  if (is.null(newdata) || isTRUE(fit$frame[["latres_off"]])) return(NULL)
  bk <- frame_latres(fit$frame)[[rspec$resp_name]]
  if (is.null(bk)) return(NULL)
  latres_refuse_newdata(fit, newdata)
  ac <- bk[["aux"]][["ac"]]
  n <- nrow(newdata)
  env <- rspec$formula_env
  # R would otherwise find t() or g() in the search path and hand a function
  # to the term's own checks
  want <- c(all.vars(ac[["time_expr"]]), ac[["gr_vars"]])
  miss <- setdiff(want, names(newdata))
  if (length(miss)) {
    frm_stop(ac[["label"]], " draws new latent residuals by group and time, ",
             "so newdata needs ", paste0("`", miss, "`", collapse = ", "),
             call. = FALSE)
  }
  gv <- if (is.null(ac[["gr_expr"]])) NULL else {
    autocor_gr_value(ac, newdata, env)
  }
  gidx <- if (is.null(gv)) rep(1L, n) else as.integer(factor(gv))
  ti <- autocor_time_index(ac, newdata, gidx, env, n)
  tidx <- ti$idx
  if (identical(ac[["struct"]], "unstr")) {
    tidx <- match(ti$levels[ti$idx], ac[["time_levels"]])
    if (anyNA(tidx)) {
      frm_stop("unstr(): time point '", ti$levels[ti$idx][is.na(tidx)][1L],
               "' of newdata is not one the model was fitted on, so the ",
               "correlation matrix has no row for it", call. = FALSE)
    }
  }
  aux <- latres_aux(ac, gidx, tidx, n, FALSE)
  if (aux[["kind"]] %in% c("ar1", "arma1")) aux[["steps"]] <- latres_steps(aux)
  list(theta_idx = bk[["theta_idx"]], blk = list(aux = aux))
}

#' One draw of the latent residuals of the new rows at a replicate's theta.
#'
#' @noRd
latres_new_draw <- function(fs, lat) {
  th <- fs$estimates[["theta"]][lat[["theta_idx"]]]
  as.numeric(latres_expand(latres_draw(th, lat[["blk"]]), lat[["blk"]], th))
}

#' A latent term's parameters under brms's variable names at the block's
#' theta: `ar[1]`, `ma[1]`, `cosy`, `cortime__t1__t2`, `sderr`, the last with
#' brms's response suffix in a multivariate model, as `sigma_y` carries it.
#'
#' @noRd
latres_variables <- function(fit, bk, theta) {
  a <- bk[["aux"]]
  nat <- latres_natural_coefs(theta[-1L], a)
  if (length(fit$spec$responses) > 1L) {
    rs <- brms_stan_name(a[["ac"]][["resp"]])
    names(nat) <- sub("^(ar|ma|cosy|cortime)", paste0("\\1_", rs), names(nat))
    return(c(nat, stats::setNames(exp(theta[1L]), paste0("sderr_", rs))))
  }
  c(nat, sderr = exp(theta[1L]))
}

#' The rows the summary reports for the residual and latent terms, in brms's
#' order: a term's coefficients, then its `sderr`.
#'
#' @noRd
ac_report_rows <- function(fit) {
  tr <- autocor_trans_rows(fit)
  lrs <- frame_latres(fit$frame)
  if (!length(lrs)) return(tr)
  sdr <- sdr_of(fit)
  Vfull <- sdr$cov.fixed
  pos <- which(rownames(Vfull) == "theta")
  th <- fit$estimates[["theta"]]
  rows <- list()
  for (bk in lrs) {
    ti <- bk[["theta_idx"]]
    nr <- cs_natural_rows(bk, th[ti])
    V <- Vfull[pos[ti], pos[ti], drop = FALSE]
    se <- vapply(seq_along(nr$term), function(k) {
      cs_natural_se(nr$grad[k, ], V)
    }, 0)
    se[!is.finite(se) | latres_block_limit(fit, bk, V)$na] <- NA_real_
    o <- c(seq_along(nr$term)[-1L], 1L)
    rows[[length(rows) + 1L]] <- data.frame(
      block = bk[["term_label"]], term = nr$term[o], type = nr$type[o],
      est_t = nr$est[o], se_t = se[o], stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, c(list(tr), rows))
  rownames(out) <- NULL
  out
}

#' A fit whose predictor of `mu` on new rows carries one fresh draw of the
#' latent residuals, for a method that evaluates a whole draw at a time through
#' the fit's predictor (frmtmb.sample's `posterior_epred()` and
#' `posterior_predict()`). `lat` is `latres_new_spec()`'s, or `NULL`, and the
#' fit comes back as it was. The shift is `lp_eta_design()`'s own `eta_shift`,
#' so no reader of the predictor needs to know.
#'
#' @noRd
latres_shifted <- function(fit, lat, resp) {
  if (is.null(lat) || isTRUE(fit$frame[["latres_off"]])) return(fit)
  fit[["eta_shift"]] <- list(key = linpred_key(resp, "mu"),
                             h = latres_new_draw(fit, lat))
  fit
}

#' A fit that a grid path reads, with its latent residuals off: new rows draw
#' none and are not refused, brms's `incl_autocor = FALSE`.
#'
#' @noRd
autocor_cond_strip_latent <- function(fit) {
  fit$frame[["latres_off"]] <- TRUE
  fit
}

# ------------------------------------------------- what a fit tells the user

# On a true ridge the SE of log sderr was 1.4 to 192 or NaN; on a
# well-determined strong AR, 0.03 to 0.07 (dev/review-latres-m35b.R, -m3c.R). So
# nearness to 1 alone says nothing: ar 0.99922 on 2000 rows is a number the data
# give.
latres_se_flat <- 1

#' Which of a latent block's rows (`sderr`, then its coefficients) sit at a
#' limit of their range, where a standard error would read as knowledge the data
#' do not give. `se` is the standard errors of the block's `theta`, NaN when not
#' known. A row is at its limit when the value is near the end of its range and
#' the fit's numbers say that end is not determined: the SE of log `sderr` is
#' flat, or the row's own SE is not finite. `sderr` is near 0 relative to
#' `scale`, the response's sd on an identity link, so a fit in small units is
#' not at a limit.
#'
#' `at` names the rows at a limit. `na` also holds the term's coefficients when
#' `sderr` is at 0, which the data do not show but which are at no limit.
#'
#' @noRd
latres_at_limit <- function(theta, a, se = NULL, scale = 1) {
  tol <- 1e-3
  se <- se %||% rep(NaN, length(theta))
  flat <- !is.finite(se[1L]) || se[1L] > latres_se_flat
  cf <- latres_natural_coefs(theta[-1L], a)
  near <- if (a[["kind"]] == "unstr") {
    logical(length(cf))
  } else if (identical(a[["kind"]], "cosy")) {
    cf < tol | cf > 1 - tol
  } else {
    abs(cf) > 1 - tol
  }
  at_cf <- near & (flat | !is.finite(se[-1L]))
  # ar and ma cancel to white noise, where each is free along the other
  if (identical(a[["kind"]], "arma1") && any(at_cf) &&
        abs(cf[1L] + cf[2L]) < 0.01) {
    at_cf[] <- TRUE
  }
  zero <- exp(theta[1L]) < tol * scale && flat
  # a coefficient on the ridge takes sderr with it, as sderr -> 0 there
  ridge <- a[["kind"]] %in% c("ar1", "arma1") && at_cf[1L] && flat
  at <- c(zero || ridge, at_cf)
  na <- at
  if (zero) na[] <- TRUE
  list(at = unname(at), na = unname(na), se = unname(se))
}

#' `latres_at_limit()` for a fit's block. `V` is the covariance of the block's
#' theta; the sdreport is read only when a value is near an end, so a fit inside
#' its range pays nothing.
#'
#' @noRd
latres_block_limit <- function(fit, bk, V = NULL) {
  th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
  a <- bk[["aux"]]
  r <- a[["ac"]][["resp"]]
  link <- fit$spec$responses[[r]]$family[["links"]][["mu"]][["name"]] %||%
    "identity"
  scale <- 1
  if (identical(link, "identity")) {
    scale <- stats::sd(as.numeric(fit$frame[["y"]][[r]]), na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) scale <- 1
  }
  blind <- latres_at_limit(th, a, NULL, scale)
  if (!any(blind$at)) return(blind)
  if (is.null(V)) {
    V <- tryCatch({
      Vf <- sdr_of(fit)$cov.fixed
      p <- which(rownames(Vf) == "theta")[bk[["theta_idx"]]]
      Vf[p, p, drop = FALSE]
    }, error = keep_missing(NULL))
  }
  if (is.null(V)) return(blind)
  v <- diag(as.matrix(V))
  se <- sqrt(pmax(v, 0))
  se[!is.finite(v)] <- NaN
  latres_at_limit(th, a, se, scale)
}

#' The sd of one row's latent residual, at the block's theta: `sderr`'s scale of
#' a stationary process, which is one number whatever the row's position.
#'
#' @noRd
latres_marginal_sd <- function(theta, a) {
  s <- exp(theta[1L])
  cf <- latres_coefs(theta[-1L], a)
  switch(a[["kind"]],
    ar1 = s * sqrt(1 + theta[2L]^2),
    ma1 = s * sqrt(1 + cf$ma^2),
    arma1 = s * sqrt((1 + cf$ma^2 + 2 * cf$ar * cf$ma) * (1 + theta[2L]^2)),
    s)
}

#' brms's names for a latent block's `theta`, in order, for a message.
#'
#' @noRd
latres_theta_names <- function(bk) {
  a <- bk[["aux"]]
  np <- length(bk[["theta_idx"]]) - 1L
  c("sderr", names(latres_natural_coefs(numeric(np), a)))
}

#' Tell the user what a latent block's limit means, at the level of a singular
#' fit (`se_boundary_act()`).
#'
#' @noRd
latres_limit_report <- function(fit, control) {
  for (bk in frame_latres(fit$frame)) {
    th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
    a <- bk[["aux"]]
    state <- latres_block_limit(fit, bk)
    lim <- state$at
    if (!any(lim)) next
    nm <- latres_theta_names(bk)
    val <- c(exp(th[1L]), latres_natural_coefs(th[-1L], a))
    # a stop short of the maximum silences the generic check, so a row beside
    # the limit with no standard error is said here
    other <- !state$na & !is.finite(state$se)
    # beside a chain at its limit a row can read NA with a finite, large standard
    # error, and no sentence named it
    unsaid <- state$na & !lim
    # a coefficient at its limit reads as the limit, which is what it means
    shown <- ifelse(seq_along(val) == 1L, signif(val, 2L), round(val, 3L))
    chain <- a[["kind"]] %in% c("ar1", "arma1") && lim[1L] && lim[2L]
    cancel <- identical(a[["kind"]], "arma1") && lim[2L] && lim[3L] &&
      val[2L] * val[3L] < 0
    msg <- paste0(
      "Boundary (singular) fit: ",
      paste0(nm[lim], " = ", shown[lim], collapse = ", "),
      " in ", bk[["term_label"]],
      if (sum(lim) > 1L) " are at the limits of their ranges. " else {
        " is at the limit of its range. "
      },
      if (cancel) {
        paste0("Their roots cancel, so the data show no correlation in the ",
               "latent and neither coefficient is known. ")
      } else if (chain && val[2L] > 0) {
        paste0("The data read as cosy() or (1 | group): the marginal sd of ",
               "the latent, ", format(latres_marginal_sd(th, a), digits = 3),
               ", is the number they identify, and neither of the two is ",
               "known. ")
      } else if (chain) {
        paste0("The data read as one size per group whose sign alternates by ",
               "row, with a marginal sd of ",
               format(latres_marginal_sd(th, a), digits = 3),
               ", and neither of the two is known. ")
      } else if (lim[1L]) {
        paste0("The data show no latent residual, so the term's coefficients ",
               "are not known. ")
      },
      if (any(other)) {
        paste0(paste(nm[other], collapse = ", "),
               " has no standard error either. ")
      },
      if (chain && any(unsaid)) {
        paste0(paste(nm[unsaid], collapse = ", "),
               " is not known there either. ")
      },
      "Est.Error and the interval are NA in summary() and confint_varcorr(); ",
      "diagnose() lists the limit")
    se_boundary_act(control$check_se %||% "warning", msg)
  }
  invisible(NULL)
}

#' The rows `latres_at_limit()` marks, in the shape of `diagnose()$singular`.
#'
#' @noRd
latres_singular_rows <- function(fit) {
  out <- NULL
  for (bk in frame_latres(fit$frame)) {
    th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
    a <- bk[["aux"]]
    lim <- latres_block_limit(fit, bk)$at
    if (!any(lim)) next
    val <- c(exp(th[1L]), latres_natural_coefs(th[-1L], a))
    out <- rbind(out, data.frame(block = bk[["term_label"]],
                                 term = latres_theta_names(bk)[lim],
                                 value = unname(val[lim])))
  }
  out
}

#' The derivatives of each row's log density in its linear predictor at the fit,
#' by differences of the family's own density, so a family needs no slot for
#' them. The step is a share of the latent's sd, and `d2 * var` and
#' `d3 * sd^3` do not read the response's units.
#'
#' @noRd
latres_row_derivs <- function(fit, bk) {
  a <- bk[["aux"]]
  r <- a[["ac"]][["resp"]]
  fam <- fit$spec$responses[[r]]$family
  dp <- eval_dpars(fit)[[r]]
  lk <- fit$frame[["linpreds"]][[linpred_key(r, "mu")]][["link"]]
  eta <- lk$linkfun(as.numeric(dp[["mu"]]))
  y <- fit$frame[["y"]][[r]]
  av <- fit$frame[["aterm_values"]][[r]]
  ex <- fit_extras(fit)
  ll <- function(e) {
    dp[["mu"]] <- lk$linkinv(e)
    as.numeric(row_lpdf(fam, y, y, dp, av, ex))
  }
  sds <- latres_marginal_sd(fit$estimates[["theta"]][bk[["theta_idx"]]], a)
  h <- 1e-3 * sds
  l0 <- ll(eta)
  lp <- ll(eta + h)
  ln <- ll(eta - h)
  list(d2 = -(lp - 2 * l0 + ln) / h^2,
       d3 = (ll(eta + 2 * h) - 2 * lp + 2 * ln - ll(eta - 2 * h)) / (2 * h^3),
       var = sds^2)
}

#' The curvature of each row's log density in its linear predictor.
#'
#' @noRd
latres_row_info <- function(fit, bk) pmax(latres_row_derivs(fit, bk)$d2, 0)

# M1's condition, chosen on 57 designs against the exact MLE (dev/latres-r3.R
# and dev/latres-findings.md, "R3"). No number here has the response's units: a
# row's share of its latent is d2 v / (1 + d2 v), and the third derivative is
# taken at the posterior sd of the latent. A log density quadratic in eta is
# exact whatever the share: lognormal reads 6e-9 here and a poisson of mean 1,
# 0.09.
latres_cubic_floor <- 1e-3
latres_share_floor <- 0.18
latres_share_disp <- 0.3
latres_info_disp <- 2.5
latres_trials_floor <- 3

#' The first-order Laplace correction of a fit at outer parameters `par`:
#' `-1/8 sum(h4 v^2) + 1/8 sum(h3 h3 v v v) + 1/12 sum(h3 h3 v^3)`, with `h3`, `h4`
#' the row's derivatives of `-log density` in `eta` and `v = A H^-1 A'` the
#' covariance of the latents at the inner mode. Read at the fit it is the exact
#' value's distance from Laplace: median 0.340 to 0.012 on designs no side had
#' seen (dev/reviews/2026-11-02-latres.md, RF1).
#'
#' @noRd
latres_tk <- function(fit, par, bk, map) {
  env <- fit$obj$env
  if (!is.finite(suppressWarnings(fit$obj$fn(par)))) return(NA_real_)
  full <- env$last.par
  at <- fit_set_outer(fit, par, map)
  at$estimates[["b"]] <- as.numeric(full[env$random])
  a <- bk[["aux"]]
  r <- a[["ac"]][["resp"]]
  fam <- at$spec$responses[[r]]$family
  dp <- eval_dpars(at)[[r]]
  lk <- at$frame[["linpreds"]][[linpred_key(r, "mu")]][["link"]]
  eta <- lk$linkfun(as.numeric(dp[["mu"]]))
  y <- at$frame[["y"]][[r]]
  av <- at$frame[["aterm_values"]][[r]]
  ex <- fit_extras(at)
  ll <- function(e) {
    dp[["mu"]] <- lk$linkinv(e)
    as.numeric(row_lpdf(fam, y, y, dp, av, ex))
  }
  h <- 0.02 * latres_marginal_sd(at$estimates[["theta"]][bk[["theta_idx"]]], a)
  l0 <- ll(eta)
  lp <- ll(eta + h)
  ln <- ll(eta - h)
  lp2 <- ll(eta + 2 * h)
  ln2 <- ll(eta - 2 * h)
  h3 <- -(lp2 - 2 * lp + 2 * ln - ln2) / (2 * h^3)
  h4 <- -(lp2 - 4 * lp + 6 * l0 - 4 * ln + ln2) / h^4
  pj <- latres_jacobian(bk, at$estimates[["theta"]][bk[["theta_idx"]]])
  A <- Matrix::sparseMatrix(i = pj$i, j = pj$j, x = pj$x,
                            dims = c(a[["n"]], length(env$random)))
  H <- env$spHess(full, random = TRUE)
  # TMB fills one template's values, and Matrix keeps a factor in its slot, so
  # Cholesky() on it hands back the factor of an earlier point: a residual of
  # 0.037 here, where C read 0.97 for 1.10
  H@factors <- list()
  Ch <- Matrix::Cholesky(H, LDL = FALSE)
  v <- as.matrix(A %*% Matrix::solve(Ch, Matrix::t(A)))
  dv <- diag(v)
  -sum(h4 * dv^2) / 8 + sum(outer(h3 * dv, h3 * dv) * v) / 8 +
    sum(outer(h3, h3) * v^3) / 12
}

#' How far that correction moves log `sderr` and `ar[1]`: `V grad(C)`, the
#' gradient by differences of 0.05 in each outer parameter but `beta`. `Inf` where
#' the fit has no covariance to read it by. The tape's last point goes back on
#' exit, or a reader of `env$last.par` reads the last step: 1.5 nats on "the
#' readers agree with the tape".
#'
#' @noRd
latres_shift_tk <- function(fit, bk) {
  par <- fit$opt$par
  map <- outer_par_map(fit)
  env <- fit$obj$env
  lp0 <- env$last.par
  lb0 <- env$last.par.best
  on.exit({
    env$last.par <- lp0
    env$last.par.best <- lb0
  })
  V <- tryCatch(as.matrix(sdr_of(fit)$cov.fixed), error = function(e) NULL)
  if (is.null(V) || nrow(V) != length(par)) return(Inf)
  g <- numeric(length(par))
  for (k in which(names(par) != "beta")) {
    pp <- par
    pp[k] <- pp[k] + 0.05
    pm <- par
    pm[k] <- pm[k] - 0.05
    g[k] <- (latres_tk(fit, pp, bk, map) - latres_tk(fit, pm, bk, map)) / 0.1
  }
  g[!is.finite(g)] <- 0
  th <- which(names(par) == "theta")[bk[["theta_idx"]][1:2]]
  if (anyNA(th) || !all(is.finite(V[th, ]))) return(Inf)
  d <- as.numeric(V[th, ] %*% g)
  t2 <- par[th[2]]
  max(abs(d[1]),
      abs((t2 + d[2]) / sqrt(1 + (t2 + d[2])^2) - t2 / sqrt(1 + t2^2)))
}

#' Why a fit's Laplace approximation is biased, in a clause, or `NULL`.
#'
#' @noRd
latres_laplace_reason <- function(fit, bk) {
  r <- bk[["aux"]][["ac"]][["resp"]]
  fam <- fit$spec$responses[[r]]$family
  # above the limit the correction's dense v is the cost the limit avoids, so
  # the row-level condition stays there
  if (fit$frame[["n_obs"]] > getOption("frmtmb.latres_se_rows", 3000L)) {
    return(latres_share_reason(fit, bk, fam, r))
  }
  # its error is not local: the shift is 0.02 to 0.06 with the exact MLE far off
  # along a flat ridge, so the correction cannot read it
  if (identical(fam[["type"]], "ordinal")) return("each row is one category")
  s <- tryCatch(latres_shift_tk(fit, bk), error = function(e) Inf)
  if (s <= 0.15) return(NULL)
  if (!is.finite(s)) {
    return("the fit has no curvature to read its first-order correction by")
  }
  paste0("the first-order correction moves log(sderr) or ar[1] by ",
         format(s, digits = 2), ", above 0.15")
}

#' The row-level condition above the row limit.
#'
#' @noRd
latres_share_reason <- function(fit, bk, fam, r) {
  g <- latres_row_derivs(fit, bk)
  i2 <- pmax(g$d2, 0)
  cubic <- mean(abs(g$d3) / (i2 + 1 / g$var)^1.5, na.rm = TRUE)
  if (!isTRUE(cubic >= latres_cubic_floor)) return(NULL)
  tr <- fit$frame[["aterm_values"]][[r]][["trials"]]
  if (identical(fam[["family"]], "bernoulli")) return("each row is binary")
  if (identical(fam[["type"]], "ordinal")) return("each row is one category")
  if (!is.null(tr) && stats::median(tr) <= latres_trials_floor) {
    return(paste0("a row has at most ", stats::median(tr), " trials"))
  }
  share <- mean(i2 * g$var / (1 + i2 * g$var), na.rm = TRUE)
  shown <- format(share, digits = 2)
  if (share < latres_share_floor) {
    return(paste0("a row determines ", shown, " of its latent, below ",
                  latres_share_floor))
  }
  if (!length(latres_log_dpars(fit$spec$responses[[r]]$family))) return(NULL)
  link <- fam[["links"]][["mu"]][["name"]] %||% "identity"
  plain <- !identical(link, "identity") &&
    is.null(latres_mu_range(fam, FALSE))
  if (share < latres_share_disp ||
        (plain && mean(i2, na.rm = TRUE) < latres_info_disp)) {
    return(paste0("the latent trades off against the family's dispersion, ",
                  "with a row determining ", shown, " of its latent"))
  }
  NULL
}

#' The warning of a fit whose rows are too thin for the Laplace approximation.
#' It is about the estimates, so `check_se` does not reach it; its class does:
#' `suppressWarnings(expr, classes = "frmtmb_latres_laplace")`.
#'
#' @noRd
latres_laplace_report <- function(fit) {
  for (bk in frame_latres(fit$frame)) {
    why <- tryCatch(latres_laplace_reason(fit, bk), error = keep_missing(NULL))
    if (is.null(why)) next
    frm_warning(
      "The Laplace approximation of ", bk[["term_label"]], " is biased on ",
      "this response: ", why, ". On bernoulli, 20 groups of 10, sderr reads ",
      "0.39 where the exact MLE is 1.10, and on beta_binomial() 0.26 where it ",
      "is 0.65, with code 0 and no other sign. frm_sample() samples the same ",
      "model and agrees with brms; ?frmtmb-autocor has the table",
      call. = FALSE,
      class = "frmtmb_latres_laplace")
  }
  invisible(NULL)
}

#' The response's mean over a gaussian latent on the linear predictor of `mu`,
#' row by row, by Gauss-Hermite: `E[mean(linkinv(eta + Z))]`, `Z ~ N(0, var)`.
#'
#' `var` is the row's total latent variance, whatever carries it (a latent
#' residual on a new row, a `me()` value on a new row), so a model with several
#' has one rule: their variances add, as the terms are independent given the
#' parameters. A row at 0 reads the plain mean. `dpars` are the family's dpars
#' on these rows, and `mu` among them is replaced at each node, so a
#' zero-inflated or ordinal mean, or a truncated one, integrates as a whole.
#' `aterms` are the rows' addition terms. `mean_fn` is `response_mean()`, or
#' another reader of the same arguments for a caller that wants `mu` alone.
#'
#' brms's `posterior_epred()` reads one draw of the latent per posterior draw,
#' whose mean is this integral; at the latent's 0 a log link is off by a factor
#' of 1.95 at a marginal sd of 1.15 and 14.8 at 2.31.
#'
#' Accuracy of 41 nodes: under 1e-8 relative on a log link up to an sd of 4, and
#' 1.5e-7 on a logit at 2.31, whose poles sit near the real line.
#'
#' @noRd
response_mean_gh <- function(eta, var, fam, dpars, aterms = list(),
                             mean_fn = response_mean) {
  sds <- sqrt(rep_len(var, length(eta)))
  linkinv <- fam[["links"]][["mu"]][["linkinv"]]
  gh <- gh_std_normal(41L)
  total <- 0
  lost <- 0
  for (k in seq_along(gh$x)) {
    dpars[["mu"]] <- linkinv(eta + sds * gh$x[k])
    v <- mean_fn(fam, dpars, aterms)
    bad <- !is.finite(v)
    if (any(bad)) {
      v[bad] <- 0
      lost <- lost + gh$w[k] * bad
    }
    total <- total + gh$w[k] * v
  }
  # a family's mean can fail at a node of weight 1e-20, a mu of 5000 under
  # trunc(ub = 100): that node says nothing, and a cell that lost real weight
  # has no answer, which NaN says
  out <- total / (1 - lost)
  out[lost >= latres_gh_lost] <- NaN
  out
}

#' `fitted()` on new rows of a fit with a latent residual: `response_mean_gh()`
#' at the row's marginal variance. `out` is the value at the latent's 0,
#' returned as it is where the integral is that value (a mean linear in an
#' identity link, or a grid path with the term dropped). The predictor reads a
#' fit's `eta_shift`, so `fit_fd_se()`'s step moves this mean too.
#'
#' @noRd
latres_fitted_mean <- function(fit, newdata, re_formula, resp, dpar,
                               allow_new_levels, out) {
  rn <- resp %||% names(fit$spec$responses)[1L]
  bk <- frame_latres(fit$frame)[[rn]]
  if (is.null(bk) || isTRUE(fit$frame[["latres_off"]])) return(out)
  rspec <- fit$spec$responses[[rn]]
  fam <- rspec$family
  link <- fam[["links"]][["mu"]][["name"]] %||% "identity"
  if (identical(link, "identity") && mean_is_mu(fam) && !has_trunc(rspec)) {
    return(out)
  }
  s <- latres_marginal_sd(fit$estimates[["theta"]][bk[["theta_idx"]]],
                          bk[["aux"]])
  eta <- as.numeric(frm_linpred(fit, newdata = newdata, type = "link",
                                dpar = "mu", resp = rn, re_formula = re_formula,
                                allow_new_levels = allow_new_levels))
  dp <- dpars_natural(fit, rspec, newdata, re_formula, allow_new_levels)
  reader <- if (identical(dpar, "mu")) {
    function(fam, dpars, aterms) dpars[["mu"]]
  } else if (length(dim(out))) {
    # category probabilities are no response_mean(): the fit's own path reads
    # them, at the node's predictor
    linkfun <- fam[["links"]][["mu"]][["linkfun"]]
    # eta holds fit_fd_se()'s chain step, and this list replaces it: without it
    # d out / d eta reads 0, which drops a gp()'s kriging variance
    key <- linpred_key(rn, "mu")
    h0 <- if (identical(fit[["eta_shift"]][["key"]], key)) {
      fit[["eta_shift"]][["h"]]
    } else {
      0
    }
    function(fam, dpars, aterms) {
      fit[["eta_shift"]] <- list(key = key,
                                 h = linkfun(dpars[["mu"]]) - eta + h0)
      frm_linpred(fit, newdata = newdata, type = "response", dpar = dpar,
                  resp = resp, re_formula = re_formula,
                  allow_new_levels = allow_new_levels)
    }
  } else {
    response_mean
  }
  total <- response_mean_gh(eta, s^2, fam, dp,
                            aterms_for_newdata(rspec, newdata),
                            mean_fn = reader)
  attributes(total) <- attributes(out)
  total
}

#' Whether a standard error of a row of this fit would solve the whole joint
#' precision densely, above the row count where that stops being a cost a method
#' can pay unasked: 6 to 12 s at 2000 rows, 80 s at 4000, 984 s and 491 MB of
#' saved cache at 8000, and 80 GB at 100000, where the fit itself takes 30 s.
#' `options(frmtmb.latres_se_rows = )` moves the limit.
#'
#' @noRd
latres_skip_se <- function(fit, newdata = NULL) {
  # rows of a grid read beta and theta alone, 0.02 to 0.31 s at 4000 rows
  if (!is.null(newdata) || !length(frame_latres(fit$frame))) return(FALSE)
  limit <- getOption("frmtmb.latres_se_rows", 3000L)
  if (fit$frame[["n_obs"]] <= limit) return(FALSE)
  frm_message("The standard errors of the fitted rows are skipped (NA): a ",
              "latent residual over ", fit$frame[["n_obs"]], " rows would ",
              "need the whole joint covariance as a dense matrix, about ",
              format(8 * length(fit$estimates[["b"]])^2 / 1e9, digits = 2),
              " GB. options(frmtmb.latres_se_rows = ) moves the limit of ",
              limit, " rows; frm_linpred() gives the values at no cost",
              class = "frmtmb_latres_se_skipped")
  TRUE
}

#' `predict()`'s parameter draw holds a latent block's own parameters at their
#' estimates. Its `err` is `exp(theta_1)` times unit-scale latents, so a draw of
#' `theta_1` with a standard error of 10 multiplies them by `exp(30)` before the
#' inverse link: 3.7e304 on a poisson fit of 60 rows whose largest count is 10.
#' A random effect holds `b` on its natural scale, and its draw cannot do this.
#' The latents are still drawn, given the estimates.
#'
#' @noRd
latres_held_pos <- function(fit, map) {
  bks <- frame_latres(fit$frame)
  if (!length(bks)) return(integer(0))
  th <- which(map$comp == "theta")
  th[unlist(lapply(bks, `[[`, "theta_idx"), use.names = FALSE)]
}

#' The forms of a latent residual `frm()` refuses before the fit, where
#' `frm_sample()` fits them. `frm_sample()` builds its objective with
#' `dry_run = "objective"`, which does not reach this.
#'
#' `cov = FALSE` after the inverse link: on benign designs 54 of 120 fits read
#' true, and no number from the fit says which: Beta at a margin of 6.6 reads
#' true and Gamma(log) at 23.6 reads ar -0.86 for 0.47
#' (dev/review-latres-m4b.R). A link that does not hold the mean under `cov =
#' TRUE` has such a number, and `latres_refuse_margin()` reads it after the fit.
#'
#' @noRd
latres_refuse_laplace <- function(spec, frame) {
  for (bk in frame_latres(frame)) {
    a <- bk[["aux"]]
    resp <- spec$responses[[a[["ac"]][["resp"]]]]
    fam <- resp$family
    lab <- a[["ac"]][["label"]]
    # two latent vectors in one inner problem, which no measurement has read
    if (isTRUE(resp$aterms[["mi"]])) {
      frm_stop(lab, ": mi() on the response cannot be combined with a latent ",
               "residual on frm(). The imputed values and the latents would ",
               "share one inner problem, and that Laplace approximation has ",
               "not been measured against an exact value. frm_sample() ",
               "samples the model", call. = FALSE)
    }
    if (!isTRUE(a[["on_response"]])) next
    frm_stop(
      lab, " on ", fam[["family"]], "(link = '",
      fam[["links"]][["mu"]][["name"]] %||% "identity", "'): the AR part adds ",
      "to the mean after the inverse link, and no number from the fit says ",
      "whether the Laplace approximation reads that likelihood. Against the ",
      "exact marginal, 54 of 120 benign fits read true: Beta at phi 50 does, ",
      "and Gamma(log) at a mean of 20 reads ar -0.86 for 0.47. frm_sample() ",
      "samples this model, and brms's own chains struggle on it (up to 3507 ",
      "of 4000 divergent). cov = TRUE keeps the term on the link scale, where ",
      "frm() reads it", call. = FALSE)
  }
  invisible(NULL)
}

#' Which of the fit's named parameters `latres_limit_report()` will speak for.
#'
#' @noRd
latres_limit_par <- function(fit, nms) {
  out <- logical(length(nms))
  for (bk in frame_latres(fit$frame)) {
    ti <- bk[["theta_idx"]]
    said <- ti[latres_block_limit(fit, bk)$na]
    if (length(said)) out <- out | nms %in% paste0("theta_", said)
  }
  out
}

#' A message with brms's names in place of a latent block's `theta_k`, which
#' name nothing a user wrote.
#'
#' @noRd
latres_par_text <- function(fit, s) {
  for (bk in frame_latres(fit$frame)) {
    nm <- latres_theta_names(bk)
    ti <- bk[["theta_idx"]]
    for (k in seq_along(ti)) {
      s <- gsub(paste0("\\btheta_", ti[k], "\\b"),
                paste0(nm[k], " (in ", bk[["term_label"]], ")"), s, perl = TRUE)
    }
  }
  s
}

# Where a dispersion's log stops being a number. A latent per row stands in for
# the row's mean, so a dispersion the data do not ask for runs off, and past
# about 25 the differences of lgamma() in beta_binomial() and Beta() are
# rounding: logLik +6.5e17 and +1.9e30. A box at 15 on a constant dispersion
# gives the logLik the uncapped fit gave where that one stayed sane, within 1e-3
# but on seeds where the uncapped fit was already rounding, and fewer fits stop
# with a nonzero code than NaN on the tape or no cap: 2 of 30 on beta_binomial
# where those give 13 and 9 (dev/latres-r1.R). A dispersion with a predictor has
# no box to hold it, so the tape reads NaN above 25.
latres_disp_box <- 15
latres_disp_cap <- 25

#' The dispersions of each response with a latent residual: its dpars on a log
#' link, other than `mu`.
#'
#' @noRd
latres_disp_dpars <- function(frame, spec) {
  out <- list()
  if (!is.finite(latres_disp_cap)) return(out)
  for (r in names(frame_latres(frame))) {
    d <- latres_log_dpars(spec$responses[[r]]$family)
    if (length(d)) out[[r]] <- d
  }
  out
}

#' A family's dpars on a log link, other than `mu`: its dispersions.
#'
#' @noRd
latres_log_dpars <- function(fam) {
  links <- fam[["links"]]
  d <- setdiff(names(links), "mu")
  d[vapply(links[d], function(l) identical(l[["name"]], "log"), NA)]
}

#' A dispersion unchanged below the cap, and NaN above it, on the tape too.
#'
#' @noRd
latres_disp_guard <- function(v) v + nan_unless_pos(exp(latres_disp_cap) - v)

#' The box on each constant dispersion of a response with a latent residual,
#' over the outer parameter vector `nm`. A tighter box the user wrote stays.
#'
#' @noRd
latres_disp_upper <- function(frame, spec, upper, map) {
  if (!is.finite(latres_disp_box)) return(upper)
  disp <- latres_disp_dpars(frame, spec)
  fixed <- frame[["betad_fixed_idx"]]
  pos <- which(map$comp == "betad")
  free <- setdiff(seq_along(frame[["par_template"]][["betad"]]), fixed)
  for (r in names(disp)) {
    for (d in disp[[r]]) {
      lp <- frame[["linpreds"]][[linpred_key(r, d)]]
      if (is.null(lp) || !identical(lp[["par"]], "betad") ||
            !isTRUE(ncol(lp[["X"]]) == 1L)) next
      k <- pos[match(lp[["idx"]], free)]
      k <- k[!is.na(k)]
      upper[k] <- pmin(upper[k], latres_disp_box)
    }
  }
  upper
}

# The weight response_mean_gh() may drop at a cell: the total there is at most
# this share of a bounded mean off, and an unbounded one has no failing node
# before its tail overflows, past any weight the rule holds.
latres_gh_lost <- 1e-6

# ---- the margin of a link that does not hold the mean

#' The smallest distance of the fit's `eta` from the link's edge, in marginal
#' sds of the latent, or `NULL` for a link that holds the mean.
#'
#' @noRd
latres_margin <- function(fit, bk) {
  a <- bk[["aux"]]
  r <- a[["ac"]][["resp"]]
  range <- latres_mu_range(fit$spec$responses[[r]]$family, FALSE)
  if (is.null(range)) return(NULL)
  th <- fit$estimates[["theta"]][bk[["theta_idx"]]]
  eta <- as.numeric(frm_linpred(fit, type = "link", dpar = "mu", resp = r))
  d <- if (identical(range, "unit")) pmin(eta, 1 - eta) else eta
  min(d, na.rm = TRUE) / latres_marginal_sd(th, a)
}

# Below this, by each fit's own margin, 14 of 150 fits read true, a start reads
# NaN on 0.42 of them and log sderr covers at 0.08. From 2 up, 0.56 read true
# with no NaN, as on links that hold the mean (dev/review-latres-m4bm.R, 680
# fits).
latres_margin_floor <- 2

#' The refusal of a fit whose mean stands too near its link's edge, after the
#' fit because only the fit knows the margin. `frm_sample()` does not reach it.
#'
#' @noRd
latres_refuse_margin <- function(fit) {
  for (bk in frame_latres(fit$frame)) {
    m <- tryCatch(latres_margin(fit, bk), error = keep_missing(NULL))
    if (!isTRUE(m < latres_margin_floor)) next
    a <- bk[["aux"]]
    fam <- fit$spec$responses[[a[["ac"]][["resp"]]]][["family"]]
    frm_stop(
      a[["ac"]][["label"]], " on ", fam[["family"]], "(link = '",
      fam[["links"]][["mu"]][["name"]] %||% "identity", "'): the mean stands ",
      format(m, digits = 2), " marginal sds of the latent from the edge of ",
      "its link, below ", latres_margin_floor, ", so a latent can carry it ",
      "out of the family's support. There, against the exact marginal, 14 of ",
      "150 fits read true and log sderr covers at 0.08. frm_sample() samples ",
      "this model, and a link that holds the mean is what frm() reads",
      call. = FALSE)
  }
  invisible(NULL)
}

#' A fit whose dispersion reached its box says so, naming the term.
#'
#' @noRd
latres_disp_report <- function(fit) {
  disp <- latres_disp_dpars(fit$frame, fit$spec)
  for (r in names(disp)) {
    dp <- eval_dpars(fit)[[r]]
    lab <- frame_latres(fit$frame)[[r]][["aux"]][["ac"]][["label"]]
    for (d in disp[[r]]) {
      big <- suppressWarnings(max(log(dp[[d]]), na.rm = TRUE))
      if (!is.finite(big) || big < latres_disp_box - 0.01) next
      frm_warning(
        "The dispersion `", d, "` beside ", lab, " ran to log(", d, ") = ",
        format(big, digits = 3), ": the latent leaves the rows nothing to ask ",
        "of it, and from a log of about 25 the objective is rounding, so it is ",
        "held. Its value is a limit and not an estimate", call. = FALSE,
        class = "frmtmb_latres_dispersion")
    }
  }
  invisible(NULL)
}
