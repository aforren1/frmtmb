# REVIEW script 03, claim 1: the row density identity and an
# INDEPENDENT reference for the cells.
#
# Own design and own seeds (5701). Unequal group lengths, interior time
# gaps, rows shuffled, p = 2 and q = 2 where the term allows it.
#
# Three instruments, and the review says which is which:
#
#  (1) IDENTITY, tape form. For a model with NO random effect
#      `fit$obj$fn(p)` IS the taped negative log likelihood, so the row
#      sum must equal `-obj$fn(p)`. The internal parameter vector is
#      rebuilt from the draw by inverting `parList()`, and the inversion
#      is checked by a round trip on `last.par.best` before it is used.
#  (2) IDENTITY, R form. `build_objective(frame)(estimates)` is the R
#      composition RTMB tapes, and for a random-effect model it is the
#      JOINT nll, so the row sum must equal it minus each block's own
#      prior. Used where obj$fn is a Laplace marginal instead.
#  (3) MEASUREMENT, and the one that is not frmtmb's own code:
#      brms 2.23.0's generated model block for `arma()`, transliterated
#      here from `make_stancode()` output read in this review, driven by
#      the draws matrix's OWN brms-scale columns. The aterm cases add
#      their own closed forms (weights multiply, cens swaps the density
#      for a tail probability, trunc divides by the retained mass).
#
#   Rscript dev/arcovsample-rev-03-rows.R

LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))
cat("frmtmb ", format(packageVersion("frmtmb")), " frmtmb.sample ",
    format(packageVersion("frmtmb.sample")), " (lane build)\n", sep = "")

SEED <- 5701L

# ---- brms 2.23.0's model block, transliterated -----------------------
# for (n in 1:N) {
#   mu[n] += Err[n, 1:Kma] * ma;
#   err[n] = Y[n] - mu[n];
#   for (i in 1:J_lag[n]) Err[n + 1, i] = err[n + 1 - i];
#   mu[n] += Err[n, 1:Kar] * ar;
# }
# Written independently of the package's own test helper, from
# make_stancode() output printed in this review; `Err` is N x max_lag
# there and the extra row here is never read.
rev_mu <- function(mu, Y, J_lag, ar, ma) {
  N <- length(Y)
  ml <- max(length(ar), length(ma))
  Err <- matrix(0, N + 1L, ml)
  err <- numeric(N)
  for (n in seq_len(N)) {
    if (length(ma)) mu[n] <- mu[n] + sum(Err[n, seq_along(ma)] * ma)
    err[n] <- Y[n] - mu[n]
    if (J_lag[n] > 0L) {
      for (i in seq_len(J_lag[n])) Err[n + 1L, i] <- err[n + 1L - i]
    }
    if (length(ar)) mu[n] <- mu[n] + sum(Err[n, seq_along(ar)] * ar)
  }
  mu
}
# brms's data_ac: J_lag[n] counts how many of the previous max_lag rows
# belong to the same group as row n + 1, in the (gr, time) sort order.
rev_jlag <- function(g, max_lag) {
  N <- length(g)
  J <- integer(N)
  for (n in seq_len(N - 1L)) {
    ind <- n:max(1L, n + 1L - max_lag)
    J[n] <- sum(g[ind] == g[n + 1L])
  }
  J
}

# ---- data: unequal group lengths, interior gaps, shuffled ------------
mk_data <- function(seed = SEED) {
  set.seed(seed)
  ng <- 6L
  # unequal lengths BY CONSTRUCTION, 4 to 11 rows, then interior time
  # points removed so that a lag in rows is not a lag in time
  lens <- c(4L, 7L, 11L, 5L, 9L, 6L)
  rows <- do.call(rbind, lapply(seq_len(ng), function(i) {
    tt <- seq_len(lens[i] + 2L)
    if (lens[i] >= 5L) tt <- tt[-c(3L, 5L)] else tt <- tt[-2L]
    data.frame(g = factor(i, levels = seq_len(ng)), t = tt[seq_len(lens[i])])
  }))
  n <- nrow(rows)
  rows$x <- stats::rnorm(n)
  rows$z <- stats::runif(n, -1, 1)
  u <- stats::rnorm(ng, 0, 0.5)
  e <- unlist(lapply(split(seq_len(n), rows$g), function(r) {
    as.numeric(stats::arima.sim(list(ar = c(0.4, 0.2), ma = c(0.3, 0.15)),
                                length(r), sd = 0.5))
  }))
  rows$y <- 0.8 + 0.6 * rows$x + u[as.integer(rows$g)] + e
  rows$yp <- rows$y + 4
  rows$w <- stats::runif(n, 0.4, 2.2)
  rows$cc <- sample(c(0L, 0L, 0L, 1L, -1L), n, TRUE)
  rows$y2 <- -0.3 + 0.4 * rows$x + stats::rnorm(n, 0, 0.7)
  rows[sample(n), ]
}
dd <- mk_data()
cat("N = ", nrow(dd), "  group sizes: ",
    paste(as.integer(table(dd$g)), collapse = " "), "  seed ", SEED,
    "\n", sep = "")

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- the inverse of parList(), checked by a round trip ---------------
par_vec <- function(fit, est) {
  p0 <- fit$obj$env$last.par.best
  # round trip first: if parList() and this inversion disagree on the
  # fit's OWN vector the instrument is broken and nothing below counts
  rt <- p0
  l0 <- fit$obj$env$parList(p0)
  for (nm in unique(names(p0))) rt[names(p0) == nm] <- as.numeric(l0[[nm]])
  stopifnot(identical(unname(rt), unname(p0)))
  p <- p0
  for (nm in unique(names(p0))) p[names(p0) == nm] <- as.numeric(est[[nm]])
  p
}

re_prior_of <- function(sh) {
  tot <- 0
  for (bk in sh$frame[["re_blocks"]] %||% list()) {
    f <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll
    tot <- tot + as.numeric(f(sh$estimates[["b"]][bk$b_idx],
                              sh$estimates[["theta"]][bk$theta_idx], bk))
  }
  tot
}

# ---- one model: sample, then run all three instruments ---------------
run_one <- function(label, form, family, data = dd, iter = 400L,
                    cells = NULL, has_re = FALSE) {
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, family = family, data = data, chains = 1,
               iter = iter, refresh = 0, seed = 11)))
  fit <- frmtmb.sample:::draws_base_fit(ds)
  ll <- log_lik(ds)
  nd <- nrow(ds$draws)
  idx <- frmtmb.sample:::draws_par_index(fit)
  r_obj <- rep(NA_real_, nd)
  r_bld <- numeric(nd)
  for (k in seq_len(nd)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k, idx)
    nllR <- as.numeric(frmtmb::build_objective(sh$frame)(sh$estimates))
    r_bld[k] <- sum(ll[k, ]) - (-nllR - re_prior_of(sh))
    if (!has_re) {
      p <- par_vec(fit, sh$estimates)
      r_obj[k] <- sum(ll[k, ]) - (-as.numeric(fit$obj$fn(p)))
    }
  }
  cat("\n---- ", label, " ----\n", sep = "")
  cat("  draws x cols: ", nd, " x ", ncol(ll), "; any non-finite cell: ",
      any(!is.finite(ll)), "\n", sep = "")
  cat("  (2) max |rowsum - (-build_objective + re_prior)| : ",
      format(max(abs(r_bld)), digits = 12), "  rel ",
      format(max(abs(r_bld) / abs(rowSums(ll))), digits = 4), "\n",
      sep = "")
  if (!has_re) {
    cat("  (1) max |rowsum - (-obj$fn)| TAPE                 : ",
        format(max(abs(r_obj)), digits = 12), "  rel ",
        format(max(abs(r_obj) / abs(rowSums(ll))), digits = 4), "\n",
        sep = "")
  } else {
    cat("  (1) obj$fn is a LAPLACE marginal here, not compared\n")
  }
  if (!is.null(cells)) {
    dm <- ds$draws
    dev <- vapply(seq_len(nd), function(k) {
      max(abs(ll[k, ] - cells(dm[k, ], data)))
    }, numeric(1))
    cat("  (3) max |cell - brms transliteration| over draws  : ",
        format(max(dev), digits = 6), "  cell sd ",
        format(stats::sd(ll), digits = 4), "\n", sep = "")
    # the shift must MATTER: the same reference with ar = ma = 0
    cat("      spread of rowSums across draws: ",
        format(stats::sd(rowSums(ll)), digits = 6), "\n", sep = "")
  }
  lo <- suppressWarnings(loo(ds))
  cat("  elpd_loo ", format(lo$estimates["elpd_loo", "Estimate"],
                            digits = 9),
      "  finite pointwise: ", all(is.finite(lo$pointwise[, "elpd_loo"])),
      "\n", sep = "")
  invisible(list(ds = ds, ll = ll))
}

ord <- order(dd$g, dd$t)
gs <- as.integer(dd$g)[ord]

# mu0 from the draws matrix's own brms-scale columns, plus a check that
# this really is frmtmb's unshifted linear predictor
lin0 <- function(v, data, re = FALSE) {
  m <- v[["b_Intercept"]] + v[["b_x"]] * data$x
  if (re) {
    rn <- grep("^r_g\\[", names(v), value = TRUE)
    lev <- sub("^r_g\\[([^,]+),.*$", "\\1", rn)
    m <- m + unname(v[rn])[match(as.character(data$g), lev)]
  }
  m
}
cells_gauss <- function(p, q, re = FALSE, sig = NULL) {
  function(v, data) {
    ar <- if (p > 0L) unname(v[paste0("thetaac_", seq_len(p))]) else numeric(0)
    ma <- if (q > 0L) {
      unname(v[paste0("thetaac_", p + seq_len(q))])
    } else numeric(0)
    J <- rev_jlag(gs, max(p, q))
    mu <- rev_mu(lin0(v, data, re)[ord], data$y[ord], J, ar, ma)
    s <- if (is.null(sig)) rep(v[["sigma"]], nrow(data)) else sig(v, data)[ord]
    out <- numeric(nrow(data))
    out[ord] <- stats::dnorm(data$y[ord], mu, s, log = TRUE)
    out
  }
}

res <- list()

# 1. gaussian arma(2,2), no random effect: all three instruments
res$arma22 <- run_one("gaussian arma(2,2), p = q = 2, no RE",
                      bf(y ~ x + arma(t, g, p = 2, q = 2)), gaussian(),
                      cells = cells_gauss(2L, 2L))

# 2. distributional sigma
res$sigma <- run_one("gaussian ar(2), sigma ~ z",
                     bf(y ~ x + ar(t, g, p = 2), sigma ~ z), gaussian(),
                     cells = function(v, data) {
                       ar <- unname(v[c("thetaac_1", "thetaac_2")])
                       J <- rev_jlag(gs, 2L)
                       mu <- rev_mu(lin0(v, data)[ord], data$y[ord], J, ar,
                                    numeric(0))
                       s <- exp(v[["b_sigma_Intercept"]] +
                                  v[["b_sigma_z"]] * data$z)[ord]
                       out <- numeric(nrow(data))
                       out[ord] <- stats::dnorm(data$y[ord], mu, s,
                                                log = TRUE)
                       out
                     })

# 3. weights(): the column is w * the density
res$weights <- run_one("gaussian ma(2), weights(w)",
                       bf(y | weights(w) ~ x + ma(t, g, q = 2)),
                       gaussian(),
                       cells = function(v, data) {
                         ma <- unname(v[c("thetaac_1", "thetaac_2")])
                         J <- rev_jlag(gs, 2L)
                         mu <- rev_mu(lin0(v, data)[ord], data$y[ord], J,
                                      numeric(0), ma)
                         out <- numeric(nrow(data))
                         out[ord] <- stats::dnorm(data$y[ord], mu,
                                                  v[["sigma"]],
                                                  log = TRUE)
                         out * data$w
                       })

# 4. cens(): +1 right censored -> log S, -1 left -> log F
res$cens <- run_one("gaussian arma(2,2), cens(cc)",
                    bf(y | cens(cc) ~ x + arma(t, g, p = 2, q = 2)),
                    gaussian(),
                    cells = function(v, data) {
                      ar <- unname(v[c("thetaac_1", "thetaac_2")])
                      ma <- unname(v[c("thetaac_3", "thetaac_4")])
                      J <- rev_jlag(gs, 2L)
                      mu <- rev_mu(lin0(v, data)[ord], data$y[ord], J, ar,
                                   ma)
                      s <- v[["sigma"]]
                      yy <- data$y[ord]
                      cc <- data$cc[ord]
                      o <- stats::dnorm(yy, mu, s, log = TRUE)
                      o[cc == 1L] <- stats::pnorm(yy[cc == 1L],
                                                  mu[cc == 1L], s,
                                                  lower.tail = FALSE,
                                                  log.p = TRUE)
                      o[cc == -1L] <- stats::pnorm(yy[cc == -1L],
                                                   mu[cc == -1L], s,
                                                   log.p = TRUE)
                      out <- numeric(nrow(data))
                      out[ord] <- o
                      out
                    })

# 5. trunc(): the density divided by the retained mass
res$trunc <- run_one("gaussian ar(2), trunc(lb = 0) on y + 4",
                     bf(yp | trunc(lb = 0) ~ x + ar(t, g, p = 2)),
                     gaussian(),
                     cells = function(v, data) {
                       ar <- unname(v[c("thetaac_1", "thetaac_2")])
                       J <- rev_jlag(gs, 2L)
                       mu0 <- v[["b_Intercept"]] + v[["b_x"]] * data$x
                       mu <- rev_mu(mu0[ord], data$yp[ord], J, ar,
                                    numeric(0))
                       s <- v[["sigma"]]
                       o <- stats::dnorm(data$yp[ord], mu, s, log = TRUE) -
                         stats::pnorm(0, mu, s, lower.tail = FALSE,
                                      log.p = TRUE)
                       out <- numeric(nrow(data))
                       out[ord] <- o
                       out
                     })

# 6. student() with nu fitted, plus (1 | g)
res$stu <- run_one("student arma(2,2) + (1 | g), nu fitted",
                   bf(y ~ x + arma(t, g, p = 2, q = 2) + (1 | g)),
                   student(), has_re = TRUE,
                   cells = function(v, data) {
                     ar <- unname(v[c("thetaac_1", "thetaac_2")])
                     ma <- unname(v[c("thetaac_3", "thetaac_4")])
                     J <- rev_jlag(gs, 2L)
                     mu <- rev_mu(lin0(v, data, re = TRUE)[ord],
                                  data$y[ord], J, ar, ma)
                     s <- v[["sigma"]]
                     out <- numeric(nrow(data))
                     out[ord] <- stats::dt((data$y[ord] - mu) / s,
                                           df = v[["nu"]], log = TRUE) -
                       log(s)
                     out
                   })

# 7. set_rescor(TRUE), ar(1) on both responses, gaussian and student.
#    The cell reference is mvtnorm, which is outside frmtmb.
run_rescor <- function(label, family) {
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + ar(t, g)) + bf(y2 ~ x + ar(t, g)) +
                 set_rescor(TRUE), family = family, data = dd,
               chains = 1, iter = 400L, refresh = 0, seed = 11)))
  fit <- frmtmb.sample:::draws_base_fit(ds)
  ll <- log_lik(ds)
  nd <- nrow(ds$draws)
  idx <- frmtmb.sample:::draws_par_index(fit)
  r_obj <- numeric(nd); r_mvt <- numeric(nd)
  J <- rev_jlag(gs, 1L)
  for (k in seq_len(nd)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k, idx)
    p <- par_vec(fit, sh$estimates)
    r_obj[k] <- sum(ll[k, ]) - (-as.numeric(fit$obj$fn(p)))
    # mvtnorm at the SAME mu, sigma and C, with mu from the review's own
    # transliteration rather than from eval_dpars()
    v <- ds$draws[k, ]
    C <- frmtmb:::us_chol_cor(sh$estimates[["thetar"]], 2L)
    m1 <- rev_mu((v[["b_y_Intercept"]] + v[["b_y_x"]] * dd$x)[ord],
                 dd$y[ord], J, v[["thetaac_1"]], numeric(0))
    m2 <- rev_mu((v[["b_y2_Intercept"]] + v[["b_y2_x"]] * dd$x)[ord],
                 dd$y2[ord], J, v[["thetaac_2"]], numeric(0))
    s1 <- v[["sigma_y"]]; s2 <- v[["sigma_y2"]]
    nu <- if ("nu" %in% names(v)) v[["nu"]] else NULL
    S <- diag(c(s1, s2)) %*% C %*% diag(c(s1, s2))
    ref <- vapply(seq_along(ord), function(i) {
      xx <- c(dd$y[ord][i], dd$y2[ord][i]); mm <- c(m1[i], m2[i])
      if (is.null(nu)) {
        mvtnorm::dmvnorm(xx, mean = mm, sigma = S, log = TRUE)
      } else {
        suppressWarnings(mvtnorm::dmvt(xx, delta = mm, sigma = S,
                                       df = nu, log = TRUE))
      }
    }, numeric(1))
    o <- numeric(nrow(dd)); o[ord] <- ref
    r_mvt[k] <- max(abs(ll[k, ] - o))
  }
  cat("\n---- ", label, " ----\n", sep = "")
  cat("  draws x cols: ", nd, " x ", ncol(ll), "; non-finite cell: ",
      any(!is.finite(ll)), "\n", sep = "")
  cat("  (1) max |rowsum - (-obj$fn)| TAPE : ",
      format(max(abs(r_obj)), digits = 12), "  rel ",
      format(max(abs(r_obj) / abs(rowSums(ll))), digits = 4), "\n",
      sep = "")
  cat("  (3) max |cell - mvtnorm| over draws: ",
      format(max(r_mvt), digits = 6), "  cell sd ",
      format(stats::sd(ll), digits = 4), "\n", sep = "")
  if ("nu" %in% colnames(ds$draws)) {
    cat("      nu range: ", format(min(ds$draws[, "nu"]), digits = 4),
        " to ", format(max(ds$draws[, "nu"]), digits = 4), "\n", sep = "")
  }
  lo <- suppressWarnings(loo(ds))
  cat("  elpd_loo ", format(lo$estimates["elpd_loo", "Estimate"],
                            digits = 9), "\n", sep = "")
  invisible(ds)
}
res$rc_g <- run_rescor("rescor + ar(1) on both, gaussian", gaussian())
res$rc_t <- run_rescor("rescor + ar(1) on both, student", student())

cat("\nDONE\n")
