# Review, cross-lane cases no lane tested (rellib-r14):
#   x1  mo() from a minimum of 2 and of -2 with ar(cov = TRUE) under
#       frm_sample(): simplex size, names, chart order, log_lik() against
#       the dense conditional with brms's mo() mean, loo(); a user prior
#       on sigma (innovation reading) and on simo together
#   x2  mo() new-data refusals through log_lik(newdata =),
#       posterior_predict(newdata =) and predict() on autocor models;
#       valid new data through log_lik(newdata =) against the dense
#       conditional of the new groups
#   x3  sar() prior tables, lag and error, both routes, on a weights
#       matrix whose smallest eigenvalue is above -1; validate_prior();
#       the tables passed back to frm() and frm_sample(); brms 2.23.0's
#       default_prior() and stancode() bounds on the same data
#   x4  the seam's jacobian slot with acprior's class coefficients: a
#       bounds-only prior on errorsar under frm_sample(), then log_lik()
#       against the dense conditional of the error SAR
#   x5  brms-40 and brms-41's stated workarounds in frmtmb
#   Rscript dev/review-rel076-cross.R <arm> > dev/review-rel076-out/cross-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "\n")
cond <- function(y, m, S) {
  vapply(seq_along(y), function(i) {
    o <- seq_along(y)[-i]
    A <- solve(S[o, o, drop = FALSE])
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% (y[o] - m[o]))
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    stats::dnorm(y[i], cm, sqrt(cv), log = TRUE)
  }, 0)
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
try_msg <- function(lab, expr) {
  r <- tryCatch({
    v <- suppressWarnings(expr)
    paste("RAN, dim", paste(dim(v) %||% length(v), collapse = "x"))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-34s %s\n", lab, substr(gsub("\n", " ", r), 1, 230)))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
sim_ar <- function(ng, nt, phi, sd) {
  unlist(lapply(seq_len(ng), function(g) {
    as.numeric(stats::arima.sim(list(ar = phi), nt, sd = sd))
  }))
}
mo_ar_data <- function(lo) {
  set.seed(20261029)
  d <- expand.grid(t = 1:5, g = 1:24)
  d$x <- stats::rnorm(nrow(d))
  d$k <- sample(lo:(lo + 3), nrow(d), replace = TRUE)
  d$y <- 1 + 0.5 * d$x + c(0, 0.5, 0.8, 1.4)[d$k - lo + 1] +
    sim_ar(24, 5, 0.6, 0.6)
  d
}
check_mo_ar <- function(ds, d) {
  m <- as.matrix(ds$draws)
  cat("  variables:", paste(colnames(m), collapse = " "), "\n")
  sx <- grep("^simo_mok1\\[", colnames(m))
  D <- length(sx)
  cat("  simplex weights:", D, "for values", min(d$k), "to", max(d$k),
      "| rows sum to 1:", isTRUE(all.equal(rowSums(m[, sx]),
                                            rep(1, nrow(m)))),
      "| ar in (-1, 1):", all(abs(m[, "ar[1]"]) < 1), "\n")
  ll <- log_lik(ds)
  worst <- 0
  for (s in seq_len(nrow(m))) {
    cs <- c(0, cumsum(m[s, sx]))
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x +
      m[s, "bsp_mok"] * D * cs[d$k - min(d$k) + 1]
    ref <- numeric(nrow(d))
    for (g in unique(d$g)) {
      r <- which(d$g == g)
      S <- m[s, "sigma"]^2 * m[s, "ar[1]"]^abs(outer(d$t[r], d$t[r], "-"))
      ref[r] <- cond(d$y[r], mu[r], S)
    }
    worst <- max(worst, rel(ll[s, ], ref))
  }
  cat("  log_lik", paste(dim(ll), collapse = "x"), "against the dense",
      "conditional (brms's mo() mean): max relative", format(worst), "\n")
  lo <- suppressWarnings(loo(ds))
  cat("  loo elpd", format(lo$estimates["elpd_loo", "Estimate"]),
      "pareto k > 0.7:", sum(lo$diagnostics$pareto_k > 0.7), "\n")
  invisible(m)
}

if (arm == "x1" || arm == "x2") {
  fo <- y ~ x + mo(k) + ar(t, g, cov = TRUE)
  for (lo in c(2L, -2L)) {
    d <- mo_ar_data(lo)
    cat("\n== x1: mo(k) with k in", lo, "to", lo + 3, "+ ar(cov = TRUE)\n")
    f <- frm(fo, data = d)
    cat("  frm(): converged", isTRUE(f$opt$convergence == 0), " logLik",
        format(as.numeric(logLik(f))), "\n")
    dp <- default_prior(fo, data = d, route = "sample")
    print(dp[, c("prior", "class", "coef", "lb", "ub")], row.names = FALSE)
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fo, data = d, chains = 2, iter = 600, seed = 4,
                 refresh = 0)))
    m <- check_mo_ar(ds, d)
    if (lo == 2L) {
      pr <- c(set_prior("normal(0, 0.5)", class = "sigma"),
              set_prior("dirichlet(2)", class = "simo", coef = "mok1"))
      ds2 <- suppressWarnings(suppressMessages(
        frm_sample(fo, data = d, prior = pr, chains = 2, iter = 600,
                   seed = 4, refresh = 0)))
      cat("  -- user priors on sigma and simo:\n")
      check_mo_ar(ds2, d)
      fm <- frm(fo, data = d, prior = pr)
      cat("  frm(prior =) MAP converged:", isTRUE(fm$opt$convergence == 0),
          "\n")
    }
    if (arm == "x2") {
      cat("\n== x2: new data refusals, k fitted on", lo, "to", lo + 3, "\n")
      nd0 <- d[d$g %in% 1:3, ]
      bad <- list(
        fractional = within(nd0, k[1] <- lo + 0.5),
        missing = within(nd0, k[1] <- NA),
        factor = within(nd0, k <- factor(k)),
        character = within(nd0, k <- as.character(k)),
        `below minimum` = within(nd0, k[1] <- lo - 1L),
        `above maximum` = within(nd0, k[1] <- lo + 4L))
      for (nm in names(bad)) {
        cat(" ", nm, ":\n")
        try_msg("log_lik(draws, newdata)", log_lik(ds, newdata = bad[[nm]]))
        try_msg("posterior_predict(draws, newdata)",
                posterior_predict(ds, newdata = bad[[nm]]))
        try_msg("posterior_epred(draws, newdata)",
                posterior_epred(ds, newdata = bad[[nm]]))
        try_msg("predict(frm fit, newdata)", predict(f, newdata = bad[[nm]]))
      }
      cat("  valid new data (groups 1 to 3):\n")
      try_msg("posterior_predict(draws, newdata)",
              posterior_predict(ds, newdata = nd0))
      lln <- tryCatch(log_lik(ds, newdata = nd0), error = function(e) e)
      if (inherits(lln, "error")) {
        cat("  log_lik(newdata) ERROR:", conditionMessage(lln), "\n")
      } else {
        sx <- grep("^simo_mok1\\[", colnames(m))
        D <- length(sx); worst <- 0
        for (s in seq_len(nrow(m))) {
          cs <- c(0, cumsum(m[s, sx]))
          mu <- m[s, "b_Intercept"] + m[s, "b_x"] * nd0$x +
            m[s, "bsp_mok"] * D * cs[nd0$k - min(d$k) + 1]
          ref <- numeric(nrow(nd0))
          for (g in unique(nd0$g)) {
            r <- which(nd0$g == g)
            S <- m[s, "sigma"]^2 *
              m[s, "ar[1]"]^abs(outer(nd0$t[r], nd0$t[r], "-"))
            ref[r] <- cond(nd0$y[r], mu[r], S)
          }
          worst <- max(worst, rel(lln[s, ], ref))
        }
        cat("  log_lik(newdata)", paste(dim(lln), collapse = "x"),
            "against the dense conditional of the new groups: max",
            "relative", format(worst), "\n")
        llf <- log_lik(ds)
        cat("  equals the in-sample columns of those rows:",
            isTRUE(all.equal(unname(lln), unname(llf[, d$g %in% 1:3]),
                             tolerance = 1e-12)), "\n")
      }
    }
  }
}

sar_w <- function() {
  # a ring of 20 with chords, row-standardized: odd cycles, so the
  # smallest eigenvalue is above -1
  n <- 20
  A <- matrix(0, n, n)
  for (i in 1:n) {
    j <- i %% n + 1; A[i, j] <- A[j, i] <- 1
    j <- (i + 1) %% n + 1; A[i, j] <- A[j, i] <- 1
  }
  A[1, 10] <- A[10, 1] <- 1
  A / rowSums(A)
}

if (arm == "x3" || arm == "x4") {
  suppressMessages(library(frmtmb.spatial))
  cat("frmtmb.spatial", format(packageVersion("frmtmb.spatial")), "\n")
  W <- sar_w()
  ev <- Re(eigen(W, only.values = TRUE)$values)
  cat("W: n", nrow(W), "eigenvalues from", format(min(ev)), "to",
      format(max(ev)), "\n")
  n <- nrow(W)
  set.seed(3)
  d <- data.frame(x = stats::rnorm(n * 1))
  d$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * d$x +
                           stats::rnorm(n, 0, 0.7)))
  d2 <- list(W = W)
}

if (arm == "x3") {
  for (ty in c("lag", "error")) {
    fo <- stats::as.formula(paste0("y ~ x + sar(W, type = '", ty, "')"))
    cl <- if (ty == "lag") "lagsar" else "errorsar"
    cat("\n== x3:", ty, "SAR\n")
    for (route in c("fit", "sample")) {
      dp <- frmtmb::default_prior(fo, data = d, data2 = d2, route = route)
      cat("  route", route, ": rows", paste(dp$class, collapse = " "),
          "|", cl, "prior", dp$prior[dp$class == cl], "lb",
          dp$lb[dp$class == cl], "ub", dp$ub[dp$class == cl], "\n")
      vp <- tryCatch(validate_prior(dp, fo, data = d, data2 = d2),
                     error = function(e) e)
      cat("  validate_prior(table):", if (inherits(vp, "error"))
        paste("ERROR", conditionMessage(vp)) else
          paste(nrow(vp), "rows; source", paste(unique(vp$source),
                                                 collapse = ",")), "\n")
      if (route == "fit") {
        f0 <- frm(fo, data = d, data2 = d2)
        f1 <- tryCatch(frm(fo, data = d, data2 = d2, prior = dp),
                       error = function(e) e)
        cat("  frm(prior = fit table):", if (inherits(f1, "error"))
          paste("ERROR", conditionMessage(f1)) else
            paste("logLik equal", isTRUE(all.equal(
              as.numeric(logLik(f0)), as.numeric(logLik(f1)),
              tolerance = 1e-10)), "| estimate", format(coef(f1)[[1]] %||%
                                                          NA)), "\n")
        cat("  frm estimate of", cl, ":",
            format(summary(f0)$autocor %||% NA), "\n")
      } else {
        ds <- tryCatch(suppressWarnings(suppressMessages(
          frm_sample(fo, data = d, data2 = d2, prior = dp, chains = 2,
                     iter = 500, seed = 3, refresh = 0))),
          error = function(e) e)
        if (inherits(ds, "error")) {
          cat("  frm_sample(prior = sample table) ERROR:",
              conditionMessage(ds), "\n")
        } else {
          m <- as.matrix(ds$draws)
          cat("  frm_sample(prior = sample table):", cl, "draws from",
              format(min(m[, cl])), "to", format(max(m[, cl])), "\n")
        }
      }
    }
    # brms 2.23.0 on the same data, no compile
    bp <- tryCatch({
      stopifnot(requireNamespace("brms", quietly = TRUE))
      fb <- stats::as.formula(paste0("y ~ x + sar(W, type = '", ty, "')"))
      g <- brms::default_prior(fb, data = d, data2 = d2)
      sc <- brms::stancode(fb, data = d, data2 = d2)
      list(g = g, sc = grep(cl, strsplit(sc, "\n")[[1]], value = TRUE))
    }, error = function(e) e)
    if (inherits(bp, "error")) {
      cat("  brms ERROR:", conditionMessage(bp), "\n")
    } else {
      cat("  brms default_prior", cl, "row: prior", bp$g$prior[bp$g$class == cl],
          "lb", bp$g$lb[bp$g$class == cl], "ub", bp$g$ub[bp$g$class == cl],
          "\n")
      cat("  brms stancode lines:\n"); cat(paste("    ", bp$sc), sep = "\n")
      sd <- brms::standata(stats::as.formula(paste0(
        "y ~ x + sar(W, type = '", ty, "')")), data = d, data2 = d2)
      cat("  brms standata min_eigenMsar", format(sd$min_eigenMsar, digits = 17),
          "max_eigenMsar", format(sd$max_eigenMsar, digits = 17), "\n")
    }
    # the bounds a user types for a row-standardized W: the ends -1 and 1
    fo2 <- fo
    dpf <- frmtmb::default_prior(fo, data = d, data2 = d2)
    cat("  frmtmb's table ends, 17 digits: lb",
        format(dpf$lb[dpf$class == cl], digits = 17), "ub",
        format(dpf$ub[dpf$class == cl], digits = 17), "\n")
    for (b in list(c(NA, 1), c(0, 1), c(0, NA))) {
      pr <- set_prior("", class = cl, lb = b[1], ub = b[2])
      r <- tryCatch({frm(fo2, data = d, data2 = d2, prior = pr); "RAN"},
                    error = function(e) paste("ERROR", conditionMessage(e)))
      cat("  frm(", cl, "lb =", b[1], "ub =", b[2], "):", r, "\n")
    }
  }
}

if (arm == "x4") {
  fo <- y ~ x + sar(W, type = "error")
  f <- frm(fo, data = d, data2 = d2)
  fr <- f$frame
  ac <- fr[["autocor"]][[1]]
  cc <- tryCatch(frmtmb:::autocor_class_coefs(fr), error = function(e) e)
  cat("\n== x4: autocor_class_coefs(frame) for errorsar:\n")
  utils::str(cc, max.level = 2)
  cat("  seam entry has a jacobian:", !is.null(ac[["jacobian"]]) ||
        !is.null(frmtmb:::cs_entry), "\n")
  pr <- set_prior("", class = "errorsar", lb = 0)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fo, data = d, data2 = d2, prior = pr, chains = 2,
               iter = 600, seed = 3, refresh = 0)))
  m <- as.matrix(ds$draws)
  cat("  variables:", paste(colnames(m), collapse = " "), "\n")
  cat("  errorsar draws from", format(min(m[, "errorsar"])), "to",
      format(max(m[, "errorsar"])), "(lb = 0)\n")
  n <- nrow(d)
  ll <- log_lik(ds)
  worst <- 0
  for (s in seq_len(nrow(m))) {
    IA <- diag(n) - m[s, "errorsar"] * W
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x
    S <- m[s, "sigma"]^2 * solve(crossprod(IA))
    worst <- max(worst, rel(ll[s, ], cond(d$y, mu, S)))
  }
  cat("  log_lik against the dense conditional of the error SAR: max",
      "relative", format(worst), "\n")
  lo <- suppressWarnings(loo(ds))
  cat("  loo elpd", format(lo$estimates["elpd_loo", "Estimate"]), "\n")
  # a bound at the end of the interval passes back (bounds nothing)
  pr2 <- set_prior("", class = "errorsar", lb = -1, ub = 1)
  f2 <- tryCatch(frm(fo, data = d, data2 = d2, prior = pr2),
                 error = function(e) e)
  cat("  frm(errorsar lb = -1, ub = 1):", if (inherits(f2, "error"))
    paste("ERROR", conditionMessage(f2)) else
      paste("logLik equal", isTRUE(all.equal(as.numeric(logLik(f)),
                                             as.numeric(logLik(f2))))),
    "\n")
  pr3 <- set_prior("", class = "errorsar", lb = -1.5)
  f3 <- tryCatch(frm(fo, data = d, data2 = d2, prior = pr3),
                 error = function(e) e)
  cat("  frm(errorsar lb = -1.5):", if (inherits(f3, "error"))
    paste("ERROR", conditionMessage(f3)) else "RAN", "\n")
  pr4 <- set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1")
  f4 <- tryCatch(frm(fo, data = d, data2 = d2, prior = pr4),
                 error = function(e) e)
  cat("  frm(class theta thetaac_1):", if (inherits(f4, "error"))
    paste("ERROR", conditionMessage(f4)) else "RAN", "\n")
}

if (arm == "x5") {
  cat("\n== x5: brms-40, attr(x, 'min') is not read\n")
  set.seed(5)
  d <- data.frame(x = c(1, 2, 3, 4)[sample(4, 200, TRUE)])
  d$y <- 0.4 * d$x + stats::rnorm(200)
  d1 <- d; attr(d1$x, "min") <- 0
  f0 <- frm(y ~ mo(x), data = d)
  f1 <- frm(y ~ mo(x), data = d1)
  cat("  logLik equal with and without the attribute:",
      isTRUE(all.equal(as.numeric(logLik(f0)), as.numeric(logLik(f1)))),
      "| simplex length", length(f1$estimates[["zeta1"]] %||% NA), "\n")
  nd <- data.frame(x = 1:4)
  cat("  fitted(newdata = 1:4) equal:",
      isTRUE(all.equal(fitted(f0, newdata = nd), fitted(f1, newdata = nd))),
      "\n")
  cat("\n== x5: brms-41, mv subset() codes new data per response\n")
  set.seed(5)
  n <- 300
  d <- data.frame(k = sample(1:4, n, TRUE))
  d$y <- 0.3 * d$k + stats::rnorm(n)
  d$y2 <- -0.2 * d$k + stats::rnorm(n)
  d$s1 <- d$k > 1
  d$s2 <- TRUE
  fo <- bf(y | subset(s1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)) +
    set_rescor(FALSE)
  f <- frm(fo, data = d)
  nd <- data.frame(k = 2:4, s1 = TRUE, s2 = TRUE)
  for (r in c("y", "y2")) {
    p <- fitted(f, newdata = nd, resp = r)[, 1]
    # the fit's own rows with those values: the same model mean
    fr <- fitted(f, resp = r)[, 1]
    kk <- d$k[if (r == "y") d$s1 else d$s2]
    own <- tapply(fr, kk, mean)[as.character(2:4)]
    cat("  resp", r, "fitted(newdata k = 2:4):", signif(p, 6),
        "| fit's rows at 2:4:", signif(own, 6), "| equal:",
        isTRUE(all.equal(unname(p), unname(own))), "\n")
  }
  nd1 <- data.frame(k = 1, s1 = TRUE, s2 = TRUE)
  p1 <- tryCatch(fitted(f, newdata = nd1, resp = "y"), error = function(e) e)
  cat("  predict(resp y, k = 1, below y's fitted minimum 2):",
      if (inherits(p1, "error")) paste("ERROR", conditionMessage(p1)) else
        "RAN (not refused)", "\n")
}
if (arm == "x6") {
  # the chart (acprior) and the mo() softmax (mocode) in one draw matrix:
  # a bounds-only prior on ar puts ar through autocor_bound_chart(), and
  # chart_draws_back() runs before mo_draws_to_chart()
  d <- mo_ar_data(2L)
  fo <- y ~ x + mo(k) + ar(t, g, cov = TRUE)
  cat("\n== x6: mo(k) + ar(cov = TRUE) with set_prior('', class = 'ar',",
      "lb = 0.2): the chart and the simplex together\n")
  pr <- set_prior("", class = "ar", lb = 0.2)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fo, data = d, prior = pr, chains = 2, iter = 600, seed = 4,
               refresh = 0)))
  m <- check_mo_ar(ds, d)
  cat("  ar draws from", format(min(m[, "ar[1]"])), "to",
      format(max(m[, "ar[1]"])), "(lb = 0.2)\n")
  # thetaac_1 is the internal parameter of ar: ar = th / sqrt(1 + th^2)
  th <- m[, "thetaac_1"]
  cat("  thetaac_1 column is ar's internal parameter:",
      isTRUE(all.equal(unname(th / sqrt(1 + th^2)), unname(m[, "ar[1]"]))),
      "\n")
}
cat("CROSS DONE", arm, "\n")
