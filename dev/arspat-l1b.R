# Lane arspat: the L1 identity of each lifted combination against brms
# 2.23.0's own Stan program, and its fit against brms's optimum.
#
#   Rscript dev/arspat-l1.R [case ...]      (from the worktree root)
#
# A joint case (se()) reads brms's log_prob at frmtmb's estimate and at two
# points stepped 0.05 off it in every coordinate (dev/equiv-lib.R eq_l1()):
# the constant lp - ours must be 0 at all three, and brms's gradient must
# vanish at the estimate. A missing-value case reads frmtmb's logLik, which
# integrates the missing values out, against brms's program integrated
# over Ymi: by its Laplace value on gaussian (exact, the program is
# quadratic in Ymi), and by quadrature over each missing value on student
# (the groups are independent, so the integral factors by group and each
# factor is one-dimensional). Version b of dev/arspat-l1.R: moderate nu on
# the student cases, a real (1 | h) effect, a smaller unstr() design, and
# the sar() cases. Output: dev/arspat-out/l1b-all.txt
lib <- Sys.getenv("ARSPAT_LIB", "C:/Users/adf44/source/r/wt-arspat-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R")) sys.source(f, envir = eq_env)

# tdf: the residual of each group is scaled by one sqrt(tdf / chi^2_tdf),
# a multivariate t draw, so that a student() fit has a moderate nu. With
# gaussian data nu runs to about 1e9, where the two programs part by about
# 1e-5 nats (dev/arspat-out/l1-all.txt).
ac_data <- function(ng = 12, nt = 6, seed = 29, tdf = NULL, re = 0) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  t = rep(seq_len(nt), ng))
  n <- nrow(d)
  d$x <- rnorm(n)
  d$x2 <- rnorm(n)
  d$h <- factor(rep(seq_len(ng / 2), each = 2 * nt))
  e <- unlist(lapply(seq_len(ng), function(i) {
    as.numeric(arima.sim(list(ar = 0.5), nt, sd = 0.7))
  }))
  if (!is.null(tdf)) e <- e * rep(sqrt(tdf / rchisq(ng, tdf)), each = nt)
  d$s <- runif(n, 0.2, 0.6)
  d$y <- 1 + 0.7 * d$x + e + rnorm(n, 0, d$s) + rnorm(ng / 2)[d$h] * re
  d
}

# a symmetric 4-nearest-neighbor graph, row-standardized: real
# eigenvalues, which brms takes
sar_W <- function(n) {
  set.seed(5)
  xy <- cbind(runif(n), runif(n))
  D <- as.matrix(stats::dist(xy))
  diag(D) <- Inf
  K <- t(apply(D, 1, function(r) as.numeric(rank(r) <= 4)))
  K <- pmax(K, t(K))
  K / rowSums(K)
}

cases <- list(
  se_ar = list(bf = "y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)"),
  se_ar_het = list(bf = "bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE), sigma ~ x2)"),
  se_ma = list(bf = "y | se(s, sigma = TRUE) ~ x + ma(t, g, cov = TRUE)"),
  se_arma = list(bf = "y | se(s, sigma = TRUE) ~ x + arma(t, g, cov = TRUE)"),
  se_cosy = list(bf = "y | se(s, sigma = TRUE) ~ x + cosy(t, g)"),
  se_unstr = list(bf = "y | se(s, sigma = TRUE) ~ x + unstr(t, g)",
                  drop = c(3, 10, 11, 40), nt = 4, ng = 24),
  se_ar_re = list(bf = "y | se(s, sigma = TRUE) ~ x + (1 | h) + ar(t, g, cov = TRUE)",
                  re = 0.8),
  se_ar_student = list(bf = "y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)",
                       family = "student()", tdf = 4, ng = 30),
  se_cosy_student = list(bf = "y | se(s, sigma = TRUE) ~ x + cosy(t, g)",
                         family = "student()", tdf = 4, ng = 30),
  mi_ar = list(bf = "y | mi() ~ x + ar(t, g, cov = TRUE)", mi = c(3, 20, 41),
               marg = "laplace"),
  mi_cosy = list(bf = "y | mi() ~ x + cosy(t, g)", mi = c(3, 20, 41),
                 marg = "laplace"),
  mi_unstr = list(bf = "y | mi() ~ x + unstr(t, g)", mi = c(3, 20, 41),
                  marg = "laplace"),
  mi_ar_student = list(bf = "y | mi() ~ x + ar(t, g, cov = TRUE)",
                       mi = c(3, 20, 41), family = "student()",
                       marg = "quad", tdf = 4, ng = 30),
  mi_fcor = list(bf = "y | mi() ~ x + fcor(M)", mi = c(3, 20, 41),
                 marg = "laplace", data2 = "M"),
  mi_sar_error = list(bf = "y | mi() ~ x + sar(W, type = 'error')",
                      mi = c(3, 20, 41), marg = "laplace", data2 = "W"),
  mi_sar_lag = list(bf = "y | mi() ~ x + sar(W, type = 'lag')",
                    mi = c(3, 20, 41), marg = "laplace", data2 = "W")
)

run_case <- function(nm) {
  cs <- cases[[nm]]
  d <- ac_data(ng = cs$ng %||% 12, nt = cs$nt %||% 6, tdf = cs$tdf,
               re = cs$re %||% 0)
  if (!is.null(cs$drop)) d <- d[-cs$drop, ]
  if (!is.null(cs$mi)) d$y[cs$mi] <- NA
  data2 <- if (identical(cs$data2, "M")) {
    list(M = 0.4^abs(outer(seq_len(nrow(d)), seq_len(nrow(d)), "-")))
  } else if (identical(cs$data2, "W")) {
    list(W = sar_W(nrow(d)))
  }
  fam_txt <- cs$family %||% "gaussian()"
  bf_frm <- eval(parse(text = if (startsWith(cs$bf, "bf(")) cs$bf else
    paste0("bf(", cs$bf, ")")), envir = asNamespace("frmtmb"))
  bf_brms <- eval(parse(text = if (startsWith(cs$bf, "bf(")) cs$bf else
    paste0("bf(", cs$bf, ")")), envir = asNamespace("brms"))
  fam_frm <- eval(parse(text = fam_txt), envir = asNamespace("frmtmb"))
  fam_brms <- eval(parse(text = fam_txt), envir = asNamespace("brms"))
  args_frm <- list(bf_frm, family = fam_frm, data = d)
  if (!is.null(data2)) args_frm$data2 <- data2
  fit <- do.call(frm, args_frm)
  cat("frmtmb logLik", format(as.numeric(logLik(fit)), digits = 12),
      "convergence", fit$opt$convergence, "\n")
  spec <- list(brms_args = if (!is.null(data2)) list(data2 = data2) else
    list())
  prior <- eq_env$eq_brms_prior(spec, bf_brms, fam_brms, d)
  joint <- length(fit$obj$env$random) > 0L
  refit <- function(k) {
    a <- args_frm
    a$control <- frmtmb_control(optimizer = eq_env$eq_step_optimizer(
      c(1, -1)[k]), restarts = 0, check_se = "ignore")
    do.call(frm, a)
  }
  if (is.null(cs$marg)) {
    l1 <- eq_env$eq_l1(spec, bf_brms, fam_brms, d, fit, prior, joint,
                       refit = refit)
    cat(sprintf("L1 const at estimate %.3e, grad %.3e\n", l1$at$const,
                l1$at$grad))
    for (k in seq_along(l1$extra)) {
      e <- l1$extra[[k]]
      if (!is.null(e$error)) cat("point", k, "ERROR", e$error, "\n") else
        cat(sprintf("L1 const at point %d %.3e\n", k, e$const))
    }
    cat("lp", format(l1$at$lp, digits = 12), "ours",
        format(l1$at$ours, digits = 12), "\n")
    if (!joint) {
      # the fit against brms's own maximum: flat priors and no Jacobian,
      # so rstan::optimizing() finds the maximum likelihood point
      l4 <- eq_env$eq_l4_opt(spec, bf_brms, fam_brms, d, prior, fit)
      cat(sprintf("L4 brms max %.10f frmtmb max %.10f diff %.3e (%d of 6 starts)\n",
                  l4$z[["brms_max"]], l4$z[["frm_max"]], l4$z[["diff"]],
                  l4$starts_ok))
    }
    return(invisible(l1))
  }
  # the marginal over the missing values, brms's program on every row
  args <- eq_env$eq_brms_args(spec, d, prior)
  code <- do.call(brms::make_stancode, c(list(bf_brms, family = fam_brms),
                                         args))
  sdat <- suppressMessages(do.call(brms::make_standata,
                                   c(list(bf_brms, family = fam_brms), args)))
  rtab <- do.call(eq_env$brms_ranef_table,
                  c(list(bf_brms, d, fam_brms, prior), spec$brms_args))
  sf <- suppressMessages(rstan::sampling(eq_env$brms_stan_model(code),
                                         data = sdat, chains = 0))
  # every Ymi, including those that left frm()'s objective, which
  # brms_inner_index() leaves out on purpose
  ymi_index <- function(pars) {
    b <- pars
    attr(b, "mi_gone") <- NULL
    b[["Ymi"]] <- pars[["Ymi"]] + 0.0123
    which(rstan::unconstrain_pars(sf, b) != rstan::unconstrain_pars(sf, pars))
  }
  one <- function(f) {
    pars <- eq_env$stan_pars_from_fit(f, sdat, code, rtab)
    ours <- -f$obj$fn(f$opt$par)
    ix <- ymi_index(pars)
    if (cs$marg == "laplace") {
      u0 <- rstan::unconstrain_pars(sf, pars)
      g <- function(u) rstan::grad_log_prob(sf, u,
                                            adjust_transform = FALSE)[ix]
      h <- 0.01
      H <- vapply(seq_along(ix), function(j) {
        up <- u0; dn <- u0
        up[ix[j]] <- up[ix[j]] + h
        dn[ix[j]] <- dn[ix[j]] - h
        (g(up) - g(dn)) / (2 * h)
      }, numeric(length(ix)))
      H <- (H + t(H)) / 2
      u <- u0
      u[ix] <- u[ix] - solve(H, g(u0))
      lp <- rstan::log_prob(sf, u, adjust_transform = FALSE)
      ld <- determinant(-H, logarithm = TRUE)
      val <- lp + length(ix) / 2 * log(2 * pi) - 0.5 * as.numeric(ld$modulus)
      return(c(stan = val, ours = ours, const = val - ours))
    }
    # one-dimensional quadrature per missing value, each in its own group
    u0 <- rstan::unconstrain_pars(sf, pars)
    lp0 <- rstan::log_prob(sf, u0, adjust_transform = FALSE)
    tot <- lp0
    for (j in ix) {
      fj <- function(v) vapply(v, function(vv) {
        u <- u0
        u[j] <- vv
        exp(rstan::log_prob(sf, u, adjust_transform = FALSE) - lp0)
      }, 0)
      I <- stats::integrate(fj, -Inf, Inf, rel.tol = 1e-12)
      tot <- tot + log(I$value)
    }
    c(stan = tot, ours = ours, const = tot - ours)
  }
  pts <- list(one(fit))
  for (k in 1:2) pts[[k + 1L]] <- one(refit(k))
  for (k in seq_along(pts)) {
    cat(sprintf("marginal point %d: stan %.10f ours %.10f const %.3e\n",
                k, pts[[k]][["stan"]], pts[[k]][["ours"]],
                pts[[k]][["const"]]))
  }
  invisible(pts)
}

which_cases <- commandArgs(trailingOnly = TRUE)
if (!length(which_cases)) which_cases <- names(cases)
for (nm in which_cases) {
  cat("\n######## ", nm, "\n", sep = "")
  r <- tryCatch(run_case(nm), error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n")
  })
}
