# Within-group residual correlation: brms's "R-side" autocorrelation
# terms ar(), ma(), arma(), cosy() and unstr().
#
# WHAT THIS IS. A random-effect (G-side) structure such as
# `ar1(week + 0 | subj)` gives every subject a latent AR(1) curve and
# leaves the residual iid. An R-side term instead makes the RESIDUAL of
# one group a single multivariate draw,
#
#   y_g ~ N(mu_g, D_g R(phi) D_g),   D_g = diag(sigma_i, i in g),
#
# with R a correlation matrix over the group's time points. Nothing is
# added to the linear predictor and no new latent vector appears; the
# likelihood itself changes shape. That is what nlme::gls(correlation =
# corAR1()) fits, and what brms fits under `cov = TRUE`.
#
# PARAMETERIZATION OF sigma. R is UNIT-DIAGONAL here, so `sigma` is the
# MARGINAL residual standard deviation, as it is in nlme and as it is
# everywhere else in this package (sigma(), pearson residuals, se(),
# every distributional model). brms diverges: its
# `cholesky_cor_ar1(ar, n)` returns `chol(rho^|i-j| / (1 - rho^2))` and
# `cholesky_cor_ma1` returns `chol(...)` with `1 + ma^2` on the
# diagonal, so under `ar()`/`ma()`/`arma()` brms's `sigma` is the
# INNOVATION sd. The two agree after
#
#   sigma_marginal = sigma_innovation / sqrt(1 - phi^2)      (AR(1))
#   sigma_marginal = sigma_innovation * sqrt(1 + theta^2)    (MA(1))
#
# and the correlation parameters themselves agree exactly. brms's
# `cosy()` and `unstr()` are already unit-diagonal, so those match
# without conversion.
#
# LAGS COME FROM THE TIME LEVELS. The correlation between two rows is
# read off their positions in the GLOBAL set of time levels, so a group
# missing week 3 gets cor(week2, week4) = rho^2, not rho. brms indexes
# `ar`/`ma`/`arma`/`cosy` by position WITHIN the group (its Stan code
# takes `chol_cor[1:nobs[i], 1:nobs[i]]`, with no time index at all;
# only `unstr` carries `Jtime_tg`), so on ragged data brms treats a
# missing row as no gap. The gap-aware reading is nlme's
# (`corAR1(form = ~ week | subj)`) and is what the validation suite
# pins against gls(); on complete balanced groups the two coincide.
#
# COST. The density is evaluated one PATTERN at a time, a pattern being
# a distinct set of present time levels. Groups sharing a pattern share
# one k x k correlation submatrix and are done in a single vectorized
# `dmvnorm` over a k x G matrix of standardized residuals, so a
# balanced design costs exactly one on-tape Cholesky per evaluation.
# Ragged data costs one per distinct pattern. Nothing ever builds an
# n x n matrix.

#' Within-group residual correlation (R-side autocorrelation)
#'
#' `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` are written as terms
#' of the model formula, next to the fixed and random effects, and make
#' the residuals of one group a single correlated draw instead of
#' independent ones:
#'
#' \deqn{y_g \sim N(\mu_g,\; D_g R D_g), \qquad
#'       D_g = \mathrm{diag}(\sigma_i, i \in g).}
#'
#' `R` is a unit-diagonal correlation matrix over the time points, so
#' `sigma` keeps its usual meaning - the marginal residual standard
#' deviation - and a `sigma ~ ...` distributional model enters through
#' the diagonal. Nothing is added to the linear predictor, so
#' [fitted()], [predict.frmtmb_fit()] and `se.fit` are unchanged; what
#' changes is the likelihood. This is the model `nlme::gls(correlation
#' = corAR1())` fits, and the one brms fits under `cov = TRUE`.
#'
#' @section Structures:
#' \describe{
#'   \item{`ar(time, gr, p = 1, cov = TRUE)`}{Stationary AR(`p`)
#'     correlation. `p = 1` gives \eqn{R_{ij} = \rho^{|i-j|}}.}
#'   \item{`ma(time, gr, q = 1, cov = TRUE)`}{Invertible MA(`q`):
#'     correlation dies after lag `q`.}
#'   \item{`arma(time, gr, p = 1, q = 1, cov = TRUE)`}{Stationary,
#'     invertible ARMA(`p`, `q`).}
#'   \item{`cosy(time, gr)`}{Compound symmetry - one correlation shared
#'     by every pair. Equivalent in fit to a random intercept per group,
#'     but the correlation may also be negative (down to
#'     `-1 / (d - 1)`), which a variance component cannot express.}
#'   \item{`unstr(time, gr)`}{One free correlation per pair of time
#'     levels, through the same Cholesky parameterization the `us()`
#'     random-effect structure uses.}
#' }
#' The AR and MA coefficients are estimated through partial
#' autocorrelations (the Monahan/Jones transform `nlme::corARMA` also
#' uses), so every parameter value is a stationary and invertible
#' process and the optimizer cannot leave the parameter space. The
#' ARMA autocorrelation function is exact, not a truncated MA(\eqn{\infty})
#' expansion.
#'
#' @section Arguments:
#' The argument order is brms's, so the FIRST positional argument is
#' `time` and the second is `gr`: write `cosy(gr = subj)`, not
#' `cosy(subj)`.
#' \describe{
#'   \item{`time`}{The variable whose levels index the correlation. Omit
#'     it (brms's `time = NA`) to use each row's position within its
#'     group.}
#'   \item{`gr`}{The grouping variable the residual factorizes over.
#'     Omit it to treat the whole data set as one series.}
#'   \item{`p`, `q`}{Autoregressive and moving-average orders.}
#'   \item{`cov`}{Must be `TRUE`. brms's default `cov = FALSE` is a
#'     different likelihood (a residual regression that conditions on
#'     each group's first rows), which is not implemented; the call is
#'     refused rather than silently reinterpreted.}
#' }
#'
#' @section Families:
#' `gaussian()` and `student()` only - the two families with a real
#' residual, and exactly the two brms treats this way (`student()` gets
#' the multivariate-t analog, with one `nu` per group, so a predicted
#' `nu ~ ...` is refused). brms accepts the same spelling for other
#' families but fits a different model there: a latent Gaussian AR
#' process added to the linear predictor. That model has a spelling of
#' its own here - a random effect over the time factor, `+
#' ar1(factor(week) + 0 | subj)`, or `toep()` / `us()` for a freer lag
#' structure - and the refusal names it.
#'
#' @section Time points, gaps and ragged groups:
#' The lag between two rows is the distance between their positions in
#' the GLOBAL set of time levels, so a group missing week 3 gets
#' `cor(week2, week4) = rho^2`. That is `nlme`'s reading
#' (`corAR1(form = ~ week | subj)`) and what the agreement tests pin
#' down; brms instead indexes `ar()`/`ma()`/`arma()`/`cosy()` by
#' position within the group and treats a missing row as no gap (only
#' its `unstr()` carries a time index). On complete balanced groups the
#' two coincide. Time levels that are whole numbers but not consecutive
#' warn, because the lag is then not the one the labels suggest.
#'
#' Groups of different sizes are handled by construction: the density is
#' evaluated one PATTERN at a time (a pattern being a distinct set of
#' present time levels), so a balanced design costs one Cholesky per
#' gradient evaluation and ragged data costs one per distinct pattern.
#' Duplicate `(gr, time)` pairs are refused, as they are in brms.
#'
#' @section What cannot be combined with it:
#' The likelihood is a joint density over each group, so it no longer
#' factorizes into per-row contributions. `weights()`, `cens()`,
#' `trunc()`, `se()` and `mi()` on the same response are therefore
#' refused, as are `rescor = TRUE`, mixture families,
#' `quadrature = TRUE`, `frm_simulate()` and
#' `residuals(type = "osa")`. brms refuses the same core set. Random
#' effects ARE allowed and are the point of the feature: the marginal
#' likelihood is a Laplace approximation over the modes with the
#' correlated residual density inside, which reproduces
#' `nlme::lme(random = ~ 1 | subj, correlation = corAR1())`.
#'
#' @section Where the parameters appear:
#' On the internal (unconstrained) scale they are the `thetaac_*` rows
#' of [confint.frmtmb_fit()] and of `vcov(full = TRUE)`. On the natural
#' scale, `summary()` prints them under "Within-group residual
#' correlation" and [confint_varcorr()] reports one row per parameter
#' under brms's names - `ar[1]`, `ma[1]`, `cosy`,
#' `cortime__<t1>__<t2>` - with a delta-method interval (Fisher-z for
#' the bounded ones, the identity for the coefficients of a
#' higher-order AR/MA process). [hypothesis()] sees them as `ar1`,
#' `ma1`, `cosy` and `cortime__<t1>__<t2>`. [autocor_matrix()] returns
#' the fitted `R`. They are NOT part of [VarCorr()], which reports
#' random-effect blocks. `set_prior()` cannot target them yet.
#'
#' @section Divergence from brms:
#' brms parameterizes `ar()`, `ma()` and `arma()` by the INNOVATION
#' standard deviation (its `cholesky_cor_ar1()` divides by
#' `1 - ar^2`), while `cosy()` and `unstr()` use the marginal one. Here
#' every structure uses the marginal `sigma`, so that `sigma()`, pearson
#' residuals and a `sigma ~ x` model mean one thing throughout. The
#' correlation parameters agree with brms exactly; the scales relate by
#' `sigma_marginal = sigma_innovation / sqrt(1 - phi^2)` for AR(1) and
#' `sigma_marginal = sigma_innovation * sqrt(1 + theta^2)` for MA(1).
#' brms also limits `cov = TRUE` to order one ("Covariance formulation
#' of ARMA structures is only possible for effects of maximal order
#' one"); higher `p` and `q` are supported here.
#'
#' @return `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` are formula
#'   terms, not free-standing functions: `bf()` reads them at parse
#'   time, and the value they contribute is the fitted autocorrelation
#'   block of the model, reachable through [autocor_matrix()],
#'   [VarCorr()] and [confint_varcorr()]. This page itself documents
#'   the term grammar and returns nothing.
#' @name frmtmb-autocor
#' @seealso [autocor_matrix()] for the fitted correlation matrix,
#'   [confint_varcorr()] for natural-scale intervals, and
#'   `vignette("brms-migration")` for the porting notes.
#' @examples
#' set.seed(1)
#' d <- expand.grid(week = 1:5, subj = factor(1:30))
#' d$x <- rnorm(150)
#' e <- as.vector(vapply(1:30, function(i) {
#'   as.vector(stats::filter(rnorm(5), 0.6, "recursive"))
#' }, numeric(5)))
#' d$y <- 1 + 0.5 * d$x + e
#'
#' fit <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
#'            data = d)
#' summary(fit)
#' autocor_matrix(fit)
#'
#' # compound symmetry, and the unstructured correlation over the five
#' # weeks
#' frm(bf(y ~ x + cosy(week, subj)) + gaussian(), data = d)
#' frm(bf(y ~ x + unstr(week, subj)) + gaussian(), data = d)
#'
#' # a random intercept alongside the correlated residual is allowed
#' frm(bf(y ~ x + (1 | subj) + ar(week, subj, cov = TRUE)) + gaussian(),
#'     data = d)
#' @family autocorrelation
NULL

# The brms spellings this file understands.
autocor_structs <- c("ar", "ma", "arma", "cosy", "unstr")

# Largest number of time levels a residual block may carry. The density
# factorizes a dense d x d Cholesky on every tape evaluation, so this
# caps the per-gradient cost; it is also the size at which a "one group
# holding the whole data set" mistake stops being cheap.
autocor_max_dim <- 300L

# unstr() estimates d (d - 1) / 2 free correlations, so its cap is much
# lower than the density's: 50 levels is already 1225 parameters.
autocor_max_unstr <- 50L

#' Evaluate an `ar()`/`ma()`/`arma()` tuning argument (`p`, `q`, `cov`).
#'
#' @noRd
autocor_arg <- function(expr, nm, env, fn) {
  val <- tryCatch(eval(expr, env), error = function(e) {
    stop(fn, "(): the argument ", nm, " = ", deparse1(expr),
         " could not be evaluated in the formula environment: ",
         conditionMessage(e), call. = FALSE)
  })
  if (nm == "cov") {
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop(fn, "(): cov must be TRUE or FALSE", call. = FALSE)
    }
    return(isTRUE(val))
  }
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val) ||
      val < 0 || val != trunc(val)) {
    stop(fn, "(): ", nm, " must be a single non-negative whole number ",
         "(got ", deparse1(expr), ")", call. = FALSE)
  }
  as.integer(val)
}

#' brms's `ar(time, gr, p, cov)` and relatives, parsed into the spec
#' entry the frame builds a residual block from.
#'
#' The argument ORDER is brms's: the first positional argument is
#' `time`, the second is `gr`. `cosy(g)` therefore names a time
#' variable, not a group - the same trap brms has - so the duplicate
#' check below points at `cosy(gr = g)`.
#'
#' @noRd
parse_autocor_call <- function(tm, env) {
  fn <- as.character(tm[[1L]])[1L]
  argn <- switch(fn,
                 ar = c("time", "gr", "p", "cov"),
                 ma = c("time", "gr", "q", "cov"),
                 arma = c("time", "gr", "p", "q", "cov"),
                 c("time", "gr"))
  a <- match_special_args(tm, argn, fn)
  na_arg <- function(e) is.null(e) || identical(e, quote(NA)) ||
    (is.logical(e) && length(e) == 1L && is.na(e))
  time_expr <- if (na_arg(a$time)) NULL else a$time
  gr_expr <- if (na_arg(a$gr)) NULL else a$gr
  if (fn == "unstr" && (is.null(time_expr) || is.null(gr_expr))) {
    stop("unstr() needs both a time variable and a grouping variable: ",
         "unstr(week, subj)", call. = FALSE)
  }
  p <- if (fn %in% c("ar", "arma")) {
    if (is.null(a$p)) 1L else autocor_arg(a$p, "p", env, fn)
  } else 0L
  q <- if (fn %in% c("ma", "arma")) {
    if (is.null(a$q)) 1L else autocor_arg(a$q, "q", env, fn)
  } else 0L
  if (fn %in% c("ar", "ma", "arma")) {
    if (p + q < 1L) {
      stop(fn, "(): at least one of p and q must be greater than zero",
           call. = FALSE)
    }
    cov <- if (is.null(a$cov)) FALSE else {
      autocor_arg(a$cov, "cov", env, fn)
    }
    if (!cov) {
      # brms's cov = FALSE is a different likelihood, not a different
      # implementation of the same one: it adds ar * (y - mu) at the
      # previous row to the linear predictor and keeps a univariate
      # normal density, which conditions on the first p rows instead of
      # giving them their stationary distribution. Fitting the
      # covariance form under that spelling would silently disagree
      # with brms; refusing says so.
      stop(fn, "(): only the residual-covariance formulation is ",
           "implemented, so the call needs cov = TRUE: ", fn, "(",
           if (!is.null(time_expr)) paste0(deparse1(time_expr), ", ") else "",
           if (!is.null(gr_expr)) paste0(deparse1(gr_expr), ", ") else "",
           "cov = TRUE). brms's default cov = FALSE is the ",
           "residual-regression form, a different likelihood (it ",
           "conditions on the first observations of each group rather ",
           "than giving them their stationary distribution). cov = ",
           "TRUE is the marginal multivariate-normal residual that ",
           "nlme::gls(correlation = corAR1()) fits", call. = FALSE)
    }
  }
  list(fn = fn, struct = fn, time_expr = time_expr, gr_expr = gr_expr,
       p = p, q = q, label = deparse1(tm))
}

#' Free correlations of a `d x d` unstructured matrix.
#'
#' Its own function because `d * (d - 1L) %/% 2L` does NOT say this:
#' `%/%` binds tighter than `*`, so that expression is
#' `d * ((d - 1) %/% 2)`, which is right for odd `d` and wrong for even
#' `d` - a block that then silently carries too few parameters.
#'
#' @noRd
autocor_n_cor <- function(d) as.integer(d * (d - 1L) / 2L)

#' Number of `thetaac` entries a residual block needs.
#'
#' @noRd
autocor_npar <- function(ac) {
  switch(ac[["struct"]],
         ar = ac[["p"]], ma = ac[["q"]], arma = ac[["p"]] + ac[["q"]],
           cosy = 1L,
         unstr = autocor_n_cor(ac[["d"]]))
}

#' Starting `thetaac`: white noise for the ARMA family and the identity
#' for unstr; `cosy` starts where its own homcs-style transform puts a
#' zero, as the `homcs` random-effect block does.
#'
#' @noRd
autocor_start <- function(ac) numeric(autocor_npar(ac))

# ------------------------------------------------------- ACF machinery
#
# All AD-safe: arithmetic and one small linear solve, no branching on
# parameter values, no base matrix()/c().

#' Unconstrained reals -> partial autocorrelations in (-1, 1). The same
#' bounded transform the ar1() and toep() random-effect structures use,
#' so a rho reads the same wherever it appears.
#'
#' @noRd
autocor_pacf <- function(th) th / sqrt(1 + th^2)

#' Partial autocorrelations -> AR (or MA) coefficients, by the
#' Levinson-Durbin recursion. Every input maps onto a stationary AR
#' polynomial (an invertible MA polynomial), which is the
#' Monahan/Jones parameterization nlme's corARMA also uses, so the
#' optimizer cannot step outside the parameter space.
#'
#' @noRd
autocor_levinson <- function(pac) {
  "c" <- RTMB::ADoverload("c")
  phi <- pac[1]
  k <- 1L
  while (k < length(pac)) {
    k <- k + 1L
    pk <- pac[k]
    prev <- seq_len(k - 1L)
    phi <- c(phi[prev] - pk * phi[rev(prev)], pk)
  }
  phi
}

#' Autocorrelation function of a stationary ARMA(p, q) process at lags
#' `0..lagmax`, exactly (Brockwell & Davis 3.3.1): the psi weights of
#' the MA(inf) expansion up to lag q, then the `m x m` linear system
#' with `m = max(p, q) + 1` for gamma(0..m-1), then the AR recursion for
#' the rest. No truncation and no infinite sum, so the tape carries the
#' exact ACF; verified against `stats::ARMAacf()` to 1e-16.
#'
#' @noRd
autocor_arma_acf <- function(phi, th, lagmax) {
  "[<-" <- RTMB::ADoverload("[<-")
  p <- length(phi)
  q <- length(th)
  zero <- if (p) phi[1] * 0 else th[1] * 0
  m <- max(p, q) + 1L
  # psi weights psi_0 .. psi_q (psi[j + 1] is psi_j)
  psi <- rep(zero, q + 1L)
  psi[1] <- zero + 1
  if (q) {
    for (j in seq_len(q)) {
      s <- th[j]
      if (p) {
        for (i in seq_len(min(j, p))) s <- s + phi[i] * psi[j - i + 1L]
      }
      psi[j + 1L] <- s
    }
  }
  # gamma(k) - sum_j phi_j gamma(|k - j|) = sum_{j >= k} theta_j psi_{j-k}
  Av <- rep(zero, m * m)
  rhs <- rep(zero, m)
  for (k in 0:(m - 1L)) {
    Av[k * m + k + 1L] <- Av[k * m + k + 1L] + 1
    if (p) {
      for (j in seq_len(p)) {
        pos <- abs(k - j) * m + k + 1L
        Av[pos] <- Av[pos] - phi[j]
      }
    }
    s <- if (k == 0L) psi[1L] else zero
    if (q) {
      for (j in seq_len(q)) {
        if (j >= k) s <- s + th[j] * psi[j - k + 1L]
      }
    }
    rhs[k + 1L] <- s
  }
  # RTMB::solve, not base's: the S4 advector method is not imported
  g <- RTMB::solve(RTMB::matrix(Av, m, m)) %*% RTMB::matrix(rhs, m, 1)
  acv <- rep(zero, lagmax + 1L)
  for (k in 0:min(m - 1L, lagmax)) acv[k + 1L] <- g[k + 1L, 1]
  if (lagmax > m - 1L) {
    for (k in m:lagmax) {
      s <- zero
      if (p) for (j in seq_len(p)) s <- s + phi[j] * acv[k - j + 1L]
      acv[k + 1L] <- s
    }
  }
  acv / acv[1L]
}

#' Banded correlation matrix from an autocorrelation vector (lag 0
#' first). The gather keeps the advector class that `matrix()` would
#' strip.
#'
#' @noRd
autocor_band <- function(acv, d) {
  M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
  RTMB::matrix(acv[as.vector(M)], d, d)
}

#' The `d x d` unit-diagonal correlation matrix of a residual block at
#' one `thetaac` segment. AD-safe; also correct on plain numerics, which
#' is what the reporting and simulation paths feed it.
#'
#' @noRd
autocor_cor <- function(theta, ac) {
  d <- ac[["d"]]
  if (ac[["struct"]] == "unstr") return(us_chol_cor(theta, d))
  if (ac[["struct"]] == "cosy") {
    # bounded on (-1/(d-1), 1), the widest range that keeps every
    # sub-block positive definite - nlme's corCompSymm range. brms
    # bounds cosy on [0, 1] instead, so a negative estimate here has no
    # brms counterpart.
    a <- 1 / (d - 1)
    rho <- -a + (1 + a) / (1 + exp(-theta[1]))
    return(diag(d) * (1 - rho) + rho)
  }
  pac <- autocor_pacf(theta)
  phi <- if (ac[["p"]]) autocor_levinson(pac[seq_len(ac[["p"]])]) else pac[0]
  mth <- if (ac[["q"]]) autocor_levinson(pac[ac[["p"]] +
                                         seq_len(ac[["q"]])]) else pac[0]
  autocor_band(autocor_arma_acf(phi, mth, d - 1L), d)
}

#' Natural-scale parameters of a residual block, named as brms names
#' them (`ar[1]`, `ma[1]`, `cosy`, `cortime__<t1>__<t2>`).
#'
#' @noRd
autocor_natural <- function(theta, ac) {
  if (ac[["struct"]] == "unstr") {
    C <- us_chol_cor(theta, ac[["d"]])
    pr <- which(lower.tri(C), arr.ind = TRUE)
    return(stats::setNames(
      C[lower.tri(C)],
      paste0("cortime__", ac[["time_levels"]][pr[, 2]], "__",
             ac[["time_levels"]][pr[, 1]])))
  }
  if (ac[["struct"]] == "cosy") {
    a <- 1 / (ac[["d"]] - 1)
    return(c(cosy = -a + (1 + a) / (1 + exp(-theta[1]))))
  }
  pac <- autocor_pacf(theta)
  out <- numeric(0)
  if (ac[["p"]]) {
    out <- c(out, stats::setNames(autocor_levinson(pac[seq_len(ac[["p"]])]),
                                  paste0("ar[", seq_len(ac[["p"]]), "]")))
  }
  if (ac[["q"]]) {
    out <- c(out, stats::setNames(
      autocor_levinson(pac[ac[["p"]] + seq_len(ac[["q"]])]),
      paste0("ma[", seq_len(ac[["q"]]), "]")))
  }
  out
}

#' Which reporting transform each natural-scale parameter of a block
#' takes: Fisher-z for anything that lives inside (-1, 1), and the
#' identity for the AR/MA coefficients of a higher-order process, which
#' do not (a stationary AR(2) has phi_1 in (-2, 2)).
#'
#' @noRd
autocor_types <- function(ac) {
  if (ac[["struct"]] == "unstr") {
    return(rep("cor", autocor_n_cor(ac[["d"]])))
  }
  if (ac[["struct"]] == "cosy") return("cor")
  c(rep(if (ac[["p"]] <= 1L) "cor" else "raw", ac[["p"]]),
    rep(if (ac[["q"]] <= 1L) "cor" else "raw", ac[["q"]]))
}

# ---------------------------------------------------------- assembly
#
# Everything below runs off the tape, once, at frame assembly.

#' Everything a residual correlation term cannot be combined with.
#'
#' Each refusal is a likelihood that stops factorizing over rows once
#' the residual is a joint density, so accepting the pair would fit a
#' different model in silence. brms refuses the same set from
#' `stan_log_lik_gaussian_time()`, which rejects `weights()`, `cens()`
#' and `trunc()` outright ("Invalid addition arguments for this model")
#' and `rescor` from `frame_ac()` ("Explicit covariance terms cannot be
#' modeled when 'rescor' is estimated at the same time").
#'
#' @noRd
check_autocor_response <- function(resp, spec, av, yv) {
  ac <- resp$autocor
  fam <- resp$family[["family"]]
  if (!fam %in% c("gaussian", "student")) {
    stop(ac[["fn"]], "(): a residual correlation needs a family with real ",
         "residuals, so gaussian() or student(); '", fam, "' has none. ",
         "brms accepts the same call for other families but fits a ",
         "different model there - a latent gaussian AR process added to ",
         "the linear predictor - which is spelled as a random effect ",
         "here: replace ", ac[["label"]], " with ar1(factor(",
         if (is.null(ac[["time_expr"]])) "time" else
           deparse1(ac[["time_expr"]]),
         ") + 0 | ",
         if (is.null(ac[["gr_expr"]])) "group" else deparse1(ac[["gr_expr"]]),
         "), or toep()/us() for a freer lag structure", call. = FALSE)
  }
  if (is.matrix(yv)) {
    stop(ac[["fn"]], "(): the response must be a numeric vector",
         call. = FALSE)
  }
  if (isTRUE(spec$rescor)) {
    stop(ac[["fn"]], "(): residual correlation within a response cannot be ",
         "combined with rescor = TRUE. Both describe the residual ",
         "covariance - one across time, one across responses - and the ",
         "joint structure is their Kronecker product, which is not ",
         "implemented. brms refuses the same pair. Fit rescor = FALSE, ",
         "or drop the time term", call. = FALSE)
  }
  if (!is.null(resp$family[["mix"]])) {
    stop(ac[["fn"]], "(): a mixture likelihood has no single residual to ",
         "correlate", call. = FALSE)
  }
  bad <- c(if (!is.null(av[["weights"]])) "weights()",
           if (!is.null(av[["cens"]])) "cens()",
           if (!is.null(av[["trunc_lb"]]) ||
             !is.null(av[["trunc_ub"]])) "trunc()",
           if (!is.null(av[["se"]])) "se()",
           if (isTRUE(resp$aterms[["mi"]])) "mi()")
  if (length(bad)) {
    stop(ac[["fn"]], "(): ", paste(bad, collapse = ", "),
         " cannot be combined with a residual correlation term. The ",
         "likelihood is a joint density over each group, so it no ",
         "longer factorizes into per-row contributions that a ",
         "frequency weight could repeat, a censoring indicator could ",
         "replace with a tail probability, a truncation bound could ",
         "renormalize, or a known standard error could be added to. ",
         "brms refuses weights(), cens() and trunc() here for the same ",
         "reason", call. = FALSE)
  }
  if (fam == "student") {
    # brms's student_t_time_*_lpdf takes `real nu`: the multivariate t
    # has one shape for the whole group, so a row-varying nu has no
    # multivariate counterpart
    nud <- resp$dpars[["nu"]]
    predicted <- !is.null(nud) && is.null(nud$constant) &&
      (length(nud$re %||% list()) || length(nud$smooth %||% list()) ||
         !identical(deparse1(reformulas::RHSForm(nud$fixed)), "1"))
    if (predicted) {
      stop(ac[["fn"]], "(): student() with a residual correlation needs a ",
           "constant nu; a predicted 'nu ~ ...' has no multivariate-t ",
           "counterpart, because the group shares one shape parameter",
           call. = FALSE)
    }
  }
  ac
}

#' Time-level index of every row plus the level labels.
#'
#' With no `time` variable the index is the row's POSITION within its
#' group (brms's default when `time = NA`); with one, it is the row's
#' position in the global, sorted set of levels, so a gap in a group is
#' a gap in the lag.
#'
#' @noRd
autocor_time_index <- function(ac, mf, gidx, env, n) {
  if (is.null(ac[["time_expr"]])) {
    idx <- integer(n)
    for (g in split(seq_len(n), gidx)) idx[g] <- seq_along(g)
    return(list(idx = idx, levels = as.character(seq_len(max(idx)))))
  }
  tv <- eval(ac[["time_expr"]], mf, env)
  if (anyNA(tv)) {
    stop(ac[["fn"]], "(): the time variable '", deparse1(ac[["time_expr"]]),
         "' has missing values", call. = FALSE)
  }
  if (is.factor(tv)) {
    lv <- levels(droplevels(tv))
    idx <- match(as.character(tv), lv)
  } else if (is.numeric(tv)) {
    lv <- sort(unique(as.numeric(tv)))
    idx <- match(as.numeric(tv), lv)
    lv <- format(lv, trim = TRUE)
  } else {
    lv <- sort(unique(as.character(tv)))
    idx <- match(as.character(tv), lv)
  }
  list(idx = idx, levels = lv)
}

#' Non-consecutive whole-number time levels mean the lag counted here
#' (one step per level) is not the lag the labels suggest.
#'
#' Unlike the `ar1()` random-effect warning this one CAN be acted on:
#' the gap is only invisible because the missing level is missing from
#' every group. Say so rather than let cor(t6, t10) be fitted as rho.
#'
#' @noRd
autocor_warn_gaps <- function(ac, lv) {
  pos <- suppressWarnings(as.numeric(lv))
  if (anyNA(pos) || any(pos != trunc(pos))) return(invisible(NULL))
  gap <- which(abs(diff(pos)) != 1)
  if (!length(gap)) return(invisible(NULL))
  i <- gap[1L]
  warning(ac[["fn"]], "(): the time levels present are whole numbers but not ",
          "consecutive ('", lv[i], "' is followed by '", lv[i + 1L],
          "'), and the lag is counted in LEVELS, so that gap counts as ",
          "a single step. Pad the level set (a row per missing time, ",
          "with an NA response and na.action = na.exclude, is enough) ",
          "or use ou(num_factor(", deparse1(ac[["time_expr"]]),
          ") + 0 | ...) for a continuous-position structure",
          call. = FALSE)
  invisible(NULL)
}

#' Turn one parsed autocorrelation term plus the model frame into the
#' residual block the objective evaluates.
#'
#' The returned entry carries the structure, the number of time levels,
#' and one entry per PATTERN (a distinct set of present time levels):
#' the flat row-index vector of a `k x G` matrix of residuals, and the
#' gather that lifts the pattern's `k x k` correlation submatrix out of
#' the full `d x d` one.
#'
#' @noRd
autocor_block <- function(ac, resp, mf, env, n) {
  fn <- ac[["fn"]]
  gidx <- if (is.null(ac[["gr_expr"]])) {
    rep(1L, n)
  } else {
    gv <- eval(ac[["gr_expr"]], mf, env)
    if (anyNA(gv)) {
      stop(fn, "(): the residual factorizes over gr = ",
           deparse1(ac[["gr_expr"]]),
           ", so every row needs a group; that variable has ",
           sum(is.na(gv)), " missing value(s)", call. = FALSE)
    }
    as.integer(factor(gv))
  }
  ti <- autocor_time_index(ac, mf, gidx, env, n)
  ac[["d"]] <- length(ti$levels)
  ac[["time_levels"]] <- ti$levels
  ac[["group_levels"]] <- if (is.null(ac[["gr_expr"]])) "1" else {
    levels(factor(eval(ac[["gr_expr"]], mf, env)))
  }
  if (ac[["d"]] < 2L) {
    stop(fn, "(): the residual correlation needs at least 2 time ",
         "points; '",
         if (is.null(ac[["time_expr"]])) "(row order)" else
           deparse1(ac[["time_expr"]]), "' gives ", ac[["d"]], call. = FALSE)
  }
  if (ac[["d"]] > autocor_max_dim) {
    stop(fn, "(): ", ac[["d"]], " time points would build a dense ", ac[["d"]],
         " x ", ac[["d"]], " residual covariance on every gradient ",
         "evaluation (cap ", autocor_max_dim,
         "). Name a grouping variable - ", fn, "(",
         if (!is.null(ac[["time_expr"]])) paste0(deparse1(ac[["time_expr"]]),
           ", ")
         else "", "gr = subject) - so the density factorizes over ",
         "groups", call. = FALSE)
  }
  if (ac[["struct"]] == "unstr" && ac[["d"]] > autocor_max_unstr) {
    stop("unstr(): ", ac[["d"]], " time points means ",
         autocor_n_cor(ac[["d"]]), " free correlations (cap ",
         autocor_max_unstr, " levels). Use ar(), arma() or cosy() for a ",
         "structure that does not grow with the number of time points",
         call. = FALSE)
  }
  key <- paste(gidx, ti$idx, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- key[duplicated(key)][1L]
    parts <- strsplit(dup, "\r", fixed = TRUE)[[1L]]
    stop(fn, "(): time points within groups must be unique; group '",
         ac[["group_levels"]][as.integer(parts[1L])], "' has ",
         sum(key == dup), " rows at time '",
         ti$levels[as.integer(parts[2L])], "'. ",
         if (is.null(ac[["gr_expr"]])) {
           paste0("The first argument of ", fn,
                  "() is the TIME variable and the second is the ",
                  "grouping variable, so ", fn, "(g) reads g as time; ",
                  "write ", fn, "(gr = g) if g is the group.")
         } else {
           "Aggregate the repeated rows, or add the replicate to the grouping variable."
         }, call. = FALSE)
  }
  if (!is.null(ac[["time_expr"]])) autocor_warn_gaps(ac, ti$levels)
  # one pattern per distinct set of present time levels
  by_g <- split(seq_len(n), gidx)
  pkeys <- vapply(by_g, function(rows) {
    paste(sort(ti$idx[rows]), collapse = ",")
  }, "")
  ac[["patterns"]] <- lapply(split(seq_along(by_g), pkeys), function(gs) {
    tset <- sort(ti$idx[by_g[[gs[1L]]]])
    k <- length(tset)
    rows <- vapply(gs, function(g) {
      r <- by_g[[g]]
      r[order(ti$idx[r])]
    }, integer(k))
    list(k = k, G = length(gs), rows = as.vector(rows),
         gather = as.vector(outer(tset, tset, function(a, b) {
           (b - 1L) * ac[["d"]] + a
         })))
  })
  names(ac[["patterns"]]) <- NULL
  ac[["npar"]] <- autocor_npar(ac)
  ac[["resp"]] <- resp
  ac[["n_groups"]] <- length(by_g)
  ac
}

# ---------------------------------------------------------- likelihood

#' Log-density of one response under a residual correlation block.
#'
#' Standardizing by `sigma` first is what keeps this vectorized when
#' `sigma` varies by row: `z_g = (y_g - mu_g) / sigma_g` has correlation
#' matrix R whatever the diagonal is, so one `dmvnorm` covers every
#' group sharing a pattern and the scale re-enters as `-sum(log sigma)`.
#'
#' `nu` (student) is a scalar: the multivariate t has one shape for the
#' whole group, which is also why a predicted `nu` is refused upstream.
#'
#' @noRd
autocor_loglik <- function(z, R, ac, log_sigma_sum, nu = NULL) {
  ll <- -log_sigma_sum
  for (pt in ac[["patterns"]]) {
    k <- pt$k
    Rp <- if (k == ac[["d"]]) R else {
      RTMB::matrix(as.vector(R)[pt$gather], k, k)
    }
    zs <- z[pt$rows]
    if (is.null(nu)) {
      ll <- ll + if (k == 1L) {
        sum(RTMB::dnorm(zs, 0, 1, log = TRUE))
      } else {
        sum(RTMB::dmvnorm(t(RTMB::matrix(zs, k, pt$G)), 0, Rp,
                          log = TRUE))
      }
      next
    }
    # multivariate Student-t with scale matrix R (brms's
    # multi_student_t_lpdf): log|R| comes from the gaussian density at
    # the origin, which avoids needing a Cholesky of an advector matrix
    Z <- RTMB::matrix(zs, k, pt$G)
    l0 <- RTMB::dmvnorm(0 * Z[, 1], 0, Rp, log = TRUE)
    half_ldet <- -(l0 + 0.5 * k * log(2 * pi))
    tZ <- t(Z)
    qv <- as.vector(((tZ %*% RTMB::solve(Rp)) * tZ) %*% rep(1, k))
    ll <- ll + pt$G * (lgamma((nu + k) / 2) - lgamma(nu / 2) -
                         0.5 * k * log(nu * pi) - half_ldet) -
      0.5 * (nu + k) * sum(log(1 + qv / nu))
  }
  ll
}

# ---------------------------------------------------------- reporting

#' `varcorr_trans_rows()` for the residual blocks: one row per natural
#' scale parameter, on the scale its interval is built on (Fisher-z for
#' the bounded ones, the identity for higher-order AR/MA coefficients),
#' with delta-method standard errors from the `thetaac` block of the
#' outer covariance.
#'
#' @noRd
autocor_trans_rows <- function(fit) {
  acs <- fit$frame[["autocor"]]
  if (!length(acs)) return(NULL)
  sdr <- sdr_of(fit)
  Vfull <- sdr$cov.fixed
  pos <- which(rownames(Vfull) == "thetaac")
  th <- fit$estimates[["thetaac"]]
  rows <- list()
  for (ac in acs) {
    t0 <- th[ac[["theta_idx"]]]
    Vth <- Vfull[pos[ac[["theta_idx"]]], pos[ac[["theta_idx"]]], drop = FALSE]
    types <- autocor_types(ac)
    gfun <- function(tt) {
      v <- unname(autocor_natural(tt, ac))
      ifelse(types == "cor",
             atanh(pmin(pmax(v, -varcorr_cor_clamp), varcorr_cor_clamp)),
             v)
    }
    est0 <- unname(autocor_natural(t0, ac))
    bd <- types == "cor" & abs(est0) >= varcorr_cor_clamp
    est_t <- ifelse(types == "cor", atanh(est0), est0)
    g0 <- gfun(t0)
    J <- vapply(seq_along(t0), function(i) {
      h <- 1e-5 * max(abs(t0[i]), 1)
      tp <- t0; tp[i] <- tp[i] + h
      tm <- t0; tm[i] <- tm[i] - h
      (gfun(tp) - gfun(tm)) / (2 * h)
    }, numeric(length(g0)))
    J <- matrix(J, nrow = length(g0))
    se_t <- sqrt(pmax(diag(J %*% Vth %*% t(J)), 0))
    se_t[bd | !is.finite(est_t) | !is.finite(se_t)] <- NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      block = ac[["block_label"]], term = names(autocor_natural(t0, ac)),
      type = types, est_t = est_t, se_t = se_t,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Estimated within-group residual correlation matrix
#'
#' The correlation matrix `R` of the `ar()`, `ma()`, `arma()`, `cosy()`
#' or `unstr()` term of a fit, over the model's time levels. The
#' residual covariance of a group is `D R D` restricted to the time
#' points that group has, with `D` the diagonal matrix of that group's
#' `sigma` values, so `R` plus [sigma()] (or `predict(dpar = "sigma")`)
#' describes the whole residual structure.
#'
#' @param fit A `frmtmb_fit`.
#' @param resp Response name, for a multivariate fit.
#' @return A correlation matrix with the time levels as dimnames, or
#'   `NULL` when the fit has no residual correlation term.
#' @seealso [rescor_matrix()] for the ACROSS-response residual
#'   correlation of a `rescor = TRUE` fit, which is a different
#'   structure and cannot be combined with this one.
#' @examples
#' set.seed(1)
#' d <- expand.grid(week = 1:5, subj = factor(1:25))
#' e <- as.vector(apply(matrix(rnorm(125), 5, 25), 2, function(z) {
#'   as.vector(stats::filter(z, 0.6, "recursive"))
#' }))
#' d$x <- rnorm(125)
#' d$y <- 1 + 0.5 * d$x + e
#' fit <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
#'            data = d)
#' autocor_matrix(fit)
#' @export
autocor_matrix <- function(fit, resp = NULL) {
  acs <- fit$frame[["autocor"]]
  if (!length(acs)) return(NULL)
  ac <- if (is.null(resp)) acs[[1L]] else acs[[resp]]
  if (is.null(ac)) {
    stop("autocor_matrix(): no residual correlation term for response '",
         resp, "'", call. = FALSE)
  }
  R <- autocor_cor(fit$estimates[["thetaac"]][ac[["theta_idx"]]], ac)
  R <- as.matrix(R)
  dimnames(R) <- list(ac[["time_levels"]], ac[["time_levels"]])
  R
}

#' Per-group residual draws for `simulate()`: the group's correlation
#' submatrix, Cholesky-factorized off the tape, applied to standard
#' normal (or scaled-t) innovations and scattered back to the rows.
#'
#' @noRd
autocor_draw_resid <- function(ac, R, sigma, n, nu = NULL) {
  e <- numeric(n)
  for (pt in ac[["patterns"]]) {
    k <- pt$k
    Rp <- if (k == ac[["d"]]) R else {
      matrix(as.vector(R)[pt$gather], k, k)
    }
    L <- t(chol(Rp))
    Z <- L %*% matrix(stats::rnorm(k * pt$G), k, pt$G)
    if (!is.null(nu)) {
      # multivariate t: one chi-square scaling per GROUP, not per row
      Z <- Z * rep(sqrt(nu / stats::rchisq(pt$G, nu)), each = k)
    }
    e[pt$rows] <- as.vector(Z)
  }
  e * sigma
}
