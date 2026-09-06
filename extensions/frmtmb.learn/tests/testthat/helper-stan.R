# The identity tier: every family checked against an independent Stan
# program of the same model, at the same point.
#
# WHY AN IDENTITY AND NOT AN AGREEMENT. A recursion vectorized across
# subjects and masked for unequal trial counts is exactly the kind of
# code that can be wrong in a way that still converges and still
# recovers roughly the right parameters. A second implementation, in
# another language and another automatic-differentiation system, written
# from the published equations rather than from this package's sources,
# either agrees to machine precision or does not.
#
# WHAT MAKES THE DIFFERENCE EXACTLY ZERO. Each program declares exactly
# the parameters frmtmb estimates, on the scales frmtmb estimates them
# on, so the map is the identity and the comparison carries no Jacobian.
# The random-effect standard deviation rides along as DATA, which keeps
# log_prob() free of a constrained parameter. frm() carries no default
# prior, the programs use `target +=` with `_lpmf` and `_lpdf` so they
# drop no normalizing constant, and every declared parameter is
# unbounded so `adjust_transform = FALSE` has no Jacobian to discard.
#
# None of hBayesDM's code is read or copied. The equations are the
# published models and the programs below are written from them.

# Compiling one program costs about a minute and the program is a
# function of the model only, never of the data, so the cache is content
# addressed on the Stan code. The rstan version joins the key because a
# DSO built by one rstan is not loadable by another.
ln_stan_cache_dir <- function() {
  dir <- Sys.getenv("FRMTMB_STAN_CACHE", "")
  if (!nzchar(dir)) {
    dir <- testthat::test_path("..", "..", "..", "..", "dev", "stan-cache")
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  normalizePath(dir, winslash = "/", mustWork = TRUE)
}

# md5 over a file rather than over a string keeps this on base R: a test
# suite must not gain a dependency to hash one character vector.
ln_stan_cache_key <- function(code) {
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))), f)
  unname(tools::md5sum(f))
}

# Programs already in hand in THIS session. Not an optimization: a
# stanfit model object compiled in this session and then re-read from
# its own RDS comes back with a DSO that will not initialize, so the
# first call for a program works and every later call in the same
# session fails. Reading an RDS written by an EARLIER session is fine.
.ln_stan_models <- new.env(parent = emptyenv())

ln_stan_model <- function(code) {
  key <- ln_stan_cache_key(code)
  hit <- .ln_stan_models[[key]]
  if (!is.null(hit)) return(hit)
  path <- file.path(ln_stan_cache_dir(), paste0(key, ".rds"))
  if (file.exists(path)) {
    mod <- try(readRDS(path), silent = TRUE)
    if (!inherits(mod, "try-error")) {
      .ln_stan_models[[key]] <- mod
      return(mod)
    }
  }
  mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(mod, path)
  .ln_stan_models[[key]] <- mod
  mod
}

# Gated the way the brms comparison tier is: outside R CMD check BOTH
# variables are needed to opt in.
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
skip_unless_stan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
    testthat::skip("set FRMTMB_BRMS_FIT_TESTS=true to run the Stan tier")
  }
}

# The block, rebuilt from the data rather than read off the fit, so that
# a mistake in this package's own frame assembly shows up as a mismatch
# instead of cancelling on both sides.
ln_stan_block <- function(data, id = "id", trial = "trial") {
  gv <- factor(data[[id]])
  tv <- data[[trial]]
  rows <- lapply(split(seq_len(nrow(data)), gv), function(r) r[order(tv[r])])
  len <- lengths(rows)
  nt <- max(len)
  idx <- t(matrix(vapply(rows, function(r) c(r, rep(r[1L], nt - length(r))),
                         integer(nt)), nrow = nt))
  mask <- t(matrix(vapply(len, function(l) as.numeric(seq_len(nt) <= l),
                          numeric(nt)), nrow = nt))
  list(S = length(rows), T = nt, idx = idx, mask = mask,
       subj = as.integer(gv))
}

# The data every program in this file takes, plus whatever the family
# adds. `X` is rebuilt with model.matrix() rather than read off the
# frame, for the same reason the block is.
ln_stan_data <- function(fit, data, form, extra = list()) {
  blk <- ln_stan_block(data)
  X <- unname(stats::model.matrix(form, data))
  sd_u <- sqrt(unname(frmtmb::VarCorr(fit)[[1L]])[1L, 1L])
  c(list(N = nrow(data), S = blk$S, T = blk$T, K = ncol(X),
         idx = blk$idx, mask = blk$mask, subj = blk$subj,
         choice = as.integer(data$choice), X = X, sd_u = sd_u),
    extra)
}

# frmtmb's full parameter vector, split the way every program below
# declares it: `beta` is the PRIMARY distributional parameter's fixed
# effects, `betad` is every other parameter's in declaration order, and
# `b` is the subject deviations themselves rather than a standardized
# version of them.
#
# as.array() on all three is load-bearing for the one-parameter case. A
# length-one R numeric reaches rstan with no `dim`, and a program that
# declares `vector[1] bd` refuses it outright: "dims declared=(1); dims
# found=()". bandit2arm_delta() is the only family here with a single
# non-primary parameter, so it was the only one to meet this, and it
# ERRORED rather than disagreeing, which is why its row was missing from
# the identity table instead of showing a bad residual.
ln_stan_pars <- function(par) {
  list(b = as.array(unname(par[names(par) == "beta"])),
       bd = as.array(unname(par[names(par) == "betad"])),
       u = as.array(unname(par[names(par) == "b"])))
}

# The two checks, at one point.
#
#   A  Stan's log_prob at frmtmb's estimates equals frmtmb's joint log
#      density there, up to a constant this asserts is zero.
#   B  the gradient with respect to the subject effects vanishes, which
#      is what catches a wrong map: frmtmb puts them at their
#      conditional modes.
#   C  the gradient with respect to the POPULATION parameters equals
#      frmtmb's own joint gradient there. It is NOT asserted against
#      zero, and could not be: maximum likelihood makes the MARGINAL
#      likelihood stationary, not the joint one, so this number is far
#      from zero. Comparing it against frmtmb's own gradient instead is
#      an identity of the same kind as the value check, at no extra
#      cost, and it catches a wrong fixed-effect map that a
#      single-point value identity can miss.
#
#      obj$env$f(par, order = 1) is the gradient of the joint NEGATIVE
#      log density, aligned with par, so Stan's gradient of the log
#      density is its negation. The Stan positions are in declaration
#      order, b then bd then u, so the complement of the u block lines
#      up with beta then betad without further bookkeeping.
#
# Setting FRMTMB_LEARN_LP_LOG to a file path makes each check append its
# numbers there. The identity table in dev/learn-findings.md is that
# file rather than a transcription of it, so the published residuals
# cannot drift from the ones the suite asserts.
ln_lp_check <- function(fit, code, sdat, par = NULL, tol = 1e-6,
                        tol_grad = 1e-4, check_grad = TRUE, label = NA) {
  mod <- ln_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  par <- if (is.null(par)) fit$obj$env$last.par.best else par
  pars <- ln_stan_pars(par)
  up <- rstan::unconstrain_pars(sf, pars)
  lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
  ours <- -fit$obj$env$f(par)
  testthat::expect_lt(abs(lp - ours), tol * max(1, abs(ours)))
  gin <- NA_real_
  gall <- NA_real_
  gout <- NA_real_
  if (check_grad) {
    gr <- rstan::grad_log_prob(sf, up, adjust_transform = FALSE)
    bump <- pars
    bump$u <- bump$u + 1
    inner <- which(rstan::unconstrain_pars(sf, bump) != up)
    testthat::expect_gt(length(inner), 0)
    gin <- max(abs(gr[inner]))
    gall <- max(abs(gr))
    testthat::expect_lt(gin, tol_grad)
    # check C: the population block against frmtmb's own joint gradient
    pop <- setdiff(seq_along(up), inner)
    gj <- fit$obj$env$f(par, order = 1)
    ours_g <- -gj[names(par) %in% c("beta", "betad")]
    testthat::expect_length(pop, length(ours_g))
    gout <- max(abs(gr[pop] - ours_g))
    testthat::expect_lt(gout, 1e-6 * max(1, max(abs(ours_g))))
  }
  logf <- Sys.getenv("FRMTMB_LEARN_LP_LOG", "")
  if (nzchar(logf)) {
    cat(sprintf("%s|%.10f|%.10f|%.3e|%.3e|%.3e|%.3e|%d|%d\n", label,
                lp, ours, lp - ours, gin, gall, gout, sdat$S, sdat$N),
        file = logf, append = TRUE)
  }
  invisible(list(lp = lp, ours = ours, const = lp - ours, g_inner = gin,
                 g_outer = gout))
}
