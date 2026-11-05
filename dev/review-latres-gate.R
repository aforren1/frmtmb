# Short check 2, item 3: the exact change, tried in-process on the final build.
# latres_shift() is the body proposed for R/autocor-latent.R, and
# latres_laplace_reason_new() replaces latres_laplace_reason() up to the row limit.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
ns <- asNamespace("frmtmb")
old_reason <- get("latres_laplace_reason", ns)

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
  l0 <- ll(eta); lp <- ll(eta + h); ln <- ll(eta - h); lp2 <- ll(eta + 2 * h); ln2 <- ll(eta - 2 * h)
  h3 <- -(lp2 - 2 * lp + 2 * ln - ln2) / (2 * h^3)
  h4 <- -(lp2 - 4 * lp + 6 * l0 - 4 * ln + ln2) / h^4
  pj <- latres_jacobian(bk, at$estimates[["theta"]][bk[["theta_idx"]]])
  A <- Matrix::sparseMatrix(i = pj$i, j = pj$j, x = pj$x, dims = c(a[["n"]], length(env$random)))
  Ch <- Matrix::Cholesky(env$spHess(full, random = TRUE), LDL = FALSE)
  v <- as.matrix(A %*% Matrix::solve(Ch, Matrix::t(A)))
  dv <- diag(v)
  -sum(h4 * dv^2) / 8 + sum(outer(h3 * dv, h3 * dv) * v) / 8 + sum(outer(h3, h3) * v^3) / 12
}
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
    pp <- par; pp[k] <- pp[k] + 0.05
    pm <- par; pm[k] <- pm[k] - 0.05
    g[k] <- (latres_tk(fit, pp, bk, map) - latres_tk(fit, pm, bk, map)) / 0.1
  }
  g[!is.finite(g)] <- 0
  th <- which(names(par) == "theta")[bk[["theta_idx"]][1:2]]
  if (anyNA(th) || !all(is.finite(V[th, ]))) return(Inf)
  d <- as.numeric(V[th, ] %*% g)
  t2 <- par[th[2]]
  max(abs(d[1]), abs((t2 + d[2]) / sqrt(1 + (t2 + d[2])^2) - t2 / sqrt(1 + t2^2)))
}
new_reason <- function(fit, bk) {
  r <- bk[["aux"]][["ac"]][["resp"]]
  fam <- fit$spec$responses[[r]]$family
  if (fit$frame[["n_obs"]] > getOption("frmtmb.latres_se_rows", 3000L)) return(old_reason(fit, bk))
  if (identical(fam[["type"]], "ordinal")) return("each row is one category")
  s <- tryCatch(latres_shift_tk(fit, bk), error = function(e) Inf)
  if (s <= 0.15) return(NULL)
  if (!is.finite(s)) return("the fit has no curvature to read its first-order correction by")
  paste0("the first-order correction moves log(sderr) or ar[1] by ", format(s, digits = 2), ", above 0.15")
}
for (nm in c("latres_tk", "latres_shift_tk", "new_reason")) {
  f <- get(nm)
  environment(f) <- ns
  assign(nm, f)
}
environment(latres_tk) <- ns
sh_env <- new.env(parent = ns)
sh_env$latres_tk <- latres_tk
environment(latres_shift_tk) <- sh_env
rn_env <- new.env(parent = ns)
rn_env$latres_shift_tk <- latres_shift_tk
rn_env$old_reason <- old_reason
environment(new_reason) <- rn_env
assignInNamespace("latres_laplace_reason", new_reason, ns = "frmtmb")
if (nzchar(Sys.getenv("REV_GATE_TEST"))) {
  suppressPackageStartupMessages(library(testthat))
  Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
  res <- test_file("tests/testthat/test-latres.R", package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  cat("tests", sum(df$nb), "fail", sum(df$failed), "err", sum(df$error), "
")
  print(unique(df$test[df$failed > 0 | df$error]))
  for (r in res) for (e in r$results) if (inherits(e, "expectation_failure") || inherits(e, "expectation_error")) cat(substr(r$test, 1, 50), "::", substr(gsub("
", " ", conditionMessage(e)), 1, 150), "
")
  quit("no")
}
run <- function(lab, sim, form, fam, seeds = 1:6, G = 20, K = 10) {
  w <- 0L; msg <- ""; ts <- 0
  for (s in seeds) {
    set.seed(s)
    d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
    d <- d[order(d$g, d$t), ]
    n <- nrow(d)
    d$x <- rnorm(n)
    d$e <- as.vector(replicate(G, arima.sim(list(ar = .6), K, sd = 1)))
    d$tr <- 10L
    d$y <- sim(d)
    wr <- character(0)
    t0 <- proc.time()[[3]]
    withCallingHandlers(frm(form, data = d, family = fam), warning = function(c) {
      if (inherits(c, "frmtmb_latres_laplace")) wr <<- c(wr, conditionMessage(c))
      invokeRestart("muffleWarning")
    }, message = function(c) invokeRestart("muffleMessage"))
    ts <- ts + proc.time()[[3]] - t0
    if (length(wr)) { w <- w + 1L; msg <- wr[1] }
  }
  cat(sprintf("%-20s warned %d of %d, %.2f s a fit: %s\n", lab, w, length(seeds), ts / length(seeds), substr(msg, 1, 150)))
}
run("bernoulli", function(d) rbinom(nrow(d), 1, plogis(.2 + .5 * d$x + d$e)), y ~ x + ar(t, g, cov = TRUE), bernoulli())
run("poisson mean 20", function(d) rpois(nrow(d), exp(3 + .3 * d$x + .5 * d$e)), y ~ x + ar(t, g, cov = TRUE), poisson())
run("cumulative", function(d) factor(findInterval(.5 * d$x + d$e + rlogis(nrow(d)), c(-1, 0, 1)) + 1L, ordered = TRUE), y ~ x + ar(t, g, cov = TRUE), cumulative())
run("lognormal, sigma 2", function(d) rlnorm(nrow(d), .3 * d$x + d$e, 2), y ~ x + ar(t, g, cov = TRUE), lognormal())
run("skew_normal, units 10", function(d) 10 * (.3 * d$x + d$e + rnorm(nrow(d), 0, .5)), y ~ x + ar(t, g, cov = TRUE), skew_normal())
run("Beta, phi 4", function(d) { m <- plogis(.5 * d$e); pmin(pmax(rbeta(nrow(d), 4 * m, 4 * (1 - m)), 1e-5), 1 - 1e-5) }, y ~ x + ar(t, g, cov = TRUE), Beta())
