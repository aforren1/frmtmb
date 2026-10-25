p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.sample/tests/testthat/test-autocor-seam-draws.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("# The draws of a model with a residual correlation term of core's
# registry carry brms's name for its parameter. frmtmb.spatial's sar()
# reports lagsar or errorsar; a term registered here, through the
# exported frmtmb_register_autocor() alone, reports lagsar the same way,
# so the test needs no other extension.",
"# The draws of a model with a residual correlation term of core's
# registry carry the name the term gives its parameter, read off the
# registry (frmtmb_autocor_pars()), and a bounded parameter gets brms's
# flat default prior on its interval. frmtmb.spatial's sar() reports
# lagsar or errorsar; the terms registered here, through the exported
# frmtmb_register_autocor() alone, report lagsar or a name frmtmb.sample
# has never seen, so the test needs no other extension.")
rep1("reg_lagterm <- function() {
  frmtmb::frmtmb_register_autocor(",
"reg_lagterm <- function(par = \"lagsar\") {
  frmtmb::frmtmb_register_autocor(")
rep1("    natural = function(theta, aux) {
      c(lagsar = aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1])))
    },",
"    natural = function(theta, aux) {
      stats::setNames(aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1])),
                      par)
    },
    bounds = function(aux) cbind(lower = aux$lo, upper = aux$hi),")
rep1("      r <- y - rho * as.vector(aux$W %*% y) - mu
      n <- length(y)
      lsig <- if (length(sigma) == 1L) n * log(sigma) else sum(log(sigma))
      list(u = r / sigma, logjac = sum(log1p(-rho * aux$ev)) - lsig)",
"      r <- y - rho * as.vector(aux$W %*% y) - mu
      # sigma_varies = FALSE: sigma is one value
      list(u = r / sigma,
           logjac = sum(log1p(-rho * aux$ev)) - length(y) * log(sigma))")
rep1("      as.vector(solve(diag(nrow(aux$W)) - rho * aux$W, mu + sigma * u))",
"      solve(diag(nrow(aux$W)) - rho * aux$W, mu + sigma * u)")
rep1("  # the absent case: a model without the term has no such column
  f0 <- frm(y ~ x, data = cs$d)
  expect_false(\"lagsar\" %in% frmtmb.sample:::draws_cov_cols(f0)$names)
})",
"  # the absent case: a model without the term has no such column
  f0 <- frm(y ~ x, data = cs$d)
  expect_false(\"lagsar\" %in% frmtmb.sample:::draws_cov_cols(f0)$names)
})

test_that(\"a parameter that is not called lagsar is a column under its own
           name\", {
  frmtmb::frmtmb_local_registries()
  reg_lagterm(\"toyphi\")
  cs <- lag_case()
  fit <- frm(y ~ x + lagtest(W), data = cs$d, data2 = list(W = cs$W))
  dv <- frmtmb.sample:::draws_cov_cols(fit)
  expect_identical(dv$names, \"toyphi\")
  expect_true(\"thetaac\" %in% names(dv$cols))
  # and its default prior is the flat one on its interval, by its name
  dp <- frmtmb.sample:::default_priors_for(fit)
  cls <- vapply(unclass(dp), function(s) s$class %||% \"\", \"\")
  expect_true(\"toyphi\" %in% cls)
  expect_false(\"lagsar\" %in% cls)
})")
txt <- paste0(txt, "

test_that(\"the default prior on a bounded parameter is flat on its
           interval, as brms's on rho\", {
  skip_sampler()
  frmtmb::frmtmb_local_registries()
  # a term whose likelihood does not depend on its parameter: the
  # draws of the parameter are then its prior. brms's prior on lagsar
  # and errorsar is flat on the interval; a prior flat on the internal
  # parameter instead puts the draws at the ends.
  lo <- -0.4
  hi <- 0.9
  frmtmb::frmtmb_register_autocor(
    \"flatphi\", args = \"V\",
    build = function(spec, rows, n_all, lookup) list(n = length(rows)),
    npar = 1L, start = function(aux) 0,
    natural = function(theta, aux) {
      c(toyphi = lo + (hi - lo) / (1 + exp(-theta[1])))
    },
    bounds = function(aux) cbind(lower = lo, upper = hi),
    whiten = function(y, mu, sigma, theta, aux) {
      list(u = (y - mu) / sigma, logjac = -length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) mu + sigma * u,
    sigma_varies = FALSE)
  set.seed(3)
  d <- data.frame(x = stats::rnorm(30))
  d$y <- 1 + d$x + stats::rnorm(30)
  fit <- allow_warnings(
    frm(y ~ x + flatphi(V), data = d, data2 = list(V = 1)),
    \"Standard errors are not available\")
  msgs <- character()
  ds <- withCallingHandlers(
    frm_sample(fit, chains = 4, iter = 3000, warmup = 500, refresh = 0,
               seed = 11),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart(\"muffleMessage\")
    })
  # the announcement names the parameter and its prior
  expect_true(any(grepl(\"toyphi +uniform\\(-0.4, 0.9\\)  \\[flat on its \",
                        msgs)))
  z <- as.matrix(ds)[, \"toyphi\"]
  expect_true(all(z > lo & z < hi))
  # each decile of the interval holds a tenth of the draws, within four
  # standard errors at the draws' effective size
  ess <- posterior::ess_basic(z)
  se <- sqrt(0.1 * 0.9 / ess)
  frac <- tabulate(findInterval(z, lo + (hi - lo) * (1:9) / 10) + 1L, 10L) /
    length(z)
  expect_lt(max(abs(frac - 0.1)) / se, 4)
  # and the mean is the interval's middle, within four standard errors
  expect_lt(abs(mean(z) - (lo + hi) / 2) / (sd(z) / sqrt(ess)), 4)
})
")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
