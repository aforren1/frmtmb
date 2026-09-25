# gr(g, by = f) on draws: ranef(), coef() and VarCorr() are brms's own
# output on the same draws, compared with identical() through the shim
# test-brms-output.R builds (brms's brm(empty = TRUE) object given the
# sampler's draws under the names frmtmb.sample puts on them).
#
# frmtmb fits a by-split term as one block per by-level, and VarCorr()
# names its coefficients per by-level (Intercept:fa). brms's ranef() has
# one column per coefficient over EVERY level of g, which is the layout
# draws_ranef_layout() has to merge back.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

grby_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    skip_if_not_installed("brms")
    skip_if_not_installed("posterior")
    if (is.null(cache)) {
      set.seed(11)
      dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
      dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
      u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
      dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] +
                             u[dd$g, 2] * dd$x, 1)
      fit <- frm(bf(y ~ x + (1 + x | gr(g, by = f))), family = gaussian(),
                 data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 400, refresh = 0, seed = 3)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

grby_shim <- function(cs, M) {
  b <- suppressWarnings(suppressMessages(
    brms::brm(y ~ x + (1 + x | gr(g, by = f)), data = cs$dd, empty = TRUE,
              backend = "rstan")))
  sf <- cs$ds$stanfit
  nc <- sf@sim$chains
  n <- nrow(M) %/% nc
  sim <- sf@sim
  sim$samples <- lapply(seq_len(nc), function(ch) {
    rows <- (ch - 1L) * n + seq_len(n)
    stats::setNames(lapply(seq_len(ncol(M)), function(j) unname(M[rows, j])),
                    colnames(M))
  })
  sim$iter <- n
  sim$warmup <- 0
  sim$warmup2 <- rep(0, nc)
  sim$n_save <- rep(n, nc)
  sim$permutation <- lapply(seq_len(nc), function(i) seq_len(n))
  sim$pars_oi <- sim$fnames_oi <- colnames(M)
  sim$dims_oi <- stats::setNames(rep(list(integer(0)), ncol(M)),
                                 colnames(M))
  sim$n_flatnames <- ncol(M)
  sf@sim <- sim
  b$fit <- sf
  b
}

# brms stores r_ draws every level of one coefficient before the next;
# frmtmb stores them level by level
grby_brms_order <- function(M) {
  r <- grep("^r_", colnames(M))
  nm <- colnames(M)[r]
  lev <- as.integer(sub("^r_g[[]([0-9]+),.*$", "\\1", nm))
  cf <- sub("^r_g[[][0-9]+,(.*)[]]$", "\\1", nm)
  o <- order(match(cf, unique(cf)), lev)
  M[, r] <- M[, r[o]]
  colnames(M)[r] <- nm[o]
  M
}

test_that("ranef() and coef() of a by-split term are brms's", {
  cs <- grby_case()
  ds <- cs$ds
  # every level of g, one column per coefficient
  dn <- dimnames(ranef(ds)$g)
  expect_identical(as.vector(dn[[1L]]), as.character(1:16))
  expect_identical(dn[[3L]], c("Intercept", "x"))
  shb <- grby_shim(cs, grby_brms_order(ds$draws))
  expect_identical(ranef(ds), brms:::ranef.brmsfit(shb))
  expect_identical(ranef(ds, FALSE), brms:::ranef.brmsfit(shb, FALSE))
  expect_identical(coef(ds), brms:::coef.brmsfit(shb))
})

test_that("VarCorr() of a by-split term is brms's", {
  cs <- grby_case()
  ds <- cs$ds
  vd <- VarCorr(ds, summary = FALSE)
  expect_identical(colnames(vd$g$sd),
                   c("Intercept:fa", "x:fa", "Intercept:fb", "x:fb"))
  M <- cbind(ds$draws,
             "sd_g__Intercept:fa" = unname(vd$g$sd[, "Intercept:fa"]),
             "sd_g__x:fa" = unname(vd$g$sd[, "x:fa"]),
             "sd_g__Intercept:fb" = unname(vd$g$sd[, "Intercept:fb"]),
             "sd_g__x:fb" = unname(vd$g$sd[, "x:fb"]),
             "cor_g__Intercept:fa__x:fa" = vd$g$cor[, "Intercept:fa", "x:fa"],
             "cor_g__Intercept:fb__x:fb" = vd$g$cor[, "Intercept:fb", "x:fb"])
  sv <- grby_shim(cs, M)
  expect_identical(VarCorr(ds), brms:::VarCorr.brmsfit(sv))
})
