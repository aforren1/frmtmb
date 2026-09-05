# The agreement tier for the drift-diffusion families: each family's
# simulator against that family's OWN log-density, the same contract
# frmtmb's tests/testthat/test-simulate-density.R holds the built-in
# families to.
#
# WHY it is worth its own file here. The extension's other tests check
# the density against RWiener and the exported ddm_simulate() against
# the generative process, but the seam frmtmb actually uses is neither
# of those: it is the family's `sim` SLOT, reached by frm_simulate()
# and by simulate() on a fit. Nothing tied that slot to the density
# beside it.
#
# The shape of the check is different from a scalar family's, because a
# choice-RT density is DEFECTIVE: `lpdf` at a row gives the joint
# density of (time, boundary), so it integrates to the probability of
# that boundary rather than to one. The `sim` slot draws a time
# CONDITIONAL on the boundary the row already observed. So the
# reference is the defective density renormalized on its own boundary,
# and the two boundaries' masses are asserted to sum to one, which is
# itself the check that no probability has gone missing.

skip_on_cran()

# The tolerance argument is the one made in the core tier: the seed is
# fixed, so a false alarm would be permanent rather than flaky, and at
# this alpha the chi-square critical value is around 40 while a
# simulator drawing from the wrong distribution lands in the hundreds.
DDM_ALPHA <- 1e-6
DDM_N <- 4000L
DDM_CELLS <- 8L
DDM_SEED <- 20250905L

#' Cell probabilities from a family's own log-density on a bracket,
#' discretized with trapezoid weights. Returns the cut points, each
#' cell's probability CONDITIONAL on the bracket, and the total
#' (defective) mass the bracket carries.
#' Per-row family DATA the likelihood needs but no addition term
#' supplies (gddm's grid index arithmetic). Assembly builds it once from
#' the response; a density evaluated on a grid of candidate responses
#' has to build it the same way, or it is scoring the wrong rows.
ddm_at <- function(fam, y, at) {
  at <- lapply(at, function(v) rep(v, length.out = length(y)))
  if (!is.null(fam[["aterm_data"]])) {
    at <- c(at, fam[["aterm_data"]](y, at))
  }
  at
}

ddm_cells <- function(fam, dp, at, lo, hi, B = DDM_CELLS, m = 8001L) {
  y <- seq(lo, hi, length.out = m)
  d <- exp(as.numeric(fam[["lpdf"]](
    y, lapply(dp, function(v) rep(v, length(y))), ddm_at(fam, y, at))))
  h <- diff(y)
  w <- c(h[1L], h[-1L] + h[-length(h)], h[length(h)]) / 2
  p <- d * w
  mass <- sum(p)
  p <- p / mass
  cp <- cumsum(p)
  idx <- unique(pmin(findInterval(seq_len(B - 1L) / B, cp) + 1L, m - 1L))
  cuts <- unique(y[idx])
  cell <- findInterval(y, cuts) + 1L
  list(cuts = cuts, mass = mass,
       p = as.numeric(tapply(p, factor(cell,
                                       levels = seq_len(length(cuts) + 1L)),
                             sum, default = 0)))
}

expect_ddm_gof <- function(draws, cells, label) {
  n <- length(draws)
  obs <- tabulate(findInterval(draws, cells[["cuts"]]) + 1L,
                  length(cells[["p"]]))
  ex <- cells[["p"]] * n
  ok <- !is.na(ex) & ex >= 5
  obs <- obs[ok]
  ex <- ex[ok]
  df <- length(obs) - 1L
  expect_gte(df, 1L, label = paste0(label, ": usable cells"))
  expect_lt(sum((obs - ex)^2 / ex), stats::qchisq(1 - DDM_ALPHA, df),
            label = paste0(label, " chi-square (df ", df, ")"))
}

# ---------------------------------------------------------------- wiener

test_that("wiener: the sim slot draws from the density's own conditional", {
  skip_if_not_installed("RWiener")
  fam <- wiener()
  dp <- list(mu = 1.1, bs = 1.5, ndt = 0.2, bias = 0.45)
  # the family derives its ndt link from the response, so the density
  # under test is the finalized one, the family the likelihood scores
  fin <- fam[["family_finalize"]](fam, c(0.3, 1.0), list(dec = c(0, 1)))
  lo <- dp[["ndt"]] + 1e-6
  hi <- dp[["ndt"]] + 30
  mass <- 0
  for (b in c(0, 1)) {
    at <- list(dec = b)
    cells <- ddm_cells(fin, dp, at, lo, hi)
    mass <- mass + cells[["mass"]]
    set.seed(DDM_SEED + b)
    draws <- fam[["sim"]](lapply(dp, function(v) rep(v, DDM_N)),
                          list(dec = rep(b, DDM_N)), DDM_N)
    lab <- paste0("wiener[boundary ", b, "]")
    expect_true(all(draws > dp[["ndt"]]),
                label = paste0(lab, " above the non-decision time"))
    expect_ddm_gof(draws, cells, lab)
  }
  # nothing has gone missing between the two boundaries
  expect_equal(mass, 1, tolerance = 1e-4)
})

test_that("wiener: frm_simulate() reaches the family seam", {
  skip_if_not_installed("RWiener")
  set.seed(4)
  n <- 2L * DDM_N
  d0 <- ddm_simulate(n, mu = 1.1, bs = 1.5, ndt = 0.2, bias = 0.45)
  dd <- data.frame(rt = d0$rt, up = as.integer(d0$upper))
  dp <- list(mu = 1.1, bs = 1.5, ndt = 0.2, bias = 0.45)
  s <- frmtmb::frm_simulate(rt | dec(up) ~ 1, dd, family = wiener(),
                            newparams = list(Intercept = dp[["mu"]],
                                             bs = dp[["bs"]],
                                             ndt = dp[["ndt"]],
                                             bias = dp[["bias"]]),
                            nsim = 1L, seed = DDM_SEED)
  v <- s[["sim_1"]]
  expect_equal(length(v), n)
  expect_true(all(v > dp[["ndt"]]))
  # each row's draw is conditional on THAT row's boundary, so the two
  # boundaries are tested apart
  fin <- wiener()[["family_finalize"]](wiener(), dd$rt, list(dec = dd$up))
  for (b in c(0L, 1L)) {
    idx <- which(dd$up == b)
    cells <- ddm_cells(fin, dp, list(dec = b), dp[["ndt"]] + 1e-6,
                       dp[["ndt"]] + 30)
    expect_ddm_gof(v[idx], cells, paste0("wiener/frm_simulate[", b, "]"))
  }
})

test_that("wiener: simulate() on a fit reaches the same seam", {
  skip_if_not_installed("RWiener")
  set.seed(5)
  n <- 400L
  d0 <- ddm_simulate(n, mu = 1.2, bs = 1.4, ndt = 0.25, bias = 0.5)
  dd <- data.frame(rt = d0$rt, up = as.integer(d0$upper))
  fit <- frmtmb::frm(rt | dec(up) ~ 1, family = wiener(), data = dd)
  sm <- stats::simulate(fit, nsim = 3L, seed = DDM_SEED)
  expect_equal(dim(sm), c(n, 3L))
  expect_false(anyNA(sm))
  # a redrawn time still clears the fitted non-decision time, which is
  # the one bound a wiener draw can never violate
  expect_true(all(as.matrix(sm) > 0))
})

# ------------------------------------------------------------------ gddm

test_that("gddm: the sim slot draws from the density's own conditional", {
  set.seed(6)
  g <- gddm_simulate(600L, coh = rep(c(0.2, 0.5), 300L))
  av <- list(dec = as.numeric(g$upper), vint1 = as.numeric(g$cond))
  if (!is.null(g$coh)) av$vreal1 <- g$coh
  # the modeled time window is fixed at finalize time and the density
  # is not defined outside it, so it is named rather than inherited
  # from whatever the largest simulated response happened to be
  fam <- gddm(control = gddm_control(t_max = 6))
  # gddm installs its density AND its simulator in family_finalize, so
  # neither exists on the family as written
  expect_null(fam[["sim"]])
  fin <- fam[["family_finalize"]](fam, g$rt, av)
  expect_false(is.null(fin[["sim"]]))

  dp <- list(mu = 1.0, bs = 1.0, bias = 0.5, ndt = 0.2)
  dp <- dp[fin[["dpars"]]]
  lo <- dp[["ndt"]] + 1e-6
  hi <- fin[["gddm"]][["ctl"]][["t_max"]] - 1e-6
  mass <- 0
  cond <- 1
  for (b in c(0, 1)) {
    at <- list(dec = b, vint1 = cond)
    cells <- ddm_cells(fin, dp, at, lo, hi)
    mass <- mass + cells[["mass"]]
    set.seed(DDM_SEED + b)
    ydum <- rep(dp[["ndt"]] + 0.5, DDM_N)
    draws <- fin[["sim"]](lapply(dp, function(v) rep(v, DDM_N)),
                          ddm_at(fin, ydum, at), DDM_N)
    lab <- paste0("gddm[boundary ", b, "]")
    expect_false(anyNA(draws), label = paste0(lab, " no NA draws"))
    expect_ddm_gof(draws, cells, lab)
  }
  expect_equal(mass, 1, tolerance = 1e-2)
})

test_that("gddm: frm_simulate() reaches the family seam", {
  # This is the regression. gddm's simulator is installed by
  # family_finalize, and frm_simulate() used to read the family as
  # WRITTEN, before assembly had finalized it, and refuse the model
  # for having no simulator, while simulate() on a fit of the same
  # model worked.
  set.seed(7)
  g <- gddm_simulate(500L, coh = rep(c(0.2, 0.5), 250L))
  gd <- data.frame(rt = g$rt, up = as.integer(g$upper),
                   cond = as.integer(g$cond))
  if (!is.null(g$coh)) gd$coh <- g$coh
  s <- frmtmb::frm_simulate(rt | dec(up) + vint(cond) ~ 1, gd,
                            family = gddm(),
                            newparams = list(Intercept = 1.0, bs = 1.0,
                                             bias = 0.5, ndt = 0.2),
                            nsim = 1L, seed = DDM_SEED)
  expect_equal(nrow(s), nrow(gd))
  expect_false(anyNA(s[["sim_1"]]))
  expect_true(all(s[["sim_1"]] > 0.2))
})

# ------------------------------------------------------------------- lba

test_that("lba: the sim slot draws from the density's own conditional", {
  n_acc <- 2L
  fam <- lba(n_acc)
  dp <- list(v1 = 1.6, v2 = 1.0, A = 0.5, k = 0.4, ndt = 0.2)
  fin <- fam[["family_finalize"]](fam, c(0.5, 0.9), list(vint1 = c(1, 2)))
  lo <- dp[["ndt"]] + 1e-6
  hi <- dp[["ndt"]] + 12
  mass <- 0
  for (w in seq_len(n_acc)) {
    at <- list(vint1 = w)
    cells <- ddm_cells(fin, dp, at, lo, hi)
    mass <- mass + cells[["mass"]]
    set.seed(DDM_SEED + w)
    draws <- fin[["sim"]](lapply(dp, function(v) rep(v, DDM_N)),
                          list(vint1 = rep(w, DDM_N)), DDM_N)
    lab <- paste0("lba[winner ", w, "]")
    expect_true(all(draws > dp[["ndt"]]),
                label = paste0(lab, " above the non-decision time"))
    expect_ddm_gof(draws, cells, lab)
  }
  # the race is over the accumulators, so their winning probabilities
  # are a distribution
  expect_equal(mass, 1, tolerance = 1e-3)
})

# --------------------------------------------------------------- coverage

test_that("every family this package defines declares a simulator", {
  # the extension's families all draw through the family seam, so none
  # of them should be refused by frm_simulate() or simulate()
  fams <- list(wiener = wiener(), lba = lba(2L), gddm = gddm())
  for (nm in names(fams)) {
    f <- fams[[nm]]
    # gddm derives its simulator from the data, so ask the family it
    # would become rather than the family as written
    if (!is.null(f[["family_finalize"]]) && is.null(f[["sim"]])) {
      g <- gddm_simulate(20L)
      av <- list(dec = as.numeric(g$upper), vint1 = as.numeric(g$cond))
      if (!is.null(g$coh)) av$vreal1 <- g$coh
      f <- f[["family_finalize"]](f, g$rt, av)
    }
    expect_false(is.null(f[["sim"]]),
                 label = paste0(nm, " declares a simulator"))
    expect_null(f[["sim_refusal"]],
                label = paste0(nm, " has no refusal to state"))
  }
})
