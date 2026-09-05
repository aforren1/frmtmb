## gddm() against PyDDM, which is the canonical implementation of the
## model this family fits (Shinn, Lam and Murray 2020).
##
## WHY THIS TIER EXISTS. Everything else that checks gddm() checks it
## against itself. test-gddm-solver.R compares the flux with this
## package's own Navarro-Fuss density, which is only defined for a
## CONSTANT drift and FIXED bounds, so it says nothing about leak,
## collapse, coherence, start-point variability or lapse.
## test-gddm-gradients.R proves the derivative matches the value it is
## the derivative of, which is true of a wrong density too.
## test-gddm-recovery.R is circular by construction, because
## gddm_simulate() draws from the solver's own density on purpose. This
## file is the only place a generalized component is compared with a
## number that did not come out of this package.
##
## WHAT IS COMPARED. Three things per case: the defective density at
## both walls on a common time grid, read through the family's own
## lpdf; the two boundary masses; and the log-likelihood of a small
## fixed dataset whose response times fall BETWEEN grid nodes, so the
## likelihood's interpolation is exercised rather than bypassed.
##
## WHAT IS NOT COMPARED. The leading edge. Within a few steps of zero
## both solvers spread a little mass where the true first-passage
## density is exponentially small, for the same structural reason, so a
## comparison there would agree about an artefact. The fixture starts
## at a decision time of 0.2 s, which is where gddm_control()'s help
## page already begins claiming accuracy. Also not compared: anything
## outside the ten cases below, and any grid other than the ones named
## here.
##
## NO PYTHON RUNS HERE. The fixture is frozen, generated once by
## dev/gddm-pyddm-reference.py, and the metadata assertions below make
## a regenerated fixture that moved fail loudly instead of quietly
## redefining what is correct.

fx <- function(nm) {
  testthat::test_path("fixtures", paste0("gddm-pyddm-", nm, ".csv"))
}
gr_dens_tab <- utils::read.csv(fx("density"))
gr_ll_tab <- utils::read.csv(fx("loglik"))
gr_meta_tab <- utils::read.csv(fx("meta"), colClasses = "character")

gr_meta <- function(case, key) {
  v <- gr_meta_tab$value[gr_meta_tab$case == case & gr_meta_tab$key == key]
  if (!length(v)) NA_character_ else v[[1L]]
}
gr_num <- function(case, key) as.numeric(gr_meta(case, key))

## ------------------------------------------------------------------
## The parameter translation, written out once
## ------------------------------------------------------------------
##
## gddm() and PyDDM parameterize the same model differently, and three
## of the differences are sign or reciprocal traps. This table is the
## contract; dev/gddm-pyddm-reference.py holds the same table on the
## Python side and the fixture records the gddm() numbers, so the two
## halves cannot drift apart without the metadata check below failing.
##
##   quantity              gddm()                  PyDDM
##   ------------------------------------------------------------------
##   noise                 dx = a(x,t) dt + dW     NoiseConstant(noise=1)
##   boundary              +/- B(t), B(0) = bs/2   B = bs/2
##   start point           bias in (0,1), the      ICPoint(x0 =
##                         fraction of the           (bs/2)(2 bias - 1))
##                         separation above the
##                         lower boundary
##   constant drift        a = mu                  DriftConstant(drift = mu)
##   leak                  a = mu - leak x         DriftLinear(drift = mu,
##                                                   x = -leak, t = 0)
##   exponential collapse  B = (bs/2) exp(-t/tau)  BoundCollapsingExponential(
##                                                   B = bs/2, tau = 1/tau)
##   linear collapse       B = (bs/2)(1 -          BoundCollapsingLinear(
##                           kappa t/t_max)          B = bs/2,
##                                                   t = (bs/2) kappa/t_max)
##   coherence drift       a = sign(C) mu          DriftConstant(drift =
##                           (|C|/cmax)^alpha        that number)
##   start variability     uniform, half width     ICRange(sz = sz bs/2)
##                           sz of the half
##                           separation
##   non-decision time     ndt                     OverlayNonDecision(
##                                                   nondectime = ndt)
##
## THE LEAK SIGN. gddm()'s `leak` is the paper's l: POSITIVE is leaky
## integration, which pulls the accumulator back toward zero. PyDDM's
## DriftLinear takes the coefficient of x directly, which is the
## NEGATIVE of that. R/gddm.R:404 says so; the `leak` and `unstable`
## cases are what would catch it being got wrong, one for each sign.
##
## THE BOUND PARAMETERIZATION. gddm() names the SEPARATION `bs` and
## puts the walls at plus and minus bs/2, so that a separation fitted
## here is comparable with one from wiener(). PyDDM names the HALF
## separation B. And gddm()'s exponential time constant `tau` is the
## RECIPROCAL of PyDDM's, which parameterizes the same bound by a rate
## (R/gddm.R:523).
##
## THE LAPSE. gddm() renormalizes each condition's pair of densities
## and THEN mixes, adding a density of 0.5 lapse / t_max at each wall.
## PyDDM's OverlayUniformMixture adds 0.5 c / len(t_domain) of MASS to
## each bin of a defective density, which is a density of
## 0.5 c / (T_dur + dt): one bin wider. The generator checks the two
## formulas against each other and the fixture carries gddm()'s, so
## `lapse_formula_max_abs_diff` below is the evidence that the
## difference is only the window and not the mixture.
##
## THE RENORMALIZATION. gddm() divides by the realized mass, so its
## likelihood is conditional on a response inside the window; PyDDM
## keeps the undecided mass. The fixture is renormalized to match and
## records the undecided mass it dropped, which for these cases runs
## from 6.5e-5 to 6.8e-3.

gr_spec <- list(
  constant = list(args = list(),
                  par = c(mu = 1.2, bs = 2, bias = 0.5, ndt = 0.25)),
  bias_ndt = list(args = list(),
                  par = c(mu = 0.8, bs = 2, bias = 0.65, ndt = 0.40)),
  leak = list(args = list(drift = list(gddm_drift_constant(),
                                       gddm_drift_leak())),
              par = c(mu = 1.5, leak = 1.2, bs = 2, bias = 0.5,
                      ndt = 0.25)),
  unstable = list(args = list(drift = list(gddm_drift_constant(),
                                           gddm_drift_leak())),
                  par = c(mu = 0.6, leak = -0.8, bs = 2, bias = 0.5,
                          ndt = 0.25)),
  exp_collapse = list(args = list(bound = gddm_bound_exponential()),
                      par = c(mu = 1.2, bs = 2, tau = 1.5, bias = 0.5,
                              ndt = 0.25)),
  lin_collapse = list(args = list(bound = gddm_bound_linear()),
                      par = c(mu = 1.2, bs = 2, kappa = 0.6, bias = 0.5,
                              ndt = 0.25)),
  coh_high = list(args = list(drift = gddm_drift_coherence(cmax = 0.512)),
                  par = c(mu = 2, alpha = 1.3, bs = 2, bias = 0.5,
                          ndt = 0.25),
                  cmax = 0.512, coh = 0.512),
  coh_low = list(args = list(drift = gddm_drift_coherence(cmax = 0.512)),
                 par = c(mu = 2, alpha = 1.3, bs = 2, bias = 0.5,
                         ndt = 0.25),
                 cmax = 0.512, coh = 0.256),
  sz_uniform = list(args = list(start = gddm_start_uniform()),
                    par = c(mu = 1, bs = 2, bias = 0.5, sz = 0.2,
                            ndt = 0.25)),
  lapse = list(args = list(lapse = "uniform"),
               par = c(mu = 1.2, bs = 2, bias = 0.5, ndt = 0.25,
                       lapse = 0.05)))

GR_CASES <- names(gr_spec)
GR_TMAX <- 3

## Tolerances. Every one is about twice a MEASURED discrepancy, each
## measured with the rule the assertion below actually uses, and all of
## them are listed in dev/gddm-reference-findings.md. At the shipped
## grid the worst measured log-density difference is 1.62e-3
## (bias_ndt) and the worst log-likelihood difference 2.9e-3 (leak).
GR_TOL_DLOG <- 2.5e-3
GR_TOL_LOGLIK <- 5e-3

## The boundary-mass tolerance is measured with gr_mass()'s TRAPEZOID
## rule, which is the rule the assertion uses, and the worst is 1.53e-4
## on coh_low. 3e-4 is 1.96 times that. It was widened from 2e-4
## rather than left alone: 2e-4 was justified against 7.1e-5, which is
## the worst under the RECTANGLE rule, and against the shipped rule
## that left only 1.31x, the thinnest margin in the file.
##
## Why the two rules differ by about 1e-4, which is worth knowing
## before reading the number as a solver disagreement: both sides of
## this comparison are rectangle sums. gd_densities()'s renormalize
## divides by (sum(pu) + sum(pl)) * dt, and PyDDM's Solution.prob() is
## a plain sum of per-bin masses. The trapezoid rule is kept anyway,
## because it is the mass of the piecewise-linear density the
## likelihood actually interpolates and the rule gd_draw() samples
## from, and because it is the right rule for the lapse case, whose
## density is nonzero at both ends of the window. The cost of keeping
## it is that roughly 1e-4 of the number below is quadrature
## convention rather than a difference between the two solvers.
GR_TOL_MASS <- 3e-4

## exp_collapse is the one case where the tolerance is set by the
## REFERENCE rather than by gddm(). It is the only case with no closed
## form AND a moving bound, so PyDDM solves it with backward Euler on a
## grid-snapped bound: first order in dt AND in dx. Its own convergence
## residual, recorded in the fixture, is 2.6e-3.
##
## The 3.25e-3 gddm() sits from the fixture is that residual, and it
## was measured rather than argued by analogy. Solving the case one
## level finer than the fixture and forming the first-order Richardson
## limit, gddm()'s distance from PyDDM is:
##
##   gddm grid        vs fixture   vs one finer   vs PyDDM's limit
##   dt .01  ny 201     3.250e-3      1.965e-3        8.080e-4
##   dt .005 ny 401     2.832e-3      1.547e-3        2.779e-4
##   dt .0025 ny 801    2.764e-3      1.459e-3        1.689e-4
##
## So gddm() does not sit at a floor: against where PyDDM is GOING it
## converges, and at the shipped grid it is 8.1e-4 away, inside the
## 7.3e-4 to 1.6e-3 band the other nine cases occupy. The whole
## remaining 2.8e-3 is the reference's own discretization error. The
## measured 3.25e-3 gets about 1.5 times, as everything else does; the
## tolerance is wide because the fixture is wrong by 2.6e-3, not
## because gddm() is. What it costs is discrimination: a gddm() error
## on this case below about 2e-3 would be invisible here, which is why
## the convergence tier below also asserts on this case.
GR_TOL_DLOG_EXP <- 5e-3

## ------------------------------------------------------------------
## Driving the family
## ------------------------------------------------------------------

## The finalized family, its per-row index data and its density, the way
## frm() would assemble them. gddm() installs the density AND the grid
## in family_finalize(), so a family as written has neither: the object
## the density is read off has to be the finalized one.
gr_finalize <- function(case, y, upper, dt, ny) {
  sp <- gr_spec[[case]]
  p <- as.list(sp$par)
  ctl <- gddm_control(t_max = GR_TMAX, dt = dt, ny = as.integer(ny),
                      # named rather than taken from the data, so the
                      # non-decision-time shift kernel is the same
                      # length whatever response times are passed in
                      max_ndt = p$ndt)
  fam <- do.call(gddm, c(sp$args, list(control = ctl)))
  at <- list(dec = as.numeric(upper),
             vint1 = rep(1, length(y)))
  if (!is.null(sp$coh)) at$vreal1 <- rep(sp$coh, length(y))
  fin <- fam[["family_finalize"]](fam, y, at)
  at <- c(at, fin[["aterm_data"]](y, at))
  dp <- lapply(fin[["dpars"]], function(nm) rep(p[[nm]], length(y)))
  names(dp) <- fin[["dpars"]]
  list(fam = fin, at = at, dp = dp)
}

gr_lpdf <- function(case, y, upper, dt = 0.01, ny = 201L) {
  f <- gr_finalize(case, y, upper, dt, ny)
  exp(as.numeric(f$fam[["lpdf"]](y, f$dp, f$at)))
}

## The mass at each wall, from the same density the likelihood reads.
## The TRAPEZOID rule, not the sum: the likelihood interpolates
## linearly between nodes, so the trapezoid rule is the mass of the
## function it actually evaluates, and it is the rule gd_draw() already
## uses to turn that density into draws. The difference matters only
## for the lapse case, whose density is nonzero at both ends of the
## window.
gr_mass <- function(case, dt = 0.01, ny = 201L) {
  p <- as.list(gr_spec[[case]]$par)
  y <- rep(p$ndt + 0.5, 2L)
  f <- gr_finalize(case, y, c(1, 0), dt, ny)
  bag <- f$fam[["gddm"]]
  v <- as.numeric(gd_densities(f$dp, bag$comp, bag$ctl, f$at[[".gddm"]]))
  nb <- bag$ctl$nt + 1L
  tr <- function(z) (sum(z) - (z[[1L]] + z[[nb]]) / 2) * dt
  c(up = tr(v[seq_len(nb)]), lo = tr(v[nb + seq_len(nb)]))
}

## What the fixture's boundary mass means for a case with a lapse:
## gddm() renormalizes to unit mass and then mixes, so each wall gains
## half the lapse. The fixture stores the PRE-lapse probabilities.
gr_expected_mass <- function(case) {
  p <- as.list(gr_spec[[case]]$par)
  lp <- if (is.null(p$lapse)) 0 else p$lapse
  c(up = (1 - lp) * gr_num(case, "prob_upper") + lp / 2,
    lo = (1 - lp) * gr_num(case, "prob_lower") + lp / 2)
}

## Everything a case needs from the density, in ONE call. Both walls and
## the fixed dataset's response times are all condition 1 with the same
## parameters, so the family solves once and gathers; asking three times
## would triple the cost of the file for nothing. Cached because the
## curve, the log-likelihood and the mass are asserted in three separate
## tests off the same solve.
gr_cache <- new.env(parent = emptyenv())

gr_eval <- function(case) {
  hit <- gr_cache[[case]]
  if (!is.null(hit)) return(hit)
  d <- gr_dens_tab[gr_dens_tab$case == case, ]
  l <- gr_ll_tab[gr_ll_tab$case == case, ]
  n <- nrow(d)
  m <- nrow(l)
  got <- gr_lpdf(case, c(d$t, d$t, l$rt),
                 c(rep(1, n), rep(0, n), l$upper))
  out <- list(d = d, l = l, up = got[seq_len(n)],
              lo = got[n + seq_len(n)], ll = got[2L * n + seq_len(m)])
  gr_cache[[case]] <- out
  out
}

## Worst absolute log-density difference over a case's grid, both walls.
gr_worst_dlog <- function(case, dt = 0.01, ny = 201L) {
  if (dt == 0.01 && ny == 201L) {
    e <- gr_eval(case)
    return(max(abs(log(e$up) - log(e$d$dens_upper)),
               abs(log(e$lo) - log(e$d$dens_lower))))
  }
  d <- gr_dens_tab[gr_dens_tab$case == case, ]
  n <- nrow(d)
  got <- gr_lpdf(case, c(d$t, d$t), c(rep(1, n), rep(0, n)), dt, ny)
  max(abs(log(got[seq_len(n)]) - log(d$dens_upper)),
      abs(log(got[n + seq_len(n)]) - log(d$dens_lower)))
}

# ------------------------------------------------------------------
# The fixture is the one this test was written against
# ------------------------------------------------------------------

test_that("the fixture records the environment it was built in", {
  ## A regenerated fixture is allowed; a regenerated fixture that
  ## quietly moved is not. These are the assertions that make
  ## regeneration a visible act.
  expect_identical(gr_meta("_fixture", "pyddm"), "0.9.0")
  expect_identical(gr_meta("_fixture", "numpy"), "2.5.2")
  expect_identical(gr_meta("_fixture", "scipy"), "1.18.1")
  expect_match(gr_meta("_fixture", "python"), "^3[.]13[.]")
  expect_identical(gr_meta("_fixture", "n_cases"), "10")
  expect_equal(gr_num("_fixture", "T_dur"), GR_TMAX)
  expect_equal(gr_num("_fixture", "grid_step"), 0.02)
  expect_equal(gr_num("_fixture", "dens_floor"), 1e-6)
  expect_setequal(setdiff(unique(gr_dens_tab$case), "_fixture"), GR_CASES)
  expect_setequal(unique(gr_ll_tab$case), GR_CASES)

  for (cs in GR_CASES) {
    ## the parameters the test builds the model with ARE the parameters
    ## the reference was solved at, one by one
    for (nm in names(gr_spec[[cs]]$par)) {
      expect_equal(gr_num(cs, paste0("par_", nm)),
                   unname(gr_spec[[cs]]$par[[nm]]),
                   label = paste0(cs, " par_", nm))
    }
    expect_equal(gr_num(cs, "T_dur"), GR_TMAX, label = paste0(cs, " T_dur"))
    ## the reference's own convergence residual is smaller than the
    ## tolerance granted, so the tolerance is measuring gddm() and not
    ## the reference. exp_collapse is the exception and is named.
    tol <- if (cs == "exp_collapse") GR_TOL_DLOG_EXP else GR_TOL_DLOG
    expect_lt(gr_num(cs, "ref_residual_dlog"), tol)
  }
  ## the coherence covariate and its scale, which live outside `par`
  expect_equal(gr_num("coh_high", "par_coh"), 0.512)
  expect_equal(gr_num("coh_low", "par_coh"), 0.256)
  expect_equal(gr_num("coh_high", "par_cmax"), 0.512)
  ## and the drift the nonlinearity is supposed to produce
  expect_equal(gr_num("coh_high", "drift_value"), 2)
  expect_equal(gr_num("coh_low", "drift_value"), 2 * 0.5^1.3)
})

test_that("the fixture's payload hashes to what the generator wrote", {
  ## Rebuilt from the PARSED numbers rather than from the file's bytes,
  ## so a checkout that rewrote the line endings still hashes the same.
  skip_if_not(exists("sha256sum", where = asNamespace("tools")),
              "tools::sha256sum() needs R >= 4.5")
  e <- function(x) sprintf("%.12e", x)
  txt <- paste0(
    paste(c("case,t,dens_upper,dens_lower",
            paste0(gr_dens_tab$case, ",", e(gr_dens_tab$t), ",",
                   e(gr_dens_tab$dens_upper), ",",
                   e(gr_dens_tab$dens_lower))), collapse = "\n"), "\n",
    paste(c("case,rt,upper,dens",
            paste0(gr_ll_tab$case, ",", e(gr_ll_tab$rt), ",",
                   gr_ll_tab$upper, ",", e(gr_ll_tab$dens))),
          collapse = "\n"), "\n")
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  writeBin(charToRaw(txt), f)
  expect_identical(unname(tools::sha256sum(f)),
                   gr_meta("_fixture", "content_sha256"))
})

test_that("PyDDM's uniform lapse mixture is the one gddm() applies", {
  ## The generator built the lapse two ways, once with PyDDM's own
  ## OverlayUniformMixture and once by hand, and they agreed to machine
  ## precision. So the only difference left between the two packages'
  ## lapse is the window it is spread over, which is stated rather than
  ## absorbed: PyDDM uses T_dur + dt and gddm() uses t_max.
  expect_lt(gr_num("lapse", "lapse_formula_max_abs_diff"), 1e-12)
  expect_equal(gr_num("lapse", "lapse_gddm_density"), 0.5 * 0.05 / GR_TMAX)
  expect_lt(abs(gr_num("lapse", "lapse_pyddm_density") -
                  gr_num("lapse", "lapse_gddm_density")), 2e-6)
})

# ------------------------------------------------------------------
# The densities themselves
# ------------------------------------------------------------------

test_that("every generalized component reproduces PyDDM's density", {
  for (cs in GR_CASES) {
    tol <- if (cs == "exp_collapse") GR_TOL_DLOG_EXP else GR_TOL_DLOG
    expect_lt(gr_worst_dlog(cs), tol, label = paste0("dlog[", cs, "]"))
  }
})

test_that("the two boundary masses match PyDDM's", {
  ## Absolute, not relative. A wall that carries two percent of the
  ## mass would otherwise be held to a tolerance fifty times tighter
  ## than the wall that carries the rest, which is not what a
  ## probability tolerance means.
  for (cs in GR_CASES) {
    got <- gr_mass(cs)
    want <- gr_expected_mass(cs)
    expect_lt(abs(got[["up"]] - want[["up"]]), GR_TOL_MASS,
              label = paste0("mass up [", cs, "]"))
    expect_lt(abs(got[["lo"]] - want[["lo"]]), GR_TOL_MASS,
              label = paste0("mass lo [", cs, "]"))
    ## and the density is a proper one once the window is conditioned
    ## on: renormalize divides by the realized mass, so the two walls
    ## have to carry all of it between them. Not to the last bit,
    ## because renormalize uses the plain sum while this uses the
    ## trapezoid rule, and the two differ by half the density at the
    ## ends of the window.
    expect_lt(abs(got[["up"]] + got[["lo"]] - 1), 3e-4,
              label = paste0("total mass [", cs, "]"))
  }
})

test_that("the fixed dataset's log-likelihood matches PyDDM's", {
  ## The response times here are odd multiples of 0.005 while the grid
  ## is a multiple of 0.01, so every one of them falls BETWEEN two
  ## nodes and the density comes back through the likelihood's own
  ## linear interpolation. The curve test above reads nodes; this one
  ## reads the code path a fit uses.
  for (cs in GR_CASES) {
    e <- gr_eval(cs)
    expect_false(any(abs(e$l$rt / 0.01 - round(e$l$rt / 0.01)) < 1e-9),
                 label = paste0("off-node response times [", cs, "]"))
    ## absolute: the tolerance is a statement about the log-likelihood
    ## of six observations, not a fraction of whatever it happens to be
    expect_lt(abs(sum(log(e$ll)) - gr_num(cs, "loglik")), GR_TOL_LOGLIK,
              label = paste0("loglik[", cs, "]"))
  }
})

# ------------------------------------------------------------------
# The agreement is a limit, not a coincidence
# ------------------------------------------------------------------

test_that("refining gddm's grid moves it toward PyDDM, not merely near it", {
  ## A tolerance on its own would be satisfied by a solver that happened
  ## to land inside it. What says the reference is the thing gddm() is
  ## converging TO is that the discrepancy falls as the grid is refined.
  ## Three cases: the calibration case, a moving bound where PyDDM has
  ## a closed form, and the case whose sequence flattens because the
  ## FIXTURE stops improving rather than because gddm() does.
  skip_on_cran()
  grids <- list(c(0.02, 101), c(0.01, 201), c(0.005, 401))
  for (cs in c("constant", "lin_collapse", "exp_collapse")) {
    e <- vapply(grids, function(g) gr_worst_dlog(cs, g[[1L]], g[[2L]]),
                numeric(1))
    expect_lt(e[[2L]], e[[1L]],
              label = paste0("default beats coarse [", cs, "]"))
    expect_lt(e[[3L]], e[[2L]],
              label = paste0("fine beats default [", cs, "]"))
  }
  ## and the two that converge do so at the second order both schemes
  ## claim, quartering rather than halving
  for (cs in c("constant", "lin_collapse")) {
    e <- vapply(grids, function(g) gr_worst_dlog(cs, g[[1L]], g[[2L]]),
                numeric(1))
    expect_gt(e[[1L]] / e[[2L]], 3, label = paste0("second order [", cs, "]"))
    expect_gt(e[[2L]] / e[[3L]], 3, label = paste0("second order [", cs, "]"))
  }
})

# ------------------------------------------------------------------
# A third reference, so that two solvers cannot be wrong together
# ------------------------------------------------------------------

## Euler-Maruyama on the SDE itself. Absorption is only ever noticed at
## a monitoring time, so a simulated first passage is LATE, and late by
## order sqrt(h) rather than h (Broadie, Glasserman and Kou 1997). Two
## gaps, h and 2h, are driven by the SAME Brownian increments, so the
## difference the extrapolation leans on is nearly noiseless and the
## extrapolated value carries one level's Monte Carlo error instead of
## the sum of two. Measured, that pairing is worth a factor of four.
gr_euler <- function(h, n, mu, B0, tau, x0, tmax, seed) {
  set.seed(seed)
  sh <- sqrt(h)
  xf <- rep(x0, n); xc <- rep(x0, n)
  tf <- rep(NA_real_, n); tc <- rep(NA_real_, n)
  wf <- rep(NA_integer_, n); wc <- rep(NA_integer_, n)
  lf <- rep(TRUE, n); lc <- rep(TRUE, n)
  acc <- numeric(n)                 # the coarse step's half-built increment
  live <- rep(TRUE, n)
  for (k in seq_len(floor(tmax / h))) {
    idx <- which(live)
    if (!length(idx)) break
    dw <- sh * stats::rnorm(length(idx))
    b <- B0 * exp(-(k * h) / tau)
    j <- idx[lf[idx]]
    if (length(j)) {
      xf[j] <- xf[j] + mu * h + dw[lf[idx]]
      hu <- xf[j] >= b
      hit <- hu | xf[j] <= -b
      if (any(hit)) {
        tf[j[hit]] <- k * h
        wf[j[hit]] <- as.integer(hu[hit])
        lf[j[hit]] <- FALSE
      }
    }
    jc <- idx[lc[idx]]
    if (length(jc)) {
      acc[jc] <- acc[jc] + dw[lc[idx]]
      if (k %% 2L == 0L) {
        xc[jc] <- xc[jc] + mu * 2 * h + acc[jc]
        acc[jc] <- 0
        hu <- xc[jc] >= b
        hit <- hu | xc[jc] <= -b
        if (any(hit)) {
          tc[jc[hit]] <- k * h
          wc[jc[hit]] <- as.integer(hu[hit])
          lc[jc[hit]] <- FALSE
        }
      }
    }
    live <- lf | lc
  }
  list(tf = tf, wf = wf, tc = tc, wc = wc)
}

test_that("a path simulation says the two grid solvers do not share a bias", {
  ## Both gddm() and PyDDM march a Fokker-Planck equation on a grid, so
  ## agreeing with each other does not by itself rule out an error they
  ## make in the same direction. This runs the stochastic differential
  ## equation instead, which shares no machinery with either, on the
  ## case where a shared bias would be least surprising: the
  ## exponentially collapsing bound, where PyDDM cannot use
  ## Crank-Nicolson at all and gddm() only can because it changes
  ## variable to hold the walls still.
  skip_on_cran()
  p <- as.list(gr_spec$exp_collapse$par)
  n <- 80000L
  h <- 0.001
  r <- gr_euler(h, n, mu = p$mu, B0 = p$bs / 2, tau = p$tau, x0 = 0,
                tmax = GR_TMAX, seed = 20250905L)
  both <- !is.na(r$tf) & !is.na(r$tc)
  ## everything absorbs well inside the window under a collapsing bound,
  ## so the two resolutions decide the same paths
  expect_gt(mean(both), 0.999)

  ## Richardson in sqrt(h), with the standard error taken from the
  ## PAIRED combination rather than from two independent means
  rr <- 1 / (sqrt(2) - 1)
  rich <- function(af, ac) {
    z <- af + rr * (af - ac)
    c(est = mean(z), se = stats::sd(z) / sqrt(length(z)))
  }
  p_up <- rich(as.numeric(r$wf[both] == 1L), as.numeric(r$wc[both] == 1L))
  m_all <- rich(r$tf[both], r$tc[both])

  ## The reference's own answers for the same two functionals. The
  ## unconditional mean is the two conditional means the fixture stores,
  ## weighted by the two boundary probabilities.
  ref_p <- gr_num("exp_collapse", "prob_upper")
  ref_m <- ref_p * gr_num("exp_collapse", "mean_dt_upper") +
    (1 - ref_p) * gr_num("exp_collapse", "mean_dt_lower")

  ## gddm()'s own answers, from the density the likelihood reads
  f <- gr_finalize("exp_collapse", rep(p$ndt + 0.5, 2L), c(1, 0),
                   0.005, 401L)
  bag <- f$fam[["gddm"]]
  v <- as.numeric(gd_densities(f$dp, bag$comp, bag$ctl, f$at[[".gddm"]]))
  nb <- bag$ctl$nt + 1L
  tg <- seq(0, bag$ctl$t_max, by = bag$ctl$dt) - p$ndt
  vu <- v[seq_len(nb)]
  vl <- v[nb + seq_len(nb)]
  gd_p <- sum(vu) / (sum(vu) + sum(vl))
  gd_m <- sum(tg * (vu + vl)) / sum(vu + vl)

  ## First, the two grid solvers agree with each other far more closely
  ## than the path simulation can resolve. That is the premise; if it
  ## failed, the comparison below would be measuring something else.
  expect_lt(abs(gd_p - ref_p), 1e-3)
  expect_lt(abs(gd_m - ref_m), 1e-3)

  ## Then the path simulation agrees with both, inside four standard
  ## errors. Measured at n = 200000 the gaps are 1.2 and 1.5 standard
  ## errors, and across four independent seeds, 800000 paths, nothing
  ## larger than about 1e-3 in P(upper) is resolvable and nothing is
  ## detected: it is the MAGNITUDES that carry the no-shared-bias
  ## claim. Which SIDE of the reference a run lands on carries nothing
  ## and is not read here. The two functionals are anti-correlated by
  ## construction, because under a positive drift more upper hits means
  ## faster decisions, so opposite directions is the null expectation
  ## rather than evidence, and the observed signs flip with the seed.
  ## The seed is fixed, so this is a deterministic assertion; single
  ## seeds have been seen to reach 1.84 standard errors, so four is
  ## real headroom for an arithmetic change but not lavish.
  for (nm in c("pyddm", "gddm")) {
    pv <- if (nm == "pyddm") ref_p else gd_p
    mv <- if (nm == "pyddm") ref_m else gd_m
    expect_lt(abs(p_up[["est"]] - pv), 4 * p_up[["se"]],
              label = paste0("P(upper) vs ", nm))
    expect_lt(abs(m_all[["est"]] - mv), 4 * m_all[["se"]],
              label = paste0("mean decision time vs ", nm))
  }
  ## and the check has to have some resolution to be worth making: a
  ## shared bias of a few parts in a thousand would have shown
  expect_lt(p_up[["se"]], 3e-3)
  expect_lt(m_all[["se"]], 3e-3)
})
