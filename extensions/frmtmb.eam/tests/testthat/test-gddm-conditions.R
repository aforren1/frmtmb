## The within-condition contract, and its enforcement.
##
## gd_densities() solves the Fokker-Planck equation once per condition
## and reads every parameter at that condition's FIRST ROW. A parameter
## that varies inside a condition therefore never reaches the
## likelihood, and through frmtmb.eam 0.7.0 that was silent: the model
## fitted, no warning was raised, and adding 100 to a covariate on 118
## of 120 rows left the objective BITWISE identical at the same
## parameter vector, for `mu`, `ndt`, `bs` and `bias` alike. See
## dev/gddm-findings.md and dev/gddm-scripts/.
##
## Every refusal below is a model that FITTED on 0.7.0, so this file
## fails on the unfixed package rather than erroring on an absent
## symbol.

gc_toy <- function(n = 120L, seed = 5L) {
  set.seed(seed)
  d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                     control = gddm_control(t_max = 2))
  d$cond <- rep(1:2, length.out = n)
  d$x <- rnorm(n)
  ## a subject factor that CROSSES the condition index, which is the
  ## random-effect form of the same defect
  d$s <- factor(rep(1:4, length.out = n))
  d
}
gc_small <- function() gddm_control(t_max = 2, dt = 0.05, ny = 51L)

## Nothing here needs a taped objective: the check runs at frame
## assembly, so dry_run = "frame" exercises it and skips the solve.
gc_frame <- function(spec, data, family = gddm(control = gc_small())) {
  frm(spec, family = family, data = data, dry_run = "frame")
}

test_that("a dpar varying within a condition is refused, and named", {
  d <- gc_toy()
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "the parameter `mu` is not constant within every condition",
    fixed = TRUE)
  ## the variable is named too, because the remedy is to put it in the
  ## index and gddm_conditions() takes variable names
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "`x` varies inside condition 1", fixed = TRUE)
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "gddm_conditions()", fixed = TRUE)
})

test_that("every parameter of the family is covered, not just mu", {
  ## the loop over dpars is the part that can fail open by hitting
  ## `next` on every name, so each of the four is constructed
  d <- gc_toy()
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ x,
                bias = 0.5), d),
    "the parameter `ndt` is not constant", fixed = TRUE)
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ x, ndt ~ 1,
                bias = 0.5), d),
    "the parameter `bs` is not constant", fixed = TRUE)
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1,
                bias ~ x), d),
    "the parameter `bias` is not constant", fixed = TRUE)
  ## and a parameter a component brings, rather than one the base
  ## family always has
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, tau ~ x, ndt ~ 1,
                bias = 0.5), d,
             gddm(bound = gddm_bound_exponential(),
                  control = gc_small())),
    "the parameter `tau` is not constant", fixed = TRUE)
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1,
                lapse ~ x, bias = 0.5), d,
             gddm(lapse = "uniform", control = gc_small())),
    "the parameter `lapse` is not constant", fixed = TRUE)
})

test_that("a random effect whose grouping crosses the index is refused", {
  ## The Z half of the same defect. On 0.7.0 this model fitted and
  ## reported a between-subject standard deviation of 7.6e-11, because
  ## only the subjects sitting on a condition's first row ever reached
  ## the density: the other deviations stayed at exactly zero.
  d <- gc_toy()
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "`s` varies inside condition", fixed = TRUE)
  ## and it is the drift that carries it, so the drift is what is named
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "the parameter `mu` is not constant", fixed = TRUE)
})

test_that("the message names the user's own condition label", {
  ## the index is not renumbered for the message: a user who labelled
  ## their conditions 10 and 20 is told 20, not 2
  d <- gc_toy()
  d$cond <- ifelse(seq_len(nrow(d)) <= 60L, 10L, 20L)
  d$x <- 0
  d$x[[nrow(d)]] <- 1        # varies inside the SECOND condition only
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "`x` varies inside condition 20, and inside 1 of 2 conditions",
    fixed = TRUE)
})

test_that("an index that separates the varying term is accepted", {
  ## The half that makes the refusal worth having. A check that fires
  ## on a correct model is worse than no check, so every refusal above
  ## is paired with the model it becomes once the index is right.
  d <- gc_toy()
  d$x <- rep(c(0, 1), length.out = nrow(d))
  d$cond <- gddm_conditions(d, cond, x)
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                bias = 0.5), d))
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ x, ndt ~ x,
                bias ~ x), d))
  d2 <- gc_toy()
  d2$cond <- gddm_conditions(d2, cond, s)
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
                bias = 0.5), d2))
  ## an intercept-only model has nothing that can vary
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1,
                bias = 0.5), gc_toy()))
})

test_that("a computed design column is not refused for its own rounding", {
  ## poly() orthogonalizes over the whole column, so two rows built
  ## from bitwise identical inputs come back a few hundred ulp apart.
  ## An exact comparison, which is the shape gd_check_response() uses
  ## on the vreal() covariates, refuses this correct model; the check
  ## therefore carries a tolerance, scaled by the two entries being
  ## compared with a floor under it for the entries near zero.
  d <- gc_toy()
  d$z <- rep(c(0, 0.128, 0.256, 0.512), length.out = nrow(d))
  d$cond <- gddm_conditions(d, z)
  P <- stats::poly(d$z, 2)
  ref <- match(d$cond, d$cond)
  dev <- abs(P - P[ref, , drop = FALSE])
  ## the rounding is real: this is the measurement the tolerance exists
  ## for, expressed against the tolerance the guard itself applies to
  ## each pair rather than as an absolute number
  expect_gt(max(dev), 0)
  tol <- 1e-8 * pmax(abs(P), abs(P[ref, , drop = FALSE])) +
    gd_col_floor(as.vector(P))
  expect_lt(max(dev / tol), 1)
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ poly(z, 2), bs ~ 1, ndt ~ 1,
                bias = 0.5), d))
  ## and a step above the tolerance is still refused, so the tolerance
  ## has not opened the guard
  d$z[[2L]] <- d$z[[2L]] + 1e-4 * max(abs(d$z))
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ poly(z, 2), bs ~ 1, ndt ~ 1,
                bias = 0.5), d),
    "the parameter `mu` is not constant", fixed = TRUE)
})

test_that("a condition of one row cannot vary and is not refused", {
  ## the degenerate end of the index: one solve per row is the most
  ## resolution the design can carry, and the first row IS every row
  d <- gc_toy()
  d$cond <- seq_len(nrow(d))
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x + s, bs ~ x, ndt ~ x,
                bias ~ x), d))
})

test_that("the refusal is gddm's alone and leaves other families be", {
  ## the check runs on every frame frmtmb assembles once this package
  ## is loaded, so a family that has no conditions must be untouched
  ## however the data vary
  d <- gc_toy()
  d$win <- d$upper + 1L
  expect_no_error(
    frm(bf(rt ~ x + s), family = gaussian(), data = d,
        dry_run = "frame"))
  expect_no_error(
    frm(bf(rt | dec(upper) ~ x + (1 | s)), family = wiener(),
        data = d, dry_run = "frame"))
  expect_no_error(
    frm(bf(rt | vint(win) ~ x + (1 | s)), family = lba(2), data = d,
        dry_run = "frame"))
})

test_that("the vreal() check this one is modelled on still fires", {
  ## the older, exact check on the drift covariate. It is the precedent
  ## and the two must not drift apart: this pins that a covariate
  ## varying inside a condition is still refused by its own message.
  d <- gc_toy()
  d$coh <- rep(c(0.3, 0.3, 0.7, 0.7), length.out = nrow(d))
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) + vreal(coh) ~ 1, bs ~ 1,
                ndt ~ 1, bias = 0.5), d,
             gddm(drift = gddm_drift_coherence(),
                  control = gc_small())),
    "vreal1 is not constant within every condition", fixed = TRUE)
  ## and the same covariate aligned with the index is accepted
  d$coh <- ifelse(d$cond == 1L, 0.3, 0.7)
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) + vreal(coh) ~ 1, bs ~ 1,
                ndt ~ 1, bias = 0.5), d,
             gddm(drift = gddm_drift_coherence(),
                  control = gc_small())))
})

test_that("the check reads the condition under either spelling", {
  ## vint() numbers positionally: alongside dec() the condition is the
  ## FIRST vint() value, and inside vint(upper, cond) it is the second.
  ## A check that read a fixed slot would compare the boundary
  ## indicator instead, which varies inside every condition by design
  ## and would refuse every model.
  d <- gc_toy()
  ## set outright: the draw at this seed lands on one boundary only,
  ## and dec() needs a two-level factor. Nothing here fits, so which
  ## row took which boundary does not matter.
  d$resp <- factor(rep(c("lower", "upper"), length.out = nrow(d)),
                   levels = c("lower", "upper"))
  expect_no_error(
    frm(bf(rt | dec(resp) + vint(cond) ~ 1, bs ~ 1, ndt ~ 1,
           bias = 0.5),
        family = gddm(control = gc_small()), data = d,
        dry_run = "frame"))
  expect_error(
    frm(bf(rt | dec(resp) + vint(cond) ~ x, bs ~ 1, ndt ~ 1,
           bias = 0.5),
        family = gddm(control = gc_small()), data = d,
        dry_run = "frame"),
    "the parameter `mu` is not constant", fixed = TRUE)
})

test_that("a non-finite entry cannot open the tolerance", {
  ## The tolerance is scaled by the column's largest entry. Taken over
  ## the whole column that maximum is Inf the moment one entry is, the
  ## tolerance is Inf, and nothing is ever refused again: the guard
  ## fails open on one infinite covariate. The scale is therefore taken
  ## over the FINITE entries, and the infinite ones are compared
  ## exactly. Constructed on the helper, because a column of Inf does
  ## not survive as far as a fitted model.
  v <- c(1, 2, Inf, Inf)
  gi <- c(1L, 1L, 2L, 2L)
  first <- c(1L, 3L)
  expect_identical(gd_varying_groups(v, gi, first), 1L)
  ## one Inf against one finite value inside a condition is a
  ## difference, not a rounding
  expect_identical(gd_varying_groups(c(1, Inf, 5, 5), gi, first), 1L)
  ## and a column that really is constant is still constant
  expect_length(gd_varying_groups(c(1, 1, Inf, Inf), gi, first), 0L)
  ## two missing values in one condition are the same datum
  expect_length(gd_varying_groups(c(NA, NA, 5, 5), gi, first), 0L)
  expect_identical(gd_varying_groups(c(NA, 2, 5, 5), gi, first), 1L)
})

test_that("the tolerance does not widen with the column's largest entry", {
  ## The first shipped form was 1e-8 times the COLUMN MAXIMUM, applied
  ## to every pair in the column, so a column whose dynamic range
  ## exceeds 1e8 carried a band wider than its own small entries.
  ## Review construction: an intertemporal-choice delay of 1, 2 and 3
  ## seconds beside ten years. The band was 3.15 s against a
  ## within-condition spread of 2 s, the model was ACCEPTED, and two
  ## data sets differing on 116 of 120 rows gave a bitwise equal
  ## objective with `mu.delay` fitted from 4 rows of 120.
  gi <- rep(1:2, each = 30L)
  first <- c(1L, 31L)
  wide <- c(rep(c(1, 2, 3), each = 10L), rep(3.15e8, 30L))
  expect_identical(gd_varying_groups(wide, gi, first), 1L)
  ## the old rule, written out, so the flip is pinned rather than
  ## remembered: it accepted exactly this column
  ref <- first[gi]
  expect_length(
    which(abs(wide - wide[ref]) > 1e-8 * max(abs(wide))), 0L)
  ## and the answer must not depend on whether the user centered
  expect_identical(gd_varying_groups(wide - mean(wide), gi, first), 1L)

  ## A covariate measured in small units sits beside an intercept of 1
  ## in the same design. A tolerance floored at 1, which is the shape
  ## frmtmb.ode uses, would be 1e-8 absolute and would accept a
  ## covariate whose whole range is 1e-9.
  g2 <- rep(1:2, each = 2L)
  f2 <- c(1L, 3L)
  expect_identical(gd_varying_groups(c(0, 1e-9, 2e-9, 2e-9), g2, f2),
                   1L)
  ## the floor is only a floor: it never exceeds the column, and it is
  ## exactly zero where there is nothing to carry
  expect_gt(gd_col_floor(c(1, 1, 1)), 0)
  expect_lt(gd_col_floor(c(1, 1, 1)), 1)
  expect_identical(gd_col_floor(c(0, 0)), 0)
  expect_identical(gd_col_floor(factor(c("a", "b"))), 0)
})

test_that("every offending variable is named in one message", {
  ## A refusal that stops at the first offender makes a user with three
  ## missing variables pay three round trips, each a whole frame
  ## assembly, and the complete list is what gddm_conditions() takes.
  d <- gc_toy()
  d$w <- rnorm(nrow(d))
  d$blk <- factor(rep(1:4, length.out = nrow(d)))
  e <- expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ x + w, bs ~ blk, ndt ~ 1,
                bias = 0.5), d))
  m <- conditionMessage(e)
  expect_match(m, "the parameters `mu` and `bs` are not constant",
               fixed = TRUE)
  for (v in c("`blk` varies inside condition",
              "`w` varies inside condition",
              "`x` varies inside condition")) {
    expect_match(m, v, fixed = TRUE)
  }
  ## and it hands back the call that fixes it, with every variable in it
  expect_match(m, "gddm_conditions(data, blk, w, x)", fixed = TRUE)
  ## one offender still reads as one
  expect_match(
    conditionMessage(expect_error(
      gc_frame(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                  bias = 0.5), d))),
    "the parameter `mu` is not constant", fixed = TRUE)
})

test_that("gddm_simulate() refuses a parameter it would half-read", {
  ## The same defect in an exported function: gddm_simulate() solves
  ## once per distinct `coh` and reads every parameter at that value's
  ## first trial, so a per-trial parameter vector was accepted and half
  ## of it ignored. On 0.7.0, mu = c(-2.5 x 200, +2.5 x 200) at coh = 0
  ## drew BOTH halves from -2.5, with no error and no warning.
  n <- 40L
  mu <- c(rep(-2.5, n / 2L), rep(2.5, n / 2L))
  ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
  expect_error(
    gddm_simulate(n, mu = mu, bs = 2, ndt = 0.2, coh = 0,
                  control = ctl),
    "`mu` varies between trials that share a coherence", fixed = TRUE)
  ## it names every offending parameter, not just the first
  expect_match(
    conditionMessage(expect_error(
      gddm_simulate(n, mu = mu, bs = rep(c(1, 3), each = n / 2L),
                    ndt = 0.2, coh = 0, control = ctl))),
    "`mu` and `bs` vary between trials", fixed = TRUE)
  ## and the uses that must keep working do
  expect_no_error(
    gddm_simulate(n, mu = 1.5, bs = 2, ndt = 0.2, control = ctl))
  expect_no_error(
    gddm_simulate(n, mu = mu, bs = 2, ndt = 0.2,
                  coh = rep(0:1, each = n / 2L), control = ctl))
})

test_that("a nonlinear formula is covered from both sides", {
  ## Under nl = TRUE the primary parameter has no right-hand side of
  ## its own: its covariates are in the nonlinear BODY, and each
  ## nonlinear parameter carries its own formula. A check that read
  ## only the right-hand sides would pass the first, and one that read
  ## only the body would pass the second.
  d <- gc_toy()
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ a * exp(-b * x), a ~ 1, b ~ 1,
                bs ~ 1, ndt ~ 1, bias = 0.5, nl = TRUE), d),
    "the parameter `mu` is not constant", fixed = TRUE)
  expect_error(
    gc_frame(bf(rt | vint(upper, cond) ~ a * 1, a ~ x, bs ~ 1,
                ndt ~ 1, bias = 0.5, nl = TRUE), d),
    "the parameter `a` is not constant", fixed = TRUE)
  ## and with the covariate in the index the same model is accepted
  d$cond <- gddm_conditions(d, cond, x)
  expect_no_error(
    gc_frame(bf(rt | vint(upper, cond) ~ a * exp(-b * x), a ~ 1, b ~ 1,
                bs ~ 1, ndt ~ 1, bias = 0.5, nl = TRUE), d))
})
