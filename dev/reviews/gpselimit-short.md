# Short check of lane gpselimit, fix round 1, 2026-11-17

Review: `dev/reviews/2026-11-17-gpselimit.md`. Lane record:
`dev/gpselimit-findings.md` sections 0, 1 and 9. Scripts:
`dev/reviews/gpselimit-short-a.R` (sweep, plateau, `by = f`, nl2),
`gpselimit-short-probe.R` (held-parameter cases, doubted stop, order),
`gpselimit-short-cmp.R` (base against lane). Logs in
`dev/reviews/gpselimit-short-out/` (local). Arms: `lane` =
wt-gpselimit-lib over rellib-r22, `base` = rellib-r22. Every `lib:` line
names the arm's frmtmb. Two R processes at most; no package code edited;
no git operation beyond `git status` and `git show`.

## Verdict: SHIP

The four fixes do what the review asked. One LOW gap in the FIX 3 gate
(code-read only, below) can be filed instead of fixed.

## 1. FIX 1, the outward probe: confirmed

`a-lane.txt`, `a-base.txt`, `cmp.txt`:

- The 5 plateau stops (rows 6, 12, 18, 24, 30, started at 3e-5 of the
  spacing): lane gives one `frmtmb_gp_basis` warning with "stopped short
  of its maximum", then the SE warning, and loses `theta_1=flat,
  theta_2=flat`. Base: the plain warning and "boundary". Log-likelihoods
  -40.581497, -31.070729, -51.344372, -87.886060, -20.108328, the same on
  both builds.
- The 13 true white rows (1-5, 9, 10, 14, 21, 22, 25-27): lane gives the
  boundary message, `theta_1=boundary, theta_2=boundary`.
- `identical()` base against lane in estimates, `opt$par`,
  log-likelihood, every SE and the lost reasons: all 30 sweep rows (17
  resolved, 13 white). On the 5 plateau stops only the lost reasons
  differ, as intended.

Can the probe mislead with the rest held at the fit? With the rest held,
each probed point is a valid parameter point, so a gain the probe finds
is a real point of higher likelihood: the held rest can hide a gain, but
it cannot make one. Measured both ways (`probe-lane.txt`):

- True white fits beside another free parameter, rows 1, 3, 21:
  `sigma ~ z`, `sigma ~ x` (the gp covariate) and a second `gp(z, k =
  12)`. 7 of the 9 end white (row 1 `sigma ~ x` and `gp2` resolve). All 7
  are certified, none short. The outward ridge is a cost at t = -2
  (+3.35e-07 to +1.82e-06) and t = -4 (+1.87e-05 to +1.01e-04). Refits
  from lscale 1 and 3 spacings return to the same log-likelihood within
  7e-9 (14 of 14).
- The reverse, plateau starts at 3e-5 on rows 6 and 12 with the same
  three extra parameters: 6 of 6 are called short (lane: warning "stopped
  short", "flat"; base: plain warning, "boundary"). Outward at t = -4:
  -2.45e-05 to -6.79e-05; at t = -8: -0.072 to -0.200. The held rest did
  not hide the gain.

## 2. FIX 2, gp(by = f): confirmed

`a-*.txt`, rows 1, 3, 9, 21, 25 (seeds 9600 + row), `f` alternating:

| row | lane: white theta NaN | base: white theta NaN | lane condition |
|---|---|---|---|
| 1 | 4 of 4 | 2 of 4 | message |
| 3 | 4 of 4 | 2 of 4 | message |
| 9 | 4 of 4 | 2 of 4 | message |
| 21 | 4 of 4 | 0 of 4 | message |
| 25 | 4 of 4 | 2 of 4 | message |

Lane loses all four as "boundary"; the intercept and sigma SEs stay
finite (2 of 2 on every row). Estimates, `opt$par` and log-likelihood are
identical to base on all 5; only SEs and lost reasons differ.

Note, not the lane's: row 1 stops at nlminb code 1 "singular
convergence (7)" with no convergence warning on either build (sub-GP
`f = a` at 8e-20 of the spacing, ridge exactly 0 at t = -8 to 2). The
FIX 3 gate reads the convergence check's verdict, so this fit gets the
message. That is the convergence check's call, not the lane's.

## 3. FIX 3, gp_end_warn() after check_convergence(): confirmed

- The review's nl2 seeds 2 and 4 (b0 -0.8, k 20) now warn: lane order
  `m:frmtmb_nl_restart > w:frmtmb_conv_code > [w:frmtmb_conv_gradient] >
  w:frmtmb_gp_basis` (plain warning), the same as base. On a nonlinear
  body the order did not change: frm() calls `gp_end_warn()` after
  `nl_fit_probed()`, which was already after the convergence check.
- A non-nl doubted stop (rows 3 and 21 from a moved start, `iter.max =
  2`, code 1): certified ridge, and the lane gives the warning, not the
  message. Order: base `gp_basis > conv_code > conv_gradient`, lane
  `conv_code > conv_gradient > gp_basis`.
- What else moved: on a non-nl fit the gp() conditions (nugget, Laplace
  and basis) now follow everything from sdreport() through
  `latres_laplace_report()`: the convergence warnings, separation,
  `frmtmb_converged_at_max`, "Hessian is not positive definite",
  `fit_end_checks()` (family checks, ps() knot span) and the latres
  reports. They still precede `mi_laplace_note()` and the SE check.
- Tests: no test reads condition order near gp(). Grep of every core and
  extension test file that fits a gp() for positional reads (`[[1]]`,
  `[1]`, `head()`, `identical(w, ...)`) of collected conditions, and for
  `expect_snapshot`: none. The handlers in test-gp-start.R,
  test-gp-field.R, test-gp-white-limit.R and frmtmb.sample's
  test-brms-draws-order.R collect by class or by `grepl()`. No `_snaps`
  directory exists.
- LOW gap, code-read only: the gate is the class prefix `frmtmb_conv_`.
  `check_convergence()` returns early, before it sets `conv_classes`,
  when separation is found and the code is not 0 (`frmtmb_separation`
  replaces the code warning there). A certified white field on such a fit
  would get the boundary message beside the separation warning. Not
  constructed. Suggest a backlog line, or set `conv_classes` on that
  return too.

## 4. FIX 4, R/par-template.R: no other class dropped

- `nl_fit_off_saddle()`: `hold_gp` holds only a `frmtmb_gp_basis`
  message. A `frmtmb_gp_basis` warning returns from it and is held by
  `warning = hold`, as before. The held gp message does not match
  `nl_flat_warning()` (text patterns "Standard errors are not
  available|are not identified", class `frmtmb_limit_fit`) or the
  `from_flat` grep, so it does not turn on a ridge restart. It is
  released on both paths that return the first fit, and dropped with the
  rest of `held` where a restart returns `best`.
- `nl_fit_probed()`: the two filters add only `"frmtmb_gp_basis"` to the
  `inherits()` list. Only a condition of that class is newly dropped; the
  only signals of that class are `gp_basis_message()`'s two conditions
  (grep of R/). The replacing fits run through `run(st, FALSE)`,
  unannounced, so frm()'s `gp_end_warn()` speaks for the returned fit
  once.

## 5. Tests

| file | build | result |
|---|---|---|
| test-gp-white-limit.R | lane | pass=42 fail=0 err=0 skip=0 warn=0 |
| test-gp-white-limit.R | rellib-r22 | pass=24 fail=14 err=2 skip=0 warn=5 |
| test-gp-by.R | lane | pass=71 fail=0 err=0 skip=2 warn=0 |
| test-gp-field.R | lane | pass=22 fail=0 err=0 skip=0 warn=0 |

The base failures are behavioral, not only missing symbols: the `by = f`
test fails on `all(is.nan(se[wh]))`, the "boundary" reason and the
message; the plateau test on "stopped short of its maximum"; the
boundary test on "white limit of its basis". (`test-wl-*.txt`,
`test-by-lane.txt`, `test-field-lane.txt`.)
