# brms default priors against frm(prior =): the measurement

Lane: `wt-brms-priors`, worktree
`C:/Users/adf44/source/r/frmtmb-wt-brms-priors`, branch off `564e185`.
This is the "Follow-on, priors" section of
`dev/brms-likelihood-tests.md` carried out. It decides nothing and
changes no code under `R/`. The maintainer decides; this file is the
numbers the decision needs.

## The headline, before the detail

The plan predicted one finding: that frmtmb puts a standard-deviation
prior on the link scale where brms puts it on the natural scale with a
Jacobian, and that the two therefore disagree by a term that moves the
mode. **The prediction is out of date. Measuring it corrected the
premise and turned up three larger problems in its place.**

The correction first, then the three problems.

1. **`class = "sd"` already agrees with brms.** `R/priors.R:1307`
   evaluates an `"sd"`-scaled entry at `exp(theta)` and adds `theta`
   as the log-Jacobian, which is exactly brms's placement under
   `adjust_transform = TRUE`. The two differ by `log(2)` per
   parameter, the half-t renormalizer brms writes and frmtmb does not,
   and a constant does not move a mode. `man/set_prior.Rd:59` already
   documents the natural-scale placement. Nothing to fix here.

2. **brms's `get_prior()` defaults reach `frm(prior =)` as nothing at
   all.** Every row `get_prior()` writes carries `source == "default"`,
   and `as_priorlist()` (`R/priors.R:460`) drops exactly those. Handing
   a default table to `frm(prior = )` emits one message and fits an
   unpenalized model: `fit$prior` is `NULL`. Measured on all seven
   shapes. This is the finding a user meets first.

3. **The link-scale placement that does survive is on the
   INTERCEPT, and the mismatch there is not a Jacobian.** brms centers
   its design matrix and puts the `Intercept` prior on the intercept at
   the MEAN of the predictors; frmtmb puts it on the intercept at
   `x = 0`. Same density, different argument. For a mean-zero
   predictor the gap is numerically nothing; for `sleepstudy`, where
   `mean(Days) == 4.5`, it is a different prior. On `Reaction ~ Days`
   the same prior string moves the SLOPE by 3.5e-05 in brms and by
   -0.084 in frmtmb, which is 0.069 standard errors. It is the only
   one of these that biases a regression coefficient, and it is the
   largest effect this lane measured.

4. **The one class where the packages genuinely disagree about
   placement is a distributional parameter without a linear predictor**
   (`sigma`, `sigma1`, `shape`, ...). Those are refused by name rather
   than mistranslated, so the disagreement is visible rather than
   silent. But the nearest spelling frmtmb's own error message
   suggests sits on the link scale where brms is on the natural scale,
   and it captures between -0.5% and 35% of what brms's prior intends.
   The natural placement is implemented, unreachable from
   `frm(prior = )`, and reproduces brms's mode to nine figures.

## Environment

- Windows 11, R 4.6.1, rstan 2.32.7, StanHeaders 2.39.1, brms 2.23.0,
  RTMB 1.9, frmtmb 0.50.0 via `pkgload::load_all`.
- rstan cannot compile without the C++17 repair the flat-prior tier
  documents: `R_MAKEVARS_USER` points at a file whose `CXX17FLAGS`
  ends in `-std=gnu++17`. Same fault, same repair, recorded at
  `dev/brms-likelihood-tests.md`.
- Compiled programs under `dev/stan-cache-priors`, keyed the same way
  `helper-brms.R` keys `dev/stan-cache`. Both are in `.gitignore`.

## The shapes

Data and seeds are the flat-prior tier's wherever the tier has the row,
so a number here can be read against a number there.

| id | shape | data | tier row |
| --- | --- | --- | --- |
| S1 | gaussian, `Reaction ~ Days + (1 \| Subject)` | `sleepstudy` | C0's data |
| S2 | `y ~ x + z`, `sigma ~ x` | seed 11, n = 150 | 1 |
| S3 | `Reaction ~ Days + (Days \| Subject)` | `sleepstudy` | C0 |
| S4 | `cumulative(y ~ x)` | seed 5, n = 300 | 12a |
| S5 | nonlinear `y ~ a * exp(-b * x)` | seed 7, n = 120 | 5 |
| S6 | `mixture(gaussian, gaussian)`, `theta1 ~ x` | seed 37, n = 400 | 17 |
| S7 | `Reaction ~ Days`, no random effect | `sleepstudy` | new here |
| S2b | S2 with n = 1500 | seed 11 | 1, resized |
| S5b | S5 with n = 1200 | seed 7 | 5, resized |

And the whole result on one screen. "honored" counts the live
`get_prior()` rows `frm(prior = )` accepts; "lost" counts the rest.
"check A" is `log_prob` minus frmtmb's penalized objective at frmtmb's
estimates, on the set both sides carry.

| id | live rows | honored | lost | check A, AT=FALSE | check A, AT=TRUE | what drives the residual |
| --- | --- | --- | --- | --- | --- | --- |
| S1 | 3 | 2 | 1 (`sigma`) | -2.6971 | +4.3354 | sd Jacobian, log 2, centering |
| S2 | 2 | 2 | 0 | +0.0033 | +0.0033 | centering only |
| S2b | 2 | 2 | 0 | +0.0003 | +0.0003 | centering only |
| S3 | 4 | 3 | 1 (`sigma`) | -3.3956 | +4.8478 | sd Jacobian, log 2, centering, LKJ coordinate |
| S4 | 1 | 0 | 1 (thresholds) | 0 | +0.3985 | nothing on frmtmb's side to compare |
| S5 | 1 | 0 | 1 (`sigma`) | 0 | -2.0032 | nothing on frmtmb's side to compare |
| S6 | 6 | 2 | 4 | 0 | +1.2949 | **exact agreement on the honored set** |
| S7 | 2 | 1 | 1 (`sigma`) | +0.2223 | +4.0820 | centering only |

The `nat` column of the later tables is the one that matters for the
decision: on S5 and S5b it is `log(2)` exactly, and on S7 it leaves
only the centering behind.

## The translation surface

What `frm(prior = )` does with each row brms's `get_prior()` writes.
Measured by handing each row through `frm()` on its own, with `source`
flipped to `"user"` so that the default-drop above is out of the way
and the row's own fate is what shows.

| shape | row | class | frm(prior =) |
| --- | --- | --- | --- |
| S1 | `student_t(3, 288.7, 59.3)` | `Intercept` | honored, on the UNCENTERED intercept |
| S1 | `student_t(3, 0, 59.3)` | `sd` | honored, natural scale with the Jacobian |
| S1 | `student_t(3, 0, 59.3)` | `sigma` | **refused: class** |
| S2 | `student_t(3, 0.8, 2.5)` | `Intercept` | honored, on the UNCENTERED intercept |
| S2 | `student_t(3, 0, 2.5)` | `Intercept`, dpar `sigma` | honored, link scale, as brms |
| S3 | `lkj(1)` | `cor` | honored, on frmtmb's own parameterization |
| S3 | `student_t(3, 288.7, 59.3)` | `Intercept` | honored, on the UNCENTERED intercept |
| S3 | `student_t(3, 0, 59.3)` | `sd` | honored, natural scale with the Jacobian |
| S3 | `student_t(3, 0, 59.3)` | `sigma` | **refused: class** |
| S4 | `student_t(3, 0, 2.5)` | `Intercept` | **refused: no target** |
| S5 | `student_t(3, 0, 2.5)` | `sigma` | **refused: class** |
| S6 | `student_t(3, 0, 2.5)` | `sigma1` | **refused: class** |
| S6 | `student_t(3, 0, 2.5)` | `sigma2` | **refused: class** |
| S6 | `logistic(0, 1)` | `theta2` | **refused: class** |
| S6 | `student_t(3, -0.3, 2.5)` | `Intercept`, dpar `mu1` | honored |
| S6 | `student_t(3, -0.3, 2.5)` | `Intercept`, dpar `mu2` | honored |
| S6 | `logistic(0, 1)` | `Intercept`, dpar `theta1` | **refused: distribution** |

Four distinct refusals, and they are not the same kind of thing.

**refused: class.** `check_brms_prior_class()` (`R/priors.R:507`)
accepts `b`, `Intercept`, `sd`, `cor` and refuses everything else by
name. The message for a distributional parameter is the one this lane
was sent to measure:

> A brms prior with class = "sigma" (student_t(3, 0, 59.3)) has no
> faithful frmtmb spelling. frmtmb's classes are b, Intercept, sd, cor
> and theta. If this names a distributional parameter, its frmtmb
> spelling is class = "Intercept", dpar = "sigma", and that density
> sits on the LINK scale where brms puts it on sigma itself. Write the
> prior you mean with set_prior() directly

The message is accurate and the refusal is the right call. Its one
defect is that it fires for `theta2` as well, and there suggests
`class = "Intercept", dpar = "theta2"`, which does not exist: the
special-cased hint tests `identical(cls, "theta")` and brms spells a
mixture proportion `theta2`. A user following that advice gets
`Prior target not found`.

**refused: no target.** S4's ordinal `Intercept` row. The class is one
of the four the translator accepts, so it passes
`check_brms_prior_class()`, and then `resolve_priorlist()` finds
nothing to attach it to and stops with `Prior target not found
(class=Intercept)`. brms's default prior on ordinal thresholds has no
frmtmb spelling at all: `get_prior()`'s own table offers no slot for
them either.

**refused: distribution.** `parse_prior_dist()` (`R/priors.R:533`)
knows `normal`, `student_t`, `cauchy`, `exponential` and `lkj`.
brms's default on a mixture's `theta` is `logistic(0, 1)`, and its
defaults on `simo` simplexes and on `shape`/`phi` are `dirichlet` and
`gamma`. None of those parse.

**dropped as a default.** Every row above, before any of the four,
when the table is passed whole. This is the one that matters most,
because it is silent apart from one message and the fit succeeds.

## Where each class is actually evaluated

`prior_logdens()` (`R/priors.R:1307`) is the whole story:

```r
prior_logdens <- function(x, dist, scale) {
  jac <- 0
  if (identical(scale, "sd")) {
    jac <- x          # theta = log sd; add the Jacobian
    x <- exp(x)
  }
  ...
  prior_base_logdens(x, dist) + jac
}
```

`resolve_priorlist()` decides which entries get `scale = "sd"`.
Measured by resolving each shape's prior list against its fit and
reading the entries back:

| brms row | frmtmb entry | scale | equals brms? |
| --- | --- | --- | --- |
| `Intercept` (mu) | `beta[1]` | internal | same density, different argument (centering) |
| `Intercept`, dpar with a predictor | `betad[i]` | internal | YES, brms is on the link scale too |
| `sd` | `theta[i]` | **sd** | YES, up to `log(2)` per parameter |
| `cor` | `theta[j..k]` | internal, LKJ | same density on the correlation, different unconstrained parameterization |
| dpar class, nearest spelling | `betad[i]` | internal | NO, brms is on the natural scale |
| dpar class, `natural = TRUE` | `betad[i]` | **sd** | YES, up to `log(2)` per parameter |

`natural = TRUE` is an internal field on a resolved spec. `set_prior()`
has no argument for it and only `frm_sample()`'s default-prior builder
sets it, so **the placement that matches brms for a distributional
parameter exists in the code and is unreachable from `frm(prior = )`.**
This lane sets it by hand to measure it.

## How each shape was measured

Three Stan programs per shape, from the same `get_prior()` table:

- **flat**: every `prior` cell blanked. The baseline.
- **honored**: only the rows `frm(prior = )` accepts.
- **full**: brms's defaults, unchanged.

and three frmtmb fits:

- **hon**: the honored rows, re-spelled through `set_prior()` so that
  the default-drop above does not eat them.
- **link**: hon plus the refused dpar rows in the nearest spelling
  frmtmb's own refusal message names, `class = "Intercept"` with
  `dpar = <class>`, which sits on the link scale.
- **nat**: the same rows with the internal `natural = TRUE` flag set by
  hand, which is brms's placement.

Subtracting the flat program's `log_prob` at the SAME point isolates
each side's prior contribution without assuming the flat-prior tier's
result. frmtmb's own penalty is separated the only way it can be:
`fit$obj` has the prior taped into it and `logLik()` reports the
penalized value, so a second tape of the same model with no prior is
evaluated at the penalized fit's parameters.

## Check A, and the attribution

`frm - stan` is frmtmb's log prior minus the Stan program's, at
frmtmb's estimates. **Every residual is attributed to a named class
with nothing left over: the largest unattributed remainder over the
seven shapes is 1e-8, and most are 0 to machine precision.**

| shape | set | frm log prior | stan log prior | frm - stan | unattributed |
| --- | --- | --- | --- | --- | --- |
| S1 | hon | -7.03864998 | -9.73573525 | +2.69708527 | 0.00000000 |
| S1 | link | -12.12437948 | -14.29928480 | +2.17490532 | 0.00000000 |
| S1 | nat | -8.86354672 | -14.29990914 | +5.43636242 | 0.00000000 |
| S2 | hon | -3.84319726 | -3.83988718 | -0.00331008 | 0.00000000 |
| S2b | hon | -3.83928921 | -3.83903801 | -0.00025120 | 0.00000000 |
| S3 | hon | -11.30284344 | -14.69839383 | +3.39555038 | 0.00000000 |
| S3 | link | -16.38833018 | -19.20875820 | +2.82042802 | 0.00000000 |
| S3 | nat | -13.26517844 | -19.20919179 | +5.94401335 | 0.00000000 |
| S4 | hon | 0 | 0 | 0 | 0 |
| S5 | hon | 0 | 0 | 0 | 0 |
| S5 | link | -2.30452263 | -1.22597847 | -1.07854416 | 0.00000000 |
| S5 | nat | -3.91814319 | -1.22598900 | -2.69215419 | 0.00000000 |
| S5b | link | -2.26672797 | -1.22645415 | -1.04027381 | 0.00000000 |
| S5b | nat | -3.81163547 | -1.22645548 | -2.58518000 | 0.00000000 |
| S6 | hon | -4.79474787 | -4.79474787 | **0.00000000** | 0.00000000 |
| S6 | link | -8.62999146 | -8.98684079 | +0.35684932 | 0.00000000 |
| S6 | nat | -8.88585482 | -8.98657649 | +0.10072167 | 0.00000000 |
| S7 | hon | -5.32456812 | -5.10222613 | -0.22234198 | 0.00000000 |
| S7 | link | -10.41090412 | -9.87946086 | -0.53144326 | 0.00000000 |
| S7 | nat | -6.93502484 | -9.88004414 | +2.94501930 | 0.00000000 |

The attribution, class by class. Each row is one Stan `lprior`
statement against the frmtmb entry that answers it.

| class | frmtmb writes | brms writes | difference | closed form |
| --- | --- | --- | --- | --- |
| `sd` | `st(sd) + log(sd)` | `st(sd) + k log 2` | `sum log(sd) - k log 2` | exact |
| dpar, `natural` | `st(x) + log(x)` | `st(x) + log 2` | `log(x) - log 2` | exact |
| dpar, link | `st(log x)` | `st(x) + log 2` | not a constant, not a Jacobian | - |
| `Intercept` | `st(b0)` | `st(b0 + m'b)` | `st(b0) - st(b0 + m'b)` | exact |
| `cor` = `lkj(e)` | `-log 2 + (e + (d-1)/2) log(1-r^2)` | `-log 2` (d = 2, e = 1) | `(e + (d-1)/2) log(1-r^2)` | exact |

Worked, on S3, which carries all four at once:

| statement | at | brms | frmtmb | difference |
| --- | --- | --- | --- | --- |
| `student_t(3, 288.7, 59.3)` Intercept | 298.995 vs 251.985 | -5.10349146 | -5.32399662 | -0.22050516 |
| `student_t(3, 0, 59.3)` sd (two) | 24.906, 5.985 | -8.90175519 | -5.28363185 | +3.61812333 |
| `lkj(1)` cor | rho = 0.0371157 | -0.69314718 | -0.69521497 | -0.00206779 |
| total | | -14.69839383 | -11.30284344 | +3.39555038 |

`+3.61812333` is `log(24.9060881) + log(5.9852937) - 2 log 2` to eight
figures: **frmtmb's `sd` prior IS brms's, plus the Jacobian brms adds
under `adjust_transform = TRUE`, minus the renormalizer.**
`-0.00206779` is `1.5 * log(1 - 0.0371157^2)`, the Jacobian between
Stan's Cholesky coordinate and frmtmb's, to eight figures. And
`-0.69314718` is `log(1/2)` exactly, because `lkj(1)` on a 2x2
Cholesky factor is a constant.

### Is the residual the sum of the log-Jacobian terms? No

The plan asked exactly that, and the answer is worth stating plainly
because it is nearly yes and the gap is where the findings are. The
residual is

    sum of log-Jacobians for the PRIORED transformed parameters
  - one log(2) per lower-bounded element
  + the intercept centering term
  + the LKJ coordinate change

and the first line is not the same as `jac`, the Jacobian sum Stan
adds under `adjust_transform = TRUE`. `jac` runs over EVERY
constrained parameter, priored or not, because a Jacobian belongs to
the transform rather than to the prior; frmtmb carries a Jacobian only
where a prior asked for one. On S1's honored set the two differ by
`log(sigma)`, since `sigma` is constrained in Stan and unpriored on
both sides: `jac` is 7.0324881 and the priored part is
`log(sd) = 3.6019897`.

So "the residual is the Jacobian sum" is true only for a model in
which every constrained parameter carries a prior, both sides
parameterize it the same way, and no prior sits on a centered
quantity. S5 is the only shape here that meets all three, and there
the residual is exactly `log(2)`.

## Check B, and which density frmtmb maximized

`lpF - pen` and `lpT - pen` are check A under the two settings; `jac`
is the Jacobian sum Stan adds over every constrained parameter,
priored or not.

| shape | set | lpF - pen | lpT - pen | jac | max\|grad\| F | max\|grad\| T | z block |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | hon | -2.69708527 | +4.33540280 | 7.0324881 | 16.81 | 16.81 | 2.1e-14 |
| S1 | link | -2.17490532 | +4.85757873 | 7.0324840 | 17.14 | 16.81 | 1.2e-14 |
| S1 | nat | -5.43636242 | +1.59804932 | 7.0344117 | 17.8 | 16.8 | 1.3e-14 |
| S2 | hon | +0.00331008 | +0.00331008 | 0 | 0.01364 | 0.01364 | - |
| S2b | hon | +0.00025120 | +0.00025120 | 0 | 0.002053 | 0.002053 | - |
| S3 | hon | -3.39555038 | +4.84783462 | 8.2433850 | 28.81 | 27.81 | 1.5e-14 |
| S3 | link | -2.82042802 | +5.42295377 | 8.2433818 | 29.05 | 28.05 | 1.1e-14 |
| S3 | nat | -5.94401335 | +2.30066861 | 8.2446820 | 29.78 | 28.78 | 9.0e-15 |
| S4 | hon | 0.00000000 | +0.39852080 | 0.3985208 | 1.6e-4 | 1 | - |
| S5 | hon | -1.4e-14 | -2.00317470 | -2.0031747 | 1.6e-6 | 1 | - |
| S5 | link | +1.07854416 | -0.92316236 | -2.0017065 | 0.3557 | 0.6443 | - |
| **S5** | **nat** | +2.69215419 | **+0.69314718** | -1.9990070 | 1 | **3.5e-05** | - |
| S5b | link | +1.04027381 | -0.85203259 | -1.8923064 | 0.3438 | 0.6562 | - |
| **S5b** | **nat** | +2.58518000 | **+0.69314718** | -1.8920328 | 1 | **2.9e-05** | - |
| S6 | hon | -0.00000000 | +1.29485696 | 1.2948570 | 3.3e-4 | 1 | - |
| S6 | link | -0.35684932 | +0.93811689 | 1.2949662 | 0.3865 | 1.001 | - |
| S6 | nat | -0.10072167 | +1.20046582 | 1.3011875 | 1 | 1 | - |
| S7 | hon | +0.22234198 | +4.08201467 | 3.8596727 | 0.05557 | 0.9996 | - |
| S7 | link | +0.53144327 | +4.39111205 | 3.8596688 | 0.7025 | 0.2975 | - |
| S7 | nat | -2.94501930 | +0.91547367 | 3.8604930 | 0.9999 | **0.0556** | - |

**S5 is the row that answers the question.** It is the only shape with
no confounder: its one default prior is the half-t on `sigma`, `sigma`
is its only constrained parameter, and brms does not center a nonlinear
predictor, so the intercept-centering difference is absent.

- with no prior at all, `lpF - pen` is 1e-14 and the AT=FALSE gradient
  is 1.6e-06: frmtmb maximizes the density Stan reports with
  `adjust_transform = FALSE`, which is the flat tier's result.
- with the prior in frmtmb's LINK spelling, neither gradient vanishes
  (0.3557 and 0.6443). They differ by exactly 1, the derivative of the
  one log transform, so frmtmb's optimum sits strictly between the two
  Stan optima.
- with the prior on the NATURAL scale, `lpT - pen` is `0.69314718`,
  which is `log(2)` to eight figures, and the AT=TRUE gradient is
  3.5e-05. **frmtmb then maximizes exactly the density brms samples**,
  and the only difference left is the half-t renormalizer, which is a
  constant.

The joint shapes S1 and S3 say the same thing less cleanly, because
the outer gradient of a joint density is not zero at a marginal MAP;
the flat tier's check C says so and asserts only the z block. That
block is 1e-14 in every variant here: **a prior on the outer
parameters leaves the inner problem alone, so the conditional modes
are still exactly Stan's.**

**S7's `nat` row isolates the centering's gradient signature.** With
the sigma prior on the natural scale, every constrained parameter in
the program is priored and carried the same way on both sides, so the
only thing left holding the AT=TRUE gradient off zero is the
`Intercept` argument. It measures 0.0556, and the same 0.0556 appears
in S7's `hon` row, where no sigma prior exists on either side. The
centering is therefore a term of its own, separable from the placement
question and not repaired by settling it.

S4 is worth reading twice: the two settings differ by `0.39852080`,
the Jacobian of the `ordered[nthres]` threshold vector, and the
AT=TRUE gradient is exactly 1. frmtmb carries no prior at all there,
so the whole of brms's default is what separates them.

## The consequence: where the two modes land

Stan's mode is found on the unconstrained scale with the program's own
analytic gradient, under both settings. `adjust_transform = TRUE` is
the density brms's posterior is a posterior of, so its mode is the
target.

**On the two shapes where the question is answerable** (fixed effects,
so Stan's optimum is a mode and not a degenerate joint):

| shape | n | quantity | frmtmb MLE | frmtmb link | frmtmb natural | brms mode | link miss | in SEs | shift captured |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S5 | 120 | sigma | 0.134906316 | 0.135104527 | 0.135469736 | 0.135469737 | 3.65e-04 | 0.042 | 35.2% |
| S5b | 1200 | sigma | 0.150702504 | 0.150723778 | 0.150765019 | 0.150765019 | 4.12e-05 | 0.013 | 34.0% |

- `b_a` and `b_b` agree to nine figures in every column: the placement
  moves the dispersion parameter and nothing else.
- **the natural spelling reproduces brms's mode to nine figures.**
- the link spelling captures about a third of the shift brms's prior
  intends, at both sample sizes. The absolute miss falls roughly like
  1/n (3.65e-04 at n = 120, 4.12e-05 at n = 1200) and is a few
  hundredths of a standard error either way; the FRACTION does not
  fall, because it is a property of the placement and not of n.

**On the shapes with a random effect the question cannot be answered
this way, and that is itself a measurement.** Stan's joint optimum
runs the standard deviations away: on S3 it gives `sd_1 = (110.2,
66.1)` under `adjust_transform = TRUE` against frmtmb's marginal
`(24.91, 5.99)`, and `sd_1 = 130.3` on S1 against `36.67`. Maximizing
over the group-level effects is not maximizing the marginal
likelihood, which is why the flat tier compares densities at a point
rather than optimizers. What CAN be measured there is frmtmb against
itself, which isolates the placement exactly:

| shape | quantity | no prior | link | natural | natural - link | in SEs | link captures |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S1 | sigma | 30.891909 | 30.8917888 | 30.9556683 | +0.0638795 | 0.037 | -0.19% |
| S3 | sigma | 25.5425536 | 25.5424444 | 25.6100129 | +0.0675685 | 0.045 | -0.16% |

`link captures` is the link spelling's shift as a fraction of the
natural spelling's, and on these two shapes it is **negative and
negligible**: the link placement moves sigma the wrong way by about a
ten-thousandth of what brms's prior intends.

The reason is worth stating, because it generalizes. brms's default
scales are derived from the spread of the RESPONSE: `sleepstudy` gets
`student_t(3, 0, 59.3)` on a sigma near 31. On the natural scale that
is a real, if weak, prior. On the link scale the same numbers put a
half-t of scale 59.3 over a `log sigma` of 3.43, which is 0.06 scale
units from the location and therefore flat. **A prior calibrated to a
natural-scale quantity becomes vacuous when it is placed on the link
scale, and brms's defaults are all calibrated that way.** The 35% on
S5 (scale 2.5, log sigma near -2) and the 0.2% on S1 (scale 59.3, log
sigma near 3.4) are the same effect at two calibrations.

**The intercept centering, measured where it is structural.** In S2
and S2b the predictors are mean-zero by construction, so the gap is
sampling noise in `colMeans(X)` and falls like 1/sqrt(n): the residual
is -0.00331 at n = 150 and -0.00025 at n = 1500, and the mode moves by
1.6e-04 and 2.0e-06. On `sleepstudy`, where `mean(Days) = 4.5`, it is
structural and does not fall with n. The measured density difference
is -0.2205 nats on S3, -0.2118 on S1 and -0.2223 on S7.

S7 is `Reaction ~ Days` with the random effect dropped, added for one
reason: it is the only shape where the centering is structural AND
Stan's optimum is a real mode, so the consequence can be measured
rather than argued.

| quantity | no prior | frmtmb, brms's prior | brms | frmtmb - brms | in SEs |
| --- | --- | --- | --- | --- | --- |
| Days | 10.46724855 | 10.38284256 | 10.46728315 | -0.08444059 | **0.069** |
| Intercept (centered) | 298.5079014 | 298.6627025 | 298.461943 | +0.2007595 | 0.057 |
| sigma (natural placement) | 47.4498179 | 47.48875616 | 47.48791992 | +0.00083624 | 0.00034 |

Standard errors are frmtmb's own, on the scale each quantity lives on:
`Days` from `summary()` directly (1.2313), `sigma` by the delta method
from the log scale (47.45 * 0.0527045), and the centered intercept from
`sigma / sqrt(n)`, which is what the variance of a fitted mean at the
predictor mean reduces to in this design (47.45 / sqrt(180) = 3.537).

Read the first row twice. brms's default `Intercept` prior moves the
SLOPE by 3.5e-05, which is nothing, because the intercept it
constrains is the one at the mean of `Days` and that is orthogonal to
the slope by construction. The same prior string, honored by
`frm(prior = )`, moves the slope by -0.0844, because the intercept it
constrains is the one at `Days = 0`, which in this design is strongly
correlated with the slope. **The centering difference is the only one
of the three that biases a regression coefficient rather than a
dispersion parameter**, and at 0.069 standard errors on 180 rows it is
the largest single effect this lane measured.

The sigma row of the same table is the counterpart: under the natural
placement frmtmb's sigma sits 3.4e-04 standard errors from brms's, and
what is left of even that is the centering feeding through the joint
optimization.

## S6, the mixture: three ways to lose a prior in one model

`mixture(gaussian, gaussian)` with `theta1 ~ x` has six live defaults
and keeps two.

Both columns are the statement's value at frmtmb's honored-set
estimates, so they answer "what is this prior worth at frmtmb's
answer" rather than "at whose optimum".

| row | fate | brms | frmtmb |
| --- | --- | --- | --- |
| `student_t(3, -0.3, 2.5)` Intercept mu1 | honored | -1.95356109 | -1.95356109 |
| `student_t(3, -0.3, 2.5)` Intercept mu2 | honored | -2.84118678 | -2.84118678 |
| `student_t(3, 0, 2.5)` sigma1 | refused: class | -1.33270914 | 0 |
| `student_t(3, 0, 2.5)` sigma2 | refused: class | -1.31113713 | 0 |
| `logistic(0, 1)` Intercept theta1 | refused: distribution | -1.54829744 | 0 |
| `logistic(0, 1)` theta2 | refused: class | 0 (brms writes no statement either) | 0 |

The two honored rows agree **exactly**: `frm - stan` is 0.00000000,
because a mixture's `y ~ 1` has no predictor on either mu, so there is
no centering to disagree about. That is the cleanest confirmation that
the `Intercept` gap is centering and nothing else.

The two sigmas behave exactly as row 5's single sigma. Under the
natural spelling, `frm - stan` is `log(sigma1) - log(2)` and
`log(sigma2) - log(2)` to eight figures: -0.66925066 and -0.77773524.

The logistic on `Intercept_theta1` is the third kind of loss and no
placement decision touches it: `parse_prior_dist()` knows five
densities and `logistic` is not one. It is worth -1.548 nats, and on
its own it moves `Intercept_theta1` from brms's 0.80907215 to
frmtmb's 0.81437589, a gap of 5.3e-03 that the sigma placement cannot
close. `theta2` is a slot in `get_prior()`'s table that generates no
Stan statement at all, so refusing it costs nothing.

## Recommendation

The plan asked for one decision, about `sd`. The measurement says that
one needs no decision and puts four others in its place. They are
independent: each can be taken or refused on its own. Every cost below
is a file and a line, from a sweep of `tests/`, `extensions/*/tests/`,
`R/`, `man/`, `vignettes/`, `README.md`, `NEWS.md` and `SPEC.md`.

### D0. The `sd` placement: nothing to decide

frmtmb's `class = "sd"` is already brms's placement, documented at
`man/set_prior.Rd:59`, and now pinned by `test-brms-priors.R`. The
only difference is `log(2)` per parameter, a constant that moves no
mode.

Cost of doing nothing: none. Cost of "fixing" it toward the link
scale: a real regression, and `tests/testthat/test-setprior.R:58-67`
already fails it.

### D1. brms's defaults dropped in silence

`R/priors.R:457-459` states the reason for the drop: "a row brms
filled in itself is brms's default, and frmtmb chooses its own (see
the Default priors section of `frm_sample()`); keeping both would
apply two densities to one parameter."

**That reason is true on the `frm_sample()` path and false on the
`frm()` path**, where frmtmb applies no defaults at all. There is
nothing for a brms default to collide with, so all the drop achieves
on this path is that `frm(prior = get_prior(...))` fits an unpenalized
model.

#### The row the USER wrote is dropped too

The drop is keyed on `source`, and `source` records who BUILT the row,
not who wrote the density in it. brms sets it when `get_prior()`
returns and does not update it when the user edits the `prior` cell in
place, which is the ordinary brms workflow:

```r
gp <- get_prior(...); gp$prior[i] <- "normal(0, 20)"; brm(prior = gp)
```

Measured on S1, editing the `sd` row to `normal(0, 20)`:

| what | value |
| --- | --- |
| the edited row's `source`, in brms's own frame | still `default` |
| brms's `lprior` line for `sd_1` | becomes `normal_lpdf(sd_1 \| 0, 20)` |
| frmtmb's `fit$prior` | `NULL` |
| logLik, edited table | -897.039321503 |
| logLik, no prior at all | -897.039321503 |
| logLik, `set_prior("normal(0, 20)", class = "sd")` | -898.925661682 |

So brms honors the edit and frmtmb drops it, silently, and the message
the user is shown misattributes it:

> Translating a brms prior: dropped 3 row(s) brms had filled in as its
> own defaults.

One of those three rows was written by the user. The user-facing
statement of D1 is therefore not "brms's defaults are ignored" but
**"`frm(prior = )` drops rows by the `source` column, and `source` does
not track authorship after `get_prior()` returns"**, which is a
correctness bug rather than an ergonomics gap.

Nothing in the user's brms experience warns them either: `brm()` and
`brm(prior = get_prior(...))` generate **byte-identical** Stan code
(verified on S1), so passing the default table back is genuinely a
no-op there, and there is no reason to suspect the row matters.

Pinned, as today's wrong behavior and labelled as such, at
`tests/testthat/test-brms-priors.R:82-119`.

| option | what changes | cost |
| --- | --- | --- |
| **D1a (recommended)** honor `source == "default"` on the `frm()` path, keep dropping on the `frm_sample()` path | one condition in `as_priorlist()` (`R/priors.R:460`), which has to be told which path called it; all six call sites are known (`R/fit.R:440`, `R/simulate-new.R:476`, `R/par-template.R:146`, `R/priors.R:1685`, `frmtmb.sample/R/sample.R:778` and `:1217`) | three tests invert: `test-brms-priors.R:54-80`, `test-brms-priors.R:82-119` (the edited row) and `test-prior-compat.R:437-441`, all of which assert the drop. `NEWS.md:512-513` and `R/priors.R:646-666` ("`frm()` ... is flat in every slot until a prior is set") are rewritten. **D1a is only safe after D2**, because a default table carries a `sigma` row that `check_brms_prior_class()` still refuses, so honoring the table would turn a silent no-op into an error |
| D1b keep the drop, raise `message()` to `warning()` | one line at `R/priors.R:482` | the message test at `test-prior-compat.R:440` still passes; the user still hand-copies every row |
| D1c document only | one paragraph | the surprise stays, and the README's "a prior object brms itself built is translated" (`README.md:59-62`) stays misleading for the object a brms user actually has |

### D2. The placement for a distributional parameter

This is the question the plan meant to ask, relocated from `sd` to
`sigma`. The refusal itself is right and should stay. What is wrong is
the advice it gives (`R/priors.R:515-519`):

> its frmtmb spelling is class = "Intercept", dpar = "sigma", and that
> density sits on the LINK scale where brms puts it on sigma itself

Accurate, and a user who follows it gets a prior that captures -0.5%
to 35% of what brms's does.

| option | what changes | cost |
| --- | --- | --- |
| **D2a (recommended)** expose `natural` on `set_prior()`, and route the brms dpar classes to it instead of refusing them | `set_prior()`'s signature (`R/priors.R:249`) and its spec list (`R/priors.R:304-306`, the one line that must gain the field); `check_brms_prior_class()` (`R/priors.R:507-525`) routes instead of stopping. `print.frmtmb_priorlist()` **already renders** `scale=natural` at `R/priors.R:618`, and `resolve_priorlist()` already reads the flag at `R/priors.R:1191` | `test-prior-compat.R:424-426` greps `"LINK scale"` out of the refusal and goes; `frmtmb.sample/tests/.../test-sample-direct.R:220` asserts `expect_null(sg$natural)` and breaks if the new argument defaults to `FALSE` rather than absent; `R/priors.R:1185-1190` ("a `set_prior()` spec never has the field") becomes false; `vignettes/brms-migration.Rmd:63-65`, `SPEC.md:429-431`, `vignettes/inputs.Rmd:187` and `NEWS.md:507-512` all say "refused by name" and are rewritten; `natural_dpar_prior()` (`frmtmb.sample/R/sample.R:524-528`) becomes redundant. **No existing frmtmb spelling changes meaning** |
| D2a' also flip `class = "Intercept"`/`"b"` with a `dpar` to the natural scale | `R/priors.R:1191`'s `sc <-` line | **do not.** `dpar` includes location parameters. `tests/testthat/test-simulate-ergonomics.R:118-119` writes `normal(log(0.6), 1e-9)` on `dpar = "sigma"`, deliberately on the link scale; on the natural scale every draw is negative and `frm_simulate()` rejects them all, which errors at line 120. That is the ONE hard failure. `frmtmb.ddm/tests/testthat/test-surface.R:133` writes `set_prior("normal(0, 0.1)", class = "b", dpar = "mu")` and the change moves the resolved entry from `internal` to `sd` and `mu.cond` from `+0.265` to `-0.330`, a SIGN FLIP that the file did not catch, because it asserted only `abs(pen) < abs(unpen)` and `abs(-0.330) < abs(0.917)` holds. Damage a suite cannot see is worse than damage that stops it, so this strengthens the refusal rather than weakening it. The assertion is now bounded on both sides (`test-surface.R:134-142`), which is worth having whatever is decided about D2 |
| D2b keep the refusal, rewrite the message so it stops offering the link spelling as an equivalent | `R/priors.R:515-519` and the grep at `test-prior-compat.R:426` | free and honest, but the user is left told what they cannot do and not what they can |
| D2c expose `natural`, leave the brms class refused | `set_prior()` only | a ported script still stops on the row, but the user can then write a prior that means what brms means |

D2a is recommended because the placement is already implemented,
already tested for `sd`, already what `frm_sample()`'s own defaults
use (`frmtmb.sample/R/sample.R:561`), and already rendered by
`print()`. Exposing it is a smaller change than the documentation it
lets us delete. D2b is the minimum that is not misleading.

### D3. The intercept centering

The largest measured effect and the hardest to fix. brms's `Intercept`
prior constrains `b0 + colMeans(X)'b`; frmtmb's constrains `b0`. Not a
scale question, and no Jacobian repairs it.

| option | what changes | cost |
| --- | --- | --- |
| D3a make `class = "Intercept"` address the centered intercept | `prior_logdens()` indexes parameter SLOTS (`R/priors.R:1777-1786`); a density on a linear combination of slots is a new kind of entry with its own gradient | structural, and it changes what every existing `class = "Intercept"` prior means. Not costed here |
| **D3b (recommended for now)** document it with the measured size | `R/priors.R:17-20` and `man/set_prior.Rd:56-57` ("Link scale."), `vignettes/brms-migration.Rmd` | none, but a ported brms script keeps a slope bias of 0.069 SE on data like `sleepstudy` |
| D3c add an opt-in `center = TRUE` | the same machinery as D3a | structural, but nothing changes meaning under the user |

D3b is recommended only because D3a is a larger piece of work than
this lane can cost out, and because the gap is invisible today. The
number belongs in the migration vignette: on `Reaction ~ Days`, the
same prior string moves the slope by 3.5e-05 in brms and by -0.084 in
frmtmb.

### D4. The rows with no spelling at all

Feature gaps rather than placement choices, each separate:

- **ordinal thresholds.** `class = "Intercept"` is accepted and finds
  no target (`R/priors.R:1091`), and frmtmb's own `get_prior()` emits
  no threshold row: its five loops (`R/priors.R:742-823`) read design
  columns, RE blocks, autocor, `rescor` and `theta`, and
  `theta_components` is closed to `c("theta", "thetaac", "thetar")`
  (`R/priors.R:957-961`). Every ordinal model brms builds carries this
  default, so it is the likeliest of the four to be met. **One route
  already reaches it and nobody is told about it**: `resolve_priors()`
  enumerates every `par_template` component except `b` and `miss`
  (`R/priors.R:1711-1771`), so
  `frm(prior = list(tau_raw = prior_normal(0, 5)))` resolves, with
  `scale = "internal"`. Whether that is the right scale is a second
  question, because `cumulative` and `sratio` store
  `(tau_1, log increments)` while `cratio` and `acat` store the
  thresholds themselves.
- **densities.** `logistic`, `dirichlet` and `gamma` are all brms
  defaults and none parses (`R/priors.R:548-573`). `logistic` and
  `gamma` are one `switch` arm each; `dirichlet` needs a simplex
  target, which frmtmb has none of.
- **brms's mixture `theta`.** Correctly refused, and **the refusal
  misfires**: the special-cased hint tests `identical(cls, "theta")`
  (`R/priors.R:509`) and brms spells it `theta2`, so the generic dpar
  branch fires and tells the user to write
  `class = "Intercept", dpar = "theta2"`, which then fails with
  `Prior target not found`. A one-line fix to the condition, worth
  taking whatever is decided about D2.
- **`simo`.** No simplex prior of any kind, recorded already at
  `R/families.R:4159-4168`.

### What this lane did not do

Nothing under `R/` was touched. The `natural = TRUE` measurement sets
an internal field on a resolved spec by hand inside the test helper,
precisely so that measuring the placement did not require changing it.

The new test pins today's behavior, and the pins are observations of
the package rather than of the helper: `bp_classify_rows()` puts the
class verdict through `check_brms_prior_class()`, the density verdict
through `parse_prior_dist()`, and the target question through
`as_priorlist()`, so a decision that widens any of the three has to
move them. Under D1 the tests at lines 54 and 82 invert; under D3a the
centering assertions change; and under D2a, measured against a working
build of it, **18 of the 61 assertions fail** (15 failures and 3
errors) across 6 of the 11 blocks, starting with the status vectors at
lines 141 and 163 and running through the `link` branch of the row-5
test at 202-223. That is the intended cost, not collateral damage.

**Corrected, 2026-09-05.** This paragraph previously read "whichever
way each decision goes, `test-brms-priors.R` is where it fails first
and says what changed", and named "the status vectors and the `link`
branch of the row-5 test" as what D2a would move. Both were false as
written. `bp_classify_rows()` decided the class and density verdicts
from its own hardcoded lists (`bp_classes()`, `bp_dists()`) instead of
from frmtmb's, so no D2 decision could move a status: the review
measured the file passing 49 of 49 unchanged against a working D2a
build. The lists are gone and the guarantee above is now a
measurement.

## Timings

Windows 11, R 4.6.1, three R processes competing for one toolchain, so
a compile that costs about 100 s alone cost two to four minutes here.

| step | wall |
| --- | --- |
| S1 (3 programs) | 521 s |
| S2 (2 distinct programs) | 299 s |
| S2b (2) | 425 s |
| S3 (3) | 981 s |
| S4 (2) | 505 s |
| S5 (2) | 325 s |
| S5b (0 new programs) | 21 s |
| S6 (3) | 830 s |
| S7 (3) | 723 s |
| the frmtmb side alone, all shapes, no Stan | 34 s |
| `test-brms-priors.R`, empty cache, 15 programs | 3000 s |
| `test-brms-priors.R`, warm cache | 52 s |
| `test-brms-priors.R`, warm cache, 61 assertions (punch round) | 17 to 27 s over three runs |

The cold figure is the one CI's first run pays, and it was measured
with three other lanes compiling on the same machine; alone it would
be closer to the flat-prior tier's 1482 s for 24 programs. Four of the
six flat programs are byte-identical to ones that tier already
compiles (verified by cache key, not by eye), so a shared cache pays
for them once.

S5b is worth a line: a ten-fold larger data set produced
byte-identical Stan programs, because brms rounds its data-dependent
default scales and both samples gave `student_t(3, 0, 2.5)`. The
content-addressed cache then served every program from S5.

## Deliberate omissions

- **No mode distance for the random-effect shapes.** Stan's optimum
  over the joint density is not the marginal mode, and the numbers say
  so: `sd_1 = (110.2, 66.1)` against frmtmb's `(24.91, 5.99)`. The
  frmtmb-against-itself comparison is reported instead, which isolates
  the placement without the estimand difference.
- **No `simo`, `car`, `gp` or autocorrelation row.** The autocor
  classes already use the natural placement with an exact Jacobian
  (`R/priors.R:1635-1651`), so they sit on the same side of D0 as
  `sd`; measuring them would add compiles and no decision. `simo` has
  no prior of any kind.
- **No multivariate row.** The flat-prior tier's translator has no
  rule for `Lrescor` yet, so the row stops before any prior is
  reached.
- **`natural` is not exposed.** Setting the internal field by hand in
  the helper is a measurement; exposing it is D2a and belongs to the
  maintainer.
- **The `theta2` message defect is reported, not fixed.** It is one
  line under `R/`, and this lane changes nothing under `R/`.
- **No larger data set for S1, S3, S4 or S6.** The density identities
  are constants by construction and do not move with n; the mode
  distance is either degenerate (S1, S3) or absent (S4, which carries
  no frmtmb prior at all). S2b and S5b cover the two shapes where n
  changes the answer.

## Punch round, 2026-09-05

The eight items of `dev/review-brms-priors.md`'s punch list, worked and
measured. Nothing under `R/` changed here either; the one edit outside
this lane's own files is the ddm assertion of item 8.

Environment: private library `.../scratchpad/psp-lib`, frmtmb 0.50.0
installed from this worktree with `frmtmb.ddm` 0.2.0 on top; Stan cache
`.../scratchpad/psp-stan-cache`, seeded by COPYING the lane's 19
programs so that the lane's own cache stayed read-only;
`R_MAKEVARS_USER` at `.../scratchpad/psp-makevars`. One `test_file()`
per process, `FRMTMB_BRMS_FIT_TESTS=true`.

### 1. The helper asks the package (must fix)

`bp_classes()` and `bp_dists()` are gone from
`tests/testthat/helper-brms-priors.R`. In their place:

| verdict | who decides it now | helper line |
| --- | --- | --- |
| `"refused: class"` | `frmtmb:::check_brms_prior_class()` | `:28-31` |
| `"refused: distribution"` | `frmtmb:::parse_prior_dist()` | `:33-36` |
| `"refused: no target"` | `frmtmb:::as_priorlist()` on the row with `source` flipped to `"user"`, then `resolve_prior_input()` | `:44-49`, `:86-91` |

`as_priorlist()` is the path `frm(prior = )` itself takes, so the target
question follows whatever the class gate decides, including a routing
decision that does not exist today. `bp_frm_prior()` (`:105-116`) falls
back to the same translator for a class `set_prior()` has no name for:
no such row exists today, and one exists under any decision that routes
a distributional class instead of refusing it. `bp_shape()`'s dpar
probe reads the `dist_ok` column rather than a list of its own
(`:304`).

Neutral against today's package, checked before anything else: the
status vectors of S1/S3, S2, S4, S5, S6 and S7 are identical to the
pinned ones, and the file runs 61 passed, 0 failed, 0 error warm
(17.2 to 27.3 s over three runs, no new compiles; the spread is other
lanes competing for the machine).

**The proof.** The review's D2a-A shape was rebuilt on a scratch copy
of the core (`.../scratchpad/psp-d2a`, the worktree minus `dev/` and
`extensions/`): `set_prior()` gains `natural = FALSE` and writes the
field only when TRUE, `check_brms_prior_class()` refuses only `theta*`,
and `as_priorlist()` routes a distributional class to
`class = "Intercept"`, `dpar = <class>`, `natural = TRUE`. Installed to
`.../scratchpad/psp-d2a-lib`. Against it the file gives **35 passed, 15
failed, 3 errors: 18 of the 61 assertions move**, across 6 of the 11
blocks:

| line | block | what moved |
| --- | --- | --- |
| 141 | every default row's fate | S1/S3's `sigma` row: `refused: class` becomes `honored` |
| 163 | same | S6's `sigma1` and `sigma2` likewise |
| 202 | row 5 | S5's only default is now honored |
| 203, 204, 207 | row 5 | the honored set carries the sigma prior, so check A's residual is `log(2)` rather than 0 and the AT=FALSE gradient is 118 rather than 1.6e-06 |
| 218, 223 | row 5 | the `link` branch: nothing is refused, so the dpar probe set is empty and `r$link` is `NULL` |
| 256 | row 5, mode | the same cause, at `r$nat` |
| 321, 334 | row C, `sd` | the residual gains the sigma placement |
| 379, 390 | row C, `cor` | the honored program is now the full one, so the half-t count changes |
| 452, 494 | S7 | the status vector and the centering residual |
| 552, 553, 561 | row 17 | the mixture's honored entries are no longer two `beta` rows |

Under the same build the file used to pass 49 of 49, which is the
review's measurement and the reason this was the item it insisted on.

### 2. The two claims that rested on it (must fix)

Rewritten in "What this lane did not do", with the correction stated
rather than quietly applied.

### 3. S6's row count (must fix)

`5 | 2 | 3` was wrong in three places. Fixed in the summary table
(`6 | 2 | 4`), in the S6 section's opening sentence ("six live
defaults"), and in `test-brms-priors.R`'s own comment on the mixture,
which said "three of five" where four of six are refused: `sigma1`,
`sigma2` and `theta2` by class, `Intercept`/`theta1` by density.

### 4. The D2a-prime evidence (should fix)

Both halves re-measured here rather than taken from the review.

| build | `test-simulate-ergonomics.R` | ddm `mu.cond`, no prior | with the prior | entry scale |
| --- | --- | --- | --- | --- |
| control | 41 passed, 0 failed | 0.9165381216 | **+0.2654012893** | `internal` |
| D2a-prime | 38 passed, **1 error at :120** | 0.9165381216 | **-0.3301564534** | `sd` |

One hard failure, not two, and on the ddm surface a silent sign flip:
`abs(-0.330) < abs(0.917)` holds, so the assertion as written passed
under the change. The D2a' cell of the D2 table now says exactly that,
and why it makes the refusal stronger rather than weaker.

### 5. The edited row (should fix)

Reported under D1 as "The row the USER wrote is dropped too", with the
three logLik values and brms's own `lprior` line, and pinned at
`tests/testthat/test-brms-priors.R:82-119`. Every assertion there
records today's wrong behavior and says so, so D1a flips a named
expectation. The `brm()` against `brm(prior = get_prior(...))`
byte-identity was re-checked here and holds.

### 6. The theta2 message (should fix)

`test-brms-priors.R:175-176` now reads the refusal text itself, so the
one-line fix to `identical(cls, "theta")` has to come past an
expectation that names the spelling the message currently invents.

### 7. S7 in the suite (nice to have)

`test-brms-priors.R:439-498`. Cost zero compiles: all three of S7's
programs were already in the cache, verified by key before the test was
written. It pins the shape of the headline number rather than the
number alone.

| quantity | measured here |
| --- | --- |
| slope, no prior | 10.46724855 |
| slope, brms's `Intercept` prior through `frm(prior = )` | 10.38302871 |
| shift | -0.08421984 |
| `SE(Days)` | 1.231295522 |
| shift in SEs | **0.0684** |
| brms's own shift, same prior, AT=TRUE mode against the flat program | 1.21e-05, which is 9.8e-06 SE |

The 0.0684 is pinned with a 2 percent tolerance: the lane reported
0.069 and the review 0.0685 by two other routes, and the spread between
the three is optimizer tolerance. The centering identity is pinned
exactly, `raw + mean(Days) * slope` against the centered intercept, and
the density residual equals the centering term to 1e-10.

### 8. The ddm assertion (nice to have)

`extensions/frmtmb.ddm/tests/testthat/test-surface.R:134-142`, the one
edit outside this lane's files. `abs(pen) < abs(unpen)` is replaced by
the ratio `pen / unpen` bounded below by 0 and above by 1, so a shrink
past zero fails. Verified both ways: 47 passed and 0 failed against the
control core, and against the D2a-prime build 46 passed and **1
failed**, naming `shrink` at -0.36. The assertion it replaces passed
under both.

### Verification

| run | result |
| --- | --- |
| `test-brms-priors.R`, warm, control | **61 passed**, 0 failed, 0 error, 0 skipped, 27.3 s on the last run, no new compiles |
| `test-brms-priors.R`, warm, D2a-A | 35 passed, 15 failed, 3 error |
| `test-message-uniqueness.R` | 6 passed, 0 failed |
| `test-bracket-access.R` | 7 passed, 0 failed |
| ddm `test-surface.R`, control core | 47 passed, 0 failed |
| ddm `test-surface.R`, D2a-prime core | 46 passed, 1 failed |
| `test-simulate-ergonomics.R`, control | 41 passed, 0 failed |
| `test-simulate-ergonomics.R`, D2a-prime | 38 passed, 1 error |

Scope, by `git diff --name-only 564e185` plus untracked: `.gitignore`,
`dev/brms-likelihood-tests.md`,
`extensions/frmtmb.ddm/tests/testthat/test-surface.R`,
`dev/brms-priors-findings.md`, `dev/review-brms-priors.md`,
`tests/testthat/helper-brms-priors.R` and
`tests/testthat/test-brms-priors.R`. Nothing under `R/`. The main
checkout is at `5dfdd84` with an empty `git status`, unmoved by this
round. Nothing was committed.
