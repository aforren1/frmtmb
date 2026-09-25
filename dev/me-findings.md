# Lane me: brms `me()` noise-free terms

Worktree `/home/user/wt/me`, branch `lane/me`, base `91cbf355`
(frmtmb 0.63.0). Private library `/opt/rlib/lane-me`; base measured
against `/opt/rlib/base`. Machine: Linux cloud container, R 4.5.3,
brms 2.23.0, rstan 2.32.7.

## 1. What brms means, read from brms 2.23.0

Read from the installed namespace (`brms:::me`, `frame_me`,
`get_uni_me`, `data_Xme`, `stan_Xme`, `prior_Xme`, `rename_Xme`,
`stan_sp`, `prepare_predictions_sp`) and from `make_stancode()`:

- `me(x, sdx, gr = NULL)`: `Xn ~ normal(Xme, noise)` for every latent
  value, and `Xme = meanme + sdme * zme`, `zme ~ std_normal()`. One
  latent value per row, or per level of `gr` (then `x` and `sdx` must
  be constant within a level; brms keeps the first value).
- One latent vector per DISTINCT call for the whole model: shared by
  every dpar, nlpar and response that uses it. A variable used in two
  different calls is refused: "Variable 'x' is used in different calls
  to 'me'. Associated calls are: ...".
- Several calls with the same `gr` are correlated by default
  (`Lme`, `lkj(1)`), switched off by `+ set_mecor(FALSE)`. brms's
  `bf(..., mecor = )` does not exist in 2.23.0 (brms reads it as a
  fixed dpar and refuses it).
- Names: coefficient `bsp_<prefix>_<rename(term)>` (`bsp_mexsx`,
  `bsp_mexsx:mezsz`, `bsp_sigma_mexsx`), hyperparameters
  `meanme_<rename(me<x>)>`, `sdme_...`, `corme__mex__mez`
  (`corme_<gr>__...` with a group), latent values `Xme_mex[i]`.
  Prior classes `meanme`, `sdme` (coef `mex`), `corme`; the special
  coefficient is class `b` coef `mexsx`.
- brms's summary does not print meanme/sdme/corme; `variables()` has
  them.
- New data: brms draws `Xme ~ normal(x_new, sdx_new)`. In-sample with
  the default `save_pars(latent = FALSE)` it does the same, with a
  warning; with latent values saved it uses their draws.
- A factor multiplier (`me(x, sx):f`) reaches brms's Csp machinery,
  where brms 2.23.0 emits both columns against the same `Csp_1`
  (`/tmp/lanes/me/me-stancode2.R`, last block). frmtmb refuses factor
  multipliers, as it already does for `mo()` and `mi()`.

## 2. What changed

- `R/me.R` (new): parsing of a term that calls `me()` (expanded with
  `terms()`, so `me(x, sx) * z` is `me(x, sx)`, `z`, `me(x, sx):z`),
  the model-wide registry with brms's duplicate refusal, the frame
  builder (brms's value checks and messages), the placeholder columns,
  `me_latent()` (the density, used on and off the tape), new-data
  values, the natural-scale hyperparameter table, `set_mecor()`
  (exported) and the `?frmtmb-me` page.
- Parameterization: the latent values are appended to the `miss`
  component (so every path that treats `mi()` latent values as inner
  treats these the same way), and three outer components are added:
  `meanme` (natural scale, named `meanme_<coef>`, which is brms's name
  and scale), `logsdme`, `thetame` (the `us_chol_cor()` parameters of
  each correlated group). The latent density is on the latent values
  directly (centered), which is why brms's `zme` differs from them by a
  linear map with log-Jacobian `N * (sum(log(sdme)) +
  sum(log(diag(L))))`.
- `R/parse.R`: the `me()` branch in `parse_linpred()`, the nonlinear
  body refusal, `spec$mecor` (`spec_mecor()`).
- `R/bf.R`: `+ set_mecor()` on `bf()` and `mvbf()` formulas; `mvbf()`
  carries an mv-level setting through a further `+ bf()`.
- `R/frame.R`: model-frame variables, the me frame after the `mi()`
  responses, placeholder columns, template components, `frame$me`.
- `R/objective.R`: `me_latent()` once per evaluation when the model
  has `me()` terms; `lp_eta_fixed()` adds the `me()` columns.
- `R/predict.R`: in-sample columns at the latent modes; new data at
  the observed value.
- `R/brms-names.R`: `bsp_` for `me()` columns; `Xme_<coef>[i]` labels
  (level names with `gr =`) for the latent values of sampled draws.
- `R/confint.R`: `variables()` / `hypothesis()` reach `meanme_`,
  `sdme_`, `corme__`; `anova()` refuses fits whose `me()` calls differ.
- `R/methods-fit.R`: `summary()$me` and a "Noise-free Terms (me())"
  section.
- `R/priors.R`: classes `meanme`, `sdme`, `corme` refused by name;
  brms's default `lkj(1)` corme row dropped (flat).
- `R/fit.R`, `R/importance.R`, `R/sandwich.R`: quadrature, importance
  and `vcov(cluster =)` refuse `me()` fits by name. The sandwich
  guard read `mi_map` only, so without the added condition it would
  have passed `me()` fits silently although its Rd says they are
  refused; the mutant run in section 5 shows it returning a matrix.
- `R/compat.R`: feature `me()` (kind special) and its rules.
- `R/conditional-effects.R`: the noisy variables are plottable.

## 3. Validation

### 3.1 Closed-form marginal likelihood (exact Laplace case)

`dev/me-closed-form.R`, output `dev/me-closed-form.txt`. Gaussian
response, linear `me()` terms: each row (or group) of `(y, x_obs, ...)`
is multivariate normal once the latent values are integrated out,
written independently with `mvtnorm::dmvnorm()`. Evaluated at frmtmb's
estimates (an identity at a point):

    one term               diff -5.684e-14 rel 3.070e-16
    two terms, corme       diff  0.000e+00 rel 0.000e+00
    two terms, mecor FALSE diff  1.137e-13 rel 4.068e-16
    gr = g, 25 levels      diff -1.421e-14 rel 1.648e-16

ML optimum, one term, `optim(BFGS, reltol = 1e-14)` from a perturbed
start on the closed form: maxima differ by 1.54e-09 (optim is higher,
frmtmb's nlminb stops that short), coefficients differ by at most
1.97e-05 standard errors. Seeds 31, 32, 33 (as in
`tests/testthat/test-me.R`).

### 3.2 brms's own log density

`dev/me-brms-lp.R`, output `dev/me-brms-lp.txt`. brms's Stan program
with flat priors, compiled, `rstan::log_prob(adjust_transform = FALSE)`
at frmtmb's estimates and latent modes translated into brms's
parameters (`stan_pars_from_fit()` gained `meanme`, `sdme`, `zme`,
`Lme` rules in `tests/testthat/helper-brms.R`), against frmtmb's joint
density plus the map's log-Jacobian. Seed 23, n = 80.

    shape            diff        rel        max|grad zme|
    cor_interaction  -8.527e-14  3.654e-16  1.091e-11
    mecor_false      -2.842e-14  1.148e-16  3.997e-15
    sigma_dpar        0.000e+00  0.000e+00  3.997e-15
    gr_level         -2.842e-14  1.517e-16  5.329e-15
    poisson          -2.842e-14  1.285e-16  7.697e-12

The zero gradient on brms's `zme` block says the conditional modes are
brms's modes too. The first three shapes are gated test row 23 of
`test-brms-likelihood.R` (run: 9 pass, 0 fail, with
`FRMTMB_BRMS_FIT_TESTS=true`, `/opt/rlib/stan` first on the library
path). Before the reviewer's `/opt/rlib/stan` library existed, the
programs were compiled with `-DTBB_INTERFACE_NEW` through a lane-local
`R_MAKEVARS_USER` (`/tmp/lanes/me/me-Makevars`); nothing in the
repository depends on that.

### 3.3 Joint density written from brms's Stan code (no Stan)

`test-me.R`, "interactions, dpars and non-gaussian families": for
`me(x, sx) * me(z, sz)`, `me(x, sx) * w`, `me(x, sx):w:z`, poisson and
`sigma ~ me(x, sx)`, minus the inner objective at the optimum equals an
R implementation of brms's joint density to a relative 1e-10.

### 3.4 Structure against `make_standata()`

`test-me.R`, last block: `Xn_k`, `noise_k`, `Jme_2` equal brms's
standata, and the special-term coefficient names equal `get_prior()`'s
class `b` coefs.

## 4. Downstream methods

Verified by fits in `/tmp/lanes/me/me-try2.R`, `me-try4.R` and in
`test-me.R`; recorded in `R/compat.R`.

| method | behavior |
|---|---|
| `fixef`, `vcov`, `coef`, `summary` | `mexsx` rows; summary adds `meanme_`, `sdme_`, `corme__` with delta-method SEs |
| `variables`, `hypothesis` | brms names, including the hyperparameters on brms's scale |
| `confint` | internal scale (`meanme_mex`, `logsdme_mex`, `thetame_1`), as for `theta`; `method = "profile"` works |
| `fitted`, `predict`, `residuals`, `simulate` | in-sample at the latent modes; new data at the observed `x`, missing `x` refused |
| `conditional_effects`, `emmeans` | work through new data (observed-value convention) |
| `ranef`, `VarCorr` | unaffected; latent values are not random-effect blocks (brms's `ranef()` has no `Xme` either) |
| `residuals(type = "osa")` | runs |
| `anova` | refuses fits whose `me()` calls differ (new) |
| `prior_summary`, `get_prior`, `set_prior` | class `b` coef `mexsx` honored; `meanme`, `sdme`, `corme` refused by name |
| `frm_sample()` | runs; draws labeled `Xme_mex[i]`, `meanme_mex`, `logsdme_mex` |
| REML, `profile = TRUE`, `autoscale`, `sparse_x`, `weights()` | run; profile, autoscale, sparse_x give the ML fit (logLik differences 3.7e-09, 0, 0) |
| quadrature, importance, `vcov(cluster =)` | refused by name |
| `mvbf`, `rescor`, `mi()` together, nlpar formula, mixture | run |

## 5. Tests

`tests/testthat/test-me.R` (new, 70 expectations, 12 s). Against the
base build: pass=1 err=10 (every block dies on `could not find function
"me"`), the weak form. The behavioral pin: a mutant build with the
`vcov(cluster =)` guard and the `bsp_` naming for `me()` removed
(installed to a throwaway library, sources restored) gives pass=66
fail=4: `vcov(cluster = )` returned a covariance matrix for a `me()`
fit, and the coefficient was named `b_mexsx`.

## 6. Decisions and what is not done

- **Priors on meanme, sdme, corme are refused**, like `simo`: frmtmb
  has no slot for them, and the ML fit has none. Implementing them
  means a natural-scale density with a Jacobian on `logsdme` and an LKJ
  on `thetame`; it is the next step if wanted.
- **New-data predictions use the observed value.** brms draws
  `N(x, sdx)`; the means agree for a linear term, but frmtmb's
  prediction interval does not include the measurement noise. This is
  recorded in `R/compat.R` (`me()` x `predict` conditional), in
  `?frmtmb-me` and in the vignette.
- **In-sample predictions use the latent modes**, which is brms with
  `save_pars(latent = TRUE)`, not brms's default.
- **REML with a `me()` term in mu** runs but is not an exact restricted
  likelihood: the slope and the latent value it multiplies are both
  integrated, and their product is not gaussian. Marked conditional.
  The same argument applies to a `mi(x)` predictor under REML, whose
  compat rule says "works"; not changed here, flagged for review.
- **`AIC()` across fits with different `me()` calls is not refused.**
  `anova()` is; AIC has no comparison hook. Documented in `?frmtmb-me`.
- Refused by name: `me()` in a group-level term, in a nonlinear body,
  inside another call (`I(me(x, sx)^2)`, `s(me(x, sx))`), in one
  interaction with `mo()`/`mi()`, with a factor multiplier, `gr =` that
  is not a bare name, conflicting `set_mecor()` settings between
  responses.
- brms suite: `brm:100` flipped to pass (generator
  `/tmp/lanes/shared/brmsport-gen.R` run; its output for
  `test-brms-suite-brm.R` equals the hand edit, and the only other file
  it touched is the frmtmb.sample copy of `helper-brms-suite.R`). The
  `brmsfit-methods:423/425` rows stay "cannot transfer": they read
  posterior draws or the stanfit, and fixture 3 keeps `Age` for
  `me(Age, AgeSD)` because brms's `Trt` is a factor multiplier; its
  `changed` note now says so. `dev/brmsport-ledger.tsv` is a record of
  an earlier run and still shows the old row.
- `SPEC.md` section 6 lists `me()` as deferred in the historical
  milestone plan; left unchanged, like the other items there that have
  since shipped.

## 7. Speed

`dev/me-timing.R`, outputs `dev/me-timing-base.txt` and
`dev/me-timing-lane.txt`, on a model without `me()` terms (n = 3000,
`(1 + x | g)`, `sigma ~ z`). Value and gradient sum at the optimum are
bitwise identical between base and lane (`f = 4307.0738699271124`,
`sum|grad| = 165.45392255284997`, 17 digits printed), so the tape is
unchanged. Clock, min of 5 whole fits: base 0.505 s, lane 0.465 s, with
the arithmetic control at 0.067 s and 0.047 s in the same processes;
the machine was loaded (1 GB free during the base run), so this says
"no slower" and nothing finer.

On a model with `me()` terms the added per-evaluation work is one
vectorized `dnorm` over the latent values, one `dnorm` or `dmvnorm` for
the latent model and one index per term; no R loop over rows.

## NEWS entry

* `me(x, sdx)` and `me(x, sdx, gr = g)` noise-free predictors, with
  brms's meaning: the latent values are integrated by the Laplace
  approximation (exact for a gaussian response with linear `me()`
  terms), several terms are correlated unless `set_mecor(FALSE)` is
  added, and the names are brms's (`bsp_mexsx`, `meanme_mex`,
  `sdme_mex`, `corme__mex__mez`, `Xme_mex[i]`). The log density agrees
  with brms 2.23.0's Stan program to 1e-13 on five shapes. New data
  use the observed value of the noisy variable. brms's prior classes
  `meanme`, `sdme` and `corme` are refused by name. `vcov(cluster = )`
  now refuses `me()` fits, as its documentation already said, and
  `anova()` refuses fits whose `me()` calls differ.
