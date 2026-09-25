# Prior rows keyed by response, dpar and nlpar: lane `wt-mvprior`

Worktree `C:\Users\adf44\source\r\frmtmb-wt-mvprior`, base `51bfaa4e`
(frmtmb 0.62.0, frmtmb.sample 0.10.0). The BASE arm reads the round's
read-only build `C:/Users/adf44/source/r/rellib-r3`; the LANE arm reads
`C:/Users/adf44/source/r/mvprior-lib` first. brms 2.23.0, R 4.6.1,
RTMB 2.0. Library order and arm switch: `dev/mvprior-prelude.R`.

**Answer first.** brms 2.23.0 keys every `b`, `Intercept`, `sd`,
distributional-parameter and residual-correlation prior row by a
PREFIX: the response in a multivariate model, the `dpar` of each
location where a family has several (categorical, multinomial,
mixture), and the `nlpar` on a nonlinear location. A specification
applies only where its prefix EQUALS the row's, an empty field
included. So a prefix-less prior is refused wherever the model has no
prefix-less row, and a prefix is refused where the row carries none.
frmtmb broadcast such priors to every response, every category or
component, or every nonlinear parameter. In three shapes it gave a
silent wrong answer: class `b` without `nlpar` on `a ~ 1 + z` went to
`a_z` only (and an `lb`-only `b` bounded `a_z`); class `b` on `y ~ 1`
applied nothing while `prior_summary()` listed it; and class `b` on a
`cs()` ordinal model missed the category-specific coefficients (all of
them on `ord ~ cs(x)`, where brms applies it to `bcs`). It now refuses
the prefix-less spellings with a classed `frmtmb_error`, reaches the
`cs()` coefficients as brms does, and `default_prior()` lists the rows
with brms's prefixes. frmtmb applies the same rule to its own families
with several locations (`mixture_mvn()`, `lca()`, `hmm()`, `lba()`,
`rdm()`), where brms has no counterpart. By the user's decisions of
2026-09-24 three more brms refusals are matched and two frmtmb
spellings are kept (below). On 131 specifications on the same data,
base frmtmb agrees with brms on 71 and the lane on 128, the 3 left
being the two spellings the user kept; on 36 `cs()` spellings, base
disagrees on 12 and the lane on 0. Of 324 recorded quantities, 317 are
`identical()` across builds; the 7 that are not are explained below, and
one of them is a real movement: `frm_sample()` draws on a mixture with
`dpar` priors (punch round 2 correction). The ported ledger is
unchanged at 250 of 494.

## The user's decisions, 2026-09-24

Held after the first pass; decided by the user on 2026-09-24.

| spelling | brms 2.23.0 | decision | lane |
| --- | --- | --- | --- |
| class `b` where the predictor has no slope (`y ~ 1`; `bf(y2 ~ 1)` with `resp`) | REFUSE (`b`, `b_y2`) | refuse, as brms | refused; the fit stops, so `prior_summary()` never lists it |
| `rescor` with `resp` | REFUSE (`Lrescor_y1`) | refuse, as brms | refused: "Drop resp" |
| `resp = "y"` on a univariate model, any class | REFUSE (`b_y`, `sd_y`, `ar_y`, `sigma_y`) | refuse, as brms | refused: "resp applies only to a multivariate model" |
| `cor` with `resp` | REFUSE (`L_y1`) | keep, frmtmb extension | accepted; documented in `?set_prior` |
| `Intercept` with `nlpar` | REFUSE (`Intercept_a`, `Intercept_ynl_a`) | keep, frmtmb extension | accepted; documented; the nonlinear refusal names it |

## How to reproduce every number here

| what | script | output under `dev/mvprior-log/` |
| --- | --- | --- |
| the models and specifications both sides run | `dev/mvprior-cases.R` | |
| what brms does | `dev/mvprior-brms-probe.R` | `brms-probe.txt` |
| what frmtmb does, per arm | `MVPRIOR_ARM=base\|lane dev/mvprior-probe.R` | `probe-base.txt`, `probe-lane.txt` |
| verdicts side by side | `dev/mvprior-compare.R` | `compare.txt` |
| brms's defaults on the sample tests' data | `dev/mvprior-brms-sample-defaults.R` | `brms-sample-defaults.txt` |
| sampling defaults per arm | `dev/mvprior-sample-defaults.R` | `sample-defaults-base.txt`, `-lane.txt` |
| brms's rows for a nonlinear sigma | `dev/mvprior-brms-nlf.R` | printed; quoted below |
| the reviewer's base problems | `dev/mvprior-filed.R` | `filed.txt` |
| install, lane library only | `dev/mvprior-install.R core sample` | `install-*.txt` |
| one test file, one process (optional trace) | `dev/mvprior-run-tests.R` | |
| one test file, readable failures | `dev/mvprior-onefile.R` | |
| call sites the suites reach | `MVPRIOR_TRACE=...` with the two drivers below | `detect3-*.log`, `detect3-*-trace.tsv` |
| nothing else moved | `dev/mvprior-bitwise.R`, `dev/mvprior-bitwise-inspect.R` | `bitwise-*.txt`, `bitwise-compare.txt` |
| the release tiers | `dev/mvprior-run-suite.ps1`, `dev/mvprior-run-gated.ps1` | `suite.log`, `gated.log` |
| `R CMD check --as-cran` | `dev/mvprior-run-check.ps1` | `check.log`, `<pkg>-testthat.Rout` |
| the ported ledger | `dev/brmsport-record.sh`, `dev/brmsport-ledger.R`, `dev/brmsport-gen.R` | `brmsport-*.txt`, `brmsport-ledger-before.tsv` |
| class b and cs(), per arm and brms | `dev/mvprior-cs.R` | `cs-brms.txt`, `cs-base.txt`, `cs-lane.txt` |
| every several-location family, per arm | `dev/mvprior-extfam.R` | `extfam-base.txt`, `extfam-lane.txt` |
| frm_sample() draws per arm (the reviewer's probe, rerun) | `dev/mvprior-sample-draws.R`, `-compare.R` | `sample-draws-compare.txt` |
| the library holds this tree | `dev/mvprior-libmatch.R` | `libmatch.txt`, `libmatch-control.txt` |

## brms 2.23.0, per shape

`stancode()` with the prior, and `default_prior()`, on `mvprior_data()`
(seed 2309). REFUSE is always "The following priors do not correspond
to any model parameter: ... Function 'default_prior' might be helpful
to you." Full output: `brms-probe.txt`; frmtmb per arm beside it in
`compare.txt`.

| model | brms REFUSES | brms accepts |
| --- | --- | --- |
| `bf(y1 ~ x + (1 \| g)) + bf(y2 ~ x + (1 \| g))` | `b`, `b coef = x`, `Intercept`, `sigma`, `sd`, bounds-only `b`, all without `resp`; `resp = "nosuch"` | the same with `resp`: `normal_lpdf(b_y1 \| 0, 5)` |
| `bf(y1 ~ x, sigma ~ z) + bf(y2 ~ x)` | `b` and `Intercept` with `dpar = sigma` and no `resp` | `b_sigma_y1`, `Intercept_sigma_y1` |
| gaussian + poisson | `b`, `Intercept`, `sigma` without `resp`; `sigma resp = cnt` | `resp` where the slot exists |
| `mvbind(y1, y2) ~ x + (1 + x \| g)`, rescor | `b`, `Intercept`, `sigma` without `resp`; `rescor resp = y1`; `cor resp = y1` | `rescor`, `cor`, `cor group = g` without `resp` |
| two `student()` | `nu` without `resp` | `nu resp = y1` |
| `ar()`, `ma()`, `cosy()`, `unstr()` in one response | each class without `resp`, and with the other response | on the response that has it |
| `mi()` model | `b`, `b coef = mixm`, `Intercept`, `sigma` without `resp` | `b resp = y` reaches `b_y` and `bsp_y` |
| `bf(cat ~ x)`, categorical | `b`, `b coef = x`, `Intercept` without `dpar`; `dpar = nosuch` | `b dpar = mub`: `b_mub`; `Intercept dpar = mub` |
| `bf(cat ~ x + (1 \| g))` | `sd`, `sd group = g` | `sd dpar = mub`: `sd_1`; rows per `mub`, `muc` |
| `bf(ym ~ x)`, two gaussians mixed | `b`, `b coef = x`, `Intercept`, `sd` without `dpar`; class `sigma` | `b dpar = mu1`, `Intercept dpar = mu2`, `sigma1` |
| `bf(y1 ~ x) + bf(cat ~ x, categorical())` | `b resp = cat`, `Intercept resp = cat`; `b dpar = mub` without `resp` | `b dpar = mub resp = cat`: `b_mub_cat` |
| `ynl ~ a * exp(-b * ax)`, `a ~ 1 + (1 \| g)` | `b`, `b coef = Intercept`, `Intercept`, bounds-only `b`, without `nlpar`; `Intercept nlpar = a` | `b nlpar = a`, `sigma` |
| the same with `a ~ 1 + z` | `b`, bounds-only `b` without `nlpar` | `b coef = z nlpar = a`: `b_a[2]` |
| `bf(y1 ~ x) + bf(ynl ~ ...)`, nonlinear | `b nlpar = a` without `resp`; `b resp = ynl` without `nlpar`; `Intercept nlpar = a` with or without `resp` | `b nlpar = a resp = ynl` |
| `y ~ 1 + (1 \| g), sigma ~ 1` | `b`, bounds-only `b`, `b dpar = sigma`, `sd resp = y` | `Intercept`, `sd` |
| `bf(y1 ~ x) + bf(y2 ~ 1)` | `b resp = y2` | `b resp = y1`, `Intercept resp = y2` |
| univariate `y ~ x`, with or without `ar()` | `b`, `Intercept`, `sigma`, `ar` with `resp = "y"` | every prefix-less class |

`default_prior()` for `bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1,
b ~ 1)` lists the `a` and `b` rows by `nlpar` with an EMPTY `dpar`
(`dev/mvprior-brms-nlf.R`); base frmtmb filled both columns, a row
`set_prior()` refuses to spell.

## What changed

**`R/priors.R`**

- `resp_missing_refusal()` with `resp_keyed_prior_classes` (`b`,
  `Intercept`, `ar`, `ma`, `cosy`, `cortime`; a distributional class is
  stored as `Intercept`). Called first in the resolver loop, before
  `dpar_shape_refusal()`, which on base judged a resp-less `sigma`
  spelling against whichever response it met first and gave a wrong
  reason. Its advice is found by RESOLVING the same specification with
  each response's `resp`, so every call it suggests resolves (punch
  minor 1: reading `prior_table()` missed a `dpar` the table did not list,
  and `Intercept nlpar =`).
- `multi_location_family()` and `lp_prior_dpar()`: a location gets the
  empty `dpar` only when its family has one location. `categorical()`,
  `multinomial()`, `mixture()` and a family with the same `mix` slot (the
  latent-class family) have several, and brms names each. Read by
  `prior_table()`, `block_dpar()`, `block_col_prefix()` (so `sd` rows and
  `sd` scope follow) and `target_coefs()`. `primary_dpars` is no longer
  read as "the location": it also holds every location of such a family
  and the nonlinear parameters REML integrates, and that reading caused
  both the broadcast and the nonlinear bug (punch M1).
- `target_coefs()`: a specification with no `dpar` or `nlpar` reaches
  neither a nonlinear parameter nor one of several locations, and is
  refused naming the `nlpar` or `dpar` values that exist. Class `b` that
  reaches no coefficient is refused ("no population-level slope").
- A `resp` naming no response is refused first, naming the responses. A
  `resp` on a univariate model is refused, and so is `rescor` with a
  `resp`. `sd_spec_reach()` no longer normalizes a univariate `resp`,
  because that input cannot reach it now.
- `dpar_shape_refusal()` treats a location as having a predictor, so
  `Intercept dpar = mub` on `cat ~ 1` stays accepted, as in brms.
- `prior_table()` no longer offers a natural-scale class for a location
  dpar, which a several-location family would otherwise have listed.
- The nonlinear refusal names `Intercept` with `nlpar` as a spelling that
  works, matching `?set_prior` (punch minor 2).
- `?set_prior` and `?default_prior`: the prefix rule for multivariate,
  several-location and nonlinear models; `rescor` takes no `resp`; the
  two frmtmb spellings brms refuses. `?default_prior` said the `rescor`
  row is addressed by `resp`, which was wrong.

**`extensions/frmtmb.sample`**

- `default_priors_for()` and `default_sd_prefixes()` write each
  location's `dpar` through `default_lp_dpar()`, which mirrors
  `lp_prior_dpar()`. Without it `frm_sample()` fails on every categorical
  and mixture model, because its own defaults meet the new refusal. The
  rows now match brms's (`sample-defaults-lane.txt` against
  `brms-sample-defaults.txt`: Intercept and sd per `mub` and `muc` at
  `student_t(3, 0, 2.5)`; per `mu1` and `mu2` at `student_t(3, 1.3, 3.3)`
  and `student_t(3, 0, 3.3)`).
- `announce_default_priors()` labels each line with its `dpar`, `nlpar`
  and `resp` (punch minor 3; it printed `Intercept` with no response on
  a multivariate model).
- `evidence-ratio.R`: `er_slot_spelling()` gives a full `set_prior()`
  call for a flat coefficient, with `resp`, `dpar`, `nlpar` and `coef` as
  the model needs (punch minor 3).

NEWS for both packages under a development-version heading, BREAKING.

## Call sites

Found by running with a trace on `frm_stop()` that logs every firing of
this lane's refusals, and of the `sd` refusal a several-location family
now meets, INCLUDING the ones a test swallows (`dev/mvprior-run-tests.R`).

| file | what | fix |
| --- | --- | --- |
| `tests/testthat/test-prior-compat.R`, "resp picks one response of a multivariate model" | asserted the broadcast: class `b` with no `resp` reaching both responses | asserts brms's refusal |
| `tests/testthat/test-nlf.R`, "post-processing follows the nonlinear parameter a body names" | asserted base's rows for `nlf(sigma ~ a + b * z)`: `a` and `b` in the `dpar` column | asserts brms's layout: in `nlpar`, not in `dpar` |

The second was found by the traced suite FAILING, not by a firing: it
is a table layout, not a refusal. The other firings are refusals the
tests expect (`test-prior-sd-scope.R`, `test-setprior.R`,
`test-prior-compat.R`'s `sd` and `resp = "zzz"` cases) and this lane's
own tests. The trace also showed `resp = "zzz"` reaching the new
no-slope refusal, described as "the location predictor of response zzz
has no slope". An unknown `resp` is now refused first, naming the
responses, and a test pins it.

Static search of every `.Rd`, `.Rmd` and `R/` file in core and all seven
extensions for a prior next to a multivariate, nonlinear, categorical,
multinomial, mixture or latent-class model: no example or vignette
passes a spelling that is now refused. `?frm_simulate` and the frmtmb
vignette use prefix-less priors on univariate models, which brms
accepts.

## Tests seen failing on the base build

| file | base | lane |
| --- | --- | --- |
| `tests/testthat/test-prior-mv-resp.R` (new) | 15 pass, 19 fail, 2 error | 36 pass |
| `tests/testthat/test-prior-location-dpar.R` (new; extended in punch round 2) | 11 pass, 38 fail, 3 error | 55 pass |
| `extensions/frmtmb.latent/tests/testthat/test-prior-dpar.R` (new, round 2) | 2 pass, 8 fail | 10 pass |
| `extensions/frmtmb.eam/tests/testthat/test-prior-dpar.R` (new, round 2) | 1 pass, 6 fail | 7 pass |
| `tests/testthat/test-prior-compat.R` | 193 pass, 1 fail | 194 pass |
| `tests/testthat/test-nlf.R` | 74 pass, 1 fail | 75 pass |
| `extensions/frmtmb.sample/tests/testthat/test-default-priors-brms.R` | 15 pass, 4 fail | 19 pass |
| `extensions/frmtmb.sample/tests/testthat/test-evidence-ratio.R` | 33 pass, 1 fail | 34 pass |

Every erroring assertion on base is the LAST in its block, so nothing
after it goes unrun there. The accepted-spelling blocks pass on both
arms by design: they guard against a refusal reaching further than
brms's. The evidence-ratio block samples (one chain, 400 iterations) and
asserts the advice through `hypothesis()`, not through the helper.

## What moved, and what did not

A fit that used a refused prior does not move: it stops before fitting.

`dev/mvprior-bitwise.R`: 17 designs under ML, `REML = TRUE` and
`control(profile = TRUE)`, saving `vcov()`, `fixef()`, `logLik()`,
`summary()` and `predict()` at seed 7, plus frm_sample()'s resolved
default priors and their negative log prior at the estimate on ten
designs. The univariate designs carry prefix-less priors; the nonlinear
one `nlpar` priors; the multivariate ones priors with `resp`, rescor,
and none; the categorical and mixture ones priors with `dpar`, and none.
**237 of 243 identical.** The 6 that differ:

- `cat_dpar` under ML, REML and profile: base REFUSES its
  `sd dpar = "mub"` prior (brms accepts it, `sd_1`), so base has no
  value; the lane fits.
- `sample defaults cat_dpar`, entries and nlp: the same refusal on base.
- `sample defaults mix_dpar`, entries: the SAME entries in a different
  order (`dev/mvprior-bitwise-inspect.R`: "same entries, order aside:
  TRUE"). The order moved because the defaults now carry `dpar`, so a
  user prior with `dpar` replaces its default instead of stacking over
  it. **Corrected in punch round 2:** this is not "no movement". The
  first pass called the negative log prior bitwise identical; the
  reviewer measured a 1 ulp difference at the estimate on its own
  design, and a reordered sum can differ in the last bit. And the
  draws move: rerunning the reviewer's `r2-sample.R` on the final build
  (`dev/mvprior-sample-draws.R`, `sample-draws-compare.txt`, one chain,
  300 iterations, seed 11) gives `frm_sample()` draws that differ by up
  to **8.741103** on a mixture with `dpar` priors, and `prior_summary()`
  lists 4 rows where base listed 6, because base listed two defaults
  the user's priors override. With no user prior, and with a partial
  one (`mix2_partial`), the draws are identical. NEWS says so.

On `cat_noprior`, `cat_re_noprior` and `mix_noprior` the resolved
defaults are identical, order included, so sampling those models does
not move. Two designs refuse identically on both arms and are recorded
by message: an `lb` on class `b` under REML or profile, and a mixture
under REML or profile.

## The ported ledger

Rebuilt on the final build with `FRMTMB_PORT_ROOT` set to this worktree
and `FRMTMB_PORT_LIB` to the lane library. No stale verdict, no stop.
Bin 1: **250 of 494 before and after**; frmtmb.sample 94 of 165 both
times. The regenerated files are content-identical to HEAD: the scripts
write CRLF on this box, they were normalized to LF, and `git diff` is
empty for all of them (punch minor 8). `git status` may still flag them
from a stale index stat, which is not a change.

## Tiers

These are the punch-round-2 tiers. The lane library was built at 15:10
on 2026-09-24 (frmtmb `R/frmtmb.rdb` 15:10:34, frmtmb.sample 15:10:44,
frmtmb.latent 15:10:48, frmtmb.eam 15:10:56). Every covered source
predates it (`R/priors.R` 15:01, `sample.R` 15:02, `lca.R` and `hmm.R`
15:07). After the eighth library loss (18:41) it was checked function
by function against the tree (`dev/mvprior-libmatch.R`, `libmatch.txt`):
frmtmb 1143 identical, 0 differ; frmtmb.sample 211, 0; frmtmb.latent 62,
0; frmtmb.eam 155, 0. The same script pointed at rellib-r3
(`libmatch-control.txt`) reports the 6 and 4 functions this lane
changed, so it can see a stale install. Counts read off the logs by awk
over the RESULT lines.

```
suite  (suite.log, 18:10)  SUITE ran 283 of 283  pass 16525 fail 0 err 0 skip 176
gated  (gated.log, 15:47)  GATED ran 39 of 39    pass 3283  fail 2 err 0 skip 0
check  (check.log, 20:08)  frmtmb 1 NOTE | frmtmb.sample OK | frmtmb.latent 1 NOTE | frmtmb.eam 1 NOTE
       frmtmb-testthat.Rout          FAIL 0 | SKIP 169 | PASS 10969
       frmtmb.sample-testthat.Rout   FAIL 0 | SKIP 11  | PASS 1815
       frmtmb.latent-testthat.Rout   FAIL 0 | SKIP 4   | PASS 351
       frmtmb.eam-testthat.Rout      FAIL 0 | SKIP 6   | PASS 1648
```

Both tier logs predate the 18:39 shutdown and postdate the build. The
check was rerun after the restore, because its first run was cut at
18:37. Each of the three NOTEs is this machine's environmental "V8
unavailable" on the HTML manual. The check trees were deleted after
their `testthat.Rout` files were copied out. The gated tier skips nothing (every RESULT line reads `skip=0`)
and ran with `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`,
`FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ` and `FRMTMB_STAN_CACHE` set,
as the release driver sets them. Its two failures are
`test-drmtmb-agreement.R`, "skew_normal: same likelihood; frmtmb's
default start stalls", IDENTICAL on the base build (`drmtmb-base.txt`,
`drmtmb-lane.txt`): the block pins the stall 0.62.0 fixed.

The traced runs (`detect*`) are an instrument, not a tier: wrapping
`frm_stop()` changes the call `test-conditions.R` reads back. That file
passes untraced on both arms. The last traced suite (`detect4-*`,
283 files) fires the lane's refusals only in the tests that expect them
(`test-prior-compat.R`, `test-prior-sd-scope.R`, `test-setprior.R`) and
in this lane's own tests, including the new frmtmb.latent and
frmtmb.eam files; no other call site.

A frmtmb.sample test writes `tests/testthat/Rplots.pdf` through an
unguarded device (wt-simnewdata's fix); it is deleted before each check.

## Punch round 2

**cs() (M1).** brms 2.23.0 puts class `b` on the category-specific
coefficients of `cs()` on `acat()`, `sratio()` and `cratio()`
(`normal_lpdf(to_vector(bcs) | ...)`), and `coef = "x"` names the term
`cs(x)` (`bcs[1]`). Measured on all three families and three formulas,
12 spellings each (`dev/mvprior-cs.R`, `cs-brms.txt`). Base missed the
cs() coefficients: on `ord ~ cs(x)` class `b` applied nothing, on
`ord ~ z + cs(x)` it reached `z` only, and `coef = "x"` or `"w"` was
refused; 12 of 36 verdicts disagreed with brms, and 6 of the accepts
were silently partial or empty. The first pass of punch round 1 then
refused `ord ~ cs(x)` with a false "no slope ... brms refuses it too".
The lane now reaches every cs() coefficient under class `b` and one
term's under its `coef`, lists a `b` row per cs() term, and counts cs()
terms as slopes: 0 of 36 disagree (`cs-lane.txt`). A MAP fit with a
class `b` prior on a cs() model moves, by design, because the prior now
applies. Row order differs from brms's (brms sorts; frmtmb keeps design
order), which no resolution reads.

**The rule, applied to every family with several locations (M2, minor
3).** `multi_location()` now reads the response: several non-nonlinear
location dpars, a `mix` slot, or a categorical or multinomial family.
Measured per family on both arms (`dev/mvprior-extfam.R`,
`extfam-base.txt`, `extfam-lane.txt`):

| family | package | locations | base, class b with no dpar | lane |
| --- | --- | --- | --- | --- |
| `lca(K = 3)` | frmtmb.latent | theta1, theta2 | reached `theta1_x theta2_x` | refused; rows per dpar |
| `hmm(K = 2)` | frmtmb.latent | mu1, mu2 | reached `mu1_x mu2_x`; one row, empty dpar | refused; rows per dpar |
| `lba(3)` | frmtmb.eam | v1, v2, v3 | reached every `v*_x`; one row | refused; rows per dpar |
| `rdm(3)` | frmtmb.eam | v1, v2, v3 | reached every `v*_x`; one row | refused; rows per dpar |
| `mixture_mvn(K = 2, D = 2)` | frmtmb | mu1d1 ... mu2d2 | reached all four | refused; rows per dpar |

Checked and NOT changed, because each has one location dpar:
`wiener_gng()` (mu), `gddm()` (mu, the default), and the frmtmb.learn
families (each names one `primary`: alpha, Arew, tau, ...). No other
extension family declares several primary dpars (`grep
primary_dpars extensions/*/R`). Refusals on a family brms does not have
say "That is frmtmb's rule for a family whose location is several
distributional parameters"; "as in brms" is kept for categorical,
multinomial and `mixture()` only. `?lca` and `?hmm` show the `dpar`
spelling; NEWS entries in frmtmb, frmtmb.latent and frmtmb.eam are
BREAKING. frmtmb.sample's `default_lp_dpar()` mirrors the rule.

**Messages (minor 5).** One location is named as one ("the
distributional parameter muyes"). On a multivariate model the resp
advice offers resp plus dpar calls for a response whose location is
several dpars (`dpar = "mub", resp = "cat"`), and a test runs each
suggested call.

**Bitwise, round 2.** The battery now also fits `lca()`, `hmm()`,
`lba()`, `rdm()` and `mixture_mvn()` with `dpar` priors, a cs() model
with no prior, and one with an `Intercept` prior. **317 of 324
identical.** The 7 that differ: the 6 above (`cat_dpar` refused on base,
the `mix_dpar` reorder), and `hmm_dpar ML` summary, whose `call` and
`family` slots differ only in the environments of the closures they
carry: `deparse()` of both is identical, and every numeric quantity of
that fit (vcov, fixef, logLik, predict) is identical.

**Filed (item 6)** in `dev/test-backlog.md`, as the reviewer reported
them.

## The machine crash, 2026-09-24 09:24

The first punch-round tier runs were cut by a machine crash. Their logs
(`detect2-*`) are void as tiers and quoted nowhere as counts. They were
read for call sites only: `test-nlf.R`'s layout assertion came from
there, and a `test-sampling-ported.R` error ("ordinal posterior_epred
has brms's draws x obs x category shape") did not reproduce on a rerun
of that file (213 pass, `sampling-ported-lane.txt`). The lane library
was checked against the sources by mtime and reinstalled from the tree
before every tier quoted below.

## Found and NOT fixed

1. **`lb` on class `b` under `REML = TRUE` or `profile = TRUE`** is
   refused as "Unknown parameter(s) in bounds: x", because the
   coefficient is integrated out. Pre-existing and identical on base;
   the message does not say why.
2. **Five base problems the punch reviewer found**, reproduced and filed
   in `dev/test-backlog.md` ("Filed by wt-mvprior after punch round 1"):
   three `bf()` summed with `+`; `cumulative()` in a multivariate model;
   `me()`; `0 + Intercept`; `student()` with `set_rescor(TRUE)`. brms
   accepts all five (`filed.txt`).
3. **The drmTMB agreement block that pins the skew_normal stall** fails
   on base and lane alike, because 0.62.0 fixed the stall; its own
   comment says to flip it. Not this lane's.
4. **The latent-class family** now needs `dpar` on its location priors,
   because it carries a mixture's `mix` slot. brms has no such family;
   the rule follows the mixture it is built as. The traced suite reaches
   no `frmtmb.latent` test or example that passes such a prior.
