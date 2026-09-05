# Review: one simplex per mo() term (wt-mo-terms)

Reviewer's independent verification. Lane branch `wt-mo-terms` in
`C:\Users\adf44\source\r\frmtmb-wt-mo-terms`, uncommitted on top of
564e185. Main (`C:\Users\adf44\source\r\frmtmb`) at 5dfdd84, clean, and
untouched by this review.

Method: the worktree installed into a private library
`...\scratchpad\rm-lib`; a `git archive 564e185` extraction installed
into `...\scratchpad\rm-lib-base`; brms 2.23.0, rstan 2.32.7,
StanHeaders 2.39.1 read from the user library.

## 1. Enumeration order: CONFIRMED

brms side reproduced with `make_standata()` and, for the term labels
themselves, `brms:::brmsterms(bf(f))$dpars$mu$sp`. Data: `inc` 0:3
(D=3), `w` 0:2 (D=2), `v` 0:4 (D=4), numeric `z`, `z2`.

| formula | brms sp labels, in order | Jmo |
|---|---|---|
| `y ~ mo(inc) + z` | mo(inc) | 3 |
| `y ~ mo(inc) * z` | mo(inc), mo(inc):z | 3,3 |
| `y ~ z * mo(inc)` | mo(inc), mo(inc):z | 3,3 |
| `y ~ mo(inc) + mo(inc):z` | mo(inc), mo(inc):z | 3,3 |
| `y ~ mo(w) + mo(inc):z` | mo(w), mo(inc):z | 2,3 |
| `y ~ mo(inc) + mo(inc):z + mo(w)` | mo(inc), mo(w), mo(inc):z | 3,2,3 |
| `y ~ mo(inc) * mo(w)` | mo(inc), mo(w), mo(inc):mo(w) | 3,2,3,2 |
| **8th (mine)** `y ~ mo(v):z + mo(inc) * z2 + mo(w)` | mo(inc), mo(w), mo(v):z, mo(inc):z2 | **3,2,4,3** |

Every one of the lane's seven rows reproduces exactly, `Imo` and `Jmo`
both. The eighth is the sharper test: three mo() variables with three
different category counts, an interaction written FIRST, and a `*` in
the middle. brms answers 3,2,4,3, both main effects ahead of both
interactions and, within each block, the written order. The lane's two
rules hold.

frmtmb after the change (`frm(..., dry_run = "frame")`, read from
`$linpreds[["y.mu"]]$mo`):

| formula | frmtmb labels | zeta | D |
|---|---|---|---|
| `y ~ mo(inc) * z` | moinc, moinc:z | 1,2 | 3,3 |
| `y ~ z * mo(inc)` | moinc, moinc:z | 1,2 | 3,3 |
| `y ~ mo(inc) + mo(inc):z` | moinc, moinc:z | 1,2 | 3,3 |
| `y ~ mo(w) + mo(inc):z` | mow, moinc:z | 1,2 | 2,3 |
| `y ~ mo(inc) + mo(inc):z + mo(w)` | moinc, mow, moinc:z | 1,2,3 | 3,2,3 |
| 8th | moinc, mow, mov:z, moinc:z2 | 1,2,3,4 | **3,2,4,3** |

`D` equals `Jmo` elementwise in every case, and the label sequence is
brms's sp sequence in brms's b-coefficient spelling. `mo(inc) * mo(w)`
is refused by frmtmb ("the special on exactly one side of ':' or '*'"),
as claimed.

The pristine 564e185 tree disagrees on four of the seven, and the
divergence is worse than a permutation, because simplexes are ALIASED:

| formula | 564e185 labels | zeta | D |
|---|---|---|---|
| `y ~ mo(inc) * z` | moinc:z, moinc | 1,**1** | 3,3 |
| `y ~ mo(inc) + mo(inc):z` | moinc, moinc:z | 1,**1** | 3,3 |
| `y ~ mo(inc) + mo(inc):z + mo(w)` | moinc, moinc:z, mow | 1,1,2 | 3,3,2 |
| 8th | mov:z, moinc:z2, moinc, mow | 1,2,2,3 | **4,3,3,2** |

against brms's 3,2,4,3 on the eighth. Not one position agrees.

## 2. Identity on the tier's data: CONFIRMED, 0.834911 nats and two parameters

`y ~ mo(inc) * z` on row 3b's data (`set.seed(3)`, n = 300), fitted in
each installed tree:

| | logLik | `df` | zeta components |
|---|---|---|---|
| 564e185 | -419.908892845 | 7 | `zeta1` |
| wt-mo-terms | -419.073982034 | 9 | `zeta1`, `zeta2` |

Gap **0.834910811 nats**, matching the reported 0.834911, and exactly
**two more** free parameters (the second simplex's `D - 1 = 2` softmax
coordinates). AIC moves the other way, 853.818 to 856.148, which is the
honest cost of the two parameters and is what the NEWS entry warns about.

Hand computation of the two-simplex log-likelihood at the lane's own
estimates (two independent `step()` maps, `b_main * m1 + b_int * m2 * z`)
reproduces `logLik(fit)` to **3.4e-13**. Substituting `m1` for `m2` (the
old sharing) at those same estimates gives -421.107837, worse, as it
must be. The lane's model is exactly what it says it is.

Observation, not a defect. On this data the true interaction is zero, so
the second simplex is barely identified and the optimizer runs it to the
boundary: `zeta2` raw (-3.97, 17.43), simplex (2.7e-08, 5.0e-10, 1.0),
with `summary()` printing standard errors of 44075 and 4902 and
`confint()` returning (-86389, 86381). Nothing is broken; brms would show
the same flat direction as a flat posterior. But a user who fits
`mo(x) * z` with no real interaction now gets a printed table with a
meaningless row in it where before there was none. The NEWS says standard
errors move; it does not say they can become uninformative. Worth one
clause.

## 3. The tier rows: CONFIRMED, every count and every constant

Run one row per process (`test_file(desc = ...)` under
`pkgload::load_all`), `FRMTMB_BRMS_FIT_TESTS=true`, a private
`FRMTMB_STAN_CACHE`, `R_MAKEVARS_USER` with `CXX17FLAGS` ending in
`-std=gnu++17`:

| row | pass | fail | error | skip |
|---|---|---|---|---|
| the translator round-trips through Stan's constraints | 11 | 0 | 0 | 0 |
| the simplex and group-level rules round-trip | 30 | 0 | 0 | 0 |
| row 2, `y ~ mo(inc) + z` | 2 | 0 | 0 | 0 |
| row 3, `y ~ mo(inc):z` | 2 | 0 | 0 | 0 |
| row 3b, `y ~ mo(inc) * z` | 6 | 0 | 0 | 0 |
| check C, `y ~ mo(inc):z + (1 or g)` | 3 | 0 | 0 | 0 |

Every count is the lane's. The one-process-per-row rule is real: the
helper documents the stanfit DSO failure at `helper-brms.R:85-95` and I
did not try to defeat it.

The constants are exact, not merely inside tolerance. Measured
independently on a formula of my own (section 4), `lp - logLik` came back
1.38629436112 against `2 * lgamma(3)` = 1.38629436112, a residual of
**-1.1e-13** on a `1e-6 * abs(logLik)` tolerance. Two simplexes, two
`lgamma(3)`s; the single-simplex rows carry one.

## 4. The bsp ordering claim: CONFIRMED, and it was silent

The lane says the pre-change order would have handed brms's `bsp` the
wrong coefficient. Most orderings are caught by accident, because a wrong
`Jmo` gives Stan a wrong-length simplex and it errors. To make the
failure SILENT the two mo() variables must have the same category count,
so I built one:

    d$inc and d$q both 0:3 (D = 3);   y ~ mo(inc):z + mo(q)

brms orders these `mo(q)` (main effect) then `mo(inc):z`, `Jmo = 3,3`.
Both trees fit the SAME model here (two different variables, so no
sharing is involved), which isolates the ordering from the simplex-count
change. The tier's own `stan_pars_from_fit()` (the lane's copy, applied
to both trees) then produces:

| | translated `bsp` | `simo_1` | `simo_2` |
|---|---|---|---|
| lane | 0.638550, 0.576264 | 0.6932, 0.1167, 0.1901 | 0.0547, 0.7607, 0.1846 |
| 564e185 | **0.576264, 0.638550** | **0.0547, 0.7607, 0.1846** | **0.6932, 0.1167, 0.1901** |

Every one swapped, no length mismatch, no error. `brms_lp_check()`:

* lane: `measured_const` 1.38629436112, residual **-1.1e-13**, PASS.
* 564e185 order: **fails by 272.70 nats**.

So the claim is right, and understated: the old order was not just wrong,
it was wrong in a way nothing would have caught.

## 5. test-v18.R's rewrite: CONFIRMED

Both reference NLLs re-implemented from the diff and optimized on the
file's own data (`set.seed(11)`, n = 500):

| reference | max logLik | matches |
|---|---|---|
| OLD, one shared simplex (`p[5:7]` used for both terms) | -537.81025239354 | the 564e185 fit, to **3.3e-11** |
| NEW, two simplexes (`p[5:7]`, `p[8:10]`) | -535.12111601704 | the lane fit, to **7.4e-11** |

The lane's "7e-11" is exact. The old reference WAS the old model, so
rewriting it was mandatory, not cosmetic. `expect_gt(logLik(fit),
-op_shared$value)` is right in direction: the shared model sets
`zeta2 = zeta1` and is a nested submodel, so its optimum can only be
lower. It passes on the lane by **2.689 nats** and fails on 564e185
(margin -3.3e-11, i.e. the two are the same fit there), which is the
correct behavior for an assertion meant to pin the change.

**Defect, minor.** The comment introducing it, `test-v18.R:36-37`, says
the bound holds "on data whose two shapes differ". The data at
`test-v18.R:9` is

    y <- 1 + (1.5 + 0.8 * z) * cz_t[inc + 1] + 0.3 * z + rnorm(n, 0, 0.7)

`cz_t` multiplies both `1.5` and `0.8 * z`, so the main effect's shape
and the interaction's are IDENTICAL by construction. The assertion still
holds, but for the generic reason (two extra free parameters cannot fit
worse), not the stated one. A comment that explains the wrong "why".


## 6. Consumers: CONFIRMED, and the checks discriminate

Fitted on data where the two shapes genuinely differ, so a consumer that
reused one simplex for both terms would be caught:

    sA <- c(0, 0.60, 0.75, 1);  sB <- c(0, 0.10, 0.80, 1)
    y  <- 1 + 1.5 * sA[inc+1] + 0.3 * z + 0.9 * sB[inc+1] * z + N(0, 0.6)

The fit recovers two clearly different ladders,
`c1 = (0, 2.0774, 2.3139, 3)` for `moinc` and
`c2 = (0, 0.5823, 2.2810, 3)` for `moinc:z`. Hand model:
`b0 + bz*z + b_mo*c1[inc+1] + b_moz*c2[inc+1]*z`.

| consumer | result |
|---|---|
| `fitted()` | max abs difference from hand **8.9e-16** |
| `predict(newdata)`, rows presenting categories 3,0,2,1,3,0 | **exactly 0** on all six |
| `predict()` on an ordered-factor `mo(incf)` fit, newdata levels 3,2,1,0 | **exactly 0**; same logLik as the numeric-coded fit |
| `conditional_effects("inc")` | matches hand at all four categories, difference **0**, with `z` held at `mean(z) = 0.1237` so the interaction ladder is actually in play |
| `conditional_effects("z")` | see below |
| `simulate(nsim = 500)` | `max abs(rowMeans - eta)` = 0.104 against `4 * MC se` = 0.107; mean per-row sd 0.5977 vs sigma 0.5981 |

`conditional_effects("z")` is the sharpest of these. It draws a straight
line in `z` with `inc` held at category 2. Regressing that line:

* intercept 2.11416541, which implies a `c1` value of 2.3138865, exactly
  `c1[3]`;
* slope 1.0016062, which implies a `c2` value of 2.2810013, exactly
  `c2[3]`.

Had `ce_lp_vars()` reused the main effect's simplex for the interaction,
the slope would have been 1.0123038. It is not. So the claim that
`R/conditional-effects.R` needed no change because it already reads
`mi$zeta` per entry is verified numerically, not just by inspection.

`simulate()` is likewise discriminating: the same row means sit 1.54
away from the shared-simplex eta and 0.104 from the two-simplex one.

Two cosmetic notes, both pre-existing and neither the lane's doing:
`conditional_effects()` on a numeric `mo()` variable builds a 100-point
CONTINUOUS grid over `[0, 3]` and rounds each grid point to the nearest
category, so the plot is a four-step staircase drawn with 100 points
rather than four points. Values at the categories are exact.

## 7. Names and surfaces: CONFIRMED, with one wrong hint

Component spellings are untouched. `par_template()` returns `beta`,
`betad`, `zeta1`, `zeta2`: the same `zeta<j>` spelling as before, one
more of them. `par_template_names()` still spells the coordinates
`zeta1_1`, `zeta1_2`, `zeta2_1`, `zeta2_2`; `confint(parm = "zeta2_1")`
resolves; and `start =` built from `par_template()` with both `zeta1` and
`zeta2` set is accepted and reaches the same optimum.

`variables()` returns `Intercept | z | moinc | moinc:z |
sigma_Intercept | sigma`, with no `zeta` row. It had none on 564e185
either, so "unchanged in spelling" is true, but only because
`variables()` never listed simplex coordinates. Pre-existing; the change
doubles the number of parameters it does not mention.

`get_prior()` emits classes `b` and `Intercept` only, no `zeta` and no
`simo`, on both trees. Its `b` row order does move with the change:
`moinc` now precedes `moinc:z`, where 564e185 had them the other way.
That is user-visible and is what the NEWS entry describes.

Translation of a brms prior object:

* `brms::get_prior()` output, whose `simo` rows are `dirichlet(1)` with
  `source == "default"`: translates, four default rows dropped with the
  usual message, and the fit reaches -419.073982, the same optimum.
* `brms::prior(dirichlet(2), class = "simo", coef = "moinc1")`, with
  `source == "user"`: refused, by name.

**Defect, minor, pre-existing but now on an advertised path.** The
refusal message (`R/priors.R:507-525`, else branch) is generic and its
hint is wrong for `simo`:

> frmtmb's classes are b, Intercept, sd, cor and theta. If this names a
> distributional parameter, its frmtmb spelling is class = "Intercept",
> dpar = "simo", and that density sits on the LINK scale where brms puts
> it on simo itself.

`simo` is not a distributional parameter and `class = "Intercept", dpar =
"simo"` is not something a caller can write. The sentence also
undercounts the vocabulary: `set_prior()`'s `match.arg` at
`R/priors.R:264-265` takes ten classes (`b`, `Intercept`, `sd`, `cor`,
`theta`, `ar`, `ma`, `cosy`, `cortime`, `rescor`) and the message names
five. The lane's findings say the row is "refused by name, with a message
telling the caller to write what they mean"; the first and last clauses
are true, the middle of the message is misdirection.

One more, not the lane's doing and not worth blocking: a caller who edits
`get_prior()` output in place (`gp$prior[gp$class == "simo"] <-
"dirichlet(2)"`) leaves `source` at `"default"`, and `as_priorlist()`
(`R/priors.R:459-462`) drops the row by source, silently. Reported here
only because it is the natural way to reach the simo path and it does not
reach it.
## 8. mo() in a non-mu dpar: the fit is right, the lane's claim about the translator is NOT

`bf(y ~ x, sigma ~ mo(q)) + gaussian()` fits, and the numbering is what
the lane says: the simplex goes into the ONE global sequence as `zeta1`,
while brms declares it `simo_sigma_1`. With `mo()` in both dpars,
`y.mu`'s takes `zeta1` and `y.sigma`'s takes `zeta2`, one global run.

The lane's omissions bullet then says of the tier's translator: "a model
with mo() in two dpars would collide. Pre-existing, untouched, and now at
least loud: the rule errors instead of silently building a wrong-length
simplex."

**That last clause is false, and I reproduced it.** The new guard
(`helper-brms.R:436-441`) fires only on `is.null(zeta)`, i.e. only when
frmtmb has FEWER zetas than the largest trailing integer brms declares.
The dpar collision is not that case. Running the lane's own
`stan_pars_from_fit()`:

| case | brms declares | translated |
|---|---|---|
| A. `sigma ~ mo(q)` only | `simo_sigma_1` | `zeta1` -- **correct, by accident** (the trailing integer happens to be 1 and frmtmb's only simplex is `zeta1`) |
| B. `y ~ mo(inc)+x`, `sigma ~ mo(q)`, both D=3 | `simo_1`, `simo_sigma_1` | both get `zeta1`. `simo_sigma_1` should be `zeta2` (0.3760, 0.4141, 0.2099); it silently receives mu's (0.4850, 0.2772, 0.2378). **No error, no warning, right length.** |
| C. same but `sigma ~ mo(w)`, D=2 | `simo_1` (len 3), `simo_sigma_1` (len 2) | `simo_sigma_1` receives mu's length-**3** simplex. **A wrong-length simplex, built silently** -- precisely the outcome the guard's own comment says it prevents. |

The guard is still worth having: it is what catches the PRE-change
sharing (brms declares `simo_2`, frmtmb has only `zeta1`), which is a
real regression tripwire. But the comment above it,

> both packages enumerate one simplex per special term, in terms() order,
> so simo_<j> is zeta<j>

is true only for `mu`. For any other dpar brms writes `simo_<dpar>_<j>`
and the trailing-integer rule aliases it onto the mu sequence. Case A
passing is luck, not correctness. The honest fix is three lines: refuse
a `simo` name that carries a dpar suffix, since no tier row covers one.

The FIT is not affected by any of this. This is a test-helper defect and
a findings-document overstatement, not shipped behavior.

## Extra, and the largest finding: `zeta<j>` is NOT `simo_<j>` for every family

The branch's headline naming claim is stated without qualification in
three shipped places:

* `R/fit.R:236-238` / `man/frm.Rd:277-279`: "`zeta<j>` is brms's
  `simo_<j>`, and it belongs to the term brms names `<label><1>`".
* `NEWS.md`, the new development-version entry: "so `zeta<j>` is brms's
  `simo_<j>`".
* `R/frame.R:1899-1903`, the comment above the sorted loop.

It is false for every family that declares `extra_pars`.

`assemble_frame()` initializes `extras <- list()` at `R/frame.R:1126`
and then, at `R/frame.R:1420`, REPLACES it with the family's extra
parameters:

    extras <- resp$family[["extra_pars"]](y[[resp$resp_name]], av)

That runs in the per-response loop that ends at `R/frame.R:1422`, before
Phase 1 opens at `:1432` and long before the mo loop at `:1904`. So by
the time the mo loop evaluates `paste0("zeta", length(extras) + 1L)`,
`extras` is already non-empty, and the numbering starts past the family's
parameters. Measured, `y ~ mo(inc) * z`:

| family | par_template | mo term to zeta |
|---|---|---|
| gaussian | beta, betad, **zeta1**, **zeta2** | moinc to zeta1, moinc:z to zeta2 |
| cumulative | beta, tau_raw, **zeta2**, **zeta3** | moinc to zeta2, moinc:z to zeta3 |
| sratio | beta, tau_raw, **zeta2**, **zeta3** | same |
| cratio | beta, tau_raw, **zeta2**, **zeta3** | same |
| acat | beta, tau_raw, **zeta2**, **zeta3** | same |
| cox | beta, sbhaz_raw, **zeta2**, **zeta3** | same |

brms calls them `simo_1` and `simo_2` in all six. `frm_compat()` lists
`cumulative x mo()` (and sratio, cratio, acat) as conditional, i.e.
supported, so this is reachable by a user following the documentation.
Six families declare `extra_pars`: the four ordinal ones
(`R/families.R:2189, 2436, 2480, 2547`), `fam_cox` (`:3835`, whose
`sbhaz_raw` takes slot 1, measured above) and `mixture_mvn` (`:3292`).
Every one of them shifts the mo() simplex index.

Three things follow.

1. **The numbering itself is pre-existing, not a regression.** 564e185
   also puts the first ordinal simplex at `zeta2` (I ran it). What is new
   is that the branch now DOCUMENTS an equality that does not hold.

2. **The tier's new guard earns its keep here, and the findings do not
   say so.** Running the lane's own `stan_pars_from_fit()` on a fitted
   `bf(yo ~ mo(inc) * z) + cumulative()`:

       brms simo pars : simo_1, simo_2
       frmtmb         : beta, tau_raw, zeta2, zeta3
       TRANSLATOR ERROR: brms declares simo_1 but frmtmb has no zeta1

   Loud, by name, exactly as intended. Note this is the case the guard
   actually catches -- not the two-dpar collision the findings credit it
   with, which it does not catch (section 8).

3. **The sentence needs a scope, or the numbering needs a base.** Either
   is defensible. Qualifying the docs is one line and changes no
   behavior: say the mo() simplexes are numbered in brms's special-term
   order, and that `zeta<j>` is `simo_<j>` when the family contributes no
   parameters of its own, with ordinal thresholds shifting the index.
   Renumbering so mo() simplexes always start at 1 is the better end
   state but it is a separate change with its own compatibility cost
   (`start =`, `par_template()` and `confint()` spellings all move for
   ordinal fits), and it is not what this branch was asked for.

The branch is not wrong to leave the numbering alone. It is wrong to
promise, in `?frm` and in NEWS, an equality it only delivers for families
without extra parameters.

## Extra: is the sort rule right, or right on seven examples?

`mo_terms_in_brms_order()` sorts by interaction order and relies on
`order()` being stable to keep the written order within one order. Every
formula the lane checked has AT MOST ONE mo() interaction, so stability
among order-2 terms is never exercised by the test suite. I built six
more, all against `brms:::brmsterms(...)$dpars$mu$sp` and `Jmo`:

| formula | brms sp order | frmtmb D | match |
|---|---|---|---|
| `y ~ mo(inc):z + mo(inc) * z2` | mo(inc), mo(inc):z, mo(inc):z2 | 3,3,3 | yes |
| `y ~ mo(w):z + z2 * mo(inc) + mo(v)` | mo(inc), mo(v), mo(w):z, mo(inc):z2 | 3,4,2,3 | yes |
| `y ~ mo(inc) * z + mo(w) * z2` | mo(inc), mo(w), mo(inc):z, mo(w):z2 | 3,2,3,2 | yes |
| `y ~ mo(w):z2 + mo(inc):z + mo(v)` | mo(v), mo(w):z2, mo(inc):z | 4,2,3 | yes |
| `y ~ z2 * mo(v) + mo(inc):z + mo(w)` | mo(v), mo(w), mo(v):z2, mo(inc):z | 4,2,4,3 | yes |
| `y ~ mo(v) + mo(w) * z + mo(inc)` | mo(v), mo(w), mo(inc), mo(w):z | 4,2,3,2 | yes |

Rows 2, 3, 5 and 6 have two order-2 terms and so do exercise the
stability the lane is relying on, and rows 2, 5 and 6 also have the
`*`-implied main effect arriving out of written order. The rule holds in
all of them. I could not construct a counterexample. The sort is right,
not lucky.

`y ~ mo(inc) * z + mo(w) * z2` is the cheapest of these to add and it is
the one the current tests are missing: it is the only shape where getting
the within-order tie wrong would change the answer.

## Extra: mo() and mi() in one formula share brms's bsp

The lane's omissions list says mi() interactions are not re-sorted
because "they carry no simplex, so nothing is numbered against brms"
(`dev/mo-terms-findings.md:159-161`, and the same words at
`test-mo-terms.R:155`). That reason does not hold. brms puts mo() and
mi() terms in ONE `bsp` vector, in one terms() order across both:

    y ~ mo(inc) * z + mi(x) * z
    brms bsp order  : mo(inc) | mi(x) | mo(inc):z | mi(x):z
    frmtmb X columns: moinc   | moinc:z | mix:z    | mix

Three of the four positions disagree, all scalars, so the tier's `bsp`
rule (`fe[setdiff(names(fe), known)]`, `helper-brms.R:391-395`) would
hand brms the wrong coefficient in exactly the way section 4 shows for
mo(). It is latent, not live: there is no `mi(` anywhere in
`test-brms-likelihood.R` or `test-brms-agreement.R`, and the branch made
it no worse than 564e185, which had `moinc:z` first as well.

So the OMISSION is fine; the REASON given for it is wrong, and it is
wrong in a way that would mislead whoever writes the first mo()+mi()
tier row. Say instead that mi() terms carry no simplex, that their bsp
position is not yet asserted against brms anywhere, and that sorting
them is the obvious follow-up when a tier row needs it.

## Extra: man/frm.Rd does not match its roxygen source

The claim is "a Monotonic effects section on ?frm (R/fit.R roxygen
only)". `man/frm.Rd` carries a hand edit that roxygen will silently
revert. Roxygenising a full copy of the working tree with roxygen2 8.1.0
(the version DESCRIPTION pins in `Config/roxygen2/version`) changes
exactly one line in the whole package:

    -coordinates in a \verb{zeta<j>} component of \code{coef(fit, "template")}, and
    +coordinates in a \verb{zeta<j>} component of \code{\link[=par_template]{par_template()}}, and

`R/fit.R:232` writes `[par_template()]`; `man/frm.Rd:273` says
`coef(fit, "template")`. The two other `[par_template()]` references in
the same roxygen block (`R/fit.R:63` and `:126`) DO render as links at
`man/frm.Rd:73` and `:140`, which is what makes this a hand edit rather
than a roxygen quirk. No other Rd file differs.

Consequences: the rendered help loses a cross-reference that the rest of
the page has, and the next `document()` run reverts the text without
anyone noticing it was ever different. `R CMD check` does not catch it,
which is why as-cran came back clean.

## 9. Scope

`git diff --name-only 564e185` plus untracked, in the worktree:

    NEWS.md
    R/fit.R                                (roxygen only, all 31 added lines are #')
    R/frame.R                              (the only behavior change)
    R/parse.R                              (comment only, 3 lines)
    man/frm.Rd
    tests/testthat/helper-brms.R           (+7: the simo guard)
    tests/testthat/test-brms-likelihood.R  (row 3b, and the header's exemption paragraph)
    tests/testthat/test-v15.R              (one comment)
    tests/testthat/test-v18.R              (the reference NLL)
    ?? dev/mo-terms-findings.md
    ?? tests/testthat/test-mo-terms.R

That is the claimed scope, with one addition the claim does not name:
`tests/testthat/helper-brms.R` is a third tier file, not one of "the two
named hunks". It is a seven-line guard and it is the right change; it
just is not in the list.

NEWS is under `# frmtmb (development version)` at the top of the file and
its cross-references check out: `NEWS.md:2510` is the v0.18 entry with
the false "(brms convention)" parenthetical, and `NEWS.md:151-156` is the
0.50.0 entry that withdrew the claim but not the behavior. The new entry
describes both correctly.

Main (`C:\Users\adf44\source\r\frmtmb`) is at 5dfdd84 with a clean
working tree, unmoved and untouched by this review.

Loose ends in the findings document itself: it cites
`dev/mo-probe/probe1.R`, `probe2.R` (line 6), `dev/mo-makevars` (line
111) and `dev/mo-stan-cache` (line 112). None of the four exists in the
worktree, tracked, untracked or ignored. Either commit them or drop the
citations; a findings document that points at files nobody can open is
worse than one that says "reproduced with make_standata()".

### The core suite, one file per process

`pkgload::load_all` per file, `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS`
unset, both trees. Totals reconciled from well-formed rows only, with a
conflict check across duplicate runs (zero conflicts on either tree).

| | files | pass | fail | error | skip |
|---|---|---|---|---|---|
| lane, claimed | 105 | 5473 | 0 | 0 | 22 |
| **lane, measured** | **105** | **5473** | **0** | **0** | **22** |
| 564e185, claimed | 104 | 5423 | 0 | 0 | 21 |
| **564e185, measured** | **104** | **5423** | **0** | **0** | **21** |

Exact, every column, on both trees. The claim that exactly three files
moved also holds: diffing the two per-file tables, the ONLY files whose
counts differ are

    test-brms-likelihood.R   lane 0,0,0,19   base 4,0,0,18
    test-mo-terms.R          lane 53,0,0,0   base absent
    test-v18.R               lane 27,0,0,0   base 26,0,0,0

and nothing else in 104 shared files moved by a single expectation.
`test-v15.R` is 41 on both, confirming it changed by comment only.

The four files, measured:

| file | measured | lane's claim |
|---|---|---|
| `test-mo-terms.R` | 53, 0, 0, 0 | 53 passes, new |
| `test-v18.R` | 27, 0, 0, 0 | 26 to 27 |
| `test-v15.R` | 41, 0, 0, 0 | unmoved at 41 |
| `test-brms-likelihood.R` | 0, 0, 0, 19 | 4 to 0 passes, 18 to 19 skips |

`test-brms-likelihood.R` reporting zero passes and nineteen skips in a
plain run is the visible consequence of moving row 3b's guard from
`skip_unless_brms()` to `skip_unless_brms_fit()`, exactly as described.
With the env var set that row is 6 passes (section 3).

### as-cran

`R CMD build` then `R CMD check --as-cran --no-manual` with
`_R_CHECK_CRAN_INCOMING_=false` and quarto's bundled pandoc on PATH,
against a copy of the working tree:

    Status: OK

No ERRORs, no WARNINGs, no NOTEs. The suite runs inside it
(`checking tests ... [21m] OK`) and the vignettes re-build clean, so the
behavior change breaks no vignette that fits a monotonic model. Note
this passes despite the `man/frm.Rd` drift in the punch list, which is
the point of that item: `R CMD check` cannot see it.

## Verdict: GO WITH FIXES

The decision is right and the implementation is right. brms does build
one simplex per special term and does enumerate them in `stats::terms()`
order; I reproduced that on the lane's seven formulas, on an eighth of my
own with three mo() variables at three different category counts, and on
six more built specifically to exercise the tie-breaking the sort depends
on. frmtmb after the change agrees with brms in every one. `R/frame.R` is
the only behavior change, it is small, and the claim that
`R/objective.R`, `R/predict.R` and `R/conditional-effects.R` needed
nothing holds both by inspection (`objective.R:159` and `predict.R:48`
index `mi$zeta` per entry) and numerically.

Nothing on the punch list changes a fitted number. But four items state
something that is not true, and the first of them is shipped
documentation.

### Punch list

**P1. `R/fit.R:236-238`, `man/frm.Rd:277-279`, `NEWS.md` (new entry),
`R/frame.R:1899-1903`: "`zeta<j>` is brms's `simo_<j>`" is false for
six families.** `R/frame.R:1420` replaces `extras` with the family's
`extra_pars` before Phase 1, so for `cumulative`, `sratio`, `cratio` and
`acat` the first monotonic simplex is `zeta2`, not `zeta1`
(`tau_raw` takes slot 1); brms still calls it `simo_1`. Measured on
`y ~ mo(inc) * z`: gaussian gives `zeta1, zeta2`, all four ordinal
families give `zeta2, zeta3`, and so does `cox` (`sbhaz_raw` in slot 1).
`frm_compat()` advertises
`cumulative x mo()` as supported, so a user can reach this by following
the docs. The numbering is pre-existing (564e185 does the same); the
unconditional promise is new. Give the sentence a scope, or renumber as
a separate change. Six families declare `extra_pars` and all six shift
the index: the four ordinal ones, `fam_cox` and `mixture_mvn`.

**P2. `tests/testthat/helper-brms.R:436-439` and
`dev/mo-terms-findings.md:169-171`: the guard is credited with the wrong
save.** The findings say a mo() in two dpars "errors instead of silently
building a wrong-length simplex". It does not: the guard fires only on
`is.null(zeta)`, and in the dpar case `zeta1` exists, so `simo_sigma_1`
silently receives mu's simplex, and with different category counts
silently builds a WRONG-LENGTH one. Reproduced. The case the guard
really catches is P1's ordinal one, where it errors by name ("brms
declares simo_1 but frmtmb has no zeta1") -- worth saying, because it is
the tripwire that makes P1 safe. Fix the sentence in both places, and
either refuse a dpar-suffixed `simo` name in three lines or state
plainly that the trailing-integer rule is mu-only.

**P3. `man/frm.Rd:273` is hand-edited and out of sync with
`R/fit.R:232`.** Roxygen 8.1.0 (the version `Config/roxygen2/version`
pins) regenerates that one line and no other line in the package:
`R/fit.R:232` writes `[par_template()]`, the Rd says
`coef(fit, "template")`. The two other `[par_template()]` references in
the same block DO render as links, at `man/frm.Rd:73` and `:140`. The
next `document()` reverts the text silently. `R CMD check` does not catch
it.

**P4. `tests/testthat/test-v18.R:37` describes data the test does not
have.** It says the bound holds "on data whose two shapes differ";
`test-v18.R:9` builds `y <- 1 + (1.5 + 0.8 * z) * cz_t[inc + 1] + 0.3 *
z`, in which `cz_t` is the shape of BOTH terms. The assertion passes, by
2.689 nats, but for the generic nesting reason, not the stated one. The
same sentence at `test-mo-terms.R:51` IS correct, because
`mo_terms_data()` really does build two shapes. Keep it only where true.

**P5. `dev/mo-terms-findings.md:159-161` and `test-mo-terms.R:155`: the
reason for not sorting mi() is wrong.** "They carry no simplex, so
nothing is numbered against brms" is false: brms puts mo() and mi() terms
in ONE `bsp`, interleaved in terms() order. For
`y ~ mo(inc) * z + mi(x) * z` brms wants
`mo(inc) | mi(x) | mo(inc):z | mi(x):z`; frmtmb reports
`moinc | moinc:z | mix:z | mix`. Three of four positions disagree. The
omission is still right (no `mi(` appears in either brms tier file); the
justification is not.

**P6. `dev/mo-terms-findings.md:6, 111, 112` cite four files that do not
exist**: `dev/mo-probe/probe1.R`, `probe2.R`, `dev/mo-makevars`,
`dev/mo-stan-cache`, absent tracked, untracked and ignored alike. Commit
them or drop the citations.

**P7. The scope line undercounts the tier files.**
`tests/testthat/helper-brms.R` is a third changed tier file, not one of
"the two named hunks". The change is right; name it.

**P8. Missing test, one line.** Every formula in `test-mo-terms.R` has at
most ONE mo() interaction, so `order()`'s stability among order-2 terms,
which the whole sort rests on, is never exercised.
`y ~ mo(inc) * z + mo(w) * z2` covers it: brms gives
`mo(inc), mo(w), mo(inc):z, mo(w):z2`, `Jmo = 3,2,3,2`, and the lane's
frame agrees (I ran it).

**P9. `R/priors.R:507-525`, pre-existing, now on an advertised path.**
The brms `simo` prior IS refused by name, as claimed, but the hint tells
the caller to write `class = "Intercept", dpar = "simo"`, meaningless for
a simplex, and names five of the ten classes `set_prior()` accepts at
`R/priors.R:264-265`. `?frm` now documents the simo refusal, so the
message is reachable by following the docs.

**P10. NEWS understates one user-visible consequence.** When the
interaction is near zero the second simplex is unidentified and ML runs
it to the boundary: on the tier's own row 3b data `summary()` prints
standard errors of 44075 and 4902 for `zeta2` and `confint()` returns
(-86389, 86381). The entry says standard errors move; one clause saying
they can become uninformative would save a bug report.

### The dated dev records

`dev/reviews/2026-09-05-brms-likelihood-identity.md` and
`dev/brms-vignette-audit.md` should stay exactly as they are: both are
dated snapshots, both are internally consistent as of their dates, and
the second review's numbers (-419.908893 npar 7, -419.073982 npar 9) are
the ones I reproduced. `dev/brms-likelihood-tests.md` is the exception
and should gain a pointer: it is not a dated record but a live PLAN with
a "Status:" header, that header still says "the exemption list it creates
has exactly that one entry" while the test file it governs now says the
list is empty, and its own body at `:126-129` already said "today that
list is empty". One line in the Status paragraph, and one at the head of
the `### Divergence` section at `:253` saying it was resolved on
wt-mo-terms, is enough.

### Edits I made to the worktree

Four, all comment or generated-text only, no behavior and no expectation
counts moved. Re-ran the three affected test files afterwards: 27, 53,
and 0/19, identical to before.

1. `tests/testthat/test-v18.R:36-40` -- replaced the comment that said
   the bound holds "on data whose two shapes differ" with one that says
   what this data actually is (both terms share `cz_t`, so the gain is
   the generic two-parameter slack) and points at `test-mo-terms.R` for
   the case where the shapes really differ. P4.
2. `tests/testthat/helper-brms.R:436-447` -- extended the `simo_` comment
   to say the `simo_<j>` = `zeta<j>` rule holds only for mu and only for
   a family with no extra parameters, that ordinal and cox fits start at
   `zeta2` and ARE caught by the guard below, and that a non-mu dpar is
   NOT caught and no row exercises it. P1/P2. The guard's code is
   untouched.
3. `man/frm.Rd:273` -- `\code{coef(fit, "template")}` to
   `\code{\link[=par_template]{par_template()}}`, which is what roxygen
   8.1.0 generates from `R/fit.R:232`. Verified: roxygenising a copy of
   the working tree now changes no Rd file at all. P3.
4. `tests/testthat/test-mo-terms.R:154-160` -- replaced "they carry no
   simplex, so nothing is numbered against brms" with the true reason,
   naming the shared `bsp` and the order brms wants, and saying the mi()
   list is left unsorted because no tier row uses `mi()`. P5.

Everything else on the punch list is left for the lane: P1's shipped
sentences in `R/fit.R`, `man/frm.Rd`, `NEWS.md` and `R/frame.R` are a
scope decision the maintainer should make, not a comment fix; P2's
optional three-line guard changes test behavior; P6 through P10 are the
lane's own findings document, a missing test, and NEWS wording.

# Punch re-check, 2026-09-05

Second pass, after the lane applied the punch list. Same rules and the
same private library. Everything below is checked against the DIFF and
re-run, not against the lane's prose. Main is still 5dfdd84 with a clean
tree and was not touched.

## Every punch item, against the diff

| item | where | verdict |
|---|---|---|
| P1 numbering promise made conditional | `R/fit.R:250-256`, `man/frm.Rd:286-291`, `NEWS.md`, `R/frame.R:1898-1903` | **done** |
| P1 translator resolves by position | `helper-brms.R:200-214`, `:451-471` | **done** |
| P1 ordinal test | `test-mo-terms.R:176-207` | **done** |
| P2 dpar suffix refused by name | `helper-brms.R:451-458`, test at `test-mo-terms.R:209-226` | **done** |
| P3 roxygenise idempotent | whole `man/` | **done** |
| P4 tier plan records closure | `dev/brms-likelihood-tests.md:9-16`, `:260-265` | **done** |
| P5 mi() reason corrected | `test-mo-terms.R:154-160`, findings `:185-190` | **done** |
| P6 stale probe citations removed | `dev/mo-terms-findings.md` | **done** |
| P7 Scope names helper-brms.R | findings `:219-226` | **done** |
| P8 two-interaction formula added | `test-mo-terms.R:120-125` | **done** |
| P9 simo prior branch, ten classes | `R/priors.R:514-526` | **done** |
| P10 NEWS near-zero clause | `NEWS.md` | **done** |

Detail on the three that mattered.

**P1.** The frame's CODE is unchanged -- still
`paste0("zeta", length(extras) + 1L)` -- which is the right call: the
lane documented the offset instead of renumbering, so no ordinal fit's
`start =` or `par_template()` spelling moves. The comment now says so in
as many words: the correspondence is carried by "the POSITION ... not
the zeta number, which starts past whatever the family already put in
`extras`". The roxygen adds a paragraph naming ordinal thresholds and
`cox()`'s baseline hazard and ending "Read the terms in order rather
than parsing the number", NEWS says "What lines up with brms is the
ORDER, not the `zeta<j>` number", and `man/frm.Rd:286-291` carries the
same text. Roxygenising a full copy of the working tree now changes
**no Rd file at all** (it changed exactly one before), so P3 is closed
by construction and the two are consistent.

The translator is genuinely rewritten, not re-commented:
`brms_mo_terms_of(fit, "mu")` returns the frame's mo list for one
predictor and `simo_<j>` resolves through `mo[[j]][["zeta"]]`, with a
range check that names the shortfall. The new `cumulative()` test
asserts `zeta2`/`zeta3` and that `simo_1` is `zeta2`'s simplex and
`simo_2` is `zeta3`'s.

**P2.** The guard is `if (!grepl("^simo_[0-9]+$", nm))`, so any
dpar-suffixed name is refused before the integer is read. I re-ran my
own three-case probe against the new helper:

| case | before | now |
|---|---|---|
| A. `sigma ~ mo(q)` only | translated, correct BY LUCK | refused by name |
| B. mo() in both dpars, same D | silently got mu's simplex | refused by name |
| C. mo() in both dpars, different D | silently built a WRONG-LENGTH simplex | refused by name |

Refusing A as well as B and C is the right call and goes past what I
asked for: A only ever worked because the trailing integer happened to
be 1.

**P4.** Both heads record closure. The Status paragraph now reads "had
exactly that one entry. UPDATE, branch `wt-mo-terms`: that entry is
closed", and the `### Divergence` section opens with "RESOLVED on branch
`wt-mo-terms`" and says the rest is kept as the argument for which side
was right. That is exactly the "dated record plus a pointer" shape I
recommended.

## Runs

A harness fault of MINE first, so the numbers below are trustworthy: my
`tier-one.R` loaded the package with `export_all = FALSE`, which made
row 12 error on `ord_tau_from_raw`, an internal. That is my scaffolding,
not the lane's code. I fixed it and re-ran every tier row.

Core, one file per process, `NOT_CRAN=true`, no `FRMTMB_BRMS_FIT_TESTS`:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| `test-mo-terms.R` | **74** | 0 | 0 | 0 |
| `test-v18.R` | 27 | 0 | 0 | 0 |
| `test-v15.R` | 41 | 0 | 0 | 0 |
| `test-message-uniqueness.R` | 6 | 0 | 0 | 0 |
| `test-bracket-access.R` | 7 | 0 | 0 | 0 |
| `test-ordinal.R` | 12 | 0 | 0 | 0 |
| `test-setprior.R` | 27 | 0 | 0 | 0 |
| `test-prior-compat.R` | 105 | 0 | 0 | 0 |

74 is the lane's claim exactly. The other seven are identical to what I
measured for the same files in the first round, so `test-mo-terms.R` is
the only file that moved and it moved by **+21**.

Tier, `FRMTMB_BRMS_FIT_TESTS=true`, one row per process, corrected
harness:

| row | pass | fail | error |
|---|---|---|---|
| translator round-trips through Stan's constraints | 11 | 0 | 0 |
| simplex and group-level rules round-trip | 30 | 0 | 0 |
| row 2, `y ~ mo(inc) + z` | 2 | 0 | 0 |
| row 3, `y ~ mo(inc):z` | 2 | 0 | 0 |
| row 3b, `y ~ mo(inc) * z` | 6 | 0 | 0 |
| check C, `y ~ mo(inc):z + (1 | g)` | 3 | 0 | 0 |
| row 12, ordinal families | 10 | 0 | 0 |

Every count unchanged by the translator rewrite, and row 12 is green
once my harness stopped hiding internals.

## The ordinal question: where the check lives, and whether it is enough

**There is no tier row that combines an ordinal family with `mo()`.**
Row 12 is `y ~ x` across cumulative/sratio/cratio/acat with no monotonic
term; grepping `mo(` in `test-brms-likelihood.R` returns only rows 2, 3,
3b and check C, all gaussian. The ordinal-with-`mo()` check lives ONLY in
`test-mo-terms.R:176-207`.

That test is STRUCTURAL. It builds the Stan code and data, calls
`stan_pars_from_fit()`, and asserts the mapping and the simplex lengths.
It never calls `brms_lp_check()`, so it never reaches Stan's `log_prob`.
It pins the thing that was broken -- the `simo_<j>` to `zeta<j+1>`
mapping -- and it would fail loudly if that regressed. What it cannot
see is whether the COMPOSITION is right: ordinal threshold centering and
monotonic placeholder columns are each pinned against Stan separately
(row 12, and rows 2/3/3b), never together.

So I ran the missing check myself, on that test's own data and model:

    brms_lp_check(bf(yo ~ mo(inc) * z), cumulative(), dd, fit,
                  const = 2 * lgamma(3))

    Stan lp        -262.693761971
    frmtmb logLik  -264.080056332
    measured const    1.38629436112   (2 * lgamma(3) = 1.38629436112)
    RESIDUAL          3.77476e-15
    max |grad|        2.52878e-06

**The identity holds.** So the answer to "is it enough" is: enough for
correctness today and enough to catch a regression in the mapping, but
NOT enough as a guard on the composition, and the composition is exactly
what a future change to ordinal centering or to the monotonic column
placeholder would break. Promoting it to a tier row costs one line and
no compile time that has not already been paid -- the program is in the
cache now. I would take that trade; it is the only ordinal+`mo()`
log-density evidence in existence and right now it lives in my scratch
directory rather than in the suite.

## What is stale in dev/mo-terms-findings.md

The document was partly updated (it describes `brms_mo_terms_of()` and
the `cumulative()` test at `:127-131`), but five things predate the
punch round:

1. **`:237-240`, the core-suite table.** Still says the branch is 105
   files / **5473** passes / 22 skips. `test-mo-terms.R` went 53 to 74
   and the other seven files I re-measured are unmoved, so the branch
   total is **5494**. Off by 21.
2. **`:242`.** "`test-mo-terms.R`, new, 53 passes" -- now 74, and the
   sentence does not mention the two tests the punch round added.
3. **`:213-216`.** Says `dev/brms-likelihood-tests.md` and the other two
   records "still describe the divergence as live" and "were not in the
   change list". The punch round edited that file (P4), so both halves
   are now false for it. The other two records are still correct.
4. **`:174-181`, as-cran.** "Status: OK ... [34m]" was measured before
   this round. `R/priors.R`, `R/fit.R`, `man/frm.Rd`, `NEWS.md`,
   `helper-brms.R` and two test files have changed since, and it has not
   been re-run. See the residual item below.
5. **`:147-154`, the tier table.** Predates the translator rewrite. The
   counts happen to still be right (I re-measured all six), but there is
   no row for the ordinal translation the round added, and the sentence
   at **`:156`**, "The last row is the flip", points at check C when the
   flip is row 3b, the fifth row. That mis-reference is older than this
   round but still wrong.

## Updated verdict: GO

The three substantive findings from the first pass are closed, and
closed properly rather than papered over: the numbering promise is
conditional and truthful, the translator maps by position and refuses
what it cannot map, and the mi() reasoning is corrected. Every count I
was given reproduces. Nothing on the list below changes a fitted number
or blocks a merge.

### Remaining items

* **R1. `dev/mo-terms-findings.md:237-240, :242, :213-216, :174-181,
  :147-156`** -- the five stale passages above. The core-suite total
  (5473, should be 5494) is the one a reader would act on.
* **R2. `tests/testthat/test-brms-likelihood.R`** -- no ordinal+`mo()`
  row. The structural check at `test-mo-terms.R:176-207` is good but
  stops short of `log_prob`; I verified the identity holds at a residual
  of 3.8e-15, so promoting it is cheap and safe.
* **R3. as-cran not re-run this round.** Low risk: roxygenise is
  idempotent, `man/par_template.Rd` exists so the `\link[=par_template]`
  target resolves, and the same link form already passed as-cran twice
  at 564e185. But the claim at findings `:174-181` is currently
  unverified for the tree as it now stands, and one `R CMD check` would
  settle it.
* **R4. `dev/mo-terms-findings.md:156`** -- "The last row is the flip"
  names the wrong row. Pre-existing, one word to fix.

No edits by me this round. The four comment fixes I made in the first
pass all survived the lane's rewrite, except the `helper-brms.R` comment,
which the lane replaced with a better one attached to code that now
actually does what it says.
