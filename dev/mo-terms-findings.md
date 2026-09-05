# One simplex per mo() term (wt-mo-terms)

## Ground truth: how brms enumerates monotonic simplexes

Probed with brms 2.23.0, `make_standata()` and `get_prior()`, from
throwaway scripts under the scratchpad that are not kept. `Imo` is the
number of SIMPLEXES, not of terms; `Jmo[i]` is simplex i's dimension.
`brms:::brmsterms(bf(f))$dpars$mu$sp` reads the term labels directly and
is the cheaper check.

| formula | brms sp term labels | Imo | Jmo | simo coefs |
|---|---|---|---|---|
| `y ~ mo(inc) + z` | `mo(inc)` | 1 | 3 | `moinc1` |
| `y ~ mo(inc) * z` | `mo(inc)`, `mo(inc):z` | 2 | 3, 3 | `moinc1`, `moinc:z1` |
| `y ~ z * mo(inc)` | `mo(inc)`, `mo(inc):z` | 2 | 3, 3 | same |
| `y ~ mo(inc) + mo(inc):z` | `mo(inc)`, `mo(inc):z` | 2 | 3, 3 | same |
| `y ~ mo(w) + mo(inc):z` | `mo(w)`, `mo(inc):z` | 2 | 2, 3 | `mow1`, `moinc:z1` |
| `y ~ mo(inc) + mo(inc):z + mo(w)` | `mo(inc)`, `mo(w)`, `mo(inc):z` | 3 | 3, 2, 3 | `moinc1`, `mow1`, `moinc:z1` |
| `y ~ mo(inc) * mo(w)` | `mo(inc)`, `mo(w)`, `mo(inc):mo(w)` | 4 | 3, 2, 3, 2 | `moinc1`, `mow1`, `moinc:mow1`, `moinc:mow2` |

Two rules, both confirmed by the `Jmo`/`Xmo` fingerprints (variables with
different category counts disambiguate which term owns which simplex):

1. Terms are enumerated in `terms()` order: every main effect before any
   interaction, formula order within one interaction order. brms does not
   keep the written order.
2. Within one term, one simplex per mo() variable, in the order the
   variables appear in the term. Only reachable for `mo(x):mo(w)`, which
   frmtmb refuses outright.

frmtmb before this branch, from `frm(..., dry_run = "frame")`:

| formula | frmtmb mo column order |
|---|---|
| `y ~ mo(inc) * z` | `moinc:z`, `moinc` (interaction FIRST) |
| `y ~ mo(w) + mo(inc):z` | `mow`, `moinc:z` (agrees) |
| `y ~ mo(inc) + mo(inc):z + mo(w)` | `moinc`, `moinc:z`, `mow` (disagrees) |

`R/parse.R:896-899` appends the `mult` entry before the `mult = NULL`
main effect it implies, and never reorders by interaction order. So the
mo list has to be sorted before the frame numbers the simplexes.

`y ~ mo(inc) * z * u` is refused by frmtmb (parse.R: the special must be
on exactly one side of `:` / `*`), so no three-way case exists to order.
brms 2.23.0 itself errors on `y ~ mo(inc) * g` for a 3-level factor `g`.

## Prior vocabulary

brms: class `simo`, coef `<b coefficient name><index within term>`, so
`moinc1` for the main effect and `moinc:z1` for the interaction, prior
`dirichlet(1)`.

frmtmb `set_prior()` does not accept class `simo` at all
(`R/priors.R:264-265` match.arg list), and nothing in frmtmb emits a
`simo` row. There is no simplex prior to rename.

## What changed

`R/frame.R`

* `mo_terms_in_brms_order()` and `mo_interaction_order()`, new, sort the
  parsed mo list by interaction order before the frame walks it.
  `order()` is stable, so terms of equal order keep the written order,
  which is what `terms(keep.order = FALSE)` does. frmtmb allows only a
  single non-special multiplier, so the orders in play are 1 and 2; the
  count is written generically so a deeper multiplier would still sort
  right.
* The `mo_zetas` cache is gone. Each term takes
  `paste0("zeta", length(extras) + 1L)` as before. Nothing but zeta
  enters `extras` before the mo loop (`bcs` comes after), so that
  expression is the mo index only when `extras` was empty at that
  point. It is NOT: `R/frame.R:1420` replaces `extras` with the
  family's `extra_pars` earlier in the same response loop, so the four
  ordinal families, `cox()` and `mixture_mvn` all push the first
  monotonic simplex to `zeta2`. What lines up with brms is the
  POSITION in the mo list, not the number in the name. The numbering
  is unchanged from 564e185; only the documentation of it is new.
* Sorting the list also reorders the monotonic X columns, so
  `y ~ mo(x) * z` now reports `mox` before `mox:z`. That is what makes
  the tier's `bsp` rule right: it reads the non-X coefficients out of
  `fixef()` in frmtmb's order and hands them to brms's `bsp`, which is
  in brms's special-term order.

`R/parse.R`: the comment claiming the sharing was corrected. No code
change; the parser still emits the interaction first and the frame
no longer cares.

`R/objective.R`, `R/predict.R`, `R/conditional-effects.R`: NO change
needed. `lp_eta_fixed()`, `mo_col_values()`, `patch_mo_cols()`,
`pred_design()` and `ce_lp_vars()` all loop over `lp[["mo"]]` and read
`mi$zeta` per entry, so they were already per-term. The sharing lived
entirely in the frame, as an aliasing of the zeta NAME across entries.

## Prior naming

No rename was possible, because frmtmb has no simplex prior.
`set_prior()`'s class vocabulary (`R/priors.R:264`) is b, Intercept,
sd, cor, theta, ar, ma, cosy, cortime, rescor; `get_prior()` emits no
row for a `zeta<j>`; and `check_brms_prior_class()` refuses a brms
`simo` row with an explicit prior by name, with a message telling the
caller to write what they mean. A brms `get_prior()` object still
translates, because its `simo` rows are flat defaults and the
translator drops those before the class check.

The hint in that refusal named five of the ten classes `set_prior()`
takes and told the caller to write `class = "Intercept", dpar = "simo"`,
which means nothing for a simplex. `?frm` now documents the refusal, so
the message became reachable by following the docs; it gained a `simo`
branch saying there is no slot to carry the row into, and the generic
branch now lists all ten classes.

What the branch does deliver is the correspondence a prior would need:
brms's `simo_<j>` is the j-th entry of the frame's mo list, whose
`label` is brms's b coefficient name, so brms's simo coef is
`paste0(label, 1)` (frmtmb permits one mo() per term, so the index
inside the term is always 1). `test-mo-terms.R` asserts that against
`brms::get_prior()`. The `zeta` NUMBER is deliberately not part of the
correspondence, for the reason above.

Component names keep the existing scheme, unnamed inside, so
`par_template_names()` keeps spelling them `zeta1_1`, `zeta1_2` and
`confint()`, `variables()` and `start =` are unchanged. Renumbering so
monotonic simplexes always start at 1 would be the better end state,
but it moves those spellings for every ordinal and cox fit and is a
separate change.

The tier's translator was rewritten to match: `simo_<j>` resolves
through `brms_mo_terms_of(fit, "mu")[[j]]$zeta` rather than through
`paste0("zeta", j)`, which is what lets an ordinal fit translate at
all. `test-mo-terms.R` covers a `cumulative()` fit whose simplexes are
`zeta2` and `zeta3`.

## brms log-density tier

rstan 2.32.7 with StanHeaders 2.39.1 needs `R_MAKEVARS_USER` pointing at
a file whose `CXX17FLAGS` ends in `-std=gnu++17`, pointed at by
`R_MAKEVARS_USER`, and a `FRMTMB_STAN_CACHE` directory. Both are local
scaffolding, written under the scratchpad and not kept in the tree.

Run each row in ITS OWN PROCESS. `testthat::test_file()` re-sources
`helper-brms.R`, which resets the in-session `.brms_stan_models`
environment; a stanfit RDS written earlier in the SAME session then
fails to reload with a null DSO pointer, which the helper documents at
`helper-brms.R:85-95`. Five `test_file()` calls in one process made row
2 error on a program the translator test had just compiled.

| row | result | admitted constant |
|---|---|---|
| the translator round-trips through Stan's constraints | 11 pass | none |
| the simplex and group-level rules round-trip | 30 pass | none |
| row 2, `y ~ mo(inc) + z` | 2 pass | `lgamma(3)` |
| row 3, `y ~ mo(inc):z` | 2 pass | `lgamma(3)` |
| row 3, `y ~ mo(inc) * z` | 6 pass | `2 * lgamma(3)` |
| check C, `y ~ mo(inc):z + (1 | g)` | 3 pass | `lgamma(3)` |

Re-run unchanged after the translator was rewritten to map by position
(2, 2, 6, 3, 11, 30), from a cold cache, one row per process.

The last row is the flip. It used to assert the structural difference
(`Imo == 2` against a single `zeta1`); it now asserts `zeta1`, `zeta2`
and then runs `brms_lp_check()`, so checks A and B both pass with the
flat Dirichlet admitted once per simplex. The tier's exemption list is
empty, and the header comment says so.

## frm_compat()

No change. The resolved table has 106 rows naming `mo()`, and none
claims a shared simplex. The two whose note matches "shar" are the
`gr_cov` and `gr_prec` rows, whose "Terms sharing an |ID| key" is about
merging grouping blocks. The nearest rows stay true as written:
`mo() x kind:family` says the interaction multiplier must be a single
numeric column "because the simplex carries one coefficient", which is
about the multiplier and not about sharing a shape; `mo() x mi_pred()`
still refuses `mo():mo()` and `mo():mi()`. `mo() x prior` is still
"untested", which is honest: there is no simplex prior to test.

## as-cran

`R CMD build` then `R CMD check --as-cran --no-manual` with
`_R_CHECK_CRAN_INCOMING_=false` and the quarto tools on PATH:
**Status: OK**, no ERRORs, no WARNINGs, no NOTEs. The suite runs inside
it (`checking tests ... [34m] OK`) and the vignettes re-build clean
(356s), so the behavior change did not break any vignette that fits a
monotonic model.

## Deliberate omissions

* `mi()` interaction terms are NOT re-sorted, and the reason is NOT
  that they are unconstrained by brms. brms puts mo() and mi() terms in
  ONE `bsp`, interleaved in a single terms() order across both, so for
  `y ~ mo(inc) * z + mi(x) * z` it wants
  `mo(inc) | mi(x) | mo(inc):z | mi(x):z` where frmtmb reports
  `moinc | moinc:z | mix:z | mix`: three of four positions disagree.
  The omission is safe only because nothing asserts it - there is no
  `mi(` in `test-brms-likelihood.R` or in `test-brms-agreement.R`, and
  564e185 had the same disagreement. Sorting the mi() list the same way
  is the obvious follow-up for whoever writes the first mo()-plus-mi()
  tier row, and it must come before that row, not after.
* No `simo` prior class was added. That is a feature (a Dirichlet
  penalty in the objective, a new `set_prior()` class, `get_prior()`
  rows), not the naming alignment this branch was asked for.
* A `mo()` term in a non-mu dpar still numbers into the same global
  `zeta<j>` sequence, while brms spells it `simo_<dpar>_<j>`. The fit
  is right; it is the tier's translator that had no answer, and reading
  the trailing integer alone silently handed such a name mu's simplex,
  at mu's length. It now refuses a `simo` name carrying a dpar suffix
  BY NAME, with a test. Matching brms's per-dpar numbering is a
  separate decision and was not taken: no tier row fits mo() outside
  mu.
* `mo(x):mo(w)` gets two simplexes in brms. frmtmb refuses mo():mo()
  outright (`R/parse.R`), so the within-term rule has nothing to apply
  to and none was written.
* `dev/brms-likelihood-tests.md`, `dev/brms-vignette-audit.md` and
  `dev/reviews/2026-09-05-brms-likelihood-identity.md` still describe
  the divergence as live. They are dated records of the investigation
  and were not in the change list; NEWS carries the correction.

## Scope

Three files in the brms tier changed, not two:

* `tests/testthat/test-brms-likelihood.R`, two hunks - the header
  paragraph that describes the exemption list, and the row 3b block.
* `tests/testthat/helper-brms.R`, the `simo_` translator rule, plus the
  `brms_mo_terms_of()` accessor it needs. The sibling lane owns
  everything else in both files and none of it was touched.

Outside the tier: `R/frame.R` (the mo section and two new sort
helpers), `R/parse.R` (a comment), `R/fit.R` (roxygen), `R/priors.R`
(the brms prior-class hint), `man/frm.Rd` (generated), `NEWS.md`,
`tests/testthat/test-v15.R` (a comment), `tests/testthat/test-v18.R`,
and the new `tests/testthat/test-mo-terms.R`.

## Core suite

One file per process under `pkgload::load_all` with `NOT_CRAN=true`,
run twice: once against a pristine copy of the branch point (564e185)
extracted with `git archive`, once against the working tree.

|  | files | pass | fail | error | skip |
|---|---|---|---|---|---|
| 564e185 | 104 | 5423 | 0 | 0 | 21 |
| this branch | 105 | 5473 | 0 | 0 | 22 |

Three files moved, and only three:

* `test-mo-terms.R`, new, 53 passes. The simplex-count contract, the
  enumeration order against `Imo`/`Jmo`/`Xmo_j`, brms's simo coef
  names, the likelihood ordering against the shared-simplex submodel,
  and a `mo()`-with-`mi()` two-way interaction model.
* `test-v18.R`, 26 -> 27. Its reference NLL used one simplex for both
  monotonic columns, which WAS the old model, so it had to be rewritten
  with two; the direct-ML identity still holds to 7e-11. The
  `expect_length(grep("^zeta", ...), 1L)` became `2L` (same count), and
  one assertion was added: the old shared-simplex optimum, fitted
  directly, is a strict lower bound on this data.
* `test-brms-likelihood.R`, 4 -> 0 passes and 18 -> 19 skips WITHOUT
  `FRMTMB_BRMS_FIT_TESTS`. Row 3b used to assert a structural
  difference, which needed only brms and so ran here; it now runs
  `brms_lp_check()`, which needs Stan, so its guard moved from
  `skip_unless_brms()` to `skip_unless_brms_fit()` and it skips in a
  plain run. With the env set it is 6 passes.

`test-v15.R` changed by a comment only and its count is unmoved (41).

