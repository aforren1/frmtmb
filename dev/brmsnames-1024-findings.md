# Coefficient names that differ from brms: lane `brmsnames`, 2026-10-24

Worktree `C:\Users\adf44\source\r\frmtmb-wt-brmsnames`, branch
`wt-brmsnames`, base 18b536fc (frmtmb 0.74.0). The BEFORE arm reads
`C:/Users/adf44/source/r/rellib-r12`; the AFTER arm ("lane") reads
`C:/Users/adf44/source/r/wt-brmsnames-lib` first, which holds core and
every extension built from this worktree. brms 2.23.0, rstan 2.32.7.

This record is `dev/brmsnames-1024-findings.md`, not
`dev/brmsnames-findings.md` as the brief names it: that path is the
tracked, 857-line record of the earlier lane of the same name (item
2.6c, frmtmb 0.58.0), and writing over it would lose that record. Every
script and log of this lane is prefixed `dev/brmsnames-1024-` for the
same reason; the earlier lane's `dev/brmsnames-*.R` are untouched.

The blocks between `BEGIN` and `END` markers are pasted from the logs
by `dev/brmsnames-1024-fill.R`; the logs themselves are local
(`dev/brmsnames-1024-log/`, gitignored).

## How to reproduce every number here

| what | script | output |
| --- | --- | --- |
| the cases and their data (seed 1024) | `dev/brmsnames-1024-cases.R` | |
| brms's names, from a real brmsfit | `dev/brmsnames-1024-brms.R <case> ...` | `dev/brmsnames-1024-out/brms-<n>.rds` |
| frmtmb's names, per arm | `dev/brmsnames-1024-frm.R base\|lane <case> ...` | `dev/brmsnames-1024-out/frm-<arm>-<n>.rds` |
| the name table | `dev/brmsnames-1024-table.R` | `table.md` |
| the raw names of every case | `dev/brmsnames-1024-dump.R` | `dump.txt` |
| names of every fit of a suite, per arm | `dev/brmsnames-1024-suitenames.sh` (runs `dev/brmsnames-1024-suitenames.R`) | `dev/brmsnames-1024-suite/<arm>/<pkg>/` |
| the before/after diff of those names | `dev/brmsnames-1024-suitediff.R [pkg]` | `suitediff.md`, `suitediff-sample.md` |
| the ported brms suite, record mode, per arm | `dev/brmsnames-1024-port.sh base\|lane` | `dev/brmsnames-1024-port-out/<arm>/` |
| which ported rows flip | `dev/brmsnames-1024-portdiff.R` | `portdiff.md` |
| brms's rename of odd columns | `dev/brmsnames-1024-probe3.R` | `probe3.txt` |
| the parser's protected names | `dev/brmsnames-1024-probe0.R` | |
| mo() coding, brms and frmtmb | `dev/brmsnames-1024-probe1.R`, `-probe2.R` | `probe1.txt`, `probe2.txt` |
| drop1()/step() and re_formula, per arm | `dev/brmsnames-1024-bugfix.R` | `bugfix-<arm>.txt` |
| prior table row order against brms, per arm | `dev/brmsnames-1024-order.R` | `order-<arm>.txt` |
| the alias machinery, before and after | `git grep` (command in "Fix round 1") | `alias-grep.txt` |
| the record's generated blocks | `dev/brmsnames-1024-fill.R` | this file |

brms's names are never recalled. `dev/brmsnames-1024-brms.R` runs
`standata()`, `stancode()` and `default_prior()`, then builds a real
brmsfit the way `help("brm")` documents it (example `fit8`): an empty
brmsfit (`brm(empty = TRUE)`), a stanfit of brms's own Stan code and
standata (one chain, 40 iterations, seed 1), and `brms::rename_pars()`.
`variables()`, `fixef()`, `vcov()`, `summary()`, `coef()`, `ranef()`,
`VarCorr()`, `conditional_effects()` and `hypothesis()` are then brms's
own methods on that fit. The compiled programs are cached in
`dev/stan-cache` under the key the test suite's `brms_stan_model()` uses.

## What was wrong

1. **The alias.** `parse_linpred()` (`R/parse.R`) renames every call
   whose function has the name of a covariance structure, so that
   `reformulas::splitForm()` does not drop the term
   (`dev/rtmb-pitfalls.md` item 9). The names are `cs_specials`:
   `cs_names()` less `smooth`, `us_t` and `diag_t`, plus `rr` and
   `propto` (`dev/brmsnames-1024-probe0.R`): `us`, `diag`, `homdiag`,
   `ar1`, `cs`, `ou`, `toep`, `hetar1`, `homcs`, `homtoep`, `exp`,
   `gau`, `mat`, `gr_cov`, `gr_prec`, `gp`, `hsgp`, `equalto`, plus any
   registered structure (`car`, `spde` with frmtmb.spatial). Only a
   name that R finds as a function is renamed, so in practice `exp()`,
   `diag()`, and `cs()` with brms attached. The alias `.frm_exp(x)`
   stayed in the fixed formula and in the bar terms after splitForm,
   so `model.matrix()` named the column `.frm_exp(x)` and every name
   built from the column carried it: `.frm_expx` in `fixef()`, `vcov()`,
   `summary()`, `coef()`, `ranef()`, `VarCorr()`; `b_.frm_expx`,
   `sd_g__.frm_expx`, `cor_g__Intercept__.frm_expx` in `variables()`
   and the draws; `.frm_exp(x)` in `confint()` and `default_prior()`.
   Since 0.15.0 (65fed080, 2026-08-31), when the alias was added.
   `cs()` inside a bar already had its own repair of the alias
   (`check_cs_in_bar()`, `R/frame.R`).
2. **The prior table's coef column.** `prior_coef_brms()`
   (`R/priors.R`, since 0.65.0, 19714a68) turned only `(Intercept)`
   and a smooth's `s(x).fx1` into brms's spelling, so every other
   transformed column was listed under its `model.matrix()` name
   (`log(xp)`, `I(x^2)`, `poly(x, 2)1`, `factor(k)2`), where brms lists
   its `rename()` (`logxp`, `IxE2`, `polyx21`, `factork2`). The
   resolver took the column name, the column name without parentheses,
   or that spelling. So `coef = "IxE2"`, `"polyx21"`,
   `"polyx2rawEQTRUE1"` and `"expxD2"` were refused with "Prior target
   not found", and `coef = "logxp"` and `"factork2"` were taken (they
   are the column name with its parentheses dropped) while the table
   listed `log(xp)` and `factor(k)2` (`dump.txt`, `prior_ok`).
3. **The fixef() order of a nonlinear parameter.** `brms_fixef_rows()`
   (`R/brms-shapes.R`) puts every intercept first, as brms does for a
   centered predictor, but brms does not center a nonlinear parameter's
   formula, so its intercept is one element of its `b` vector. Found
   by this lane's table (cases 31 and 39), not in the backlog.

## What changed

- `R/parse.R`: new `unsub_specials()`; `parse_linpred()` undoes the
  alias in `sf$fixedFormula`, `sf$reTrmFormulas` and `sf$reTrmAddArgs`
  right after `splitForm()`. The alias still protects the term from
  splitForm and never leaves the parser. With the term as written, the
  design column, the model frame column and every name built from them
  are R's and brms's.
- `R/priors.R`: `prior_coef_brms()` is `brms_rename()` of the design
  column (the smooth columns keep their `sx_1` rule), so the table lists
  brms's names, sorts them as brms sorts them, and the resolver takes
  them. The column name still resolves (`prior_coef_match()` is
  unchanged). The `@param coef` of `set_prior()` says so.
- `R/brms-shapes.R`: `brms_fixef_rows()` keeps the intercept of a
  nonlinear parameter in its place.
- Tests: `tests/testthat/test-brms-names-specials.R` (core) and
  `extensions/frmtmb.sample/tests/testthat/test-draws-names-specials.R`;
  their counts on each arm are in "The tests, seen to fail on 0.74.0".
- NEWS: BREAKING bullets and two bug fixes in core, two bullets in
  frmtmb.sample.

Fix round 1 (review `dev/reviews/2026-10-24-brmsnames.md`) added, per
finding (section "Fix round 1" below has the evidence):

- F1: no code; NEWS bug-fix bullets and two tests for `drop1()`/`step()`
  and `re_formula`, which the alias fix repairs.
- F2: `prior_table()` sorts its rows as brms's `default_prior()` does,
  and `prior_rows_sd_order()` (the filled table of `validate_prior()`
  and `prior_summary()`) as brms's `validate_prior()` does; class
  `"theta"` stays last. Chosen over stating the limit: two lines each,
  and the review found nothing that reads a prior table row by
  position, which this lane's runs of every suite confirm.
- F3: `prior_coef_match_brms()` in `R/priors.R` lets class `"sd"`
  (`sd_spec_reach()`, `sd_coef_shared()`) and a `cs()` term take the
  design column's spelling; `prior_spec_slot()` keys a coef by
  `brms_rename()`, so the two spellings are one slot.
- F4: tests of the old spelling through `set_prior()`, `prior_string()`,
  `validate_prior()`, `frm()` and `frm_sample()`, each checking the
  coefficient a tight prior moved.
- F6: `brms_lp_prefix()` (`R/brms-names.R`) writes a nonlinear
  parameter's prefix as brms does, resp then nlpar (`yn_a`), and
  `brms_fixef_rows()` orders a multivariate model's rows response by
  response, as brms declares them (found on the same case).
- F5: `fixef.frmtmb_draws()` (frmtmb.sample `R/methods-draws.R`) orders
  its rows by core's exported `brms_fixef_rows()`, the fit's order.
- F8: `sub_specials()` binds nothing; the `prot` environment, the
  `shared` argument of `parse_linpred()` and the `.frm_cs` branch of
  `check_cs_in_bar()` are gone; comments and `dev/rtmb-pitfalls.md`
  item 9 updated.
- F7: filed below ("Defects found and not fixed", item 1), not fixed.
- F9: left to consolidation.

frmtmb.sample's `R/` changed in one place (F5). No extension reads a
`.frm_` name (`grep -rn '\.frm_' extensions/*/R/` finds none);
frmtmb.sample takes its names from core's `brms_coef_table()` and
`brms_par_labels()`, so its draws follow the fix (see the sample test
and the sample suite diff below).

### What brms is right about, and the one place frmtmb keeps its own names

brms is not clearly wrong on any row this lane changed. `confint()` is
the one accessor that keeps a different vocabulary, on purpose: brms
has no `confint()` method (brms 2.23.0's NAMESPACE has none), and
frmtmb's rows are the internal parameters on their internal scale
(`sigma_(Intercept)` is log sigma, `theta_1` a log sd), with the row
names of `vcov(full = TRUE)`, as `?confint.frmtmb_fit` documents.
Renaming those rows `sigma` or `sd_g__x` would name a value that is
not the one in the row. The fix removes the alias there too: the row
of `exp(x)` is `exp(x)`, the design column R writes, and brms's
`b_expx` and the bare `expx` address it (`parm =`), which the new test
pins.

## The name table

Each term of `dev/brmsnames-1024-cases.R` fitted by brms 2.23.0 and
by both frmtmb arms. Table 1 is `fixef()`; Table 2 compares every
interface the brief lists. A "DIFF" that remains on the lane is
explained after the table.

<!-- BEGIN table.md -->


Table 1. fixef() row names.

| # | case | brms 2.23.0 | frmtmb 0.74.0 | lane |
|---|---|---|---|---|
| 1 | `exp` | Intercept expx | Intercept .frm_expx | Intercept expx |
| 2 | `log` | Intercept logxp | Intercept logxp | Intercept logxp |
| 3 | `sqrt` | Intercept sqrtxp | Intercept sqrtxp | Intercept sqrtxp |
| 4 | `Isq` | Intercept IxE2 | Intercept IxE2 | Intercept IxE2 |
| 5 | `poly` | Intercept polyx21 polyx22 | Intercept polyx21 polyx22 | Intercept polyx21 polyx22 |
| 6 | `polyraw` | Intercept polyx2rawEQTRUE1 polyx2rawEQTRUE2 | Intercept polyx2rawEQTRUE1 polyx2rawEQTRUE2 | Intercept polyx2rawEQTRUE1 polyx2rawEQTRUE2 |
| 7 | `scale` | Intercept scalex | Intercept scalex | Intercept scalex |
| 8 | `factor` | Intercept factork2 factork3 factork4 | Intercept factork2 factork3 factork4 | Intercept factork2 factork3 factork4 |
| 9 | `exp:z` | Intercept expx:z | Intercept .frm_expx:z | Intercept expx:z |
| 10 | `exp*f` | Intercept expx fb fc expx:fb expx:fc | Intercept .frm_expx fb fc .frm_expx:fb .frm_expx:fc | Intercept expx fb fc expx:fb expx:fc |
| 11 | `log:factor` | Intercept logxp:factork1 logxp:factork2 logxp:factork3 logxp:factork4 | Intercept logxp:factork1 logxp:factork2 logxp:factork3 logxp:factork4 | Intercept logxp:factork1 logxp:factork2 logxp:factork3 logxp:factork4 |
| 12 | `poly*z` | Intercept polyx21 polyx22 z polyx21:z polyx22:z | Intercept polyx21 polyx22 z polyx21:z polyx22:z | Intercept polyx21 polyx22 z polyx21:z polyx22:z |
| 13 | `I(exp)` | Intercept Iexpx | Intercept I.frm_expx | Intercept Iexpx |
| 14 | `exp+log+sqrt` | Intercept expx logxp sqrtxp | Intercept .frm_expx logxp sqrtxp | Intercept expx logxp sqrtxp |
| 15 | `exp(x/2)` | Intercept expxD2 | Intercept .frm_expxD2 | Intercept expxD2 |
| 16 | `s` | Intercept sx_1 | Intercept sx_1 | Intercept sx_1 |
| 17 | `s+exp` | Intercept expz sxp_1 | Intercept .frm_expz sxp_1 | Intercept expz sxp_1 |
| 18 | `s(exp)` | Intercept sexpx_1 | (error) | (error) |
| 19 | `t2` | Intercept t2xz_1 t2xz_2 t2xz_3 | Intercept t2xz_1 t2xz_2 t2xz_3 | Intercept t2xz_1 t2xz_2 t2xz_3 |
| 20 | `mo` | Intercept mok | Intercept mok | Intercept mok |
| 21 | `mo:z` | Intercept mok:z | Intercept mok:z | Intercept mok:z |
| 22 | `mo+exp` | Intercept expx mok | Intercept .frm_expx mok | Intercept expx mok |
| 23 | `me` | Intercept mexosdx | Intercept mexosdx | Intercept mexosdx |
| 24 | `me+exp` | Intercept expz mexosdx | Intercept .frm_expz mexosdx | Intercept expz mexosdx |
| 25 | `mi` | y_Intercept xm_Intercept y_expx xm_z y_mixm | y_Intercept xm_Intercept y_.frm_expx xm_z y_mixm | y_Intercept xm_Intercept y_expx xm_z y_mixm |
| 26 | `cs` | Intercept[1] Intercept[2] Intercept[3] z x[1] x[2] x[3] | Intercept[1] Intercept[2] Intercept[3] z x[1] x[2] x[3] | Intercept[1] Intercept[2] Intercept[3] z x[1] x[2] x[3] |
| 27 | `cs(exp)` | Intercept[1] Intercept[2] Intercept[3] z expx[1] expx[2] expx[3] | Intercept[1] Intercept[2] Intercept[3] z expx[1] expx[2] expx[3] | Intercept[1] Intercept[2] Intercept[3] z expx[1] expx[2] expx[3] |
| 28 | `offset` | Intercept x | Intercept x | Intercept x |
| 29 | `offset+exp` | Intercept expx | Intercept .frm_expx | Intercept expx |
| 30 | `nl` | a_Intercept b_Intercept | a_Intercept b_Intercept | a_Intercept b_Intercept |
| 31 | `nl a~exp` | a_Intercept a_expz b_Intercept | a_Intercept b_Intercept a_.frm_expz | a_Intercept a_expz b_Intercept |
| 32 | `sigma~exp` | Intercept sigma_Intercept x sigma_expz | Intercept sigma_Intercept x sigma_.frm_expz | Intercept sigma_Intercept x sigma_expz |
| 33 | `re exp` | Intercept expx | Intercept .frm_expx | Intercept expx |
| 34 | `re log` | Intercept z | Intercept z | Intercept z |
| 35 | `mv exp` | y_Intercept y2_Intercept y_expx y2_logxp | y_Intercept y2_Intercept y_.frm_expx y2_logxp | y_Intercept y2_Intercept y_expx y2_logxp |
| 36 | `re exp ||` | Intercept expx | Intercept .frm_expx | Intercept expx |
| 37 | `re mo` | Intercept mok | Intercept mok | Intercept mok |
| 38 | `offset(exp)` | Intercept x | Intercept x | Intercept x |
| 39 | `nl a~exp sigma~z` | sigma_Intercept a_Intercept a_expz b_Intercept b_z sigma_z | a_Intercept b_Intercept sigma_Intercept a_.frm_expz b_z sigma_z | sigma_Intercept a_Intercept a_expz b_Intercept b_z sigma_z |
| 40 | `mv nl` | sigma_yn_Intercept y_Intercept yn_a_Intercept yn_a_expz yn_b_Intercept sigma_yn_z y_x y_expz | a_yn_Intercept b_yn_Intercept y_Intercept sigma_yn_Intercept a_yn_.frm_expz y_x y_.frm_expz sigma_yn_z | sigma_yn_Intercept y_Intercept yn_a_Intercept yn_a_expz yn_b_Intercept sigma_yn_z y_x y_expz |

Table 2. Each interface against brms 2.23.0, base / lane.

| # | case | fixef | vcov | summary | variables | default_prior | set_prior | hypothesis | draws | coef / ranef / VarCorr | conditional_effects |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 2 | `log` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 3 | `sqrt` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 4 | `Isq` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 5 | `poly` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 6 | `polyraw` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 7 | `scale` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 8 | `factor` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | n/a / n/a |
| 9 | `exp:z` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | DIFF / DIFF |
| 10 | `exp*f` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 11 | `log:factor` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | n/a / n/a |
| 12 | `poly*z` | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 13 | `I(exp)` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 14 | `exp+log+sqrt` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 15 | `exp(x/2)` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 16 | `s` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / DIFF | ok / ok |
| 17 | `s+exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / DIFF | ok / ok |
| 18 | `s(exp)` | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR | ERR / ERR |
| 19 | `t2` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / DIFF | DIFF / DIFF |
| 20 | `mo` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / DIFF | ok / ok | ok / ok |
| 21 | `mo:z` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / DIFF | ok / ok | DIFF / DIFF |
| 22 | `mo+exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / DIFF | ok / ok | ok / ok |
| 23 | `me` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 24 | `me+exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 25 | `mi` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | ok / ok | DIFF / DIFF |
| 26 | `cs` | ok / ok | ok / ok | ok / ok | order / order | ok / ok | ok / ok | ok / ok | order / order | n/a / n/a | DIFF / DIFF |
| 27 | `cs(exp)` | ok / ok | ok / ok | ok / ok | order / order | ok / ok | ok / ok | ok / ok | order / order | n/a / n/a | DIFF / DIFF |
| 28 | `offset` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | n/a / n/a | ok / ok |
| 29 | `offset+exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | n/a / n/a | ok / ok |
| 30 | `nl` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 31 | `nl a~exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok | ok / ok |
| 32 | `sigma~exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | n/a / n/a | DIFF / DIFF |
| 33 | `re exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok |
| 34 | `re log` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | order / order | ok / ok |
| 35 | `mv exp` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | ok / ok | DIFF / DIFF |
| 36 | `re exp ||` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / ok | ok / ok |
| 37 | `re mo` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | DIFF / DIFF | ok / ok | ok / ok |
| 38 | `offset(exp)` | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok | ok / ok |
| 39 | `nl a~exp sigma~z` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | n/a / n/a | ok / ok |
| 40 | `mv nl` | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | DIFF / ok | DIFF / ok | DIFF / ok | DIFF / order | DIFF / DIFF | DIFF / DIFF |

Tally over 40 cases (arm interface verdict: count).

| arm | interface | verdict | cases |
|---|---|---|---|
| base | ce | DIFF | 9 |
| base | ce | n/a | 2 |
| base | ce | ok | 28 |
| lane | ce | DIFF | 9 |
| lane | ce | n/a | 2 |
| lane | ce | ok | 28 |
| base | draws | DIFF | 21 |
| base | draws | ok | 16 |
| base | draws | order | 2 |
| lane | draws | DIFF | 4 |
| lane | draws | ok | 28 |
| lane | draws | order | 7 |
| base | fixef | DIFF | 18 |
| base | fixef | ok | 21 |
| lane | fixef | ok | 39 |
| base | hyp | DIFF | 18 |
| base | hyp | ok | 21 |
| lane | hyp | ok | 39 |
| base | prior | DIFF | 27 |
| base | prior | ok | 12 |
| lane | prior | ok | 39 |
| base | re | DIFF | 6 |
| base | re | n/a | 6 |
| base | re | ok | 26 |
| base | re | order | 1 |
| lane | re | DIFF | 4 |
| lane | re | n/a | 6 |
| lane | re | ok | 28 |
| lane | re | order | 1 |
| base | set_prior | DIFF | 27 |
| base | set_prior | ok | 12 |
| lane | set_prior | ok | 39 |
| base | summary | DIFF | 18 |
| base | summary | ok | 21 |
| lane | summary | ok | 39 |
| base | variables | DIFF | 18 |
| base | variables | ok | 19 |
| base | variables | order | 2 |
| lane | variables | ok | 32 |
| lane | variables | order | 7 |
| base | vcov | DIFF | 18 |
| base | vcov | ok | 21 |
| lane | vcov | ok | 39 |

Cases with `.frm_` in any recorded output: base 18 of 40 , lane 0 of 40 

<!-- END table.md -->

### Every lane cell that is not "ok", and why it stays

- **draws DIFF on `mo`, `mo:z`, `mo+exp`, `re mo` (cases 20 to 22,
  37): a model defect, not a name.** brms subtracts the minimum of an
  integer `mo()` predictor (`get_mo_values()`, brms
  `R/formula-sp.R` lines 603 to 627: `x <- x - min_value`), so `k` in
  1 to 4 has `Jmo = 3` and a simplex `simo_mok1[1..3]`. frmtmb codes
  the integers as they are (`R/frame.R`, `codes <- as.integer(v)`,
  `D_mo <- max(codes)`), so `D = 4` and the draws hold
  `simo_mok1[1..4]`. The extra cell is a category no row is in: the
  fit warns that `(Intercept)`, `mok` and `zeta1_3` have no standard
  error, and `mo(k - 1)` gives the same log likelihood with one
  parameter fewer (`probe1.txt`, `probe2.txt` below). Out of this
  lane's scope (the model, its newdata coding and its tests would
  change); filed below. It is in both arms.
- **variables and draws "order" on `mi`, `cs`, `cs(exp)`, `sigma~exp`,
  `mv exp`, `nl a~exp sigma~z`, `mv nl`.** Same names, another order:
  brms puts
  every centered intercept first, then the `b` vectors predictor by
  predictor. Already recorded as found and not fixed by the earlier
  lane (`dev/brmsnames-findings.md`, "Found and NOT fixed", item 1):
  frmtmb.sample keeps the template order, which every positional read
  of the draws relies on. `fixef()`, `vcov()` and `summary()` are in
  brms's order on every case now.
- **coef / ranef / VarCorr on `s`, `s+exp`, `t2` (DIFF) and `re log`
  (order).** On a smooth-only model brms's `coef()` and `ranef()` stop
  ("No group-level effects detected"), and frmtmb's `ranef()` returns
  the smooth's penalized coefficients, `s(x).1` to `s(x).8`, because a
  smooth is a random-effect block there. On `re log`, brms's `coef()`
  lists the group-level terms first (`Intercept logxp z`,
  `brms:::coef.brmsfit` appends the population-only terms) and
  frmtmb's lme4-shaped `coef()` (a data frame per group) lists the
  population-level terms first (`Intercept z logxp`). On `mv nl`
  (case 40) brms's `coef()` lists every population-level coefficient
  of the model under the group, and frmtmb's those of the predictor
  that has the group-level term (`yn_a_Intercept yn_a_expz`), under
  the same names. None is a name defect; each is the shape frmtmb chose
  for these methods.
- **conditional_effects DIFF on 8 cases.** The NAMES of the effects are
  not built from coefficient names, and the alias never reached them;
  the difference is which effects are chosen by default: brms gives
  `x:z` alone for `exp(x):z` and frmtmb `x z x:z`; brms adds the
  interaction for `t2(x, z)` and `mo(k):z`; brms gives `z x` for
  `cs(x)` and frmtmb `z:cats__`; brms adds `z` for `sigma ~ exp(z)`;
  and brms prefixes the response on a multivariate model
  (`y.y_x`). Both arms agree. Filed below.
- **`s(exp)` (case 18) ERR on both arms.** `s(exp(x))` is refused, "The
  model uses `exp(x)`, which is not a column of `data`": a smooth of an
  expression is not supported. brms names it `sexpx_1` and
  `sds_sexpx_1`. Filed below.

The table's `default_prior` and `set_prior` columns leave out the rows
of the classes frmtmb refuses by name, with their reason, in
`brms_prior_class_refusal()`: `sds` (a smooth is an `sd` block with
`group = "s(x)"` here), `meanme` and `sdme` (no prior slot for me()'s
latent mean and sd). These are prior-vocabulary decisions of earlier
rounds, not names, and are unchanged.

brms's rename agrees with `brms_rename()` on the odd columns the suite
diff met (`probe3.txt`):

```text
<!-- BEGIN probe3.txt -->

y ~ splines::ns(x2, 3)     brms: splines::nsx231 splines::nsx232 splines::nsx233
                           lane: splines::nsx231 splines::nsx232 splines::nsx233  identical=TRUE
y ~ c3                     brms: c30.2050.402 c30.4020.599 c30.5990.796 c30.7960.994
                           lane: c30.2050.402 c30.4020.599 c30.5990.796 c30.7960.994  identical=TRUE
y ~ f                      brms: fcMd fe
                           lane: fcMd fe  identical=TRUE
y ~ I(x * kx)              brms: IxMUkx
                           lane: IxMUkx  identical=TRUE
y ~ x:log(x2)              brms: x:logx2
                           lane: x:logx2  identical=TRUE
y ~ log(abs(x) + 1)        brms: logabsxP1
                           lane: logabsxP1  identical=TRUE
y ~ I(x^2) + x2            brms: IxE2 x2
                           lane: IxE2 x2  identical=TRUE

<!-- END probe3.txt -->
```

## The tests, seen to fail on 0.74.0

`tests/testthat/test-brms-names-specials.R` computes brms's names at
test time from `brms::standata()` and `brms::default_prior()` (gate
`skip_unless_brms()`, no compile), for 15 formulas: `exp(x)`,
`exp(x):z`, `exp(x) * f`, `I(exp(x))`, `exp(x / 2)`,
`exp(x) + log(xp) + sqrt(xp)`, `I(x^2) + scale(z)`,
`poly(z, 2, raw = TRUE) + factor(k)`, `sigma ~ exp(z)`, a nonlinear
`a ~ exp(z)`, `(1 + exp(x) | g)`, `(1 + exp(x) || g)`, a multivariate
`exp(x)` and `log(xp)`, a Poisson `exp(x) + offset(log(e))`, and a
multivariate nonlinear model (fix round 1). It
checks `fixef()`, `vcov()`, `summary()`, `variables()`, `hypothesis()`
on every brms name, the prior table against brms's, `set_prior()` with
every brms coef through `validate_prior()`, and `VarCorr()`, `coef()`,
`ranef()` and `hypothesis(class = "sd")` of the group-level `exp(x)`.
Two blocks need no brms: no `.frm_` in any name a user sees (fixef,
variables, confint, coef, the prior table, `fixef(flatten = TRUE)`),
and the confint() rows and their brms aliases. The nonlinear order is
pinned against brms's `fixef()` of case 39, stored in the test with
its script. The frmtmb.sample test pins the draws of
`y ~ exp(x) + (1 + exp(x) | g)` against brms's `variables()` of case
33, stored the same way.

The `log`, `sqrt`, `I()`, `poly()`, `scale()` and `factor()` terms are
the absent case for the alias: their fixef() names were already
brms's on 0.74.0 (Table 1), and only their prior rows change.

```text
<!-- BEGIN seen-fail.txt -->

== newtest-base
lib: frmtmb 0.74.0 C:/Users/adf44/source/r/rellib-r12; frmtmb 0.74.0 C:/Users/adf44/source/r/rellib-r12
DETAIL FAIL [population-level names are brms's in every accessor] Expected `fx` to have the same values as `want`. | Actual: "Intercept", ".frm_expx" | Expected: "Intercept", "expx" | Needs: ".frm_expx" | Absent: "expx"
DETAIL FAIL [population-level names are brms's in every accessor] Expected bf(y ~ exp(x)) to be TRUE. | Differences: | `actual`:   FALSE | `expected`: TRUE  | 
DETAIL ERROR [population-level names are brms's in every accessor] <frmtmb_error/error/condition> | Error: Some parameters cannot be found in the model:  | 'b_expx' | variables() lists every name; brms's default class = "b" puts b_ before e
DETAIL FAIL [the prior table lists brms's coefficients and takes them] Expected `bns_key(fp)` to have the same values as `bns_key(bp)`. | Actual: "b|.frm_exp(x)||||" | Expected: "b|expx||||" | Needs: "b|.frm_exp(x)||||" | Absent: "b|expx|||
DETAIL ERROR [the prior table lists brms's coefficients and takes them] <frmtmb_error/error/condition> | Error: Prior target not found (class=b, coef=expx)
DETAIL FAIL [group-level names of exp(x) are brms's] Expected bf(y ~ exp(x) + (1 + exp(x) | g)) to be TRUE. | Differences: | `actual`:   FALSE | `expected`: TRUE  | 
DETAIL FAIL [group-level names of exp(x) are brms's] Expected `rownames(VarCorr(fit)$g$sd)` to have the same values as `sdc`. | Actual: "Intercept", ".frm_expx" | Expected: "Intercept", "expx" | Needs: ".frm_expx" | Absent: "expx"
DETAIL FAIL [group-level names of exp(x) are brms's] Expected `colnames(coef(fit)$g)` to have the same values as `c("Intercept", "expx")`. | Actual: "Intercept", ".frm_expx" | Expected: "Intercept", "expx" | Needs: ".frm_expx" | Absent: "ex
DETAIL FAIL [group-level names of exp(x) are brms's] Expected `if (length(dim(rf)) == 3L) dimnames(rf)[[3]] else colnames(rf)` to have the same values as `sdc`. | Actual: "Intercept", ".frm_expx" | Expected: "Intercept", "expx" | Needs: ".f
DETAIL ERROR [group-level names of exp(x) are brms's] <frmtmb_error/error/condition> | Error: Some parameters cannot be found in the model:  | 'sd_g__expx' | variables() lists every name; brms's default class = "b" puts b_ before each name,
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x)) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x):z) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x) * f) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ I(exp(x))) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x / 2)) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x) + log(xp) + sqrt(xp)) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ x, sigma ~ exp(z)) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ a + b * x, a ~ exp(z), b ~ 1, nl = TRUE) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x) + (1 + exp(x) | g)) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(y ~ exp(x) + (1 + exp(x) || g)) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected mvbf(bf(y ~ exp(x)), bf(y2 ~ log(xp)), rescor = FALSE) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected bf(yc ~ exp(x) + offset(log(e))) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [no placeholder reaches a name a user sees] Expected mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + exp(z), b ~ 1, sigma ~ z, nl = TRUE), bf(y ~ x + exp(z)), rescor = FALSE) to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FA
DETAIL FAIL [confint() names the design column as written] Expected `rownames(confint(fit))[1:6]` to be identical to `c("(Intercept)", "exp(x)", "fb", "fc", "exp(x):fb", "exp(x):fc")`. | Differences: |     actual           | expected       
DETAIL ERROR [confint() names the design column as written] <frmtmb_error/error/condition> | Error: Unknown parameter(s) in confint(parm =): b_expx. Available: (Intercept), .frm_exp(x), fb, fc, .frm_exp(x):fb, .frm_exp(x):fc, sigma_(Interce
DETAIL FAIL [a nonlinear parameter's intercept keeps its place, as in brms] Expected `rownames(fixef(fit))` to be identical to `want`. | Differences: |     actual            | expected              | [1] "a_Intercept"     - "sigma_Intercept
DETAIL FAIL [a nonlinear parameter's intercept keeps its place, as in brms] Expected `rownames(vcov(fit))` to be identical to `want`. | Differences: |     actual            | expected              | [1] "a_Intercept"     - "sigma_Intercept"
DETAIL FAIL [a nonlinear parameter's intercept keeps its place, as in brms] Expected `rownames(summary(fit)$fixed)` to be identical to `want`. | Differences: |     actual            | expected              | [1] "a_Intercept"     - "sigma_I
DETAIL FAIL [a multivariate nonlinear parameter is named as brms names it] Expected `rownames(fixef(fit))` to be identical to `want`. | Differences: |     actual               | expected                 | [1] "a_yn_Intercept"     - "sigma_y
DETAIL FAIL [a multivariate nonlinear parameter is named as brms names it] Expected `all(c(paste0("b_", want), "sd_g__yn_a_Intercept") %in% variables(fit))` to be TRUE. | Differences: | `actual`:   FALSE | `expected`: TRUE  | 
DETAIL ERROR [a multivariate nonlinear parameter is named as brms names it] <frmtmb_error/error/condition> | Error: Some parameters cannot be found in the model:  | 'b_yn_a_Intercept' | variables() lists every name; brms's default class = "
DETAIL FAIL [drop1() and step() see an exp() term] Expected `rownames(d1)` to be identical to `c("<none>", "exp(x)", "z")`. | Differences: | `actual`:   "<none>" ".frm_exp(x)" "z" | `expected`: "<none>" "exp(x)"      "z" | 
DETAIL FAIL [drop1() and step() see an exp() term] Expected `d1["exp(x)", "Df"]` to equal 1. | Differences: |   `actual`:  NA | `expected`: 1.0 | 
DETAIL FAIL [drop1() and step() see an exp() term] Expected `abs(d1["exp(x)", "AIC"] - red)` < `1e-06 * abs(red)`. | Actual comparison: NA >= 0.0006
DETAIL FAIL [drop1() and step() see an exp() term] Expected `rownames(fixef(s))` to be identical to `c("Intercept", "expx")`. | Differences: | `actual`:   "Intercept" ".frm_expx" "z" | `expected`: "Intercept" "expx"          | 
DETAIL ERROR [re_formula can name a group-level exp() term] <frmtmb_error/error/condition> | Error: fitted(): the re_formula term (1 + exp(x) | g) matches no group-level term of this fit. A term is kept when its grouping factor is one the f
DETAIL FAIL [set_prior(coef =) takes the design column and brms's name] Expected poly(x, 2)1 to be identical to `pr[2]`. | Differences: | `actual`:   "poly(x, 2)1" | `expected`: "polyx21"     | 
DETAIL ERROR [set_prior(coef =) takes the design column and brms's name] <frmtmb_error/error/condition> | Error: Prior target not found (class=b, coef=polyx21)
DETAIL FAIL [a smooth column's two spellings are one prior slot] Expected `validate_prior(...)` to throw a error.
DETAIL FAIL [a smooth column's two spellings are one prior slot] Expected `allow_boundary(...)` to throw a error.
DETAIL ERROR [class sd and cs() take the design column's spelling too] <frmtmb_error/error/condition> | Error: No random-effect SDs match class=sd, coef=I(x^2), group=g. Per-coefficient standard deviations: group "g": "Intercept", "IxE2". A
DETAIL ERROR [a second nonlinear response may call exp()] Error in `.frm_exp(z)`: could not find function ".frm_exp"
RESULT test-brms-names-specials.R pass=8 fail=33 err=9 skip=0 warn=0
== newtest-lane
lib: frmtmb 0.74.0 C:/Users/adf44/source/r/wt-brmsnames-lib; frmtmb 0.74.0 C:/Users/adf44/source/r/wt-brmsnames-lib
RESULT test-brms-names-specials.R pass=257 fail=0 err=0 skip=0 warn=0
== newtest-sample-base
lib: frmtmb.sample 0.21.1 C:/Users/adf44/source/r/rellib-r12; frmtmb 0.74.0 C:/Users/adf44/source/r/rellib-r12
DETAIL FAIL [draws of exp(x) are named as brms names them] Expected `grep("^(b|sd|cor)_", v, value = TRUE)` to be identical to `want`. | Differences: |     actual                        | expected                     | [1] "b_Intercept"    
DETAIL FAIL [draws of exp(x) are named as brms names them] Expected `any(grepl(".frm_", v, fixed = TRUE))` to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [draws of exp(x) are named as brms names them] Expected `rownames(fixef(ds))` to be identical to `c("Intercept", "expx")`. | Differences: | `actual`:   "Intercept" ".frm_expx" | `expected`: "Intercept" "expx"      | 
DETAIL ERROR [draws of exp(x) are named as brms names them] <frmtmb_error/error/condition> | Error: Some parameters cannot be found in the model:  | 'b_expx' | variables() lists every name; brms's default class = "b" puts b_ before each nam
DETAIL ERROR [frm_sample() takes a coef in either spelling] <frmtmb_error/error/condition> | Error: Prior target not found (class=b, coef=IzE2)
DETAIL FAIL [fixef() of the draws is in brms's order] Expected `rownames(fixef(fit))` to be identical to `want`. | Differences: | `actual`:   "Intercept" "sigma_Intercept" "x" "sigma_.frm_expz" | `expected`: "Intercept" "sigma_Intercept" "x
DETAIL FAIL [fixef() of the draws is in brms's order] Expected `rownames(fixef(ds))` to be identical to `want`. | Differences: | `actual`:   "Intercept" "x"               "sigma_Intercept" "sigma_.frm_expz" | `expected`: "Intercept" "sigma_
DETAIL FAIL [fixef() of the draws is in brms's order] Expected `colnames(fixef(ds, summary = FALSE))` to be identical to `want`. | Differences: | `actual`:   "Intercept" "x"               "sigma_Intercept" "sigma_.frm_expz" | `expected`: "I
RESULT test-draws-names-specials.R pass=1 fail=6 err=2 skip=0 warn=0
== newtest-sample-lane
lib: frmtmb.sample 0.21.1 C:/Users/adf44/source/r/wt-brmsnames-lib; frmtmb 0.74.0 C:/Users/adf44/source/r/wt-brmsnames-lib
RESULT test-draws-names-specials.R pass=9 fail=0 err=0 skip=0 warn=0

<!-- END seen-fail.txt -->
```

## Names over the whole suite, before and after

`dev/brmsnames-1024-suitenames.sh` runs every test file of a package,
one process per file, with `frm()` traced in frmtmb's namespace (so
the fits frmtmb makes inside `update()` and refits count too), and
records for each fit, in call order: the internal coefficient names,
brms's names (`brms_coef_table()`), `variables()`, the outer parameter
labels (the rows of `confint()` and `vcov(full = TRUE)`) and the coef
column of the prior table. None of these computes standard errors, so
the hook changes nothing a test sees; the test counts of every file
are the same on both arms. For frmtmb.sample the hook is on
`frm_sample()` and records the draws' `variables()`. The lane's two new
test files are not paired (on 0.74.0 their blocks stop at the first
failure, so the fits do not line up by call order).

Every changed name must be accounted for by one of the lane's
changes: `alias` (the base name holds `.frm_` and the lane's is the
same without it), `prior` (the lane's coef is `brms_rename()` of the
base's), `alias+prior` (both), `prior (row order)` (the table's rows
moved to brms's order and, as a whole, are those renames), and
`mv-nl prefix` (a brms name with the same parts in brms's order,
`b_a_y_Intercept` to `b_y_a_Intercept`). No name on the lane holds
`.frm_`. In the core suite one file fits an `exp()` term
(`test-v15.R`, `exp(as.numeric(pos)/25)`) and one a multivariate
nonlinear model (`test-nlf.R`); every other change is a prior row.
Nothing is UNACCOUNTED. (A first rerun showed one: in
`test-diagnostics-ux.R` a fit whose prior table has no rows makes
`prior_table()` fail on both arms, and the new sort changed the error
message; the sort now skips an empty table, and both arms fail
alike.)

Core:

<!-- BEGIN suitediff.md -->


Test files paired: 274 

Files on one arm only: test-brms-names-specials.rds 

Fits paired: 3753 

Names compared, per field: internal 13570, brms 13570, variables 182102, outer 18663, prior 30365 

Files whose fit count differs: none 

Files whose test counts differ between arms: none 

Names holding `.frm_` on the lane: 0 

Changed names: 4743 

| field | alias | alias+prior | mv-nl prefix | prior | prior (row order) |
|---| ---|---|---|---|--- |
| brms | 1 | 0 | 2 | 0 | 0 |
| internal | 1 | 0 | 0 | 0 | 0 |
| outer | 1 | 0 | 0 | 0 | 0 |
| prior | 0 | 1 | 0 | 354 | 4380 |
| variables | 1 | 0 | 2 | 0 | 0 |

Distinct changes (406), with the files they occur in:

| why | field | base | lane | files |
|---|---|---|---|---|
| alias | brms | `b_.frm_expas.numericposD25` | `b_expas.numericposD25` | test-v15 |
| alias | internal | `.frm_exp(as.numeric(pos)/25)` | `exp(as.numeric(pos)/25)` | test-v15 |
| alias | outer | `.frm_exp(as.numeric(pos)/25)` | `exp(as.numeric(pos)/25)` | test-v15 |
| alias | variables | `b_.frm_expas.numericposD25` | `b_expas.numericposD25` | test-v15 |
| alias+prior | prior | `.frm_exp(as.numeric(pos)/25)` | `expas.numericposD25` | test-v15 |
| mv-nl prefix | brms | `b_a_y_Intercept` | `b_y_a_Intercept` | test-nlf |
| mv-nl prefix | brms | `b_b_y_Intercept` | `b_y_b_Intercept` | test-nlf |
| mv-nl prefix | variables | `b_a_y_Intercept` | `b_y_a_Intercept` | test-nlf |
| mv-nl prefix | variables | `b_b_y_Intercept` | `b_y_b_Intercept` | test-nlf |
| prior | prior | `c3(0.000157,0.201]` | `c30.0001570.201` | test-tmb-examples |
| prior | prior | `c3(0.201,0.401]` | `c30.2010.401` | test-tmb-examples |
| prior | prior | `c3(0.401,0.6]` | `c30.4010.6` | test-tmb-examples |
| prior | prior | `c3(0.6,0.8]` | `c30.60.8` | test-tmb-examples |
| prior | prior | `c3(0.8,1]` | `c30.81` | test-tmb-examples |
| prior | prior | `fc-d` | `fcMd` | test-brms-names |
| prior | prior | `I(x * kx)` | `IxMUkx` | test-aterm-expr |
| prior | prior | `I(x^1)` | `IxE1` | test-se-tier |
| prior | prior | `I(x^2)` | `IxE2` | test-brms-names, test-curve-illcond, test-methods-audit, test-refit-data, test-se-tier |
| prior | prior | `I(x^3)` | `IxE3` | test-curve-illcond, test-se-tier |
| prior | prior | `I(x^4)` | `IxE4` | test-curve-illcond, test-se-tier |
| prior | prior | `I(x^5)` | `IxE5` | test-curve-illcond, test-se-tier |
| prior | prior | `I(x^6)` | `IxE6` | test-curve-illcond |
| prior | prior | `I(x^7)` | `IxE7` | test-curve-illcond |
| prior | prior | `I(x^8)` | `IxE8` | test-curve-illcond |
| prior | prior | `I(x2^2)` | `Ix2E2` | test-refcat |
| prior | prior | `income^4` | `incomeE4` | test-case-studies |
| prior | prior | `income^5` | `incomeE5` | test-case-studies |
| prior | prior | `lf(x)` | `lfx` | test-brms-parity-defects |
| prior | prior | `log(abs(z) + 1)` | `logabszP1` | test-ce-parity, test-emmeans |
| prior | prior | `log(x^2)` | `logxE2` | test-adefects |
| prior | prior | `log(x2)` | `logx2` | test-refcat, test-refit-data, test-refit-refcat |
| prior | prior | `poly(x, 2)1` | `polyx21` | test-ce-parity, test-confint-anova, test-edgecases, test-refcat, test-sparsex |
| prior | prior | `poly(x, 2)2` | `polyx22` | test-ce-parity, test-confint-anova, test-edgecases, test-refcat, test-sparsex |
| prior | prior | `poly(x, 2, raw = TRUE)1` | `polyx2rawEQTRUE1` | test-edgecases |
| prior | prior | `poly(x, 2, raw = TRUE)2` | `polyx2rawEQTRUE2` | test-edgecases |
| prior | prior | `poly(x, 8)1` | `polyx81` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)2` | `polyx82` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)3` | `polyx83` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)4` | `polyx84` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)5` | `polyx85` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)6` | `polyx86` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)7` | `polyx87` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x, 8)8` | `polyx88` | test-optrap, test-se-flags, test-se-reparam |
| prior | prior | `poly(x2, 2)1` | `polyx221` | test-refit-data |
| prior | prior | `poly(x2, 2)2` | `polyx222` | test-refit-data |
| prior | prior | `poly(z, 2)1` | `polyz21` | test-emmeans |
| prior | prior | `poly(z, 2)2` | `polyz22` | test-emmeans |
| prior | prior | `poly(zz, 2)1` | `polyzz21` | test-emmeans |
| prior | prior | `poly(zz, 2)2` | `polyzz22` | test-emmeans |
| prior | prior | `scale(x)` | `scalex` | test-edgecases, test-refcat, test-refit-data |
| prior | prior | `scale(z)` | `scalez` | test-emmeans |
| prior | prior | `scale(z):fb` | `scalez:fb` | test-emmeans |
| prior | prior | `scale(z):fc` | `scalez:fc` | test-emmeans |
| prior | prior | `splines::ns(x2, 3)1` | `splines::nsx231` | test-refit-data, test-refit-refcat |
| prior | prior | `splines::ns(x2, 3)2` | `splines::nsx232` | test-refit-data, test-refit-refcat |
| prior | prior | `splines::ns(x2, 3)3` | `splines::nsx233` | test-refit-data, test-refit-refcat |
| prior | prior | `x:log(x2)` | `x:logx2` | test-refit-data |
| prior (row order) | prior | `` | `x` | test-adefects, test-arma-na-newdata, test-aterm-expr, test-autocor-cond, test-autocor, test-brms-names, test-brms-parity-defects, test-brms-shapes-punch2, test-cdfs, test-ce-levels, test-compat, test-corr2-alias, test-curve, test-custom-family, test-dirichlet-multinomial, test-dpar-equate, test-get-prior-route, test-glmm-poisson, test-gr-by, test-group-labels, test-guardrail1, test-hurdle-cum-thres-cs, test-hurdle-nb-zoib, test-importance, test-interop, test-lkj, test-method-residue, test-methods, test-multiple-methods, test-multiple-pooling, test-multivariate, test-mv-gaps, test-newdata-aterms, test-nl, test-nlf, test-optima, test-ordinal-disc-thres, test-ordinal-mixture, test-ordinal-thres-names, test-parity-integration, test-portability, test-pp-check-types, test-predfix, test-predict-re-uncertainty, test-prior-compat, test-prior-sd-coef, test-prior-sd-scope, test-prior-support, test-priors-autocor-classes, test-priors-bounds-grcov, test-quadrature-defects, test-re-formula-partial, test-refcat, test-review-fixes, test-review-v29, test-se-check, test-se-flags, test-se-tier, test-setprior, test-simplex-families, test-simplex-mvbf, test-simplex-rdirichlet, test-student-share, test-tre, test-us-sigma, test-v11, test-v12 |
| prior (row order) | prior | `` | `Intercept` | test-adefects, test-aliased-grouping, test-allfit-start, test-api-spellings, test-arg-refusal, test-autocor-cond, test-autoscale, test-backlog, test-bcm-data-analysis, test-bcm-esp, test-bcm-mpt, test-bcm-psychophysics, test-bcm-retention, test-boot, test-brms-families, test-brms-formula-priors, test-brms-names, test-brms-parity-defects, test-brms-port, test-brms-shapes-punch, test-brms-shapes-punch2, test-case-studies, test-cdfs, test-ce-bands, test-ce-facets, test-ce-levels, test-ce-options, test-ce-parity, test-cens-windows, test-conditional-smooths, test-confint-anova, test-covstruct-seam, test-curve-difference, test-curve-new-levels, test-curve-span, test-curve, test-data2, test-dates, test-diagnostics-ux, test-dirichlet-multinomial, test-double-bar, test-dry-run, test-edgecases, test-effects, test-emmeans, test-families-extreme, test-fit-end-hook, test-functional, test-get-prior-route, test-glmm-gaussian, test-glmm-poisson, test-gr-by, test-group-labels, test-guardrail1, test-hurdle-nb-zoib, test-id-kron, test-importance, test-influence-plot, test-input-validation, test-interop, test-lkj, test-map, test-me, test-method-residue, test-methods-audit, test-methods, test-mm, test-mo-factor-re, test-multiple-methods, test-multiple-pooling, test-multivariate, test-mv-gaps, test-naming-collisions, test-negbinomial2, test-nl-body-vars, test-nl-lexical, test-nl-rtmb-scope, test-nl, test-nlf, test-offset-grid, test-open-issues, test-optima, test-optrap, test-optstart, test-par-template, test-pooled-anova, test-portability, test-pp-check-types, test-predfix, test-predict-lp-basis, test-predict-newdata, test-predict-re-uncertainty, test-prior-compat, test-prior-sd-coef, test-prior-sd-scope, test-priors-bounds-grcov, test-ps, test-quadrature-defects, test-re-formula-partial, test-refit-data, test-reml-shapes, test-review-fixes, test-review-v25, test-review-v28, test-review-v29, test-rl-example, test-sandwich, test-scale-contract, test-se-check, test-se-conv, test-se-exact-hessian, test-se-flags, test-se-lost-re-predict, test-se-reparam-reml, test-se-reparam, test-se-tier, test-simplex-families, test-simulate-ergonomics, test-simulate-newdata, test-smooth-population, test-smooths, test-sparsex, test-special-names, test-student-share, test-sugar, test-tabular-inputs, test-tmb-examples, test-tre, test-update-pool, test-update-prior-drop, test-v07, test-v11, test-v12, test-v14, test-v15, test-v17, test-v18, test-v19, test-verbose, test-weights-unbounded, test-zi |
| prior (row order) | prior | `` | `grpb` | test-autocor |
| prior (row order) | prior | `` | `gpx` | test-brms-likelihood, test-brms-parity-defects, test-covstruct-floor, test-curve-deriv, test-curve-difference, test-curve, test-emmeans, test-fd-chain, test-gp-by, test-gp-multidim, test-se-check, test-v18 |
| prior (row order) | prior | `` | `age` | test-brms-port, test-update-prior-drop |
| prior (row order) | prior | `` | `traity1` | test-case-studies, test-id-kron, test-review-v29 |
| prior (row order) | prior | `` | `traity2` | test-case-studies, test-id-kron, test-review-v29 |
| prior (row order) | prior | `` | `sx_1` | test-case-studies, test-conditional-smooths |
| prior (row order) | prior | `` | `time` | test-case-studies |
| prior (row order) | prior | `` | `mix` | test-case-studies, test-ce-bands, test-ce-parity, test-review-fixes, test-subset-rate, test-v16, test-v17, test-v18 |
| prior (row order) | prior | `` | `fb` | test-ce-options, test-double-bar |
| prior (row order) | prior | `` | `mixidxEQg1` | test-ce-parity, test-subset-rate |
| prior (row order) | prior | `` | `w` | test-ce-parity, test-nl, test-nlf, test-subset-rate |
| prior (row order) | prior | `` | `z` | test-ce-parity, test-families-extreme, test-hurdle-cum-thres-cs, test-mv-gaps, test-nl, test-portability, test-prior-compat, test-review-v29, test-simplex-families, test-student-share, test-update-pool |
| prior (row order) | prior | `` | `2` | test-compat, test-ordinal-mixture, test-parity-integration |
| prior (row order) | prior | `` | `sz_1` | test-conditional-smooths, test-smooths |
| prior (row order) | prior | `` | `tim4` | test-covstruct |
| prior (row order) | prior | `` | `t` | test-curve-deriv, test-curve-new-levels, test-optrap |
| prior (row order) | prior | `` | `fa` | test-double-bar, test-open-issues |
| prior (row order) | prior | `` | `fc` | test-double-bar, test-formula-cmc |
| prior (row order) | prior | `` | `tim6` | test-edgecases |
| prior (row order) | prior | `` | `fC` | test-emmeans |
| prior (row order) | prior | `` | `fB` | test-emmeans |
| prior (row order) | prior | `` | `gc` | test-formula-cmc |
| prior (row order) | prior | `` | `Days` | test-glmm-gaussian, test-v14 |
| prior (row order) | prior | `` | `gpxfa` | test-gp-by |
| prior (row order) | prior | `` | `gpxw` | test-gp-by |
| prior (row order) | prior | `` | `gpxfIntercept` | test-gp-by |
| prior (row order) | prior | `` | `gpxzfax` | test-gp-by |
| prior (row order) | prior | `` | `gpxzfa` | test-gp-by |
| prior (row order) | prior | `` | `gpx1x2x1` | test-gp-multidim |
| prior (row order) | prior | `` | `gpx1x2` | test-gp-multidim |
| prior (row order) | prior | `` | `t4` | test-lkj |
| prior (row order) | prior | `` | `mixm` | test-mi-bsp-names |
| prior (row order) | prior | `` | `mixm:z` | test-mi-bsp-names |
| prior (row order) | prior | `` | `mmcc1c2` | test-mm |
| prior (row order) | prior | `` | `mox` | test-mo-factor-re |
| prior (row order) | prior | `` | `mox:z` | test-mo-factor-re |
| prior (row order) | prior | `` | `mix:z` | test-mo-terms, test-v16, test-v18 |
| prior (row order) | prior | `` | `moinc` | test-mo-terms |
| prior (row order) | prior | `` | `moinc:z` | test-mo-terms |
| prior (row order) | prior | `` | `moinc1` | test-mo-terms |
| prior (row order) | prior | `` | `1` | test-mv-gaps, test-ordinal-mixture |
| prior (row order) | prior | `` | `g4` | test-naming-collisions, test-skew-normal-start |
| prior (row order) | prior | `` | `group` | test-nl |
| prior (row order) | prior | `` | `f2` | test-offset-grid |
| prior (row order) | prior | `` | `3` | test-ordinal-mixture, test-parity-integration |
| prior (row order) | prior | `` | `xs` | test-predfix |
| prior (row order) | prior | `` | `sex` | test-ps |
| prior (row order) | prior | `` | `logx2` | test-refit-data |
| prior (row order) | prior | `` | `f9` | test-se-check, test-se-tier |
| prior (row order) | prior | `` | `g9` | test-se-check |
| prior (row order) | prior | `` | `z2` | test-se-flags |
| prior (row order) | prior | `` | `x2` | test-simulate-newdata |
| prior (row order) | prior | `` | `g2` | test-skew-normal-start |
| prior (row order) | prior | `` | `g1` | test-skew-normal-start |
| prior (row order) | prior | `` | `v` | test-special-names |
| prior (row order) | prior | `` | `mo` | test-special-names |
| prior (row order) | prior | `` | `mi` | test-special-names |
| prior (row order) | prior | `` | `ar` | test-special-names |
| prior (row order) | prior | `` | `ma` | test-special-names |
| prior (row order) | prior | `` | `arma` | test-special-names |
| prior (row order) | prior | `` | `cosy` | test-special-names |
| prior (row order) | prior | `` | `unstr` | test-special-names |
| prior (row order) | prior | `` | `mixidxEQg2` | test-subset-rate |
| prior (row order) | prior | `` | `tim1` | test-v11 |
| prior (row order) | prior | `` | `tim5` | test-v15 |
| prior (row order) | prior | `` | `spp2` | test-v16, test-v17 |
| prior (row order) | prior | `` | `spp3` | test-v16, test-v17 |
| prior (row order) | prior | `` | `faca` | test-weights-unbounded |
| prior (row order) | prior | `1` | `z` | test-compat, test-ordinal-mixture |
| prior (row order) | prior | `1` | `` | test-mv-gaps, test-ordinal-mixture, test-parity-integration |
| prior (row order) | prior | `1` | `3` | test-ordinal-mixture |
| prior (row order) | prior | `2` | `` | test-compat, test-mv-gaps, test-ordinal-mixture, test-parity-integration |
| prior (row order) | prior | `2` | `1` | test-ordinal-mixture |
| prior (row order) | prior | `3` | `1` | test-compat, test-ordinal-mixture |
| prior (row order) | prior | `3` | `z` | test-ordinal-mixture |
| prior (row order) | prior | `3` | `2` | test-ordinal-mixture |
| prior (row order) | prior | `3` | `` | test-parity-integration |
| prior (row order) | prior | `age` | `Intercept` | test-brms-port, test-update-prior-drop |
| prior (row order) | prior | `ar` | `` | test-special-names |
| prior (row order) | prior | `ar` | `Intercept` | test-special-names |
| prior (row order) | prior | `arma` | `` | test-special-names |
| prior (row order) | prior | `arma` | `Intercept` | test-special-names |
| prior (row order) | prior | `c1` | `` | test-mm |
| prior (row order) | prior | `cosy` | `` | test-special-names |
| prior (row order) | prior | `cosy` | `Intercept` | test-special-names |
| prior (row order) | prior | `Days` | `` | test-boot, test-ce-bands, test-confint-anova, test-glmm-gaussian, test-interop, test-lkj, test-methods, test-se-check, test-se-tier, test-sugar, test-v07, test-v11, test-v14, test-v19, test-verbose |
| prior (row order) | prior | `Days` | `Intercept` | test-glmm-gaussian, test-v14 |
| prior (row order) | prior | `f1` | `` | test-se-check, test-se-tier |
| prior (row order) | prior | `f10` | `f1` | test-se-check, test-se-tier |
| prior (row order) | prior | `f11` | `f10` | test-se-check |
| prior (row order) | prior | `f12` | `f11` | test-se-check |
| prior (row order) | prior | `f13` | `f12` | test-se-check |
| prior (row order) | prior | `f14` | `f13` | test-se-check |
| prior (row order) | prior | `f15` | `f14` | test-se-check |
| prior (row order) | prior | `f16` | `f15` | test-se-check |
| prior (row order) | prior | `f17` | `f16` | test-se-check |
| prior (row order) | prior | `f18` | `f17` | test-se-check |
| prior (row order) | prior | `f19` | `f18` | test-se-check |
| prior (row order) | prior | `f2` | `Intercept` | test-offset-grid |
| prior (row order) | prior | `f2` | `f19` | test-se-check |
| prior (row order) | prior | `f2` | `f10` | test-se-check, test-se-tier |
| prior (row order) | prior | `f20` | `f2` | test-se-check |
| prior (row order) | prior | `f21` | `f20` | test-se-check |
| prior (row order) | prior | `f22` | `f21` | test-se-check |
| prior (row order) | prior | `f23` | `f22` | test-se-check |
| prior (row order) | prior | `f24` | `f23` | test-se-check |
| prior (row order) | prior | `f25` | `f24` | test-se-check |
| prior (row order) | prior | `f26` | `f25` | test-se-check |
| prior (row order) | prior | `f27` | `f26` | test-se-check |
| prior (row order) | prior | `f28` | `f27` | test-se-check |
| prior (row order) | prior | `f29` | `f28` | test-se-check |
| prior (row order) | prior | `f3` | `f29` | test-se-check |
| prior (row order) | prior | `f3` | `f20` | test-se-check |
| prior (row order) | prior | `f3` | `f2` | test-se-check, test-se-tier |
| prior (row order) | prior | `f30` | `f3` | test-se-check |
| prior (row order) | prior | `f31` | `f30` | test-se-check |
| prior (row order) | prior | `f32` | `f31` | test-se-check |
| prior (row order) | prior | `f33` | `f32` | test-se-check |
| prior (row order) | prior | `f34` | `f33` | test-se-check |
| prior (row order) | prior | `f35` | `f34` | test-se-check |
| prior (row order) | prior | `f36` | `f35` | test-se-check |
| prior (row order) | prior | `f37` | `f36` | test-se-check |
| prior (row order) | prior | `f38` | `f37` | test-se-check |
| prior (row order) | prior | `f39` | `f38` | test-se-check |
| prior (row order) | prior | `f4` | `f39` | test-se-check |
| prior (row order) | prior | `f4` | `f3` | test-se-check, test-se-tier |
| prior (row order) | prior | `f40` | `f4` | test-se-check |
| prior (row order) | prior | `f41` | `f40` | test-se-check |
| prior (row order) | prior | `f42` | `f41` | test-se-check |
| prior (row order) | prior | `f43` | `f42` | test-se-check |
| prior (row order) | prior | `f44` | `f43` | test-se-check |
| prior (row order) | prior | `f45` | `f44` | test-se-check |
| prior (row order) | prior | `f46` | `f45` | test-se-check |
| prior (row order) | prior | `f47` | `f46` | test-se-check |
| prior (row order) | prior | `f48` | `f47` | test-se-check |
| prior (row order) | prior | `f49` | `f48` | test-se-check |
| prior (row order) | prior | `f5` | `f49` | test-se-check |
| prior (row order) | prior | `f5` | `f4` | test-se-check, test-se-tier |
| prior (row order) | prior | `f50` | `f5` | test-se-check |
| prior (row order) | prior | `f51` | `f50` | test-se-check |
| prior (row order) | prior | `f52` | `f51` | test-se-check |
| prior (row order) | prior | `f53` | `f52` | test-se-check |
| prior (row order) | prior | `f54` | `f53` | test-se-check |
| prior (row order) | prior | `f55` | `f54` | test-se-check |
| prior (row order) | prior | `f56` | `f55` | test-se-check |
| prior (row order) | prior | `f57` | `f56` | test-se-check |
| prior (row order) | prior | `f58` | `f57` | test-se-check |
| prior (row order) | prior | `f59` | `f58` | test-se-check |
| prior (row order) | prior | `f6` | `f59` | test-se-check |
| prior (row order) | prior | `f6` | `f5` | test-se-check, test-se-tier |
| prior (row order) | prior | `f60` | `f6` | test-se-check |
| prior (row order) | prior | `f7` | `f60` | test-se-check |
| prior (row order) | prior | `f7` | `f6` | test-se-check, test-se-tier |
| prior (row order) | prior | `f8` | `f7` | test-se-check, test-se-tier |
| prior (row order) | prior | `f9` | `f8` | test-se-check, test-se-tier |
| prior (row order) | prior | `fa` | `fb` | test-double-bar, test-open-issues |
| prior (row order) | prior | `fa:x` | `fb:x` | test-double-bar |
| prior (row order) | prior | `faca` | `` | test-weights-unbounded |
| prior (row order) | prior | `facb` | `Intercept` | test-weights-unbounded |
| prior (row order) | prior | `fb` | `Intercept` | test-ce-options, test-double-bar |
| prior (row order) | prior | `fb` | `fc` | test-double-bar, test-mo-factor-re, test-open-issues |
| prior (row order) | prior | `fb` | `` | test-double-bar, test-formula-cmc |
| prior (row order) | prior | `fb` | `x` | test-double-bar |
| prior (row order) | prior | `fb` | `fb:mox` | test-mo-factor-re |
| prior (row order) | prior | `fB` | `Intercept` | test-emmeans |
| prior (row order) | prior | `fB` | `fC` | test-emmeans |
| prior (row order) | prior | `fb:mox` | `fc` | test-mo-factor-re |
| prior (row order) | prior | `fb:x` | `fc:x` | test-double-bar |
| prior (row order) | prior | `fc` | `` | test-double-bar, test-open-issues |
| prior (row order) | prior | `fc` | `z` | test-double-bar |
| prior (row order) | prior | `fc` | `x` | test-double-bar |
| prior (row order) | prior | `fc` | `Intercept` | test-double-bar |
| prior (row order) | prior | `fc` | `fb` | test-double-bar, test-formula-cmc |
| prior (row order) | prior | `fc` | `mox` | test-mo-factor-re |
| prior (row order) | prior | `fc` | `fc:mox` | test-mo-factor-re |
| prior (row order) | prior | `fC` | `fB` | test-emmeans |
| prior (row order) | prior | `fC` | `` | test-emmeans |
| prior (row order) | prior | `fc:mox` | `mox` | test-mo-factor-re |
| prior (row order) | prior | `fc:x` | `` | test-double-bar |
| prior (row order) | prior | `fcq` | `` | test-mm |
| prior (row order) | prior | `g1` | `` | test-se-check, test-skew-normal-start |
| prior (row order) | prior | `g1` | `g2` | test-skew-normal-start |
| prior (row order) | prior | `g10` | `g1` | test-se-check |
| prior (row order) | prior | `g11` | `g10` | test-se-check |
| prior (row order) | prior | `g12` | `g11` | test-se-check |
| prior (row order) | prior | `g13` | `g12` | test-se-check |
| prior (row order) | prior | `g14` | `g13` | test-se-check |
| prior (row order) | prior | `g15` | `g14` | test-se-check |
| prior (row order) | prior | `g16` | `g15` | test-se-check |
| prior (row order) | prior | `g17` | `g16` | test-se-check |
| prior (row order) | prior | `g18` | `g17` | test-se-check |
| prior (row order) | prior | `g19` | `g18` | test-se-check |
| prior (row order) | prior | `g2` | `Intercept` | test-naming-collisions |
| prior (row order) | prior | `g2` | `g19` | test-se-check |
| prior (row order) | prior | `g2` | `g3` | test-skew-normal-start |
| prior (row order) | prior | `g2` | `g1` | test-skew-normal-start |
| prior (row order) | prior | `g20` | `g2` | test-se-check |
| prior (row order) | prior | `g3` | `g2` | test-naming-collisions, test-skew-normal-start |
| prior (row order) | prior | `g3` | `g20` | test-se-check |
| prior (row order) | prior | `g3` | `g4` | test-skew-normal-start |
| prior (row order) | prior | `g4` | `g3` | test-naming-collisions, test-se-check, test-skew-normal-start |
| prior (row order) | prior | `g4` | `` | test-skew-normal-start |
| prior (row order) | prior | `g5` | `g4` | test-se-check |
| prior (row order) | prior | `g6` | `g5` | test-se-check |
| prior (row order) | prior | `g7` | `g6` | test-se-check |
| prior (row order) | prior | `g8` | `g7` | test-se-check |
| prior (row order) | prior | `g9` | `g8` | test-se-check |
| prior (row order) | prior | `ga` | `` | test-ps |
| prior (row order) | prior | `gb` | `` | test-formula-cmc |
| prior (row order) | prior | `gc` | `gb` | test-formula-cmc |
| prior (row order) | prior | `gpx` | `` | test-brms-likelihood, test-brms-parity-defects, test-covstruct-floor, test-curve-deriv, test-curve-difference, test-curve, test-emmeans, test-gp-by, test-gp-multidim, test-se-check, test-v18 |
| prior (row order) | prior | `gpx` | `Intercept` | test-fd-chain |
| prior (row order) | prior | `gpx1x2` | `` | test-gp-multidim |
| prior (row order) | prior | `gpx1x2x1` | `gpx1x2x2` | test-gp-multidim |
| prior (row order) | prior | `gpx1x2x2` | `` | test-gp-multidim |
| prior (row order) | prior | `gpxfa` | `gpxfb` | test-gp-by |
| prior (row order) | prior | `gpxfb` | `gpxfc` | test-gp-by |
| prior (row order) | prior | `gpxfc` | `` | test-gp-by |
| prior (row order) | prior | `gpxfIntercept` | `gpxfb` | test-gp-by |
| prior (row order) | prior | `gpxw` | `` | test-gp-by |
| prior (row order) | prior | `gpxzfa` | `gpxzfb` | test-gp-by |
| prior (row order) | prior | `gpxzfax` | `gpxzfaz` | test-gp-by |
| prior (row order) | prior | `gpxzfaz` | `gpxzfbx` | test-gp-by |
| prior (row order) | prior | `gpxzfb` | `gpxzfc` | test-gp-by |
| prior (row order) | prior | `gpxzfbx` | `gpxzfbz` | test-gp-by |
| prior (row order) | prior | `gpxzfbz` | `gpxzfcx` | test-gp-by |
| prior (row order) | prior | `gpxzfc` | `` | test-gp-by |
| prior (row order) | prior | `gpxzfcx` | `gpxzfcz` | test-gp-by |
| prior (row order) | prior | `gpxzfcz` | `` | test-gp-by |
| prior (row order) | prior | `group` | `Intercept` | test-methods, test-nl |
| prior (row order) | prior | `grpa` | `` | test-autocor |
| prior (row order) | prior | `grpb` | `grpa` | test-autocor |
| prior (row order) | prior | `Intercept` | `` | test-adefects, test-aliased-grouping, test-allfit-start, test-api-spellings, test-arg-refusal, test-autocor-cond, test-autoscale, test-backlog, test-bcm-data-analysis, test-bcm-esp, test-bcm-mpt, test-bcm-retention, test-boot, test-brms-families, test-brms-formula-priors, test-brms-names, test-brms-parity-defects, test-brms-port, test-brms-shapes-punch, test-brms-shapes-punch2, test-case-studies, test-cdfs, test-ce-bands, test-ce-facets, test-ce-levels, test-ce-options, test-ce-parity, test-conditional-smooths, test-confint-anova, test-covstruct-seam, test-curve-deriv, test-curve-difference, test-curve-new-levels, test-curve-span, test-curve, test-data2, test-dates, test-diagnostics-ux, test-double-bar, test-dry-run, test-edgecases, test-effects, test-emmeans, test-families-extreme, test-fd-chain, test-fit-end-hook, test-functional, test-get-prior-route, test-glmm-gaussian, test-glmm-poisson, test-gr-by, test-group-labels, test-guardrail1, test-id-kron, test-importance, test-influence-plot, test-input-validation, test-interop, test-lkj, test-map, test-me, test-method-residue, test-methods-audit, test-methods, test-mm, test-mo-factor-re, test-multiple-methods, test-multiple-pooling, test-multivariate, test-mv-gaps, test-naming-collisions, test-negbinomial2, test-nl-body-vars, test-nl-lexical, test-nl-rtmb-scope, test-nl, test-nlf, test-offset-grid, test-open-issues, test-optima, test-optrap, test-optstart, test-par-template, test-pooled-anova, test-portability, test-pp-check-types, test-predfix, test-predict-lp-basis, test-predict-newdata, test-predict-re-uncertainty, test-prior-compat, test-prior-sd-coef, test-prior-sd-scope, test-priors-bounds-grcov, test-ps, test-quadrature-defects, test-re-formula-partial, test-refit-data, test-reml-shapes, test-review-fixes, test-review-v25, test-review-v28, test-review-v29, test-rl-example, test-sandwich, test-scale-contract, test-se-check, test-se-conv, test-se-exact-hessian, test-se-flags, test-se-lost-re-predict, test-se-reparam-reml, test-se-reparam, test-se-tier, test-simplex-families, test-simulate-ergonomics, test-simulate-newdata, test-smooth-population, test-smooths, test-sparsex, test-special-names, test-student-share, test-sugar, test-tabular-inputs, test-tmb-examples, test-tre, test-update-pool, test-update-prior-drop, test-us-sigma, test-v11, test-v12, test-v14, test-v15, test-v17, test-v18, test-v19, test-verbose, test-weights-unbounded, test-zi |
| prior (row order) | prior | `Intercept` | `xbig` | test-autoscale |
| prior (row order) | prior | `Intercept` | `xz` | test-autoscale |
| prior (row order) | prior | `Intercept` | `xc` | test-bcm-psychophysics |
| prior (row order) | prior | `Intercept` | `Days` | test-boot, test-ce-bands, test-confint-anova, test-glmm-gaussian, test-interop, test-lkj, test-methods, test-se-check, test-se-tier, test-sugar, test-v07, test-v11, test-v14, test-v19, test-verbose |
| prior (row order) | prior | `Intercept` | `x` | test-brms-formula-priors, test-brms-parity-defects, test-cens-windows, test-dirichlet-multinomial, test-double-bar, test-famgaps, test-gr-by, test-group-labels, test-methods, test-nlf, test-open-issues, test-predict-newdata, test-prior-compat, test-prior-sd-scope, test-scale-contract, test-setprior, test-simplex-families, test-sparsex, test-student-share, test-tre, test-xbeta-zibb-hurdle-cum |
| prior (row order) | prior | `Intercept` | `z` | test-brms-names, test-brms-shapes, test-distributional, test-hurdle-nb-zoib, test-method-residue, test-student-share |
| prior (row order) | prior | `Intercept` | `time` | test-case-studies |
| prior (row order) | prior | `Intercept` | `fb` | test-double-bar, test-mo-factor-re, test-open-issues |
| prior (row order) | prior | `Intercept` | `fa:x` | test-double-bar |
| prior (row order) | prior | `Intercept` | `fc` | test-double-bar |
| prior (row order) | prior | `Intercept` | `fa` | test-double-bar, test-open-issues |
| prior (row order) | prior | `Intercept` | `t2` | test-lkj |
| prior (row order) | prior | `Intercept` | `group` | test-methods, test-nl |
| prior (row order) | prior | `Intercept` | `c1` | test-mm |
| prior (row order) | prior | `Intercept` | `fcq` | test-mm |
| prior (row order) | prior | `Intercept` | `ga` | test-ps |
| prior (row order) | prior | `Intercept` | `v` | test-special-names |
| prior (row order) | prior | `Intercept` | `mo` | test-special-names |
| prior (row order) | prior | `Intercept` | `mi` | test-special-names |
| prior (row order) | prior | `Intercept` | `ar` | test-special-names |
| prior (row order) | prior | `Intercept` | `ma` | test-special-names |
| prior (row order) | prior | `Intercept` | `arma` | test-special-names |
| prior (row order) | prior | `Intercept` | `cosy` | test-special-names |
| prior (row order) | prior | `Intercept` | `unstr` | test-special-names |
| prior (row order) | prior | `Intercept` | `xs` | test-weights-unbounded |
| prior (row order) | prior | `Intercept` | `facb` | test-weights-unbounded |
| prior (row order) | prior | `logx2` | `Intercept` | test-refit-data |
| prior (row order) | prior | `ma` | `` | test-special-names |
| prior (row order) | prior | `ma` | `Intercept` | test-special-names |
| prior (row order) | prior | `mi` | `` | test-special-names |
| prior (row order) | prior | `mi` | `Intercept` | test-special-names |
| prior (row order) | prior | `mix` | `z` | test-case-studies, test-mo-terms, test-review-fixes, test-v16, test-v17, test-v18 |
| prior (row order) | prior | `mix` | `` | test-ce-bands, test-subset-rate |
| prior (row order) | prior | `mix` | `w` | test-ce-parity |
| prior (row order) | prior | `mix:z` | `` | test-mo-terms, test-v16, test-v18 |
| prior (row order) | prior | `mixidxEQg1` | `` | test-ce-parity, test-subset-rate |
| prior (row order) | prior | `mixidxEQg2` | `` | test-subset-rate |
| prior (row order) | prior | `mixm` | `z` | test-mi-bsp-names |
| prior (row order) | prior | `mixm:z` | `` | test-mi-bsp-names |
| prior (row order) | prior | `mixmi` | `z` | test-brms-parity-defects |
| prior (row order) | prior | `mmcc1c2` | `Intercept` | test-mm |
| prior (row order) | prior | `mo` | `` | test-special-names |
| prior (row order) | prior | `mo` | `Intercept` | test-special-names |
| prior (row order) | prior | `moinc` | `` | test-mo-terms |
| prior (row order) | prior | `moinc:z` | `` | test-mo-terms |
| prior (row order) | prior | `moinc1` | `mow` | test-mo-terms |
| prior (row order) | prior | `mow` | `` | test-mo-terms |
| prior (row order) | prior | `mox` | `Intercept` | test-mo-factor-re |
| prior (row order) | prior | `mox` | `mox:fb` | test-mo-factor-re |
| prior (row order) | prior | `mox` | `` | test-mo-factor-re |
| prior (row order) | prior | `mox:fb` | `mox:fc` | test-mo-factor-re |
| prior (row order) | prior | `mox:fc` | `` | test-mo-factor-re |
| prior (row order) | prior | `mox:w1` | `mox:z1` | test-mo-factor-re |
| prior (row order) | prior | `mox:z` | `` | test-mo-factor-re |
| prior (row order) | prior | `mox:z1` | `mox:w1` | test-mo-factor-re |
| prior (row order) | prior | `sex` | `Intercept` | test-ps |
| prior (row order) | prior | `spp1` | `` | test-v16, test-v17 |
| prior (row order) | prior | `spp2` | `` | test-v16, test-v17 |
| prior (row order) | prior | `spp3` | `spp1` | test-v16, test-v17 |
| prior (row order) | prior | `sx_1` | `` | test-case-studies |
| prior (row order) | prior | `sx_1` | `Intercept` | test-conditional-smooths |
| prior (row order) | prior | `sz_1` | `` | test-conditional-smooths, test-smooths |
| prior (row order) | prior | `sz_1` | `x` | test-nlf |
| prior (row order) | prior | `t` | `Intercept` | test-curve-deriv, test-curve-new-levels, test-optrap |
| prior (row order) | prior | `t1` | `` | test-lkj |
| prior (row order) | prior | `t2` | `t3` | test-lkj |
| prior (row order) | prior | `t2` | `t1` | test-lkj |
| prior (row order) | prior | `t3` | `t4` | test-lkj |
| prior (row order) | prior | `t3` | `t2` | test-lkj |
| prior (row order) | prior | `t4` | `` | test-lkj |
| prior (row order) | prior | `t4` | `t3` | test-lkj |
| prior (row order) | prior | `tim1` | `` | test-covstruct, test-edgecases, test-v15 |
| prior (row order) | prior | `tim1` | `tim2` | test-v11 |
| prior (row order) | prior | `tim2` | `tim1` | test-covstruct, test-edgecases, test-v15 |
| prior (row order) | prior | `tim2` | `tim3` | test-v11 |
| prior (row order) | prior | `tim3` | `tim2` | test-covstruct, test-edgecases, test-v15 |
| prior (row order) | prior | `tim3` | `tim4` | test-v11 |
| prior (row order) | prior | `tim4` | `tim3` | test-covstruct, test-edgecases, test-v15 |
| prior (row order) | prior | `tim4` | `` | test-v11 |
| prior (row order) | prior | `tim5` | `tim4` | test-edgecases, test-v15 |
| prior (row order) | prior | `tim6` | `tim5` | test-edgecases |
| prior (row order) | prior | `time` | `` | test-case-studies |
| prior (row order) | prior | `traity1` | `` | test-case-studies, test-id-kron, test-review-v29 |
| prior (row order) | prior | `traity1` | `traity2` | test-case-studies, test-id-kron, test-review-v29 |
| prior (row order) | prior | `traity2` | `` | test-case-studies, test-id-kron, test-review-v29 |
| prior (row order) | prior | `unstr` | `` | test-special-names |
| prior (row order) | prior | `unstr` | `Intercept` | test-special-names |
| prior (row order) | prior | `v` | `` | test-special-names |
| prior (row order) | prior | `v` | `Intercept` | test-special-names |
| prior (row order) | prior | `w` | `` | test-ce-parity, test-nl, test-nlf, test-subset-rate |
| prior (row order) | prior | `x` | `` | test-adefects, test-arma-na-newdata, test-aterm-expr, test-autocor-cond, test-autocor, test-brms-formula-priors, test-brms-names, test-brms-parity-defects, test-brms-shapes-punch2, test-cdfs, test-cens-windows, test-compat, test-corr2-alias, test-custom-family, test-dirichlet-multinomial, test-double-bar, test-dpar-equate, test-get-prior-route, test-glmm-poisson, test-gr-by, test-group-labels, test-hurdle-cum-thres-cs, test-hurdle-nb-zoib, test-importance, test-interop, test-method-residue, test-methods, test-multivariate, test-mv-gaps, test-newdata-aterms, test-nl, test-nlf, test-open-issues, test-optima, test-ordinal-disc-thres, test-ordinal-mixture, test-ordinal-thres-names, test-parity-integration, test-portability, test-pp-check-types, test-predfix, test-predict-re-uncertainty, test-prior-compat, test-prior-sd-scope, test-prior-support, test-priors-autocor-classes, test-quadrature-defects, test-refcat, test-review-fixes, test-review-v29, test-setprior, test-simplex-families, test-simplex-mvbf, test-simplex-rdirichlet, test-student-share, test-tre |
| prior (row order) | prior | `x` | `Intercept` | test-brms-names, test-brms-parity-defects, test-ce-levels, test-curve, test-famgaps, test-gr-by, test-group-labels, test-guardrail1, test-lkj, test-methods, test-multiple-methods, test-multiple-pooling, test-nl, test-nlf, test-portability, test-predfix, test-predict-newdata, test-prior-compat, test-prior-sd-coef, test-prior-sd-scope, test-priors-bounds-grcov, test-re-formula-partial, test-review-v29, test-scale-contract, test-se-check, test-se-flags, test-se-tier, test-setprior, test-sparsex, test-student-share, test-tre, test-us-sigma, test-v11, test-v12, test-xbeta-zibb-hurdle-cum |
| prior (row order) | prior | `x` | `3` | test-compat, test-ordinal-mixture |
| prior (row order) | prior | `x` | `x:fb` | test-double-bar |
| prior (row order) | prior | `x` | `fb` | test-double-bar |
| prior (row order) | prior | `x` | `2` | test-mv-gaps |
| prior (row order) | prior | `x` | `sz_1` | test-nlf |
| prior (row order) | prior | `x` | `z` | test-ordinal-disc-thres, test-simplex-families |
| prior (row order) | prior | `x` | `1` | test-ordinal-mixture, test-parity-integration |
| prior (row order) | prior | `x:fb` | `x:fc` | test-double-bar |
| prior (row order) | prior | `x:fc` | `` | test-double-bar |
| prior (row order) | prior | `x1` | `Intercept` | test-simulate-newdata |
| prior (row order) | prior | `x2` | `x1` | test-simulate-newdata |
| prior (row order) | prior | `xbig` | `Intercept` | test-autoscale |
| prior (row order) | prior | `xc` | `` | test-bcm-psychophysics |
| prior (row order) | prior | `xs` | `` | test-predfix, test-weights-unbounded |
| prior (row order) | prior | `xz` | `Intercept` | test-autoscale |
| prior (row order) | prior | `z` | `Intercept` | test-brms-names, test-brms-shapes, test-ce-parity, test-distributional, test-method-residue, test-update-pool |
| prior (row order) | prior | `z` | `mixmi` | test-brms-parity-defects |
| prior (row order) | prior | `z` | `` | test-case-studies, test-double-bar, test-families-extreme, test-hurdle-cum-thres-cs, test-hurdle-nb-zoib, test-mi-bsp-names, test-mv-gaps, test-nl, test-ordinal-disc-thres, test-ordinal-mixture, test-prior-compat, test-review-fixes, test-simplex-families, test-student-share, test-v16, test-v17, test-v18 |
| prior (row order) | prior | `z` | `x` | test-compat, test-nl, test-ordinal-mixture, test-portability, test-review-v29, test-simplex-families, test-student-share |
| prior (row order) | prior | `z` | `mix` | test-mo-terms |
| prior (row order) | prior | `z1` | `x` | test-se-flags |
| prior (row order) | prior | `z2` | `z1` | test-se-flags |

<!-- END suitediff.md -->

frmtmb.sample:

<!-- BEGIN suitediff-sample.md -->


Test files paired: 65 

Files on one arm only: test-draws-names-specials.rds 

Fits paired: 637 

Names compared, per field: internal 2086, brms 2086, variables 4640, outer 3006, prior 5130 

Files whose fit count differs: none 

Files whose test counts differ between arms: none 

Names holding `.frm_` on the lane: 0 

Changed names: 1175 

| field | prior (row order) |
|---| --- |
| prior | 1175 |

Distinct changes (38), with the files they occur in:

| why | field | base | lane | files |
|---|---|---|---|---|
| prior (row order) | prior | `` | `fb` | test-add-criterion, test-adefects, test-bracket-access, test-brms-output, test-brms-pins, test-brms-shapes-draws, test-brms-suite-brmsformula, test-brms-suite-brmsterms, test-brms-suite-families, test-brms-suite-helper-copy, test-brms-suite-methods, test-compat-preflight, test-conditional-effects-draws, test-conditions-census, test-conditions, test-cox-draws, test-default-priors-bounds, test-default-priors-brms, test-default-priors-dpar, test-default-priors-fit, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-formrobust-draws, test-generic-collision, test-gev-tmp-xi, test-gp-by-draws, test-gr-by-draws, test-incl-thres-draws, test-laplace-draws, test-linpred-stack, test-loo-predict, test-loo, test-message-uniqueness, test-mo-simplex-draws, test-nb2dm-draws, test-newdata-aterms-draws, test-ordinal-disc-thres-draws, test-ordinal-draws-names, test-ordinal-mixture-draws, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-prior-route, test-prior-update, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-scale, test-simplex-draws, test-simplex-family-notes, test-simplex-mvbf-draws, test-simulators, test-sratio-draws, test-stan-control, test-subset-rate-draws, test-thres-draws, test-tmbstan-build-guard, test-weights-unbounded, test-xbeta-zibb-hurdle-draws |
| prior (row order) | prior | `` | `x` | test-add-criterion, test-adefects, test-bracket-access, test-brms-output, test-brms-pins, test-brms-shapes-draws, test-brms-suite-brmsformula, test-brms-suite-brmsterms, test-brms-suite-families, test-brms-suite-helper-copy, test-brms-suite-methods, test-compat-preflight, test-conditional-effects-draws, test-conditions-census, test-conditions, test-cox-draws, test-default-priors-bounds, test-default-priors-brms, test-default-priors-dpar, test-default-priors-fit, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-formrobust-draws, test-generic-collision, test-gev-tmp-xi, test-gp-by-draws, test-gr-by-draws, test-incl-thres-draws, test-laplace-draws, test-linpred-stack, test-loo-predict, test-loo, test-message-uniqueness, test-mo-simplex-draws, test-nb2dm-draws, test-newdata-aterms-draws, test-ordinal-disc-thres-draws, test-ordinal-draws-names, test-ordinal-mixture-draws, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-prior-route, test-prior-update, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-scale, test-simplex-draws, test-simplex-family-notes, test-simplex-mvbf-draws, test-simulators, test-sratio-draws, test-stan-control, test-subset-rate-draws, test-thres-draws, test-tmbstan-build-guard, test-weights-unbounded, test-xbeta-zibb-hurdle-draws |
| prior (row order) | prior | `` | `Intercept` | test-adefects, test-brms-pins, test-brms-shapes-draws, test-compat-preflight, test-conditional-effects-draws, test-default-priors-bounds, test-default-priors-brms, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-laplace-draws, test-loo-predict, test-loo, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-stan-control, test-weights-unbounded |
| prior (row order) | prior | `` | `z` | test-brms-pins, test-brms-shapes-draws |
| prior (row order) | prior | `` | `gpx` | test-default-priors-brms, test-gp-by-draws, test-reparam |
| prior (row order) | prior | `` | `gpxfa` | test-gp-by-draws |
| prior (row order) | prior | `` | `mix` | test-laplace-draws, test-loo |
| prior (row order) | prior | `` | `2` | test-ordinal-mixture-draws |
| prior (row order) | prior | `` | `1` | test-ordinal-mixture-draws |
| prior (row order) | prior | `` | `3` | test-ordinal-mixture-draws |
| prior (row order) | prior | `1` | `z` | test-ordinal-mixture-draws |
| prior (row order) | prior | `1` | `` | test-ordinal-mixture-draws |
| prior (row order) | prior | `2` | `` | test-ordinal-mixture-draws |
| prior (row order) | prior | `2` | `z` | test-ordinal-mixture-draws |
| prior (row order) | prior | `2` | `x` | test-ordinal-mixture-draws |
| prior (row order) | prior | `3` | `1` | test-ordinal-mixture-draws |
| prior (row order) | prior | `3` | `` | test-ordinal-mixture-draws |
| prior (row order) | prior | `fb` | `x` | test-add-criterion, test-adefects, test-bracket-access, test-brms-output, test-brms-pins, test-brms-shapes-draws, test-brms-suite-brmsformula, test-brms-suite-brmsterms, test-brms-suite-families, test-brms-suite-helper-copy, test-brms-suite-methods, test-compat-preflight, test-conditional-effects-draws, test-conditions-census, test-conditions, test-cox-draws, test-default-priors-bounds, test-default-priors-brms, test-default-priors-dpar, test-default-priors-fit, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-formrobust-draws, test-generic-collision, test-gev-tmp-xi, test-gp-by-draws, test-gr-by-draws, test-incl-thres-draws, test-laplace-draws, test-linpred-stack, test-loo-predict, test-loo, test-message-uniqueness, test-mo-simplex-draws, test-nb2dm-draws, test-newdata-aterms-draws, test-ordinal-disc-thres-draws, test-ordinal-draws-names, test-ordinal-mixture-draws, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-prior-route, test-prior-update, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-scale, test-simplex-draws, test-simplex-family-notes, test-simplex-mvbf-draws, test-simulators, test-sratio-draws, test-stan-control, test-subset-rate-draws, test-thres-draws, test-tmbstan-build-guard, test-weights-unbounded, test-xbeta-zibb-hurdle-draws |
| prior (row order) | prior | `gpx` | `` | test-default-priors-brms, test-gp-by-draws, test-reparam |
| prior (row order) | prior | `gpxfa` | `gpxfb` | test-gp-by-draws |
| prior (row order) | prior | `gpxfb` | `` | test-gp-by-draws |
| prior (row order) | prior | `Intercept` | `` | test-add-criterion, test-adefects, test-bracket-access, test-brms-output, test-brms-pins, test-brms-shapes-draws, test-brms-suite-brmsformula, test-brms-suite-brmsterms, test-brms-suite-families, test-brms-suite-helper-copy, test-brms-suite-methods, test-compat-preflight, test-conditional-effects-draws, test-conditions-census, test-conditions, test-cox-draws, test-default-priors-bounds, test-default-priors-brms, test-default-priors-dpar, test-default-priors-fit, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-formrobust-draws, test-generic-collision, test-gev-tmp-xi, test-gp-by-draws, test-gr-by-draws, test-incl-thres-draws, test-laplace-draws, test-linpred-stack, test-loo-predict, test-loo, test-message-uniqueness, test-mo-simplex-draws, test-nb2dm-draws, test-newdata-aterms-draws, test-ordinal-disc-thres-draws, test-ordinal-draws-names, test-ordinal-mixture-draws, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-prior-route, test-prior-update, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-scale, test-simplex-draws, test-simplex-family-notes, test-simplex-mvbf-draws, test-simulators, test-sratio-draws, test-stan-control, test-subset-rate-draws, test-thres-draws, test-tmbstan-build-guard, test-weights-unbounded, test-xbeta-zibb-hurdle-draws |
| prior (row order) | prior | `Intercept` | `x` | test-default-priors-brms, test-reparam |
| prior (row order) | prior | `Intercept` | `t2` | test-sampling-ported |
| prior (row order) | prior | `mix` | `` | test-laplace-draws |
| prior (row order) | prior | `mix` | `z` | test-loo |
| prior (row order) | prior | `t2` | `t3` | test-sampling-ported |
| prior (row order) | prior | `t3` | `t4` | test-sampling-ported |
| prior (row order) | prior | `t4` | `` | test-sampling-ported |
| prior (row order) | prior | `w` | `Intercept` | test-brms-pins |
| prior (row order) | prior | `x` | `` | test-add-criterion, test-adefects, test-bracket-access, test-brms-output, test-brms-pins, test-brms-shapes-draws, test-brms-suite-brmsformula, test-brms-suite-brmsterms, test-brms-suite-families, test-brms-suite-helper-copy, test-brms-suite-methods, test-compat-preflight, test-conditional-effects-draws, test-conditions-census, test-conditions, test-cox-draws, test-default-priors-bounds, test-default-priors-brms, test-default-priors-dpar, test-default-priors-fit, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-formrobust-draws, test-generic-collision, test-gev-tmp-xi, test-gp-by-draws, test-gr-by-draws, test-incl-thres-draws, test-laplace-draws, test-linpred-stack, test-loo-predict, test-loo, test-message-uniqueness, test-mo-simplex-draws, test-nb2dm-draws, test-newdata-aterms-draws, test-ordinal-disc-thres-draws, test-ordinal-draws-names, test-ordinal-mixture-draws, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-prior-route, test-prior-update, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-scale, test-simplex-draws, test-simplex-family-notes, test-simplex-mvbf-draws, test-simulators, test-sratio-draws, test-stan-control, test-subset-rate-draws, test-thres-draws, test-tmbstan-build-guard, test-weights-unbounded, test-xbeta-zibb-hurdle-draws |
| prior (row order) | prior | `x` | `Intercept` | test-add-criterion, test-adefects, test-bracket-access, test-brms-output, test-brms-pins, test-brms-shapes-draws, test-brms-suite-brmsformula, test-brms-suite-brmsterms, test-brms-suite-families, test-brms-suite-helper-copy, test-brms-suite-methods, test-compat-preflight, test-conditional-effects-draws, test-conditions-census, test-conditions, test-cox-draws, test-default-priors-bounds, test-default-priors-brms, test-default-priors-dpar, test-default-priors-fit, test-draws-methods, test-draws-no-stanfit, test-draws-sd-names, test-draws-spellings, test-evidence-ratio, test-formrobust-draws, test-generic-collision, test-gev-tmp-xi, test-gp-by-draws, test-gr-by-draws, test-incl-thres-draws, test-laplace-draws, test-linpred-stack, test-loo-predict, test-loo, test-message-uniqueness, test-mo-simplex-draws, test-nb2dm-draws, test-newdata-aterms-draws, test-ordinal-disc-thres-draws, test-ordinal-draws-names, test-ordinal-mixture-draws, test-parallel-chains, test-postfit-draws, test-ppcheck-loo, test-predfix-new-levels, test-prior-route, test-prior-update, test-re-formula-draws, test-reparam, test-sample-direct, test-sample-seed, test-sampling-ported, test-scale, test-simplex-draws, test-simplex-family-notes, test-simplex-mvbf-draws, test-simulators, test-sratio-draws, test-stan-control, test-subset-rate-draws, test-thres-draws, test-tmbstan-build-guard, test-weights-unbounded, test-xbeta-zibb-hurdle-draws |
| prior (row order) | prior | `x` | `w` | test-brms-pins |
| prior (row order) | prior | `x` | `3` | test-ordinal-mixture-draws |
| prior (row order) | prior | `x` | `2` | test-ordinal-mixture-draws |
| prior (row order) | prior | `z` | `x` | test-brms-pins, test-ordinal-mixture-draws |
| prior (row order) | prior | `z` | `Intercept` | test-brms-shapes-draws |
| prior (row order) | prior | `z` | `` | test-loo |

<!-- END suitediff-sample.md -->

## The ported brms suite

`dev/brmsnames-1024-port.sh` runs the 15 generated files in record
mode (`dev/brmsport-run.R` with a third argument), one process per
file, on each arm, into an ignored directory (`dev/brmsport-record.sh`
writes into the tracked `dev/brmsport-log/`). The lane arm's library
is `wt-brmsnames-lib`, which `dev/brmsport-run.R` puts first; that it
holds this build was checked separately (`find.package("frmtmb")` and
`exists("unsub_specials")` from that library). No row flips, and no
row changes its message. Few lines of the generated files call
`exp()`, `poly()`, `I()`, `log()` or `scale()` at all (the count per
file, `grep -c`, is below). `dev/brmsport-verdicts.tsv` is not edited.

```text
<!-- BEGIN port-grep.txt -->

tests/testthat/test-brms-suite-brm.R:0
tests/testthat/test-brms-suite-brmsfit-helpers.R:0
tests/testthat/test-brms-suite-brmsformula.R:0
tests/testthat/test-brms-suite-brmsterms.R:0
tests/testthat/test-brms-suite-data-helpers.R:0
tests/testthat/test-brms-suite-emmeans.R:0
tests/testthat/test-brms-suite-families.R:0
tests/testthat/test-brms-suite-methods.R:3
tests/testthat/test-brms-suite-priors.R:0
tests/testthat/test-brms-suite-standata.R:2
extensions/frmtmb.sample/tests/testthat/test-brms-suite-brmsformula.R:0
extensions/frmtmb.sample/tests/testthat/test-brms-suite-brmsterms.R:0
extensions/frmtmb.sample/tests/testthat/test-brms-suite-families.R:0
extensions/frmtmb.sample/tests/testthat/test-brms-suite-helper-copy.R:0
extensions/frmtmb.sample/tests/testthat/test-brms-suite-methods.R:2
total lines: 7

<!-- END port-grep.txt -->
```

<!-- BEGIN portdiff.md -->


Library lines, base:  

Library lines, lane:  

Assertions recorded: base 598 , lane 598 

Held: base 485 , lane 485 

On one arm only: 0 

Rows whose outcome flips: 0 


Rows with the same outcome and another message: 0 


<!-- END portdiff.md -->

## Every suite, one file per process

`dev/release/run-par.sh` with `FRMTMB_RELLIB` set to the lane library,
over every test file of core and the eight extensions, in two tiers:
plain (`NOT_CRAN=true`) and gated (also `FRMTMB_BRMS_FIT_TESTS`,
`FRMTMB_DRMTMB_FIT_TESTS` and `FRMTMB_FUZZ`, as `run-gated.ps1` sets
them). The gated tier's only skips are the scale tier, which has its
own runner (`run-scale.ps1`). These runs are after fix round 1. Totals
by `dev/brmsnames-1024-parsum.R`.

```text
<!-- BEGIN brmsnames-par-summary.txt -->

== tier plain : BRMSNAMES-PLAIN ran 443 of 443 
             pkg files  pass fail err skip warn
          frmtmb   275 19178    0   0  201    0
 frmtmb.coupling    11   551    0   0    5    0
      frmtmb.eam    32  1838    0   0    5    0
   frmtmb.latent    12   679    0   0    2    0
    frmtmb.learn    15   430    0   0   13    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    66  3033    0   0    9    0
  frmtmb.spatial     7   328    0   0    2    0
 frmtmb.survival    14   405    0   0    0    0
files with fail, err or warn: 0 
core library of each file:
C:/Users/adf44/source/r/wt-brmsnames-lib 
                                     443 
files with skips: 45 

== tier gated : BRMSNAMES-GATED ran 443 of 443 
             pkg files  pass fail err skip warn
          frmtmb   275 21740    0   0    1    0
 frmtmb.coupling    11   551    0   0    5    0
      frmtmb.eam    32  1844    0   0    3    0
   frmtmb.latent    12   679    0   0    2    0
    frmtmb.learn    15   501    0   0    2    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    66  3303    0   0    1    0
  frmtmb.spatial     7   391    0   0    0    0
 frmtmb.survival    14   405    0   0    0    0
files with fail, err or warn: 0 
core library of each file:
C:/Users/adf44/source/r/wt-brmsnames-lib 
                                     443 
files with skips: 7 


<!-- END brmsnames-par-summary.txt -->
```

## Every suite in one session, three orders

`dev/release/run-sessions.sh` with `FRMTMB_RELLIB` set to the lane
library: a package's whole suite in one R session, forward, reverse
and with every extension preloaded, in the plain tier (`splain`) and
the gated tier with the drmTMB and fuzz gates on (`sgated`). Before
fix round 1 all nine packages ran this way, 27 sessions per tier, all
fail 0, err 0, warn 0, every LEAK line tagged. After fix round 1, which
the block below shows, core and frmtmb.sample ran again, the two
packages the fixes reach (every extension's suite ran one file per
process above). Summary by `dev/brmsnames-1024-sesssum.sh`.

```text
<!-- BEGIN brmsnames-sessions-summary.txt -->

== splain
SESSION frmtmb forward files=265/275 pass=19178 fail=0 err=0 skip=201 warn=0
SESSION frmtmb reverse files=265/275 pass=19278 fail=0 err=0 skip=201 warn=0
SESSION frmtmb preload files=265/275 pass=19372 fail=0 err=0 skip=203 warn=0
SESSION frmtmb.sample forward files=61/66 pass=3033 fail=0 err=0 skip=9 warn=0
SESSION frmtmb.sample reverse files=61/66 pass=3033 fail=0 err=0 skip=9 warn=0
SESSION frmtmb.sample preload files=61/66 pass=3033 fail=0 err=0 skip=9 warn=0
BRMSNAMES-SPLAIN ran 6 of 6
LEAK lines: 30, untagged: 0
DETAIL FAIL/ERROR lines: 0
lib lines of the sessions:
      6 lib: frmtmb 0.74.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.coupling 0.6.1 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.eam 0.13.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.latent 0.7.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.learn 0.7.3 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.ode 0.7.2 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.sample 0.21.1 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.spatial 0.1.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.survival 0.1.0 C:/Users/adf44/source/r/wt-brmsnames-lib

== sgated
SESSION frmtmb forward files=275/275 pass=21740 fail=0 err=0 skip=1 warn=0
SESSION frmtmb reverse files=275/275 pass=21840 fail=0 err=0 skip=1 warn=0
SESSION frmtmb preload files=275/275 pass=21934 fail=0 err=0 skip=3 warn=0
SESSION frmtmb.sample forward files=66/66 pass=3303 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample reverse files=66/66 pass=3303 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample preload files=66/66 pass=3303 fail=0 err=0 skip=1 warn=0
BRMSNAMES-SGATED ran 6 of 6
LEAK lines: 135, untagged: 0
DETAIL FAIL/ERROR lines: 0
lib lines of the sessions:
      6 lib: frmtmb 0.74.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.coupling 0.6.1 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.eam 0.13.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.latent 0.7.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.learn 0.7.3 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.ode 0.7.2 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.sample 0.21.1 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.spatial 0.1.0 C:/Users/adf44/source/r/wt-brmsnames-lib
      6 lib: frmtmb.survival 0.1.0 C:/Users/adf44/source/r/wt-brmsnames-lib


<!-- END brmsnames-sessions-summary.txt -->
```

## R CMD check --as-cran

`dev/brmsnames-1024-check.sh frmtmb` and `frmtmb.sample`, each built
and checked in `dev/brmsnames-1024-check/<pkg>/`, with
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH, and
`R_LIBS` the lane library, then rellib-r12, then the user library.

```text
<!-- BEGIN check-summary.txt -->

== frmtmb
* checking HTML version of manual ... [26s] NOTE
Skipping checking math rendering: package 'V8' unavailable
tests: [ FAIL 0 | WARN 0 | SKIP 228 | PASS 19090 ]
Status: 1 NOTE

== frmtmb.sample
tests: [ FAIL 0 | WARN 0 | SKIP 17 | PASS 3019 ]
Status: OK


<!-- END check-summary.txt -->
```

## Fix round 1

The review (`dev/reviews/2026-10-24-brmsnames.md`) gave FIX. The code
of each finding is listed under "What changed"; the evidence is here.
Every count of the suites, sessions and checks above is from the runs
after these fixes.

**F1, drop1()/step() and re_formula.** `dev/brmsnames-1024-bugfix.R`,
the review's constructions (seeds 3 and 77). On 0.74.0 the `exp(x)`
row of `drop1()` has `Df 0` and the full model's AIC, and `step()`
keeps `z`; on the lane the row has `Df 1` and the AIC of `y ~ z`
refitted, and `step()` drops `z`. On 0.74.0 the `re_formula` naming the
fit's own term is refused; on the lane it gives the full prediction,
which is also that of the same model on a column `xx = exp(x)`. Both
are pinned in `test-brms-names-specials.R` ("drop1() and step() see an
exp() term", "re_formula can name a group-level exp() term"), and both
fail on rellib-r12 (`seen-fail.txt` above).

```text
<!-- BEGIN bugfix-base.txt -->

arm base C:/Users/adf44/source/r/rellib-r12 
Single term deletions

Model: y ~ exp(x) + z

            Df    AIC
<none>         607.84
.frm_exp(x)  0 607.84
z            1 605.84
AIC of y ~ z refitted: 625.68897 
AIC of y ~ exp(x) refitted: 605.84179 
step() keeps: Intercept .frm_expx z 

fitted(), re_formula = NULL:             1.226339906 1.745194277 1.149327590 
fitted(), re_formula = ~ (1 + exp(x) | g): fitted(): the re_formula term (1 + exp(x) | g) matches no group-level term of this fit. A term is kept when its grouping factor is one the fit has and its columns are among that term's columns. The fit's group-level terms are: (1 + .frm_exp(x) | g). brms drops such a term silently, which turns a misspelled grouping factor into a population-level prediction; use re_formula = NA to ask for that 
the same model on a column xx = exp(x):  1.226339906 1.745194277 1.149327590 

<!-- END bugfix-base.txt -->
```

```text
<!-- BEGIN bugfix-lane.txt -->

arm lane C:/Users/adf44/source/r/wt-brmsnames-lib 
Single term deletions

Model: y ~ exp(x) + z

       Df    AIC
<none>    607.84
exp(x)  1 625.69
z       1 605.84
AIC of y ~ z refitted: 625.68897 
AIC of y ~ exp(x) refitted: 605.84179 
step() keeps: Intercept expx 

fitted(), re_formula = NULL:             1.226339906 1.745194277 1.149327590 
fitted(), re_formula = ~ (1 + exp(x) | g): 1.226339906 1.745194277 1.149327590 
the same model on a column xx = exp(x):  1.226339906 1.745194277 1.149327590 

<!-- END bugfix-lane.txt -->
```

**F2, prior table order.** `dev/brmsnames-1024-order.R` compares the
row order of `default_prior()` and of `validate_prior()` with brms's on
five models (two dpars, group-level terms in both, a nonlinear model),
on the rows both packages have. On 0.74.0 none agrees; on the lane all
ten do. One test did read the order: `test-simplex-mvbf.R` pinned the
`lncor` rows of two `logistic_normal()` responses as `Y`, `V`, the
formula's order. brms 2.23.0 lists `V` first (it sorts by response;
`dev/brmsnames-1024-lncor.R`), and the test now pins brms's order. The
suite-wide name diff below lists these moves as "prior (row order)".

```text
<!-- BEGIN lncor.txt -->

   class resp
1  lncor    V
10 lncor    Y

<!-- END lncor.txt -->
```

```text
<!-- BEGIN order-base.txt -->

arm base C:/Users/adf44/source/r/rellib-r12 
y ~ exp(x) + z + (1 + z | g)                                 default_prior same order: FALSE; validate_prior: FALSE
  brms  : b///// b/z//// cor///// cor//g/// Intercept///// sd///// sd//g/// sd/Intercept/g/// sd/z/g/// sigma///// 
  frmtmb: Intercept///// b///// b/z//// sigma///// sd///// sd//g/// sd/Intercept/g/// sd/z/g/// cor///// cor//g/// 
  brms v  : b///// b/z//// Intercept///// sd///// sd//g/// sd/Intercept/g/// sd/z/g/// sigma///// 
  frmtmb v: Intercept///// b///// b/z//// sigma///// sd///// sd//g/// sd/Intercept/g/// sd/z/g/// 
y ~ x + f, sigma ~ z                                         default_prior same order: FALSE; validate_prior: FALSE
  brms  : b///// b/fb//// b/fc//// b/x//// Intercept///// b////sigma/ b/z///sigma/ Intercept////sigma/ 
  frmtmb: Intercept///// b///// b/fb//// b/fc//// b/x//// Intercept////sigma/ b////sigma/ b/z///sigma/ 
  brms v  : b///// b/fb//// b/fc//// b/x//// b////sigma/ b/z///sigma/ Intercept///// Intercept////sigma/ 
  frmtmb v: Intercept///// b///// b/fb//// b/fc//// b/x//// Intercept////sigma/ b////sigma/ b/z///sigma/ 
y ~ x + (1 | g), sigma ~ z + (1 | g)                         default_prior same order: FALSE; validate_prior: FALSE
  brms  : b///// b/x//// Intercept///// sd///// sd//g/// sd/Intercept/g/// b////sigma/ b/z///sigma/ Intercept////sigma/ sd////sigma/ sd//g//sigma/ sd/Intercept/g//sigma/ 
  frmtmb: Intercept///// b///// b/x//// Intercept////sigma/ b////sigma/ b/z///sigma/ sd///// sd//g/// sd/Intercept/g/// sd////sigma/ sd//g//sigma/ sd/Intercept/g//sigma/ 
  brms v  : b///// b/x//// b////sigma/ b/z///sigma/ Intercept///// Intercept////sigma/ sd///// sd////sigma/ sd//g/// sd/Intercept/g/// sd//g//sigma/ sd/Intercept/g//sigma/ 
  frmtmb v: Intercept///// b///// b/x//// Intercept////sigma/ b////sigma/ b/z///sigma/ sd///// sd////sigma/ sd//g/// sd/Intercept/g/// sd//g//sigma/ sd/Intercept/g//sigma/ 
yc ~ log(xp) + (1 | g) + (0 + x | f)                         default_prior same order: FALSE; validate_prior: FALSE
  brms  : b///// Intercept///// sd///// sd//f/// sd/x/f/// sd//g/// sd/Intercept/g/// 
  frmtmb: Intercept///// b///// sd///// sd//f/// sd/x/f/// sd//g/// sd/Intercept/g/// 
  brms v  : b///// Intercept///// sd///// sd//f/// sd/x/f/// sd//g/// sd/Intercept/g/// 
  frmtmb v: Intercept///// b///// sd///// sd//f/// sd/x/f/// sd//g/// sd/Intercept/g/// 
y ~ a + b * x, a ~ exp(z) + (1 | g), b ~ 1 + z, nl = TRUE    default_prior same order: FALSE; validate_prior: FALSE
  brms  : sigma///// b/////a b/Intercept////a sd/////a sd//g///a sd/Intercept/g///a b/////b b/Intercept////b b/z////b 
  frmtmb: b/////a b/Intercept////a b/////b b/Intercept////b b/z////b sigma///// sd/////a sd//g///a sd/Intercept/g///a 
  brms v  : b/////a b/Intercept////a b/////b b/Intercept////b b/z////b sd/////a sd//g///a sd/Intercept/g///a sigma///// 
  frmtmb v: b/////a b/Intercept////a b/////b b/Intercept////b b/z////b sigma///// sd/////a sd//g///a sd/Intercept/g///a 

<!-- END order-base.txt -->
```

```text
<!-- BEGIN order-lane.txt -->

arm lane C:/Users/adf44/source/r/wt-brmsnames-lib 
y ~ exp(x) + z + (1 + z | g)                                 default_prior same order: TRUE; validate_prior: TRUE
y ~ x + f, sigma ~ z                                         default_prior same order: TRUE; validate_prior: TRUE
y ~ x + (1 | g), sigma ~ z + (1 | g)                         default_prior same order: TRUE; validate_prior: TRUE
yc ~ log(xp) + (1 | g) + (0 + x | f)                         default_prior same order: TRUE; validate_prior: TRUE
y ~ a + b * x, a ~ exp(z) + (1 | g), b ~ 1 + z, nl = TRUE    default_prior same order: TRUE; validate_prior: TRUE

<!-- END order-lane.txt -->
```

**F3, F4, the design column's spelling.** Pinned by "set_prior(coef =)
takes the design column and brms's name" (four columns, both
spellings, through `set_prior()`, `prior_string()`, `validate_prior()`
and `frm()`, each asserting which coefficient a `normal(3, 0.001)`
moved, plus the duplicate refusal of the two spellings together),
"class sd and cs() take the design column's spelling too" (`I(x^2)`
and `poly(x, 2)2` as sd coefs, asserting which sd a `normal(2, 0.001)`
moved; `cs(I(x^2))` in `validate_prior()`), and frmtmb.sample's
"frm_sample() takes a coef in either spelling" (the draws of `b_IzE2`
at 3, and of no other coefficient).

**F6, the multivariate nonlinear prefix.** Case 40 of the name table
(`mvbf(bf(yn ~ a * exp(b * x), a ~ 1 + exp(z) + (1 | g), b ~ 1,
sigma ~ z, nl = TRUE), bf(y ~ x + exp(z)))`) against brms's brmsfit:
every interface is "ok" on the lane, including `fixef()` order, which
needed the response-by-response order too. Pinned in "a multivariate
nonlinear parameter is named as brms names it" (stored brms order, and
`hypothesis()` on `yn_a_Intercept`, `yn_a_expz`, `yn_b_Intercept` and
`sd`'s `g__yn_a_Intercept`), and through `bns_cases`, which compares
fixef() with brms's `standata()` names at test time. The suite-wide
diff accounts for the renamed prefixes as `mv-nl prefix`.

**F5, fixef() of the draws.** frmtmb.sample's `fixef.frmtmb_draws()`
orders its rows by `brms_fixef_rows()` of the draws' fit; the draws
themselves keep their columns. Pinned in "fixef() of the draws is in
brms's order" against brms's stored `fixef()` of case 32.

**F8, the dead alias machinery.** Every reference to the alias, the
`prot` environment and the `shared` argument on 0.74.0 and on the lane
(`git grep -n '\.frm_\|\bprot\b\|shared_env\|, shared)\|shared = NULL'
18b536fc -- R extensions/*/R`, and `grep -rn` of the same on the
worktree, the `.frm_obs`, `.frm_mo`, `.frm_mod` and `.frm_ps_` names of
other code left out). After `unsub_specials()` nothing evaluates the
alias, so the bindings had no reader; the two `assign()`s, `prot`, the
`shared` argument and the `.frm_cs` branch are removed, and every suite
passes without them (below). `shared_env` stays as the response's
`formula_env`, an empty child of the formula environment as before, so
a fit carries the same environment.

```text
<!-- BEGIN alias-grep.txt -->

== 18b536fc (0.74.0): every reference to the alias, the prot env or the shared arg
18b536fc:R/frame.R:695:      # after frmtmb), parse_linpred() protects the call as .frm_cs()
18b536fc:R/frame.R:697:      if (!calls_function(lhs, "cs") && !calls_function(lhs, ".frm_cs")) {
18b536fc:R/frame.R:700:      re$bar <- str2lang(gsub(".frm_cs(", "cs(", deparse1(re$bar),
18b536fc:R/parse.R:1733:parse_linpred <- function(rhs_form, env, shared = NULL) {
18b536fc:R/parse.R:1989:  # .frm_<name> aliases bound in a child environment.
18b536fc:R/parse.R:1990:  prot <- new.env(parent = env)
18b536fc:R/parse.R:1998:        pn <- paste0(".frm_", nm)
18b536fc:R/parse.R:2001:          assign(pn, fn, envir = prot)
18b536fc:R/parse.R:2019:  env_lp <- if (protected) prot else env
18b536fc:R/parse.R:2475:  shared_env <- new.env(parent = env)
18b536fc:R/parse.R:2741:    parse_linpred(reformulas::RHSForm(f, as.form = TRUE), env, shared_env)
18b536fc:R/parse.R:2786:                        environment(pf) %||% env, shared_env)
18b536fc:R/parse.R:2874:    formula_env = shared_env

== lane
R/parse.R:1714:#' Undo `parse_linpred()`'s `.frm_<name>` aliases of the special names
R/parse.R:1724:    if (startsWith(nm, ".frm_") && substring(nm, 6L) %in% specials) {
R/parse.R:2009:  # .frm_<name> for splitForm alone: unsub_specials() undoes it right
R/parse.R:2019:          e[[1]] <- as.name(paste0(".frm_", as.character(hd)))
R/parse.R:2049:  # model.matrix(), whose column `.frm_exp(x)` became the coefficient
R/parse.R:2050:  # `.frm_expx` where brms has `expx`.
R/parse.R:2498:  # evaluates. It held the parser's .frm_<name> aliases until those were
R/parse.R:2501:  shared_env <- new.env(parent = env)
R/parse.R:2900:    formula_env = shared_env

<!-- END alias-grep.txt -->
```

## Short check (SHIP), items R1 and R2

**R1, a smooth column's two spellings.** `prior_spec_slot()` now keys a
coef by `prior_coef_brms()` instead of `brms_rename()`, so
`s(x).fx1` and `sx_1` are one slot (`brms_rename()` gave `sx.fx1`).
`dev/brmsnames-1024-smoothdup.R` (seed 31): on 0.74.0 both priors are
accepted and the first is dropped without a word; on the lane the pair
is refused as a duplicate. Pinned in "a smooth column's two spellings
are one prior slot", which fails on rellib-r12 (2 failures). NEWS has a
bug-fix bullet.

```text
<!-- BEGIN smoothdup-base.txt -->

arm base C:/Users/adf44/source/r/rellib-r12 
validate_prior accepts both; user rows:
             prior class coef
 normal(-3, 0.001)     b sx_1
frm accepts both; fixef sx_1 = -3 

<!-- END smoothdup-base.txt -->
```

```text
<!-- BEGIN smoothdup-lane.txt -->

arm lane C:/Users/adf44/source/r/wt-brmsnames-lib 
validate_prior: Duplicated prior specifications are not allowed: 'b_sx_1' is given 2 times. Write one specification per slot, with its density and both bounds in the same call, e.g. set_prior("normal(0, 1)", class = "b", lb = 0) 
frm: Duplicated prior specifications are not allowed: 'b_sx_1' is given 2 times. Write one specification per slot, with its density and both bounds in the same call, e.g. set_prior("normal(0, 1)", class = "b", lb = 0) 

<!-- END smoothdup-lane.txt -->
```

**R2, a second nonlinear response that calls `exp()`.** The review's
case 15 (`dev/brmsnames-1024-mvnl2.R`, the review's data, seed 2410)
stops on 0.74.0 with 'could not find function ".frm_exp"' and fits on
the lane. Pinned in "a second nonlinear response may call exp()" on
this lane's data, which errors the same way on rellib-r12. NEWS has a
bug-fix bullet.

```text
<!-- BEGIN mvnl2-base.txt -->

arm base C:/Users/adf44/source/r/rellib-r12 
frm: could not find function ".frm_exp" 

<!-- END mvnl2-base.txt -->
```

```text
<!-- BEGIN mvnl2-lane.txt -->

arm lane C:/Users/adf44/source/r/wt-brmsnames-lib 
fits; logLik -287.5422161 
fixef: sigma_yn_Intercept sigma_y_Intercept yn_a_Intercept yn_a_z yn_b_Intercept sigma_yn_z y_c_Intercept y_c_expz y_d_Intercept sigma_y_x 

<!-- END mvnl2-lane.txt -->
```

The test files rerun after these two items, one file per process,
plain and gated (`NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`), on the lane
library:

```text
<!-- BEGIN short-tests.txt -->

== shortplain
RESULT test-brms-names-specials.R pass=257 fail=0 err=0 skip=0 warn=0
RESULT test-setprior.R pass=27 fail=0 err=0 skip=0 warn=0
RESULT test-brms-formula-priors.R pass=210 fail=0 err=0 skip=0 warn=0
RESULT test-brms-priors.R pass=0 fail=0 err=0 skip=12 warn=0
RESULT test-brms-suite-priors.R pass=0 fail=0 err=0 skip=0 warn=0
RESULT test-dirichlet-prior.R pass=24 fail=0 err=0 skip=0 warn=0
RESULT test-get-prior-route.R pass=80 fail=0 err=0 skip=0 warn=0
RESULT test-prior-compat.R pass=199 fail=0 err=0 skip=0 warn=0
RESULT test-prior-location-dpar.R pass=55 fail=0 err=0 skip=0 warn=0
RESULT test-prior-mv-resp.R pass=36 fail=0 err=0 skip=0 warn=0
RESULT test-prior-sd-coef.R pass=37 fail=0 err=0 skip=0 warn=0
RESULT test-prior-sd-scope.R pass=25 fail=0 err=0 skip=0 warn=0
RESULT test-prior-support.R pass=150 fail=0 err=0 skip=0 warn=0
RESULT test-priors-autocor-classes.R pass=62 fail=0 err=0 skip=0 warn=0
RESULT test-priors-bounds-grcov.R pass=50 fail=0 err=0 skip=0 warn=0
RESULT test-update-prior-drop.R pass=29 fail=0 err=0 skip=0 warn=0
RESULT test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
RESULT test-brms-parity-defects.R pass=172 fail=0 err=0 skip=0 warn=0
RESULT test-draws-names-specials.R pass=9 fail=0 err=0 skip=0 warn=0
RESULT test-default-priors-brms.R pass=40 fail=0 err=0 skip=0 warn=0
BRMSNAMES-SHORTPLAIN ran 20 of 20
lib lines: 20 of 20 name the lane library
== shortgated
RESULT test-brms-names-specials.R pass=257 fail=0 err=0 skip=0 warn=0
RESULT test-setprior.R pass=27 fail=0 err=0 skip=0 warn=0
RESULT test-brms-formula-priors.R pass=210 fail=0 err=0 skip=0 warn=0
RESULT test-brms-priors.R pass=103 fail=0 err=0 skip=0 warn=0
RESULT test-brms-suite-priors.R pass=36 fail=0 err=0 skip=0 warn=0
RESULT test-dirichlet-prior.R pass=24 fail=0 err=0 skip=0 warn=0
RESULT test-get-prior-route.R pass=80 fail=0 err=0 skip=0 warn=0
RESULT test-prior-compat.R pass=199 fail=0 err=0 skip=0 warn=0
RESULT test-prior-location-dpar.R pass=55 fail=0 err=0 skip=0 warn=0
RESULT test-prior-mv-resp.R pass=36 fail=0 err=0 skip=0 warn=0
RESULT test-prior-sd-coef.R pass=37 fail=0 err=0 skip=0 warn=0
RESULT test-prior-sd-scope.R pass=25 fail=0 err=0 skip=0 warn=0
RESULT test-prior-support.R pass=150 fail=0 err=0 skip=0 warn=0
RESULT test-priors-autocor-classes.R pass=62 fail=0 err=0 skip=0 warn=0
RESULT test-priors-bounds-grcov.R pass=50 fail=0 err=0 skip=0 warn=0
RESULT test-update-prior-drop.R pass=29 fail=0 err=0 skip=0 warn=0
RESULT test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
RESULT test-brms-parity-defects.R pass=172 fail=0 err=0 skip=0 warn=0
RESULT test-draws-names-specials.R pass=9 fail=0 err=0 skip=0 warn=0
RESULT test-default-priors-brms.R pass=40 fail=0 err=0 skip=0 warn=0
BRMSNAMES-SHORTGATED ran 20 of 20
lib lines: 20 of 20 name the lane library

<!-- END short-tests.txt -->
```

## What I did not do

- **OpenBLAS.** The change touches names, a formula's call heads and
  the order of `fixef()` rows; no arithmetic. The lane-env rule asks
  for the OpenBLAS pass on a change to numerical code or an exactness
  test, and this is neither, so it was not run.
- **The variable order of `variables()` and the draws** (see above):
  a positional contract of frmtmb.sample, recorded by the earlier lane.
- **conditional_effects() default effects, `s(exp(x))`, the mo()
  coding, the smooth `ranef()`, the `coef()` order**: filed below.
- **`frm_compat()` rows.** The brief lists "compat rows". No compat
  rule or feature names a coefficient, and nothing in `R/compat.R`
  reads a design column, so there was nothing to change; the core and
  extension suites, which test the compat tables, pass.

## Defects found and not fixed

1. **BACKLOG, next round (MAJOR): `mo()` of an integer that does not
   start at 0 has one simplex cell too many, per unit of the minimum.**
   brms 2.23.0 subtracts the minimum (`get_mo_values()`, brms
   `R/formula-sp.R` 603 to 627), so values 1 to 4 have `Jmo = 3`;
   frmtmb codes `D = max(x)` (`R/frame.R`, `codes <- as.integer(v)`,
   `D_mo <- max(codes)`), since 0.15.0 (65fed080, 2026-08-31, the
   commit that added `mo()`), and refuses negative integers, which brms
   takes. This lane's probe (seed 1024, `probe1.txt`, `probe2.txt`
   below): `mo(k)` with 1 to 4 has the log likelihood of `mo(k - 1)`
   with one more parameter, and no standard error for `(Intercept)`,
   `mok` or `zeta1_3`; the draws hold `simo_mok1[1..4]` where brms has
   `[1..3]`. The review's numbers (`dev/review-brmsnames-mo.R`, seed 5,
   n = 300, its finding 7), all with brms's `Jmo = 3` and `Xmo_1` in
   0 to 3:

   | term | frmtmb logLik (df) | Intercept (SE) | b (SE) |
   |---|---|---|---|
   | `mo(k0)`, 0 to 3 (brms's model) | -427.181869 (5) | 2.0834 (0.1138) | 1.5159 (0.0547) |
   | `mo(k)`, 1 to 4 | -427.181869 (6) | 1.8520 (NaN) | 1.1948 (NaN) |
   | `mo(k2)`, 2 to 5 | -427.181869 (7) | 1.9612 (NaN) | 0.9340 (NaN) |

   Point predictions agree (fitted() at k = 1 to 4 to 1e-6); the ML
   coefficients are not identified (NaN SEs, with a warning naming the
   cause); 3 of 4 predictions lose their SEs; and `frm_sample()`
   samples another model than brms with only ESS warnings: Intercept
   0.912 (sd 1.330) against brms's 2.098 (0.117), `b_mok` 1.426
   (0.336) against 1.510 (0.053), and a cell `simo_mok1[1]` with
   posterior mean 0.176 on a category no row is in (2 chains, 1000
   iterations, seed 11). Integer scales that start at 1 (Likert items)
   are the common input for `mo()`. The fix: subtract the minimum,
   store it for newdata as brms does (`attr(x, "min")`), and take
   negative integers.

```text
<!-- BEGIN probe1.txt -->

[1] 3

 1  2  3  4 
48 52 63 37 
[1] 3 2 1 2 0 1

<!-- END probe1.txt -->
```

```text
<!-- BEGIN probe2.txt -->

Warning message:
Standard errors are not available for 3 of 6 parameters. (Intercept), mok, zeta1_3 (simplex of mok): the likelihood is flat along them, or too nearly flat for its curvature to be measured, at the estimates, so the data do not determine them there (a standard deviation at 0 does this, and so do a mo() simplex split by a category no row is in and parameters that enter only through a combination). The other standard errors are kept, and a prediction that moves along these directions gets none; summary() and diagnose() list them. See the 'Convergence problems' section of vignette('diagnostics') 
'log Lik.' -287.1117 (df=6)
'log Lik.' -287.1117 (df=5)
           Estimate Est.Error Q2.5 Q97.5
Intercept 1.9136116       NaN  NaN   NaN
mok       0.1366851       NaN  NaN   NaN
           Estimate  Est.Error       Q2.5     Q97.5
Intercept 1.9544419 0.14675647 1.66680447 2.2420793
mok0      0.1686364 0.07414547 0.02331391 0.3139588
List of 1
 $ :List of 10
  ..$ expr     : symbol k
  ..$ codes    : int [1:200] 4 3 2 3 1 2 2 2 1 3 ...
  ..$ D        : int 4
  ..$ levels   : NULL
  ..$ zeta     : chr "zeta1"
  ..$ col      : int 2
  ..$ label    : chr "mok"
  ..$ mult     : NULL
  ..$ mult_expr: NULL
  ..$ mm_col   : chr NA
[1] "b_Intercept" "bsp_mok"     "sigma"      

<!-- END probe2.txt -->
```

2. **`s(exp(x))` is refused** ("The model uses `exp(x)`, which is not a
   column of `data`"); brms fits it (`sexpx_1`, `sds_sexpx_1`). The
   smooth's frame variables are taken as names (`dpar_frame_rhs()`,
   `as.name(tm)` over the term strings).
3. **conditional_effects() picks other default effects than brms** on
   `exp(x):z`, `t2(x, z)`, `mo(k):z`, `cs(x)`, `sigma ~ exp(z)`, `mi()`
   and multivariate models (Table 2; brms's lists are in `dump.txt`).
4. **`ranef()` of a smooth-only model** returns the smooth's penalized
   coefficients under the internal names `s(x).1`, ...; brms stops,
   and names those coefficients `s_sx_1[1]`, ... in `variables()`.
5. **`variables()` of an ML fit with `mo()`** has no `simo_` entries;
   brms's has `simo_mok1[1]`, ..., and frmtmb.sample's draws have them.
6. **The column order of `variables()` and of the draws** is
   frmtmb's template order, not brms's (every centered intercept
   first, then response by response); recorded by the earlier lane
   (`dev/brmsnames-findings.md`, "Found and NOT fixed", item 1). The
   summaries now follow brms: `fixef()` of the fit and of the draws.

## Versions and floors

- frmtmb: a minor bump (BREAKING name and order changes).
- frmtmb.sample: one code change (F5, `fixef()` of the draws), two new
  tests and NEWS bullets: a patch bump. Its new tests fail on frmtmb
  0.74.0 (seen above), and F5 calls `brms_fixef_rows()`, which it
  already imported. The floor is left to consolidation (review
  finding 9): raise it to the frmtmb release that carries this change,
  or skip the names tests on an older core.
- No other extension changes; every extension suite passes on this
  core (above).
