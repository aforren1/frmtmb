# Lane dblbar: every `||` term as brms reads it

Round opened 2026-10-13, base d2ef6312 (frmtmb 0.70.0), branch
`wt-dblbar`, library `C:/Users/adf44/source/r/wt-dblbar-lib`. Backlog
item r2 of "Filed by lane moparity after punch round 1"
(`dev/test-backlog.md`). Nothing is committed.

## What brms does (read from its source, then run)

brms 2.23.0 `split_re_terms()` (`R/formula-re.R`, lines 290 to 370 of
the tarball in `dev/brms-suite/`) reads `(lhs || g)` as
`(lhs | gr(g, cor = FALSE))`: "||-syntax overwrites the 'cor'
argument", then the mid section becomes `|`. The left-hand side stays
ONE term, coded by `get_model_matrix()` with the formula's contrasts,
and `frame_re()` takes its coefficients from the column names. An
`|ID|` key is only applied when `cor` is TRUE. Every claim below was
then RUN on brms 2.23.0 (`standata()`, `default_prior()`, `stancode()`,
no compile), not recalled.

## What changed

- `expand_double_verts()` (`R/parse.R`) makes every `||` term ONE
  `diag(lhs | rhs)` term. lme4's piecewise expansion is gone, except
  under `gr(cov = )` / `gr(prec = )` (below). Moparity's m7 special case
  for `mo()` is now the general case.
- `gr(g, cor = )` and `mm(g1, g2, cor = )` are read (`bar_read_cor()`,
  `drop_cor_arg()`): `cor = FALSE` is the `||` term, `cor = TRUE` is
  the plain term, and under `||` `cor` is ignored, as brms ignores it.
  `gr(g, cor = )` with nothing else is the factor `g`. `cor` must be
  TRUE or FALSE (`eval_spec_arg()`).
- `(x || p | g)` and `(x | p || g)` are refused in brms's words
  ("Invalid syntax used in group-level terms"). brms refuses both
  (`dev/dblbar-log/brms-probe.txt`, case idbar; `brms-id.txt`). The base
  read the first as a grouping variable `p` ("The model uses `p`, which
  is not a column of `data`") and crashed on the second.
- `(x | p | gr(g, cor = FALSE))` is refused. brms drops the key silently
  (`brms-id.txt`: two separate blocks, no `L_`), which the project
  refuses on principle (a silently ignored argument).
- `gr(cov = )` / `gr(prec = )` under `||`: the known-covariance block
  has an unstructured covariance within a level, and frmtmb has no
  diagonal one. So the term is still split into lme4's one-term pieces,
  each its own block, which is brms's model for numeric columns
  (`brms-id.txt`: `r_1_k = sd_1[k] * (Lcov_1 * z_1[k])` per column). A
  piece with several columns is refused in the frame (`within_struct`
  field, `R/frame.R`): `(1 + f || gr(g, cov = A))` and the explicit
  `diag(1 + x | gr(g, cov = A))`. Both silently became a CORRELATED
  `gr_cov` block on the base, because `gr(cov = )` overwrote the
  covariance class. Any other structure around `gr(cov = )` with
  several columns (`ar1(f + 0 | gr(g, cov = A))`,
  `homdiag(0 + x + z | gr(g, cov = A))`) is refused for the same
  reason; on the base it was silently `gr_cov` too. A ONE-column term
  in any structure is the plain known-covariance term and is fitted,
  as on the base: fix round 1 moved this check from the parser, which
  had refused `homdiag(0 + x | gr(g, cov = A))` (a correct `gr_cov/1`
  model on the base, review item 3), to the frame, where the column
  count is known.
- An uncorrelated block names no correlation: `block_has_cor()`
  (`R/covstruct.R`) is FALSE for the structures that are diagonal by
  construction (`diag`, `diag_t`, `homdiag`). `hyp_env_vals()`
  (`variables()`, hypothesis names on fits and on frmtmb.sample's
  draws) and `varcorr_trans_rows()` (`summary()$random`,
  `confint_varcorr()`) used `dim > 1`, so such a block listed `cor_`
  parameters fixed at 0: a `diag()` block, a `homdiag()` block (six for
  `homdiag(tim + 0 | g)` in the review's structs-cmp) and the `diag_t`
  block of a Student-t `||` term. `varcorr_layout()` already excluded
  `diag`/`homdiag` but not `diag_t`. Correction from the review (item
  2): on the base, frmtmb.sample's `variables(draws)` and
  `VarCorr(draws)` had no `cor_` for a diag block; only `hypothesis()`
  on the draws resolved one, at 0. The draws carry `theta_*`, not
  `sd_` names, on base and lane alike.
  brms has no cor parameter for a `cor = FALSE` term (`stancode()`: no
  `L_`, `Cor_` or `cor_` line in 31 of the 32 probe cases brms accepts;
  the 32nd, idmix, correlates only its `|p|` term;
  `dev/dblbar-corcount.R`). moparity's `(1 + mo(x) || g)`
  had `cor_g__Intercept__moxo` on the base.
- `re_term_cnms()` (`R/re-formula.R`) reads a `||` term of a
  `re_formula` with the `|` term's columns, as brms does. On the base
  `re_formula = ~ (1 + f | g)` on a `(1 + f || g)` fit dropped lme4's
  `fa` column and changed the prediction (test-double-bar.R, seen to
  fail).

Docs: `?frm` has a new section "Uncorrelated group-level terms, ||";
the mo paragraph, `?ranef`, `?frmtmb-multimembership`, three
compatibility rules (`R/compat.R`), the mm() refusal message and the
brms-migration vignette follow. NEWS has a development-version
section.

## brms comparison, case by case

Data: `dev/dblbar-data.R`, seed 20261013, n = 360, 24 levels of `g`, a
3-level factor `f`. Cases: `dev/dblbar-cases.R`. brms:
`dev/dblbar-brms-probe.R` (output `dev/dblbar-log/brms-probe.txt`).
frmtmb: `dev/dblbar-frm-probe.R base|lane` (frame, `default_prior()`,
a fit's `variables()`). Comparison: `dev/dblbar-compare.R`, table by
`dev/dblbar-table.R` (`dev/dblbar-log/table.txt`), pasted verbatim.

"Z" is the random-effect design column for column, up to the order of
the levels: each frmtmb block's per-coefficient n x levels matrices,
stacked, against brms's `Z_<id>_<k>` placed by `J_<id>` (for `mm()`,
`W_<id>_<m> * Z_<id>_<k>_<m>` summed over members), with identical
coefficient sets after brms's `rename()`. `mo()` columns have no `Z` in
brms's standata and are compared by name. "names" is the set of `sd_`
and `cor_` names (brms's rebuilt from `frame_re()`, by-levels by
`get_rnames()`'s rule). "sd prior rows" is the set of `default_prior()`
rows of class sd/cor with no coef (class, group, resp, dpar, nlpar).
"cor" is whether both correlate.

| case | spelling | brms sds | frm sds base/lane | Z base/lane | names base/lane | sd prior rows base/lane | cor base/lane |
|---|---|---|---|---|---|---|---|
| num1 | `bf(y ~ x + (1 + x \|\| g))` | 2 | 2/2 | yes/yes | yes/yes | yes/yes | yes/yes |
| num2 | `bf(y ~ x + z + (1 + x + z \|\| g))` | 3 | 3/3 | yes/yes | yes/yes | yes/yes | yes/yes |
| num0 | `bf(y ~ x + (0 + x \|\| g))` | 1 | 1/1 | yes/yes | yes/yes | yes/yes | yes/yes |
| fac1 | `bf(y ~ f + (1 + f \|\| g))` | 3 | 4/3 | NO/yes | NO/yes | yes/yes | NO/yes |
| fac0 | `bf(y ~ f + (0 + f \|\| g))` | 3 | 3/3 | yes/yes | NO/yes | yes/yes | NO/yes |
| facm1 | `bf(y ~ f + (-1 + f \|\| g))` | 3 | 3/3 | yes/yes | NO/yes | yes/yes | NO/yes |
| facb | `bf(y ~ f + (f \|\| g))` | 3 | 4/3 | NO/yes | NO/yes | yes/yes | NO/yes |
| mix | `bf(y ~ z + f + (1 + z + f \|\| g))` | 4 | 5/4 | NO/yes | NO/yes | yes/yes | NO/yes |
| inter | `bf(y ~ x * f + (1 + x * f \|\| g))` | 6 | 8/6 | NO/yes | NO/yes | yes/yes | NO/yes |
| inter2 | `bf(y ~ x + f:x + (1 + f:x \|\| g))` | 4 | 4/4 | yes/yes | NO/yes | yes/yes | NO/yes |
| ipoly | `bf(y ~ poly(x, 2) + (1 + poly(x, 2) \|\| g))` | 3 | 3/3 | yes/yes | NO/yes | yes/yes | NO/yes |
| iI | `bf(y ~ x + I(x^2) + (1 + I(x^2) \|\| g))` | 2 | 2/2 | yes/yes | yes/yes | yes/yes | yes/yes |
| grby | `bf(y ~ x + (1 + x \|\| gr(g, by = hg)))` | 4 | 4/4 | yes/yes | yes/yes | yes/yes | yes/yes |
| grbyf | `bf(y ~ f + (1 + f \|\| gr(g, by = hg)))` | 6 | 8/6 | NO/yes | NO/yes | yes/yes | NO/yes |
| grcor | `bf(y ~ f + (1 + f \|\| gr(g, cor = TRUE)))` | 3 | NA/3 | -/yes | -/yes | -/yes | -/yes |
| cross | `bf(y ~ x + (1 + x \|\| g1:g2))` | 2 | 2/2 | yes/yes | yes/yes | yes/yes | yes/yes |
| crossf | `bf(y ~ f + (1 + f \|\| g1:g2))` | 3 | 4/3 | NO/yes | NO/yes | yes/yes | NO/yes |
| nest | `bf(y ~ x + (1 + x \|\| g1/g2))` | 4 | 4/4 | yes/yes | yes/yes | NO/NO | yes/yes |
| nestf | `bf(y ~ f + (1 + f \|\| g1/g2))` | 6 | 8/6 | NO/yes | NO/yes | NO/NO | NO/yes |
| mm | `bf(y ~ x + (1 + x \|\| mm(m1, m2)))` | 2 | 2/2 | yes/yes | yes/yes | NO/NO | yes/yes |
| mmf | `bf(y ~ f + (1 + f \|\| mm(m1, m2)))` | 3 | 4/3 | NO/yes | NO/yes | NO/NO | NO/yes |
| idmix | `bf(y ~ x + f + (1 + x \| p \| g) + (1 + f \|\| g1))` | 5 | 6/5 | NO/yes | NO/yes | yes/yes | yes/yes |
| sigma | `bf(y ~ x + (1 \| g1), sigma ~ (1 + f \|\| g))` | 4 | 5/4 | NO/yes | NO/yes | yes/yes | NO/yes |
| sigman | `bf(y ~ x + (1 \| g1), sigma ~ (1 + x \|\| g))` | 3 | 3/3 | yes/yes | yes/yes | yes/yes | yes/yes |
| nl | `bf(y ~ a + b * x, a ~ 1 + (1 + f \|\| g), b ~ 1, nl = TRUE)` | 3 | 4/3 | NO/yes | NO/yes | yes/yes | NO/yes |
| mo1 | `bf(y ~ mo(xo) + (1 + mo(xo) \|\| g))` | 2 | 2/2 | yes/yes | NO/yes | yes/yes | NO/yes |
| twin | `bf(y ~ f + (1 + f \|\| g) + (1 + f \|\| g))` | 3 | 4/3 | NO/yes | NO/yes | yes/yes | NO/yes |
| cmc | `bf(y ~ 0 + f + (0 + f \|\| g), cmc = FALSE)` | 2 | 2/2 | yes/yes | NO/yes | yes/yes | NO/yes |
| mv | `mvbf(bf(y ~ f + (1 + f \|\| g)), bf(z ~ x + (1 + x \|\| g)), rescor = FALSE)` | 5 | 6/5 | NO/yes | NO/yes | yes/yes | NO/yes |

Generated totals (`table.txt`): lane, 29 cases compared, Z, names and
cor agree on 29, sd prior rows on 25; base, 28 compared, all three
agree on 9, sd prior rows on 24.

Cases where one side refuses:

| case | spelling | brms 2.23.0 | base | lane |
|---|---|---|---|---|
| idgr | `(1 + x \|\| gr(g, id = 'p'))` + `sigma ~ (1 \| gr(g, id = 'p'))` | one block of 3, uncorrelated | refused: gr() has no `id =` | same refusal (not this lane) |
| idbar | `(1 + x \|\| p \| g)` | "Invalid syntax used in group-level terms." | "The model uses `p`, which is not a column" | "Invalid syntax used in group-level terms: ..." |
| mof | `(mo(xo) * f \|\| g)` | errors, brms#1828 | moparity's brms 2.22.0 reading | unchanged |
| cs1, csx | `(cs(1) \|\| g)`, `(cs(1 + x) \|\| g)` (acat) | uncorrelated cs block | refused: no group-level cs() | same refusal (not this lane) |
| dup | `(1 + f \|\| g) + (1 \| g)` | "Duplicated group-level effects" | refused, same | refused, same; the message now names `diag(1 + f \| g)` |

`grcor` shows brms's override: `(1 + f || gr(g, cor = TRUE))` has no
correlation in brms (no `L_` in `stancode()`), and frmtmb now agrees;
the base refused `gr(cor = )`.

## Fits

Script `dev/dblbar-fits.R base|lane` (data seed 20261013),
comparison `dev/dblbar-fits-compare.R`, output
`dev/dblbar-log/fits-compare.txt`. brms fits: `dev/dblbar-brms-fit.R`
(4 chains of 2000, sampler seed 1, rstan; max R-hat 1.002, 1.001 and
1.002 for num1, fac1 and mix).

- Same model, three spellings, in one process: `(x || g)`,
  `(x | gr(g, cor = FALSE))` and `diag(x | g)` give IDENTICAL
  estimates, logLik and predictions (`identical()`, all of
  `fitted()` at new levels with `sample_new_levels = "gaussian"`, at
  `re_formula = NA` and at `re_formula = ~ (1 | g)`), for num1, fac1
  and mix. They are one parse, so this is an identity, checked.
- Numeric `(1 + x || g)` with GAUSSIAN group effects keeps its
  estimates (a Student-t one does not; see the post-fit section): lane
  against base,
  logLik -537.089366402417 against -537.089366402416 (difference
  -1.14e-13), shared values max relative difference 1.51e-12 (6.79e3
  ulps), fitted rows 2.21e-13, new-level predictions 1.95e-12. The
  same 1.51e-12 separates the lane's one diag block from the lane's
  explicit `(1 | g) + (0 + x | g)`: the likelihood is the same and the
  random effects come in another order, so the sparse factorization
  rounds differently. Not bitwise, by that construction.
- Factor `(1 + f || g)` against brms, `(ML - posterior median) /
  posterior sd`: b_Intercept 0.037, b_fb 0.022, b_fc -0.002,
  sd_g__Intercept -0.318, sd_g__fb -0.222, sd_g__fc -0.190, sigma
  -0.119. `(1 + z + f || g)`: b within 0.023; sd_g__Intercept -0.368,
  sd_g__fb -0.228, sd_g__fc -0.263, sigma -0.178, and sd_g__z -1.037,
  where the ML estimate is at the boundary (5.03e-05) and the posterior
  median of an sd cannot be. ML standard deviations sit below posterior
  medians, as expected of ML against a posterior with a
  student_t(3, 0, 2.5) prior. The base fit of the same data is another
  model: logLik -593.358 against the lane's -597.101 (fac1), -590.838
  against -594.540 (mix), with one more sd.
- Layout against the brms fit, fac1, mix and num1: `summary()$random`
  rows identical to brms's (`sd(Intercept) sd(fb) sd(fc)`), `VarCorr()$g`
  has only `sd`, as brms's has, with brms's row names; `ranef()` and
  `coef()` have brms's level and coefficient names (frmtmb's `ranef()`
  is a level-by-coefficient matrix, brms's an array with a summary
  dimension; that is frmtmb's general layout).
- New levels (`allow_new_levels = TRUE`, `sample_new_levels =
  "gaussian"`), fac1, rows f = b and f = c: frmtmb 1.097 (se 1.259) and
  0.509 (se 1.263) against brms posterior means 1.075 (sd 1.407) and
  0.517 (sd 1.412). frmtmb's se is the block's diagonal covariance
  (sd_Intercept^2 + sd_fb^2 = 1.229^2) plus the fixed-effect part.
  `"uncertainty"` is refused for a maximum likelihood fit, as before.
- frmtmb.sample on fac1 (`dev/dblbar-smoke2.R`, 2 chains of 600,
  seed 1): `VarCorr(draws)$g$sd` means 1.122, 0.753, 0.764 against
  brms's medians 1.083, 0.736, 0.740. The draws are named
  `r_g[g01,fb]` and `theta_*`; they carry no `sd_` and no `cor_` name,
  on the base too (review item 2), so they show nothing about this
  change. `hypothesis()` on them is where the zero `cor_` was.

## Post-fit paths (`dev/dblbar-smoke.R base|lane`)

`update(formula. = ~ . + x)` and `update(newdata = )` keep the one
diag block with brms's names; `default_prior(route = "sample")` lists
`student_t(3, 0, 2.5)` on class sd and group g and no cor; frmtmb.sample's
`frm_sample()`, `VarCorr()`, `posterior_samples()` and `as_tmbstan()`
run; no `cor_` anywhere. On the base the same calls list three
`cor_g__f*` names and a fourth sd.

`refit(fit, newresp = )` (`dev/dblbar-refit.R`, new response seed 5)
keeps `diag[(Intercept),fb,fc]` and brms's names; its logLik
-611.461283386 against -611.461283382 for a direct fit of the same
data.

**A Student-t `||` term is a different model now, numeric included
(BREAKING, review item 1).** `dev/dblbar-student.R base|lane`, data
seed 20261013, `dev/dblbar-log/student.txt`:

```
base : block diag_t dim 1 (Intercept) nu 5
base : block diag_t dim 3 fa,fb,fc nu 5
base : logLik -599.18149034
base : numeric block diag_t dim 1 (Intercept)
base : numeric block diag_t dim 1 x
base : numeric logLik -538.611555525
lane : block diag_t dim 3 (Intercept),fb,fc nu 5
lane : logLik -600.00577809
lane : numeric block diag_t dim 2 (Intercept),x
lane : numeric logLik -539.120960868
```

The lane has one block, so one latent scale per level shared by the
term's columns. brms 2.23.0 does the same: `udf_1` and `dfm_1`
multiply every `r_1_k` of the factor term (`brms-id.txt`) and of the
numeric term (the review's `brms-stud.txt`). The base gave each piece
its own scale, and named three zero `cor_` on the factor case. My
first record measured the factor case only and NEWS said the numeric
fit kept its estimates; NEWS now has a BREAKING bullet and the
1.5e-12 sentence is limited to gaussian group effects. A new test in
`test-double-bar.R` asserts one `diag_t` block for both cases.
frmtmb's `nu` stays fixed (5 here), its documented choice.

Every other reader of group-level terms is reached by the suites
below: `re_formula` (`R/re-formula.R`), `ranef()`/`coef()`,
`update()`, frmtmb.sample's naming (`hyp_env_vals()`,
`varcorr_layout()`), the compatibility table (read by frmtmb.eam's
`test-family.R`), and the four extensions that build random effects
(latent, coupling, learn, eam), whose R code generates no `||` term
(grep of `extensions/*/R`). No reader outside `R/parse.R` re-parses a
`||` term: `expandDoubleVerts`, `splitForm` and `findbars` are called
in `R/parse.R` only.

## Tests changed, and why

- `test-open-issues.R` "|| over a factor ...": pinned lme4's split,
  `(f || g)` as `(1 | g) + diag(0 + f | g)` with 2 blocks and numeric
  `||` as 2 blocks. Now `(f || g)` is `diag(1 + f | g)` (identical),
  one block with `(Intercept), fb, fc`, and a different likelihood from
  the lme4 split; numeric `||` is one block with the split's estimates
  (relative tolerance 1e-8) and logLik (1e-12).
- `test-mm.R`: `(1 + c1 || mm(g1, g2))` pinned 2 blocks; now 1, identical
  to `diag()` and to `mm(cor = FALSE)`. `mm(cor = FALSE)` was asserted
  refused; it is read now, so the refusal test uses `mm(id = )` and
  `cor` inside a structure.
- `test-brms-formula-priors.R`: `(x || g) + (x || g)` pinned 2 blocks;
  now 1.
- `test-brms-parity-defects.R` "old_levels reads ONE seen group":
  needed two blocks on one factor and used `(1 + x || g)` for it; now
  spells the two blocks `(1 | g) + (0 + x | g)`.
- `test-bcm-psychophysics.R`, `test-bcm-signal-detection.R`,
  `test-bcm-model-selection.R`: addressed lme4's pieces by their term
  labels (`"0 + xc | subj"`); now address the one block and its column.
- `helper-fuzz.R`: `dbar`/`dbar_f` block counts 2 and 3 become 1; the
  known-pending entry "double-bar-factor" (the brms_re_blocks
  invariant on `(fac || g)`) is removed because the coefficient count
  now matches brms; the "DIVERGENCE" comment is corrected.
- `test-mo-factor-re.R`: comment only.
- New `test-double-bar.R` (8 tests): the 15 spellings' sd names with no
  `cor_`; one block addressed by its factor; `gr(cor = )`, `mm(cor = )`
  and `diag()` identical to `||`; not lme4's model; the `|ID|`
  refusals with their absent cases; `gr(cov = )` numeric kept, factor
  refused; `re_formula`; and (structural brms tier, no compile) 10
  spellings' Z against `brms::standata()` and sd rows against
  `brms::default_prior()`. Seen to fail on rellib-r8: pass=52 fail=57
  err=5, every one of the 8 tests failing
  (`dev/dblbar-tests/dbl-base.log`); lane pass=137 fail=0.

- `test-methods.R` "an uncorrelated slope puts two blocks on one
  factor": pinned `(1 + x || g)` as two blocks; now asserts that the
  written-out `(1 | g) + (0 + x | g)` is two blocks and `(1 + x || g)`
  is one, with the same conditional modes (relative 1e-6).
- `test-parse.R` "double-bar and diag() terms parse ...": pinned 2
  parsed terms for `(x || g)`; now 1, of class diag.

## Suites

Driver `dev/dblbar-run-par.sh` (a copy of `dev/release/run-par.sh`
writing to `dev/dblbar-tests/`), one test file per R process, through
`dev/release/run-tests.R` with `FRMTMB_RELLIB` = the lane library ahead
of rellib-r8, so every extension loads the lane's core. Summaries by
`dev/dblbar-tier-summary.R`, pasted as generated.

First full run (`lane-all`, before two fixes): 372 files, 23375 pass,
6 fail, 3 error. All nine were mine. `test-methods.R` (4 fail, 1 err)
and `test-parse.R` (1 fail) pinned the two-block split and were
rewritten (above). `test-review-v29.R` (1 fail, 1 err) and
frmtmb.sample's `test-brms-pins.R` (1 err, "rr npar needs the rank")
found a fault in my first `block_has_cor()`, which counted parameters:
`equalto()` (no parameter, a fixed correlation) lost its `cor_` name,
and `rr()`'s `npar` cannot be called without the rank. The predicate
now names the diagonal structures instead.

Final ungated run on the final build (`final-all`, NOT_CRAN=true,
`dev/dblbar-tests/final-all.summary`):

```
             pkg files  pass fail err skip warn
          frmtmb   230 16586    0   0  179    0
 frmtmb.coupling    11   542    0   0    5    0
      frmtmb.eam    29  1747    0   0    3    0
   frmtmb.latent    10   360    0   0    2    0
    frmtmb.learn    15   430    0   0   13    0
      frmtmb.ode    11   549    0   0    1    0
   frmtmb.sample    51  2587    0   0    4    0
   frmtmb.spline    15   590    0   0    1    0
all: 372 files; pass 23391 fail 0 err 0 skip 208 warn 0
not clean: none
no result lines: 0
per-file logs naming the lane core: 372 of 372
```

Gated tier (`lane-gated`, FRMTMB_BRMS_FIT_TESTS=true, the 43 files that
call a gated skip; final R code, before the last test-only edit of
`test-methods.R`, which is not gated):

```
           pkg files pass fail err skip warn
        frmtmb    36 3090    0   0   13    0
 frmtmb.sample     7  524    0   0    0    0
all: 43 files; pass 3614 fail 0 err 0 skip 13 warn 0
```

The 13 skips are `test-drmtmb-agreement.R`'s own gate; run with
FRMTMB_DRMTMB_FIT_TESTS=true it is pass=131 fail=0 err=0 skip=0
(`lane-drm`).

Ported brms suite, bin 1: the 15 `test-brms-suite-*.R` files give the
same counts on rellib-r8 and on the lane (`base-suite.log` against
`lane-gated.log`: brm 23, brmsfit-helpers 3, brmsformula 16, brmsterms
5, data-helpers 6, emmeans 11, families 84, methods 162, priors 36,
standata 87; frmtmb.sample brmsformula 16, brmsterms 4, families 84,
helper-copy 36, methods 61; all pass, none stale). No bin-1 row flips:
`dev/brmsport-verdicts.tsv` has no row with a `||` term. brms's own
`||` assertions (`tests.standata.R` 480 to 491, `tests.stancode.R` 942,
1414 to 1419, 1589 to 1592) test brms's Stan data names and Stan code,
which are not in bin 1.

OpenBLAS 0.3.26 with LAPACK (`bash dev/ciharden-openblas.sh 0.3.26
lapack`; R copy in `dev/ciharden-out/`, deleted after), lane library,
one process per file (`dev/dblbar-tests/openblas/`):
test-double-bar.R pass=137, test-open-issues.R 46, test-mm.R 118,
test-mo-factor-re.R 34 (skip 2, gated), test-glmm-gaussian.R 16,
test-brms-formula-priors.R 208; fail=0 err=0 warn=0 in every one. The
`identical()` assertions compare spellings that are one parse, so they
hold under any BLAS; the cross-structure ones use relative tolerances
(1e-8 for estimates, 1e-12 for logLik).

The last two test-only edits (the title of the test-open-issues.R test,
and test-methods.R's comparison) were rerun on the final build
(`lane-last`): test-open-issues.R pass=46, test-methods.R pass=65,
fail=0 err=0 in both.

## R CMD check

frmtmb only (no extension changed): `R CMD build` then
`R CMD check --as-cran` in `dev/dblbar-check/`, with
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH, and
`R_LIBS` = the lane library, rellib-r8, the user library
(`dev/dblbar-check/check.out`). Tests [20m] OK, vignettes [319s] OK.
`Status: 1 NOTE`, the expected one: "Skipping checking math rendering:
package 'V8' unavailable". The tarball predates the retitling of one
test in test-open-issues.R (above).

## Fuzz tier (run by the reviewer)

The lane edited `helper-fuzz.R` (block counts, the removed
known-pending entry "double-bar-factor") but did not run the fuzz
tier; the appendix's `skip=1` for test-fuzz.R is its gate. The
reviewer ran it with `FRMTMB_FUZZ=true` and the full plan of 300
(`dev/dblbar-rev-fuzz.R`, `dev/dblbar-rev-tests/fuzz-lane.log` and
`fuzz-base.log`): lane pass=2 fail=0; base with the lane's helper
pass=1 fail=1, the `real_new` expectation, because without the
pending entry the base's `(fac || g)` coefficient count is a new
finding. The edit holds and is load-bearing.

## Fix round 1 (review `dev/reviews/2026-10-13-dblbar.md`)

1. NEWS: a BREAKING bullet for the Student-t `||` term (one latent
   scale per level, shared by its columns), numeric case measured
   above; the 1.5e-12 sentence limited to gaussian group effects.
2. NEWS: the `cor_` bullet corrected (`hypothesis()` on the draws, not
   their variables) and `homdiag()` named.
3. `homdiag(0 + x | gr(g, cov = A))` fitted again: the structure check
   moved from the parser to the frame and fires only on several
   columns (`within_struct`). Tested in `test-double-bar.R`; seen to
   fail on the pre-fix lane build (`dev/dblbar-tests/fix1-pre.log`:
   err=2, this and item 6).
4. Fuzz recorded above.
5. Docs: `?frm`'s `||` section says `mm(cor = FALSE)` is the `||` term
   over the membership grouping, and names the shared Student-t scale;
   the doubled "and" in `?frmtmb-multimembership` is gone; NEWS
   sentences are 25 words or fewer (`dev/dblbar-wordcount.R`).
6. The duplicate-effect message names the term as written (above);
   the `term_label` messages are filed.
7. The two pre-existing items are filed above; `?frmtmb-student-re` is
   corrected.

Reruns on the fixed build, one file per process, lane core in every
log. `fix1` (`dev/dblbar-tests/fix1.summary`): the 30 core files that
use `||`, `gr(cov/prec)`, `homdiag`, Student-t blocks or the duplicate
message, plus all 51 of frmtmb.sample's:

```
           pkg files pass fail err skip warn
        frmtmb    30 3670    0   0    2    0
 frmtmb.sample    51 2587    0   0    4    0
all: 81 files; pass 6257 fail 0 err 0 skip 6 warn 0
not clean: none
```

test-double-bar.R is pass=150 there. `fix1-gated`
(FRMTMB_BRMS_FIT_TESTS=true, frmtmb.sample's 7 gated files and
test-mo-factor-re.R): 8 files, pass 568, fail 0, err 0, skip 0, warn 0.

The Rd changed (`frm.Rd`, `frmtmb-student-re.Rd`,
`frmtmb-multimembership.Rd`), so `R CMD check --as-cran` was run again
on a fresh tarball, with the same environment as before
(`dev/dblbar-check/check.out`): tests [11m] OK, vignettes [240s] OK,
`Status: 1 NOTE` (V8 math rendering). The full 372-file tier was not
rerun. The fix touches the `gr(cov = )` branch of the parser, one
frame check, the duplicate message and docs, and every file that
reaches those paths is in `fix1`.

## Fix round 2 (review section "Round 2")

1. **Mixed spellings in the duplicate message.** The written-text map
   was keyed by the rewritten bar, so `(1 + f || g)` beside
   `(1 + f | gr(g, cor = FALSE))` or beside `diag(1 + f | g)` named
   one term twice. `expand_double_verts()` now keeps one queue per
   bar, in the order the terms are written, and records an explicit
   `diag()` term as well (NA, so it keeps its `diag(...)` name).
   `parse_linpred()` takes the next entry for each diag bar. Both
   cases are in `test-double-bar.R`; seen to fail on the round-1
   build (`dev/dblbar-tests/fix2-pre.log`: err=1, "appears in both
   (1 + f | gr(g, cor = FALSE)) and (1 + f | gr(g, cor = FALSE))").
2. Nits: the `double_bar` rule text in `R/compat.R` says
   `(x | mm(g1, g2, cor = FALSE))` is the `||` term over the
   membership grouping. A refused `||` term over `gr(cov = )` is named
   as written ("(1 + f || gr(g, cov = A)) has a term with 3 columns"),
   through the same queue, which now also takes each lme4 piece. A
   refused structure keeps its arguments (`within_label`):
   "rr(0 + f | gr(g, cov = A), d = 1) has 3 columns". Both are tested.
3. Filed below: frmtmb.sample's `seed =`.

Reruns on the round-2 build (`fix2.summary`, `fix2b.log`), one file
per process, lane core in every log: test-double-bar.R pass=154,
test-parse.R 44, test-arg-refusal.R 114, test-compat-register.R 103,
test-compat.R 735, test-brms-formula-priors.R 208,
test-mo-factor-re.R 34 (skip 2, gated), test-gr-by.R 71,
test-covstruct.R 7, test-mm.R 118, frmtmb.ode test-compat.R 15,
frmtmb.sample test-compat-preflight.R 86, frmtmb.eam test-family.R 61
(it reads the compatibility table), test-ce-levels.R 101,
test-predict-re-uncertainty.R 78, test-review-v29.R 138,
test-open-issues.R 46; fail=0 err=0 warn=0 in all. No Rd changed in
this round (`man/` diff as in round 1), so R CMD check was not rerun.

## Release

A minor bump: a breaking change in what `||` fits, two new readings
(`gr(cor = )`, `mm(cor = )`) and new refusals. No extension needs a new
floor on frmtmb: none calls a changed function directly, and their
suites pass on the lane's core as installed.

## Not done, and defects found but not fixed

- **sd priors per coefficient.** brms's `default_prior()` lists an sd
  row per coefficient (`sd fb g`), and brms honors
  `set_prior(class = "sd", coef = "fb", group = "g")`. frmtmb lists the
  class and group rows only, for every term, and refuses coef on class
  sd ("class = "sd" addresses a whole random-effect BLOCK"), pointing at
  class "theta". Pre-existing and the same on the base for `|` terms
  (`dev/dblbar-log/smoke-base.txt`). A `||` term is a diag block whose
  sds are separate parameters, so it could be honored there; it is a
  change to `R/priors.R`'s sd addressing for every term and is left to
  a priors lane.
- **Nested and multi-membership prior labels.** frmtmb's
  `default_prior()` labels the groups of `(1 + x | g1/g2)` `g2:g1` and
  `g1`, and of `(1 + x | mm(m1, m2))` `mm(m1, m2)`; brms writes `g1:g2`
  and `mmm1m2`, and `set_prior(group = "g1:g2")` or
  `group = "mmm1m2"` is refused ("No random-effect SDs match"). Same
  for `|` terms on the base (`dev/dblbar-preexist.R`,
  `dev/dblbar-log/preexist.txt`). This is the 4 of 29 sd-prior-row
  mismatches above.
- **Order of nested groups.** `variables()` lists `sd_g1:g2__*` before
  `sd_g1__*`; brms lists `g1` first. Pre-existing for `|` terms.
- `gr(g, id = )` is refused (frmtmb spells it `|ID|`); group-level
  `cs()` is refused; `re_formula = ~ diag(x | g)` names no group-level
  term (brms has no `diag()`), so it reads as `NA`. All pre-existing.
- A diagonal within-level covariance for `gr(cov = )` / `gr(prec = )`
  would let a factor sit in a `||` term over a known matrix. Not built:
  it is a new structure through every reader of `gr_cov`.
- Messages that name a `||` term by its rewritten form. Fix round 1
  carries the written text of each `||` and `cor = FALSE` term
  (`written_text`, set from `expand_double_verts()`'s map), and the
  duplicate-effect message now names `(1 + f || g)` or
  `(1 + f | gr(g, cor = FALSE))` as written (tested; seen to fail on the
  pre-fix lane build, `fix1-pre.log`). Not done: every message and
  table that reads a block's `term_label` still says
  `1 + f | g [diag]`, for example the class "cor" prior refusal and
  `confint_varcorr()$block`. `term_label` is also the key that
  `ranef()[[...]]`, summary rows and prior targeting match on, so
  changing it is not cosmetic; it is left as filed.
- **brms shares one Student-t latent scale across all the Student-t
  terms of a grouping factor** (pre-existing, the review's item 6). In
  `(1 | gr(g, dist = "student")) + (0 + x | gr(g, dist = "student"))`
  brms reuses `dfm_1` for the second term; frmtmb gives each block its
  own scale. `?frmtmb-student-re` said "brms builds it the same way",
  true only within a block; fix round 1 corrects the sentence and
  says to write one term to share the scale. The model difference
  itself is not fixed.
- **frmtmb.sample's `frm_sample(seed = )` does not fix every chain**
  (pre-existing, against frmtmb.sample; review round 2).
  `frm_sample(seed = 3)` with 2 chains gives different draws on two
  calls in one process unless `set.seed()` comes first; rellib-r8 does
  the same (`dev/dblbar-rev2-seed2.R`, `dblbar-rev2-log/seed2-base.txt`
  and `seed2-lane.txt`: "same seed identical: FALSE, with set.seed
  before each: TRUE"). One chain on `(1 | g)` reproduces. `seed =`
  should fix every chain's inits. Not fixed here.
- **frmtmb.sample's draws carry `theta_*`, not brms's `sd_` names**
  (pre-existing). `VarCorr(draws)` and `hypothesis()` on the draws
  derive the sds; `variables(draws)` and `as_draws_df()` do not list
  them.

## Overlaps with other lanes

R/parse.R (`expand_double_verts()` and helpers, the `gr()` branch),
R/frame.R (one refusal in the component builder), R/confint.R
(`varcorr_trans_rows()`, `hyp_env_vals()`), R/covstruct.R (new
`block_has_cor()`), R/methods-fit.R (`varcorr_layout()`, ranef doc),
R/re-formula.R, R/compat.R (three rule texts), R/fit.R (doc), NEWS.md,
vignettes/brms-migration.Rmd.

## Appendix: every test file's result (generated by dev/dblbar-appendix.R)

### final-all

```
frmtmb test-adefects.R pass=85 fail=0 err=0 skip=0 warn=0
frmtmb test-aliased-grouping.R pass=33 fail=0 err=0 skip=0 warn=0
frmtmb test-allfit-start.R pass=14 fail=0 err=0 skip=0 warn=0
frmtmb test-api-spellings.R pass=36 fail=0 err=0 skip=0 warn=0
frmtmb test-arg-refusal.R pass=114 fail=0 err=0 skip=0 warn=0
frmtmb test-arma-na-newdata.R pass=14 fail=0 err=0 skip=0 warn=0
frmtmb test-aterm-expr.R pass=69 fail=0 err=0 skip=0 warn=0
frmtmb test-autocor-cond.R pass=57 fail=0 err=0 skip=0 warn=0
frmtmb test-autocor.R pass=192 fail=0 err=0 skip=0 warn=0
frmtmb test-autoscale.R pass=46 fail=0 err=0 skip=0 warn=0
frmtmb test-backlog.R pass=35 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-bart.R pass=8 fail=0 err=0 skip=1 warn=0
frmtmb test-bcm-binomial.R pass=17 fail=0 err=0 skip=6 warn=0
frmtmb test-bcm-data-analysis.R pass=20 fail=0 err=0 skip=6 warn=0
frmtmb test-bcm-esp.R pass=5 fail=0 err=0 skip=2 warn=0
frmtmb test-bcm-gaussian.R pass=11 fail=0 err=0 skip=3 warn=0
frmtmb test-bcm-gcm.R pass=9 fail=0 err=0 skip=2 warn=0
frmtmb test-bcm-latent-mixtures.R pass=30 fail=0 err=0 skip=5 warn=0
frmtmb test-bcm-model-selection.R pass=13 fail=0 err=0 skip=3 warn=0
frmtmb test-bcm-mpt.R pass=15 fail=0 err=0 skip=2 warn=0
frmtmb test-bcm-psychophysics.R pass=7 fail=0 err=0 skip=2 warn=0
frmtmb test-bcm-retention.R pass=15 fail=0 err=0 skip=3 warn=0
frmtmb test-bcm-signal-detection.R pass=6 fail=0 err=0 skip=3 warn=0
frmtmb test-bcm-simple.R pass=8 fail=0 err=0 skip=1 warn=0
frmtmb test-bernoulli-coding.R pass=47 fail=0 err=0 skip=0 warn=0
frmtmb test-boot.R pass=49 fail=0 err=0 skip=0 warn=0
frmtmb test-bracket-access.R pass=33 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-agreement.R pass=172 fail=0 err=0 skip=2 warn=0
frmtmb test-brms-api-formrobust.R pass=26 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-families.R pass=1106 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-formula-priors.R pass=208 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-likelihood.R pass=17 fail=0 err=0 skip=52 warn=0
frmtmb test-brms-methods.R pass=18 fail=0 err=0 skip=47 warn=0
frmtmb test-brms-names.R pass=138 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-parity-defects.R pass=172 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-port.R pass=21 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-postfit.R pass=76 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-priors.R pass=0 fail=0 err=0 skip=12 warn=0
frmtmb test-brms-shapes-punch.R pass=55 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-shapes-punch2.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-shapes.R pass=72 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brm.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brmsfit-helpers.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brmsformula.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brmsterms.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-data-helpers.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-emmeans.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-families.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-methods.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-priors.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-standata.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-utilities.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb test-car-spde.R pass=138 fail=0 err=0 skip=0 warn=0
frmtmb test-case-studies.R pass=31 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-aterms.R pass=19 fail=0 err=0 skip=1 warn=0
frmtmb test-ce-bands.R pass=169 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-facets.R pass=92 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-levels.R pass=101 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-options.R pass=32 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-parity.R pass=30 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-plot.R pass=73 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-predict-median.R pass=11 fail=0 err=0 skip=1 warn=0
frmtmb test-cens-lccdf.R pass=27 fail=0 err=0 skip=0 warn=0
frmtmb test-cens-trunc.R pass=69 fail=0 err=0 skip=0 warn=0
frmtmb test-ci-siblings.R pass=14 fail=0 err=0 skip=0 warn=0
frmtmb test-compat-register.R pass=103 fail=0 err=0 skip=0 warn=0
frmtmb test-compat.R pass=735 fail=0 err=0 skip=0 warn=0
frmtmb test-conditional-smooths.R pass=39 fail=0 err=0 skip=0 warn=0
frmtmb test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb test-conditions.R pass=150 fail=0 err=0 skip=0 warn=0
frmtmb test-confint-anova.R pass=35 fail=0 err=0 skip=0 warn=0
frmtmb test-core-boundary.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb test-covstruct-floor.R pass=940 fail=0 err=0 skip=0 warn=0
frmtmb test-covstruct.R pass=7 fail=0 err=0 skip=0 warn=0
frmtmb test-cox-chart.R pass=17 fail=0 err=0 skip=0 warn=0
frmtmb test-cs-bar-brms.R pass=9 fail=0 err=0 skip=0 warn=0
frmtmb test-cs-factor.R pass=84 fail=0 err=0 skip=0 warn=0
frmtmb test-cumulative-cs.R pass=52 fail=0 err=0 skip=0 warn=0
frmtmb test-custom-family.R pass=126 fail=0 err=0 skip=0 warn=0
frmtmb test-data2.R pass=31 fail=0 err=0 skip=0 warn=0
frmtmb test-dates.R pass=13 fail=0 err=0 skip=0 warn=0
frmtmb test-diagnostics-ux.R pass=130 fail=0 err=0 skip=0 warn=0
frmtmb test-dirichlet-prior.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb test-disc-no-fixed.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb test-distributional.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb test-double-bar.R pass=137 fail=0 err=0 skip=0 warn=0
frmtmb test-dpar-equate.R pass=36 fail=0 err=0 skip=0 warn=0
frmtmb test-dpar-eta-api.R pass=123 fail=0 err=0 skip=0 warn=0
frmtmb test-drmtmb-agreement.R pass=0 fail=0 err=0 skip=13 warn=0
frmtmb test-dry-run.R pass=23 fail=0 err=0 skip=0 warn=0
frmtmb test-edgecases.R pass=64 fail=0 err=0 skip=0 warn=0
frmtmb test-effects.R pass=54 fail=0 err=0 skip=0 warn=0
frmtmb test-emmeans.R pass=91 fail=0 err=0 skip=1 warn=0
frmtmb test-famgaps.R pass=100 fail=0 err=0 skip=0 warn=0
frmtmb test-families.R pass=252 fail=0 err=0 skip=0 warn=0
frmtmb test-family-list.R pass=7 fail=0 err=0 skip=0 warn=0
frmtmb test-fd-chain.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb test-fit-end-hook.R pass=18 fail=0 err=0 skip=0 warn=0
frmtmb test-formula-cmc.R pass=49 fail=0 err=0 skip=0 warn=0
frmtmb test-formula-dot.R pass=26 fail=0 err=0 skip=0 warn=0
frmtmb test-frame.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb test-functional.R pass=3 fail=0 err=0 skip=0 warn=0
frmtmb test-fuzz.R pass=0 fail=0 err=0 skip=1 warn=0
frmtmb test-generic-collision.R pass=73 fail=0 err=0 skip=0 warn=0
frmtmb test-get-prior-route.R pass=79 fail=0 err=0 skip=0 warn=0
frmtmb test-glmm-binomial.R pass=10 fail=0 err=0 skip=0 warn=0
frmtmb test-glmm-gaussian.R pass=16 fail=0 err=0 skip=0 warn=0
frmtmb test-glmm-poisson.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb test-gp-by.R pass=71 fail=0 err=0 skip=2 warn=0
frmtmb test-gp-multidim.R pass=67 fail=0 err=0 skip=0 warn=0
frmtmb test-gr-by.R pass=71 fail=0 err=0 skip=0 warn=0
frmtmb test-grad-verdict.R pass=74 fail=0 err=0 skip=0 warn=0
frmtmb test-habit-prep.R pass=16 fail=0 err=0 skip=0 warn=0
frmtmb test-hurdle-cum-thres-cs.R pass=50 fail=0 err=0 skip=0 warn=0
frmtmb test-hurdle-nb-zoib.R pass=85 fail=0 err=0 skip=0 warn=0
frmtmb test-id-kron.R pass=61 fail=0 err=0 skip=0 warn=0
frmtmb test-importance.R pass=226 fail=0 err=0 skip=0 warn=0
frmtmb test-influence-plot.R pass=21 fail=0 err=0 skip=0 warn=0
frmtmb test-input-validation.R pass=43 fail=0 err=0 skip=0 warn=0
frmtmb test-interop.R pass=45 fail=0 err=0 skip=0 warn=0
frmtmb test-lkj.R pass=212 fail=0 err=0 skip=0 warn=0
frmtmb test-map.R pass=7 fail=0 err=0 skip=0 warn=0
frmtmb test-me.R pass=70 fail=0 err=0 skip=0 warn=0
frmtmb test-mean-fn.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb test-method-residue.R pass=64 fail=0 err=0 skip=0 warn=0
frmtmb test-methods-audit.R pass=58 fail=0 err=0 skip=0 warn=0
frmtmb test-methods.R pass=65 fail=0 err=0 skip=0 warn=0
frmtmb test-mi-bsp-names.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb test-mixture-reference.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb test-mixture-start.R pass=5 fail=0 err=0 skip=0 warn=0
frmtmb test-mm.R pass=118 fail=0 err=0 skip=0 warn=0
frmtmb test-mo-factor-re.R pass=34 fail=0 err=0 skip=2 warn=0
frmtmb test-mo-terms.R pass=74 fail=0 err=0 skip=0 warn=0
frmtmb test-multinomial.R pass=17 fail=0 err=0 skip=0 warn=0
frmtmb test-multiple-methods.R pass=39 fail=0 err=0 skip=0 warn=0
frmtmb test-multiple-pooling.R pass=33 fail=0 err=0 skip=0 warn=0
frmtmb test-multivariate.R pass=23 fail=0 err=0 skip=0 warn=0
frmtmb test-mv-gaps.R pass=58 fail=0 err=0 skip=0 warn=0
frmtmb test-mvn-mixture.R pass=277 fail=0 err=0 skip=0 warn=0
frmtmb test-naming-collisions.R pass=29 fail=0 err=0 skip=0 warn=0
frmtmb test-nl-body-vars.R pass=23 fail=0 err=0 skip=0 warn=0
frmtmb test-nl-lexical.R pass=31 fail=0 err=0 skip=0 warn=0
frmtmb test-nl-rtmb-scope.R pass=105 fail=0 err=0 skip=0 warn=0
frmtmb test-nl.R pass=85 fail=0 err=0 skip=0 warn=0
frmtmb test-nlf.R pass=77 fail=0 err=0 skip=0 warn=0
frmtmb test-nonfinite-gradient.R pass=5 fail=0 err=0 skip=0 warn=0
frmtmb test-numerical-robustness.R pass=723 fail=0 err=0 skip=0 warn=0
frmtmb test-offset-grid.R pass=27 fail=0 err=0 skip=1 warn=0
frmtmb test-open-issues.R pass=46 fail=0 err=0 skip=0 warn=0
frmtmb test-optima.R pass=107 fail=0 err=0 skip=0 warn=0
frmtmb test-optrap.R pass=145 fail=0 err=0 skip=0 warn=0
frmtmb test-ordinal-disc-thres.R pass=179 fail=0 err=0 skip=0 warn=0
frmtmb test-ordinal-fitted.R pass=129 fail=0 err=0 skip=0 warn=0
frmtmb test-ordinal-mixture-brms.R pass=0 fail=0 err=0 skip=3 warn=0
frmtmb test-ordinal-mixture.R pass=99 fail=0 err=0 skip=0 warn=0
frmtmb test-ordinal-thres-names.R pass=46 fail=0 err=0 skip=0 warn=0
frmtmb test-ordinal.R pass=99 fail=0 err=0 skip=0 warn=0
frmtmb test-osa-inference.R pass=34 fail=0 err=0 skip=0 warn=0
frmtmb test-par-template.R pass=51 fail=0 err=0 skip=0 warn=0
frmtmb test-parity-integration.R pass=8 fail=0 err=0 skip=0 warn=0
frmtmb test-parse.R pass=44 fail=0 err=0 skip=0 warn=0
frmtmb test-perf.R pass=10 fail=0 err=0 skip=0 warn=0
frmtmb test-pooled-anova.R pass=64 fail=0 err=0 skip=0 warn=0
frmtmb test-portability.R pass=96 fail=0 err=0 skip=0 warn=0
frmtmb test-pp-check-types.R pass=198 fail=0 err=0 skip=0 warn=0
frmtmb test-predfix.R pass=92 fail=0 err=0 skip=0 warn=0
frmtmb test-predict-lp-basis.R pass=30 fail=0 err=0 skip=0 warn=0
frmtmb test-predict-newdata.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb test-predict-re-uncertainty.R pass=78 fail=0 err=0 skip=0 warn=0
frmtmb test-prior-compat.R pass=199 fail=0 err=0 skip=0 warn=0
frmtmb test-prior-location-dpar.R pass=55 fail=0 err=0 skip=0 warn=0
frmtmb test-prior-mv-resp.R pass=36 fail=0 err=0 skip=0 warn=0
frmtmb test-prior-sd-scope.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb test-priors-autocor-classes.R pass=62 fail=0 err=0 skip=0 warn=0
frmtmb test-priors-bounds-grcov.R pass=50 fail=0 err=0 skip=0 warn=0
frmtmb test-ps.R pass=71 fail=0 err=0 skip=0 warn=0
frmtmb test-quadrature-defects.R pass=59 fail=0 err=0 skip=0 warn=0
frmtmb test-re-formula-partial.R pass=44 fail=0 err=0 skip=0 warn=0
frmtmb test-reml-nongaussian.R pass=9 fail=0 err=0 skip=0 warn=0
frmtmb test-reml-shapes.R pass=8 fail=0 err=0 skip=0 warn=0
frmtmb test-review-fixes.R pass=18 fail=0 err=0 skip=0 warn=0
frmtmb test-review-v25.R pass=60 fail=0 err=0 skip=0 warn=0
frmtmb test-review-v28.R pass=74 fail=0 err=0 skip=0 warn=0
frmtmb test-review-v29.R pass=138 fail=0 err=0 skip=0 warn=0
frmtmb test-rl-example.R pass=96 fail=0 err=0 skip=2 warn=0
frmtmb test-rsv-intercept.R pass=56 fail=0 err=0 skip=0 warn=0
frmtmb test-sandwich.R pass=58 fail=0 err=0 skip=0 warn=0
frmtmb test-scale-contract.R pass=37 fail=0 err=0 skip=0 warn=0
frmtmb test-se-check.R pass=127 fail=0 err=0 skip=0 warn=0
frmtmb test-se-flags.R pass=68 fail=0 err=0 skip=0 warn=0
frmtmb test-se-lost-re-predict.R pass=20 fail=0 err=0 skip=0 warn=0
frmtmb test-se-tier.R pass=91 fail=0 err=0 skip=0 warn=0
frmtmb test-setprior.R pass=27 fail=0 err=0 skip=0 warn=0
frmtmb test-simulate-density.R pass=493 fail=0 err=0 skip=0 warn=0
frmtmb test-simulate-ergonomics.R pass=52 fail=0 err=0 skip=0 warn=0
frmtmb test-simulate-newdata.R pass=57 fail=0 err=0 skip=0 warn=0
frmtmb test-skew-normal-start.R pass=70 fail=0 err=0 skip=0 warn=0
frmtmb test-smooth-population.R pass=40 fail=0 err=0 skip=0 warn=0
frmtmb test-smooths.R pass=60 fail=0 err=0 skip=0 warn=0
frmtmb test-sparsear1.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb test-sparsex.R pass=64 fail=0 err=0 skip=0 warn=0
frmtmb test-special-names.R pass=37 fail=0 err=0 skip=0 warn=0
frmtmb test-spectral.R pass=113 fail=0 err=0 skip=0 warn=0
frmtmb test-sratio-thresholds.R pass=29 fail=0 err=0 skip=0 warn=0
frmtmb test-structure.R pass=135 fail=0 err=0 skip=0 warn=0
frmtmb test-subset-rate.R pass=86 fail=0 err=0 skip=0 warn=0
frmtmb test-sugar.R pass=45 fail=0 err=0 skip=0 warn=0
frmtmb test-tabular-inputs.R pass=9 fail=0 err=0 skip=0 warn=0
frmtmb test-thres-refit.R pass=118 fail=0 err=0 skip=0 warn=0
frmtmb test-thres.R pass=67 fail=0 err=0 skip=0 warn=0
frmtmb test-tmb-examples.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb test-tre.R pass=102 fail=0 err=0 skip=0 warn=0
frmtmb test-trunc-postfit.R pass=35 fail=0 err=0 skip=0 warn=0
frmtmb test-unpinned-seams.R pass=30 fail=0 err=0 skip=0 warn=0
frmtmb test-update-pool.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb test-update-prior-drop.R pass=29 fail=0 err=0 skip=0 warn=0
frmtmb test-v07.R pass=25 fail=0 err=0 skip=0 warn=0
frmtmb test-v11.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb test-v12.R pass=26 fail=0 err=0 skip=0 warn=0
frmtmb test-v14.R pass=81 fail=0 err=0 skip=0 warn=0
frmtmb test-v15.R pass=41 fail=0 err=0 skip=0 warn=0
frmtmb test-v16.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb test-v17.R pass=37 fail=0 err=0 skip=0 warn=0
frmtmb test-v18.R pass=27 fail=0 err=0 skip=0 warn=0
frmtmb test-v19.R pass=23 fail=0 err=0 skip=0 warn=0
frmtmb test-verbose.R pass=40 fail=0 err=0 skip=0 warn=0
frmtmb test-vignette-wiener.R pass=21 fail=0 err=0 skip=0 warn=0
frmtmb test-xbeta-zibb-hurdle-cum.R pass=606 fail=0 err=0 skip=0 warn=0
frmtmb test-zi.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-coherence.R pass=63 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-conditions.R pass=5 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-cross-pairs.R pass=52 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-cross-spectrum.R pass=96 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-cross-wishart.R pass=211 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-epochs.R pass=41 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-message-uniqueness.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.coupling test-scale.R pass=0 fail=0 err=0 skip=5 warn=0
frmtmb.coupling test-surface.R pass=58 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-brms-parity.R pass=18 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-conditions.R pass=8 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-contaminant.R pass=32 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-defects.R pass=60 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-density.R pass=132 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-extension-api.R pass=14 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-family.R pass=61 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-gddm-conditions.R pass=58 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-gddm-family.R pass=77 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-gddm-gradients.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-gddm-recovery.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-gddm-reference.R pass=166 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-gddm-solver.R pass=94 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-lba.R pass=109 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-message-uniqueness.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-moments.R pass=29 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-ndt-bound.R pass=79 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-ndt-seam.R pass=83 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-prior-dpar.R pass=7 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-rdm-gng.R pass=221 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-sampling.R pass=100 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-scale.R pass=0 fail=0 err=0 skip=3 warn=0
frmtmb.eam test-simulate-density.R pass=69 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-surface.R pass=50 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-units.R pass=32 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-variability.R pass=140 fail=0 err=0 skip=0 warn=0
frmtmb.eam test-wiener-cdf.R pass=42 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-conditions.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-hmm-starts.R pass=84 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-hmm.R pass=119 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-lca.R pass=105 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-message-uniqueness.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-prior-dpar.R pass=10 fail=0 err=0 skip=0 warn=0
frmtmb.latent test-scale.R pass=0 fail=0 err=0 skip=2 warn=0
frmtmb.latent test-structure-latent.R pass=20 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-clean-session.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-conditions.R pass=8 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-counterfactual.R pass=69 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-engine.R pass=20 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-factorization.R pass=35 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-families.R pass=45 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-message-uniqueness.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-reference.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-rlddm-ndt.R pass=53 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-scale.R pass=0 fail=0 err=0 skip=2 warn=0
frmtmb.learn test-session.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb.learn test-stan-identity.R pass=0 fail=0 err=0 skip=11 warn=0
frmtmb.learn test-surface.R pass=146 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-compat.R pass=15 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-conditions.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-lincmt.R pass=65 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-ode-events.R pass=142 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-ode-nlf.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-ode-tv.R pass=27 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-ode.R pass=72 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-records.R pass=204 fail=0 err=0 skip=0 warn=0
frmtmb.ode test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
frmtmb.sample test-add-criterion.R pass=41 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-adefects.R pass=20 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-output.R pass=21 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-pins.R pass=334 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-shapes-draws.R pass=78 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-brmsformula.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-brmsterms.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-families.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-helper-copy.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-methods.R pass=0 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-compat-preflight.R pass=86 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-conditional-effects-draws.R pass=63 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-conditions.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-cox-draws.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-default-priors-brms.R pass=25 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-draws-methods.R pass=160 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-draws-no-stanfit.R pass=21 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-draws-spellings.R pass=97 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-evidence-ratio.R pass=36 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-formrobust-draws.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-generic-collision.R pass=70 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-gp-by-draws.R pass=20 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-gr-by-draws.R pass=7 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-incl-thres-draws.R pass=43 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-laplace-draws.R pass=92 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-loo.R pass=102 fail=0 err=0 skip=2 warn=0
frmtmb.sample test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-mo-simplex-draws.R pass=40 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-ordinal-disc-thres-draws.R pass=44 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-ordinal-draws-names.R pass=31 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-ordinal-mixture-draws.R pass=34 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-parallel-chains.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-postfit-draws.R pass=68 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-ppcheck-loo.R pass=27 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-predfix-new-levels.R pass=62 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-prior-route.R pass=9 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-prior-update.R pass=7 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-re-formula-draws.R pass=22 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-reparam.R pass=371 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-sample-direct.R pass=137 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-sampling-ported.R pass=206 fail=0 err=0 skip=1 warn=0
frmtmb.sample test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
frmtmb.sample test-simulators.R pass=54 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-sratio-draws.R pass=2 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-stan-control.R pass=47 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-subset-rate-draws.R pass=19 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-thres-draws.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-tmbstan-build-guard.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-xbeta-zibb-hurdle-draws.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-bracket-access.R pass=1 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-conditions-census.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-conditions.R pass=5 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-curve.R pass=82 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-deriv.R pass=48 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-difference.R pass=71 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-frailty.R pass=30 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-gratia.R pass=13 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-message-uniqueness.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-new-levels.R pass=57 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-royston-parmar.R pass=67 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-rp-floored.R pass=81 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-scale.R pass=0 fail=0 err=0 skip=1 warn=0
frmtmb.spline test-span.R pass=73 fail=0 err=0 skip=0 warn=0
frmtmb.spline test-surface.R pass=47 fail=0 err=0 skip=0 warn=0
```

### lane-gated

```
frmtmb test-bcm-bart.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-binomial.R pass=37 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-data-analysis.R pass=46 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-esp.R pass=17 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-gaussian.R pass=23 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-gcm.R pass=17 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-latent-mixtures.R pass=52 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-model-selection.R pass=29 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-mpt.R pass=27 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-psychophysics.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-retention.R pass=32 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-signal-detection.R pass=30 fail=0 err=0 skip=0 warn=0
frmtmb test-bcm-simple.R pass=13 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-agreement.R pass=191 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-likelihood.R pass=543 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-methods.R pass=1006 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-priors.R pass=103 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brm.R pass=23 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brmsfit-helpers.R pass=3 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brmsformula.R pass=16 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-brmsterms.R pass=5 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-data-helpers.R pass=6 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-emmeans.R pass=11 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-families.R pass=84 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-methods.R pass=162 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-priors.R pass=36 fail=0 err=0 skip=0 warn=0
frmtmb test-brms-suite-standata.R pass=87 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-aterms.R pass=37 fail=0 err=0 skip=0 warn=0
frmtmb test-ce-predict-median.R pass=14 fail=0 err=0 skip=0 warn=0
frmtmb test-drmtmb-agreement.R pass=0 fail=0 err=0 skip=13 warn=0
frmtmb test-emmeans.R pass=109 fail=0 err=0 skip=0 warn=0
frmtmb test-gp-by.R pass=86 fail=0 err=0 skip=0 warn=0
frmtmb test-mo-factor-re.R pass=44 fail=0 err=0 skip=0 warn=0
frmtmb test-offset-grid.R pass=35 fail=0 err=0 skip=0 warn=0
frmtmb test-ordinal-mixture-brms.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb test-rl-example.R pass=102 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-brmsformula.R pass=16 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-brmsterms.R pass=4 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-families.R pass=84 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-helper-copy.R pass=36 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-brms-suite-methods.R pass=61 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-loo.R pass=110 fail=0 err=0 skip=0 warn=0
frmtmb.sample test-sampling-ported.R pass=213 fail=0 err=0 skip=0 warn=0
```

### lane-drm

```
frmtmb test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0
```

### lane-last

```
frmtmb test-open-issues.R pass=46 fail=0 err=0 skip=0 warn=0
frmtmb test-methods.R pass=65 fail=0 err=0 skip=0 warn=0
```
