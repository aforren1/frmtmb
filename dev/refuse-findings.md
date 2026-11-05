# Lane refuse: models brms 2.23.0 builds that frmtmb refused

Round opened 2026-10-29, base d5d0056d (frmtmb 0.78.0, rellib-r16).
Nothing is committed. Library: `wt-refuse-lib`. Every table in section 9
is generated from the logs by `dev/refuse-gen-tables.R`. Review:
`dev/reviews/2026-11-02-refuse.md`. Section 8 says what fix round 1 changed
against what this record first said.

## 1. What changed

- **D2.** `class = "sds"` is a class on `frm()` and `frm_sample()`, brms's
  word and meaning. `coef` is the term as R deparses it, one row for all
  the penalties of a `t2()` and all the levels of a by-smooth. Class `"sd"`
  no longer reaches a smooth (BREAKING): `class = "sd", group = "s(x)"` and
  `class = "sd"` on a smooth-only model are refused, naming `sds`, as brms
  refuses them (`dev/refuse-brms-sds.R`). `sds` with `group`, a by-level
  coef and `s(z,k=5)` for `s(z, k = 5)` are refused as brms refuses them.
  On the sampling route the default is brms's `student_t(3, 0, scale)` on
  class `sds`, and `prior_plan()` reads the `sds` rows.
- **D5.** `(x | p | mm(g1, g2))`, `mm(g1, g2, id = "p")` and
  `mm(g1, g2, cov = A)`, with `mmc()` and several columns, and
  `(1 + x || mm(g1, g2, cov = A))`. Different weights under one key are
  refused in both orders, as are conflicting `cor` flags: brms reads the
  first term's. A non-symmetric or non-positive-definite `A` is refused by
  name, on `gr(cov =)` too.
- **D8.** `mi()` on the response takes brms's 14 families. On `frm()` a
  value no density reads leaves the objective (BREAKING, fix round 1: the
  coordinator's decision on B1); a value a density reads stays latent on
  `frmtmb_family(mi_chart =)`'s chart: `log` (exponential, Gamma, weibull,
  frechet, inverse.gaussian, lognormal), `logit` (beta), `circle` (von_mises), `identity`
  (gaussian, student, skew_normal, exgaussian, gen_extreme_value), `no_mode`
  (asym_laplace; a mixture takes its components' chart and refuses a read
  value). `frm_sample()` keeps every value; its draws are `Ymi_<resp>[k]` on
  the response's scale (BREAKING).
- **Cheap sweep fixes.** `se()` x `skew_normal`; `mi()` x `se()`,
  `mi(sd)` x `se()`; `gr(g, id = , cor = FALSE)` when every term under the key
  says so; `cens()` x `mixture`.
- `frm_compat()` rows, section 5.

One new export, `mi_sampling()` (section 8, M5), and one new argument,
`frmtmb_family(mi_chart =)`.

## 2. D2 evidence

`dev/refuse-table-compare.R` (seed 43): for `s(x)`, `s(x) + s(z, k = 5)`,
`s(z, bs = "cr", k = 5)`, `t2(x, z)`, `s(x, by = f)`, `s(x) + (1 + w | g)`,
`sigma ~ s(z)`, an nlpar smooth and a two-response smooth, frmtmb's rows
(class, coef, group, resp, dpar, nlpar, lb) equal brms's in order, except
`lb = 0` on the `sd` and `sigma` class rows: that is the item filed at
0.78.0 ("`default_prior()` gives `sigma` no `lb = 0`") and I left it.
`test-prior-sds.R` pins the `sd` and `sds` rows on nine forms against
brms's `default_prior()`.

`prior_summary()` lists the rows as written, on both routes, for every
class: brms also lists the flat and vectorized rows. That is the design
since 0.78.0 (`prior_rows_as_written()`), not specific to a smooth;
`validate_prior()` lists the vectorized `sds` row as brms's
`prior_summary()` does. Not changed.

L1, through the instrument: p043's rerun is section 6. Variants in
section 9: brms's log density and the objective differ by one offset, the
same at all three points on each form. The offset is not a property of the
model (m1): the review reads -0.290, +0.232, -0.470 on its own data where
this lane read -0.667, -0.848, -1.97; it is brms's normalizer of the half
density plus `frm()`'s log-Jacobian at the reference point. NEWS and the compat
row say "one offset" now, with no value.

- `sds_two` (`coef = "s(z, k = 5)"` alone) reads offsets 2.13, 2.18, 2.08.
  brms keeps its default on `s(x)`, which the row does not name, and frm()
  has no default there, so this is the two defaults differing, not the
  class. `sds_two2` names both and is constant at 1.38.
- `sds_sigma2`: the translator stops ("'data' must be of a vector type"),
  with no prior too (`flat_sigma_smooth`), so a smooth in `sigma` alone is
  a translator gap. `sigma ~ s(z)` is pinned by the table and by
  `sds_entries()` only.

## 3. D5 evidence

p070 reaches L1, offset -4.15e-09, gradient 1.61e-09, once the translator
reads `y.mu:mmc(x1, x2)` (`dev/equiv-lib.R`). p087 is L1 at 5.68e-14,
gradient 1.11e-15. `mm(cov = I)` equals plain `mm()` to a ratio below 1e-7,
and a non-identity `A` is another model (`test-mm-id-cov.R`). Tables equal
brms's rows but the same `lb` item (`dev/refuse-log/table-compare.log`).

The merged block needed no change to the phase-3 Z or to prediction: a
component's Z is built from its triplets, and `mm_newdata_parts()` reads
the component's offset. Checked on `fitted(newdata =)` against the fit's own
rows (ratio below 1e-10), and on `allow_new_levels = TRUE`.

## 4. D8 evidence

brms's set, by running it (`dev/refuse-brms-mi.R`, section 9): an mi()
response builds and compiles for gaussian, student, skew_normal, beta,
gamma, weibull, exponential, frechet, inverse.gaussian, lognormal,
exgaussian, gen_extreme_value, von_mises and asym_laplace, and every other
family stops at "Argument 'mi' is not supported". An mi() predictor takes
every family with a scalar response in both packages (section 9).

**A value no density reads** (`dev/refuse-f1-unread.R`, section 8, B1).
`logLik` and `BIC` equal the complete-case fit's to at most 6.5e-13 on all 14
families, and the coefficients to at most 4.2e-7. brms's log_prob on the
observed rows equals `logLik` to at most 3.2e-13 on 13 families;
gen_extreme_value has no rule for `tmp_xi`, and is L1 on `xi ~ 1` (`mi_gev`).
`asym_laplace` reads gradient 2.07, its kink at the optimum.

**A value a density reads.** Joint density against brms's log_prob, a gaussian
reader of the family's response, largest 2.23e-8 (student, `nu` estimated
large), the rest at most 5.7e-14; the table in section 9. An mi() predictor
under other families reading a gaussian variable (`dev/refuse-mi-pred.R`) is an
L1 identity under poisson, bernoulli, binomial, negbinomial, categorical,
Gamma, beta, student, hurdle_gamma, zero_inflated_poisson, von_mises; cumulative
and sratio fit, but the translator stops on an ordinal response in a
multivariate model. m5: asym_laplace's response reading a gaussian variable
fits, and holds as a value, but its inner gradient is 0.9 at the kink, and the
review reads the Laplace error up to +0.59 per value. Pre-existing, filed.

**Laplace against `integrate()`, per value** (part B of
`dev/refuse-f1-unread.R`, the review's design, 80 rows, 8 missing, seed 9300 +
family index): gaussian 0, Gamma -0.114, exponential -0.098, weibull -0.078,
frechet +0.008, beta -0.125, exgaussian -0.079, von_mises -0.050, lognormal
-0.037, skew_normal -0.023, inverse.gaussian +0.021, gen_extreme_value +0.008,
student 0 at this fit's `nu`. NEWS and `?bf` state these.

`mi_sd` works with a chart: L1 on `lognormal` with `mi(s)`, 3.84e-13
(`mi_lognormal_sd`).

## 5. Compat rows changed

Each on its own line in `R/compat.R`, each checked by running both packages
(`dev/refuse-sweep.R`; tests in `test-refuse-pairs.R`, `test-compat.R`):

| row | was | now | evidence |
|---|---|---|---|
| `cens()` x `disc` (p010) | untested | refused | both refuse; cens on an ordinal family |
| `cens()` x `equidistant` (p011) | untested | refused | both refuse |
| `vint()` x `sum_to_zero` (p059) | untested | refused | both refuse |
| `trials()` x `cs_pred()` (p146) | works | refused | both refuse, on binomial and on cumulative |
| `s()` x `prior` (p043) | untested | works | L1 |
| `t2()` x `prior` | untested | works | L1 |
| `\|ID\|` x `mmc()` (p070) | conditional | works | L1 |
| `gr_cov` x `mmc()` (p087) | conditional | works | L1 |
| `mm()` x `\|ID\|` | refused | works | L1; one set of weights |
| `mm()` x `gr_cov` | refused | works | L1 |
| `mi()` x `group:mi_families` (new group; p169) | refused | works | L1 per family, both routes |
| `mi()` x `asym_laplace` | refused | conditional | fits unread rows; a read value is refused |
| `mi()` x `kind:family` | refused | refused, brms's words | brms refuses |
| `mi()` x `cens()` | refused | conditional | with nothing latent both fit; brms's program is broken otherwise |
| `mi()` x `trunc()` | refused | conditional | with nothing latent both fit; filed otherwise |
| `mi()` x `se()` | refused | works | L1 |
| `mixture` x `mi()` | refused | conditional | fits unread rows; a read value is refused |
| `importance` x `mi()` | refused | conditional | unread rows have left the objective |
| `quadrature` x `mi()` | refused | conditional | unread rows have left the objective |
| `se()` x `skew_normal` | refused | works | L1 |
| `cens()` x `mixture` | refused | works | L1 |
| `trunc()` x `mixture` | refused | refused, new reason | brms truncates per component |
| `thres()` x `mvbf` | refused | works | both fit; the row was wrong |

`mi()` x `group:gaussian_like` is gone, replaced by `mi_families`. The
proposal for p014 (cs_pred x nl), p049, p166, p171 and the nl sweep rows are
lane nlfix's, p013 predread's, p069 and p124 impsamp's.

## 6. Equivalence rerun

`dev/refuse-f2-final.sh`, with absolute paths in `FRMTMB_EQUIV_LIB`:
dev/refuse-equiv-out, all 172 pairs on the final build, against the baseline
frmtmb-wt-r16base/dev/r16base-equiv-out (dev/refuse-verdict-diff.R,
dev/refuse-log/verdict-diff.txt):

```
p029 mi() x frm_curve: cannot measure (L4) -> equivalent (L1+L2); l1_const 1.14e-13 grad 3.93e-14
p043 s() x prior: defect (NA) -> defect (L1); l1_const -0.667 grad 1.62e-13
p070 |ID| x mmc(): defect (NA) -> equivalent (L1); l1_const -4.15e-09 grad 1.61e-09
p087 gr_cov x mmc(): defect (NA) -> equivalent (L1); l1_const 5.68e-14 grad 1.11e-15
p169 von_mises x mi(): defect (NA) -> defect (L1); l1_const 1.14e-13 grad 1.87e-05
 - **agree with brms (the headline): 84.6, Wald [79.3, 89.9]; bootstrap parametric [79.1, 
base: - **agree with brms (the headline): 83.6, Wald [78.2, 89.0]; bootstrap parametric [7
```

The other 167 keep verdict and level.

**Correction.** Fix round 1's report said it had rerun all 172 pairs on its final
build. It had not: `FRMTMB_EQUIV_LIB="wt-refuse-lib;rellib-r16"` holds relative
paths, so every pair stopped on "no package called 'frmtmb'", and the summary was
rebuilt from the pair files of the earlier run. Those were one message word from
that build. The same slip reached this round's first try. This table is the real
run: `p169.tsv` is dated 15:48, no log says FAILED, measured in 497 s.

- p043 and p169 stay `defect` only because
  `dev/equiv-data/adjudication.tsv` (rows 17 and 8) sets the verdict by
  hand. Their per-pair files read `equivalent` and `compat row wrong` at L1.
  **Drop both rows at consolidation.** The review ran that: 85.5.
- p169 reads `compat row wrong`, and p010, p011, p059, p146 stay there,
  because `sample.tsv` freezes the status the flow compares against; the
  rows themselves are fixed.
- p029 moved because the translator reads the univariate `Ymi`. It and p169
  are the complete-case route: brms on the observed rows.

## 7. The sweep

Every `refused` rule with a brms spelling (`dev/refuse-sweep.R`, stancode
plus a cached compile, table in section 9). What brms builds and frm()
still refuses after this lane:

| pair | decision |
|---|---|
| `mi()` x `cens()` with a latent value | NEWS: brms's program reads `Y`, Inf at a missing row; `log_prob` is -Inf at a finite point (`dev/refuse-l1-iddiag.R`); `dev/upstream-bugs.md` brms-45 |
| `trunc()` x `mixture` | NEWS: brms truncates per component (`dev/refuse-log/brms-code.log`) |
| `gr(by =, cov =)`, `mm(by =, cov =)` | NEWS: no block density |
| `gr(cov =, dist = "student")`, `mm(dist = "student")` | NEWS: not a multivariate t, no density to expand |
| a read value under `mixture`, `asym_laplace` | NEWS: no single mode with curvature |
| different weights under one `mm()` key | NEWS: brms reads the first term's; brms-46 |
| conflicting `cor` under one key | NEWS: brms reads the first term's; brms-47 |
| `mi()` x `trunc()` with a latent value | filed: brms bounds `Ymi` by the window; a per-row chart composed with the family's |
| `rescor` x `mi()`, x `se()`, x `weights()` | filed: brms's programs are simple (`Yl` in the MVN; `se` as sigma; weights on the row), but `predict_rescor_draw()`, `rescor_row_loglik()` and frmtmb.sample's draws read sigma and y directly, so not cheap |
| `\|ID\|` x `gr(dist = "student")` | filed: merged blocks are assembled gaussian; frmtmb's `nu` is fixed where brms estimates it, so L1 is not direct |
| `mm(pw =)` | filed: `target += PW[j] * std_normal_lpdf(z[1, j])`, a weight per level in the block density |
| `mixture` x `se()` | filed: the compat note has the reason, `se_dpar` names one `sigma`; brms has `sigma1 = se` |
| `mo(x):mi(z)` | filed: brms builds; the parser takes the special on one side only |

Methods: `simulate()` on an `mvbf()` is refused where brms's
`posterior_predict()` runs and frmtmb's `predict()` runs
(`dev/refuse-sweep2.R`).
Filed; the rows `mvbf` x `kind:method` and `rescor` x `simulate` are right.

No brms spelling, not swept: REML, quadrature, profile, importance, ps(),
`residuals_osa`. Left to other lanes, as the brief says: `group:autocor` x
`kind:family` (latres, D6), `group:autocor` x `nl` (nlfix, D7), `nl` x
`kind:family` (nlfix), `importance` rows (impsamp).

Where the sweep agrees (both refuse): `autocor` x weights, cens, trunc, se,
rescor, mixture, a second autocor term; `cens()`/`trunc()` x categorical,
ordinal, simplex; `se()` x lognormal; `rate()`, `disc`, `thres()`, `cs()` x a
family without them; `subset()` x rescor, me(); `rescor` x cens, trunc,
poisson; `mixture` x rescor; `bar_crossing`; `0 + Intercept` x ordinal.

### Reading the variants table

- `cens_mix` is the first try, before `lcdf` returned a probability: the L1
  read -3.01, varying by point, gradient 2.47, which is how the wrong
  contract was found. `cens_mix_none` is the fix, -5.68e-14, with
  `mix_control` beside it (no censoring), 0. Both are under `order = "none"`.
- `id_diag` stops in the translator, which reads an id group as one block:
  frmtmb holds independent blocks. `dev/refuse-l1-iddiag.R` writes the map:
  offset 0, inner gradient 2.08e-10, and the fit equals the keyless fit's
  logLik (-311.4481 both). The same script: brms's `Y` at a missing row is
  Inf under `mi() + cens()`, and log_prob is -Inf at a finite point.
- `skew_se` and `skew_se_sigma` read gradient 1.5e-5 and 1.2e-5, below the
  bound 1e-3; the offset is 2e-13.

## 8. Fix round 1

The review's verdict was FIX: `dev/reviews/2026-11-02-refuse.md`. What this
record first said about `mi()` was wrong in one way: it reported a table of
Laplace errors from one shape each and called them constants.

- **B1.** `mi_ml_frame()` in `R/frame.R`, run at the top of `fit_assembled()`,
  so that a deletion refit of `influence()` takes the route too (the first
  placement, in `frm()`, left it out: `cooks.distance()` at a missing row read
  3.5e-4; now 0, and the observed rows agree with the complete-case fit to
  8.7e-10). `mi_rows_read()` decides row by row: a predictor's `mi()`, with
  `idx` read row by row; every row under `mi(sd)`, a `cov = FALSE` autocor term,
  or a structure's own likelihood. The objective masks the rest through a
  weight, the rows stay in the frame, and `logLik`'s `nobs` counts the rows it
  read, so `BIC` equals the complete-case one. `frm(dry_run = "objective")` and
  so `frm_sample()` are not routed. `frm_sample(fit)` restores them through
  `mi_sampling()`, with the read rows at their estimates. On student and the
  other pre-existing case this is BREAKING and NEWS says so.
- **M1, M2.** `check_id_mm_weights()`, `check_id_cor_flags()`; both orders in
  tests. A key's weights now follow one spelling; `(1 | p | mm(g1, g2))` beside
  the same weights still fits.
- **M3, m10.** `mi_placeholder()` is the observed mean on every chart. A mixture
  takes its components' chart when they share one, so
  `frm_sample()` holds lognormal, Gamma and beta mixtures in their support;
  `no_laplace` is a separate flag. Draws of a lognormal mixture against brms
  (`dev/refuse-f1-draws.R mix`, seed 9700, overlapping components, 4 chains of
  2000, each component held by a prior): KS of `Ymi` median 0.028, max 0.052;
  brms's own Rhat 1.09. Two well-separated components do not mix on either
  side: a missing value's two modes are the reason, brms's Rhat 5.5.
- **M4.** `mi_chart_inv()` stops half a unit in the last place short of each
  end,
  and `mi_chart_logjac()` reads `|u|`, so it does not overflow.
  `dev/refuse-f1-m4.R` (seed 10200): the 24 fits run on the unread route, on
  the complete-case fit (equal logLik in all) and on a read `mi()` with no error
  in any; 3 of the read fits report code 1, at `phi = 0.3` near each end, and
  one on the negative side and on the complete-case route also, so no longer one
  side.
- **M5.** The seam is `mi_sampling(fit, what)`: one export, because
  frmtmb.sample reads no internal, and because both halves of it need core's
  chart (`"fit"` restores the frame, `"natural"` gives `draws_natural_cols()`
  the linkinv and linkfun per column, as `gp_brms_natural()` does). Draws are
  `Ymi_<resp>[k]`, `Ymi[k]` with one response, `Yl` under `mi(sd)`, brms's names
  (`dev/refuse-f1-draws.R`). The sampler's Jacobian is unchanged. KS against
  brms's `Ymi`, 15 values, 4 chains of 2000 each: Gamma median 0.018, max 0.034;
  von_mises median 0.021, max 0.030; `posterior_epred` and `posterior_predict`
  read them, `log_lik` is the same refusal, `posterior_summary` lists them.
- **m1.** NEWS and the compat row say "one offset". **m2.** `gr_cov_matrix()` in
  brms's words. **m3.** The two guards need a latent value. **m4.** `||` over
  `mm(cov =)` splits as over `gr(cov =)`, and `gr_cov_check_within()` now has a
  test; mutant `mm_cov_arg` fails it. **m6 (withdrawn in fix round 2, N1).** `frm()` read von_mises on the
  line,
  so its read rows have no window: `dev/refuse-f1-m6.R` (seed 10201) reads
  Laplace error -0.006, -0.056, -0.041 at mean 3.14 kappa 20, mean 3.05 kappa 2,
  mean 0 kappa 2 (it was -0.80, -0.75 and small); the sampling route keeps the
  window. **m7.** `mi_chart_fun()` clamps. **m8.** The gated file has L1 tests
  on
  `mm(cov =)`, `se()` on skew_normal, `cens()` on a mixture, and the review's
  form of the identity. Mutants `gr_cov_ident_l1`, `mix_lcdf_l1`, `se_resid_l1`
  all fail it, 2 expectations each (`dev/refuse-f1-out/mut/`). **m9.** brms-45,
  brms-46, brms-47.
- **m5, filed.** `asym_laplace` as the response of a read gaussian `mi()`:
  review reads +0.21 mean, up to +0.59, per value, inner gradient 0.99.
  Pre-existing, unchanged; it could be refused the same way.
- **NITs.** Absolute tolerances in test-mi-families.R gone; NEWS within 80
  columns, "L1 constant" gone; the misnamed test renamed; `mi_chart` refused by
  name; a response with every value missing refused by name. `id =` beside
  `|p|` stays refused, and **`s(I(x^2))` is refused where brms builds it**, the
  same on rellib-r16, filed here.

## 8a. Fix round 2

The short check's verdict was FIX (`dev/reviews/2026-11-02-refuse.md`, "Short
check of fix round 1"). What round 1 said about von_mises on the line was wrong:
every latent left on `frm()` is read linearly, so periodicity does not help.
Scripts `dev/refuse-f2-*.R`, logs `dev/refuse-f2-out/`. Mutants of the new tests
(`dev/refuse-f2-mutants.R`, `mutants.log`): `one_pass`, `weight_reads`,
`gone_reads`, `no_warn`, `empty_miss`, `placeholder` all fail.

- **N1.** The window is back on `frm()`: `mi_ml_frame()` no longer maps a chart.
  `dev/refuse-f2-vm2.R`, seed 11207, the same 13 replicates: chart `circle`,
  values inside the window in every replicate (widest -2.91 to 2.23), largest
  shift from the exact maximum on the window 0.74 SE, as the short check read
  with the window restored in-process. `logLik` lands below the exact maximum by
  2.9 to 5.2, never above. Seed 11400: L1 reads again at all three shapes: mean 0
  kappa 3 const -5.68e-14, mean 3.05 kappa 3 const 0, mean 3.14 kappa 20 const
  9.38e-13, each `equivalent`. The gated file has the middle shape as a test.
  The Laplace loss near the end is in NEWS and `?bf`, measured on the window
  (`dev/refuse-f2-m6.R`, seed 10201, reader slope 0.2): -1.22 at mean 3.14 kappa
  20, -0.98 at mean 3.05 kappa 2, -0.068 at mean 0. That is the review's -0.80 and
  -0.75, which had no reader, in the same range. `dev/refuse-f1-m6.R` reads smaller
  numbers on this build and should not be used: it integrates one period on the
  line, which misses the far end of the chart. The read arm's von_mises per value
  is now -0.076 (`dev/refuse-f1-unread.R von_mises`), where the line read -0.050.
- **N2.** `mi_rows_read()` peels from the leaves to a fixed point, with
  `mi_rows_live()` for each reader's rows. `dev/refuse-f2-chain.R`, seeds 9503,
  12000 to 12003, against hand deletion, `logLik` gap and largest coefficient gap:
  Gamma chain 0, 0, 4 latent; the review's exgaussian 0, 0, log sigma -0.404, log
  beta 0.689 (read -0.003, 0.484); a chain of three 1.14e-13, 9.04e-12; a fork
  1.71e-13, 2.40e-11; a reader weight of 0 on every missing row, latent 0, gap
  -4.69e-09 against the separate complete-case fits; on 6 of 10 rows, 0, 0. A
  cycle has a reader at every step, so it stays latent, unrefused. BIC's count is
  now `mi_rows_gone()`, the rows with no density in any response: with `y` and `z`
  both gone on a row, the attr is the hand deletion's 114, not 108. `mi_values()`
  gives a reader the value's `fitted()` mean at a row that left, so `fitted(z)`
  there is not read at the placeholder, gap 1.11e-16.
- **N3.** `mi_sampling(fit, "laplace")`, no new export. `dev/refuse-f2-n3.R`, seed
  10600, sigma and beta: formula 0.360, 2.132 (was 0.910, 1.590); fit 0.355,
  2.146; observed rows, fully sampled, 0.359, 2.139. A read value with no mode is
  refused on this route by name.
- **N4.** The `id_no_cor` record keys `mm()` on `mm_group_name()`, as the terms it
  meets do: all four orders and spellings are refused; all terms saying
  `cor = FALSE` fits.
- **N5.** `mi_laplace_note()`, after the fit on the user's own call. See below.
- **N6.** `mi_values()` guards an empty `miss`: `fitted`, `predict`, `residuals`
  run on both responses (`dev/refuse-f2-chain.R`, seed 11600).
- **N7.** `check_id_mm_weights()` moved to frame assembly and compares the resolved
  matrices: `weights = W` beside `cbind(w1, w2)` fits, gap to `W` twice 0;
  `cbind(w2, w1)` refused.
- **N8.** `importance` on an all-unread fit runs, by the same mask in
  `build_importance_objective()`'s weights: `logLik` and coefficient gaps to the
  complete-case fit at one seed 0, 0. `vcov(cluster =)` on an all-unread fit
  runs, gap to the complete-case fit 0, and its refusal now fires only on a latent
  value, so the reason holds. `mi_sampling(fit, "bad")` refused by
  `frm_match_arg()`. `id = "p"` with `cor = FALSE` over another group refused in
  brms's words in all three brms refuses (`dev/refuse-f2-brms-probe.R`); brms
  builds the `|p|` spelling, because `cor = FALSE` there drops the key, and so does
  this.

### N5: the warning, calibrated

The total is the sum over read values of `log(Laplace) - log(integrate())`, one
value at a time at the estimate. It reproduces the short check's gap to its exact
maximum: Gamma 60/0.7 -5.48 against -5.477, exgaussian 60/0.05 -7.65 against
-7.652 on every row. The threshold stays at 1. Warned against shift, on the short
check's exact table (`dev/refuse-f2-laplace.R`, seed 11200 + family, 200 rows; the
total below is scaled from 30 rows):

```
gaussian read 60 slope 0.70: total NA, warned FALSE, largest exact-SE shift 0.00
gaussian read 60 slope 0.05: total NA, warned FALSE, largest exact-SE shift 0.00
gaussian read 100 slope 0.70: total NA, warned FALSE, largest exact-SE shift 0.00
gaussian read 15 slope 0.70: total NA, warned FALSE, largest exact-SE shift 0.00
exponential read 60 slope 0.70: total -8.11, warned TRUE, largest exact-SE shift 0.06
exponential read 60 slope 0.05: total -5.16, warned TRUE, largest exact-SE shift 0.12
exponential read 100 slope 0.70: total -13.4, warned TRUE, largest exact-SE shift 0.17
exponential read 15 slope 0.70: total -1.89, warned TRUE, largest exact-SE shift 0.03
Gamma read 60 slope 0.70: total -5.35, warned TRUE, largest exact-SE shift 0.50
Gamma read 60 slope 0.05: total -2.76, warned TRUE, largest exact-SE shift 0.30
Gamma read 100 slope 0.70: total -9.37, warned TRUE, largest exact-SE shift 0.96
Gamma read 15 slope 0.70: total -1.43, warned TRUE, largest exact-SE shift 0.13
beta read 60 slope 0.70: total -4.89, warned TRUE, largest exact-SE shift 0.48
beta read 60 slope 0.05: total -6.27, warned TRUE, largest exact-SE shift 0.54
beta read 100 slope 0.70: total -9.07, warned TRUE, largest exact-SE shift 0.98
beta read 15 slope 0.70: total -1.38, warned TRUE, largest exact-SE shift 0.12
lognormal read 60 slope 0.70: total -1.5, warned TRUE, largest exact-SE shift 0.19
lognormal read 60 slope 0.05: total 3.5, warned TRUE, largest exact-SE shift 1.31
lognormal read 100 slope 0.70: total -2.02, warned TRUE, largest exact-SE shift 0.44
lognormal read 15 slope 0.70: total -0.26, warned FALSE, largest exact-SE shift 0.06
von_mises read 60 slope 0.70: total -4.37, warned TRUE, largest exact-SE shift 0.27
von_mises read 60 slope 0.05: total -5.67, warned TRUE, largest exact-SE shift 0.62
von_mises read 100 slope 0.70: total -7.38, warned TRUE, largest exact-SE shift 0.70
von_mises read 15 slope 0.70: total -0.919, warned FALSE, largest exact-SE shift 0.07
student read 60 slope 0.70: total -2.12, warned TRUE, largest exact-SE shift 0.71
student read 60 slope 0.05: total -7.75, warned TRUE, largest exact-SE shift 2.14
student read 100 slope 0.70: total -1.34, warned TRUE, largest exact-SE shift 0.57
student read 15 slope 0.70: total -0.472, warned FALSE, largest exact-SE shift 0.16
exgaussian read 60 slope 0.70: total 0.0991, warned FALSE, largest exact-SE shift 0.14
exgaussian read 60 slope 0.05: total -7.41, warned TRUE, largest exact-SE shift 2.74
exgaussian read 100 slope 0.70: total -0.728, warned FALSE, largest exact-SE shift 0.13
exgaussian read 15 slope 0.70: total 0.0715, warned FALSE, largest exact-SE shift 0.04
```

The largest shift that stayed silent on this table is 0.16 SE; every shift above
1 SE warned. On the second table (`dev/refuse-f2-misc.R`, seed 11501, shift from
the observed rows) no replicate with a shift above 1 SE stayed silent in the final
run; `SILENT` would mark it. On the first run with every row integrated one did:
slope 0.7, 60 read, total -0.90, shift from the observed rows 1.34, from the
exact maximum 0.23, so a reader with information, not a Laplace error.

```
MISC slope 0.00 read 20: shift SE median 1.11 max 1.65; |total| median 5.58 min 3.08; warned 8 of 8
MISC slope 0.00 read 60: shift SE median 3.36 max 4.35; |total| median 9.07 min 4.42; warned 8 of 8
MISC slope 0.02 read 20: shift SE median 1.25 max 1.54; |total| median 5.72 min 4.19; warned 8 of 8
MISC slope 0.02 read 60: shift SE median 2.70 max 3.87; |total| median 7.35 min 3.23; warned 8 of 8
MISC slope 0.05 read 20: shift SE median 1.43 max 1.90; |total| median 5.78 min 3.30; warned 8 of 8
MISC slope 0.05 read 60: shift SE median 3.01 max 3.74; |total| median 7.36 min 4.86; warned 8 of 8
MISC slope 0.10 read 20: shift SE median 1.09 max 1.57; |total| median 5.41 min 2.57; warned 8 of 8
MISC slope 0.10 read 60: shift SE median 3.08 max 3.96; |total| median 7.73 min 6.32; warned 8 of 8
MISC slope 0.20 read 20: shift SE median 1.15 max 2.05; |total| median 5.58 min 3.54; warned 8 of 8
MISC slope 0.20 read 60: shift SE median 2.51 max 3.07; |total| median 7.19 min 5.88; warned 8 of 8
MISC slope 0.40 read 20: shift SE median 0.37 max 0.95; |total| median 1.48 min 0.51; warned 5 of 8
MISC slope 0.40 read 60: shift SE median 0.79 max 1.35; |total| median 3.09 min 0.88; warned 7 of 8
MISC slope 0.70 read 20: shift SE median 0.17 max 0.73; |total| median 0.23 min 0.00; warned 0 of 8
MISC slope 0.70 read 60: shift SE median 0.47 max 1.34; |total| median 0.87 min 0.03; warned 4 of 8
```

**It over-warns.** A family whose per-value error is nearly constant in the
parameters warns at small shifts: exponential 15 read -1.89 at 0.03 SE, Gamma 15
read -1.43 at 0.13 SE, von_mises 60 read about -4.2 at 0.27 to 0.74 SE. The review
asked that a strong reader at 0.14 SE not warn, and on exgaussian it does not, but
on these families almost any fit with about 12 read values or more will warn. A
statistic on how the total moves with the parameters would separate them. I did not
build it: it needs a second sweep at a shifted parameter, and the threshold on the
total is what was asked for.

**Cost**, `dev/refuse-f2-cost.R`, seed 11502, medians of 3 on a loaded machine:

```
COST n=200 read=60: whole fit 0.83s, the check 0.30s (medians of 3), share of the whole 36.1%; 30 rows, 104 sweeps per row; total -1.96
   every row at this tolerance: total -3.032; the sample's scaled total -1.962
COST n=200 read=120: whole fit 1.17s, the check 0.36s (medians of 3), share of the whole 30.8%; 30 rows, 86 sweeps per row; total -5.45
   every row at this tolerance: total -5.230; the sample's scaled total -5.449
COST n=1000 read=300: whole fit 3.89s, the check 1.34s (medians of 3), share of the whole 34.4%; 30 rows, 116 sweeps per row; total -32.41
   every row at this tolerance: total -29.342; the sample's scaled total -32.408
COST n=5000 read=1500: whole fit 31.04s, the check 10.64s (medians of 3), share of the whole 34.3%; 30 rows, 114 sweeps per row; total -159.14
```

It is about a third of the whole fit at every size, since both scale with the
tape. That is the row cap of 30 at about 100 sweeps per value, set after every row
cost two thirds (`dev/refuse-f2-rule.R`: a looser `integrate()` tolerance moves the
total by at most 0.04, on eight shapes). The cap's price is sampling noise on a
heterogeneous family: n = 200, 60 read, scaled -1.96 against -3.03 with every
row. So "each read row's error" holds up to 30 rows and is an even sample above.

## 8b. Not done, and defects found

- The `lb = 0` of the `sd` and `sigma` class rows: filed at 0.78.0, untouched.
- A translator rule for `sigma ~ s(z)` and for `tmp_xi`: both named above.
- `class = "sd"` on `car()`/`spde()` keeps its old spelling: brms's `sdcar`
  stays refused by name. Not in this brief.

## 9. Generated tables and run counts

### Sweep: frm() and brms 2.23.0 per spec (dev/refuse-sweep.R, seed 4040)

| spec | rule | frm() | brms |
|---|---|---|---|
| cens_mixture | cens() x mixture | fits | builds, compiles |
| grcov_diag | gr_cov with \|\| and several columns | fits | builds, compiles |
| id_diag | \|ID\| x kind:covstruct (gr(id, cor = FALSE)) | fits | builds, compiles |
| mi_pred_gaussian_student_x | mi_pred under poisson | fits | builds, compiles |
| mi_sdy_se | mi(sdy) x se() | fits | builds, compiles |
| mi_se | mi() x se() | fits | builds, compiles |
| ordinal_mvbf | ordinal x mvbf | fits | builds, compiles |
| se_family | se() x kind:family (skew_normal) | fits | builds, compiles |
| thres_mvbf | thres() x mvbf | fits | builds, compiles |
| grby_grcov | gr_by x gr_cov | refused | builds, compiles |
| id_rr_na | \|ID\| x group:student_blocks | refused | builds, compiles |
| mi_cens | mi() x cens() | refused | builds, compiles |
| mi_rescor | mi() x rescor | refused | builds, compiles |
| mi_trunc | mi() x trunc() | refused | builds, compiles |
| mixture_mi | mixture x mi() | refused | builds, compiles |
| mixture_se | mixture x se() | refused | builds, compiles |
| mm_pw | mm(pw = ) | refused | builds, compiles |
| mm_student | mm() x group:student_blocks | refused | builds, compiles |
| mmby_cov | mm(by = , cov = ) | refused | builds, compiles |
| mo_mi | mo() x mi_pred() | refused | builds, compiles |
| rescor_se | rescor x se() | refused | builds, compiles |
| rescor_weights | rescor x weights() | refused | builds, compiles |
| trunc_mixture | trunc() x mixture | refused | builds, compiles |
| autocor_aterm_cens | group:autocor x kind:aterm (cens) | refused | refuses |
| autocor_aterm_trunc | group:autocor x kind:aterm (trunc) | refused | refuses |
| autocor_aterm_weights | group:autocor x kind:aterm (weights) | refused | refuses |
| autocor_mixture | group:autocor x mixture | refused | refuses |
| autocor_rescor | group:autocor x rescor | refused | refuses |
| autocor_se_cov | ar() x se() (cov = TRUE) | refused | refuses |
| autocor_se_nocov | ar() x se() (cov = FALSE) | refused | refuses |
| autocor_two | group:autocor x kind:autocor | refused | refuses |
| bar_crossing | bar_crossing x * | refused | refuses |
| cens_categorical | cens() x categorical | refused | refuses |
| cens_disc | p010 cens() x disc | refused | refuses |
| cens_equidistant | p011 cens() x equidistant | refused | refuses |
| cens_ordinal | cens() x group:ordinal | refused | refuses |
| cens_rate | rate() x kind:family (gaussian) | refused | refuses |
| cens_simplex | cens() x group:simplex | refused | refuses |
| cs_family | cs_pred() x kind:family | refused | refuses |
| disc_gaussian | disc x kind:family | refused | refuses |
| mixture_rescor | mixture x rescor | refused | refuses |
| rescor_cens | rescor x cens() | refused | refuses |
| rescor_student_mix | rescor x kind:family (poisson) | refused | refuses |
| rescor_trunc | rescor x trunc() | refused | refuses |
| se_lognormal | se() x kind:family (lognormal) | refused | refuses |
| subset_me | subset() x me() | refused | refuses |
| subset_rescor | subset() x rescor | refused | refuses |
| thres_gaussian | thres() x kind:family | refused | refuses |
| trials_cs | p146 trials() x cs_pred() | refused | refuses |
| trials_cs2 | p146 trials() x cs_pred() on ordinal | refused | refuses |
| trunc_categorical | trunc() x categorical | refused | refuses |
| trunc_ordinal | trunc() x group:ordinal | refused | refuses |
| trunc_simplex | trunc() x group:simplex | refused | refuses |
| vint_stz | p059 vint() x sum_to_zero | refused | refuses |
| zero_intercept_ord | 0 + Intercept x group:ordinal | refused | refuses |

Counts: 55 specs; brms builds and compiles 23; of those frm() refuses 14.

### D8 per family (dev/refuse-mi-l1.R, seed 169 + family index)

| family | L1 response const | inner grad | L1 read by mi() const | inner grad | logLik - exact, per missing | max coef shift, SE |
|---|---|---|---|---|---|---|
| asym_laplace | - | - | - | - | - | - |
| beta | 3.9932502e-12 | 5.6310512e-13 | -1.6910917e-12 | 8.7008178e-13 | -0.0289447  | 0.0016 |
| exgaussian | 2.8421709e-14 | 5.1910254e-11 | 0.0000000e+00 | 1.2597035e-11 | -0.08473924  | 0.1 |
| exponential | -5.6843419e-14 | 3.3236391e-10 | -1.1368684e-13 | 2.7839064e-11 | -0.08106147  | 4.3e-12 |
| frechet | -5.6843419e-14 | 3.0172451e-09 | 2.2737368e-13 | 2.7285507e-11 | -0.08106147  | 2.2e-06 |
| Gamma | 1.4210855e-13 | 7.6576523e-12 | 1.7053026e-13 | 2.7580871e-10 | -0.02563259  | 1e-05 |
| gaussian | -1.1368684e-13 | 7.9682883e-14 | -1.1368684e-13 | 3.3306691e-16 | 0  | 3.4e-13 |
| gen_extreme_value | no translation rule for Stan parameter t | - | argument "family" is missing, with no de | - | -0.121734  | 0.092 |
| inverse.gaussian | -1.1368684e-13 | 1.7298829e-11 | 5.6843419e-14 | 3.3367020e-10 | 0.04062295  | 0.022 |
| lognormal | -1.7053026e-13 | 1.2012613e-13 | 0.0000000e+00 | 3.4209213e-10 | 7.579123e-15  | 1.8e-13 |
| skew_normal | -3.4106051e-13 | 1.1616672e-13 | -1.7053026e-13 | 1.0360601e-12 | -0.01229713  | 0.024 |
| student | 1.1368684e-13 | 6.1437288e-17 | 2.2737368e-13 | 1.1421974e-12 | -0.1645297  | 0.058 |
| von_mises | 1.4210855e-13 | 3.5602854e-11 | -1.7053026e-13 | 6.5797368e-13 | -0.06404186  | 0.012 |
| weibull | 5.6843419e-14 | 1.0255685e-10 | 0.0000000e+00 | 1.0699441e-11 | -0.08106147  | 5.2e-06 |

### brms families, mi() response (dev/refuse-brms-mi.R, seed 169)

mi-families.tsv: builds and compiles: asym_laplace, beta, exgaussian, exponential, frechet, gamma, gaussian, gen_extreme_value, inverse.gaussian, lognormal, skew_normal, student, von_mises, weibull

mi-families.tsv: refused: acat, bernoulli, beta_binomial, binomial, categorical, com_poisson, cox, cratio, cumulative, dirichlet, dirichlet_multinomial, dirichlet2, discrete_weibull, geometric, hurdle_cumulative, hurdle_gamma, hurdle_lognormal, hurdle_negbinomial, hurdle_poisson, logistic_normal, multinomial, negbinomial, negbinomial2, poisson, shifted_lognormal, sratio, wiener, xbeta, zero_inflated_asym_laplace, zero_inflated_beta, zero_inflated_beta_binomial, zero_inflated_binomial, zero_inflated_negbinomial, zero_inflated_poisson, zero_one_inflated_beta, acat, asym_laplace, bernoulli, beta, beta_binomial, binomial, categorical, com_poisson, cox, cratio, cumulative, dirichlet, dirichlet_multinomial, dirichlet2, discrete_weibull, exgaussian, exponential, frechet, gamma, gaussian, gen_extreme_value, geometric, hurdle_cumulative, hurdle_gamma, hurdle_lognormal, hurdle_negbinomial, hurdle_poisson, inverse.gaussian, logistic_normal, lognormal, multinomial, negbinomial, negbinomial2, poisson, shifted_lognormal, skew_normal, sratio, student, von_mises, weibull, wiener, xbeta, zero_inflated_asym_laplace, zero_inflated_beta, zero_inflated_beta_binomial, zero_inflated_binomial, zero_inflated_negbinomial, zero_inflated_poisson, zero_one_inflated_beta

mi-families-pred.tsv: builds and compiles: acat, asym_laplace, bernoulli, beta, beta_binomial, binomial, categorical, com_poisson, cox, cratio, cumulative, discrete_weibull, exgaussian, exponential, frechet, gamma, gaussian, gen_extreme_value, geometric, hurdle_cumulative, hurdle_gamma, hurdle_lognormal, hurdle_negbinomial, hurdle_poisson, inverse.gaussian, lognormal, negbinomial, negbinomial2, poisson, shifted_lognormal, skew_normal, sratio, student, von_mises, weibull, xbeta, zero_inflated_asym_laplace, zero_inflated_beta, zero_inflated_beta_binomial, zero_inflated_binomial, zero_inflated_negbinomial, zero_inflated_poisson, zero_one_inflated_beta

mi-families-pred.tsv: refused: dirichlet, dirichlet_multinomial, dirichlet2, logistic_normal, multinomial, wiener

### Variants (dev/refuse-variants*.sh, dev/equiv-variant.R)

| variant | verdict | l1_const | extra points | l1_grad | l1_error |
|---|---|---|---|---|---|
| cens_mix | defect | -3.01 | [1:2] -2.85 -3.19 | 2.47 |  |
| cens_mix_none | equivalent | -5.68e-14 | [1:2] -5.68e-14 0.00 | 0.000259 |  |
| flat_sigma_smooth | equivalent |  |  |  | 'data' must be of a vector type, was 'NULL' |
| id_diag | defect |  |  |  | the leading minor of order 3 is not positive |
| mi_gev | equivalent | 2.84e-14 | [1:2] 2.84e-14 1.14e-13 | 6.91e-10 |  |
| mi_lognormal_sd | equivalent | 3.84e-13 | [1:2] 3.98e-13 4.12e-13 | 7.2e-12 |  |
| mi_sd_se | equivalent | -1.14e-13 | [1:2] -5.68e-14 0.00 | 4e-15 |  |
| mi_se | equivalent | -3.98e-13 | [1:2] -1.14e-13 -2.84e-13 | 1.27e-15 |  |
| mix_control | equivalent | 0 | [1:2] 0.00 -1.71e-13 | 0.000743 |  |
| sds_beside_sd | equivalent | 2.84 | [1:2] 2.84 2.84 | 1.18e-13 |  |
| sds_by | equivalent | -1.97 | [1:2] -1.97 -1.97 | 1.4e-13 |  |
| sds_sigma | equivalent |  |  |  | the two bases do not span the same space |
| sds_sigma2 | equivalent |  |  |  | 'data' must be of a vector type, was 'NULL' |
| sds_t2 | equivalent | -0.848 | [1:2] -0.848 -0.848 | 2.13e-14 |  |
| sds_t2coef | equivalent | -0.848 | [1:2] -0.848 -0.848 | 2.13e-14 |  |
| sds_two | defect | 2.13 | [1:2] 2.18 2.08 | 1.06e-13 |  |
| sds_two2 | equivalent | 1.38 | [1:2] 1.38 1.38 | 1.69e-13 |  |
| skew_se | equivalent | 2.27e-13 | [1:2] -1.14e-13 -1.14e-13 | 1.48e-05 |  |
| skew_se_sigma | equivalent | -1.71e-13 | [1:2] 1.14e-13 -2.84e-14 | 1.23e-05 |  |

### mi() predictor under each family (dev/refuse-mi-pred.R, seed 170)

```
poisson                frm: logLik -438.7092  const 0 inner grad 3.89e-13
bernoulli              frm: logLik -277.8472  const -5.68e-14 inner grad 1.03e-15
binomial               frm: logLik -496.2902  const -1.71e-13 inner grad 3.1e-13
negbinomial            frm: logLik -475.0247  const -1.59e-12 inner grad 7.38e-12
cumulative             frm: logLik -399.9421  L1 error: subscript out of bounds
sratio                 frm: logLik -402.7065  L1 error: subscript out of bounds
categorical            frm: logLik -394.8137  const 0 inner grad 5.44e-14
Gamma                  frm: logLik -325.8776  const -1.14e-13 inner grad 1.33e-13
beta                   frm: logLik -59.9426  const -4.65e-12 inner grad 3.11e-14
student                frm: logLik -467.4926  const 1.14e-13 inner grad 5.95e-12
hurdle_gamma           frm: logLik -392.0517  const 0 inner grad 1.2e-12
zero_inflated_poisson  frm: logLik -421.5920  const -1.14e-13 inner grad 1.15e-13
von_mises              frm: logLik -300.6847  const 5.68e-14 inner grad 2.01e-13
asym_laplace           frm: logLik -445.9234  const -1.3e-08 inner grad 0.898
```

### B1 and the read arm, per family (dev/refuse-f1-unread.R)

```
family gaussian 
UNREAD logLik -166.246834 complete-case -166.246834 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const -2.84e-14 grad 2.84e-14 
READ chart=identity Laplace - exact per value: mean +0.0000 range [+0.0000, +0.0000] conv 0
READ L1 const 0 inner grad 3.89e-16 
family student 
UNREAD logLik -231.215637 complete-case -231.215637 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const -5.68e-14 grad 8.88e-06 
READ chart=identity Laplace - exact per value: mean -0.0000 range [-0.0000, -0.0000] conv 0
READ L1 const 2.23e-08 inner grad 2.22e-16 
family skew_normal 
UNREAD logLik -172.723879 complete-case -172.723879 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const -2.84e-14 grad 0.00024 
READ chart=identity Laplace - exact per value: mean -0.0227 range [-0.0552, +0.0802] conv 0
READ L1 const -2.84e-14 inner grad 6.44e-11 
family beta 
UNREAD logLik 36.150581 complete-case 36.150581 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const 3.2e-13 grad 0.000158 
READ chart=logit Laplace - exact per value: mean -0.1253 range [-0.1450, -0.1205] conv 0
READ L1 const 0 inner grad 4.11e-12 
family Gamma 
UNREAD logLik -213.993095 complete-case -213.993095 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const -8.53e-14 grad 0.000326 
READ chart=log Laplace - exact per value: mean -0.1145 range [-0.1710, -0.0631] conv 0
READ L1 const 5.68e-14 inner grad 3.29e-10 
family weibull 
UNREAD logLik -171.616881 complete-case -171.616881 gap 0 max coef gap 1.11e-16 nmiss slots 0 BIC gap 0
UNREAD L1 const 0 grad 4.41e-05 
READ chart=log Laplace - exact per value: mean -0.0783 range [-0.0946, -0.0374] conv 0
READ L1 const 0 inner grad 3.78e-10 
family exponential 
UNREAD logLik -167.303783 complete-case -167.303783 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const 8.53e-14 grad 0.000778 
READ chart=log Laplace - exact per value: mean -0.0980 range [-0.1800, -0.0068] conv 0
READ L1 const 0 inner grad 2.33e-09 
family frechet 
UNREAD logLik -200.334926 complete-case -200.334926 gap 1.42e-13 max coef gap 1.99e-08 nmiss slots 0 BIC gap -
UNREAD L1 const -1.14e-13 grad 3.83e-06 
READ chart=log Laplace - exact per value: mean +0.0082 range [-0.0711, +0.0883] conv 0
READ L1 const 5.68e-14 inner grad 2.88e-10 
family inverse.gaussian 
UNREAD logLik -145.974649 complete-case -145.974649 gap 6.54e-13 max coef gap 6.29e-08 nmiss slots 0 BIC gap -
UNREAD L1 const -2.84e-14 grad 1.35e-06 
READ chart=log Laplace - exact per value: mean +0.0211 range [-0.0770, +0.0793] conv 0
READ L1 const 0 inner grad 3.6e-11 
family lognormal 
UNREAD logLik -243.979579 complete-case -243.979579 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const -1.14e-13 grad 2.98e-11 
READ chart=log Laplace - exact per value: mean -0.0367 range [-0.1568, +0.0325] conv 0
READ L1 const 0 inner grad 9.14e-10 
family exgaussian 
UNREAD logLik -197.469053 complete-case -197.469053 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const -8.53e-14 grad 2.08e-08 
READ chart=identity Laplace - exact per value: mean -0.0793 range [-0.1550, +0.1785] conv 0
READ L1 const 2.84e-14 inner grad 9.8e-10 
family gen_extreme_value 
UNREAD logLik -190.982657 complete-case -190.982657 gap -8.23e-11 max coef gap 4.17e-07 nmiss slots 0 BIC gap 
UNREAD L1 L1 error: no translation rule for Stan parameter tmp_xi 
READ chart=identity Laplace - exact per value: mean +0.0075 range [-0.0150, +0.0279] conv 0
READ L1 L1 error: no translation rule for Stan parameter tmp_xi_y 
family von_mises 
UNREAD logLik -144.578297 complete-case -144.578297 gap -5.46e-12 max coef gap 1.54e-07 nmiss slots 0 BIC gap 
UNREAD L1 const 8.53e-14 grad 8.2e-05 
READ chart=identity Laplace - exact per value: mean -0.0505 range [-0.1397, -0.0356] conv 0
READ L1 const -2.84e-14 inner grad 4.84e-11 
family asym_laplace 
UNREAD logLik -206.213923 complete-case -206.213923 gap 0 max coef gap 0 nmiss slots 0 BIC gap 0
UNREAD L1 const 0 grad 2.07 
READ refused: frm() cannot fit mi() on `y`: another density reads its missing values, so they are integrated o
```

### Tests on the final build

```
One file per process, on the final build (R installed 14:04; only help text and one test file edited after):
core            files=292 pass=20813 fail=0 error=0 skip=216 warn=0
frmtmb.sample   files=82 pass=3802 fail=0 error=0 skip=18 warn=0
gated tier      files=28 pass=3444 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-mi-rows-read.R pass=54 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-mi-families.R pass=181 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-refuse-pairs.R pass=43 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-mm-id-cov.R pass=49 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-double-bar.R pass=154 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-mm.R pass=117 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-v16.R pass=24 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-compat.R pass=756 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-prior-sds.R pass=52 fail=0 error=0 skip=0 warn=0
RESULT frmtmb.sample test-mi-laplace-route.R pass=3 fail=0 error=0 skip=0 warn=0
RESULT frmtmb.sample test-mi-draws-brms.R pass=0 fail=0 error=0 skip=3 warn=0
RESULT frmtmb.sample test-mm-id-cov-sample.R pass=4 fail=0 error=0 skip=0 warn=0
RESULT frmtmb.sample test-sds-prior.R pass=28 fail=0 error=0 skip=0 warn=0
RESULT frmtmb.sample test-reparam.R pass=373 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-brms-likelihood.R pass=495 fail=0 error=0 skip=0 warn=0
RESULT frmtmb test-mi-mm-brms-l1.R pass=41 fail=0 error=0 skip=0 warn=0
RESULT frmtmb.sample test-mi-draws-brms.R pass=12 fail=0 error=0 skip=0 warn=0
One session, three orders (dev/release/refuse-sess-sessions; no untagged LEAK):
SESSION frmtmb forward files=282/292 pass=20813 fail=0 err=0 skip=216 warn=0
SESSION frmtmb preload files=282/292 pass=21056 fail=0 err=0 skip=218 warn=0
SESSION frmtmb reverse files=282/292 pass=20963 fail=0 err=0 skip=216 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward files=77/82 pass=3802 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.sample preload files=77/82 pass=3802 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.sample reverse files=77/82 pass=3802 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
OpenBLAS with LAPACK: dev/ciharden-log/refuse-ob26.sum, refuse-ob32.sum, 13 of 13,
fail 0; test-mi-mm-brms-l1.R's one warning there was a missing allow_warnings(),
fixed, and rerun alone (refuse-ob26l.sum, refuse-ob32l.sum):
0.3.26
frmtmb RESULT test-mi-families.R pass=181 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-refuse-pairs.R pass=43 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mm-id-cov.R pass=49 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-prior-sds.R pass=52 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mm.R pass=117 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-v16.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-brms-likelihood.R pass=495 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-compat.R pass=756 fail=0 err=0 skip=0 warn=0
frmtmb.sample RESULT test-mi-draws-brms.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb.sample RESULT test-sds-prior.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mi-rows-read.R pass=54 fail=0 err=0 skip=0 warn=0
frmtmb.sample RESULT test-mi-laplace-route.R pass=3 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mi-mm-brms-l1.R pass=41 fail=0 err=0 skip=0 warn=0
0.3.32
frmtmb RESULT test-mi-families.R pass=181 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-refuse-pairs.R pass=43 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mm-id-cov.R pass=49 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-prior-sds.R pass=52 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mm.R pass=117 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-v16.R pass=24 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-brms-likelihood.R pass=495 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-compat.R pass=756 fail=0 err=0 skip=0 warn=0
frmtmb.sample RESULT test-mi-draws-brms.R pass=12 fail=0 err=0 skip=0 warn=0
frmtmb.sample RESULT test-sds-prior.R pass=28 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mi-rows-read.R pass=54 fail=0 err=0 skip=0 warn=0
frmtmb.sample RESULT test-mi-laplace-route.R pass=3 fail=0 err=0 skip=0 warn=0
frmtmb RESULT test-mi-mm-brms-l1.R pass=41 fail=0 err=0 skip=0 warn=0
R CMD check --as-cran:
frmtmb        Status: 2 NOTEs (V8 unavailable; examples over the elapsed limit under load)
frmtmb.sample Status: 1 NOTE (V8 unavailable)
```
