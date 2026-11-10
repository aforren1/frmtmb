# Lane translator: the survey's weak and unmeasured pairs, lifted to L1

(At the 0.81.0 release the lane's brms-49 was numbered brms-51 in
`dev/upstream-bugs.md`, since brms-49 and brms-50 were taken; this
record and the adjudication row say brms-51. The review keeps its
own text.)

Round opened 2026-11-09, base e0cce906 (frmtmb 0.80.0), measured on
`rellib-r19` (every log's `lib:` line reads
`C:/Users/adf44/source/r/rellib-r19/frmtmb 0.80.0`; the lane library is
empty, because no package code changed). This lane changes `dev/` only:
no package code, no test, no record. `dev/equiv-data/` and the generated
blocks of `dev/equiv-findings.md` are untouched. Nothing is committed.
No R CMD check was run, because no package changed.

## 1. Why each pair missed L1 on the base

Measured fresh on rellib-r19 with the base instrument
(`dev/translator-measure.sh dev/translator-before-out 3 <ids>`, the
comparison by `dev/translator-compare.R`, output
`dev/translator-before-cmp.txt`). All 9 reproduce the record's measured
verdict and level (p069's record level is its adjudication's).

| id | pair | what stopped L1 |
|---|---|---|
| p143, p042 | smooth / s() x nl | the suite's `zs` rule reads the basis `Zs_<i>_<j>`; brms names a smooth in nl parameter b1 `Zs_b1_<i>_<j>`, so `as.matrix(NULL)` failed ("'data' must be of a vector type") |
| p047 | t2() x nl | the same |
| p048 | t2() x rescor | the same, with the response as prefix (`Zs_y1_1_1`) |
| p136 | gp x nl | nothing in the translator: L1 read -8.9e-4 and gradient 0.030, the nugget, so it was set aside by the decision and read at L4 |
| p044 | s() x REML | nothing: L1 held (constant 0), but the joint identity holds at ANY point and does not read the integral REML maximizes |
| p069 | us x importance | nothing: L1 on the Laplace fit holds (-1.9e-10); the mode's own estimate has no brms quantity but the posterior |
| p131 | cs_pred() x mixture | the suite's `Intercept_<dpar>` rule has no ordinal-mixture case (brms's `Intercept_mu<k>` are thresholds), and its `bcs` rule reads an unsuffixed `Kcs` and frmtmb's `bcs<j>` in number order; frm() also stops at a degenerate boundary |
| p021 | gen_extreme_value x mixture | no rule for brms's `tmp_xi<k>` |

A second instrument gap showed up once the smooths translated: check C's
inner-block pattern `brms_inner_pat` lists `zs_<i>_<j>` only, so a
prefixed `zs` had no inner index, the gradient read NA, and the first
remeasure scored the three nl smooth identities and p048 as DEFECTS
(`dev/translator-after-out`, first pass). Fixed in the same section.

## 2. What changed in the instrument

All in `dev/`; `tests/testthat/helper-brms.R` (the suite's translator) is
not touched. The new rules wrap it in `dev/equiv-lib.R`, section
"translator extensions, round r081": a parameter a rule owns is removed
from the program text the suite's translator reads (`eq_strip_pars()`,
as the helper already does for a mixture's simplex `theta`) and filled
in afterwards with its share of the log-Jacobian.

- `zs_prefixed`: a smooth outside a univariate mu (`zs_b1_1_1` in an nl
  body, `zs_y1_1_1` under mvbf) reads its basis from
  `Zs_<prefix>_<i>_<j>`; the rest is the suite's rule, including its
  orthogonality check on the basis map. `brms_inner_pat` gains the
  prefixed `zs` (with a `stopifnot()` on both spellings and on `sds_b1_1`
  and `bs_b1`, which must stay out).
- `ord_mix_intercept`: `Intercept_mu<k>` of an ordinal mixture under
  `order = "none"` (brms announces it for ordinal components) is
  component k's thresholds, read from `tau_raw<k>` through that
  component's own family; brms's centering shift enters with the ordinal
  sign. Shared thresholds (`order = "mu"`) stop with a message.
- `bcs_prefixed`: `bcs_mu<k>`, matrix[Kcs_mu<k>, nthres], from the
  `bcs<j>` vectors named in component k's linear-predictor frame.
- `tmp_xi`: brms's `scale_xi()` inverted at frmtmb's mu and sigma;
  an xi outside brms's interval stops with the interval named.
- `eq_jitter_optimizer()` and the spec field `l1_points = list(scale,
  why)`: L1 read at three fixed points (frm()'s start moved by
  `scale * sin(1.7 j + k)`, k = 0, 1, 2) instead of at the optimum, for a
  model whose optimum the identity cannot read. No gradient is scored
  there for a model without random effects (no point is an optimum);
  with random effects the inner gradient is still read. All three points
  must evaluate: one that errors fails the pair (fix round, m4). The
  record says so in the new `l1_at` column. Used by p131 (frm() stops at
  a degenerate boundary, thresholds near 3700, where Stan's gradient is
  NaN), p021 (frm()'s optimum is outside brms's parameter space, section
  4) and, through `eq_sd_collapsed()` (`dev/equiv-specs.R`), by the 17
  pairs whose latent block, or one row of it, is degenerate at the
  estimate (sections 3b and 3c): p047, p048, p076, p139 in their own
  specs, and the 13 of the `eq_degenerate` list.
- `l1_sdat_fix`: a spec hook in `eq_l1()` that changes brms's data by a
  change of variables of its program. No survey spec sets it; only
  `dev/translator-bslogdet.R` does (section 5).
- `l1_marginal` and `eq_laplace_stan()`, `eq_bs_logdet()`: the REML
  integral, section 5.
- `l1_code_fix` and `eq_gp_nugget_fix()` (`dev/equiv-specs.R`): an exact
  gp() spec reads L1 a second time on brms's program with
  `cov[n, n] += 1e-12;` replaced by `cov[n, n] += 1e-12 * square(sdgp);`.
  When that holds, the flow scores the pair "divergence by decision" at
  L1 (new branch in `eq_measure()`); before, these pairs had no scored
  L1 and were adjudicated by hand. p136's spec now samples no posterior
  (`sample = FALSE`); its L4 reading took 190 s and added nothing once
  L1 is read.
- `eq_flatten()` adds the columns `l1_marg`, `l1_marg_ok`, `l1_dec`
  (the patched program's three constants, then its gradient),
  `l1_dec_ok` and `l1_at`. `dev/equiv-summary.R`: a REML pair is weak
  evidence only when `l1_marg_ok` is not TRUE (it was weak by feature
  name alone).

New lane scripts: `dev/translator-measure.sh` (a list of pairs, forced,
into an outdir, at most N at once), `dev/translator-compare.R` (an outdir
against the record), `dev/translator-probe.R` (what brms declares and
where the translator stops), `dev/translator-inspect.R`,
`dev/translator-mutants.R`, `dev/translator-gev-check.R`,
`dev/translator-project.R` (the rates with these results, in a scratch
copy, never the record), and from the fix round
`dev/translator-audit.R` with `-audit.sh` and `-audit-summary.R`
(section 3b) and `dev/translator-bslogdet.R` (section 5).

## 3. The rules can fail: mutants

`dev/translator-mutants.R`, log `dev/translator-mutants.log` (fix-round
run, final instrument). Each mutant breaks one rule, or one t2() block,
in process and reruns the pair; each must turn L1 from pass to fail.
Constants at the three points:

| mutant | pair | L1 ok | constants | marginal |
|---|---|---|---|---|
| bcs: components swapped | p131 | FALSE | -1.39, -0.109, 0.108 | |
| thresholds: components swapped | p131 | FALSE | 7.1, 3.44, -4.28 | |
| tmp_xi: off by 0.1 | p021 | FALSE | -0.49, -0.528, -0.776 | |
| zs: coefficients reversed | p042 | FALSE | -288, -203, -381 | |
| zs: all blocks reversed | p048 | FALSE | -9.83, -2.84, -4.79 | |
| zs: block zs_y1_1_1 alone reversed | p048 | FALSE | -3.89, -0.718, -0.743 | |
| zs: block zs_y1_1_2 alone reversed | p048 | FALSE | -3.31, -1.64, -2.65 | |
| zs: block zs_y1_1_3 alone reversed | p048 | FALSE | -2.13, -0.383, -1.22 | |
| zs: block zs_b1_1_1 alone reversed | p047 | FALSE | -17.7, -22.1, -26.2 | |
| zs: block zs_b1_1_2 alone reversed | p047 | FALSE | -165, -306, -177 | |
| zs: block zs_b1_1_3 alone reversed | p047 | FALSE | -215, -178, -119 | |
| fixed points: third point unreadable | p021 | FALSE (2 points) | -5.68e-14, 5.68e-14, NA | |
| marginal read on the ML fit | p044 | FALSE | 2.84e-14, -2.84e-14, 2.84e-14 | -0.221, 0.979, 1.02 |
| control | p131 | TRUE | 0, 0, 5.68e-14 | |
| control | p021 | TRUE | -5.68e-14, 5.68e-14, 1.71e-13 | |
| control | p042 | TRUE | -1.14e-13, -1.14e-13, -1.14e-13 | |
| control | p048 | TRUE | -1.14e-13, -3.41e-13, -1.14e-13 | |
| control | p047 | TRUE | -1.14e-13, 2.27e-13, 1.14e-13 | |
| control | p044 | TRUE | 0, -5.68e-14, 0 | 1.99e-13, 1.42e-13, -4.55e-13 |

The p044 mutant shows why the joint L1 alone was weak evidence for REML:
on the ML fit the joint identity still holds to 3e-14, and only the
marginal reading fails, by up to 1.02 nats. The "third point unreadable"
mutant moves p021's third fixed point by 3 instead of 0.1, where the
tmp_xi rule stops (xi1 -0.659 outside (-0.0090, 0.0145)); under the
first-round flow the pair passed on the two points left, and now it
fails (fix round, m4).

The p048 and p047 rows were first run at the survey's points (estimate
and +/- 0.05). There the review found the single-block reversals of
p048's blocks 2 and 3 and p047's block 3 read L1 TRUE, and my first
p048 mutant reversed all three blocks, so block 1 alone made it fail
(-34.6). The rows above are at the fixed points the two specs now use.

## 3b. Audit: which L1 pairs cannot see a latent block (fix round, M1)

At the estimate a block whose SD has collapsed adds about 1e-7 to the
predictor, far below the L1 tolerance, so a wrong map of that block
passes L1. `dev/translator-audit.R` tests every latent block of every L1
pair at the pair's own three points. It reverses the block's translated
values (levels, basis coefficients or rows), which leaves the block's
standard-normal prior exactly unchanged and moves only what the
likelihood reads. A block is BLIND when brms's log density moves by less
than the L1 tolerance (1e-6 * max(1, |lp|)) at all three points.

Run on 132 pairs: every pair whose record level holds L1, the
lane's L1 pairs, and p069 (L1 on its Laplace fit). Rebuilt with
`dev/translator-audit-summary.R dev/translator-audit-out`, saved as
`dev/translator-audit-summary.txt`:

    pairs launched 132 ; with a result 132 ; no result:
    pairs without a latent block or not translated: 66
       66 no latent block; points
    pairs with latent blocks: 66 ; blocks: 75 ; untestable (one value): 2
    BLIND blocks: 5 in 4 pairs
    seen by less than 100 x tol: 0 blocks

| id | pair | blind block | max shift | tolerance | SDs at the estimate |
|---|---|---|---|---|---|
| p047 | t2() x nl | zs_b1_1_3 | 1.62e-7 | 1.68e-4 | sds_b1_1 3.67, 2.55, 1.36e-4 |
| p048 | t2() x rescor | zs_y1_1_2 | 9.40e-8 | 5.76e-4 | sds_y1_1 5.04, 4.21e-4, 2.56e-4 |
| p048 | t2() x rescor | zs_y1_1_3 | 5.54e-8 | 5.76e-4 | same |
| p076 | cumulative x smooth | zs_1_1 | 1.18e-9 | 2.71e-4 | sds_1 8.66e-5 |
| p139 | multinomial x us | z_1 | 3.12e-7 | 7.01e-4 | sd_1 3.16e-5, sd_2 0.499 |

The two untestable blocks are the empty `Ymi` of p029 and p169: every
missing value there is unread, so brms's program on the observed rows
has none. Every other block moves brms's density by at least 100 times
the tolerance.

The audit touches 2 existing L1 verdicts (p076 and p139, both
"equivalent, L1" in the record) and the lane's p047 and p048. All four
are fixed the same way: `eq_sd_collapsed()` reads L1 at the fixed points
(`scale = 0.3`), where every SD sits near its start. Re-audited there
(`dev/translator-audit-fixed-out`): 9 blocks, 0 blind; the smallest
shifts are p047 26.2, p048 2.13, p076 4.42 and p139 59.6 nats. All four
keep their verdicts (sections 4 and 6).

### 3c. The audit per row (second fix round, F1)

The audit above reversed whole blocks, so a matrix block with one live
row and one dead row read as seen. The short check found p087 (gr_cov x
mmc(), `(1 + mmc(x1, x2) | mm(g1, g2, cov = A))`): SDs 0.463 and 1.12e-6,
and the mmc() slope's row of z reversed alone still read L1 TRUE at the
survey's points (`dev/review-translator-p087.log`).
`eq_audit_mutations()` in `dev/translator-audit.R` now reverses, besides
the whole block, each row of a matrix block on its own, and, for a
gr(by = ) block, each row within each by-group's levels (`Jby_<i>`). Each
reversal permutes standard-normal entries among themselves, so the prior
stays exact.

No sampled pair has a gr(by = ) block, and none can reach L1 yet: a
constructed one (`Rscript dev/translator-audit.R bycheck <outdir>`: group
A's levels with SD 0.8, group B's with none) stops in
`rstan::unconstrain_pars()` ("sd_1 ... dims declared=(1,2); dims
found=(1)"), because the suite's translator has no rule for gr(by = )'s
sd matrix. The by-group branch is in place for when that rule exists.

The per-row audit of all 132 pairs at the points each spec read before
this round (`dev/translator-audit2-out`, summary
`dev/translator-audit2-summary.txt`):

    pairs with latent blocks: 66 ; blocks: 130 ; untestable (one value): 2
    BLIND blocks: 14 in 13 pairs
    seen by less than 100 x tol: 0 blocks

| id | pair | blind part | max shift | tolerance | at the estimate |
|---|---|---|---|---|---|
| p053 | us x emmeans | z_1 row 2 | 1.88e-5 | 3.24e-4 | sd_1 0.308, 0.251; diag(L_1) 1, 0.00172 |
| p056 | us x predict | z_1 row 2 | 3.81e-6 | 3.01e-4 | sd_1 0.418, 0.0531; diag(L_1) 1, 0.00496 |
| p060 | von_mises x \|ID\| | z_1 row 2 | 6.69e-7 | 2.25e-4 | sd_1 0.176, 0.0439; diag(L_1) 1, 0.00299 |
| p069 | us x importance | z_1 row 2 | 1.97e-5 | 3.28e-4 | sd_1 0.446, 0.0909; diag(L_1) 1, 0.00298 |
| p070 | \|ID\| x mmc() | z_1 rows 2 and 3 | 4.14e-5, 5.60e-7 | 3.30e-4 | sd_1 0.327, 0.102, 0.023; diag(L_1) 1, 0.00249, 0.00455 |
| p087 | gr_cov x mmc() | z_1 row 2 | 1.67e-11 | 3.31e-4 | sd_1 0.463, 1.12e-6; diag(L_1) 1, 0.84 |
| p102 | predict x mmc() | z_1 row 2 | 6.78e-7 | 3.29e-4 | sd_1 0.329, 0.027; diag(L_1) 1, 0.00597 |
| p106 | smooth x mmc() | z_1 row 2 | 1.61e-5 | 3.31e-4 | sd_1 0.153, 0.109; diag(L_1) 1, 0.00422 |
| p111 | sum_to_zero x mmc() | z_1 row 2 | 5.60e-7 | 2.69e-4 | sd_1 0.738, 0.0847; diag(L_1) 1, 0.00349 |
| p114 | us x hsgp | z_1 row 2 | 8.64e-5 | 3.26e-4 | sd_1 0.489, 0.141; diag(L_1) 1, 0.00538 |
| p115 | us x me() | z_1 row 2 | 3.28e-5 | 5.39e-4 | sd_1 0.695, 0.115; diag(L_1) 1, 0.00385 |
| p116 | us x mmc() | z_1 row 2 | 2.69e-5 | 2.99e-4 | sd_1 0.0361, 0.301; diag(L_1) 1, 0.0038 |
| p122 | us x cens() | z_1 row 2 | 2.80e-5 | 2.88e-4 | sd_1 0.511, 0.105; diag(L_1) 1, 0.00424 |

Only p087 is a collapsed SD. In the other 12 every SD is well away from
0, but the correlation sits at +/-1: the second diagonal entry of L_1 is
0.0017 to 0.006, so z's second row (and p070's third) carries almost no
weight, as r = diag(sd) L z. The survey's data have 10 groups, and a
correlation at the boundary is common there. The audit first printed
SDs only; it now prints diag(L) as well, and these 13 were rerun with it
before the specs changed.

The fix is the same as in 3b: all 13 read L1 at the fixed points
(`eq_degenerate` at the end of `dev/equiv-specs.R`, which wraps each
spec in `eq_sd_collapsed()`), where every SD and correlation sits near
its start. Re-audited there (`dev/translator-audit2-fixed-out`): 43
parts, 0 blind; the smallest row shift per pair is 3.51 (p111) to 41.6
(p069) nats. Each blind row reversed alone in the full flow now fails
L1, and each unmutated pair holds (`dev/translator-mutants.log`):

| pair | row reversed: constants | control: constants |
|---|---|---|
| p087 | row 2: -11.2, -7.47, -2.16 (row 1: -79.7, -53.5, -26.5) | 1.14e-13, 0, 0 |
| p053 | -5.84, -12.3, -19.6 | -1.14e-13, 5.68e-14, 0 |
| p056 | -4.55, -1.9, -3.72 | -1.14e-13, -1.14e-13, 1.14e-13 |
| p060 | -27.5, -18.8, -18.8 | 2.84e-13, 2.27e-13, -5.68e-14 |
| p069 | -24.8, -24.2, -41.6 | 0, -1.14e-13, 1.71e-13 |
| p070 | row 2: -7.97, -3.05, -3.13; row 3: -9.71, -9.27, -8.66 | -2.84e-13, 1.71e-13, -5.68e-14 |
| p102 | -5.42, -1.25, -0.432 | -1.14e-13, -1.14e-13, -1.71e-13 |
| p106 | -5.49, -3.93, -6.24 | 0, -1.71e-13, 1.71e-13 |
| p111 | -1.59, -3.51, -2.93 | 0, -5.68e-14, -5.68e-14 |
| p114 | -7.19, -4.67, -3.64 | -1.14e-13, 5.68e-14, -1.14e-13 |
| p115 | -2.82, -3.47, -11.2 | 0, 2.27e-13, 0 |
| p116 | -7.55, -5.6, -4.47 | -5.68e-14, -1.14e-13, 2.27e-13 |
| p122 | -8.33, -5.79, -4.52 | -5.68e-14, -5.68e-14, 5.68e-14 |

The final per-row audit of all 132 pairs, with every fix in place
(`dev/translator-audit3-out`, summary `dev/translator-audit3-summary.txt`):

    pairs launched 132 ; with a result 132 ; no result:
    pairs without a latent block or not translated: 66
       66 no latent block; points
    pairs with latent blocks: 66 ; blocks: 130 ; untestable (one value): 2
    BLIND blocks: 0 in 0 pairs
    seen by less than 100 x tol: 0 blocks

Existing L1 readings the two audits touch: 15 (p076 and p139 in 3b; the
13 above). All are now read at fixed points, and section 6 shows each
keeps its verdict and level. With the lane's p047 and p048, 17 pairs
carry fixed points for a degenerate block.

## 4. Each pair, before and after

Before: the record (`dev/equiv-data/results.tsv`, measured verdict, with
the adjudicated verdict where one exists), reproduced on rellib-r19 in
`dev/translator-before-out` (9 of 9 the same; `dev/translator-before-cmp.txt`).
After: the final instrument, `dev/translator-noreg3-out`
(`dev/translator-noreg3-cmp.txt`; the first round's final run is
`dev/translator-final-out`, and differs only on p047 and p048, then read
at the estimate, and on p069, read at fixed points since section 3c).
Constants are brms's log_prob minus frmtmb's at the three points; the
tolerance is 1e-6 * max(1, |ours|).

| id | pair | record: measured (adjudicated), level | after: verdict, level | L1 constants | gradient |
|---|---|---|---|---|---|
| p143 | smooth x nl | equivalent, L4 (weak) | equivalent, L1 | -8.53e-14, 2.84e-14, 5.68e-14 | 1.53e-13 |
| p042 | s() x nl | equivalent, L4 (weak) | equivalent, L1 | -1.14e-13, -1.14e-13, -1.14e-13 | 1.44e-14 |
| p047 | t2() x nl | equivalent, L4 (weak) | equivalent, L1 at fixed points | -1.14e-13, 2.27e-13, 1.14e-13 | 3.2e-14 (inner) |
| p136 | gp x nl | equivalent, L4 (weak) | divergence by decision, L1 | patched: -4.68e-10, -1.32e-10, -5.4e-12; unpatched -8.89e-4, -8.35e-4, -7.19e-4 | patched 2.74e-5; unpatched 0.0298 |
| p069 | us x importance | defect, L1+L4 (adjudicated equivalent, L4, weak) | unchanged: defect, L1+L4; L1 now at fixed points (3c) | 0, -1.14e-13, 1.71e-13 (the fit without importance; at the estimate -1.88e-10, -1.22e-10, 3.98e-12) | 6e-15 |
| p044 | s() x REML | equivalent, L1 (weak) | equivalent, L1, no longer weak | joint 0, -5.68e-14, 0; marginal 1.99e-13, 1.42e-13, -4.55e-13 | 1.3e-13 |
| p131 | cs_pred() x mixture | defect (cannot measure), L4 | equivalent, L1 at fixed points | 0, 0, 5.68e-14 | not scored (23.6: no point is an optimum) |
| p048 | t2() x rescor | defect (cannot measure), L4 | equivalent, L1 at fixed points | -1.14e-13, -3.41e-13, -1.14e-13 | 7.22e-16 (inner) |
| p021 | gen_extreme_value x mixture | defect (cannot measure), L4 | equivalent, L1 at fixed points; see below | -5.68e-14, 5.68e-14, 1.71e-13 | not scored (47.8) |

For p131, p048 and p021 the record's adjudication rows override the
measurement; they say "cannot measure" for reasons this lane removes.
p136 has no adjudication row. p044's `bs` map between the two packages'
unpenalized columns is the identity here (log|det| = 0,
`dev/translator-probe-p044.log`), so the stated marginal constant is 0.

### p131: the record's gap at the optimizing level is local optima

The L1 identity says frmtmb's and brms's densities are one function, so
the record's L4 reading (frmtmb's maximum 0.814 nats below brms's best of
6 starts) cannot be a model difference. The review measured what it is
(`dev/review-translator-p131ms.R`, log `dev/review-translator-p131ms.log`,
seeds 1000 to 1020 around frm()'s start, sd 1): frm()'s default stops
at -267.932323 with a singular-convergence warning; the best of 21
seeded starts reaches -263.413103, 4.52 nats higher and 3.71 nats above
brms's best of 6 (-267.118266); 2 of the 21 starts end at -267.1183,
brms's value to four decimals. So the two optimizers stop at different
local, degenerate points of one function (thresholds in the thousands on
both sides). frm() warns that the component ran off to a degenerate
boundary and says to compare several starts. That is why "no frmtmb
defect" stands: the default stop is a local optimum of the right
function, reported as such.

### p021: frmtmb's optimum is outside brms's parameter space

`dev/translator-gev-check.R` (log `dev/translator-gev-check.log`, seed 21
via `eq_mix(21)`): at frm()'s estimate, logLik -295.7098703163 equals a
hand-written mixture of two GEV densities, zero off each support, to
2.84e-13. Component 1 has xi1 = -0.266391 and puts 114 of 200 rows off
its own support; brms's `scale_xi()` holds xi1 in (-0.114796, 0.436034),
the interval that keeps every row inside component 1's support. That is
why brms's best of 6 `optimizing()` starts is 2.12 nats below frmtmb's
maximum (record, `l4_z` -298, -296, -2.12). On brms's parameter space
the two densities are an L1 identity at three points. So frmtmb fits the
mixture's likelihood and brms fits it on a restricted space: brms is
the side that differs, filed as brms-51 in `dev/upstream-bugs.md`. The
measured verdict "equivalent" holds for the model; the fits differ, by
brms's restriction. Per the fix-round instruction (m2), the record
keeps an adjudication row "equivalent" with brms-51 as the reason and the
2.12-nat difference in its text, so the record still shows that the
maxima differ: the exact row is in `dev/translator-adjudication-edits.tsv`
(section 8). The log was regenerated in the fix round from the current
script; the numbers are unchanged.

## 5. p044 and p069: what settles each

p044 (REML), built: `eq_laplace_stan()` integrates brms's own program
over Intercept, b, bs and zs at each of the three outer points (a central
difference Hessian of Stan's gradient, exact on a quadratic, then one
Newton step in Stan's chart) and compares it with minus frm()'s REML
objective there, less log|det| of the bs basis map. For a gaussian
response this is the REML integral, not an approximation of it, and it
holds to 4.6e-13 at all three points. The mutant (section 3) shows it
fails, by up to 1.02 nats, when the objective is not REML's. Cost: under
a second per point.

The log-determinant term (fix round, m3). In every brms spelling I
tried, brms's unpenalized smooth columns ARE frmtmb's: both packages take
them from mgcv (`smoothCon()` and `smooth2random()`), so the map is the
identity and `eq_bs_logdet()` is 0. `dev/translator-bslogdet.R` (log
`dev/translator-bslogdet.log`, p044's data, seed 44) fits each candidate
under REML and reads the marginal identity three ways: as built, with
the term dropped, and with its sign flipped.

| candidate | log abs det A | max abs(Xs - fx) | marginal as built | dropped | flipped |
|---|---|---|---|---|---|
| s(x) | 0 | 0 | 1.99e-13, 1.42e-13, -4.55e-13 | same | same |
| s(x, m = 3), Ks = 2 | 0 | 0 | -3.75e-12, 3.01e-11, -1.04e-11 | same | same |
| t2(x, z), Ks = 3 | 0 | 0 | 2.84e-14, 8.53e-14, -5.68e-14 | same | same |
| s(x, by = f), Ks = 3 | 0 | 0 | -1.08e-12, -3.92e-12, -7.96e-13 | same | same |
| s(x) + s(z), Ks = 2 | 0 | 0 | 1.14e-13, 5.12e-13, -1.71e-13 | not read (1) | not read (1) |
| s(xs), xs = 5000 + 1000 x | 0 | 0 | -4.83e-13, -1.14e-13, 7.96e-13 | same | same |
| s(x, bs = "cr") | 0 | 0 | -2.84e-14, 2.84e-14, -2.84e-14 | same | same |
| s(x, m = 3), brms's Xs times 3 | -2.19722 | 1.87 | -1.81e-12, 2.34e-13, 4.7e-12 | -2.2, -2.2, -2.2 (FAIL) | -4.39, -4.39, -4.39 (FAIL) |

The last row is constructed: brms's own program read on its own data
with Xs multiplied by 3 (the `l1_sdat_fix` hook), which is a change of
variables of bs. The joint identity is unchanged (-5.68e-14, -2.84e-14,
-8.53e-14); the integral gains log|det A| = -2 log 3 = -2.19722, and the
identity holds only with the term in, with the sign `eq_bs_logdet()`
gives it. Dropping it misses by exactly that amount at all three points,
flipping it by twice that. So the term is right, and no survey spelling
needs it.

(1) The dropped and flipped arms of s(x) + s(z), late in one R process
that had loaded many programs, stopped in rstan with "Failed
to initialize module pointer: NULL value passed for DllInfo" before any
L1 reading; the arm as built read normally, and its log-determinant is
0, so dropping or flipping it cannot change a number. A first run of the
script stopped the same way on the dropped arm of s(x, by = f), which
the final run read normally. This is rstan holding many compiled modules
in one process, not the identity.

p069 (importance), not built. The model is settled: L1 on the same spec
without `importance =` is an identity (-1.9e-10, unchanged). What is left
is whether the importance-corrected objective computes frmtmb's own
marginal, which is a frmtmb-internal question: brms has no integrated
likelihood to compare with, only a posterior, and the L4 posterior route
has no power on this design (review 2026-11-02, impsamp). The route that
would settle it is an exact marginal at three outer points (per-group
quadrature of the joint over the 2-dimensional (1 + x | g) effects, 10
groups) against the corrected objective, within its Monte Carlo error.
The review did exactly this once (`p069exact.txt`: 200, 512 and 2000
draws within 0.37 MCSE on all nine fits). Building it into the survey
would repeat that measurement with a Monte Carlo tolerance in a survey
whose other L1 tolerances are exact; it belongs in the importance test
file rather than in the survey, and it would not change the verdict.
I recommend not building it, and keeping p069 weak with the review's
reason, as the adjudication row already says.

## 6. No existing result moves

Rerun with the final instrument, after the second fix round: 76 pairs
across all four strata (`dev/translator-noreg3-out`, comparison
`dev/translator-noreg3-cmp.txt`; earlier runs are
`dev/translator-noreg-out`, 57 pairs, and `dev/translator-noreg2-out`,
67). They are the 9 pairs of section 4, the 15 existing pairs now read
at fixed points (p053, p056, p060, p070, p076, p087, p102, p106, p111,
p114, p115, p116, p122, p139; p069 is among the 9), every nl pair of the
sample (p001, p014, p031, p032, p049, p119, p127, p158, p166, p171),
every other smooth pair (p043, p062, p063, p067, p090, p098, p105, p110,
p148), every gp and hsgp pair (p022, p023, p094, p109, p026, p082, p092,
p117), every mixture pair (p005, p033, p045, p055), the cs() pairs
(p013, p015, p016, p146, p161), importance (p027, p124), rescor (p002,
p030, p041, p068, p151, p153) and controls (p008, p009, p012, p118,
p145, p150, p152, p154, p170). Every log's `lib:` line reads
rellib-r19 (76 of 76).

From the comparison's generated lines:

    same verdict and level as the record: 64 of 76
    l1_const: in both 58, bitwise equal 44, max |new - record| 4.146159e-09
    l1_grad: in both 58, bitwise equal 43, max |new - record| 8.378571e-09
    l2_metric: in both 11, bitwise equal 11, max |new - record| 0

The 12 that "differ": the 7 lane pairs whose level or verdict moved
(section 4: p021, p042, p047, p048, p131, p136, p143); p069, whose
measured verdict and level (defect, L1+L4) are the record's, and which
the script counts as different only because the record shows its
adjudicated level; and 4 described below.

The l1_const and l1_grad numbers that are not bitwise equal are the 15
pairs now read at fixed points, by design (one, p139, keeps its constant
bitwise by chance). Every one keeps its verdict and level:

| id | verdict, level (record and now) | record constant at the estimate | constants at fixed points | gradient now |
|---|---|---|---|---|
| p053 | equivalent, L1+L2 | -2.51e-11 | -1.14e-13, 5.68e-14, 0 | 4.66e-15 |
| p056 | equivalent, L1+L2 | 2.8e-11 | -1.14e-13, -1.14e-13, 1.14e-13 | 3e-15 |
| p060 | equivalent, L1 | -8.41e-11 | 2.84e-13, 2.27e-13, -5.68e-14 | 6.83e-10 |
| p069 | defect, L1+L4 (measured; adjudicated equivalent) | -1.88e-10 | 0, -1.14e-13, 1.71e-13 | 6e-15 |
| p070 | equivalent, L1 | -4.15e-9 | -2.84e-13, 1.71e-13, -5.68e-14 | 5.65e-10 |
| p076 | equivalent, L1 | 0 | -5.68e-14, -5.68e-14, 0 | 7.08e-14 |
| p087 | equivalent, L1 | 5.68e-14 | 1.14e-13, 0, 0 | 2.26e-14 |
| p102 | equivalent, L1+L2 | 1.37e-11 | -1.14e-13, -1.14e-13, -1.71e-13 | 3.01e-14 |
| p106 | equivalent, L1 | 4.31e-11 | 0, -1.71e-13, 1.71e-13 | 1.69e-14 |
| p111 | equivalent, L1 | -6.82e-11 | 0, -5.68e-14, -5.68e-14 | 2.51e-11 |
| p114 | equivalent, L1 | -1.23e-11 | -1.14e-13, 5.68e-14, -1.14e-13 | 1.75e-14 |
| p115 | equivalent, L1 | -1.15e-10 | 0, 2.27e-13, 0 | 9.21e-15 |
| p116 | equivalent, L1 | -5.16e-11 | -5.68e-14, -1.14e-13, 2.27e-13 | 5.61e-15 |
| p122 | equivalent, L1 | -7.43e-11 | -5.68e-14, -5.68e-14, 5.68e-14 | 3.17e-13 |
| p139 | equivalent, L1 | -1.14e-12 | -1.14e-12, -7.96e-13, -1.25e-12 | 8.38e-9 |

The largest move against the record, 4.1e-9 (p070), is the record's
constant at a degenerate estimate against a clean one at the fixed
points; the tolerance there is about 3e-4. The `l1_sdat_fix` hook is
inert unless a spec sets it (p044 and p048 rerun with and without it are
identical in every column but the timing, `dev/translator-hookcheck-out`).

The 4 exact gp() pairs differ by design: measured
"cannot measure" with no scored level, adjudicated "divergence by
decision", they now MEASURE "divergence by decision" at L1. Patched
constants and gradient (three points, then the gradient): p022 -2.54e-9,
-2.93e-9, -2.86e-9, 4.7e-5; p023 0.671, 0.671, 0.671 (the free constant
of its sdgp prior), 1.41e-5; p094 -2.33e-11, -7.99e-12, -3.2e-11,
3.65e-6; p109 -7.49e-10, -1.04e-9, -1.33e-9, 1.9e-5. Their adjudication
rows become redundant.

## 7. The rates, projected (not the record)

`dev/translator-project.R` builds a scratch copy of `dev/equiv-data/`,
writes each pair's row from the record with this lane's measurements
laid over it, and runs `dev/equiv-summary.R` there without `--record`.
Control: with no measurements and no adjudication edits it gives the
record's figures exactly (89.0 [84.5, 93.6]; weak evidence out 85.3; the
6 weak pairs p042, p044, p047, p069, p136, p143; L1 71 unique pairs).

The inputs are saved in the worktree, so this section regenerates:

    Rscript dev/translator-project.R dev/translator-proj-out \
      dev/translator-adjudication-edits.tsv dev/translator-noreg3-out

(log `dev/translator-proj.log`, summary
`dev/translator-proj-out/dev/proj-out/summary.md`). The control is the
same command with `dev/translator-adjudication-none.tsv` (a header only)
and no outdir (`dev/translator-projctl.log`).

With the 76 measurements of section 6 and the edits file of section 8
(p021's row replaced by "equivalent" with brms-51 as its reason; p048,
p131, p022, p023, p094, p109 removed), 4 verdicts adjudicated by hand
(p002, p021, p069, p072):

    agree with brms (the headline): 89.7, Wald [85.2, 94.2]
    no defect (agree, or a divergence by decision): 93.0, Wald [89.2, 96.8]
    sensitivity, weak evidence out: 89.0, Wald [84.4, 93.6]; the 1 pair
      (2 draws) left weak: p069
    L1+L2: 49 (57); L1: 82 (98); L4: 6 (9)

The fix-round note expected about 89.4. That was the first round's
reading with p021's old row kept, which counts p021 as "cannot measure",
not agree. Rerun that way (`dev/translator-adjudication-edits-oldp021.tsv`,
log `dev/translator-projold.log`): headline 89.4 [84.8, 93.9], no defect
92.7, weak evidence out 88.7. An "equivalent" row with agree TRUE counts
p021 as agree, so the headline is 89.7. The headline moves little
because p136 leaves "agree" for "divergence by decision" while p048,
p131 and p021 join it; the weak-evidence row moves from 85.3 to 89.0.

## 8. What needs the consolidation or the user

- The record is rerun once at consolidation. For it to show these
  results, `dev/equiv-data/adjudication.tsv` needs the edits in
  `dev/translator-adjudication-edits.tsv` (same columns; a row whose
  verdict is DROP is removed, any other row replaces the id's row).
  That file is under `dev/`'s root-level `*.tsv` ignore rule, so it is
  local only; its rows are copied here verbatim. Remove the rows of p048
  and p131 (the "cannot measure" reasons are gone: both are L1
  identities) and of p022, p023, p094 and p109 (the flow now measures
  their verdict at L1). Replace p021's row with (tab-separated, columns
  id, verdict, agree, level, proposal, reason; proposal empty):

      p021	equivalent	TRUE	L1		gen_extreme_value mixture: an L1 identity on brms's parameter space, read at three fixed points (constants -5.7e-14, 5.7e-14, 1.7e-13). The maxima differ by 2.12 nats, and brms is the side that differs (brms-51, dev/upstream-bugs.md): brms holds each component's xi inside the interval that keeps EVERY response inside that component's support, which a mixture does not need. frm()'s maximum -295.7099 equals a hand-written mixture density to 2.8e-13; its xi1 = -0.266 puts 114 of 200 rows off component 1's support, outside brms's interval (-0.115, 0.436); brms's best of 6 optimizing starts is -297.8329, 2.1230 nats below (dev/translator-findings.md, dev/translator-gev-check.R).

  This lane did not edit the record.
- p136 moves from agree (L4, weak) to "divergence by decision" (L1),
  consistent with p022, p023, p094 and p109 on the same nugget decision.
  The brief asked for "as far as the nugget decision allows"; this is
  where it allows.
- p069 stays weak (section 5); no change proposed.

## 9. What this lane did not do

- No package code, test or record changed; no frmtmb defect was found,
  so nothing was reported for a package change.
- No R CMD check and no OpenBLAS pass: no package or numerical package
  code changed.
- Gaps still open in the translator: an mv mixture's `theta_<resp>`
  (p005), mo() outside mu, `simo_<dpar>` (p030), and an unpredicted
  `nu` in a multivariate student model (p068). These pairs agree at L4
  by the optimizing route, which compares maxima and is not weak; they
  were not in this lane's list. Also `tmp_xi` with a predicted xi: brms's
  scalar `tmp_xi` cannot hold a per-row xi, and the rule stops with a
  message when xi varies by row. No sampled pair has one.
- An ordinal mixture with shared thresholds (`order = "mu"`) has no
  rule; the rule stops with a message rather than guess.
- p131's and p021's L1 has no gradient check (no point is an optimum);
  the value at three generic points, with all three required (m4) and
  the mutants of section 3 failing by 0.1 to 7 nats, is the evidence.
- The audit (section 3b) tests the blocks check C reads its gradient on.
  It does not test a block a pair has no latent vector for; such blocks
  have their outer parameters (sd, L) read directly by the value.

