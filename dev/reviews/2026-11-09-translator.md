# Adversarial review of lane translator, 2026-11-09

Base e0cce906 (frmtmb 0.80.0). Under review: the uncommitted `dev/` diff in
`frmtmb-wt-translator` (`dev/equiv-lib.R`, `equiv-flow.R`, `equiv-specs.R`,
`equiv-summary.R`, `equiv-findings.md`, `upstream-bugs.md`) and the new
`dev/translator-*` files. No package code changed. Everything here ran on
`rellib-r19` (every log's `lib:` line reads
`C:/Users/adf44/source/r/rellib-r19/frmtmb 0.80.0`). My scripts are
`dev/review-translator-*.R`, each with its `.log` beside it; remeasured
pairs are in `dev/review-translator-out/`. I formed the view below from the
code and the diff before I read `dev/translator-findings.md`.

## Verdict: FIX

The new rules are faithful changes of variables, and no rule can absorb a
model difference: each model mutant I put into brms's program fails L1
where the rule's pair can see the block it changes, and the gp decision
path cannot score a second difference "by decision". The REML marginal
identity, the fixed points, the upstream claim, the no-regression claim
and the projection all reproduce. One claim is overstated: on p048 and
p047 the survey's three points cannot see the t2() penalty blocks whose SD
collapsed to about 1e-4, so the L1 there certified one block of three (or
two of three), not the term. The lane's own fixed points close it, and
with them both pairs still hold (M1). The rest are record and evidence
items.

| | finding | what to do |
|---|---|---|
| M1 | p048 (2 of 3 t2() blocks) and p047 (1 of 3) are blind at the survey's points: a block reversed, or scaled in brms's model, still reads L1 TRUE | read those pairs at fixed points too (`l1_points` with `scale = 0.3`, joint); I measured that it holds and catches the mutants |
| m1 | p131: the L1 identity turns the record's 0.81-nat L4 gap into a statement about local optima; the findings do not say so | add the multistart numbers below |
| m2 | p021 measured "equivalent" carries no trace that the maxima differ by 2.12 nats | keep an adjudication row with brms-49 as the reason |
| m3 | `eq_bs_logdet()` is 0 on every pair, so its sign and value are untested | state it, or find a pair where it is not 0 |
| m4 | a fixed-point pair with no gradient check accepts a point that errors (`pts_ok` treats NA as pass) | require all three points when `l1_points` is set |
| m5 | evidence hygiene: the projection's adjudication edits are not saved; `translator-gev-check.log` is from another version of its script | save the edits file; regenerate the log |

## What I checked and how

### 1. The rules against model mutants

The lane's mutants (`dev/translator-mutants.R`) break the MAP. Mine
(`dev/review-translator-modmut.R`, `dev/review-translator-p048.R`,
`dev/review-translator-branches.R`) break the MODEL: brms's program is
edited by one fixed-string substitution, through `eq_l1()`'s own
`code_fix` hook, and the translator is left alone. Constants are at the
estimate and the two stepped points (or the three fixed points).

| model mutant (brms's program) | pair | rule under test | L1 ok | constants |
|---|---|---|---|---|
| smooth's penalized part x 1.01 | p042 | zs_prefixed | FALSE | 0.0143, 0.0156, 0.0181 |
| unpenalized columns dropped | p042 | zs_prefixed | FALSE | -16.8, -15.6, -18.3 |
| t2() block 1 x 1.01 | p048 | zs_prefixed | FALSE | 0.0464, 0.0469, 0.0486 |
| t2() block 2 x 1.01 | p048 | zs_prefixed | **TRUE** | 3.3e-10, 5.16e-10, 5.13e-10 |
| t2() block 3 x 1.01 | p048 | zs_prefixed | **TRUE** | 1.21e-10, 3.8e-10, 2.24e-10 |
| cs() sign of component 2 | p131 | bcs_prefixed | FALSE | -0.395, -1.96, -2.36 |
| disc2 = 1.1 | p131 | ord_mix_intercept, fixed points | FALSE | 1.97, 1.05, 0.18 |
| GEV `t = 1 + 1.02 xi x` | p021 | tmp_xi, fixed points | FALSE | -3.02, -3.33, -3.57 |
| gp length scale x 1.01, both programs | p136 | decision path | FALSE; patched FALSE | -0.0756, -0.0867, -0.0641; patched -0.0746, -0.0859, -0.063, gradient 2.22; verdict defect |
| brms nugget 1e-10 (patched: 1e-10 sdgp^2) | p136 | decision path | FALSE; patched FALSE | -0.0671, -0.199, -0.0147; patched -0.0415, -0.128, -0.0117, gradient 0.326; verdict defect |

The two TRUE rows are M1.

The rules:

- `zs_prefixed` reads `Zs_<prefix>_<i>_<j>` and keeps the suite's
  `brms_basis_map()` with its span and orthogonality checks, so a basis or
  penalty difference cannot be mapped away. It is caught wherever the
  block carries weight (p042: 0.0143 for a 1.01 scale, -16.8 for the
  unpenalized columns dropped; p048 block 1: 0.046). It is not caught
  where the block's SD is near 0 (M1).
- `ord_mix_intercept`: p131 has only cs() terms, so `X_mu<k>` is empty and
  the centering shift is never exercised by the lane's pair. With
  `mu<k> ~ x` the rule holds (constants 0, 0, 5.68e-14) and the shift with
  its sign flipped fails (-0.253, -1.59, -1.86;
  `dev/review-translator-branches.log`). The sign is right.
- `bcs_prefixed`: p131 has one cs() column, where `byrow` cannot matter.
  With `cs(x) + cs(z)` the rule holds (0, 0, 0) and the matrix read by
  column fails (-0.917, -0.459, 1.01). The layout is right.
- `tmp_xi` is brms's `scale_xi()` inverted at frmtmb's mu and sigma; a
  wrong mu or sigma moves brms's xi and shows. It cannot see the one
  model difference that lives in its domain, brms's restriction itself,
  which the lane names (m2).
- `eq_strip_pars()` changes only the name list the suite's translator
  reads; the compiled program is brms's unstripped one.
- `brms_inner_pat`: the `stopifnot()` keeps `sds_*` and `bs_*` out.

### 2. Fixed-point L1 for the mixtures (p131, p021)

The points are frm()'s own start moved by `scale * sin(1.7 j + k)`, k = 0,
1, 2: deterministic and independent of brms and of either optimum. The
constant must be 0 to 1e-6 * |ours| at all three, so a model difference is
caught unless it vanishes at all three generic points; my mutants of the
cs() sign (-0.395, -1.96, -2.36), the discrimination (1.97, 1.05, 0.18)
and the GEV shape (-3.02, -3.33, -3.57) all fail. That is enough for an analytic
difference. It cannot see a difference confined to a region the points do
not visit, which is exactly brms-49 (m2).

p131 (m1): `dev/review-translator-p131ms.R` runs frm() with an optimizer
that keeps the best of 21 seeded starts (seeds 1000 to 1020 around
frm()'s start, sd 1). frm()'s default stops at -267.932323 (singular
convergence, warned). Its own objective reaches -263.413103 from another
start, 4.52 nats higher and 3.71 above brms's best of 6 starts
(-267.118266, `dev/translator-before-out/p131.rds`). Two starts end at
-267.1183, brms's value to four decimals: an independent confirmation of
the L1 identity. So the record's 0.814-nat gap was two optimizers at
different local, degenerate points of one function (thresholds of
thousands in both). Not a frmtmb defect (frm() warns), but the findings
should say it, since "no frmtmb defect was found" otherwise reads as if
the L4 gap had vanished.

### 3. The gp decision path

`eq_gp_nugget_fix()` changes only `cov[n, n] += 1e-12;` (brms 2.23.0's
line, `dev/review-translator-out/code-p136.stan`; `gp_exp_quad_cov()` is
already scaled by sdgp^2, so `1e-12 * square(sdgp)` is frmtmb's relative
nugget). The decision path is entered only when the plain L1 fails and
scores "divergence by decision" only when L1 holds on the patched
program, which differs from brms's in that one line, and every other
level holds. A second difference must therefore show. My mutant with
brms's length scale at 1.01 times in both programs: the patched L1
fails (-0.0746, -0.0859, -0.063, gradient 2.22) and the verdict is
"defect", not "divergence by decision". A brms nugget of 1e-10 (the
patch then gives 1e-10 * sdgp^2, 100 times frmtmb's): the patched L1
fails too (-0.0415, -0.128, -0.0117, gradient 0.326), so L1 resolves the
nugget's size and the decision cannot cover a different nugget either.

### 4. brms-49 (p021)

Confirmed against brms 2.23.0's program (`make_stancode`,
`dev/review-translator-out/code-p021.stan`, lines 147-148):
`xi1 = scale_xi(tmp_xi1, Y, mu1, sigma1); xi2 = scale_xi(tmp_xi2, Y, mu2,
sigma2);`, and `scale_xi()` (`inst/chunks/fun_scale_xi.stan`) bounds xi
over ALL of `Y`. brms's own `gen_extreme_value_lpdf()` takes `log(t)` with
no support guard, which is why it needs the bound. The numbers reproduce
(`dev/review-translator-gevcheck-rerun.log`): logLik -295.7098703163
equals the hand mixture to 2.84e-13; xi1 -0.266391 against brms's
interval (-0.114796, 0.436034), 114 of 200 rows off component 1's
support; brms's best of 6 -297.8328965241, 2.1230262078 below. xi and
sigma are single-valued over the rows at the estimate
(`dev/review-translator-gevxi.log`).

### 5. REML and the no-regression claim

`eq_laplace_stan()` is exact for the gaussian quadratic; the ML mutant of
the lane fails by up to 1.02 nats while the joint holds. On a 2-D smooth
(`y ~ s(x, z)`, Ks = 2) the identity also holds (-1.14e-13, 2.56e-13,
-2.84e-14), but `eq_bs_logdet()` is again 0 (frmtmb uses mgcv's
unpenalized columns as brms does), and dropping it changes nothing (m3).
The derivation (brms's integral is |det A| times frmtmb's) is right on
paper.

The `equiv-summary.R` change is what it says: weak = agree and (L4 by
posterior, or REML without `l1_marg_ok`).

I remeasured 10 pairs with `dev/equiv-measure.R` through
`dev/translator-measure.sh` (`dev/review-translator-out`,
`dev/review-translator-cmp.txt`): p021, p022, p044, p048, p062, p109,
p119, p131, p136, p143. Every column the comparison reads (verdict,
level, l1_const, l1_grad, l1_extra, l1_marg, l1_dec, l1_points) is
identical to the lane's outdirs, and l1_const and l1_grad are bitwise
equal to the record on the 6 pairs that have them.

### 6. The projection

Rebuilt with `dev/translator-project.R` and my own edits file
(`dev/review-translator-edits.tsv`: DROP p021, p048, p131, p022, p023,
p094, p109): headline 89.7 [85.2, 94.2], no defect 93.0 [89.2, 96.8],
weak evidence out 89.0 [84.4, 93.6], weak pairs p069 only; with p021's
row kept (`dev/review-translator-edits-keep021.tsv`): 89.4, 92.7, 88.7
(`dev/review-translator-proj.log`, `-proj-keep021.log`). These match the
lane's section 7.

## Findings

### M1 (MAJOR). Two of the lifted t2() pairs are blind to their collapsed blocks

At frm()'s estimate the t2() SDs are (5.04, 4.21e-4, 2.56e-4) for p048
and (3.67, 2.55, 1.36e-4) for p047 (`dev/review-translator-sds.log`). At
the survey's points (the estimate and +/- 0.05) a collapsed block adds
about 1e-7 to the predictor (`|s_j|` 7.7e-8, 2.8e-8 for p048,
`dev/review-translator-p048.log`), far below the tolerance 1e-6 * |ours|.
So:

- the map of one block reversed reads L1 TRUE: p048 `zs_y1_1_2` (-6.23e-8,
  -9.4e-8, -7.88e-8), `zs_y1_1_3` (-1.92e-8, -5.54e-8, -3.39e-8); p047
  `zs_b1_1_3` (-7.87e-9, -1.62e-7, -5.22e-8,
  `dev/review-translator-p047survey.log`);
- the MODEL with p048's block 2 scaled by 1.01 in brms's program reads
  L1 TRUE (3.3e-10, 5.16e-10, 5.13e-10), and with block 3 scaled
  (1.21e-10, 3.8e-10, 2.24e-10, `dev/review-translator-p048.log`).

The lane's p048 mutant reversed all three blocks, and block 1 alone made
it fail (-34.6, the same number as block 1 reversed alone). So the claim
"p048 and p047 reach L1" holds for block 1 (and 2) and says nothing about
the rest of the term. This is not the translator made to agree; it is a
reading point with no power on those blocks.

Remedy, measured (`dev/review-translator-p048fix.R`): the same pairs with
`l1_points = list(scale = 0.3)` (frm()'s start, where every SD is near
1, joint, so the inner gradient is still read): p048 control
-1.14e-13, -3.41e-13, -1.14e-13, gradient 7.2e-16; `zs_y1_1_2` reversed
fails (-3.31, -1.64, -2.65), `zs_y1_1_3` reversed fails (-2.13, -0.383,
-1.22); p047 control -1.14e-13, 2.27e-13, 1.14e-13; `zs_b1_1_3` reversed
fails (-215, -178, -119). So both identities are true once read where
the blocks count. Add the fixed points to these two specs (or, better,
read every smooth or random-effect pair at one extra point with each SD
away from 0, since the same blindness applies to any pre-existing L1 pair
whose SD collapses), and add the single-block mutants to
`dev/translator-mutants.R`.

### m1 (MINOR). p131's L4 gap is local optima; say so

Numbers in section 2 above (`dev/review-translator-p131ms.log`). The
lane's section 4 and 9 should carry them; "no frmtmb defect was found"
is right, and this is why.

### m2 (MINOR). p021's record row should keep the disagreement

With the adjudication row dropped, the flow writes "equivalent, every
level agrees" at L1 alone, and the record no longer shows that frm() and
brms reach maxima 2.12 nats apart. The lane leaves the ruling to the user
and recommends "equivalent" with brms-49 as the reason; that needs a kept
row, not a dropped one. Section 7's headline (89.7) assumes the row is
dropped.

### m3 (MINOR). `eq_bs_logdet()` has never been nonzero

0 on p044 and on `y ~ s(x, z)` under REML; a mutant that returns 0
changes no number. Its sign is right by the derivation, but nothing
measures it. Say so in the findings.

### m4 (MINOR). Fixed points accept an erroring point

`pts_ok <- all(is.na(pts) | ...)`: a point whose rule stops (for p021, an
xi outside brms's interval at a larger `scale`) drops silently to fewer
points. At the optimum the gradient backs the value; with `l1_points` on
a model without random effects nothing does. Require `l1_points == 3`
when `spec$l1_points` is set. Both p131 and p021 have 3 today.

### m5 (MINOR). Evidence hygiene

- The projection's adjudication edits file is not in the worktree, so
  section 7 cannot be regenerated from what is saved. Mine reproduce it.
- `dev/translator-gev-check.log` contains `links:` output and three
  recycling warnings that the current script does not produce; my rerun of
  the same script is clean and gives the same numbers. Regenerate the log.
- The `tmp_xi` rule stops when xi is not a single value; a predicted xi
  has no rule. State it with the other open gaps.

## Not run by me

R CMD check and OpenBLAS (no package code changed). The lane's 57-pair
no-regression set beyond my 10.

FIX: M1 (MAJOR); m1, m2, m3, m4, m5 (MINOR). No BLOCKER.

## Fix round, short check (2026-11-10)

I read the diff since my review (`eq_sd_collapsed()` and four specs in
`dev/equiv-specs.R`, the `l1_sdat_fix` hook in `eq_l1()`, the stricter
`pts_ok` in `eq_measure()`), the new `dev/translator-audit*` files,
`dev/translator-bslogdet.R` and its log, and the edits files. New scratch:
`dev/review-translator-audit2.R`, `dev/review-translator-p087.R`, each with
its `.log`.

### What holds

- **Coverage of the audit.** `dev/translator-audit-ids.txt` holds 131
  ids: all 120 pairs whose record level contains L1, plus the 11 this lane
  lifted or touched. p069 makes the 132nd log. No L1 pair is missing
  (`dev/review-translator-audit2.log`).
- **p076 and p139.** Both are still equivalent at L1 on fixed points
  (`dev/translator-noreg2-cmp.txt`: p076 -5.68e-14, -5.68e-14, 0, gradient
  7.1e-14; p139 -1.14e-12, -7.96e-13, -1.25e-12, gradient 8.4e-9). On
  those points each block is seen (`dev/translator-audit-fixed-out`:
  p076 4.42 nats, p139 59.6 and 111, p047 26 to 306, p048 2.13 to 3.89).
- **The reversal test.** It is sound for a block with one SD. Reversal
  keeps |z| and so the N(0, 1) prior exactly, so only the likelihood
  can see it. The "seen by less than 100 x tol: 0 blocks" line answers
  the near-miss question for whole blocks. The next weakest block is
  p106 z_1 at 9.4e3 x tol.
- **m2.** I reran `dev/translator-project.R` with
  `dev/translator-adjudication-edits.tsv` on `dev/translator-noreg2-out`
  and got 67 pairs and 4 adjudications: headline 89.7 [85.2, 94.2], no
  defect 93.0, weak evidence out 89.0 (`dev/review-translator-proj2.log`).
  The 89.4 was the old "cannot measure" row. An adjudicated "equivalent"
  with agree TRUE counts the same as a measured one, which is how p072
  and p069 are already scored. The row keeps the 2.12-nat gap and
  brms-49 in its reason. That is consistent with the standard.
- **m3.** For `s(x, m = 3)` with brms's Xs times 3 the term is -2.19722,
  which is -2 log 3 for Ks = 2. As built the identity holds to 4.7e-12.
  Dropping the term misses by -2.2, and flipping its sign misses by
  -4.39. The sign is now measured. Scaling Xs is a valid change of
  variables, because bs has a flat prior. The "Failed to initialize
  module pointer" arms are only the `s(x) + s(z)` dropped and flipped
  arms, where the term is 0. They change nothing here, but they are an
  rstan fault in a re-sourced environment, so watch for it in the
  consolidation run.
- **m4.** The code reads `length(pts) == 2L && !anyNA(pts)` when
  `l1_points` is set. The lane's "third point unreadable" mutant now
  fails.
- **m5.** All three items are done:
  - The edits files are saved.
  - `dev/translator-gev-check.log` is clean: same numbers, no warnings,
    no `links:` lines.
  - The predicted-xi gap is listed in findings section 9.

### F1 (MAJOR). The audit misses a half-collapsed block: p087

The audit reverses a WHOLE block, so a matrix `z_<i>` with one live row
and one collapsed row reads as "seen". p087 (`gr_cov x mmc()`, `(1 +
mmc(x1, x2) | mm(g1, g2, cov = A))`) has z_1 SDs 0.463 / 1.12e-06. The
audit calls it seen at 28.3 nats, 8.5e4 x tol
(`dev/translator-audit-out/p087.tsv`). I reversed only row 2 (the mmc()
slope) after the translator (`dev/review-translator-p087.log`):

| points | mutant | L1 ok | constants |
|---|---|---|---|
| survey (estimate, +/- 0.05) | none | TRUE | 5.68e-14, -3.41e-13, -1.14e-13 |
| survey | row 2 reversed | **TRUE** | -1.3e-11, -1.06e-11, -1.68e-11 |
| survey | row 1 reversed | FALSE | -25, -22.2, -28.3 |
| fixed, scale 0.3 | none | TRUE | 1.14e-13, 0, 0 |
| fixed, scale 0.3 | row 2 reversed | FALSE | -11.2, -7.47, -2.16 |

So "re-audits 0 blind" is not true. p087's mmc() slope, and with it the
correlation in `L_1`, is not tested by its L1 at the survey's points.
The pair is still equivalent at fixed points. p087 is the only pair in
the audit with an SD below 1e-2 next to a large one in a block the audit
calls seen (`dev/review-translator-audit2.log`). The audit cannot tell
whether a `gr(by = )` block with one collapsed by-group has the same
problem, because it reverses across levels.

Fix: in `dev/translator-audit.R`, reverse each row (term) of a matrix
block separately, and for a `by` block each by-group's levels. Then put
p087 (and anything else the audit finds) on `eq_sd_collapsed()`, and
rerun the audit and the no-regression set.

FIX: 1. F1 (MAJOR), the audit misses a half-collapsed matrix block; p087's
mmc() slope row (SD 1.12e-6) is blind at the survey's points
(`dev/review-translator-p087.R`, data seed 87 via `eq_mm(87)`). M1's
blocks, m2, m3, m4 and m5 are resolved.

## Fix round 2, check of F1 (2026-11-10)

I read `eq_audit_mutations()` in `dev/translator-audit.R`, the
`eq_degenerate` list in `dev/equiv-specs.R`,
`dev/translator-audit3-summary.txt` and `dev/translator-noreg3-cmp.txt`.
New scratch: `dev/review-translator-audit3.R` and
`dev/review-translator-corrow.R`, each with its `.log`.

### The per-row reversal is prior-exact

Reversing a row of `z` (or a by-group's slice of a row) permutes iid
std_normal entries, so the prior does not change. The same holds after
the `cov = A` step, because brms's own `z` is iid and its program
applies Lcov to it. The audit also reverses Ymi, whose prior is not iid.
That can only make a block look more visible. It is moot here: 2 of the 3
Ymi blocks are untestable, and p147's Ymi_xm is seen at 1.7e4 x tol. No
audited pair has rcar or zcar.

### The audit covers every block shape in the sample

The final audit has 130 parts in 66 pairs:

| kind | parts |
|---|---|
| z | 93 |
| zs | 19 |
| zgp | 10 |
| zme | 4 |
| Ymi | 3 |
| zerr | 1 |

55 of the parts are rows. Every block shape in the sample is among them:
correlated blocks, mm() and mmc() blocks, multi-term groups, and the
shared `|ID|` blocks across dpars and responses (p003, row 2, seen at
5.1e3 x tol). Each t2() penalty is its own zs block.

### No L1 pair can hold a gr(by = ) block

No spec in `dev/equiv-specs.R` and no feature in `sample.tsv` uses
`by =` inside `gr()`. The audit's by-group branch fired on 0 parts. With
no translation rule for gr(by = ), such a pair could not reach L1 at all.

### There is no |cor| threshold, and none is needed

The trigger is the audit's outcome, not a cutoff on |cor|. `eq_degenerate`
is a hand list built from it. What answers "just above" is the margin:
the weakest part in the final audit is p048 `zs_y1_1_3` at 3.2e3 x tol,
and nothing lies under 100 x tol.

The audit judges blindness on the value alone, which errs toward calling
a part blind. In my spot checks, L1's gradient check did catch the row-2
reversal at the survey's points, but only by 2.9 to 3.6 times its 1e-3
threshold. Moving those pairs to fixed points is therefore right.

### Spot checks of two correlation pairs

Row 2 of `z_1` reversed after the translator
(`dev/review-translator-corrow.log`):

| pair | points | mutant | L1 | constants | gradient |
|---|---|---|---|---|---|
| p116 | survey | none | TRUE | -5.2e-11, -1.7e-11, 6.9e-11 | 3.8e-11 |
| p116 | survey | row 2 reversed | FALSE, by the gradient only | -2.6e-5, -2.7e-5, -2.5e-5 (under tol) | 3.57e-3 |
| p116 | fixed | none | TRUE | -5.7e-14, -1.1e-13, 2.3e-13 | 5.6e-15 |
| p116 | fixed | row 2 reversed | FALSE | -7.55, -5.6, -4.47 | 17.6 |
| p053 | survey | none | TRUE | -2.5e-11, 2.1e-10, 2.9e-11 | 5.6e-11 |
| p053 | survey | row 2 reversed | FALSE, by the gradient only | -1.7e-5, -1.5e-5, -1.9e-5 (under tol) | 2.87e-3 |
| p053 | fixed | none | TRUE | -1.1e-13, 5.7e-14, 0 | 4.7e-15 |
| p053 | fixed | row 2 reversed | FALSE | -5.84, -12.3, -19.6 | 4.49 |

p087 at the survey's points was blind on both channels (gradient 2.5e-6,
first check above). On fixed points every check fails by nats.

### Counts

The no-regression set has 76 pairs, and 64 keep their verdict and level.
The other 12 are the pairs this lane lifted or re-scored. 19 outdir rows
read fixed points: the 13 `eq_degenerate` pairs, p047, p048, p076, p139,
p131 and p021. So 15 existing L1 readings moved, as the lane claims, and
none changed verdict.

### Note

`eq_degenerate` holds this build's estimates. A later build that moves
an estimate onto a boundary needs the audit again. Run
`dev/translator-audit.sh` beside the record at consolidation. This is a
procedure note, not a finding.

SHIP. F1 is resolved: the per-row audit is prior-exact, covers every
block shape in the sample, reads 0 blind with the weakest part at
3.2e3 x tol, and my spot checks of p116 and p053 confirm both the
blindness and the fixed-point remedy.
