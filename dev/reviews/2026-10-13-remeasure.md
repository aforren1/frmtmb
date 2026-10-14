# Review of lane remeasure, 2026-10-13

Reviewer: adversarial review of `dev/remeasure-findings.md` and the
bin-2 port, base d2ef6312 (frmtmb 0.70.0), build `rellib-r8`.

## Round 3 (targeted, after fix round 2): MERGEABLE AFTER FIXES

This round checks only the coordinator's six items. Scripts and logs
are `dev/remeasure-rev3-*`. The stage directory is deleted.

### Remaining issues

1. **14 of the 62 bin-2P `brms-internal` decision rows fail the G1
   test.** The ledger listing is `dev/remeasure-rev3-led.txt`.
   - The 10 `inv_link_<ordinal>()` array rows (`distributions:300`,
     8; `:479`, 2) are the same probability map as `dcumulative()`.
     The lane's own `brms_adapt_dord()` answers them. On brms's block
     text with brms's reference helpers (`dev/remeasure-rev3-g2.R`,
     seed 300), every line holds in 32 of 32 iterations except the
     cumulative identity, which holds in 16 of 32 (the same domain
     limit as `distributions:238`). That is 9 passes and 1
     divergence.
   - `dcategorical()` (`:527`, 2) and `inv_link_categorical()`
     (`:549`, 2) are the categorical probability map. The lane files
     categorical `log_lik` and `posterior_predict` as "not ported
     (cost)", so these 4 rows are not ported, not a decision.
   - The other 48 hold up as decisions on the sample: the `link_*()`
     inverse maps and the "inverts" blocks (22), G5 `brmsterms()`
     (18), G6 (5) and G7 `insert_refcat()` (3). frmtmb has no
     probability-to-predictor map, no parsed-formula object, no
     parameter-to-matrix function and no `insert_refcat()`.
   - Effect on the rate: 534/709 = 75.3% with the 10 ordinal rows
     ported, or 525/709 = 74.0% left unported, against 525/695 =
     75.5%. The headline would be 534/823 = 64.9%.
2. **Wording of the categorical rows and "ported".** Per the ledger,
   the categorical rows `log_lik:495` and `posterior_predict:344`
   carry `ported = TRUE` but class `not-ported`. They were run through
   the existing adapter, and it refused them because it does not map
   `mu1` and `mu2`. "Not ported (cost)" is the true outcome; "tried
   with no new adapter code" is accurate but reads as a port. Likewise
   "168 ported through one thin adapter" counts every row run through
   the adapter, including 31 absent, 15 no-cdf, 3 no-quantile and 2
   categorical rows that it refused or that failed. Say "168 run
   through the adapter, 113 pass".

### Checked and sound

- **`brms-backends` (3 wiener rows).** The block compares brms's
  RWiener and rtdists backends with each other. A frmtmb shim would
  make `d1` and `d2` the same call, so the restated assertion would be
  hollow. The lane checked what can be checked: frmtmb.eam agrees with
  both backends to 5.68e-14 on the log scale, on 9 of 10 rows (one
  brms density underflows to 0). A decision is fair. Leaving these
  rows out also does not raise the rate, since frmtmb would plausibly
  pass the density row.
- **The G1 shims can fail.** I recorded the generated distributions
  test with the true adapter and with `brms_g1_d()`'s `row_lpdf()`
  multiplied by 1 + 1e-6 (`dev/remeasure-rev3-g1sayno.R`). 25 G1 rows
  hold on the true adapter. Under the perturbation the 2 value rows
  break (`distributions:4`, the `integrate()` row; `:5`, the
  `dt()/5` identity) and the 23 shape rows still hold. Only 2 of the
  24 d/r passes compare a value; the other 22 are `length()` checks
  as brms wrote them. The ledger classes them as shape, so the 24
  value passes are 16 log_lik + 2 `:446` + 1 2T + 2 G1 + 3 split
  `:238` lines. That is correct.
- **`$lcdf` returning the CDF is the convention, not a wrong answer.**
  `frmtmb_family()` documents `lcdf` as an "AD log-safe CDF ...
  returning probabilities" (`R/families.R:113`). Every family follows
  that convention: gaussian `pnorm`, poisson `ppois`, lognormal,
  exponential `1 - exp()`, weibull, inverse.gaussian `pinvgauss`, cox,
  and frmtmb.eam's rdm and wiener `exp()` of a log CDF. The only
  caller is `fam_lcdf()` in `R/objective.R` (lines 81, 84, 126 and
  142), which uses it as F: `log(Fv - Flb)`, `log(F2 - Fv)`,
  `log(Fub - Fv)`, `log(Fub - Flb)`. frmtmb.eam's contaminant mixture
  also mixes it as a probability. No caller treats it as a log CDF,
  and the lane's `pinv_gaussian()` shim goes through `row_lpdf()` of
  a censored row, which is correct. The slot name is misleading, and
  that is all.
- **Arithmetic.**
  - 389 + 23 + 113 = 525; 525/823 = 63.8%; 823 - 128 = 695;
    525/695 = 75.5%.
  - The 2P classes sum to 270 and 2T to 59.
  - The decisions total 128 = bin 1 (stan 26, brms-internal 15,
    mcmc 4, no-draws 3) + 2P (brms-internal 62, no-density 12,
    no-quantile 3, brms-backends 3).
  - The gaps total 65 = 2 + 17 + 31 + 15.

## Round 2 (after the lane's fix round 1): MERGEABLE AFTER FIXES

Seven of my eight round-1 items are resolved, and every regenerated
number reproduces. The one item that is not resolved is item 1, the
rate. The lane's new "no counterpart (decision)" bucket still holds
rows that are not decisions, and one whole-block verdict counts 7
assertions that hold as divergences. All the fixes are to the record
and the classification. None needs package code.

### Fixes required (round 2)

1. **52 of the 67 G1 rows (`no-export`) are not decisions.** The
   survey's reason reads "frmtmb exports no d/p/q/r function under
   brms's names, by decision". It cites no record. The nearest one is
   `SPEC.md:154`: numeric d/p/q/r come from gamlss.dist or extraDistr.
   That decision is about exported NAMES. The assertions do not need
   the names. The lane itself answers brms's UNEXPORTED
   `dcumulative()` through the ordinal family object
   (`brms_adapt_dord()`), and calls that a thin read. The same read
   answers the G1 rows: `dev/remeasure-rev2-g1.R` (seed 20261013, log
   `.txt`) shims brms's `d*`/`r*` over `frmtmb::row_lpdf()` and
   `frmtmb::sim_draw()`. 24 of 24 d/r assertions hold for families
   frmtmb has. These are student (integral, the `dt()/5` identity, r),
   von_mises, skew_normal, exgaussian, inverse.gaussian,
   beta_binomial and asym_laplace (d and r), plus 5 zero-inflated d
   and 4 hurdle d. On the 5 families I checked against brms's exported
   d functions (zero_inflated_poisson, hurdle_poisson, hurdle_gamma,
   skew_normal, von_mises), they agree to 1.7e-16 or better. The 67 rows split as follows:

   | G1 rows | assertions | honest class |
   |---|---|---|
   | d and r of families frmtmb has | 24 | transferable: holds through the family object (port it, or "not ported") |
   | p of families frmtmb has, no `$lcdf` (all but inverse.gaussian) | 16 | absent (gap) |
   | frechet, gen_extreme_value (family absent, famuni) | 9 | absent (gap) |
   | wiener (frmtmb.eam has the family) | 3 | not checked; probably transferable |
   | multivariate normal and t, incl. 8 validation messages | 12 | decision (no density outside a fit) |
   | q (no quantile functions) | 3 | decision, per SPEC.md |

   My round-1 review called G1 a judgment call and accepted the
   decision label (my 82.9% assumed it). I withdraw that. The lane's
   own G2 construction shows that the counterpart exists.
2. **`distributions:238` is one whole-block verdict over 8 assertion
   lines.** 7 of them hold. I ran brms's block (seed 235, the lane's
   `brms_adapt_dord()`, brms's reference helpers;
   `dev/remeasure-rev2-d238.R`, log `.txt`): every line except the
   cumulative identity holds in all 32 iterations, and the cumulative
   identity holds in 16 of 32. These are the 3 sratio, cratio and acat
   identities and all 4 `dim()` lines. Only cumulative needs ordered
   thresholds. The lane's own `bin2-probe3.txt` agrees: sratio, cratio
   and acat show 0 NaN at unsorted thresholds. Record 7 passes and 1
   divergence (domain), or state that 7 of the 8 hold. Also correct
   the reason text: frmtmb's NaN is not only "in the categories after
   it". At draw 5 of `bin2-probe3.txt`, categories 1 and 2 are NaN,
   where brms's category 1 is 0.805.
3. **`log_lik:302` (com_poisson) is filed as "brms-internal,
   decision".** After R6 it is a known parameterization divergence:
   frmtmb's `mu` is the mean, and the reference helper's `mu` is
   `lambda^(1/shape)`. Move it to divergence, in the denominator. It is
   1 row.
4. **The 40 unported G3/G4 rows are labeled "brms-internal,
   decision".** These are categorical, multinomial, dirichlet,
   ARMA/SAR/FCOR, multivariate, wiener and custom. The 46 ported rows
   of the same groups read the same brms internals, so "brms-internal"
   is not why these were left out. Cost is. Label them "not ported"
   in the denominator, or "absent" for the families frmtmb lacks
   (dirichlet, dirichlet_multinomial, logistic_normal). Keep
   "decision" for G5 to G7 (brms's parsed-formula objects,
   `insert_refcat()`) and for G2's link-array maps.
5. **Restate the rate.** As published it is 493 / 605 = 81.5%. The
   table shows the effect of each fix, one cell per assumption:

   | applied | rate |
   |---|---|
   | fix 2 only | 500 / 605 = 82.6% |
   | fixes 1 to 3, G1 rows unported | 500 / 658 = 76.0% |
   | the same, the 24 G1 d/r rows ported (they hold) | 524 / 658 = 79.6% |
   | fixes 1 to 4, G1 unported | 500 / 698 = 71.6% |
   | fixes 1 to 4, the 24 ported | 524 / 698 = 75.1% |

   With fix 2 the headline becomes 500 of 823 (60.8%). Whichever
   choice the lane makes, the record should say which rows each
   bucket holds and why.
6. Minor wording in R6: "brms's posterior mean of `exp(Intercept)` is
   7.30" should read `exp()` of the posterior mean of the Intercept
   (exp(1.988)). The number is right.

### Verified in round 2

- **Bin 2 reproduces.** I ran `brmsport-bin2-all.sh` in a copy root
  (`dev/remeasure-rev2-b2/`, since deleted;
  `dev/remeasure-rev2-b2-all.log`). The summary, ledger, verdicts,
  survey, groups, targets, adapted ids and the generated log_lik,
  posterior_predict and distributions tests are byte-identical. The
  tier is green on all 6 files (58, 2, 46, 58, 45, 1).
- **Arithmetic.**
  - 389 + 23 + 81 = 493; 493 / 823 = 59.9%; 823 - 218 = 605;
    493 / 605 = 81.5%.
  - Bin 2: 104 + 9 + 16 + 170 + 26 + 4 = 329.
  - 2P: 81 + 8 + 9 + 1 + 103 + 67 + 1 = 270.
  - Ported 101 = 81 passes + 20 that do not pass (9 absent, 1
    adapter, 1 reference, 1 brms-internal, 8 divergence).
  - Bin 1's 51 rows split 48, 2 and 1.
  - The phrasing "2P 81 of 270, 101 ported" is correct.
- **Ordinal rows.** All 12 ordinal rows (`log_lik:462` to `:478` and
  `posterior_predict:307` to `:327`) are "pass" in the ledger, and
  `distributions:446` (2 assertions) passes
  (`dev/remeasure-rev2-rows.R`, log `.txt`).
- **Say-no check.** I reran `dev/remeasure-bin2-sayno.R` in my copy
  root (`dev/remeasure-rev2-sayno.txt`). The result is identical to
  the lane's: 38 log_lik passes; 16 of 16 value passes break, and 22
  of 22 shape passes hold, at a relative 1e-6. `log_lik:32` is
  cannot transfer, class adapter. The pass split of 19 value and 85
  shape matches round 1's listing once the ordinal ports are added.
- **Vignette rerun.** `remeasure-run.sh` against rellib-r8 in a copy
  (`dev/remeasure-rev2-wt/`, deleted; `dev/remeasure-rev2-run-r8.log`,
  237 s at 6 jobs) gives byte-identical `families`, `history`,
  `plauscmp`, `plausibility` and `secmp`, and identical outcome
  tables. `gaps` is identical once the bin-2 ledger is in the copy.
  The defects table reads it. Refilling the findings changes only the
  wall-time block.
- **Families.** The instrument now fits each family: 34 of 40 in
  core, 35 with frmtmb.eam, and none named but unfitted. `multinomial`
  goes through `multinomial(K = 3)`, and the text says that brms's
  bare spelling fails (rank 26). Honest.
- **Median.** On the common set (34 models, 124 coefficients) it is
  0.047 on every build. This matches `dev/remeasure-rev-median.txt`
  (0.0470). The "whole change" claim is withdrawn.
- **R1.** The octopus-merge wording is correct, the mechanism is
  recorded, and the owner is optstart.
- **R6 (brms-24).** I confirmed that `brms:::log_lik_com_poisson()`
  and `posterior_predict_com_poisson()` call the R-side
  `dcom_poisson()` and `rdiscrete(dist = "com_poisson")`, and that
  the Stan-side mean at rate 3, nu 2 is 1.453549
  (`dev/remeasure-rev2-cp.R`, log `.txt`), as the lane's entry says.
- **Gap list.** Rank 7 is "deliberate: as brms", like ranks 13, 15, 19
  and 20. Rank 18 and the six 2T rows on several models are open work
  for postfix. R1 to R4 are ranked with owners.

Cleanup: I deleted my large trees (`dev/remeasure-rev2-b2/`,
`dev/remeasure-rev2-wt/`), the commit-build libraries in
`C:/Users/adf44/source/r/remeasure-rev-lib/` and the extracted commit
trees in the scratchpad. Scripts and logs stay as
`dev/remeasure-rev-*` and `dev/remeasure-rev2-*`.

## Round 1 verdict: MERGEABLE AFTER FIXES

The measurement machinery is sound and reproduces exactly. The fixes
are to the record: four claims are wrong or overstated, and the
bin-2 classification moves rows out of the denominator that belong in
it. No fix needs a package change. Fixes, numbered, each with its
evidence below:

1. The "87.1% of the 551 that can transfer" figure removes 95 rows the
   lane itself calls gaps. Restate it (section B1).
2. The 2P survey classifies 12 ordinal `log_lik`/`posterior_predict`
   assertions as "brms-internal, high cost". A thin read of the
   exported family object answers them. Port them or reclassify them
   as "not ported" (section B2).
3. Say what the bin-2 passes assert: 73 of 91 are shape, length,
   name or type checks. `log_lik:32` tests the adapter's own weight
   multiplication, not frmtmb (section B3).
4. R6 is wrong as stated. brms 2.23.0's FITTED `mu` is the
   COM-Poisson rate lambda, not `lambda^(1/shape)`. The R-side helper
   the lane read disagrees with brms's own Stan program (section C).
5. R1: the history is an octopus merge, so "bisected, 58b97a2c right,
   ea711734 wrong" is not a bisection. State the evidence that holds,
   and record the mechanism found here (section D).
6. The median `|z|` headline (0.0458 against 0.0488) compares two
   different coefficient sets. On the common set both builds give
   0.0470 (section E).
7. The families instrument checks names, not fits; "frm() accepts" is
   not what it measures. `multinomial` counts as covered while gap
   rank 26 says brms's spelling fails (section F).
8. Gap-list mis-assignments: rank 7 against ranks 13, 15, 19 and 20;
   R2 to R4 have no owner; six "policy" rows contradict
   `dev/brms-api-diff.md` (section G).

All reviewer scripts and logs are `dev/remeasure-rev-*` in this
worktree. Commit builds are in
`C:/Users/adf44/source/r/remeasure-rev-lib/c<sha>`. I did not edit a
lane file.

## A. Reproduction (holds)

- `dev/remeasure-run.sh` writes into its own tree (the outcome
  tables under `dev/brms-port/`), so I ran it in a copy,
  `dev/remeasure-rev-wt/dev/` (scripts, `brms-ref`, the v067, r6 and
  r7 records), against `rellib-r8`, 6 jobs. Wall time: 135 s
  (`dev/remeasure-rev-run-r8.log`). All six reducer outputs
  (`families`, `gaps`, `history`, `plauscmp`, `plausibility`, `secmp`)
  and the three outcome tables are byte-identical to the lane's.
  `remeasure-blocks.R v067 r8` followed by `remeasure-fill.R` on a
  copy of the findings gives a file that differs from the lane's only
  in the wall-time block.
- I recounted the 0.70.0 numbers from the raw per-vignette records,
  without the lane's reducers (`dev/remeasure-rev-counts.R`, log
  `.txt`). Model: raw OK 37, spell OK 40 of 42, CLEAN 37. Post: raw
  OK 65, spell OK 73 of 102, CLEAN 65. Hand: 458 of 539 non-data
  rows (28 data rows are excluded). These counts match the record.
- I reran 0.69.1 against `rellib-r7` in the same copy (239 s,
  `dev/remeasure-rev-run-r7.log`). `families`, `gaps` and
  `plausibility` are identical to the lane's r7 records, and both
  outcome tables are byte-identical after the r7 column is rewritten.
  `history.txt` differs only in the order of its build rows, because
  the r8 column was already present. Model 36/36/39, post 64/64, hand
  453: as recorded.
- Bin 2: I ran `dev/brmsport-bin2-all.sh` in a copy root
  (`dev/remeasure-rev-b2/`, which holds the package helpers, the
  `brmsport-bin2-*` files, the suite and the bin-1 ledger). The
  summary, ledger, verdicts, survey, groups, targets, adapted ids
  and the generated tests are byte-identical. The tier is green on
  all 5 files (58, 40, 58, 39, 1 passes; `dev/remeasure-rev-b2-all.log`).
- The instruments can say no. I removed the adapter's weights line:
  `log_lik` gave pass=39 fail=1 (`dev/remeasure-rev-noweights.txt`).
  I flipped the verdict of `log_lik:50` from pass to divergence: the
  result was fail=1 (`dev/remeasure-rev-flipverdict.txt`). Both
  copies are in `dev/remeasure-rev-b3/`.
- Arithmetic: 389 + 23 + 68 = 480; 823 - 272 = 551; 480 / 551 =
  87.11%. Bin 1 has 494 rows with 494 distinct ids. Bin 2 has 294
  rows whose `n` column sums to 329 (2T 59, 2P 270). The two sets
  share no id and no (file, line) pair, so no row is counted twice.
  A both-tier bin-1 row is one row.

## B. Bin-2 classification

### B1. "Cannot transfer: absent" is a gap, and it leaves the denominator

The lane's own reason text for the 17 2T rows (`loo_R2`, `loo_epred`
and the others) is "a gap, not a decision". The 2P rows for
`discrete_weibull`, `negbinomial2`, `frechet` and `gen_extreme_value`
are open gaps that the ranked list gives to lane famuni (rank 16).
Censored `von_mises` has no CDF. The bin-1 `absent` class has the
same meaning ("a recorded gap ... not a decision"). Every one of
these rows is outside the 551 denominator, so the 87.1% treats
"frmtmb lacks it" as if it were "cannot be asked of frmtmb". The bin-1
ledger never quoted a "transferable" rate; this lane introduced the
metric.

| denominator | rows | pass rate |
|---|---|---|
| as published: all "cannot transfer" out | 551 | 480 / 551 = 87.1% |
| clear gaps back in (bin 1: 2, 2T: 17, 2P: 9) | 579 | 480 / 579 = 82.9% |
| every `absent` back in (2, 17, 76) | 646 | 480 / 646 = 74.3% |

The 67 G1 rows (brms's exported `d/p/q/r` utilities) are a judgment
call. frmtmb does not export these names, which is a decision, not a
gap. Label them that way, not "absent". The fix: split "cannot
transfer" into "no counterpart (decision)" and "absent (gap)", and
quote the rate over everything except the first.

### B2. Twelve ordinal rows are portable through the exported family object

The survey (group table, G3/G4) says ordinal threshold preps "have no
such read" and calls them high cost. G2 gives as its reason that the
ordinal `$lpdf` "takes one linear predictor and the thresholds, not
brms's array". But `brms:::log_lik_cumulative(i, prep)` takes exactly
that: `mu`, `disc` and `thres`. Each frmtmb ordinal family object
exports `ord_lpdf_make(tmap)` and `ord_sim_make(tmap)`. With
`tmap = identity` they take brms's thresholds as given.
`dev/remeasure-rev-ord.R` (seed 446, the shapes of brms's block
`tests.log_lik.R:446`, increasing thresholds; log `.txt`):

    cumulative: dim 50 x 8 finite TRUE     (sratio, cratio, acat the same)
    cumulative against brms:::log_lik_cumulative: max abs diff 6.7e-13
    simcumulative: dim 50 x 8 range 1-4    (sratio, cratio, acat the same)

That covers `log_lik:462` to `:478` (6) and `posterior_predict:307`
to `:327` (6). The G2 blocks `distributions:238` (8 assertions, which
compare `dcumulative()` and the other three against a closed form)
and `:446` (2) have the same inputs and are probably portable too. I
did not run them. `categorical(K =)` has an lpdf as well; I did not
check it. Either port these rows or classify them as "not ported,
transferable", which keeps them in the denominator.

### B3. What the passes assert

`dev/remeasure-rev-2tpass.R` (log `.txt`) lists every bin-2 pass with
its assertion:

- The 23 2T passes are all `dim()`, `names()`, `dimnames()` or
  `is.numeric()` checks on draws, plus one consistency check
  (`methods:1020`, `waic2` equals `waic_pointwise`). This is how brms
  wrote them, and they are frmtmb.sample's passes. Fine, but say so.
- Of the 68 2P passes, 51 are `length()` or `dim()` checks, and the
  adapter sets `n` itself (`brms_adapt_dpars()`). All 35
  `posterior_predict` passes are of this kind. 16 compare a frmtmb
  density value against brms's exported density or base R's: a fair
  test (see below).
- `log_lik:32` (`ll == ll_gaussian * weights[1]`) tests the ADAPTER.
  frmtmb applies weights outside `row_lpdf()`, and the adapter
  reimplements that as `ll * w[i]`. The lane quotes "removing the
  weights line breaks an assertion" as evidence that the instrument
  can say no. It is evidence that this assertion tests adapter code.
  Mark it as such, or route it through frmtmb's own weighting.

Is the adapter a fair test of frmtmb, or brms against itself? It is
fair for the 16 value rows. The answer comes from `frmtmb::brmsfamily()`,
`row_lpdf()` and `sim_draw()`; no brms function computes it. The
reference side is base R or brms's EXPORTED R densities
(`dstudent_t`, `dskew_normal`, ...). One caveat follows from section C:
a brms R-side density is not always the density brms fits.

### B4. "Divergence, policy" (mostly holds)

- `LOO()` and `WAIC()` (442, 443, 1005, 1021): `R/loo.R` documents the
  refusal. This is a real decision, but brms 2.23.0 still runs these
  spellings.
- `loo(fit1, fit1)` and `waic(fit1, fit1)` (446 to 449, 1015, 1016):
  frmtmb.sample refuses them by name. But `dev/brms-api-diff.md`
  (around line 385) lists "`compare =`, `model_names =` ... on
  `loo()` / `waic()` / `LOO()` / `WAIC()` / `loo_subsample()`" as
  **Small** open work, not a decision. Cite a decision record or
  reclassify these six rows as defects. Gap rank 18 has the same
  problem.
- `loo_subsample()` (478 to 480): `dev/arcovsample-findings.md`
  records the refusal. The stated reason ("data too large ... which a
  frmtmb model is not") is weak, but it is recorded.

## C. R6: what each `mu` means (the record is wrong)

`dev/remeasure-rev-compois.R`, seed 20261013, n = 400, drawn from the
exact COM-Poisson pmf with lambda = 9 and nu = 2 (true mean 2.737,
`lambda^(1/nu)` = 3). brms 2.23.0 fits 2 chains of 1000; frmtmb 0.70.0
`compois()`. Log `dev/remeasure-rev-compois.txt`.

| | log link Intercept | exp(Intercept) | shape / nu |
|---|---|---|---|
| frmtmb `compois()` | 1.0234 | 2.7825 (the sample mean, 2.7825) | exp(0.585) = 1.79 |
| brms `com_poisson` fit | 1.988 (posterior mean) | 7.30 | 1.80 |

brms's Stan program (`brms/chunks/fun_com_poisson.stan`) is
`y * log_mu - nu * lgamma(y + 1) - log_Z_com_poisson(log_mu, nu)`,
with series terms `(k - 1) * log_mu - nu * lgamma(k)`. So the `mu` of
a brms FIT is Shmueli's rate lambda (pmf proportional to
`mu^y / (y!)^nu`), and 7.30 is that rate at nu = 1.80. brms's R side,
`brms:::dcom_poisson()` and `brms:::log_Z_com_poisson()` (which
`log_lik()`, `posterior_predict()` and `loo()` use), is
`shape * (x * log_mu - lgamma(x + 1))` with terms
`shape * (k * log_mu - lfac)`, so there `mu = lambda^(1/shape)`.
**brms 2.23.0 is internally inconsistent.** Its post-processing
evaluates a different pmf from the one it fits. frmtmb's `mu` is the
mean (`RTMB::dcompois2`).

So R6 must say: the brms fit's `mu` is lambda; brms's R helper
(the bin-2P reference of `log_lik:302`) has `mu = lambda^(1/shape)`;
frmtmb's `mu` is the mean; all three differ. "brms's own normalizer
is approximate (mass 1.00014)" describes the R helper's asymptotic
branch. The coefficients of the two fits are not comparable
(Intercept 1.02 against 1.99). This needs a "differences from brms"
note, and the brms inconsistency is worth filing upstream.

## D. R1 (confirmed; the attribution holds, the bisection wording does not)

`dev/remeasure-repro-missings.R`, rerun by me (logs
`dev/remeasure-rev-repro-<build>.txt`):

| build | fit_imp2 | fit_imp3 |
|---|---|---|
| r8 (0.70.0) | code 1, false convergence (8), all 6 SEs NaN | code 0, both residual-SD SEs withheld, `bmi_age` 6.0975 |
| ea711734 | the same | the same |
| 58b97a2c | code 0, SEs | code 0, SEs |
| ebd12b88 | code 0, SEs | code 0, SEs |
| 48c1e6c9 (0.69.0, the common parent) | code 0, SEs | code 0, SEs |
| r7 (0.69.1, 2884ccbd) | code 0, SEs | code 0, `bmi_age` 4.5195 |

The log-likelihood is -113.672565129 and -113.674154714 on every
build. The finite-difference reference is PD, with eigenvalues
0.001617 to 5.45e5. It reproduces 0.69.1's SEs, so "35% too large" is
relative to a trustworthy reference: 6.0975 / 4.5195 = 1.349 for
`bmi_age`, and 0.025832 / 0.019397 = 1.332 for `bmi_michl:age`.

**Wording:** `git log` shows d2ef6312 with the parents 2884ccbd,
58b97a2c, ea711734 and ebd12b88, and each of these has the one parent
48c1e6c9. So 58b97a2c is not an ancestor of ea711734, and "58b97a2c
is right, ea711734 is wrong" is not a bisection. The evidence that
holds is this: the common parent is clean, ea711734 alone reproduces
the failure, and the other three merge parents are clean.

**Mechanism (new; for the owner):**
`dev/remeasure-rev-r1units.R` and `dev/remeasure-rev-r1cause.R` (logs
`.txt`):

- At 0.70.0, fit_imp2 carries `par_units` 8 8 8 8 128 128 1 1. These
  are the response units that ea711734's `predictor_par_units()`
  applies to every coefficient of a response, whatever the scale of
  its column. They are the same under `autoscale = FALSE`. 0.69.1
  carries none. This explains why `autoscale = FALSE` does not
  restore the SEs.
- `autoscale_sdreport()` takes `optimHess()` in `q = par / u` at the
  default `ndeps = 1e-3`. That is a step of 0.008 in `bmi_michl:age`
  (SE 0.019, column `age * chl` up to several hundred) and 0.128 in
  the `chl` coefficients:

      q-space optimHess ndeps 0.001: eigen -0.01153 to 34780000; se NaN ...
      q-space optimHess ndeps 1e-04: eigen 1.026 to 34880000; se 6.771 4.516 0.0342 0.01939 22.32 11.76 0.1835 0.175
      raw-space optimHess ndeps 1e-3: eigen 0.001617 545000;  se 6.772 4.519 0.0342 0.01939 22.32 11.76 0.1835 0.175

  So the Hessian goes indefinite from the finite-difference step in
  those units, not from the objective. `autoscale_fit_sound()` takes
  the same Hessian. Any fit with a response unit above 1 and a column
  of large spread can be hit, so the reach is unknown. Lane optstart
  or illse should scan for it. The record names them as the nearest
  lanes; it should also name this cause.

## E. The median `|z|` headline

`dev/remeasure-rev-median.R` (log `.txt`), from the lane's
`plausibility.rds` files:

    only v067: brms_missings.8.2        only r8: brms_distreg.12.1
    as published  v067: n 130 median 0.0488   r8: n 129 median 0.0458
    common set    v067: n 124 median 0.0470   r8: n 124 median 0.0470

The record says missings.8.2 leaving "is the whole change of the
median". It is not: distreg.12.1 enters too. On the common set
nothing changed. Part of the apparent improvement comes from R1
removing a model. Report the common-set median, or say that the two
figures are not comparable.

## F. Families instrument

`families-coverage.R` counts a core family when
`as_frmtmb_family(name)` succeeds or the name is in the registry, and
an extension family when its constructor's default call passes
`as_frmtmb_family()`. Neither fits a model. The record's phrase
"whose default call frm() accepts" is not what the script measures.
The frmtmb.eam wiener probe does fit. `dev/remeasure-rev-famfit.R`
(seed 20261013, n = 300, `frm(bf(y ~ x), family =
frmtmb::brmsfamily(name))`, log `.txt`) shows that 33 of the 34
counted core families fit. `multinomial` errors with "multinomial()
needs the number of categories". Gap rank 26 lists this as open, and
the instrument still counts it. Say the count is name-level, or count
33 plus a note.

## G. Ranked gap list

- Rank 7 (`threshold =` as an `frm()` argument) is "open, postfix",
  although brms 2.23.0 fails the same call. Ranks 13, 15, 19 and 20
  have the same property and are "deliberate: as brms". Make them
  consistent.
- R2 (`loo()` on draws silently ignores `newdata`, a wrong answer with
  no message), R3 and R4 are real defects. They are not in any ranked
  list and have no lane. R2 is the most severe kind of defect. Give
  it to postfix.
- R1 has the cause above. Give it an owner (optstart or illse), not
  "no lane".
- See B4 for rank 18.

## Not checked

- The wiener fit and R5 (I read the lane's log only).
- R2 to R4 (I read `dev/remeasure-out/probe/bin2-probe1.txt`; R2's
  log supports the claim: 40 pointwise rows for 20 rows of `newdata`).
- The brms reference fits were reused, not refitted.

## Reviewer files

`dev/remeasure-rev-wt/` (a copy for the rerun, 369 MB with the Stan
cache), `dev/remeasure-rev-b2/` and `dev/remeasure-rev-b3/` (bin-2
copy roots), `dev/remeasure-rev-tables/`, and the scripts and logs
`dev/remeasure-rev-*.R` and `dev/remeasure-rev-*.txt`. Delete them
at consolidation. Libraries: `C:/Users/adf44/source/r/remeasure-rev-lib/`
holds four commit builds.
