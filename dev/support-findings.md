# Lane support: exact densities outside the support, and the xbeta reference

Round opened 2026-11-09, base e0cce906 (frmtmb 0.80.0). Brief:
`C:/Users/adf44/source/r/round-briefs/r081-support.md`. Review:
`dev/reviews/2026-11-10-support.md` (FIX: 1 blocker, 2 major, 3 minor).
This record covers both rounds; nothing is committed.

## Part A: what frmtmb reports outside a support

### The rule, and the first version withdrawn

The first version (round 1) switched on class: every evaluation on plain
doubles was exact, the tape kept the barrier. That rule is WRONG and is
withdrawn here. Importance sampling (`R/importance.R:1471-1473`,
`imp_verify()` at `:1326`) and nested quadrature (`prime()` at
`R/quadrature-nested.R:563-569`, `quad_nest_verify()` at `:600`, the
`rule_err` pass) evaluate the objective on doubles to set the offsets of
their tapes and to check the tape. With `-Inf` there, nested fits
stopped or moved silently and `importance = 500` fits stopped (review,
finding 1; numbers below). The round-1 claim "the fit is unchanged"
rested on seven fixed-effect fits only and was false for these paths.

The rule now is the inverse. Every evaluation, on the tape and on
doubles, keeps the fitted objective's arithmetic by default. Only a
reader that reports densities asks for the exact value, by calling
`row_lpdf(..., exact = TRUE)`, which sets the option
`frmtmb.exact_support` for that call (restored on exit). The families
read it through `exact_support()` and still keep the tape's arithmetic on
an advector. The switch is on class and on the option, so it resolves
while a tape is built and leaves no branch on a parameter.

Readers that ask for the exact value:

- frmtmb.sample `draws_row_loglik()` (`extensions/frmtmb.sample/R/loo.R`),
  the only source of `log_lik()` of draws, and so of `loo()`, `waic()`,
  `psis()` and every reader of that matrix.
- The ported-suite adapters, which stand in for brms's `log_lik_<fam>()`
  and `d<fam>()` (`dev/brmsport-bin2-adapters.R`: the log_lik adapter,
  the ordinal and categorical log_lik, the G1 density shims).
- The core test helper `frm_row_loglik()`
  (`tests/testthat/helper-brms-methods.R`), which composes `log_lik()`
  for the brms-method tests.

Callers that keep the objective's arithmetic (the default, `exact =
FALSE`): the objective itself (`R/objective.R:1025`), importance sampling
(`R/importance.R:1235`), nested quadrature (`R/quadrature-nested.R:534`,
`:539`), latent-residual derivatives (`R/autocor-latent.R:909`, `:963`),
OSA (`R/predict.R:6101` and the CDF tails), truncated simulation by
inversion (`sim_cdf_pair()`, `R/families-cdf.R`), `nll()` on doubles in
`nl_check_body_read()`, `nl_coef_curved()` and `nl_domain_start()`
(`R/par-template.R`), and `check_custom_family()` (`R/interop.R:74`).
No core path reports a density of an ML fit on doubles: `log_lik()` of a
fit refuses and points to draws, and `logLik()` is `-opt$objective`.

### What changed

- `R/families-extreme.R`:
  - new `exact_support()`, the option reader.
  - `exp_neg_capped()`: exact `exp(-w)` only when a reader asks, on
    doubles; the cap at `exp(700)` everywhere else.
  - new `log_surv_exp_neg()`: `log(1 - exp(-exp(-w)))` that keeps the
    far upper tail, used for a reader only.
  - new `gev_exact()` and `gev_parts_exact()`: the GEV density, CDF and
    log survivor are exact for a reader (`-Inf`; 0 or 1; 0 or `-Inf`
    outside the support). Otherwise the tape's arithmetic.
  - Frechet `lcdf`/`lccdf`: `positive_cdf_q()` at `q <= 0`; the exact
    survivor for a reader.
  - roxygen of `?gen_extreme_value`: who gets the exact value, and the
    two places the fitted objective is not the exact likelihood (a row
    with `t` below 1e-15, and a `trunc()` or `cens()` bound outside the
    support).
- `R/objective.R`: `row_lpdf(..., exact = FALSE)`.
- `R/sampling-api.R`: the `exact` argument documented.
- `R/families.R`: new `positive_cdf_q()`; `exponential()` and
  `weibull()` CDF slots use it.
- `R/families-cdf.R`: `shifted_lognormal_cdf_log()` floors `q - ndt` as
  `1e-300 + ad_pos(d - 1e-300)` (pitfall 11), so a bound below `ndt` is
  `lS = 0` with a zero gradient and a finite `lF`; `-Inf` for a reader.
  `von_mises_cdf_log()` reads a data `q` outside `[-pi, pi]` at the end
  it passed.
- frmtmb.sample `R/loo.R`: `draws_row_loglik()` asks for `exact = TRUE`.
- frmtmb.survival `R/rp-basis.R`: new `sv_pos_q()`; `R/royston-parmar.R`:
  the CDF slots read a time at or below 0 at 1 and weight it 0.

### Evidence: the blocker is fixed

The reviewer's scripts, rerun on the fixed build (logs
`dev/support-out/rev-{is,is2,nq}-lane2.txt`) against the reviewer's base
logs (`dev/review-support-{is,is2,nq}-base.log`):

- `dev/review-support-is.R` (Laplace, `importance = 500`, one-block and
  nested quadrature at seeds 1 to 3, xi 0.5 and 0.8): every line equal
  to the base's once the timings are removed (`diff` empty). Nested seed
  2 xi 0.5 is -242.9981897, beta 1.1569858, as on the base; seeds 2 and
  3 at xi 0.8 and seed 3 at xi 0.5 fit.
- `dev/review-support-is2.R` (`importance = 500` at xi 1 and 1.5, seeds
  1 to 4): equal to the base's; seeds 2 and 4 at xi 1.5 fit, -326.5836299
  and -329.4890012.
- `dev/review-support-nq.R` (the trace of `prime()`): the whole verbose
  trace equals the base's, 16 `prime()` calls with 0 non-finite offsets
  in every one, final logLik -242.998189713707, beta 1.156985773478,
  theta -0.543440065712 -3.749768780143.

Bit for bit, pinned to core 0 (`dev/support-pin.bat`, pitfall 21):
`dev/support-fitsame-re.R` prints `opt$par` and `logLik()` in hex for the
26 random-effect fits of those two scripts (24 at seeds 1 to 3 times xi
0.5, 0.8 times Laplace, importance, one-block and nested quadrature, and
`importance = 500` at xi 1.5, seeds 2 and 4). Logs
`dev/support-out/fitsame-re-{before,after}.txt`: no fit errors on either
build, and the two files are IDENTICAL apart from the `lib:` line. The
seven fixed-effect fits of `dev/support-fitsame.R` (GEV, GEV with `cens()`
and `trunc()`, a two-GEV mixture with a component that misses rows,
Frechet, Frechet with `cens()` and `trunc()`, Weibull and exponential with
`cens()`) are again IDENTICAL to the base
(`dev/support-out/fitsame-after2.txt` against `fitsame-before.txt`).

On rounding: `dev/support-isnq-probe.R 20` fits the new test's two cases
(importance at seed 2 xi 1.5; nested quadrature at seed 3 xi 0.5) on 21
copies of `y` scaled by `1 + 2e-15 k`, k = 0..20 (logs
`dev/support-out/isnq-{before,after}.txt`): 42 fits per build, 0 errors,
and the two builds' outputs are identical (logLik to 12 digits). The
warnings across the 42 fits are the four the test allows: the proposal
covers groups poorly (21), the importance correction used all its rounds
(21), the nested rule did not settle (21), a large gradient at the
optimum (12).

### Evidence: finding 2, each site by name

With the option unset, every changed function computes the base's
arithmetic: the GEV branch is `gev_exact()`, false without the option;
`exp_neg_capped()` and the Frechet survivor test the option too; the
only unconditional changes are `positive_cdf_q()`, `sv_pos_q()`, the
shifted-lognormal floor and the von Mises clamp, which change values
only at a `q` outside the support (a bound), and the shifted-lognormal
floor, which equals the old `ad_max()` to the bit for `q - ndt` above
about 1e-284. So:

- `nl_check_body_read()`, `nl_coef_curved()`, `nl_domain_start()` read
  `nll()` on doubles: the base's values. Frechet at `w = -705`
  (`dev/support-m1-probe.R`, logs `dev/support-out/m1-*.txt`): doubles
  -1.01423e+304 and tape -1.01423e+304 on both builds; the reader alone
  gets -1.50525e+306. With the option forced on for the session (the
  first version's behavior) doubles gave -1.50525e+306 against the
  tape's -1.01423e+304.
- `latres_row_derivs()` and `latres_tk()`: `row_lpdf()` with the default,
  the base's values. No GEV case with a latent residual was constructed.
- `check_custom_family()` on the GEV at test values with one row outside
  the support: the base and the fixed build both stop at "AD gradient
  disagrees with finite differences (max relative error 1e+09)" (the
  barrier's slope); the first version stopped earlier at "lpdf is not
  finite at the test values". The fixed build matches the base.
- importance sampling and nested quadrature: above.

### Evidence: the reported value is exact

`dev/support-gev349-probe.R` rebuilds brms's block of `log_lik:349`
(seed 233; mu is `exp(eta)` at that point, row `log_lik:273`) and reads
the GEV row as the adapter does, `row_lpdf(exact = TRUE)`. Logs
`dev/support-out/gev349-probe-{before,after}.txt`:

    before (rellib-r19)
    draws 25, outside the support (t <= 0) 3
    outside: frmtmb -Inf 0, finite 3, NaN 0; brms NaN 3
    outside, frmtmb: -3.497e+248 -1.285e+133  -6.223e+85
    inside: 22 draws, max rel diff frmtmb vs brms 2.74e-15
    y = 0 at (mu 1, sigma 0.1, xi 0.3) and (mu -1, sigma 1, xi 0.3): -1.00803e+150  -1.55396e+00
    after (wt-support-lib)
    draws 25, outside the support (t <= 0) 3
    outside: frmtmb -Inf 3, finite 0, NaN 0; brms NaN 3
    outside, frmtmb: -Inf -Inf -Inf
    inside: 22 draws, max rel diff frmtmb vs brms 2.74e-15
    y = 0 at (mu 1, sigma 0.1, xi 0.3) and (mu -1, sigma 1, xi 0.3):     -Inf -1.55396

The bin-2 generated files in tier mode on the fixed build with every
extension loaded (`dev/support-bin2-tier.R`, logs
`dev/support-out/tier2-{log_lik,posterior_predict}.txt`): log_lik
pass=104 fail=1, the one failure `log_lik:580` ("recorded as 'error' and
is now 'holds'", below); posterior_predict pass=65 fail=0. distributions
(92), priors (1) and methods (58) pass with fail=0
(`dev/support-out/tier-*.txt`). `log_lik:349` keeps its recorded
`divergence`.

### Bounds outside a bounded support (review finding 3)

`dev/support-poscdf-probe.R` (round 1) and `dev/support-ndt-probe.R`
(seed 402 for the shifted lognormal, as the reviewer's
`dev/review-support-trunc.R`; 403 wiener, 404 rdm, 405 royston_parmar,
406 von Mises), logs `dev/support-out/{poscdf,ndt-probe}-{before,after}.txt`.
Each line is the fit with a bound that truncates nothing against the
fit without it:

    before (rellib-r19)
    exponential      trunc(lb = -1) logLik -623.2796505 vs -498.4377178; Intercept 1.0775988 vs 0.66145906
    weibull, frechet trunc(lb = -1): NA/NaN gradient evaluation
    shifted_lognormal trunc(lb = -1), (lb = 0), (lb = -1, ub = 1e6), cens + trunc(lb = -1), cens + trunc(lb = 0): ERROR NA/NaN gradient evaluation
    slot at q = -1, 0 (ndt 0.02): value 0 0 -Inf -Inf; taped Jacobian 12 of 12 entries NaN
    wiener trunc(lb = -1), (lb = 0)          dlogLik 0  max|dbeta| 0
    rdm trunc(lb = -1), (lb = 0)             dlogLik 0  max|dbeta| 0
    royston_parmar trunc(lb = -1), (lb = 0)  ERROR NA/NaN gradient evaluation
    von_mises trunc(lb = -4)                 dlogLik -6.42  max|dbeta| 0.0693
    von_mises trunc(lb = -pi)                dlogLik 5.68e-14  max|dbeta| 4.56e-12
    after (wt-support-lib)
    exponential, weibull, frechet trunc(lb = -1): identical logLik and Intercept
    shifted_lognormal, all five forms        dlogLik 0  max|dbeta| 0
    slot at q = -1, 0 (ndt 0.02): value 0 0 -238593 -238593; Jacobian of lS 0, non-finite entries 0
    wiener, rdm                              dlogLik 0  max|dbeta| 0
    royston_parmar trunc(lb = -1), (lb = 0)  dlogLik 0  max|dbeta| 0
    von_mises trunc(lb = -4)                 dlogLik 5.68e-14  max|dbeta| 4.56e-12
    von_mises trunc(lb = -pi)                dlogLik 5.68e-14  max|dbeta| 4.56e-12

Checked and already right: gamma, lognormal, inverse.gaussian, Beta,
hurdle_gamma, hurdle_lognormal, zero_inflated_beta (5.7e-14), cox()
(reviewer's logs); wiener and rdm (above). lba, gddm and wiener_gng take
no `trunc()`. The von Mises at `lb = -pi` keeps a 5.7e-14 difference: the
full-circle mass is formed numerically, not as exactly 1.

### Minor findings

- m1: a GEV `trunc(lb = )` below the fitted lower end moves the fit by
  7.7e-12 in logLik and 1.77e-7 in a coefficient, with `cens()` 2.2e-10
  and 1.1e-6 (the reviewer's `dev/review-support-trunc-*.log`), on both
  builds: the tape's clamp leaves `F(lb)` up to about 1e-13. The help and
  NEWS now name this second case. The support's end depends on the
  parameters, so the bound cannot be resolved off the tape.
- m2: the xbeta claim of round 1's report ("frmtmb is closer to the
  256-bit reference than betareg") was wrong and is withdrawn here. The
  measured errors are in Part B: frmtmb is closer only at `y = 1`.
- m3: see "Seen to fail" below.

### Seen to fail

`dev/support-out/t-*.txt` (one process per file, `dev/release/run-tests.R`):

    test-families-extreme.R on rellib-r19:          pass=112 fail=0 err=3
    test-families-extreme.R, first version emulated: pass=120 fail=2 err=1
    test-families-extreme.R on the fixed build:      pass=124 fail=0 err=0
    test-cens-trunc.R on rellib-r19:                 pass=74 fail=2 err=1 warn=1
    test-cens-trunc.R on the fixed build:            pass=99 fail=0 err=0
    test-royston-parmar.R on rellib-r19:             pass=67 fail=0 err=1 warn=3
    test-royston-parmar.R on the fixed build:        pass=71 fail=0 err=0

On rellib-r19 the three GEV reader blocks error with "unused argument
(exact = TRUE)", the weak form; the behavioral failure on the base is the
probe above (-3.497e+248 where the value is `-Inf`). The first version
is emulated by `dev/support-emul-run.R`, which sets
`frmtmb.exact_support` for the whole session, so every evaluation on
doubles is exact, as the first version's class switch made it: the
block "off the reader, doubles keep the fitted objective's arithmetic"
fails twice (non-finite on doubles where the tape is finite), and the
importance/quadrature block errors with "`importance` could not reproduce
this model's joint log-density ... (draw 1: -Inf against -Inf)". On the
base, test-cens-trunc's new block fails on the exponential (logLik 124.84
lower, intercept 0.416 away) and errors at the Weibull's NaN gradient,
which ends the block before the shifted lognormal and the von Mises; the
probe shows those. test-royston-parmar errors at the NaN gradient.

Round 1 recorded `test-families-extreme.R` on the base as `pass=111
fail=7`; the reviewer's rerun gave `pass=112 fail=6`. The round-1 record
described an earlier form of the mixture assertion, which passes on the
base (`exp(-1e150)` is already 0). That record is superseded by the
table above.

### The survey

Read across core and the eight extensions (a read-only search agent,
then checked by this lane). "Measured" marks a value seen in a log; the
others are from reading the code. Kinds: S, the exact value is -Inf (or a
CDF of exactly 0 or 1) and a finite number comes back; U, the exact value
is positive but below what the floor lets through; I, no representable
value changes.

| Site | Kind | Families | Returned vs exact | Readers | This lane |
|---|---|---|---|---|---|
| `R/families-extreme.R`, `core_lpdf`: `- 1e8 * g[["out"]]` with the `gev_parts()` clamp at t = 1e-15 | S | GEV, GEV mixtures | -1.008e150, -5.0e7, +6.4 (measured) vs -Inf | log_lik, loo, waic of draws | exact for a reader; the objective keeps it |
| `R/families-extreme.R`, GEV `lcdf` | S | GEV, `xi > 0` below the lower end | 1.8e-14 at `xi = 10` (measured) vs 0 | cens, trunc, OSA | exact for a reader |
| `R/families-extreme.R`, GEV `lccdf` | S | GEV, `xi < 0` above the upper end | -103.6 at q = 3, `xi = -0.5` (measured) vs -Inf | right cens, OSA | exact for a reader |
| `R/families-extreme.R`, `exp_neg_capped()` | U | Frechet, GEV | -1.01e304 for `w` in (-709.8, -700) (measured) | all | exact for a reader |
| `log(-expm1(-exp(-w)))` in both `lccdf` | U | Frechet, GEV | -Inf (measured) where the value is about -w | right cens | exact for a reader |
| `R/families.R`, exponential CDF at q < 0 | S | exponential | F(-1) = -0.65, log S(-1) = 0.5 (measured); wrong fit | fit and readers | FIXED, tape and doubles |
| weibull and frechet CDF at q <= 0 | S | weibull, frechet | NaN, fit stops (measured) | fit | FIXED, tape and doubles |
| `R/families-cdf.R`, `shifted_lognormal_cdf_log()`, `ad_max(q - ndt, 1e-300)` | S | shifted_lognormal | `log(0)`, NaN gradient, fit stops (measured) | fit | FIXED, tape and doubles |
| `R/families-cdf.R`, `von_mises_cdf_log()` at q outside [-pi, pi] | S | von_mises | an arc of negative length; logLik -6.42 (measured) | fit | FIXED, tape and doubles |
| frmtmb.survival `royston-parmar.R`, CDF slots at q <= 0 | S | royston_parmar | `log(q)` NaN, fit stops (measured) | fit | FIXED, tape and doubles |
| `R/families-extreme.R`, GEV `trunc()`/`cens()` bound outside the support | S | GEV | F within about 1e-13 of 0 or 1; fit moves 7.7e-12 (reviewer) | fit | named in help and NEWS; parameter-dependent end |
| `R/families-cdf.R:32`, `log1m_exp_ad()` floor | S and U | every `cdf_log` family | an empty or unresolvable window: about -690.8 vs -Inf | interval cens, trunc windows, OSA | not changed |
| `R/families-cdf.R`, von Mises zero-length arc floor | S | von_mises | about -690.8 plus finite terms vs -Inf | left cens at -pi, right at pi | exact for a reader at q at or past the ends; floor kept |
| `R/families.R:4803`, ordinal gap floor (`.gap_floor`) | S and U | ordinal mixtures | tied thresholds: about log(1e-300) vs -Inf | all | not changed |
| `R/predict.R:6073`, OSA tails floored at -1e4 | S and U | all, OSA only | -1e4 vs -Inf | OSA (tape) | not changed: documented |
| eam `wiener-family.R:986`, decision time held at delta | S | wiener with `allow_unreachable` | about -1.8e9 vs -Inf | all | not changed |
| eam `wiener-density.R:284`, `:308` | S | wiener variability models | finite where the ndt range passes y | all | not changed |
| eam `wiener-family.R:1059`, `:1070`, `:1108` | S and U | wiener cens/trunc | F = 1e-300 where it is 0 | cens, trunc | not changed; a bound below the support fits exactly (measured) |
| eam `wiener-contaminant.R:35`, `:143` | S | wiener contaminant | -1e300 outside crange; survival held at 1e-300 | fitted rows outside crange are refused | not changed |
| eam `wiener-gng.R:642` | S | wiener_gng go rows | about -1e11 at y <= ndt | all | not changed |
| eam `lba.R:456`, `:336`; `rdm.R:289`, `:348`, `:378`; `gddm.R:1527` | S and U | lba, rdm, gddm | -690.8 at a decision time <= 0 or a density below 1e-300 | all | not changed; rdm with a bound below the support fits exactly (measured) |
| survival `royston-parmar.R`, `sv_floor_pos` | S (parameter) | royston_parmar | about -33.6 where the hazard is not positive | all; `rp_floored()` refuses such fits | not changed |
| survival `royston-parmar.R`, the caps and `sv_squeeze()` | U | royston_parmar | cumulative hazard capped at e^30; F in [5.6e-16, 1 - 5.6e-16] | all | not changed |
| discrete densities at a non-integer or negative y | S | binomial, negbinomial, com_poisson, tweedie, discrete_weibull | finite vs -Inf | none: `valid_y` refuses such y | not changed, unreachable |

Why the rows marked "not changed" stay:

- The eam and survival floors are the extensions' own composed densities.
  Each is a mixture-component device (`allow_unreachable`, the
  contaminant) or a bound the bounded `link_ndt` keeps every fitted or
  sampled row away from, so at a fit or a draw a mixture row is the same
  number either way (`exp(-1.8e9) = exp(-Inf) = 0`). A reader sees the
  floor only at a parameter value the fit cannot reach. With the reader
  switch now in place, a follow-up can make them exact for a reader
  through `exact_support()`'s option without touching their fits.
- `log1m_exp_ad()` mixes two cases doubles cannot tell apart: an empty
  window (mass 0) and a window narrower than the tails resolve.
- The ordinal gap floor: tied thresholds in an ordinal mixture; a
  follow-up.
- The OSA floor is documented in `osa_cdf_tails()`.

### Not done

- No GEV model with a latent residual was fitted to exercise
  `latres_row_derivs()`; the argument above (default arithmetic) covers
  it.
- The eam and survival floors above stay as they are.

## Part B: the xbeta value check (`log_lik:580`)

betareg 3.2-6 is in the user library, so brms's `log_lik_xbeta()` runs.

The one generated block (brms's "log_lik for the xbeta model runs without
errors", seed 568, its own `withr::local_seed()`), in record mode on the
base build rellib-r19: `dev/support-bin2-run.R frmtmb
dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R 1424 1457`, log
`dev/support-out/xbeta-run.txt`, record `dev/support-out/xbeta-rec.tsv`:

    assert frmtmb log_lik:580 pass   TRUE
    value  frmtmb log_lik:580 error  TRUE  HOLDS: 50 values (0 NA, 0 infinite), max abs diff 8.53e-14, max rel 6.43e-14

The value check now HOLDS. Nothing in frmtmb changes for xbeta.

Against the definition in 256-bit arithmetic (Rmpfr; the end masses
through the positive series of `I_x(a, b)`, because Rmpfr's `pbetaI()`
takes whole shapes only): `dev/support-xbeta-ref.R`, log
`dev/support-out/xbeta-ref.txt`, base build:

| case | frmtmb vs 256-bit, max rel | brms or betareg vs 256-bit, max rel | max abs, frmtmb | max abs, betareg |
|---|---|---|---|---|
| brms's block, n 50 | 6.55e-14 | 1.49e-14 (brms) | | |
| y = 0, n 100 | 1.03e-13 | 3.8e-15 | 2.27e-13 | 1.14e-13 |
| y = 1, n 100 | 1.32e-14 | 2.54e-14 | 4.55e-13 | 9.38e-12 |
| y = 0.3, n 100 | 1.3e-12 | 4.69e-14 | 2.39e-13 | 1.14e-13 |

frmtmb against brms on the block: max abs 8.53e-14, max rel 6.43e-14.
Section 2 is seed 581: `mu ~ U(0, 1)`, `phi` log-uniform on (0.5, 500),
`kappa` log-uniform on (1e-3, 2). betareg is closer to the reference in
relative error at the block, at `y = 0` and at `y = 0.3`; frmtmb is
closer at `y = 1`. The largest relative error, frmtmb's 1.3e-12 at
`y = 0.3`, is at a value near 0 (absolute 2.4e-13). Both sides are exact
to the precision of their arithmetic; neither is wrong.

Other rows that needed betareg: brms skips two blocks without it
(`tests.log_lik.R:569`, `tests.posterior_predict.R:215`). The second is
`posterior_predict:228` and `:231`, two length assertions with no value
check, already `pass`, because the port's adapter does not need betareg.
No core or extension test names betareg, and no bin-1 ledger row does.

## Ledger rows that move

- `log_lik:349` (bin 2P, frmtmb): verdict and class stay `divergence` /
  `domain`; the reason reads frmtmb `-Inf` against brms `NaN`, the
  classification of `log_lik:462` and `:478`. Row edited in
  `dev/brmsport-bin2-verdicts-manual-2p.tsv`. The adapter now reads
  through `row_lpdf(exact = TRUE)`, which is what gives `-Inf`.
- `log_lik:580` (bin 2 value check): from `cannot transfer` /
  `reference` to a check that holds. Its row is removed from
  `dev/brmsport-bin2-valuechecks-manual.tsv`, because
  `dev/brmsport-bin2-ledger.R` stops on a manual verdict for a check that
  holds.
- Both edits are also made by `dev/support-bin2-moves.R`, which
  `dev/brmsport-bin2-verdicts-2p.R` sources after the 0.73.0 moves;
  rerunning that script rewrote both tables byte for byte identical to
  the hand edits (`dev/support-out/verdicts-2p-run.txt`; the reviewer
  confirmed it on a scratch copy).
- Until the ledger is regenerated, the tier file still expects
  `log_lik:580` to be `"error"` and fails it ("recorded as 'error' and
  is now 'holds'"); that is betareg's installation. No other bin-2 row
  moves on the fixed build (tier runs above).

## Tests, sessions, OpenBLAS, check (round 2, the fixed build)

All nine suites, one file per process, three at a time, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true`, the lane's library first
(`dev/release/run-par.sh support-full2`, log
`dev/release/support-full2.log`). Every per-file log's `lib:` line names
`wt-support-lib` for frmtmb (509 of 509), and for frmtmb.sample (86) and
frmtmb.survival (15):

    frmtmb files=317 pass=26347 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1854 fail=0 err=0 skip=3 warn=0
    frmtmb.latent files=12 pass=687 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=86 pass=4331 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial files=8 pass=561 fail=0 err=0 skip=0 warn=0
    frmtmb.survival files=15 pass=556 fail=0 err=0 skip=0 warn=0
    SUPPORT-FULL2 ran 509 of 509

The skips are the same 29 as in round 1 (drmTMB gate 13, `test-scale.R`
14, `test-fuzz.R` 1, `test-curve-scale.R` 1). `test-importance.R`,
`test-quadrature-nested.R` and `test-quadrature-defects.R` are among the
317 core files and pass.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
log `dev/support-out/openblas-build2.txt`, 4 threads; `dev/ciharden-out/`
deleted after), logs `dev/support-out/openblas2-*.txt`:

    RESULT test-families-extreme.R pass=124 fail=0 err=0 skip=0 warn=0
    RESULT test-cens-trunc.R pass=99 fail=0 err=0 skip=0 warn=0
    RESULT test-royston-parmar.R pass=71 fail=0 err=0 skip=0 warn=0

`R CMD check --as-cran` (`dev/support-check.sh` for core,
`dev/support-check-pkg.sh <package>` for an extension; logs
`dev/support-out/00check-<package>.log`):

    frmtmb           * checking tests ... [18m] OK     Status: 1 NOTE
    frmtmb.survival  * checking tests ... [11s] OK     Status: 1 NOTE
    frmtmb.sample    * checking tests ... [85s] OK     Status: 1 NOTE

Each NOTE is the expected one: "Skipping checking math rendering:
package 'V8' unavailable". A first core check was killed during its
install step when I stopped the driver that started it, and was rerun;
the line above is the rerun.

One session per suite, three orders (`bash dev/release/run-sessions.sh
support-sess2 2`, the lane's library first, log
`dev/release/support-sess2.log`): 27 of 27 sessions, every one
`fail=0 err=0 warn=0`, 0 DETAIL lines, and no LEAK line without a
`[with load of ...]` tag. The core lines:

    SESSION frmtmb forward files=317/317 pass=26347 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb reverse files=317/317 pass=26497 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb preload files=317/317 pass=26590 fail=0 err=0 skip=17 warn=0

Round 1's runs (`dev/release/support-full.log`, `support-sess.log`, the
first core check) were on the first version and are superseded.

## Version

frmtmb: minor bump. BREAKING: `log_lik()` of draws reports `-Inf` for a
GEV row outside the support; `exponential()`, `von_mises()` with a bound
outside the support now give the right fit. frmtmb.sample: needs the
frmtmb that adds `row_lpdf(exact = )` (a new floor); BREAKING values of
`log_lik()`. frmtmb.survival: patch (a fit that stopped now fits).
DESCRIPTION versions and floors are left to the consolidation.
