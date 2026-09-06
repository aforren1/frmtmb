# Review: wt-eam-followups (frmtmb.eam 0.3.0)

Reviewer lane: rvea. Base commit 68d6782. Branch `wt-eam-followups`, uncommitted.
Private library `.../scratchpad/rvea-lib`, holding the worktree's frmtmb 0.52.0 and
frmtmb.eam 0.3.0 (both printed by every script).

**Round one verdict: PUNCH** (documentation, headed by a vignette section that
stated the opposite of the shipped feature).
**Round two verdict, 2026-09-06: MERGE.** All five punch items closed and verified
on disk and by measurement; suite 1339/0/1, as-cran 1 NOTE (V8). Four text
corrections listed at the end, none of them code, none of them blocking.

## (a) Diff hygiene

15 paths: 14 tracked modifications plus one untracked `dev/eam-followups-findings.md`.
No `dev/*.log` in the diff. The tracked logs under `dev/brms-port/results/` are
untouched.

`extensions/frmtmb.eam/NAMESPACE`: unchanged, and no `@export` tag is added or removed
anywhere in the R diff, so the empty NAMESPACE diff is consistent rather than stale.

`extensions/frmtmb.eam/R/zzz.R`: two hunks, four rows, and nothing else.
`rdm`/`cens()` and `rdm`/`trunc()` move refused -> works; `wiener_gng`/`cens()` moves
refused -> conditional; `wiener_gng`/`trunc()` stays refused with a rewritten reason.
Confirmed as the lane reports.

Reported line numbers verified: `gng_family` :460, `gng_nodes` :436, `gng_lpdf_var`
:593, `nogo_nodes` default :350, `lccdf` :525 in R/wiener-gng.R;
`ddm_nogo_lprob_var` R/wiener-cdf.R:307; `rdm` lccdf/lcdf R/rdm.R:252/:256;
any-of at R/wiener-family.R:382 and R/gddm.R:1040.

One reported line number is wrong: the `bs` init defect is at
`R/wiener-gng.R:489`, not :566 (:566 is a `@noRd` tag). Everything else lands.

## (b) The go branch is bit-identical to wiener(). CONFIRMED

Reproduced on my own data (seed 20260906, `wiener_gng_simulate(500, mu = 1.15,
bs = 1.35, ndt = 0.22, bias = 0.55, deadline = 1.35)`, 413 go / 87 no-go), scoring
each family through its own `family_finalize()` and comparing `sprintf("%a")` bit
patterns.

| data | variability | go rows | bit-unequal |
|---|---|---|---|
| all-go (413) | none, sv, sz, st, sv+sz+st | 413 | 0 for all five |
| mixed (500) | none, sv, sz, st, sv+sz+st | 413 | 0 for all five |

`identical()` is TRUE in all ten cases. The mixed case is the one the lane did not
report and it holds: the go rows' per-row contributions equal `wiener()`'s lpdf on
the go rows alone, exactly. It holds because both finalizers set
`delta = 1e-9 * lo` off the same `lo` (`gng_finalize` R/wiener-gng.R:710 takes
`min()` over go rows only; `ddm_finalize` R/wiener-family.R:419 over the rows it is
given), so a mixed data set and its go subset produce the same `delta`.

The no-variability case also holds, which is not free: `gng_lpdf()`
(R/wiener-gng.R:596) applies `ddm_floor(y - ndt, 1e-12)` where `wiener()` passes
`y - ndt` raw. Measured, `ddm_floor(x, 1e-12) == x` on 20006 positive x, so the
floor is inert. The identity rests on that; it is not structural.

## (e) ddm_cdf_ks 12 -> 4. CONFIRMED, worst-case bound exactly zero

Reproduced on my own grid of 2352 rows (a in {0.8, 1.4, 2.5, 4}, v in
{-2, -0.5, 0.5, 1, 2, 5}, w in {0.05, 0.25, 0.45, 0.5, 0.75, 0.9, 0.95}), with the
normalized time u = t/a^2 swept deliberately THROUGH the blend centre
`ddm_cdf_u0 = 0.02` and below it: u in {1e-4, 1e-3, 3e-3, 0.01, 0.015, 0.02, 0.03,
0.05, 0.1, 0.15, 0.197, 0.3, 1, 3}.

Values, every truncation K in {2, 3, 4, 6, 8, 12, 20} against K = 12: 0 bit-unequal
of 2352, max|d| = 0. Exact AD gradients (RTMB tape, d/d(v, a, w)) on 400 sampled
rows: K = 4 and K = 2 both 0 bit-unequal of 1200, max|d| = 0. No non-finite values.

Worst-case error bound relative to the old K = 12: exactly zero, on this grid.
The mechanism checks out: the largest u still carrying a small-route weight above
1e-16 is 0.15, and at that u the j = 2 term is already below 1e-22. Note K = 2 is
also bit-identical, so 4 is margin, as the comment claims.

## (c) The no-go probability

### Mass conservation. HOLDS where the model is defined, FAILS where it is not

My own quadrature: `integrate()` on `exp(lpdf)` with `dec = 1` over (0, deadline),
split at the non-decision time, `rel.tol = 1e-12`, plus `exp(lpdf)` at `dec = 0`.
Deadline 1.5, mu 1.0, bs 1.4, ndt 0.25.

| parameters | go + no-go - 1 |
|---|---|
| plain, bias 0.50 | +8.9e-16 |
| plain, bias 0.85 | +2.2e-16 |
| sz = 0.3, bias 0.50 | +2.2e-16 |
| st = 0.20 | +2.0e-10 |
| st = 0.45 | -9.9e-08 |
| sv = 0.6 | +3.9e-14 |
| sv = 2.0 | -8.1e-05 |
| all three (sv .8, sz .3, st .2) | +1.4e-09 |
| **sz = 0.5 at bias 0.85** | **-2.9e-02** |
| **sz = 0.9 at bias 0.90** | **-1.6e-01** |
| **all three, bias .85, sv 1.2, sz .5, st .4** | **-2.5e-02** |

The lane's headline "4.9e-12" is its best case, not its worst. With DEFAULT nodes I
measure 9.9e-08 at st = 0.45 and 8.1e-05 at sv = 2.0, and the wide-sz rows are off
by percents. See finding F1.

### nogo_nodes defaults, doubled and quadrupled

Change in the no-go log probability from the default `c(sv = 15, sz = 7, st = 7)`:

| case | x2 | x4 |
|---|---|---|
| sv = 0.6 | -1.6e-13 | -1.6e-13 |
| sv = 2.0 | +2.3e-04 | +2.3e-04 |
| sv = 2.0, bs = 4.0 | -7.9e-04 | -7.9e-04 |
| st = 0.45 | -2.0e-15 | -1.6e-15 |
| sz = 0.3 | -4.4e-16 | 0 |
| all three (sv .8) | -8.0e-11 | -8.0e-11 |
| all three, bs 3, sv 1.5 | +3.8e-06 | +3.8e-06 |

x2 and x4 agree, so x2 is converged and the DEFAULT carries the error. `sz` and `st`
at 7 are saturated, as documented. `sv` at 15 is not, and the Rd says so and tells
users to raise it to 31; the measured sizes match the Rd's table. Defaults are
defensible.

### The NaN, and the EMC2 st convention

The NaN reproduces at the underlying function, `ddm_nogo_lprob(1.25, 1, 1.4, w)`:
w = 0.95 gives -4.441, w = 1.00 gives -41.04, w = 1.05 gives NaN (with
`In log1p(-w) : NaNs produced`), w = -0.05 gives +0.1412, a probability above one.
The table in the `ddm_wclamp_dev()` comment is exactly right.

The EMC2 `st0` convention IS documented for users: `man/wiener_gng.Rd:263-265` says
EMC2 takes `st0` uniform on `[t0, t0 + st0]` while this family centres `st` on
`ndt`, and gives the shift (`ndt = t0 + st0/2`) and the agreement with and without
it. A user porting parameters is told. I could NOT re-run the EMC2 comparison:
EMC2 is not installed in the user library or the R site library (only in sibling
lanes' scratch libraries), so that number is unverified here.

## (d) ddm_wclamp_dev(). Tape-safe and exactly inert, but see F1

`R/wiener-cdf.R:389`: `ddm_floor(ddm_smin(d, 1 - ddm_w_edge - w), ddm_w_edge - w)`,
over `ddm_floor(x, lo) = lo + 0.5*((x-lo) + abs(x-lo))` (R/wiener-density.R:124) and
`ddm_smin(x, cap) = 0.5*(x + cap - abs(cap - x))` (:212). Branch-free: `abs()` only,
no `if`, no comparison on a parameter, so it tapes.

Exactly inert at d = 0, measured at w in {0.05, 0.25, 0.50, 0.58, 0.85, 0.95}:
`clamp(0, w)` is `0x0p+0` and `identical(., 0)` is TRUE for every one. The algebra:
`smin(0, hi)` is `0.5*(hi - |hi|)` = 0 exactly for hi > 0, then `floor(0, lo)` is
`lo + 0.5*(-lo + |-lo|)` = `lo - lo` = 0 exactly for lo < 0.

For a NON-zero in-range deviation it is not inert: at w = 0.5 the round trip moves d
by up to 5.6e-17, one ulp. That is a rounding, not a bias, and it cannot reach the
go branch, which is unclamped.

Out of range it clamps correctly: at w = 0.85, `clamp(0.2)` is 0.15 = `1 - 1e-9 - w`
and `clamp(-0.9)` is -0.85 = `1e-9 - w`.

## (f) lccdf and lcdf. ALL CONFIRMED

Note for readers: core's `lcdf` slot holds the CDF on the PROBABILITY scale despite
its name (`gaussian`'s is a bare `RTMB::pnorm()`, and `objective.R:126` takes
`log(Fv - Flb)` itself). `lccdf` is a true log. rdm() follows both conventions.

Algebra, against numerical integration of the family's OWN density summed over all
winners (rdm has no external reference here):

| case | worst rel err in S | in F |
|---|---|---|
| rdm(2) v=(3,2) | 2.4e-09 | 2.2e-14 |
| rdm(3) v=(2.5,2,1.5) | 5.5e-07 | 1.3e-14 |
| rdm(2) slow drifts | 7.3e-14 | 4.8e-10 |

(The larger S numbers are my quadrature, not theirs; the F column is the sharp test.)

Machinery, hand-assembled log likelihood at the fitted parameters vs `logLik(frm())`:

| model | n | detail | rel |
|---|---|---|---|
| rdm(2) `cens(cc, y2)` | 600 | all four codes: 324/125/87/64 none/right/left/interval | 2.4e-16 |
| rdm(2) `cens(cc)` | 600 | right only, 166 censored | 0 (exact) |
| rdm(2) `trunc(lb = 0.30)` | 389 | | 5.3e-16 |
| rdm(2) `trunc(ub = 1.20)` | 600 | not reported by the lane; works | 4.1e-16 |

wiener_gng `cens()`: 500 trials, 405 go / 95 no-go, deadline 1.4. Scored as a
go/no-go model and as an all-go model with the no-go rows right-censored at the
deadline: logLik -172.793662634772431 both ways, difference EXACTLY 0,
`identical()` TRUE, and `fixef()` identical. The last-bit claim holds.

wiener_gng refusals: left cens, interval cens and `trunc()` are all refused. The
message is core's generic "cens()/trunc() need a family with a CDF (currently:
gaussian, lognormal, ...)". It explains the mechanism but never names `wiener_gng`
and never says the omission is deliberate. See F4.

## (g) any-of adoption. Messages are adequate; the pin count is wrong

New refusals, measured:

- `wiener`: "wiener: the density needs one of `dec` or `vint1`, which nothing on
  this response supplies. Write the addition term: rt | dec(<column>) ~ ..."
- `gddm`: the same sentence under `gddm:`.
- `wiener_gng`: "the density needs `dec`" and, with no deadline, "the density needs
  `vreal1`".
- `ddm_check_variability(what =)` now names the calling family, so
  `wiener_gng(variability = "sq")` says `wiener_gng():` and not `wiener():`. A real
  improvement.

The bar in the brief (name the missing terms and the alternatives) is met. What is
lost against the old hand-written text: the factor convention ("a factor whose
second level is the upper boundary") and, for gddm, the `vint(upper, cond)` ORDER.
See F5.

The pin changes are NOT four and NOT all text-only. Eight assertions changed plus
one added, over three files: `test-defects.R` 2 regex; `test-family.R` 2 regex plus
a NEW `expect_identical(wiener()[["required_aterms"]], list(c("dec","vint1")))`;
`test-gddm-family.R` 2 regex plus 2 `expect_identical` that change the asserted
DATA STRUCTURE from `character(0)` to `list(c("dec","vint1"))` and from `"vreal1"`
to `list(c("dec","vint1"), "vreal1")`. Those two are behavior pins, not text.

## (h) The bs init defect. REAL, but I could not make it change an answer

Reproduced. `R/wiener-gng.R:489` is `bs = function(y, aterms) log(1.4)`, and the
contract at `R/families.R:185` says each `init_dpars` value "goes straight through
that dpar's `linkfun`", which this package restates at `R/rdm.R:414` ("These are
RESPONSE-scale values"). Measured:

    wiener_gng init_dpars$bs -> 0.3364722366   (eta = log() = -1.089240)
    wiener()   init_dpars$bs -> 1.5000000000   (eta = log() =  0.405465)

So the separation starts at 0.336 where 1.4 was meant, a factor 4.16 low. Confirmed.

Impact, seven designs, shipped start vs a start patched to 1.4 through a wrapped
`family_finalize`:

| design | shipped logLik / iters | fixed logLik / iters | same optimum |
|---|---|---|---|
| bs 1.4, n 400 | -182.170269 / 14 | -182.170269 / 15 | yes |
| bs 3.5, n 300 | -56.870614 / 16 | -56.870614 / 15 | yes |
| bs 5.0, n 200 | 0.391520 / 63 | 0.391520 / 58 | yes |
| bs 0.6, n 150 | 45.246987 / 10 | 45.246987 / 10 | yes |
| sz on, bs 2.0, n 400 | -246.305877 / 34 | -246.305877 / 30 | yes |
| st on, bs 2.5, n 300 | -217.959522 / 32 | -217.959522 / 29 | yes |
| bs 1.4, n 1200 | -400.6212241388 / 13 | -400.6212241388 / 17 | yes |

Every design converges to the same optimum to 1e-6 or better and the same `bs` to
four decimals. The cost is 0 to 5 iterations, and on the n = 1200 case the WRONG
start took FEWER. It is a defect against the documented contract and it should be
fixed, but on this evidence it is a one-word cleanup, not a merge gate. See F6.

## (i) sv weak identification. Claim reproduced; the attribution is not established

3000 trials, mu 1.2, bs 1.4, ndt 0.25, bias 0.5, sv 0.6, deadline 1.5, 2358 go.

Profile over `sv` with everything else AT the truth (the Rd's claim that the surface
does peak at the truth):

    sv=0.001 -1131.495  sv=0.05 -1131.281  sv=0.20 -1128.286
    sv=0.40  -1121.228  sv=0.60 -1115.913 (truth, max)  sv=0.80 -1117.210
    sv=1.20  -1146.956

And the higher-likelihood corner: at sv = 0.001, mu = 1.04, bs = 1.34 the log
likelihood is -1114.138, HIGHER than the truth by 1.775. Both halves reproduce.

What does NOT reproduce is the attribution. The Rd says "it is the design rather
than the likelihood" and contrasts wiener(), "seeing BOTH boundaries". On a
two-boundary data set from the same generative parameters I find the same corner:
truth -1300.450, and sv = 0.001 with mu = 1.04, bs = 1.36 gives -1298.098, HIGHER
by 2.353. My search fixes `ndt` and `bias` at the truth and is a grid rather than a
fit, so it does not refute the Rd's fitted numbers; it does mean the contrast is
not demonstrated by the surface. See F7.

The caution IS in the Rd (`man/wiener_gng.Rd:273`, section "What a go/no-go design
can and cannot identify"). It is NOT in the vignette. See F2.

## (k) Cross-lane

`R/ddm-shared.R`: untouched. No diff against 68d6782, absent from `git status`.

`R/zzz.R`: two hunks, 7 lines added and 7 removed, both inside
`rdm_gng_compat_rules()`, at `@@ -182,10 +182,10` and `@@ -211,10 +211,10`. The four
rows are rdm/`cens()`, rdm/`trunc()`, wiener_gng/`cens()` and wiener_gng/`trunc()`.
No `.onLoad`, no `register`, no aterm-registration line appears on either side of
the diff. Clean for a merge with wt-protocol.

One coupling to flag: the new `required_aterms = list(c("dec", "vint1"))` in
wiener(), gddm() and wiener_gng() names the registered aterm `dec` as a string. If
wt-protocol renames or re-registers that term, these three declarations and the six
message pins move with it.

## Additional diff-hygiene measurements

Roxygen IS idempotent. Copied `extensions/frmtmb.eam` to a scratch directory, ran
`roxygen2::roxygenise()` (roxygen2 8.1.0, matching `Config/roxygen2/version`), and
diffed: no change in `man/` and none in `NAMESPACE`. One pre-existing warning is
emitted, `gddm.R:169: @details Could not resolve link to topic "gd_tri_df"`, at a
line this lane does not touch (its gddm hunks are at 1014 and 1085).

`wiener_gng_simulate.Rd` and `rdm.Rd` regenerate to what is committed, so the three
changed Rd files are in step with the R sources.

## Findings, ranked

### F1. MEDIUM. The Rd and NEWS claim an identity the lane's own clamp breaks

`R/wiener-gng.R:141-144`, shipped as `man/wiener_gng.Rd:192-195`:

> The identity that ties the two branches together survives: the go density
> integrated to the deadline plus the no-go probability is still exactly one,
> because both are the same average of the same pair.

After `ddm_wclamp_dev()` the stated REASON is false. The no-go branch averages
CLAMPED start points (R/wiener-cdf.R:318) and the go branch averages unclamped ones,
deliberately, because clamping it would break the `wiener()` bit-identity
(R/wiener-cdf.R:370-374 says so). They are no longer the same average of the same
pair, and the mass is not conserved where the two disagree.

Measured (script `rvea-c-mass.R`, quadrature `rel.tol = 1e-12`):

    sz = 0.5 at bias 0.85   go 0.919017 + nogo 0.051622 = 0.970639   (-2.9e-02)
    sz = 0.9 at bias 0.90   go 0.780184 + nogo 0.063186 = 0.843370   (-1.6e-01)
    all three, biased, wide                             = 0.975118   (-2.5e-02)

And even in range "exactly one" is only true to the quadrature: with DEFAULT
`nogo_nodes` I measure -9.9e-08 at `st = 0.45` and -8.1e-05 at `sv = 2.0`.

NEWS repeats it: "The clamped value is the correct limit, not a fudge". It is the
correct limit for the branch it is applied to. The pair is defective there.

Mitigating, and why this is MEDIUM and not high: the failing region is one where the
start-point range crosses a boundary, which is not a model. `wiener()` itself is
already broken there and worse - measured, `wiener()`'s lower-boundary lpdf returns
NaN at `sz = 0.5, bias = 0.85` for every t at or below `ndt` (4 of 12 probe points),
while `wiener_gng()`'s go branch stays finite everywhere I probed. The clamp is an
improvement on a NaN. The defect is that two user-facing documents say the region is
now handled correctly.

Fix: qualify the Rd sentence and the NEWS bullet. State that the identity holds
while the `sz` range stays inside the boundaries, that outside it only the no-go
branch is clamped, and give the size of the gap.

### F2. MEDIUM. The vignette states the opposite of the shipped feature

`vignettes/ddm.Rmd:788-798`, section "What it does not do", unmodified by this lane:

> `wiener()` offers Ratcliff's three across-trial variability parameters through
> `variability =`. This family offers none of them, and the reason is the no-go
> branch. ... shipping two of three would make `variability =` mean something
> different here than it does on `wiener()`. EMC2's `DDMGNG` does carry all three,
> by integrating them numerically in compiled code.

`wiener_gng()` now offers all three under exactly that argument. The vignette is the
document a user reads first, it ships in the tarball, it renders into the pkgdown
site, and it currently tells them the feature does not exist and argues that
shipping it would be wrong.

The `sv` weak-identification caution is also absent from the vignette; it exists
only in the Rd. The vignette's go/no-go section is where a user meets the design.

Reproduce: `sed -n '786,800p' extensions/frmtmb.eam/vignettes/ddm.Rmd`.

Fix: rewrite that subsection, add the `sv` caution and the `nogo_nodes` argument,
and rebuild the site. This is the only finding I would call a merge gate on its own.

### F3. LOW. The lane's own report miscounts three things

- The `bs` init defect is at `R/wiener-gng.R:489`, not `:566`. Line 566 is a
  `@noRd` tag.
- "Four test pins changed" is eight assertions changed plus one added, over three
  files, and two of them are not text: `test-gddm-family.R:33,35` change
  `expect_identical(gddm()[["required_aterms"]], character(0))` to
  `list(c("dec","vint1"))` and `"vreal1"` to `list(c("dec","vint1"), "vreal1")`.
  Those pin a data structure, not a message.
- "go mass plus no-go probability = 1 to 4.9e-12" is a best case. See F1.

### F4. LOW. wiener_gng's censoring refusals do not name the family

Measured message for left cens, interval cens and `trunc()` on `wiener_gng`:

> cens()/trunc() need a family with a CDF (currently: gaussian, lognormal, poisson,
> exponential, weibull, inverse.gaussian, cox). The list is not closed: a family
> supplies one through the lcdf argument of frmtmb_family(), and a family that only
> ever sees RIGHT censoring may supply the log survivor function through lccdf
> instead

This is core's generic text. It never says `wiener_gng`, and it reads as "this
family forgot to supply a CDF" when the compat row and the Rd both say the omission
is a deliberate modelling decision. A user hitting it will reasonably file a bug.
The zzz.R row's "refused by name" is true of the refusal being loud, not of the
message naming the family.

Fix (optional, cheap): the family could carry a `trunc`/left-cens refusal of its own
that says why, the way `post$mean_fn` already refuses with a reason.

### F5. LOW. The any-of adoption drops two pieces of user guidance

The retired `wiener` message explained the factor convention: "where `decision` is a
factor whose second level is the upper boundary, or a 0/1 column". The retired
`gddm` message gave the argument ORDER: "vint(upper, cond) carries the same pair as
plain integers and also works, boundary first and condition second". Neither
survives in the declared refusal, and `test-gddm-family.R` dropped the pin that
asserted the ordering was mentioned (`"vint(upper, cond)"` replaced by
`"rt | dec(<column>) ~"`).

A gddm user who supplies `vint()` in the wrong order now gets nothing at the point
of failure. The brief's bar (name the missing terms and the alternatives) IS met;
this is a regression in helpfulness, not in correctness.

### F6. LOW. The bs init defect is real and harmless on this evidence

`R/wiener-gng.R:489`. See section (h): the start is 0.336 where 1.4 was meant, and
across seven designs including bs = 5.0 and variability-on models every fit reached
the same optimum, costing 0 to 5 iterations. Fix it (one word, `log(1.4)` to `1.4`)
and re-pin, but I found no evidence it can produce a wrong answer.

### F7. LOW. The sv attribution in the Rd is not demonstrated

The Rd says the weak identification "is the design rather than the likelihood" and
contrasts `wiener()` "seeing BOTH boundaries". On two-boundary data from the same
generative parameters my grid finds the same corner: sv = 0.001 beats the truth by
2.353 there, against 1.775 for the go/no-go data. My search is a grid with `ndt` and
`bias` fixed, not a fit, so this does not refute the Rd's fitted numbers. It does
mean the contrast rests on a single pair of fits. Either widen that evidence or
soften the claim.

### F8. INFO. rdm's trunc(lb) normalizer degrades in the tail

Core forms the left-truncation normalizer as `1 - F(lb)` on the probability scale
(`R/objective.R:80-84`, then `ll - log(Fub - Flb)`), so the declared `lccdf` is not
used for it and the cancellation returns. Measured on rdm(2), v = (3,2), A = 0.8,
k = 0.5, ndt = 0.15:

| lb | exact log S(lb) | core's log(1 - F(lb)) | error |
|---|---|---|---|
| 0.30 | -0.40597 | -0.40597 | 0 |
| 2.00 | -14.02283 | -14.02283 | 5.0e-11 |
| 3.00 | -21.50636 | -21.50636 | -1.1e-07 |
| 4.00 | -28.74648 | -28.74656 | -8.7e-05 |
| 5.00 | -35.83953 | -36.04365 | -0.204 |

This is a KNOWN core limitation, documented at `R/families.R:278-285` ("a
LEFT-TRUNCATED survival model ... meets the identical representability problem from
the other side"). Not this lane's bug. But the zzz.R row now says rdm/`trunc()`
"works" and cites one fit at lb = 0.30, where log S is -0.41. The row should name
the range, as the other rows in that table do.

### Out of scope, for the owner

`wiener()`'s lower-boundary density returns NaN when `sz` pushes the start point
past a boundary. Measured at mu 1.0, bs 1.4, ndt 0.25, bias 0.85, sz 0.5: the
`dec = 0` lpdf is NaN at t in {0.001, 0.1, 0.24, 0.2500001} and finite above. Same
at bias 0.90, sz 0.9. `R/wiener-density.R` is not in this lane's diff, so this
predates it, but it is the same failure mode `ddm_wclamp_dev()` was written for and
it is still live in the sibling family.

### F1, severity revised DOWN to LOW after one more measurement

I tested whether the defective-mass region can attract an optimizer. 2000 trials,
bias 0.85, true `sz` = 0.20 (range 0.75 to 0.95, inside), profiling `sz` with
everything else at the truth. The start-point range leaves the boundary at
`sz` = 0.30:

    sz=0.05  1659.029    sz=0.20  1866.687 (truth)   sz=0.29  1725.500
    sz=0.30  1694.999    sz=0.31  1661.545 OUTSIDE   sz=0.45  1576.482 OUTSIDE
    sz=0.60  1520.269    sz=0.90  1265.186 OUTSIDE   sz=0.99  1278.646 OUTSIDE

    max inside  (sz <= 0.30): sz = 0.190, ll = 1867.975
    max outside (sz >  0.30): sz = 0.305, ll = 1678.640

The interior peak beats everything outside by 189 log units and the surface is
monotone downhill across the crossing. The clamp behaves exactly as the comment
intends: it turns a NaN into a finite, heavily penalized value that an optimizer
walks away from. There is no spurious optimum and no gradient pulling a fit out of
the valid region.

So F1 is a documentation defect only. The numbers in it stand; the consequence does
not. It stays on the punch list as a wording fix, not as a numerical concern.

## (j) Suites

Every file run in a FRESH Rscript process, `testthat::test_file()` with a silent
reporter, `NOT_CRAN=true` so that `skip_on_cran()` does not hide work.

The seven named files:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-rdm-gng.R | 221 | 0 | 0 | 0 |
| test-defects.R | 59 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 |
| test-variability.R | 140 | 0 | 0 | 0 |
| test-gddm-family.R | 77 | 0 | 0 | 0 |
| test-family.R | 35 | 0 | 0 | 0 |

All seven match the lane's counts exactly. My first pass got 210 and 122 for
rdm-gng and variability because `test_file()` leaves `NOT_CRAN` unset and four and
five `skip_on_cran()` tests fired; with it set they run and reach 221 and 140.

Core `tests/testthat/test-compat.R`: 266/0/0 with frmtmb alone, **273**/0/0 with
frmtmb.eam ALSO loaded, matching the lane. The 7 extra assertions are the eam
compatibility rows, so the four rewritten rows are genuinely exercised. Worth
recording that the lane's count only reproduces with the extension loaded.

## (j) as-cran

Consolidation's invocation: `_R_CHECK_CRAN_INCOMING_=false`,
`_R_CHECK_FORCE_SUGGESTS_=false`, `RSTUDIO_PANDOC` prepended to PATH (pandoc 3.8.3
resolved), tarball built with vignettes (`R CMD build` reported
`creating vignettes ... OK`, 331423 bytes).

**Status: 1 ERROR, 1 WARNING, 2 NOTEs.** Same as the lane. Every one attributed by
measurement:

| item | check | measured cause |
|---|---|---|
| WARNING | PDF version of manual | `Rdlatex.log` verbatim: `Error in texi2dvi(...) : pdflatex is not available`. `Sys.which("pdflatex")` from R returns the empty string. |
| ERROR | PDF version of manual without index | Same call, same log, same cause. |
| NOTE | HTML version of manual | Verbatim: "Skipping checking math rendering: package 'V8' unavailable". `requireNamespace("V8")` is FALSE; V8 is absent from the user library and the R site library. |
| NOTE | non-standard things in the check directory | Found `frmtmb.eam-manual.tex`. That file is the intermediate the failed `texi2pdf()` left behind. Downstream of the WARNING/ERROR, not independent. |

The brief's premise about pdflatex is mistaken, and the evidence is worth recording
so the next reviewer does not repeat it. `frmtmb.eam-Ex.pdf` IS produced (3535
bytes), but it is the example-graphics file R's `pdf()` device writes while running
examples. It is not a LaTeX product and its presence says nothing about pdflatex.
On this machine pdflatex is genuinely absent: not on PATH, empty from
`Sys.which()`, and the check's own LaTeX log names it. The R package `tinytex` IS
installed but its TeX distribution is not, so `tinytex::install_tinytex()` would
clear the ERROR, the WARNING and the second NOTE at once.

Nothing in the package fails. `checking Rd files`, `Rd cross-references`,
`Rd contents`, `code/documentation mismatches`, `examples`,
`examples with --run-donttest`, `tests` [502s] and
`re-building of vignette outputs` [73s] are all OK.


### Extra check: AD gradients of the go/no-go likelihood

Not asked for, but it underpins every fit. 120 trials (84 go), RTMB tape jacobian
against `numDeriv::grad(method = "Richardson", r = 6)` on the summed log
likelihood, at parameters away from the truth:

| variability | max abs | max rel |
|---|---|---|
| sv | 2.8e-08 | 4.7e-09 |
| sz | 1.2e-08 | 6.0e-09 |
| st | 3.2e-07 | 3.1e-08 |
| sv+sz+st | 3.4e-07 | 3.5e-08 |

Agreement confirmed, no NaN on any tape. My numbers are looser than the lane's
4e-10; that is the finite-difference reference on a 120-row sum, not the tape.

### Whole suite, all 17 files, fresh process each, NOT_CRAN=true

| file | pass | fail | error | skip | s |
|---|---|---|---|---|---|
| test-brms-parity.R | 13 | 0 | 0 | 0 | 2.6 |
| test-defects.R | 59 | 0 | 0 | 0 | 16.9 |
| test-density.R | 132 | 0 | 0 | 0 | 0.8 |
| test-family.R | 35 | 0 | 0 | 0 | 3.3 |
| test-gddm-family.R | 77 | 0 | 0 | 0 | 3.6 |
| test-gddm-gradients.R | 22 | 0 | 0 | 0 | 84.0 |
| test-gddm-recovery.R | 28 | 0 | 0 | 0 | 266.0 |
| test-gddm-reference.R | 166 | 0 | 0 | 0 | 48.5 |
| test-gddm-solver.R | 94 | 0 | 0 | 0 | 30.4 |
| test-lba.R | 109 | 0 | 0 | 0 | 36.1 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 | 0.7 |
| test-moments.R | 29 | 0 | 0 | 0 | 2.0 |
| test-rdm-gng.R | 221 | 0 | 0 | 0 | 240.2 |
| test-sampling.R | 94 | 0 | 0 | 1 | 4.9 |
| test-simulate-density.R | 69 | 0 | 0 | 0 | 46.7 |
| test-surface.R | 47 | 0 | 0 | 0 | 8.2 |
| test-variability.R | 140 | 0 | 0 | 0 | 321.5 |
| **TOTAL** | **1339** | **0** | **0** | **1** | |

**1339 / 0 / 1**, matching the lane exactly, file for file. The single skip is in
`test-sampling.R`. Independently, `R CMD check` ran the same suite with
`NOT_CRAN=false` and reported `checking tests ... [502s] OK`.
## Verdict

**PUNCH.** Nothing here is numerically wrong. Every substantive claim I could test
reproduced, several of them on wider grids than the lane used, and two reproduced to
the last bit. The suites and as-cran are clean once the environment is discounted.

The punch is documentation. One shipped user-facing document states the opposite of
the shipped feature, and two more overstate an identity this lane's own change
broke. A user reading the vignette today is told the feature does not exist and
given an argument for why shipping it would be wrong. That should not go out.

### Punch list

1. **`vignettes/ddm.Rmd:788-798`.** Rewrite "What it does not do". It currently
   says wiener_gng offers none of `sv`, `sz`, `st`, and argues that shipping a
   subset would be wrong. Replace with what shipped, point at `nogo_nodes =`, and
   carry over the `sv` weak-identification caution that currently lives only in the
   Rd. Rebuild the vignette and the site. (F2)
2. **`R/wiener-gng.R:141-144`** (shipped as `man/wiener_gng.Rd:192-195`) **and the
   matching NEWS bullet.** Drop "still exactly one, because both are the same
   average of the same pair" and "The clamped value is the correct limit, not a
   fudge" as unqualified statements. Say that the identity holds while the `sz`
   range stays inside the boundaries, that outside it only the no-go branch is
   clamped so the pair loses up to 16 percent of its mass, and that the region is
   strictly downhill so the clamp acts as a barrier. Also stop saying "exactly one"
   where the quadrature gives 8.1e-05 at `sv` = 2.0 with default nodes. (F1)
3. **`R/wiener-gng.R:489`.** `log(1.4)` to `1.4`. `init_dpars` takes response-scale
   values, per `R/families.R:185` and this package's own `R/rdm.R:414`. Re-pin
   `test-rdm-gng.R`. Ships in the same round, not a later one: the fix is one word,
   and it costs iterations today and could cost more on data nobody has fit yet.
   (F6)
4. **`R/zzz.R`, the rdm `trunc()` row.** It says "works" on the strength of one fit
   at `lb` = 0.30. Name the range: the normalizer is core's `1 - F(lb)` on the
   probability scale, exact where `log S(lb)` is around -0.4 and 20 percent wrong by
   -35.8. Core documents the limitation at `R/families.R:278-285`; the row should
   point at it. (F8)
5. **Correct the lane's own report** before it becomes the record: the `bs` defect
   is at `:489` not `:566`; the pin change is eight assertions plus one addition
   over three files, two of them structural rather than textual; and 4.9e-12 is a
   best case, not a bound. (F3)

Optional, not blocking: a `wiener_gng`-specific refusal for left/interval censoring
and `trunc()` so users do not read core's generic "needs a family with a CDF" as an
oversight (F4); restore the factor-level convention and gddm's `vint(upper, cond)`
ordering somewhere a user meets them (F5); widen or soften the Rd's attribution of
the `sv` corner to the go/no-go design (F7).

### What I verified and would merge on

- The go branch is bit-identical to `wiener()` on my own data, all five variability
  combinations, all-go AND mixed. Ten of ten, zero bit-unequal.
- `ddm_cdf_ks` 12 to 4 is free: bit-identical values for K in 2..20 and
  bit-identical AD gradients, on 2352 rows spanning the blend centre. Worst-case
  bound exactly zero.
- Mass conservation holds to 8.9e-16 where the model is defined; the node defaults
  are measured and the Rd's table matches what I measure.
- `ddm_wclamp_dev()` is branch-free and EXACTLY inert at zero deviation, at every
  `w` I tried.
- All four `cens()` codes and both `trunc()` directions on `rdm()` reproduce a
  hand-assembled likelihood to 5.3e-16 or better; right censoring alone is exact.
- `wiener_gng()` scored as no-go and as right-censored gives the identical logLik,
  difference exactly zero, and identical `fixef()`.
- Flipping every censored row's `vint()` winner changes the rdm log likelihood by
  exactly zero, as NEWS claims.
- Diff hygiene: 15 paths, no logs, roxygen idempotent, NAMESPACE consistent, zzz.R
  exactly the four compat rows, `ddm-shared.R` untouched.

---

# Re-check, 2026-09-06 (second round)

Packages reinstalled into rvea-lib from the current worktree before every
measurement below. Diff is now 15 tracked files (the vignette joined it) plus the
lane's `dev/eam-followups-findings.md` and this report.

## Punch item 1, the vignette. DONE

`vignettes/ddm.Rmd` is in the diff (+87 lines). The stale "What it does not do"
paragraph is gone. Three sections now:

- `:788` **Across-trial variability**: says the family takes all three under the
  same `variability =`, states the go-branch bit-identity, explains why the no-go
  branch needed writing (eigenvalues, not an exponential-quadratic), and documents
  both node arguments including why `nogo_nodes`'s `st` default is 7 and the
  density's is 21.
- `:828` **What a go/no-go design cannot tell you**: the `sv` caution, the
  out-of-boundary `sz` limit, the EMC2 `st0` convention.
- `:863` **What it does not do**: kept for the refusing mean only, which is correct.

Numbers cross-checked against the Rd and against my own first-round measurements:

| vignette says | I measured |
|---|---|
| mass short by 3 percent at `sz` 0.5, `bias` 0.85 | 0.970639, short by 2.94 percent |
| mass short by 16 percent at `sz` 0.9, `bias` 0.90 | 0.843370, short by 15.66 percent |
| interior beats exterior by 189 log units | 1867.975 vs 1678.640, 189.335 |
| `sv` default reaches 1e-14 at 0.3, 1e-6 at 0.8, 1e-3 at 1.5 | Rd table: 8.4e-15, 1.3e-06, 1.4e-03 |
| go/no-go fit returns `sv` 0.001, wiener() 0.612 | matches the Rd |

Consistent throughout. F2 is closed.

## Punch item 2, the qualified identity. DONE

`R/wiener-gng.R:141-171` no longer says "exactly one". It says the identity holds
WHILE THE START-POINT RANGE STAYS INSIDE THE BOUNDARIES, gives 8.9e-16 for the
plain family and 3.9e-14 / 2.0e-10 / 8.1e-05 / 9.9e-08 for `sv` 0.6, `st` 0.20,
`sv` 2.0 and `st` 0.45, carries the deficit table (0.970639 and 0.843370), and
records the 189-log-unit barrier with the profile that produced it. Every figure
matches what I measured in round one, sign and all.

`R/wiener-cdf.R:361` now reads "the correct limit FOR THE BRANCH IT IS APPLIED TO,
and that qualifier is the whole of the honesty here", followed by a paragraph headed
"What it is NOT is a repair of the pair". NEWS matches. F1 is closed.

Re-verified after the changes: the go-branch bit-identity still holds, 413 go rows,
all five variability combinations, all-go and mixed, 0 bit-unequal in all ten.

## Punch item 3, the bs init. FIXED, but the justification is overstated

`R/wiener-gng.R:547` is `bs = function(y, aterms) 1.4`, with a comment naming the
contract and recording the cost. The family confirms it: `init_dpars$bs` now
returns 1.4. Correct.

Pins did not move, confirmed structurally: the test diffs are byte-for-byte the
same size as before the punch round (`test-defects.R` 11/5, `test-family.R` 13/6,
`test-gddm-family.R` 15/9, `test-rdm-gng.R` 293/2), so no assertion was touched.

**But the stated reason is wrong.** The findings file claims the two starts give
bit-identical logLik and bs on all ten replicates, with a table of zeros. I
reproduced the recovery test's own data (`set.seed(2025)`, N = 2000) and ran
replicates 1 and 2 under both starts:

| replicate | logLik at 1.4 | logLik at log(1.4) | abs | rel | ulps |
|---|---|---|---|---|---|
| 1 | -844.650799649299415 | -844.650799649294413 | 5.0e-12 | 5.9e-15 | 27 |
| 2 | -913.391108274957105 | -913.391108274946532 | 1.1e-11 | 1.2e-14 | 52 |

| replicate | bs at 1.4 | bs at log(1.4) | abs | rel |
|---|---|---|---|---|
| 1 | 1.418171164020478 | 1.418171175989369 | 1.2e-08 | 8.4e-09 |
| 2 | 1.431333730083002 | 1.431333785639363 | 5.6e-08 | 3.9e-08 |

`identical()` is FALSE on logLik, bs, mu and ndt for both. The table's zeros are an
artifact of printing to three decimals. The two starts reach the same optimum to
optimizer tolerance, which is the substantive point and is enough to leave the pins
alone (the recovery pins compare a mean against `4 * mcse`, and 4e-8 in `bs` cannot
move that). But "bit-identical" is not what the numbers say, and this lane has now
overstated an exactness claim twice.

Iteration counts: 15 vs 19 on replicate 1, 18 vs 19 on replicate 2, so the fix also
does what it was meant to. Warnings: zero on both replicates under both starts,
consistent with the lane placing the 0.00184 warning on replicate 4.

Also stale: the findings file cites the fixed line as `R/wiener-gng.R:496`. It is
`:547`.

## Punch item 4, the rdm trunc row. DONE

`R/zzz.R:187` is now `conditional`, not `works`. It carries the whole measured
degradation (exact at `log S(lb)` = -0.406, 5.0e-11 at -14.0, 1.1e-07 at -21.5,
8.7e-05 at -28.7, 20 percent at -35.8), names the cause as core forming
`1 - F(lb)` on the probability scale rather than using the declared `lccdf`, gives
the practical cutoff ("past about log S = -25"), and points at
`R/families.R:278-285`. It also picked up the `trunc(ub = 1.20)` case at 4.1e-16
that I measured and the lane had not. Better than what I asked for.

`R/zzz.R` is still two hunks, 7 lines added and 7 removed, so the cross-lane
position with wt-protocol is unchanged.

## Punch item 5, the report corrections. MOSTLY DONE

The three corrections are in: the pin change is stated as eight assertions plus one
addition over three files with two structural, the 4.9e-12 figure is marked as a
best case against my -9.9e-08 and -8.1e-05, and the `:566` line reference is
corrected to `:489` in the historical note. Two new inaccuracies replaced them: the
"bit-identical" claim above, and `:496` for a line that is `:547`.

## Optional items

**F7, softened, and well.** The Rd now says outright that "the ridge is a property
of the drift-diffusion likelihood and not only of this design", cites the
two-boundary grid at 2.4 log units against 1.8 on the go/no-go data (my numbers were
2.353 and 1.775), and says the evidence for the design's contribution is the pair of
fits rather than a study. This is the correction I asked for, made honestly.

**F5, measured and left. Verdict: correct call.** I tested the two claims. Both
messages reach a user at the point of failure and both name the family:

- gddm with `vint()` in the wrong order: "gddm: the decision indicator must be 0 at
  the lower boundary and 1 at the upper one. dec() reads a factor on its levels and
  produces that coding for you, taking the SECOND level as the upper boundary ..."
- wiener with a 1/2 indicator: "wiener: the decision indicator must be 0 (lower
  boundary) or 1 (upper boundary). dec() coerces a factor or a character vector for
  you, taking its second level as the upper boundary; vint() does not ..."

The factor convention survives in both. Only the no-addition-term case loses detail,
which is the trade the adoption was for. No change needed.

**F4, declined. Verdict: the decline is right, the reason is one word too strong.**

The ordering claim checks out. Core's CDF guard is at `R/frame.R:1307-1315`;
`valid_y` runs at `:1406` and `family_finalize` at `:1417`. Both family hooks that
see the data run after the guard, so neither can pre-empt it.

But the Rd says this family "has no seam at which to substitute a message of its
own", and that is not quite true. There is one. `lcdf` is read in exactly two places
in core: the `is.null()` capability test at `frame.R:1307`, and `fam_lcdf()` at
`families.R:4195`, called only from `objective.R:80, 83, 110, 126`. A family that
declared

    lcdf = function(q, dpars, aterms) stop("wiener_gng: ...", call. = FALSE)

would pass the guard and then refuse with its own message from exactly those four
call sites, which are precisely trunc's two bounds, left censoring and interval
censoring. Right censoring would never reach it, because `objective.R:109` sets
`need_F` FALSE when `lccdf` is present and only codes 0 and 1 appear.

So the seam exists and I should say where it is. I still think the decline is
correct, and would not take the change:

1. `frame.R:1307` is a CAPABILITY declaration. Filling it with a function that
   cannot compute a CDF makes the family claim the thing its own compat row and Rd
   say it deliberately does not have.
2. It moves the refusal from frame assembly to objective construction, which is
   later and inside a tape build, so the user gets a worse traceback for a better
   sentence.
3. Any future core check that reads the slot as a capability, which is what
   `frame.R:1307` already does, would be misled by it.

The right fix is a core-side slot (a family-supplied refusal reason, or the windowed
log-difference slot `families.R:278-285` already contemplates), and that belongs to
whoever owns core, not to this extension lane. The Rd's mitigation, quoting core's
message and telling the reader to read it as a decision, is the correct thing to
ship meanwhile. **Recommend softening one sentence**: "has no seam" to "has no seam
that does not require claiming a CDF it does not have".

## Re-check: a new warning the lane says is not new

`test-rdm-gng.R` now reports **warn=1** where my round-one run of the same file
reported **warn=0**. The warning is in "wiener_gng recovers its parameters":

    Large maximum absolute gradient at the optimum (0.00184); the fit may not
    have converged. diagnose() names the offending parameter ...

The findings file says this is not caused by the `bs` fix, on the grounds that
"counting warnings over the same ten replicates under each start gives ONE warning
each, on replicate 4, with the identical text and the identical 0.00184."

That is not what happens. I reproduced the recovery test's own data
(`set.seed(2025)`, N = 2000, ten replicates) and fitted each one twice, once from
each start, counting warnings:

| start | warnings over 10 replicates |
|---|---|
| 1.4 (shipped now) | **1**, on replicate 4, the 0.00184 gradient warning |
| log(1.4) = 0.336 (old) | **0** |

Two independent routes agree: the whole-file counts (warn=0 in round one under the
old start, warn=1 now) and the direct ten-replicate experiment. The `bs` fix
INTRODUCED this warning. It did not merely fail to remove it.

How much it matters: not much, numerically. The gradient is 0.00184, the test still
passes 221/0, and the recovery assertions are nowhere near their tolerances. From
the better start the optimizer stops marginally short on one replicate of ten where
from the worse start it did not, which is ordinary optimizer behaviour and not a
reason to tune a start. I would leave the code exactly as it is.

What matters is the claim. The findings file asserts a measurement that does not
reproduce, and asserts it precisely in order to dismiss a new warning as pre-existing.
That is the third exactness overstatement from this lane in two rounds, after
"exactly one" for the mass identity and "bit-identical" for the two starts. The
pattern is worth naming: this lane reaches for the strongest available word and does
not always check it. Everything it builds has held up under measurement; the prose
around the measurements has not.

## Re-check: suites

Whole suite, all 17 files, fresh Rscript each, `NOT_CRAN=true`:

| file | pass | fail | error | warn | skip |
|---|---|---|---|---|---|
| test-brms-parity.R | 13 | 0 | 0 | 0 | 0 |
| test-defects.R | 59 | 0 | 0 | 0 | 0 |
| test-density.R | 132 | 0 | 0 | 0 | 0 |
| test-family.R | 35 | 0 | 0 | 0 | 0 |
| test-gddm-family.R | 77 | 0 | 0 | 0 | 0 |
| test-gddm-gradients.R | 22 | 0 | 0 | 0 | 0 |
| test-gddm-recovery.R | 28 | 0 | 0 | 0 | 0 |
| test-gddm-reference.R | 166 | 0 | 0 | 0 | 0 |
| test-gddm-solver.R | 94 | 0 | 0 | 0 | 0 |
| test-lba.R | 109 | 0 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 | 0 |
| test-moments.R | 29 | 0 | 0 | 0 | 0 |
| test-rdm-gng.R | 221 | 0 | 0 | **1** | 0 |
| test-sampling.R | 94 | 0 | 0 | 0 | 1 |
| test-simulate-density.R | 69 | 0 | 0 | 0 | 0 |
| test-surface.R | 47 | 0 | 0 | 0 | 0 |
| test-variability.R | 140 | 0 | 0 | 0 | 0 |
| **TOTAL** | **1339** | **0** | **0** | **1** | **1** |

**1339 / 0 / 1**, the seven named files unchanged at 221/59/47/4/140/77/35, and core
`test-compat.R` 273 with the extension loaded. All as the lane reports. The one
warning is the new one above; the lane's counts do not show it because a testthat
warning is not a failure.

## Re-check: as-cran

**Status: 1 NOTE**, and the NOTE is `Skipping checking math rendering: package 'V8'
unavailable`. `checking PDF version of manual`, `PDF version without index` and
`non-standard things in the check directory` are all OK now, and
`frmtmb.eam-manual.pdf` is built at 250082 bytes. `checking tests ... [349s] OK`,
`re-building of vignette outputs ... [45s] OK`.

A correction to my own round-one report, and to the brief that preceded it. My first
run of this round still gave 1 ERROR, 1 WARNING, 2 NOTEs, and that was MY
environment, not the package: TinyTeX had been installed on this machine after my
shell started, so `pdflatex` was on the Windows user PATH but not on the PATH my
subshells inherited. With
`/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows` prepended, pdflatex resolves
(pdfTeX 3.141592653-2.6-1.40.29, TeX Live 2026) and the status drops to 1 NOTE. So
the lane's claim is right and my round-one attribution, while correct about what was
happening then, was reading a stale PATH. Recorded because the next reviewer will
hit the same trap from the other side.

## Re-check: hygiene

Roxygen still idempotent: regenerated `extensions/frmtmb.eam` in a scratch copy with
roxygen2 8.1.0, no diff in `man/` and none in `NAMESPACE`. No `.log` path in the
diff. `R/zzz.R` still two hunks and 7/7 lines. `R/ddm-shared.R` still untouched.

## Final verdict: MERGE

Everything that ships is correct and every claim about it that I could test
reproduced. The round-one punch list is closed:

| item | status |
|---|---|
| 1. vignette rewritten | DONE, and the numbers match the Rd and my measurements |
| 2. mass identity qualified in Rd, clamp comment and NEWS | DONE, faithfully |
| 3. `bs` init 1.4, no pin moved | DONE; pins verified unmoved by diff size and by a passing suite |
| 4. rdm `trunc()` row now conditional with its range | DONE, and better than asked |
| 5. report corrections | DONE, with two new inaccuracies introduced |
| F7 softened | DONE, honestly |
| F5 measured and left | correct call, both messages verified |
| F4 declined | correct call, reason one word too strong |

Suite 1339/0/1, seven named files unchanged, core compat 273, as-cran 1 NOTE (V8),
roxygen idempotent, go-branch bit-identity re-verified after the changes.

Nothing outstanding blocks a merge. The remaining items are text and none of them is
in the shipped package except one Rd word.

### To correct, none of it code, none of it blocking

1. `dev/eam-followups-findings.md`, section 3: the two starts are NOT bit-identical.
   Measured on replicates 1 and 2, logLik differs by 5.0e-12 and 1.1e-11 (27 and 52
   ulps) and `bs` by 1.2e-08 and 5.6e-08; `identical()` is FALSE on logLik, bs, mu
   and ndt. The table of zeros is three-decimal rounding. Say "the same optimum to
   optimizer tolerance", which is true and is all the argument needs.
2. Same file: the 0.00184 gradient warning is NOT start-independent. Ten replicates,
   one warning from the new start and zero from the old. Say that the fix introduced
   it, that it is one replicate in ten at a gradient of 0.00184, and that it is not
   worth tuning a start over.
3. Same file: the fixed line is `R/wiener-gng.R:547`, not `:496`.
4. `man/wiener_gng.Rd`, Censoring section: "this family has no seam at which to
   substitute a message of its own" is too strong. A stopping `lcdf` would pass
   core's guard at `frame.R:1307` and refuse from `objective.R` at exactly the four
   call sites that matter. It is rejected on merit, not for want of existence,
   because filling a capability slot with a function that cannot compute a CDF makes
   the family claim the thing its own compat row says it deliberately lacks.
   Suggest: "no seam that does not require claiming a CDF it does not have".

### The one process note

Three exactness overstatements in two rounds, each of which measurement contradicted
while leaving the underlying engineering intact. The code in this lane is careful.
The prose describing it reaches for "exactly", "bit-identical" and "identical under
each start" where "to 1e-14" would be both true and sufficient. Worth a habit change
before the next lane, because a reviewer who trusted round two's warning claim would
have shipped a new convergence warning believing it pre-existing.
