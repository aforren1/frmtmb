# Lane fitedge: fitting edge cases and the SE check's cost

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-fitedge`, library `wt-fitedge-lib`, "before" arm `rellib-r18`.
Scripts are `dev/fitedge-*.R`; logs are in `dev/fitedge-out/` (gitignored).
Every number below is pasted from a log named beside it. Status: items 1 to
4 done, item 5 not done (see its section).

## Fix round (review dev/reviews/2026-11-06-fitedge.md: FIX)

This section holds where it disagrees with the sections below it, which are
the first submission's record. Logs: `dev/fitedge-out/r2/`.

### B1: item 1 made fits worse than 0.79.1

Two changes.

1. Default start: the fit also descends from 0.79.1's start (the objective
   scan with the body's blocks not covered) and keeps the better optimum.
   This is lane arstep's `gp_second_descent()`, copied verbatim (its text and
   roxygen are identical to arstep's, `diff` empty), called at the same place
   in `fit_assembled()` with the same two lines; `gp_start_scan()` puts the
   second start in `obj$env$.frmtmb_gp_alt` after saving and restoring the
   inner state, as arstep's does. The lost descent is in `opt$evals`
   (`tally$evals`).
2. A caller's `start`: no body start at all (`gp_start_working(body =
   FALSE)`), so the fit is 0.79.1's. A restart from a saddle inside
   `nl_fit_off_saddle()` is not the caller's (`announce_start = FALSE`) and
   keeps the body start; that is the review's case M, -4.40 against 0.79.1's
   -107.94 (`gnprobe-*`).

While pinning the review's surviving `gn1` mutant I measured one Gauss-Newton
step against two with the second descent in place
(`dev/fitedge-r2-gn1rows.R`, `gn1rows.txt`, 48 rows: the lane96 `exp(a)`
body and the fresh `a * exp(b z)` and logistic bodies): one step was never
lower than two, and higher on 4 rows by 4.49, 44.10, 2.96 and 14.72 nats (E
wig 60 2, E wig 120 1, F wig 120 1, F wig 120 2). The first submission's
"the second step is what `exp(a)` needs" was never measured on its own: the
17-to-37-nat rows were a start from `y` alone, one step and no other
coefficient. So `gp_body_start()` now takes one step.

The claim holds for the pair, the one-step body start plus the second
descent, and not for one step alone. The review's short check (arm
`nosec`, `gp_second_descent()` returning the first descent;
`gpsum-fresh-nosec.txt`) compared one step alone with the reviewed two-step
build, which had no second descent: one step higher on 4 rows (by 4.494,
44.100, 2.964, 14.723) and lower on 2 (F wig 60 1 by 9.086, F wig 60 2 by
0.140); alone, one step is worse than base on 6 of 96 rows (worst 9.086),
where two steps alone were worse on 4 (worst 3.30). Shipped is the better
of two descents: 2 short rows on that set, against 6 for the two-step pair.

The two-step build's
sweep is kept in `r2/twostep/` (fresh 6 short, worst -44.10; lane96 7
short; 0 worse than base on either).

The review's 192 rows (`dev/fitedge-r2-gpfresh.R`, the review's
`gpfresh.R` with its output moved to `r2/`; base and reference are the
review's and the first submission's logs; `dev/fitedge-r2-gpsum.R`,
`gpsum-r2.txt`):

    R2 fresh  all rows 96 | short base 26 (silent 26) reviewed 10 fix 2 (silent 2)
      | worst base -270.322 fix -36.500 | fix worse than base 0, better 24
    R2 lane96 all rows 96 | short base 30 (silent 30) reviewed 10 fix 7 (silent 7)
      | worst base -290.420 fix -0.635 | fix worse than base 0, better 25

0 of 192 rows below rellib-r18. The 9 left: E wig 120 2 (-36.50) and L wig
120 2 (-5.24), both equal to base; 5 low-count poisson rows equal to the
linear spelling; C step 60 1 (-0.130) and D step 60 1 (-0.002), equal to
base. Every one is silent, as on base.

The review's `start` cases (`dev/review-fitedge-gnprobe.R lane`,
`gnprobe-1step.txt`): every case now reads 0.79.1's value or higher.
Against the review's `gnprobe-base.txt` and `gnprobe-lane.txt`:

    case                      rellib-r18   reviewed     fix
    E a * exp(b z), default     8.2413      28.1839    28.1839
    E, start b = 3             28.1839       8.2413    28.1839
    E, start b = -3             8.2413     -17.6974     8.2413
    E, start a = 50             8.2413       8.2413     8.2413
    M b * a, default         -107.9420      -4.4031    -4.4031
    M, start b = 1           -107.9420      -4.4031  -107.9420
    F logistic, start c = 8    39.9185     -85.2537    39.9185
    F, start b = 20            39.9185      37.7627    39.9185
    E, y * 5000, default     -491.6328    -493.3618  -491.6328

M with `start` b = 1 reads 0.79.1's value: the caller's start turns the
body start off, which is the rule asked for; the reviewed build read
-4.4031 there.

Cost, the shipped build against 0.79.1's path in the same build
(`gp_body_start()` patched to NULL, which also removes the second start),
interleaved, order alternating (`dev/fitedge-r2-cost.R`, `cost.txt`):

    COST2 fresh rows 96: shipped 183.6 s, 0.79.1 path 90.9 s, ratio 2.02;
      shipped above 24, below 0
    COST2 lane96 rows 96: shipped 162.6 s, 0.79.1 path 73.3 s, ratio 2.22;
      shipped above 25, below 0

The two-step build with the second descent read 2.02 and 2.42
(`r2/twostep/cost.txt`).

### M1: multiplicative bodies

Filed in `dev/test-backlog.md` ("Filed in the round opened 2026-11-06",
fitedge) with the review's numbers. NEWS now says the fix is for a body
additive in the field's parameter and puts the multiplicative case in a
Limits section. The one-step change above reaches 4 of the 6 multiplicative
rows the review listed; E wig 120 2 (-36.50) and L wig 120 2 (-5.24) stay.

### M2: merge with arstep (exact edits for consolidation)

`diff3 -m -E` of this tree's `R/gp-start.R` and arstep's over 3cd2a0c8
(`r2/merge/merged_gp-start_E.R`): 5 conflicts. `R/fit.R` merges clean (both
add the same two lines). The resolutions:

1. `gp_marginal_data()`, after `u <- max(J)`: keep fitedge's
   `if (u > max_rows) return(NULL)`, then arstep's two lines
   `ard <- if (!is.null(ar)) gp_ar_data(ar, J, u, y, X)` and
   `if (!is.null(ar) && is.null(ard)) return(NULL)`.
2. `gp_start_working()` roxygen: keep fitedge's two paragraphs (the body
   start, `tpl` and `body`), and one `@return`:
   "`list(idx, theta, alts, body)`: the positions in the template's `theta`,
   their values, the full `theta` vectors at each other local maximum the
   working marginal reached (gp_marginal_alts()), and which of the positions
   the body start set, or `NULL`."
3. Its locals: both, `alts <- list()` and `from_body <- integer(0)`.
4. Its end: `if (length(idx)) list(idx = idx, theta = theta, alts = alts,
   body = from_body)`, then fitedge's `gp_body_start()` whole (diff3 puts all
   of it inside this conflict).
5. `gp_start_scan()`'s tail, then one copy of `gp_second_descent()` (the two
   are identical) and arstep's `gp_marginal_agrees()`. The tail:

        out <- gp_start_scan_units(obj, template, frame, val, covered)
        body <- intersect(covered, w[["body"]])
        if (!length(covered) || (!length(w$alts) && !length(body))) {
          return(out)
        }
        saved <- obj_state_save(obj)
        on.exit(obj_state_restore(obj, saved), add = TRUE)
        alt <- list()
        if (length(w$alts)) {
          tpl <- template
          tpl[["theta"]][w$idx] <- w$alts[[1L]]
          if (!is.na(val(tpl))) {
            alt <- c(alt, list(outer_from_template(tpl, obj, frame)))
          }
        }
        if (length(body)) {
          keep <- setdiff(covered, body)
          scanned[["theta"]][keep] <- template[["theta"]][keep]
          p <- gp_start_scan_units(obj, scanned, frame, val, keep)$par
          if (length(p) == length(out$par) && all(is.finite(p))) {
            alt <- c(alt, list(p))
          }
        }
        if (length(alt)) out$env$.frmtmb_gp_alt <- alt
        out
      }

   (`scanned <- template` must stay where fitedge sets it, before the working
   start overwrites `template`.)

The silent one, which no marker shows: in `gp_start_working()`'s non-`mu`
branch, add `alts <- lapply(alts, function(th) c(th, m$theta))` before
`idx <- c(idx, from)`, as arstep's section 17 says; without it a `mu` block
before a body block leaves `w$alts[[1]]` shorter than `w$idx`, and the
assignment recycles it.

Also after the merge:

- `frm()`'s roxygen (`R/fit.R`, `@param start`, arstep's paragraph) says "A
  nonlinear body, a mixture and `sigma ~ gp()` start at the best of a few
  length scales". Replace that sentence with: "A field in a parameter of a
  nonlinear body starts at the marginal of the body made linear at its
  start, and the fit also descends from the best of a few length scales and
  keeps the better; a `start` of your own turns the first off. A mixture and
  `sigma ~ gp()` start at the best of a few length scales read from the
  objective." I did not edit it here: arstep rewrote the same lines.
- `gp_marginal_outer_gr()` reads `pn %in% c("beta", "betad", "theta")`, so a
  fit with arstep's `thetaac` keeps the z tape; arstep's `sigc` (a constant
  sigma) passes the structural gate and is then decided by the value check.
  Rerun item 3's test and `dev/review-fitedge-se.R fall` on the merged tree.

### Minors

- m1, mutants (`dev/fitedge-mutrun.R`, the review's runner plus `prefix`
  and `gn2`; `r2/mut-*`, `r2/mut1-*`): 

        prefix (the reviewed build)   test-gp-start: 2 fail (-3.3 on K,
                                      123.4 on the c = 8 start)
        gn2 (two steps, reviewed)     test-gp-start: 1 fail (29.3 against 73.4)
        noothers                      test-gp-start: 2 fail (0.441; 12.6 against 29.9)
        ridgekept                     test-nl-family: 1 fail (gradient warning)
        control                       test-gp-start: pass 233, before the one-step test

  `gn1` has no target now that the code takes one step; `gn2` is its
  inverse.
- m2, the ridge pass now runs only on "Optimizer did not report
  convergence" or "Large maximum absolute gradient" (`nl_stop_warning()`).
  In one process, interleaved against the pass on any held warning
  (`dev/fitedge-r2-ridge.R`, `ridge.txt`):

        RGSUM mvs rows 8: stop-keyed 4.9 s, any-warning 6.9 s, ratio 0.72;
          logLik equal 8 of 8, warnings equal 8 of 8
        RGSUM mv4 rows 8: stop-keyed 33.5 s, any-warning 32.8 s, ratio 1.02;
          logLik equal 8 of 8, warnings equal 8 of 8

  The `mv4` values equal the reviewed build's on all 8 rows (-624.650299 to
  -639.088600), each 0.50 to 5.55 above rellib-r18's
  (`dev/review-fitedge-out/nl-mv4-fit-base.tsv`). The item 2 test now reads
  the warning set, and catches `ridgekept`.
- m3, a read `mi()` predictor beside a latent residual,
  `ln ~ mi(xm) + ar(t, g, cov = TRUE)` with `xm | mi() ~ 1`, measured
  against the exact marginal (lognormal and a gaussian `xm` are jointly
  gaussian, so the missing `xm` integrate out in closed form; no frmtmb
  code; `dev/fitedge-mi-pred.R`, `r2/mi-pred-base.txt`, seeds 73001 to
  73005, 10 groups of 6, 6 missing):

        MP seed 73001 frm ll -184.72631086 | exact max -184.72631086 | diff -8.6e-10
        MP seed 73002 frm ll -198.60148558 | exact max -198.60148556 | diff -1.5e-08
        MP seed 73003 frm ll -184.53454017 | exact max -184.53454026 | diff  9.4e-08
        MP seed 73004 frm ll -198.25386894 | exact max -198.25386894 | diff -9.7e-10
        MP seed 73005 frm ll -159.35077898 | exact max -159.35077943 | diff  4.5e-07

  (positive: frm() reached a hair higher than my optimizer, where sigma ran
  to 0.) The estimates agree to 6 digits. So the shared inner problem is
  right where it can be checked, and the refusal's reason is now true: it
  names what is still unmeasured, a read value of the latent-residual
  response itself, which sits in its own density beside its row's latent,
  and says the predictor case has been measured. Non-gaussian families
  (Gamma, weibull) have no closed form; I did not measure them.
- Filed in `dev/test-backlog.md`: the `frmtmb_latres_laplace` rounding
  pairs, the lognormal false alarm, the skew_normal three-body rows.
- NEWS: bullet 2 rephrased, the `mi()` bullet BREAKING, whole-fit times
  beside the Hessian's (`dev/fitedge-fit-time.R lane 500 2`,
  `r2/ftime-500.txt`: `frm()` + `summary()` 34.20 s against 53.00 s, ratio
  1.55, logLik bitwise equal, SEs to 5.4e-10), a Limits section.

### Tests on the fix build (one file per process, gated tier, `r2/t-*.txt`)

    test-gp-start.R        pass=234 fail=0 err=0 skip=0 warn=0
    test-latres.R          pass=301 fail=0 err=0 skip=0 warn=0
    test-nl-family.R       pass=322 fail=0 err=0 skip=0 warn=0
    test-se-check.R        pass=134 fail=0 err=0 skip=0 warn=0
    test-gp-noncentered.R  pass=76  fail=0 err=0 skip=0 warn=0

Core suite, one file per process, gated tier, the fix build on all 308 `lib:`
lines (`dev/fitedge-run-par.sh fitedge-core-r2 dev/fitedge-core.jobs 2`,
`fitedge-core-r2.log`, `r2/suite-sum.txt`):

    frmtmb 308 files pass 25580 fail 0 err 0 skip 15 warn 0
    FITEDGE-CORE-R2 ran 308 of 308

The skips are other gates (drmTMB, fuzz, scale).

`R CMD check --as-cran` of core on the fix tree (`dev/fitedge-check.sh`,
log `r2/check-core.log`, check tree deleted): `Status: 1 NOTE`, the V8
math-rendering one; tests 57 min, OK.

New or changed: in test-gp-start.R "a field in a body never ends below the
start it replaced" (the review's K row and its `start` c = 8 row, each
against the same build with the body start mocked off) and "the body start
reads every coefficient, in one step"; in test-nl-family.R the item 2 test
reads the warning set. Each fails under its mutant above.

### After the short check (SHIP, three minors)

- m1: `fit_assembled()` takes `gp_body = is.null(start)`, and passes it to
  `gp_start_scan()` and to its autoscale refit; `frm()`'s `run()` passes
  `gp_body = is.null(start)` of the caller's own `start`, so a saddle
  restart keeps the body start only when the caller gave none. So
  `frm_allfit()` (`start = fit$start`), the autoscale pre-fit and
  `refit(newresp =)` skip the body start whenever the caller gave a
  `start`, as NEWS says. Pinned in test-gp-start.R ("a caller's start turns
  the body start off on every refit path", through `frm_allfit()`, with the
  absent case: no `start`, and the body start runs). Mutant `allstart`
  (`gp_body = TRUE` as the default, the old rule on those paths):
  `r2/mut2-allstart.txt`, 1 fail (`n_body` 1, not 0).
- m2: the refusal now says "(a missing predictor read through mi() has
  been, on lognormal)". New test in test-latres.R, "a read mi() predictor
  beside a latent residual is exact": `-fit$obj$fn(par)` at a point of the
  test's choosing against the dense marginal with the missing `xm`
  integrated out (seed 73001, 6 groups of 4, 3 missing), to 1e-8 relative.
  It pins behavior both builds share, so it passes on rellib-r18 too.
- m3: `check_convergence()`'s three stop warnings carry classes:
  `frmtmb_conv_code` ("Optimizer did not report convergence"),
  `frmtmb_conv_gradient` ("Large maximum absolute gradient"),
  `frmtmb_conv_gradient_nonfinite` ("The gradient at the reported optimum
  is not finite"); the messages are unchanged, and `warnings` in its
  returned list stays unnamed. `nl_stop_warning()` reads those classes and
  no text. Test in test-nl-family.R ("the ridge pass reads the stop warnings
  by class"): each class turns the pass on, the non-finite one included, an
  SE warning and the old text with no class do not, and a one-iteration
  `y ~ x + (1 | g)` fit raises `frmtmb_conv_code` and
  `frmtmb_conv_gradient`. Mutant `textmatch` (the text match):
  `r2/mut2-textmatch.txt`, 4 fails.
- The one-step claim is restated for the pair (item B1 above).

Tests (one file per process, gated tier, `r2/tests-c.txt` to
`tests-f.txt`):

    test-gp-start.R           pass=238 fail=0 err=0 skip=0 warn=0
    test-nl-family.R          pass=329 fail=0 err=0 skip=0 warn=0
    test-latres.R             pass=303 fail=0 err=0 skip=0 warn=0
    test-allfit-start.R       pass=14  fail=0 err=0 skip=0 warn=0
    test-autoscale.R          pass=50  fail=0 err=0 skip=0 warn=0
    test-conditions-census.R  pass=11  fail=0 err=0 skip=0 warn=0

(test-nl-family.R first read 327 pass and 1 fail: a one-iteration
`y ~ x` GLM raised no warning on that seed, so the class check moved to the
mixed model above.) No core suite or R CMD check after these edits.

### Not done in the fix round

- No OpenBLAS pass, and no extension suite or session rerun: the coordinator
  asked for core, and the fix round's paths are core's (the body start, the
  second descent, the ridge pass's trigger, a message).
- No merged-tree build with arstep; the M2 edits are written, not run.
- The m3 measurement has no test: the reference is a dense density that
  would need the fit's outer layout spelled by hand.

## Shared functions this lane changes (for arstep and nlstart)

- `R/gp-start.R`: `gp_start_working()` (restructured: the block renumbering
  moved above the working response, and a branch for a block in a nonlinear
  parameter), new `gp_body_start()`, new `gp_marginal_outer_gr()`,
  `gp_marginal_value_tol` and `gp_marginal_kappa_max`;
  `gp_marginal_setup()` and `gp_marginal_data()` take
  a trailing `max_rows` (default unchanged); `gp_start_scan()` passes
  `template` to `gp_start_working()`. Fix round: `gp_start_working(body =)`
  and its `body` element, `gp_start_scan(body =)` and its second start in
  `obj$env$.frmtmb_gp_alt`, `gp_second_descent()` (arstep's, verbatim).
- `R/fit.R` (fix round): `fit_assembled()` gets `gp_body = is.null(start)`
  (passed to `gp_start_scan()` and its autoscale refit; `frm()`'s `run()`
  passes the caller's) and calls `gp_second_descent()` after the first
  descent (arstep's two lines); `check_convergence()` gives its three stop
  warnings classes.
- `R/se-check.R`: `fit_outer_hessian()` (the gradient the differences read).
- `R/par-template.R`: `nl_fit_off_saddle()` (a ridge pass), new
  `nl_all_bodies()`; fix round: new `nl_stop_warning()`.
- `R/autocor-latent.R`: `latres_refuse_laplace()` (the `mi()` clause).
- `R/compat.R`: the notes (not the status) of the rows
  (`group:autocor`, `kind:aterm`) and (`ar()`, `ma()`, `arma()`,
  `kind:aterm`), which said `mi()` on the response is refused beside a latent
  residual; they now say where (`dev/fitedge-compat-edit.R`).

## Item 1: gp() in a nonlinear body

Reproduced on rellib-r18 (`dev/fitedge-gpbody-probe.R`, log
`probe1-base.txt`): `bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)` against
`y ~ z + gp(x)`: seed 1811 `2 sin(3x) + x` 1.6536 linear, -111.8512 body
(gap -113.505); seed 58001 (the ncgp review's edge2 row) -10.1336 against
-110.5492 (gap -100.416). Conv 0, no warning on both. The review's `x18.R` is
gone (its `review-rel079-x-run/` directory was gitignored and deleted), so
the 1.12-nat `sin(x)` row itself is not reconstructed; seed 1811 `sin(x)` in
my design reads gap 0.

Cause: `gp_marginal_gate()` refuses a frame with `nlpars`, and
`gp_start_working()` skipped every predictor but `mu`. A block in `a` got
only `gp_start_scan_units()`. A theta start alone carries the body to the
linear spelling's maximum (`dev/fitedge-gpbody-proto.R`, `proto1-base.txt`:
both rows reach 1.6536 and -10.1336 from the theta of `y ~ gp(x)` on y).

Fix: `gp_body_start()`, the body linearized at the template, handed to
`gp_marginal_opt()`. `y_w - eta0 + J_a a0 + sum(J_k b_k0)` is linear in the
coefficients of `a`, in every other coefficient the body reads and in the
field (rows scaled by `J_a`). On `a + b * z` this is exactly the linear
spelling's marginal. A second Gauss-Newton step runs only where the
derivatives move between steps (a body nonlinear in its coefficients). A first
version that read `y` alone (one step, no other coefficient) left `A sin n 60
seed 2` 0.441 short and the `exp(a)` rows 17 to 37 short
(`gpbody-fit-lane-gn2.tsv` is the two-step version without the early stop;
`dev/fitedge-gpbody-diag.R`, logs `gpbody-diag-dev.txt`,
`gpbody-diag-dev2.txt`).

Sweep, `dev/fitedge-gpbody-sweep.R`: 4 bodies x 4 curves (`sin(x)`,
`2 sin(3x) + x`, `exp(0.3x)`, `sign(sin(1.5x))`) x n 60, 120 x seeds 1 to 3,
`set.seed(18000 + 100 seed + n + curve index)`. Bodies: A `a + b * z`
gaussian, B `a + b * z` poisson, C `exp(a) + b * z` gaussian, D `b + a * w`
gaussian, `a ~ gp(x)` throughout. Reference: the best of the linear spelling
(A, B) and 12 theta starts on the body (`gpbody-ref-lane.tsv`; a theta start
skips the warm start, so the reference is build-free). Short: more than 1e-3
nats below. Summary `dev/fitedge-gpbody-sum.R`, log `gpbody-sum3.txt`:

    SUM base rows 96 short 30 silent 30 worst -290.420
    SUM base body A rows 24 short 6 silent 6 worst -249.915
    SUM base body B rows 24 short 9 silent 9 worst -1.096
    SUM base body C rows 24 short 9 silent 9 worst -41.500
    SUM base body D rows 24 short 6 silent 6 worst -290.420
    SUM lane rows 96 short 10 silent 10 worst -0.635
    SUM lane body A rows 24 short 0 silent 0 worst 0.000
    SUM lane body B rows 24 short 7 silent 7 worst -0.635
    SUM lane body C rows 24 short 2 silent 2 worst -0.389
    SUM lane body D rows 24 short 1 silent 1 worst -0.002

Row by row and against the linear spelling's own default
(`dev/fitedge-gpbody-cmp.R`, `gpbody-cmp3.txt`; `gpbody-lin-lane.tsv`):

    CMP body A lane better 6 worse 0
    CMP body B lane better 6 worse 2
    CMP body C lane better 8 worse 1
    CMP body D lane better 5 worse 0
    WORSE B step n 60 seed 1 base -124.8471 lane -125.2335 ref -124.8471
    WORSE B wig n 120 seed 3 base -263.4178 lane -263.7027 ref -263.4178
    WORSE C exp n 60 seed 2 base 46.9990 lane 46.6098 ref 46.9990
    LIN body A rows 24 body equals linear (1e-6 rel) 24, above 0, below 0
    LIN body B rows 24 body equals linear (1e-6 rel) 22, above 2, below 0

So the body now reads the linear spelling's value on every A and B row, and
above it on 2. The 7 B rows still short are short on the linear spelling too:
low-count poisson, the working start's limit that 0.79.0's NEWS "Limits"
already names (worst 1.77 there). The two B rows worse than base are rows
where 0.79.1's scan happened to beat the linear spelling. C `exp` n 60 seed 2
is the one row a nonlinear body loses, 0.389. What is left is silent, as the
linear spelling's own shortfall is; no detector exists off the gaussian
(dev/ncgp-findings.md, round 2).

Cost, interleaved in one process, the fix on and off (`gp_body_start()`
patched to return NULL is 0.79.1's start), order alternating by row
(`dev/fitedge-cost12.R`, `cost12.txt`), bodies A and C, 48 rows:

    COST item 1 rows 48: fix on 110.5 s, off 61.8 s, ratio 1.79;
      ll on above off 14, below 1

As wall sums from the sweep's separate processes (load differs):
A 11.9 s to 27.1 s, B 13.5 to 23.4, C 19.1 to 44.0, D 11.7 to 30.2; the
linear spelling on the same rows: A 17.0 s, B 18.1 s. 0.79.1's body was
cheap because it had no warm start; the linear spelling pays the same kind
of start.

## Item 2: a disc saddle beside a mu ridge

Reproduced (`dev/fitedge-disc-probe.R`, `disc-base.txt`): default -377.0349
code 1 with two warnings; `start` of mu only reaches -371.5414.
Cause (`dev/fitedge-disc-starts.R`, `disc-starts-base.txt`): the first fit
stops with disc on its saddle and mu on the linear limit of exp() (b1 -13832,
b2 -2.5e-5). The restart moves disc only (the saddle body); mu restarts from
its zeros and runs to the same limit. Of 16 mu starts with disc moved, every
b1 > 0 start with b2 > 0, and every b1 = 1 start, reaches -371.5414; b1 < 0
does not.

Fix: in `nl_fit_off_saddle()`, a kept restart that still holds a warning gets
one more pass that moves every body not yet moved, from the pattern (not from
the kept values, which leave disc at d2 = 1499 and a gradient warning: that
first version read -371.5414 with "Large maximum absolute gradient (0.23 at
d1)", `disc-lane.txt`). After: -371.5414, code 0, the one warning that d1 and
d2 enter only through their product, which is true (`disc` probe on the dev
arm). A model with one body never reaches the pass: every body is moved
already. So the 450 single-body nl rows of the nlfix review cannot change, by
construction.

Sweep, `dev/fitedge-nl2-sweep.R`, models with two bodies, data as
`dev/review-nlfix-d7.R` (`gen(seed * 100 + 1)`, seeds 1 to 5, slopes +-0.6):
`disc` (cumulative, mu body `exp`, `sat`, `two` beside `nlf(disc ~ exp(d1 *
d2 * z))`), `sigma` (gaussian, the same mu bodies beside `nlf(sigma ~ s0 + s1
* s2 * z)`), `mv` (mvbf of two cumulative `exp` bodies). Reference: the best
of 8 random starts (`set.seed(seed * 7 + 1)`) and both arms. Summary
`dev/fitedge-nl2-sum.R`, log `nl2-sum.txt`:

    NL2 disc base rows 30 short 20 silent 0 worst -5.282 time 24.8 s
    NL2 disc lane rows 30 short 0 silent 0 worst 0.000 time 34.2 s
    NL2 disc lane above base 20, below base 0
    NL2 sigma base rows 30 short 1 silent 0 worst -69.944 time 8.9 s
    NL2 sigma lane rows 30 short 1 silent 0 worst -69.944 time 9.8 s
    NL2 sigma lane above base 0, below base 0
    NL2 mv base rows 10 short 0 silent 0 worst 0.000 time 12.7 s
    NL2 mv lane rows 10 short 0 silent 0 worst 0.000 time 12.0 s

Every one of the 20 disc rows was warned on 0.79.1 (code 1 and two warnings
on 17). The `sigma` row is `b1 (1 - exp(-b2 x))`, 70 nats, lane nlstart's
item 1, warned on both. Times above are wall sums in separate processes.
Interleaved in one process, the ridge pass on and off (`nl_all_bodies()`
patched to return no body), `dev/fitedge-cost12.R`, `cost12.txt`:

    COST item 2 rows 30: fix on 68.4 s, off 48.1 s, ratio 1.42;
      ll on above off 20, below 0

## Item 3: the marginal's Hessian for the SE check

`fit_outer_hessian()` takes the same `fd_hessian()` differences, at the same
steps, of the marginal tape's gradient where `gp_marginal_outer_gr()` admits
the fit. Not the AD Hessian of the marginal: the differences keep `exact =
FALSE` and the noise `E` the SE check's tiers read, so every tier decides as
before.

The gate, in order: an exact `gp()` block; maximum likelihood with no prior,
no REML, quadrature, importance, profile or whitened beta; outer parameters
`beta`, `betad`, `theta` only; `gp_marginal_setup()` admits the frame (with
no row limit: `max_rows = Inf`); the bound on the condition of the
marginal's covariance at most `gp_marginal_kappa_max = 1e6`; the marginal
and the z tape read one value at the point to 1e-9 relative.

The condition gate is measured, not assumed. The first version had the value
check only, and `test-gp-start.R`'s "y's and x's units" failed 4 times on the
lane (`t-gp-start-lane.txt`: fitted SE ratio across scales 1.0143 against the
1e-2 bound). `dev/fitedge-scale-diag.R` (`scale-diag.txt`): at seed 1 the
marginal's fitted SEs were 1.0145 times the z tape's, and the value check
passed at one scale and not at another. Which is right:
`dev/fitedge-hess-mpfr.R` (`hess-mpfr.txt`) takes the exact Hessian of the
marginal written from the model statement in 200-bit Rmpfr arithmetic:

    MP y scale 1:   max rel SE error: z tape 3.134e-05, marginal 3.455e-01
    MP y scale 0.001: max rel SE error: z tape 4.480e-05, marginal 2.296e-01

So at noise 1e-3 the marginal's dense covariance has lost its digits and the
z tape has not. The forced runs (`FE_FORCE=1 dev/fitedge-se-verdict.R`, every
gate but the structural one lifted, both routes in one process; the old route
is 0.79.1's code path) bin the disagreement by the bound
`kappa = 1 + u sd^2 / sigma^2` (`dev/fitedge-sev-scan.R force`,
`sev-scan-force.txt`):

    KAPPA (1,100]        rows  28 SE max rel 4.84e-12 NaN sets differ 0
    KAPPA (100,1e+03]    rows  59 SE max rel 9.04e-11 NaN sets differ 0
    KAPPA (1e+03,1e+04]  rows  47 SE max rel 1.56e-10 NaN sets differ 0
    KAPPA (1e+04,1e+05]  rows  62 SE max rel 2.55e-09 NaN sets differ 0
    KAPPA (1e+05,1e+06]  rows  53 SE max rel 8.65e-08 NaN sets differ 0
    KAPPA (1e+06,1e+07]  rows  53 SE max rel 2.70e-06 NaN sets differ 0
    KAPPA (1e+07,1e+08]  rows  38 SE max rel 3.20e-05 NaN sets differ 0
    KAPPA (1e+08,1e+09]  rows  11 SE max rel 5.57e-04 NaN sets differ 0
    KAPPA (1e+09,1e+12]  rows  19 SE max rel 9.98e-01 NaN sets differ 5
    KAPPA (1e+12,Inf]    rows   7 SE max rel 1.20e-04 NaN sets differ 0

The bins use `1 + n / (sigma/sd)^2` from the logs; the code's own bound on
the mpfr design reads 1.309e+10 at both scales and refuses the route, and on
the 500-position timing design (noise 0.3) 1.056e+04
(`dev/fitedge-kappa-mpfr.R`, `kappa-mpfr.txt`).

The value check alone does not see this: the largest value difference is
3.1e-12 below 1e6 and 2.2e-7 above 1e9, and the seed-1 fit passed it at
6.6e-10. In `gp_marginal_outer_gr()` the bound is computed from the point:
`(sum(diag(V_field)) + max(noise)) / min(noise)`, which bounds the condition
number of `V` for any blocks.

The sets: the ncgp review's `modes2` gaussian part (168 rows), its fresh
curves (192 rows) and 17 more (`r2-se`'s 4, the 10 spellings of
"the marginal reads the tape's objective" plus `weights`, edge2's 2).

With the shipped gate (`dev/fitedge-se-verdict.R lane <set>`, then
`dev/fitedge-sev-scan.R lane`, `sev-scan-lane.txt`):

    TOTAL rows 377, route taken 246, NaN sets differ 0, warnings differ 0,
      logLik differ 0, SE max rel 2.16e-08, fitted SE max rel 3.10e-10

So the SE verdicts do not change on these sets: the same parameters lose
their standard error, the same warnings fire (compared on their first 60
characters, sorted), and every logLik is bitwise the same. The `weights` row
keeps the z tape: its marginal reads another constant (value difference
0.17), so the value check refuses it. Its Hessian was equal to 1.2e-12 in the
forced run, so this costs time only.

Time of the outer Hessian alone, both routes in one process, interleaved,
minimum of the rounds, with a control of the z tape against itself
(`dev/fitedge-hess-time.R`, data `sin(x) + 0.3 e`, seed 20261103):

    100 positions (htime-lane-100.txt): marginal 0.035 s, z tape 0.124 s,
      ratio 3.49 (blocks of 11 calls, 3 rounds; control ratio 1.118)
    500 positions (htime-lane-500.txt): marginal 2.070 s, z tape 37.220 s,
      ratio 17.98 (2 rounds; control 0.940), SEs equal to 5.4e-10
    500 positions, dev arm, earlier and quieter (htime-dev-500.txt):
      marginal 1.050 s, z tape 14.640 s, ratio 13.94 (control 0.992)

    2000 rows on 500 positions (htime-lane-2000.txt): marginal 2.220 s,
      z tape 40.310 s, ratio 18.16 (1 round; control 0.951), SEs equal to
      9.1e-10

The machine ran 11 other lanes; the two 500-position runs differ by a factor
of 2.5 in absolute time and agree on the ratio (14 to 18). The ncgp review
measured 11.3 s for the z tape and 0.87 s for the marginal.

"2000 positions" cannot be measured: `frm()` refuses an exact `gp()` over
500 positions ("gp() without k= builds a dense 2000-point covariance",
`R/frame.R`; the first 2000 run, `htime-lane-2000.txt` before the rerun,
stopped there). The 2000 point is the ncgp review's own convention
(`dev/review-ncgp-r2-prof.R`: `n` rows on `min(n, 500)` positions), where
the marginal reduces to its 500 distinct rows and a within sum of squares.
`gp_marginal_max_rows` (1000) does not bind the Hessian route
(`max_rows = Inf`): it limits the warm start's cost, and the Hessian costs
one tape and `2 * np` gradients whatever the rows.

Whole `frm()` plus `summary()` at 100 positions, interleaved in one process,
order alternating, minimum of 6 (`dev/fitedge-fit-time.R`,
`ftime-lane-100.txt`): marginal 0.86 s, z tape 0.80 s, ratio 0.93. At 100
positions the Hessian is a small part of the fit, and the difference is
inside the load's noise (the Hessian control read 1.118). The route is asked
for once per fit (`dev/fitedge-route-count.R`: 1 ask, 0.02 s).

## Item 4: latres's mi() refusal

`latres_refuse_laplace()` now refuses only when `mi_ml_frame()` leaves a
missing row of the latent residual's response that another density reads.
It still runs before the fit, at the same call site: it calls
`mi_ml_frame()` itself, which is what "after" means here. The message names
the read rows.

Why the unread case is right: `mi_ml_frame()` takes an unread missing row out
of the objective (weight 0, `mi_gone` in `build_objective()`), and its
latent stays in the chain with no observation. That is the model of a weight
of 0 on the row, and on lognormal the exact marginal of log y over the
observed rows, `N(mu, sigma^2 I + sderr^2 C[obs, obs])` per group.
Measured, `dev/fitedge-mi-latres.R`, logs `milat-base.txt`, `milat-lane.txt`:

- rellib-r18 refuses all 7 fits (5 lognormal, Gamma, Gamma read by `w2`).
- Lane, lognormal `y | mi() ~ x + ar(t, g, cov = TRUE)`, 10 groups of 6,
  8 missing, seeds 40401 to 40405: logLik against the exact marginal's maximum
  (BFGS then Nelder-Mead on the dense density, no frmtmb code): differences
  -9.161e-09, -2.956e-11, -1.002e-08, -3.098e-12, -9.462e-09 (the exact
  optimizer reached at most 1e-8 higher). The weight-0 spelling reads the same
  logLik to every printed digit and the same `beta` to 7 digits, all 5 seeds.
- Gamma (the test's data, seed 12301): -34.69797550 for both spellings, and
  every estimate equal to 8 digits.
- A value another response reads (`bf(w2 ~ mi(y))`): still refused, now with
  "where another density reads its missing values".

The read case stays refused: nobody has measured it, and I did not.

The `cov = FALSE` form (`dev/fitedge-mi-condar.R`, `mi-condar-lane.txt`):
`mi_rows_joint()` reads FALSE for a latent residual's response (it keys on
`frame$autocor`, which a latent block is not), so the test's old note that it
"does not know a latent residual" costs nothing here. lognormal
`pos | mi() ~ x + ar(t, g)` fits (ll -45.146164, code 1); Gamma(log) under
`cov = FALSE` stops at the post-link refusal and `ma()` at brms's own
refusal, both unrelated to `mi()`.

Found, not fixed: `frmtmb_latres_laplace` fires on lognormal fits here (seeds
40401, 40403, 40405: "no curvature to read its first-order correction" and
"moves log(sderr) or ar[1] by 0.4"), where the inner problem is gaussian and
the Laplace value is exact (above). It fires the same on the weight-0
spelling, so it is not from `mi()`. Item 5's gate reads finite differences of
`h3`, `h4` that are 0 in exact arithmetic, times a covariance that is near
singular where sigma and sderr trade off.

## Whole suites

One file per process, every file of core and the eight extensions, with
`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true` (the gated tier, a superset of the
ungated one), my core first (`dev/fitedge-run-par.sh fitedge-all
dev/fitedge-all.jobs 2`, log `fitedge-all.log`, sums by
`dev/fitedge-suite-sum.R`, `suite-sum.txt`). Every one of the 499 `lib:`
lines names `wt-fitedge-lib` for frmtmb; the extensions are rellib-r18's,
built from the same commit as this tree's (none changed).

                 pkg files  pass fail err skip warn
              frmtmb   308 25574    0   0   15    0
     frmtmb.coupling    11   551    0   0    5    0
          frmtmb.eam    33  1854    0   0    3    0
       frmtmb.latent    12   683    0   0    2    0
        frmtmb.learn    16   545    0   0    2    0
          frmtmb.ode    11   549    0   0    1    0
       frmtmb.sample    85  4326    0   0    1    0
      frmtmb.spatial     8   561    0   0    0    0
     frmtmb.survival    15   552    0   0    0    0
    TOTAL files 499 pass 35195 fail 0 err 0 skip 29 warn 0
    FITEDGE-ALL ran 499 of 499

The 29 skips are other gates: 15 `FRMTMB_SCALE_TESTS`, 13
`FRMTMB_DRMTMB_FIT_TESTS`, 1 `FRMTMB_FUZZ`.

One session per suite, three ways (`dev/fitedge-run-sessions.sh`, NOT_CRAN
only), core (`fitedge-sess-core.log`):

    SESSION frmtmb forward files=295/308 pass=22524 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb reverse files=295/308 pass=22674 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb preload files=295/308 pass=22767 fail=0 err=0 skip=241 warn=0

LEAK lines: forward 10, reverse 10, preload 0; untagged 0.

The eight extensions, three ways each (`fitedge-sess-ext.log`): 24 of 24
sessions, every one fail 0, err 0, warn 0; untagged LEAK lines 0 in all 24.

    coupling 11/11 files pass 551 skip 5     eam 33/33 pass 1848 skip 5
    latent 12/12 pass 683 skip 2             learn 16/16 pass 474 skip 13
    ode 11/11 pass 549 skip 1                spatial 8/8 pass 478 skip 6
    survival 15/15 pass 552 skip 0           sample 80/85 pass 3962 skip 18

(the same counts in all three orders of each package). `files=` counts the
files that recorded an expectation; the sessions ran without the brms gate,
and I did not list which 13 core files and 5 sample files recorded none.

`R CMD check --as-cran` of core (`dev/fitedge-check.sh`: built with
vignettes, checked with the manual, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
pandoc and TinyTeX on PATH, `R_LIBS` = lane library, rellib-r18, user
library; log kept as `check-core.log`, the check tree deleted):
`Status: 1 NOTE`, the expected one ("Skipping checking math rendering:
package 'V8' unavailable"). Tests 34 min, OK. No extension changed, so no
extension was checked.

## Equivalence survey pairs reached

Item 1 reaches p136 (`gp` x `nl`, `b1 * exp(b2 * x)` with `b1 ~ gp(z)`).
Item 3 reaches no sampled pair: the exact `gp()` pairs are lognormal (p094),
student (p109), a prior (p023), mvbf (p022) or a body (p136), all outside the
gate. Item 2 reaches none: every sampled nl spec has one body
(`dev/equiv-specs.R`). Item 4 reaches none (no pair holds `mi()` beside an
autocorrelation term). Remeasured, `FRMTMB_EQUIV_LIB="wt-fitedge-lib;rellib-r18"
dev/equiv-measure.R p136 dev/fitedge-equiv-out`, against the record
(`dev/fitedge-equiv-cmp.R`, `equiv-cmp.txt`):

    EQ p136 verdict equivalent -> equivalent, level L4 -> L4
    EQ    l1_const : -0.0008887333 -> -0.0008891745
    EQ    l1_grad : 0.02977121 -> 0.02977964

No verdict changes. (`l1_by_decision` reads the 0.79.1 text; it is not mine.)

## Tests, and each seen to fail on rellib-r18

Runner `dev/fitedge-run-tests.R` (`dev/release/run-tests.R` with LIB
rellib-r18; `FRMTMB_RELLIB` empty is the base arm), one file per process.

- test-gp-start.R: rellib-r18 pass 221 fail 2 err 1
  (`seenfail-gp-start-base.txt`), lane pass 229 (`t-gp-start-lane2.txt`).
  Item 1: the body row reads -98.541 below the best start. Item 3: the z
  tape's gradient is called 10 times for the Hessian where 0 is asked (the
  behavioral failure), then `gp_marginal_outer_gr` is missing.
- test-latres.R: rellib-r18 pass 297 err 1 (`seenfail-latres-base.txt`),
  lane pass 301 (`t-latres-lane.txt`). Item 4: the refusal stops the unread
  model.
- test-nl-family.R: rellib-r18 pass 320 fail 1
  (`seenfail-nl-family-base.txt`), lane pass 321 (`t-nl-family-lane.txt`).
  Item 2: -377.0 against -371.5.

The final test-gp-start.R on rellib-r18 reads the same
(`seenfail-gp-start-base2.txt`: pass 221 fail 2 err 1); the ill-conditioned
case added to the item 3 test is past the missing-symbol error there.

The guard's absent case, each in a test: item 3 keeps the z tape under a
prior, on poisson, and at noise 1e-3 (the condition gate), and counts the z
tape's gradient calls there (`> 0`); item 4 still refuses a value another
response reads; item 2 cannot reach a model with one body (by construction,
above), and the 10 `mv` and 30 `sigma` two-body fits read the same values.

OpenBLAS with LAPACK (`dev/ciharden-openblas.sh <ver> lapack`,
`dev/fitedge-openblas-run.sh <ver>`, logs `ob<ver>-*.txt`), 0.3.26 and
0.3.32 (the condition and value gates are rounding-decided branches): on
both, test-gp-start 229, test-latres 301, test-nl-family 321, test-se-check
134 pass, 0 fail, 0 error, 0 skip. `dev/ciharden-out/` is deleted.

## Item 5: the first-order Laplace correction as a fit

Not done, and not measured on fresh designs. The brief allows it only after
items 1 to 4 and only on a measurement that holds on my own fresh designs.
What stopped it:

- The measurement needs an exact MLE per design. The review's instrument
  (`dev/review-latres-tk.R`) reads the exact MLEs that `frmtmb-wt-latres`
  wrote to `dev/latres-out/m1/` and `dev/review-latres-out/`; that worktree is
  gone, and rebuilding exact marginals for non-gaussian latent residuals on
  fresh designs (quadrature over each chain) is a lane of its own.
- It changes every estimate and the logLik a latent-residual fit reports
  (`theta + shift`, Laplace + C), so every reader (`fitted()`, `vcov()`,
  `predict()`, the draws) would read a point that is not the objective's
  maximum. That needs a design decision, BREAKING in NEWS, which I did not
  want to make on an unmeasured basis.
- Time: items 1 to 4, their sweeps and the whole-suite runs used the lane's
  process budget.

Related, found and not fixed (item 4 above): the gate that item 5 would
replace, `frmtmb_latres_laplace`, fires on lognormal, where Laplace is
exact.

## Files touched

Package: `R/gp-start.R`, `R/se-check.R`, `R/par-template.R`,
`R/autocor-latent.R`, `R/compat.R` (two notes), `NEWS.md`,
`tests/testthat/test-gp-start.R`, `tests/testthat/test-latres.R`,
`tests/testthat/test-nl-family.R`. `DESCRIPTION` changed only by
roxygen 8.1.1 writing `Config/roxygen2/version` (reverted at the end; every
lane's roxygenise does this). No `man/` page changed: every new block is
`@noRd`.

Lane scripts (tracked, `dev/fitedge-*`): `-lib.R`, `-install.sh`,
`-gpbody-probe.R`, `-gpbody-proto.R`, `-gpbody-sweep.R`, `-gpbody-sum.R`,
`-gpbody-cmp.R`, `-gpbody-diag.R`, `-disc-probe.R`, `-disc-starts.R`,
`-nl2-sweep.R`, `-nl2-sum.R`, `-mi-latres.R`, `-hess-time.R`,
`-fit-time.R`, `-route-count.R`, `-hess-mpfr.R`, `-kappa-mpfr.R`,
`-scale-diag.R`, `-se-verdict.R`, `-sev-scan.R`, `-cost12.R`,
`-compat-edit.R`, `-equiv-cmp.R`, `-openblas-run.sh`, `-check.sh`,
`-run-tests.R`, `-run-session.R`, `-run-par.sh`, `-run-sessions.sh`
(the last four are `dev/release/` copies pointed at rellib-r18), and
`dev/fitedge-all.jobs`.

## Bump and floor

Minor: `frm()` now fits a model it refused (`mi()` beside a latent
residual) and changes the estimates of fits that stopped short (items 1 and
2). No exported function, argument or class changes, and no extension code
changes, so no extension needs a new floor on frmtmb for this work.
