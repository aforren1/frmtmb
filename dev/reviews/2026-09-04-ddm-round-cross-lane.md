# rv9 findings, part 2 (gddm, D, E, F)

## PER-LANE: gddm
GD1 CALL COUNTER (gd_solve wrapped in a counter, through a real frm() fit):
  2400 rows over 6 conditions -> 6 solves for the WHOLE FIT.
  Scaling: 600/3->3 ; 2400/6->6 ; 2400/3->3 ; 1200/12->12 ; 4800/6->6. Exactly 1.00/condition.
  STRONGER than claimed: after the tape is built, one objective eval and one gradient eval
  each trigger 0 further R-level solves. The solve is taped once and replayed in compiled
  code. Claim confirmed and then some.

GD2 ACCURACY TABLE at the DEFAULT grid (dt=0.01, ny=201), reproduced EXACTLY:
    decision time >= 0.20 s -> 0.0068   (lane: 0.0068)
    decision time >= 0.08 s -> 0.7753   (lane: 0.775)
    decision time >= 0.05 s -> 4.3299   (lane: 4.33)
  Grid behaviour as documented: coarse(0.02/101) 0.07604, default 0.00681,
  fine(0.005/401) 0.00204.

GD3 EARLY-TIME LIMIT on a Roitman-Shadlen-like design (6 coherences, 6000 trials, mu=2.2 bs=2.0):
    ndt=0.25 median RT 0.861 : dec times <0.05s 0.02% ; <0.10s 0.57% ; <0.20s 7.03%
    ndt=0.30 median RT 0.892 :                  0.00% ;         0.53% ;         8.20%
    ndt=0.35 median RT 0.952 :                  0.00% ;         0.48% ;         7.42%
  ~0.5% of trials sit where the log density is wrong by ~0.8 (a factor of 2.2); essentially
  none where it is wrong by 4.3; but 7-8% sit under 0.20 s, in the transition. Exposure on
  the paper's own design is real but modest, and the lane states it fairly.

GD4 WHAT THE LAPSE REMEDY COSTS: on clean data, NOTHING.
  no lapse logLik -2345.845 (5 par); lapse=uniform logLik -2345.845 (6 par); mu/alpha/bs
  identical to 4 dp. The lapse rate goes to its boundary when there are no fast
  contaminants, so it is free insurance and only pays where lapses exist.

GD5 RENORMALIZATION PIN - justified:
  renormalize=TRUE  mu=2.2744 bs(log)=-0.0474 logLik -2345.845
  renormalize=FALSE mu=2.2183 bs(log)=-0.0955 logLik -2400.531
  logLik falls 54.7; boundary separation shrinks ~5% (exp(-0.0474)=0.954 vs 0.909).
  Direction matches; my magnitude is smaller because I fitted a constant drift, not the
  leak model the lane measured.

GD6 THE LEAK: a MULTI-START problem? NO.
  3000 rows, truth mu=2.5 bs=2.0 leak=1.5 ndt=0.28. Single fit mu=2.3982 (4% low),
  leak=1.3983 (7% low), bs(log)=0.6951 vs truth 0.6931 (exact).
  frm_allfit(): nlminb -3229.967 | optim -3229.967 | bobyqa -3229.967 | nloptr_lbfgs -5448.286
  THREE independent optimizers agree on the SAME optimum. The residual drift bias is a
  property of the likelihood at this grid (discretization), NOT a poor optimum, so
  frm_allfit() does NOT resolve it. (nloptr_lbfgs is an outlier that fails badly - a core
  frm_allfit observation, not this lane's.)

GD7 ZERO-COHERENCE GRADIENT: finite, no NaN.
  AD gradient on a 3-condition design whose first condition is C=0:
    (-3.78e-16, 2.16e-15, 1.07e-15, 9.22e-16) - all finite, no NaN.
  The naive d/dalpha of C^alpha = C^alpha log(C) is NaN at C=0. Defect genuinely avoided.

GD8 THE MULTI-ALTERNATIVE REFUSAL (gddm.R:1025-1035) against the three requirements:
  names the count?           YES
  says why?                  YES (single accumulator, two absorbing boundaries; needs
                             racing accumulators or a multi-dimensional Fokker-Planck solve)
  points at the alternative? NO - AND IT BECOMES FALSE ON MERGE.
    It ends "Neither is in this package. Collapse the response to two alternatives, or fit
    each pair separately." After the lba lane lands, racing accumulators ARE in this
    package: lba(n) is exactly the named alternative. The message must point at lba(n).
    A CROSS-LANE defect neither lane could see alone.

## D. ONE STORY FOR THE DECISION INDICATOR - measured in the merged package
  wiener + dec(factor)       works     (wiener-family.R:437 reads aterms[["dec"]])
  wiener + vint(0/1)         works
  gddm   + dec(factor)       REFUSED   "the density needs vint2, which nothing on this
                                        response supplies"
  gddm   + vint(upper, cond) works     vint1 = boundary 0/1, vint2 = condition index
  lba    + dec(factor)       REFUSED   "a decision indicator has two levels ... this one has 3"
  lba    + vint(1..n)        works     vint1 = which accumulator won, counting from 1
  gddm never reads aterms[["dec"]]; lba never reads aterms[["dec"]]. Only wiener does.

WHAT A USER OF ALL THREE SEES TODAY: three spellings for "which response happened".
  dec() works on exactly one of the three families. The same 0/1 boundary is dec() OR vint1
  in wiener but vint1-only in gddm, where vint2 additionally carries a CONDITION index that
  is not a decision at all. lba overloads vint1 again with a 1..n index whose 1-based coding
  contradicts wiener/gddm's 0-based boundary. A factor response works in wiener, is refused
  by gddm for the wrong reason (a missing vint2), and is refused by lba with a message about
  "boundaries" for a family that has none.

RECOMMENDED UNIFICATION (do not implement now):
  1. Make dec() the one spelling across all three; its registered coercion returns a 0-based
     integer level index of arity 1. wiener keeps 0/1; lba reads level+1 as the accumulator;
     gddm reads it as the boundary.
  2. Move gddm's CONDITION index out of vint2 into its own registered aterm - cond() -
     because it is not a decision and does not belong in that positional vocabulary.
  3. Keep vint() working everywhere as the raw-integer escape hatch, unchanged.
  4. Then one compat row per family says "dec() works", and the refusals name one spelling.

## E. zzz.R COMPAT ROWS - the exact punch-list
Both gddm and lba ship the BASE zzz.R unchanged:
  frmtmb_register_compat(features = c(wiener = "family"), rules = ddm_compat_rules)
VERIFIED in the merged package: frm_compat() has 5005 rows; gddm appears in NONE;
lba appears in NONE. Both new families are invisible to the compatibility surface.

NEEDED in the merged R/zzz.R:
  (a) features = c(wiener = "family", gddm = "family", lba = "family")
  (b) gddm rows, matching wiener's 13-slot vocabulary:
      dec()         refused      never reads aterms[["dec"]]
      vint()        works        vint1 boundary 0/1, vint2 condition index; BOTH required
      vreal()       works        carries the per-condition covariate
      cens()        refused      no lcdf
      trunc()       refused      no lcdf
      weights()     untested
      simulate      works        gddm_simulate() draws from the solved density
      fitted        works        fam[["post"]]$mean_fn defined at gddm.R:1126
      predict       works
      residuals     conditional  response only, as wiener
      residuals_osa untested
      REML          untested
      mixture()     untested
      quadrature    refused      by core, same reason as wiener
  (c) lba rows:
      dec()         refused      an n-alternative index, not a 0/1 boundary
      vint()        works        required_aterms = "vint1", the 1..n accumulator
      cens()        refused      no lcdf
      trunc()       refused      no lcdf
      weights()     untested
      simulate      works        sim = lba_sim_rt (lba.R:260)
      fitted        REFUSED      no post$mean_fn; the race mean has no closed form
      predict       works
      residuals     conditional  response only
      residuals_osa untested
      REML          untested
      mixture()     untested
      quadrature    refused

## F. NaN AND -Inf: THE TWO FAMILIES DISAGREE
  wiener: floors, so a below-support row gives a finite very negative log density, exp()
          exactly 0 and gradient exactly 0. Silent.
  gddm  : does not floor. VERIFIED: all 12 suite warnings are nlminb's "NA/NaN function
          evaluation" (10 at test-gddm-recovery.R:49:5, 1 at :86:3, 1 more in that file).

RECOMMENDED POLICY: adopt wiener's floor everywhere.
  WHY: a NaN objective is not a value an optimizer can use. nlminb survives it only because
  it treats NaN as "step rejected" and backtracks - luck in the line search, not a contract.
  Inside mixture() a NaN propagates through the log-sum-exp and takes every other component
  with it, which is precisely the defect ddmvar already had to fix once (test-defects.R).
  Two families in one package disagreeing on this guarantees that mixture(gddm(), ...)
  reproduces the bug mixture(wiener(), ...) no longer has.
  COST TO gddm: the 12 warnings disappear, and with them a real diagnostic signal that the
  optimizer is probing a region the grid cannot represent. Recover that as an explicit
  post-fit check - count floored rows at the optimum and warn ONCE with the count - rather
  than 12 optimizer warnings a user cannot act on.
  COST TO wiener: none; it already does this.
  The floored value must be finite-and-flat (as ddm_floor gives), never -Inf: -Inf
  differentiates to NaN, which is the same failure in a new coat.
