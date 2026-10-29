## The 0.78.0 round, 2026-10-28 to 2026-10-29

`dev/round-20261029.md` is this round's record. Each lane has its
findings (`dev/equiv-findings.md`, `dev/meme-findings.md`,
`dev/rescorpp-findings.md`) and its review with the fix-round checks
(`dev/reviews/2026-10-28-<lane>.md`). Three lanes off `ac4906bd` (0.77.0
plus one test fix, pushed), merged onto `ac4906bd` in the order equiv,
meme, rescorpp:

- `equiv`: a measured equivalence survey, dev/ only. 172 frm_compat()
  pairs of the 7,387 core pairs, drawn by stratum, weighted by reach,
  each fitted by brms 2.23.0 and frmtmb and compared at four levels:
  83.6 percent reach-weighted agreement, Wald [78.2, 89.0], and 80.5
  with the weak-evidence verdicts out. Defects D1 to D8 and 53
  compat-row proposals (`dev/equiv-data/compat-proposals.tsv`), none
  applied.
- `meme`: brms's `me()` prior classes `meanme`, `sdme` and `corme` on
  both routes (BREAKING); sampling-route group rows print brms's
  `(flat)`, so an edited class row reaches the group (BREAKING);
  frm_sample() refuses `coef` with a bound for the brms classes.
- `rescorpp`: frmtmb.sample's joint `posterior_predict()` of
  `set_rescor(TRUE)` draws through core's new export
  `predict_rescor_draw()`, brms's draws x observations x responses
  layout, and `log_lik(resp =)` as brms (BREAKING).

Versions: frmtmb **0.78.0** (minor, breaking); frmtmb.sample
**0.25.0** (minor, breaking), needing frmtmb 0.78.0, measured by
`dev/rel078-floor.sh`; spatial 0.3.0, survival 0.2.0, latent 0.7.1,
eam 0.13.0, coupling 0.6.1, ode 0.7.2 and learn 0.7.3 unchanged. The
release library is `C:/Users/adf44/source/r/rellib-r16`; `rellib-r15`
keeps 0.77.0.

### Decisions made at the consolidation, for the coordinator

- **R1, the two spellings of a nonlinear intercept.** frm_sample() now
  merges `set_prior(class = "Intercept", nlpar = "a")` and brms's class
  row on `a ~ 1` by specificity, as frm() resolves the pair (measured,
  `dev/rel078-r1-frm.R`): the `Intercept` spelling's density, else the
  class row's; per side, its bound, else the class row's. The fit's own
  prior and the call's `prior =` combine the same way, so a call's class
  row does not displace a fit's more specific row, as in brms's
  `update(prior =)`.
- **Rtools last on PATH** in the per-file and one-session runners
  (`dev/ciharden-scan.sh`, `dev/release/run-sessions.sh`), as the brief
  asks after 0.77.0's bin-2 segfault; bin 1 of the ledger ran through
  `dev/rel078-record.sh`, `dev/defects-record.sh`'s loop with that PATH.
- **The equiv reviewer's scripts.** The brief left out
  `dev/review-equiv-scratch/`; the release review (finding 1) showed
  the equiv review cites its 26 scripts there, so they are copied and
  committed, as the meme and rescorpp reviewers' scripts are. Only
  `dev/review-equiv-out2/` (rerun results) stays untracked in the lane
  tree.

### What the 0.78.0 round is evidence for

**A measurement lane needs a rerun plan, not a merge.** equiv changed no
tracked file, so the merge could not show whether the other lanes moved
its verdicts; its section 6 named the pairs they reach, and the
consolidation remeasured them on the merged build into a scratch outdir
(release fix R4). Give every measurement lane such a list.

**A merge rerun into a copy shows the tested tree is the scripted one.**
The merge script takes its root and its git directory as variables;
rerun into a scratch copy, with roxygen there, and compared file by file
with the tree the tiers read, it showed 0 differences outside `dev/`.
So the release library did not have to be rebuilt after the runs, and
no run had to be repeated.

**A rule that "the later one stands" is an order dependence.** meme's
fix round 2 kept the later of two specifications that name one slot;
the reviewer's short check found that frm(), the base build and brms all
resolve the pair by specificity in either order. When two spellings
collapse into one slot, merge them by the precedence the other route
already has, and test both orders.

### The release review and its fixes

`dev/reviews/2026-10-29-release.md`: FIX, no blocker or major. The
fixes are `dev/rel078-release-fix2.R`, sourced at the end of
`dev/rel078-release-fix.R`, and edits of the merge and commit scripts:

1. The equiv review's 26 evidence scripts are copied and committed
   (the merge skips only `dev/review-*-out*/` folders).
2. `prior_summary()` of draws lists both spellings of a nonlinear
   intercept as written; the sampler reads the merged slot (release fix
   R1's printed record).
3. `bayes_R2()` and `loo_R2()` refuse a response named twice.
4. The commit script's dry run lists paths with `git status` and
   `git ls-files`, so it takes no lock (`git add -n` takes the index
   lock).
5. Nits: the skip-gate reason in the record (RTMBode is installed; it
   is missing from frmtmb.sample's Suggests, filed), Rtools last in
   `dev/release/run-par.sh`, a stale comment in the bin-2 adapter spec.
6. Finding 8's items are filed; the `me()` newdata item is MAJOR: on new
   data the latent value is the observed `x`, and the predictive sd is
   0.60 of brms's (`dev/rel078-fix2-me-newdata.R`).

**A dry run is only write-free if every command in it is.** `git add -n`
lists what it would add, and still takes the index lock. List with
commands that only read (`git status --porcelain` with optional locks
off, `git ls-files`), and check them against each other.
