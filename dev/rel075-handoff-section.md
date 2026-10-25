## The 0.75.0 round, 2026-10-24

`dev/round-20261024.md` is this round's record. Each lane has its
findings (`dev/smallfix-findings.md`, `dev/bhaz-findings.md`,
`dev/brmsnames-1024-findings.md`, `dev/sarfcor-findings.md`) and its
review with the short-check sections (`dev/reviews/2026-10-24-<lane>.md`).
Four lanes off `18b536fc` (0.74.0 plus two CI commits), merged onto
`4fd134db` (18b536fc plus one CI commit to `pkgdown.yaml` and
`dev/release/build-docs.R`) in this order:

- `smallfix`: one refusal names every package a saved fit lacks
  (`stop_extension_parts()`), `keep_missing()` in every handler that
  could swallow it, `require_fit_packages()` at the head of methods; the
  `mixture_mvn()` test debt; the Dirichlet refusal at the model
  (BREAKING) and `print(curve, digits =)`; the CI skip gate for core.
- `bhaz`: `frmtmb_register_aterm(option =)`; brms's `bhaz()` in
  frmtmb.survival, and `cox(knots =, Boundary.knots =)`.
- `brmsnames`: brms's coefficient names where the parser's alias
  leaked (`.frm_expx` is `expx`), the prior table in brms's names and
  order, `fixef()` order of nonlinear and multivariate rows (BREAKING).
- `sarfcor`: `frmtmb_register_autocor()`, core's `fcor()`,
  frmtmb.spatial's `sar()`, `frmtmb_autocor_pars()` and the `bounds =`
  prior slot (part A); `frmtmb_covstruct_info()`, `full_rank` and the
  autoscale message (part B). The user kept both parts.

Versions: frmtmb **0.75.0** (minor, breaking); frmtmb.survival
**0.2.0**; frmtmb.spatial **0.2.0**; frmtmb.sample **0.22.0** (minor:
its NEWS has breaking changes and new features); frmtmb.latent
**0.7.1** (tests only); eam 0.13.0, coupling 0.6.1, ode 0.7.2 and learn
0.7.3 unchanged. All four changed extensions need frmtmb 0.75.0,
measured by `dev/rel075-floor.sh`. The release library is
`C:/Users/adf44/source/r/rellib-r13`; `rellib-r12` keeps 0.74.0.

### Decisions made at the consolidation, confirmed by the coordinator

The coordinator confirmed the decisions below, except one: the
consolidation kept smallfix's head check in
`recover_data.frmtmb_fit()`, and the release review showed that it
broke emmeans (next section).

- **One refusal across three seams.** sarfcor's residual-term check
  joins smallfix's collector (`seam_missing_parts()`), the "session"
  sentence names the seam that registers the part, and
  `ac_block_entry()` refuses through the same builder.
- **Two integration fixes the merged tree needed**, found by the
  consolidation's own test: `autocor_matrix()` and
  `frmtmb_autocor_pars()` check the whole fit first, and
  `frmtmb_covstruct_info()` uses `keep_missing()`.
- **The consolidation's test** of a saved fit with parts from three
  absent packages, and of a real saved `car()` and `sar()` fit. The
  brief's "cox + car + sar" cannot be one model: `sar()` takes a
  gaussian or student response, and `cox()` is refused in a
  multivariate fit.
- **frmtmb.sample is 0.22.0**, a minor bump, and **frmtmb.latent
  0.7.1** with the floor 0.75.0, because a test pins 0.75.0's message.
- **The upstream entries are renumbered**: brms-36 is bhaz's
  interval-censored `cox()`, brms-37 to brms-39 are sarfcor's; the
  references in sarfcor's findings, review and scripts follow.
- **The bin-2 `sar()` and `fcor()` rows are not ported** in this pass
  (record, "The ported brms suite").

### The release review and its fixes

`dev/reviews/2026-10-24-release.md`: FIX, no blocker. The coordinator
asked for every finding but 7 (a test-count nit, filed) as a release
fix (`dev/rel075-release-fix.R`, sourced by the resolve step):

1. (major) `emmeans()` and `ref_grid()` of a fit with a missing part
   stopped with emmeans's "Perhaps a 'data' or 'params' argument is
   needed": smallfix's head check in `recover_data.frmtmb_fit()` is
   raised inside emmeans's `try()`. The head check is removed;
   `emm_basis.frmtmb_fit()` refuses by name; both calls are in the
   saved-fit test.
2. `frmtmb_covstruct_info("sar")` names a residual correlation term.
3. One `sar()` parameter, one prior slot: no internal row in
   `default_prior()`, and a prior on it is refused by its class.
   `ar()`'s same double slot is filed.
4. frmtmb.survival's README and vignette show `bhaz()`.
5. The saved-fit step's comments name both fits.
6. `class = "lagsar"` on an error SAR names the model's `errorsar`.
7. `dev/rel075-floor.sh` installs only outside the tree and keeps
   Git's bash and PATH; rerun.

### What the 0.75.0 round is evidence for

**A head check inside someone else's `try()` is lost.** smallfix put
the fit-level check at the head of every method that could meet a
missing package, `recover_data()` among them; emmeans calls that one
inside `try()` and discards what it raises. The consolidation read the
code's own comment, kept the check, and recorded that it "refuses
first": it did, and nobody saw it. A check is placed where its caller
lets it through, and a test calls the method the way users do
(`emmeans::ref_grid()`), not the method directly.

**A safety check's refusal is reported, not worked around.** When the
tool refused a deletion written with shell variables, the
consolidation reran it with literal paths. The rule is to report the
denial and leave the deletion to the user; the record files it as an
incident.

**A merge can compose two lanes into a gap neither lane had.** Each
lane's saved-fit test passed in its own tree; only a test written for
the merged tree showed `autocor_matrix()` naming the `sar()` term alone
in a fit that also lacked a structure's package. Every integration of
two seams gets a test at the consolidation, written against the merged
tree and seen to fail without the resolution (`dev/rel075-seenfail.R`).

**A merge rerun restores generated files.** `dev/rel075-merge.sh`
restores `NAMESPACE` and `man/` to the HEAD blob, so an install without
roxygen after a rerun builds the old exports; the extensions then fail
to load ("'frmtmb_register_autocor' is not exported"). Rebuild the
release library in full after every merge rerun.

**A brief's example can be impossible.** "A saved fit with cox + car +
sar" cannot be fitted; measure that first, then test the nearest
construction and say why.
