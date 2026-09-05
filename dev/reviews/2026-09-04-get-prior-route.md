# Review: get_prior() route= (branch wt-getprior, base a233c3c)

Reviewer run 2026-09-04, R 4.6.1. Core and the extension installed from
this worktree into a private review library; the user library was a read
fallback for dependencies only. Main checkout confirmed clean at a233c3c
and unmoved.

## Verdict: GO-WITH-FIXES

The defect is real, the fix is the right one, and every headline claim
reproduces. The punch list is documentation and test coverage; no
behavior needs to change before merge.

## 1. Attach-independence, reproduced

Probe run in three configurations against the review-library core.

| configuration | fit route identical before/after attach | sample route |
| --- | --- | --- |
| (a) no extension | TRUE (formula and fit) | refused, names frmtmb.sample |
| (b) worktree frmtmb.sample, `load_all` | TRUE (formula and fit) | student_t/lkj table |
| (c) INSTALLED frmtmb.sample 0.1.0, built before this change | TRUE (formula and fit) | student_t/lkj table |

Case (c) is the load-bearing one. The registration API did not change,
so a pre-existing installed extension against the new core is both
attach-independent on the fit route AND still able to answer the sample
route. No coordinated release is required.

The lane's own `dev/getprior-load-order-probe.R`, run verbatim from
`dev/`, prints `identical: TRUE` for both calls.

Refusal text, reproduced:

    route = "sample" reports the prior defaults frm_sample() applies,
    and no loaded package states them. They live in frmtmb.sample: run
    library(frmtmb.sample) and call again (install it from
    extensions/frmtmb.sample in the frmtmb repository). For the
    defaults frm() applies, use route = "fit".

## 2. Call sites and consumers

`get_prior()` has NO runtime caller in any of the five packages. Every
occurrence in `R/` across core, frmtmb.sample, frmtmb.latent,
frmtmb.ode and frmtmb.ddm is prose in a comment or roxygen block. It is
a leaf, user-facing function. The lane's weaker claim ("no core caller
reads the prior column") is true and understated.

Test-side consumers, all checked and all unaffected:

- `tests/testthat/test-v19.R:78-96` reads `$class`, `$coef`, `$group`,
  `$dpar` by name. Passes.
- `tests/testthat/test-priors-autocor-classes.R:182-196`, `$class`
  membership. Passes.
- `tests/testthat/test-prior-compat.R:171,204,373`, column reads.
  Passes.
- `tests/testthat/test-brms-agreement.R:463` wraps **brms's**
  `get_prior()`, not frmtmb's. Untouched.
- `extensions/frmtmb.ddm/tests/testthat/test-surface.R:121-128` asserts
  `is.data.frame(gp)`, which is still TRUE.
- `R/priors.R:433 as_priorlist()` gates on `inherits(x, "brmsprior")`,
  so a `frmtmb_prior_rows` falls through untouched, exactly as a plain
  data frame did.

Interop probed directly on the returned object:

| operation | result |
| --- | --- |
| `head()`, row subset `gp[i, ]` | class kept, `route` kept, label prints |
| column subset `gp[, cols]` | class kept, **`route` dropped**, no label |
| `gp$prior`, `gp[["prior"]]` | plain character vector |
| zero rows (`gp[0, ]`, `subset()`) | label plus `<0 rows>`, prints sanely |
| `rbind(gp, gp)` | class and `route` kept from the first argument |
| `merge()` | plain data.frame |
| `as.data.frame()` | class dropped, **`route` attribute retained** |
| `transform()` | class and `route` both dropped |
| `tibble::as_tibble()` | clean `tbl_df`, 9 columns |
| `knitr::kable()` | renders normally, no label row |
| `knitr::knit_print()` | dispatches to the print method, label included |
| `write.csv()` round trip | 8 rows, clean |

Nothing breaks. Two asymmetries are noted in the punch list.

The extension vignette's route chunks are `eval = FALSE`
(`extensions/frmtmb.sample/vignettes/sampling.Rmd:65`), so knitr never
renders the object there; the `knit_print` path was exercised directly
instead and is correct.

## 3. The test seam

`swap_prior_defaults()` and `require_prior_defaults()` are both absent
from the `@rawNamespace export(...)` list at `R/sampling-api.R:189-235`
and from `NAMESPACE`. Unexported, confirmed.

**The review brief's premise is wrong, and in the lane's favor.**
`extensions/frmtmb.sample/tests/testthat/test-prior-route.R` contains no
`:::` at all. It reaches nothing internal: it calls
`frmtmb::get_prior()` and reads `attr(x, "route")`. Only core's
`tests/testthat/test-get-prior-route.R:24-25` touches the seam, which is
the package that owns it. Nothing goes on the extension-test `:::`
burn-down list.

This is already the cleaner design the brief hypothesized: core proves
attach-independence with a synthetic provider through the swap; the
extension proves its own registration through the public route without
emptying anyone's registry. The synthetic provider at
`test-get-prior-route.R:16-19` is deliberately not frmtmb.sample's, so
core's test fails for the right reason whether or not the extension is
installed. Correct.

## 4. The dontrun decision

Mechanically the lane's account is exact. `R CMD check` renders
`\dontrun` content into `frmtmb-Ex.R` as `##D` comment lines (verified
in the check directory), so `.check_packages_used_in_examples` never
parses `library(frmtmb.sample)` and cannot flag it. A
`requireNamespace()` guard is live code and would be flagged, because
core cannot Suggest a package that Imports core.

I would still drop it, for prose. Reasons:

1. `\dontrun` in CRAN idiom means "running this is inappropriate" (too
   slow, needs credentials, side effects). Here running it is perfectly
   appropriate; the block exists to hide an undeclared dependency from a
   check. Using `\dontrun` to route around a check is the wrong tool
   even when it works.
2. Rendered help shows `## Not run:`, which reads to a user as "this is
   broken or dangerous". The opposite of the intent.
3. It is redundant. The new "Which route the defaults describe" section
   (`R/priors.R:646-668`) already explains the sample route in three
   paragraphs. It just never shows the literal call.
4. Nothing executes it. `--run-donttest` does not run `\dontrun`, so
   this half of the example is verified nowhere. Its behavior IS covered
   by `extensions/frmtmb.sample/tests/testthat/test-prior-route.R:18`,
   so the example adds no assurance, only maintenance surface.

Suggested replacement: delete the `\dontrun` block from the example, and
add one sentence with an inline call to the existing route section,
naming `route = "sample"` in code font. Requires a roxygen re-run, so I
left it for the lane rather than editing generated Rd.

Either choice is defensible; this is a preference with reasons, not a
blocker.

## 5. Test and check results (all reproduced)

One file per process, `pkgload::load_all`, `NOT_CRAN=true`:

| file | pass | fail | skip |
| --- | --- | --- | --- |
| test-get-prior-route.R | 29 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 |
| test-bracket-access.R | 7 | 0 | 0 |
| test-prior-compat.R | 105 | 0 | 0 |
| test-priors-autocor-classes.R | 63 | 0 | 0 |
| test-brms-agreement.R | 168 | 0 | 3 |

29 on the new file matches the lane exactly.

Full core suite, one process, `load_all`, `NOT_CRAN=true`:
**102 files, 5282 pass, 0 fail, 4 skip.** The lane reported 5272. The
count is environment-dependent (skip guards fire differently depending
on what is installed); under `R CMD check` the same suite reports 5267
pass / 9 skip. File count and zero failures agree everywhere. Not a
defect, but the lane should quote a pass count with its environment or
not quote one.

Extension suite, one process, against the review-library core:
**134 files, 881 pass, 0 fail, 0 error, 2 skip.** Matches the lane.

`R CMD check --as-cran --no-manual` with `_R_CHECK_CRAN_INCOMING_=false`,
on a freshly built tarball: **Status: OK.** No NOTEs, no WARNINGs.
Examples OK (37s), `--run-donttest` OK (43s), tests OK (337s), vignette
re-build OK (138s), "unstated dependencies in examples" OK.

## 6. Changed surface

Tracked, vs a233c3c: `NAMESPACE`, `NEWS.md`, `R/priors.R`,
`R/sampling-api.R`, `extensions/frmtmb.sample/NEWS.md`,
`extensions/frmtmb.sample/R/zzz.R`,
`extensions/frmtmb.sample/vignettes/sampling.Rmd`,
`man/frmtmb-sampling-api.Rd`, `man/get_prior.Rd`.
Untracked: `dev/getprior-load-order-probe.R`,
`extensions/frmtmb.sample/tests/testthat/test-prior-route.R`,
`tests/testthat/test-get-prior-route.R`.

Exactly the stated surface, nothing else. Main checkout clean at a233c3c.

## Punch list

1. **`R/priors.R:835-841`** (fixed by the reviewer, comment only). The
   print method's comment said "A subset keeps it, correctly" without
   qualification. Row subsets keep `route`; column subsets do not,
   because `[.data.frame` drops the attribute. Reworded to say which is
   which and why the unlabeled case is fine.
2. **`tests/testthat/test-get-prior-route.R`, near line 131.** The print
   method's no-label branch (`R/priors.R:846-849`) is reachable in
   ordinary use, via `gp[, c("prior", "class")]`, which keeps the class
   and loses the attribute, and no test covers it. Add the column-subset
   case beside the existing row-subset assertion.
3. **`R/priors.R:826-830`.** `as.data.frame()` strips the class but
   leaves a stray `route` attribute behind, so a coerced table is not
   `identical()` to the plain data frame a reader would expect. The
   extension's own test hand-rolls a `bare()` helper at
   `extensions/frmtmb.sample/tests/testthat/test-prior-route.R:47-51`
   for exactly this reason, which says the strip belongs in core. Either
   add `as.data.frame.frmtmb_prior_rows()` that drops both, or document
   the attribute as sticky. Low severity, no current consumer trips on
   it.
4. **`dev/getprior-load-order-probe.R:14-16`.** The `load_all` fallback
   resolves only when the working directory is `dev/`. From the repo
   root the expression yields
   `C:/Users/adf44/source/frmtmb-wt-getprior/extensions/frmtmb.sample`,
   one directory level short, and the branch fails. Verified. The
   fallback is dead whenever frmtmb.sample is installed, which hid it.
   Note the required cwd in the header comment or drop the arithmetic.
5. **`R/priors.R:690-699` and `man/get_prior.Rd`.** The `\dontrun{}`
   block. See section 4; prose in the existing route section is the
   better shape. Needs a roxygen re-run, so not applied here.
6. **Lane report corrections.** (a) The extension test does not use
   `:::`; nothing joins the burn-down list. (b) "survives subsetting"
   holds for rows only. (c) The 5272 figure is environment-dependent.

## Version bump

**Minor: 0.49.1 to 0.50.0.** `get_prior()` returns a different class
with a new attribute, prints a line it did not print before, and gives a
different `prior` column to anyone whose session had frmtmb.sample
loaded. The NEWS entry calls it a behavior change and it is one. A patch
bump would understate it.

## Reviewer edits

One, comment text only, no code path touched:
`C:\Users\adf44\source\r\frmtmb-wt-getprior\R\priors.R` lines 835-841,
punch-list item 1. File re-parsed clean after the edit.
