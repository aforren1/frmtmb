# Lane brmssuite, item 2.6a

The deliverable is `dev/brms-suite-audit.md`. This file records what
the lane decided, what it did not do, and why.

## What changed in the tree

No package code. One tracked file changed: `.gitignore` gains
`dev/brms-suite/`, so the 434-file, 5 MB copy of brms's source does not
enter the repository. The audit carries the fetch command and the
tarball's sha256, which is what makes the copy disposable.

Everything else is new and unreferenced by the package: the audit, ten
scripts prefixed `brmssuite-`, and their `.tsv` outputs.

## Decisions

- **The plan was NOT edited.** `dev/extension-gaps-plan.md` is a shared
  file and two other lanes are running against it. The text for 2.6b's
  `days` cell is in the final report for whoever consolidates.
- **Six bins, not three.** The plan's bin 2 says "restate the tolerance
  because brms reports a posterior". Only 59 of 329 bin-2 assertions
  are that. The other 270 assert an exact identity against a
  hand-written density or a closed-form matrix, and what they need is
  an adapter to frmtmb's internals, not a tolerance. Bin 3 splits the
  same way: Stan, MCMC, and brms-internal utilities are three different
  reasons not to port, and only the first two are in the plan's
  wording. The three top-level bins are still reported.
- **Blocks were classified, not files.** Four files straddle bins; a
  file-level label would have put 87 bin-1 validation assertions inside
  a file called `tests.standata.R` and lost them.
- **Nothing was ported.** The spot check runs brms's blocks through a
  harness that evaluates one expression at a time; it writes no test
  file and installs nothing.
- **Defects were filed, not fixed.** 15 items in section 8 of the
  audit, each with its construction. Two of them were already recorded
  elsewhere and are labelled as such.

## What I could not settle

- **`(1 | g) + (x | g)`** is accepted by frmtmb and refused by brms.
  lme4 and glmmTMB accept it and mean two independent blocks, so this
  is a policy question rather than an obvious defect. The plan's
  tiebreaker rule says follow brms; that would break a spelling lme4
  users write on purpose. Left for the user.
- **`get_prior()` reports a `theta` prior class** for a poisson GLMM
  where brms reports none. It is one row in one call and I did not
  chase what puts it there.
- **Whether `variables(ds)` on a draws object agrees with
  `variables(fit)`.** The fit side drops brms's `b_` prefix. Checking
  the draws side needs a Stan compile, which this item did not budget,
  and item 2.5f is already in that code.

## Cost

Ten R processes, none concurrent, each a few seconds except the two
that fit one poisson GLMM on `epilepsy` (236 rows, 59 patients). No
installs. No `R CMD check`: this lane changes no package code, so there
is nothing for it to check.
