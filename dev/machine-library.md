# The R user library on this machine, and how to restore it

The user library `C:/Users/adf44/AppData/Local/R/win-library/4.6` has
been destroyed three times since it was created on 2026-08-31: around
2026-09-04, around 2026-09-07, and on 2026-09-09. Each time every
package directory survives and every FILE inside it is gone. On
2026-09-09 that was 375 directories with zero `DESCRIPTION` files
between them.

## What the signature rules out

- **Not a concurrent install.** R installs by emptying one package
  directory and unpacking into it, so a collision hollows THE PACKAGE
  BEING INSTALLED, not all 375. `dev/lane-rules.md` attributes the
  earlier losses to concurrent installs and that attribution is
  probably wrong. The private-library rule it justifies is still right,
  because a collision does real damage; the stated cause is not.
- **Not `unlink(recursive = TRUE)`** from a script, which would take the
  directories with the files.

A file-level sweep that walks the tree and deletes files is what fits.

## The suspects, none excluded

1. **Windows Disk Cleanup.** `\Microsoft\Windows\DiskCleanup\SilentCleanup`
   ran at 04:52 on 2026-09-09, about forty minutes before the loss was
   found. It is triggered by the low-disk-space condition rather than a
   schedule, and `C:` had 48.8 GB free of about 951 GB, near five
   percent. The library sits under `%LOCALAPPDATA%`, which is R 4.2 and
   later's Windows default and is inside the tree these cleaners treat
   as disposable. Storage Sense is also on, with "delete temporary files
   my apps aren't using" enabled.
2. **A managed endpoint agent.** `Get-MpPreference` is unavailable on
   this machine, so Defender is policy-managed or replaced. An endpoint
   agent that quarantines files produces the same signature.
3. **Interrupted install machinery.** Added after the fifth loss. R
   installs by deleting a package directory's contents and then
   unpacking, so a hollow directory is what an install that began and
   did not finish leaves behind. This is the only one of the three
   that explains a sweep reaching `C:/Users/adf44/source/r/rellib-r3`,
   which is not under `%LOCALAPPDATA%` and is not a temporary
   directory. See the fifth loss below.

Evidence that would separate the first two is not available here. The Task
Scheduler operational log holds about four hours (1451 records spanning
01:32 to 05:33 on 2026-09-09), so it cannot reach a previous loss, and
raising its size needs rights this account does not have.

## The fix

Move the library out of `%LOCALAPPDATA%`. In `~/.Renviron`:

    R_LIBS_USER=C:/Users/adf44/R/win-library/%v

Nothing sweeps `C:/Users/adf44/R`. Disabling `SilentCleanup` or freeing
disk space treats one suspect only.

`ZZZ-canary.txt` now sits in the library. If it is gone and the package
directories are not, the loss was a file-level sweep; compare its
timestamp against
`Get-ScheduledTaskInfo -TaskPath '\Microsoft\Windows\DiskCleanup\' -TaskName SilentCleanup`.

## Restoring it

A manifest of what was installed sits beside the library as
`win-library-4.6-manifest-<date>.csv`, with a package list in
`restore-cran.txt`. Regenerate the manifest after any large install.

The 2026-09-09 restore took about ten minutes. Install nothing by hand:
list the hollow directories, take the ones CRAN has a binary for, and
install those in one call.

    LIB <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
    .libPaths(LIB)
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    dirs <- list.dirs(LIB, recursive = FALSE, full.names = FALSE)
    hollow <- dirs[!file.exists(file.path(LIB, dirs, "DESCRIPTION"))]
    hollow <- setdiff(hollow, grep("^(00LOCK|frmtmb)", hollow, value = TRUE))
    av <- rownames(available.packages(type = "binary"))
    install.packages(intersect(hollow, av), lib = LIB, type = "binary",
                     dependencies = FALSE)

That recovered 364 of 364 with no warnings. Four packages are not on
CRAN and come from r-universe, with `dependencies = FALSE` and the same
explicit `lib`:

- `RTMBode` and `RTMBp` from `https://kaskr.r-universe.dev`
- `autotest` and `pkgcheck` from `https://ropensci.r-universe.dev`

`RTMBode` is now 1.0 from r-universe. Work recorded in `dev/` before
2026-09-09 was measured against RTMBode at commit 5242257, so a number
from `frmtmb.ode` that will not reproduce may be a version difference
rather than a wrong record.

The `frmtmb` packages themselves are excluded above because they are
installed from a checkout, not from a repository.

## The sixth loss, 2026-09-17, and the first one with a dated trigger

136 of 401 packages in the user library went hollow, and `frmtmb` in
the release library `rellib-r3` went with them. `pinlib` was untouched,
as in every previous loss. The directory mtimes are 17:09:30 to
17:10:05.

**What was running, and what the session did.** The scale tier and
`R CMD check --as-cran` were running together when the user asked for a
pause. The session stopped both drivers and then killed every remaining
`Rscript.exe` and `Rterm.exe` with `taskkill /T /F`. The damage
timestamps are that minute. That is the closest this failure has come to
a named trigger, and it points at the third suspect above, an
interrupted install: `R CMD check` installs the package under check, and
a hard kill during an install leaves exactly this signature, because R
installs by emptying a directory and then unpacking into it.

**What is NOT explained.** No step of either driver installs into the
USER library, and 136 unrelated packages went with it. So a check-time
install alone does not account for the blast radius, and the cause is
still open. Recorded because a dated, attributable trigger is more than
any earlier loss had.

**What to do differently.** Do not hard-kill R processes to pause a
release. Stop the drivers and let the R processes finish the file they
are on; if one must die, prefer stopping the driver and waiting.

**The restore.** `dev/release/restore-library.R` recovered 136 of 136
with no warnings, and the four non-CRAN packages (RTMBode, RTMBp,
autotest, pkgcheck) survived. The eight frmtmb packages were reinstalled
from the checkout. Verified by FITTING, as this file requires: gaussian
`(1 | g)` logLik -304.480044337 and poisson logLik -307.276634783 on the
seed-1 design in the session scratchpad, with brms and testthat loading.

## What this costs a round

Every private library and every reference build is derived from this
one, so a loss invalidates all of them. Rebuild the round's shared
reference library after a restore, and treat any measurement taken
across the boundary as suspect.

## The StanHeaders trap, which a restore walks straight into

`rstan` 2.32.7 declares `StanHeaders (>= 2.32.0)`. The bound is open at
the top, so `install.packages()` resolves it to the newest StanHeaders
on CRAN, which is 2.39.1. rstan 2.32.7 cannot compile against 2.39
headers. Every test block that compiles a FRESH Stan program then dies
in `compileCode()` at `make: *** Error 1`, while every block served from
`FRMTMB_STAN_CACHE` passes, so the damage looks selective and does not
look like a toolchain problem.

Measured on 2026-09-09, after the restore above: 20 of 23 gated files
green, and the three failures were exactly the ones compiling something
new. `tests/testthat/test-brms-agreement.R` came in at 168 passing
against the 187 the 0.55.1 release recorded, because its two blocks at
lines 522 and 554 call `rstan::stan_model()` directly rather than
through the cache.

The remedy is the one `frmtmb.sample`'s own `check_tmbstan_build()`
names: install StanHeaders 2.32.10 to match rstan 2.32.7.

This is the same defect as the CI failure of the same day, seen from
the other side. There, public RSPM built `tmbstan` against StanHeaders
2.39 and every chain sampled a standard normal. Here, rstan cannot
build against 2.39 at all. One open version bound, two very different
symptoms, and only one of them was loud.

So: after any restore, pin StanHeaders before believing a Stan-backed
tier, and check `packageVersion("StanHeaders")` against
`packageVersion("rstan")` rather than checking that the suite is green.
A cache makes a green suite the weaker evidence.

## What was decided, 2026-09-09

The library stays under `%LOCALAPPDATA%`. The user's call, made knowing
it has been destroyed three times in nine days. `ZZZ-canary.txt` is in
place so the next loss identifies itself as a file-level sweep rather
than being rediscovered by a failing build, and the restore recipe
above recovers 364 packages in about ten minutes.

The user library also keeps StanHeaders 2.39.1. Rather than downgrade
it, the pin lives in a separate read-only library,
`C:/Users/adf44/source/r/pinlib`, which holds StanHeaders 2.32.10 and
the tarball it was built from. A lane puts it on `.libPaths()` between
its own private library and the user library. `dev/lane-rules.md` has
the ordering and the two version checks worth running before believing
a Stan-backed tier.

That split is deliberate. The user library is the one that keeps being
destroyed and the one a restore rewrites to latest, so a pin placed
there does not survive either event. A pin outside it survives both,
and it makes the version question explicit at the top of every lane
script instead of implicit in whatever the last restore happened to
install.

## The fourth and fifth losses, 2026-09-10, and why the cause is NOT settled

Two losses happened on the same day, about eight hours apart. After
the fourth, this file said the cause was settled. **That was wrong,
and the fifth loss is what shows it.** The claim is withdrawn here
rather than edited quietly, because a confident wrong diagnosis costs
more than an open question: it tells the next session to stop looking.

### The two signatures, side by side

| | fourth, 09:21 | fifth, 17:24 |
|---|---|---|
| user-library directories emptied | 136 of 375 | 74 of 375 |
| window | not measured | **6 seconds**, 17:24:44 to 17:24:50 |
| `ZZZ-canary.txt` | **survived** | **destroyed** |
| `rellib-r3`, outside `%LOCALAPPDATA%` | untouched, 8 of 8 | **`frmtmb` emptied**, at 17:24:29 |
| `pinlib`, outside `%LOCALAPPDATA%` | untouched | untouched |
| packages removed outright rather than hollowed | none noted | `Matrix`, `mgcv` |
| `SilentCleanup` | ran at 09:21:06 | not established |
| free space released | 48.8 to 110.7 GB | not measured |
| preceded by | a consolidation running | a session ending |

### Why the fourth loss did not settle anything

The argument was: a sweep that spares a file written the previous day
while removing 136 older package trees is selecting on ACCESS TIME,
which is what Disk Cleanup does. The reasoning is sound. The evidence
was one observation of one file, and a single canary surviving one
sweep cannot distinguish "spared because recently accessed" from
"spared by chance". The fifth sweep destroyed a canary of the same
age.

Worse for that theory, the fifth sweep reached **`C:/Users/adf44/source/r/rellib-r3`**,
which is not under `%LOCALAPPDATA%`, not a temporary directory, and not
anywhere Disk Cleanup sweeps. It emptied exactly one package there,
`frmtmb`, and left the other seven intact.

### What the fifth signature looks like instead

R installs a package by deleting the target directory's contents and
then unpacking into it, so **a hollow package directory is the
signature of an install that began and did not finish**. That is
already the stated reason for this project's absolute rule about
private libraries. Seventy-four directories emptied in six seconds,
one specific package emptied in a second library fifteen seconds
earlier, and two packages removed outright, is a much better match for
install machinery interrupted mid-flight than for a disk cleaner.

It is not proof. Nothing here establishes which process it was, and
the fourth loss's `SilentCleanup` timing and 62 GB of released space
remain real and unexplained by this reading. **Both mechanisms may be
present.** That is the honest state.

### What to collect next time, so the sixth loss decides it

The two theories differ in ways that are cheap to distinguish, but
only if the evidence is captured before the restore overwrites it:

- **The exact mtime spread.** A cleaner walking a tree and an
  interrupted install have different timing shapes. Record the full
  sorted list, not the range.
- **Whether anything outside `%LOCALAPPDATA%` was hit.** This is the
  single most discriminating fact and it takes one command. Disk
  Cleanup does not touch `source/r`.
- **Whether hollow directories carry a `00LOCK`.** R's installer
  leaves one; a cleaner does not.
- **What was running.** Check for an R process, a session teardown, or
  an install started by any lane in the minutes before the timestamps.
- **The Windows event log** around the timestamp, for `SilentCleanup`
  or an antimalware action, rather than inferring from free space.

### What held through both, and is therefore worth keeping

The two decisions from 2026-09-09 both survived five losses:

**The pin lives outside `%LOCALAPPDATA%`.** `pinlib` was untouched by
the fourth AND the fifth, and StanHeaders still reads 2.32.10 against
rstan 2.32.7. A pin inside the user library would have gone twice.

**The restore is a script, not a recipe to retype under pressure.**
`dev/release/restore-library.R` recovered 129 of 129 at the fourth
loss and 67 of 67 at the fifth. The frmtmb packages it deliberately
skips, because they are not on CRAN, are reinstalled from the tree:
core first, then the extensions, into both the user library and
`rellib-r3`.

Total cost of the fifth loss was under fifteen minutes to a verified
toolchain, against roughly an hour for the first.

### Verify with a fit, not with a version string

A hollow directory still answers `packageVersion()`. It has a
`DESCRIPTION` only if the sweep spared that file, and the check that
matters is whether the package LOADS and FITS. After any restore, run
something that produces a number and compare it, rather than reading a
version:

    gaussian GLMM, seed 1, 200 rows, (1 | g)   logLik -268.5330674
    wiener DDM, seed 2, 400 trials             logLik -128.3133341

Ten digits, so a future restore can be compared rather than eyeballed.

### What is still not fixed

The library remains under `%LOCALAPPDATA%` by the user's decision, so
this will recur. What has changed is that it is cheap, and that the
next session knows the cause is open rather than closed.
