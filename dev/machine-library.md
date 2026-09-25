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

## Where the disk actually goes, measured 2026-09-23

Seven losses have been tied to the low-disk condition without anyone
naming what produces it. Lane `wt-reunc` measured it while free space
fell to 4.15 GB mid-round:

| tree | size |
|---|---|
| `%LOCALAPPDATA%/Temp` | **26.37 GB** |
| all of `source/r` | 2.29 GB |
| the R library tree | 1.79 GB |

There were NO `Rtmp*` directories at all, so this is not R's doing. The
consumers are `Temp/1` at 13.57 GB, which is agent session state, plus
`s3h2g10b` 4.6 GB, `lmdzeox4` 2.27 GB and `DiagOutputDir` 1.02 GB. All
`claude` session temp is 8.45 GB, of which that lane's own session was
53.5 MB: this accumulates across sessions rather than within one run.
Nothing was deleted, because live processes and other sessions hold
state there and one wrong `Remove-Item` is the failure mode this project
keeps paying for. Recorded for whoever decides what to do about it.

A second scan the next day, at 10 GB free, broke `Temp/1` down further:
14 GB in total, of which **7.9 GB is `Temp/1/claude`**, the scratchpad
tree shared by every session on this box, plus several dozen abandoned
`Rtmp*` directories at about 110 MB each. The `*-check` trees under
`source/r` are NOT the cause: the largest is 33 MB. So the standing
consumer is accumulated session scratchpad and R session temp litter,
which no single lane owns and no lane can safely clear, because live
processes and other sessions hold state in both.

## The seventh loss, 2026-09-22, and the first with a named mechanism

Between 217 and 223 of 410 package directories in
`%LOCALAPPDATA%/R/win-library/4.6` went hollow, found independently by two
lanes. `wt-reunc` counted 217 (134 at 17:06, 83 at 17:07); `wt-skewinit`
counted 223 over 17:05:29 to 17:07:06, hollowed ALPHABETICALLY from
`abind`, which is the shape of a sweep rather than of an install. The two
counts differ because each sampled at a different moment of the sweep.
The evidence, collected before any restore:

- nothing outside `%LOCALAPPDATA%` was touched: `rellib-r3` 8 of 8, the
  lane's own library 2 of 2, `pinlib` intact;
- no `00LOCK` anywhere, so it is not an interrupted install (the shape of
  the fifth loss);
- the canary SURVIVED;
- no R process of that lane was running: its last tier log was 16:29:56;
- free space was 38 GB of about 951 GB, which is 4 percent.

**Storage Sense is on and set to run when free space is low.** Read at
18:40 the same day from
`HKCU:SOFTWARE/Microsoft/Windows/CurrentVersion/StorageSense/Parameters/StoragePolicy`:
`01 = 1` (enabled), `2048 = 0` (run when disk space is low), `04 = 1`
(clean temporary files), `128 = 30` and `256 = 30`. The disk crossed the
low-space threshold, and the loss followed within the window. That is the
first mechanism with a named agent and a setting to point at, rather than
a correlation. `wt-skewinit` also recorded that the scheduled task
`SilentCleanup` last ran at 18:26:12 the same day, result 0.

Two remedies, neither applied without the user: turn Storage Sense off
(or stop it running on low space), and keep the disk above the
threshold. The durable fix is still the one the user declined on
2026-09-09, moving the library out of `%LOCALAPPDATA%`. Restored the same
evening by another session, 223 hollow to 0; verified afterwards by
FITTING, not by counting files: gaussian `(1 | g)` logLik -292.288298118.

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

## The seventh loss, 2026-09-22, and the second with a dated trigger

223 of 410 user-library directories went hollow while lane `wt-correct`
was running `R CMD check --as-cran`. Evidence collected BEFORE the
restore, by `dev/correct-loss7-evidence.R`
(`dev/correct-log/loss7-evidence.txt`), in the order the section above
asks for it:

| | sixth, 17:09 (09-17) | seventh, 17:05 (09-22) |
|---|---|---|
| user-library directories emptied | 136 of 401 | **223 of 410** |
| window | 35 s | **97 s**, 17:05:29 to 17:07:06 |
| `ZZZ-canary.txt` | (not recorded) | **survived**, mtime 2026-09-09 |
| `rellib-r3` | `frmtmb` emptied | **untouched, 8 of 8** |
| `pinlib` | untouched | untouched (StanHeaders 2.32.10) |
| the lane's own library | (none) | **untouched, 2 of 2** |
| `00LOCK` present | (not recorded) | **none** |
| preceded by | a hard kill of every R process | **a session cut while
  `R CMD build` was installing frmtmb to build its vignettes** |

**What was running.** `dev/correct-run-check.ps1` started at 17:02 and
`check.log` stops at "installing the package (it is needed to build
vignettes)", written 17:03. The lane's session hit an account limit
about then and its processes went away; the damage timestamps are the
two minutes after. That is the third suspect again, an install
interrupted mid-flight, and it is the second loss in a row where a
dying `R CMD check` install is the only thing running. What it still
does not explain is the blast radius: that install writes to the check
directory and to no library at all, and 223 unrelated packages went
with it.

**What the seventh adds to the record.**

- The canary survived a sweep that took more than half the library,
  which kills the remaining "spared because recently accessed" reading
  of the fourth loss: nothing about the canary is special here.
- Nothing outside `%LOCALAPPDATA%` was touched this time, including a
  library that was being READ by the dying check.
- No `00LOCK` anywhere, so whatever emptied the directories was not
  R's installer part-way through its own unpack.

**The restore** (`dev/release/restore-library.R`, one call): 223 of 223
recovered, none of them off CRAN, 0 hollow afterwards. Verified by
FITTING, with values a future restore can compare against
(`dev/correct-libcheck.R`, seed 1, n = 200, 10 groups):

    gaussian (1 | g)  logLik -295.602189818
    poisson  (1 | g)  logLik -332.137876329

**What the restore changed, and needs the user.** The user library held
**RTMB 2.0**, installed by the user on 2026-09-21. `RTMB` was among the
223, and the restore reinstalled it from CRAN, whose Windows R 4.6
binary is **RTMB 1.9**. The restore script takes what CRAN has, which is right
wherever CRAN holds the version the machine had. Every one of the 223
was on CRAN, and RTMB is the only one this lane knows was ahead of its
CRAN binary. `frmtmb` itself is
excluded by the script and was not touched. Reinstalling RTMB 2.0 means
`install.packages("RTMB", repos = "https://kaskr.r-universe.dev")` into
the user library, which is the user's call, not a lane's.

## The eighth loss, 2026-09-24, at a critical-battery shutdown

Found at 19:16 after the laptop was plugged back in. Evidence collected
BEFORE the restore:

| | seventh, 17:05 (09-22) | eighth, 18:41 (09-24) |
|---|---|---|
| user-library directories emptied | 223 of 410 | **230 of 410** |
| window | 97 s | **12 s**, 18:40:57 to 18:41:09 |
| order | alphabetical from `abind` | alphabetical, `abind` to `zoo`, with some left whole |
| `ZZZ-canary.txt` | survived | **survived** |
| `rellib-r3` | untouched | **3 of 9 emptied**: `drmTMB`, `frmtmb.eam`, `frmtmb.ode`, 18:41:03 to 18:41:08 |
| lane libraries | untouched | **`predfix-lib` 2 of 5 emptied** (`frmtmb`, `frmtmb.latent`, 18:41:08); `mvprior-lib`, `phase3b-lib`, `simnewdata-lib`, `phase3a-lib`, `pinlib` untouched |
| `00LOCK` present | none | **none** |
| free space | 38 GB | **114 GB** at 19:16 |
| preceded by | a session cut during `R CMD build` | **`Kernel-Power` 524, "Critical Battery Trigger Met", at 18:38:57**; the next System event is the boot at 19:01 |

**What this rules out.** Free space was far above any low-space
threshold, so Storage Sense on low disk is not the trigger this time.
Storage Sense also does not reach `C:/Users/adf44/source/r`, and two
libraries there were emptied inside the same 12 seconds. A lane script
is unlikely: a search of every `.R`, `.ps1` and `.sh` file changed in
the worktrees since noon found no deletion that keeps directories.

**The common factor across the sixth, seventh and eighth losses** is
that processes were killed, not that the disk was low: a hard kill of
every R process (sixth), an agent session cut mid-install (seventh), and
agent sessions and R processes dying at a critical-battery shutdown
(eighth). Which process does the deleting is still unknown. Two lanes
were idle with their check trees already deleted; `wt-mvprior` and
`wt-predfix` had last written logs at 18:37 and 18:33, and the
`wt-phase3b` review at 18:37.

**The restore**, about 19:25 to 19:40:

- `dev/release/restore-library.R`: 230 of 230 recovered from CRAN
  binaries, none off CRAN (`dev/release/restore-8th.log`). It
  reinstalled RTMB as CRAN's Windows binary, 1.9, so RTMB 2.0 went back
  from `https://kaskr.r-universe.dev`, as the user approved after the
  seventh loss.
- `rellib-r3`: `drmTMB` 0.7.0 from CRAN; `frmtmb.eam` 0.10.0 and
  `frmtmb.ode` 0.6.0 from the main checkout at `cad68e21` (the 0.62.0
  release plus a CI change).
- `predfix-lib`: `frmtmb` and `frmtmb.latent` from the `wt-predfix`
  worktree as it stood.
- 0 hollow in all three afterwards. Verified by FITTING with
  `dev/correct-libcheck.R`, identical to the seventh loss's values:
  gaussian `(1 | g)` logLik -295.602189818, poisson -332.137876329.
- StanHeaders in the user library is 2.39.1 as before; the 2.32.10 pin
  in `pinlib` was untouched.
