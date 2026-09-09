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

## The two suspects, neither excluded

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

Evidence that would separate them is not available here. The Task
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

## What this costs a round

Every private library and every reference build is derived from this
one, so a loss invalidates all of them. Rebuild the round's shared
reference library after a restore, and treat any measurement taken
across the boundary as suspect.
