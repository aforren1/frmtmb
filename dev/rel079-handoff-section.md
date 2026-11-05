## The 0.79.0 round, opened 2026-10-29, consolidated 2026-11-04

`dev/round-20261104.md` is this round's record. Each lane has its findings and its
review (`dev/reviews/2026-11-02-<lane>.md`, `2026-11-03-ncgp.md`). Seven lanes off
`d5d0056d` (0.78.0, pushed), merged onto it in the order ncgp, gpgap, refuse, impsamp,
latres, nlfix, predread. The survey's headline moved from 83.6 to 89.0, Wald [84.5,
93.6], and 0 pairs read defect.

Versions: frmtmb **0.79.0**; frmtmb.sample **0.26.0**, frmtmb.learn **0.7.4** and frmtmb.latent **0.7.2**, needing
0.79.0; frmtmb.eam 0.13.1, frmtmb.survival 0.2.1. Spatial 0.3.0, coupling 0.6.1
and ode 0.7.2 keep theirs. The library is `C:/Users/adf44/source/r/rellib-r17`; `rellib-r16`
keeps 0.78.0.

### What to know next time

- `Matrix::Cholesky(env$spHess(...))` returns a stale factor: TMB refills one template and Matrix
  caches the factor in its slot. Clear `@factors` first. A green review ran this as an in-process
  script and did not compare it with the dense solve its table rested on.
- A heredoc in this shell halves backslashes, and an R `writeLines()` to a text file on Windows
  writes CRLF. Write R with the Write tool.
- A `timeout` around `dev/rel079-merge.sh` kills it between removing the copied lane files and
  writing their list. The next run reads every copy as a clash and keeps stale files. Clear
  `git ls-files --others` first, by hand, if that happens.
- The recorded survey and ledger run after the last merge; a merge rerun undoes them.
- `rel079-resolve.R` restores each tracked file no lane changed from its HEAD blob, so an edit there
  applies once on any rerun.
- A runner that copies a script at each job's start needs a copy per job. A job read a shared copy while the
  next job wrote it, and 3 of 84 jobs ran on a lane's library. A `lib:` line from another process cannot show that.
- Read the column before naming a cause. "Past 3000 rows" went into a report with no read of `n`.
- A path past about 260 characters segfaults Rscript here, and a `-e` with `>` under Git Bash writes a file named after what follows.
- First of the next round: `marginaleffects`' standard errors on group-level terms
  (`dev/test-backlog.md`, "Filed at the 0.79.0 release").

