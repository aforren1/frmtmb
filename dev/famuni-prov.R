f <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-findings.md"
md <- readLines(f)
i <- which(md == "### Provenance")
stopifnot(length(i) == 1L)
txt <- c(
"",
"Driver `dev/famuni-final.sh` (then `dev/famuni-ob.sh` and",
"`dev/famuni-sample-reach.sh`), against the final installs listed",
"first. The order of events, all on 2026-10-13:",
"",
"1. 19:48: core, eam and sample installed from the final source. The",
"   full eight-package suite (374 files, one per process), the gated",
"   tier and the drmTMB gate ran against them (`final-suite.txt`,",
"   `final-gated.txt`, `final-drmtmb.txt`). The suite had 2 errors, both",
"   frmtmb.sample tests that call `mo_sample_entries(ent)` without the",
"   `fit` this round added (`test-cox-draws.R`, `test-mo-simplex-draws.R`).",
"   The sample R CMD check found the same two (`check-run1.log`).",
"2. `fit` got a default, `fit = NULL`, the only change after step 1",
"   (`dev/famuni-fix5h.R`). Sample was reinstalled at 20:07",
"   (`final-install-sample2.txt`); the install script also reinstalled",
"   core from its unchanged source (no file under `R/`, `man/` or",
"   `tests/` newer than the 19:48 install, checked with `find`).",
"3. Every file that can reach frmtmb.sample was rerun against the 20:07",
"   build, gated: all frmtmb.sample test files and every core and",
"   extension file that names `frmtmb.sample` or `frm_sample`",
"   (`final-reach.txt`). The family test files ran under OpenBLAS 0.3.26",
"   with LAPACK routed to it (`final-openblas.txt`).",
"4. R CMD check --as-cran of core and sample ran on the final source;",
"   eam's ran after its last change (19:48). The first core run of this",
"   round is `check-run2.log` of the first report's numbering.",
"",
"The suite's gated skips are `test-drmtmb-agreement.R`'s 13, gated on",
"`FRMTMB_DRMTMB_FIT_TESTS`; run with it, that file passes",
"(`final-drmtmb.txt`).",
"")
md <- c(md[seq_len(i)], txt, md[-seq_len(i)])
writeLines(md, f)
