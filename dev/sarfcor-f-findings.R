# Fix round 1: correct dev/sarfcor-findings.md where the review found it
# wrong (F7), and add the fix round's record. Rscript <this>
p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/sarfcor-findings.md"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("Defects of brms 2.23.0 found here (none posted upstream):

1. The fcor definiteness check passes an indefinite `M` (above).",
"Defects of brms 2.23.0 found here (none posted upstream; corrected in
fix round 1, where the review showed 1 and 3 overstated):

1. The fcor definiteness check passes an indefinite `M` (above). Not a
   wrong model: brms's Stan program takes `cholesky_decompose(Mfcor)`
   in transformed data, and the fit fails there (the review's
   `fcorindef.txt`). brms fails late; frmtmb refuses at the start, with
   the same outcome. A note in `dev/upstream-bugs.md`, not an entry.")
rep1("3. For a `W` that is not row-standardized, brms's interval for rho,
   (min eigenvalue, max eigenvalue), reaches past 1 / max eigenvalue,
   where `log1m(rho * eigenW)` is the log of a negative number.
4. sar with weights, se, cens or trunc: code that does not compile.",
"3. For a `W` that is not row-standardized, brms's interval for rho,
   (min eigenvalue, max eigenvalue), reaches past 1 / max eigenvalue,
   where `log1m(rho * eigenW)` is not a number. Stan rejects such a
   proposal, so for a W with min eigenvalue <= -1 <= 1 <= max
   eigenvalue (every 0/1 adjacency) brms's support is frmtmb's
   interval: the cost is rejected proposals and initial values, not a
   different model. A note in `dev/upstream-bugs.md`, not an entry.
4. sar with weights, se, cens or trunc: code that does not compile.

Entries in `dev/upstream-bugs.md`: defect 2 is brms-37, defect 4
brms-38, and fix round 1 found brms-39 (`posterior_epred()` fails on
every error SAR fit).")
rep1("| sar rho interval | (min ev, max ev) | (1 / min ev, 1 / max ev) | brms's is the wrong interval for an unstandardized W (defect 3) and too narrow on the negative side for a row-standardized one |",
"| sar rho interval | (min ev, max ev) | (1 / min ev, 1 / max ev) | for an unstandardized W brms's effective support is the same, through rejected proposals (defect 3); for a row-standardized W brms's is narrower on the negative side, so frmtmb can estimate rho beyond it (fix round 1, F8) |")
rep1("| fcor with indefinite M | accepted (defect 1) | refused | not a covariance |",
"| fcor with indefinite M | passes the R check, fails in Stan (defect 1) | refused | not a covariance |")
rep1("| new data | refused unless `newdata2` | refused; frmtmb has no `newdata2`; `frm_linpred()` works | matches brms's default |",
"| new data | same rows taken; other rows refused unless `newdata2` | same rows taken (fix round 1, F3); other rows refused; frmtmb has no `newdata2`; `frm_linpred()` works | matches brms's default |")
rep1("| `set_prior(class = \"lagsar\")` | flat on its interval, a class of its own | refused by class; the internal parameter takes `class = \"theta\", coef = \"thetaac_1\"` | not built (A6) |",
"| `set_prior(class = \"lagsar\")` | flat on its interval, a class of its own | a class of its own, carried to the internal parameter; frmtmb.sample's default flat on its interval (fix round 1, F2) | the `bounds` slot |")
rep1("- Priors on `lagsar`/`errorsar` by brms's class. brms lists the class
  with a flat prior on its interval; frmtmb refuses `class = \"lagsar\"`
  and takes `class = \"theta\", coef = \"thetaac_1\"`. A class would need a
  slot that maps the internal parameter to the natural one with its log
  Jacobian on the tape; not built. Consequence for frmtmb.sample: it
  samples `thetaac_1` flat on the logit scale of rho, a density
  `1 / ((rho - lo)(hi - rho))` on rho, where brms's is flat on rho.",
"- Priors on `lagsar`/`errorsar` by brms's class: built in fix round 1
  (F2), through a `bounds` slot; see there.")
rep1("- Priors on `lagsar` and `errorsar` by brms's class, and frmtmb.sample's
  default for them, which is flat on the logit of rho where brms's is
  flat on rho (A6). Building it needs a new entry slot (the natural
  parameter and its log Jacobian on the tape).",
"- Priors on `lagsar` and `errorsar`: done in fix round 1 (F2). `ar`
  and `cosy` keep frmtmb.sample's default flat on the internal
  parameter, which is not brms's flat on the natural one; the
  coordinator put that on the backlog.")
txt <- paste0(txt, "\n\n", paste(readLines(
  "C:/Users/adf44/AppData/Local/Temp/12/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/fr1-section.md"),
  collapse = "\n"), "\n")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("ok\n")
