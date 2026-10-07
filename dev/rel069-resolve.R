# Consolidation edits of the 0.69.0 release tree, applied after the
# textual merge (dev/rel069-merge.sh) and before the NEWS merge
# (dev/rel069-news-merge.R). Each edit checks that its old text is
# there, so a rerun on a fresh merge either applies all of them or
# stops on the first that no longer fits.
#
#   Rscript dev/rel069-resolve.R
rel <- "C:/Users/adf44/source/r/frmtmb-wt-release"
setwd(rel)
rd <- function(f) readLines(f, warn = FALSE)
wr <- function(x, f) {
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
edit_file <- function(f, old, new, n = 1L) {
  s <- paste(rd(f), collapse = "\n")
  hits <- gregexpr(old, s, fixed = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ",
                   substr(old, 1, 70))
  s <- gsub(old, new, s, fixed = TRUE)
  wr(strsplit(s, "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}

# ---- versions --------------------------------------------------------
edit_file("DESCRIPTION", "Version: 0.68.1.9000", "Version: 0.69.0")
edit_file("codemeta.json", '"version": "0.68.1",', '"version": "0.69.0",')
edit_file("extensions/frmtmb.sample/DESCRIPTION", "Version: 0.16.0",
          "Version: 0.17.0")
edit_file("extensions/frmtmb.sample/DESCRIPTION",
          "frmtmb (>= 0.68.1.9000)", "frmtmb (>= 0.69.0)")

# ---- a comment lane setier's tier order made stale -------------------
# se_flat_tol is read by se_tier3_could_act() before either inverse is
# kept, not "only after both inverses have failed"; the line break is
# lane optima's
edit_file("R/se-check.R",
  paste0("#' along. Only read after both inverses have failed. Exact ",
         "ridges give\n#' 1e-16 to 2e-16 there and a mo() simplex ",
         "coordinate the likelihood\n#' does not read gives an exactly ",
         "zero diagonal or a decoupled row (a\n#' softmax weight run to ",
         "0 did too, before lane optima), so the value is\n#' not ",
         "delicate; it\n#' is the threshold"),
  # 80 columns: the line after "threshold" goes on to the next line
  # (the consolidation review, m1)
  paste0("#' along. se_tier3_could_act() reads it before either inverse ",
         "is kept.\n#' Exact ridges give 1e-16 to 2e-16 there and a mo() ",
         "simplex coordinate\n#' the likelihood does not read gives an ",
         "exactly zero diagonal or a\n#' decoupled row (a softmax weight ",
         "run to 0 did too, before lane\n#' optima), so the value is not ",
         "delicate; it is the threshold\n#'"))

# ---- dev/lane-rules.md: OpenBLAS passes route LAPACK (review m7) -----
x <- rd("dev/lane-rules.md")
if (!any(x == "## Added after the 0.69.0 round")) {
  i <- which(x == "## Cost, which is a real constraint")
  stopifnot(length(i) == 1L)
  add <- c(
    "## Added after the 0.69.0 round", "",
    "- **An OpenBLAS pass is `dev/ciharden-openblas.sh <ver> lapack`**, so",
    "  that LAPACK is routed to OpenBLAS too, as on the ubuntu runners.",
    "  `dev/optima-openblas.sh` and `dev/cifix-openblas.sh` route BLAS only",
    "  and keep R's reference LAPACK; lane optima's pass with the first",
    "  could not show the two failures the 0.69.0 merged tree has under",
    "  OpenBLAS LAPACK (`dev/reviews/2026-10-07-release.md`, B1).", "")
  wr(append(x, add, after = i - 2L), "dev/lane-rules.md")
  cat("edited dev/lane-rules.md\n")
}

# ---- _pkgdown.yml, which no lane edited ------------------------------
# Lane surface gave core four new reference topics and moved the draws
# pages of pp_mixture() and plot() in frmtmb.sample, and check_pkgdown()
# stopped on the index. The redirects of the retired core pages
# reference/pp_mixture.html and reference/stancode.html would now
# overwrite live core pages, so they go; the others point where the
# topic lives now. Restored from the base blob first, so a rerun applies
# the same edits.
pk <- system2("git", c("show", "4f5ea39f:_pkgdown.yml"), stdout = TRUE)
wr(sub("\r$", "", pk), "_pkgdown.yml")
edit_file("_pkgdown.yml", "  - frm_multiple\n",
          "  - frm_multiple\n  - frm_multiple-methods\n")
edit_file("_pkgdown.yml", "  - bayes_R2\n",
          "  - bayes_R2\n  - add_criterion\n  - pp_mixture\n  - stancode\n")
edit_file("_pkgdown.yml",
  paste0('  - ["reference/pp_mixture.frmtmb_draws.html", ',
         '"frmtmb.sample/reference/pp_mixture.html"]\n',
         '  - ["reference/pp_mixture.html", ',
         '"frmtmb.sample/reference/pp_mixture.html"]\n'),
  paste0('  - ["reference/pp_mixture.frmtmb_draws.html", ',
         '"frmtmb.sample/reference/pp_mixture.frmtmb_draws.html"]\n'))
edit_file("_pkgdown.yml",
  paste0('  - ["reference/stancode.html", ',
         '"frmtmb.sample/reference/frmtmb-draws-refusals.html"]\n'), "")
edit_file("_pkgdown.yml",
  paste0('  - ["reference/standata.html", ',
         '"frmtmb.sample/reference/frmtmb-draws-refusals.html"]'),
  '  - ["reference/standata.html", "reference/stancode.html"]')
edit_file("_pkgdown.yml",
  paste0('  - ["reference/plot.frmtmb_draws.html", ',
         '"frmtmb.sample/reference/frmtmb-draws-refusals.html"]'),
  paste0('  - ["reference/plot.frmtmb_draws.html", ',
         '"frmtmb.sample/reference/plot.frmtmb_draws.html"]'))

# ---- dev/lane-rules.md: the active-binding count, rechecked ----------
# on the merged build (rellib-r7), with and without brms attached: 33 in
# frmtmb and 24 in frmtmb.sample, lane surface's count (its review, m2)
edit_file("dev/lane-rules.md",
  paste0("57 exported names after lane\n  surface (2026-10-07), 33 in ",
         "frmtmb and 24 in frmtmb.sample (56 at"),
  paste0("57 exported names at 0.69.0\n  (lane surface, rechecked on the ",
         "merged build), 33 in frmtmb and 24\n  in frmtmb.sample (56 at"))

# ---- the release harness reads rellib-r7 (idempotent) ----------------
for (f in c("dev/release/run-tests.R", "dev/release/run-check.ps1",
            "dev/release/run-docs.ps1", "dev/warnleak-scan.R",
            "dev/brmsport-run.R", "dev/brmsport-guards.R",
            "dev/ciharden-run1.R", "dev/ciharden-scan.sh")) {
  x <- rd(f)
  y <- gsub("rellib-r6", "rellib-r7", x, fixed = TRUE)
  if (!identical(x, y)) { wr(y, f); cat("harness:", f, "\n") }
  stopifnot(any(grepl("rellib-r7", rd(f), fixed = TRUE)))
}

# ---- dev/test-backlog.md ---------------------------------------------
# The base blob with each closed item marked in place, the four lanes'
# own sections replaced by one "Filed at the 0.69.0 release" section
# (dev/rel069-backlog-section.md) before "## Reference".
bl <- system2("git", c("show", "4f5ea39f:dev/test-backlog.md"),
              stdout = TRUE)
bl <- sub("\r$", "", bl)
close_item <- function(x, start, note) {
  i <- which(startsWith(x, start))
  if (length(i) != 1L) stop("backlog item not unique: ", start)
  j <- i + 1L
  while (j <= length(x) && startsWith(x[j], "  ")) j <- j + 1L
  wrapped <- strwrap(note, width = 78, initial = "  ", prefix = "  ")
  append(x, wrapped, after = j - 1L)
}
cl <- list(
  c("- **vigport defect 2: `stancode()`",
    "CLOSED at 0.69.0 (lane surface): refusals by name on a fit; pp_mixture() computed at the estimates."),
  c("- **vigport defect 4: `update()` keeps a prior",
    "CLOSED at 0.69.0 (lane surface): update() drops it as brms does, with a message."),
  c("- **vigport defects 8 and 9, owned by lane nanse**",
    "CLOSED at 0.68.0 (lane nanse); the mo() SE losses are gone at 0.69.0 (lane optima, 0 of 400 fits)."),
  c("- **fixes, consolidation item 1: a linear ridge",
    "CLOSED at 0.68.0 (lane nanse); measured closed by lane setier (dev/nanse-item3.R)."),
  c("- **`fixef()` and `summary()` report NaN standard errors without a",
    "CLOSED at 0.68.0 (lane nanse); measured closed by lane setier (dev/setier-sx.R, 0 of 140)."),
  c("- **`fitted()` on an ordinal fit whose `disc` predictor has no fixed",
    "CLOSED at 0.69.0 (lane surface)."),
  c("- **`frm_sample(fit)` on an exact `y ~ gp(x)` fit does not move**",
    "CLOSED at 0.69.0 (lane surface): the field's log density was +Inf at an underflowed sd; the variance is floored."),
  c("- **`frm_sample()` repeats a warning the frame build gives**",
    "CLOSED at 0.69.0 (lane surface): not a defect; the six warnings were six different ones, and a guard test pins the one."),
  c("- **`cs()` on a cumulative component of an ordinal mixture often stops",
    "CLOSED at 0.69.0 (lane optima): no NaN-gradient error on 40 fits; the non-convergence that remains is correct."),
  c("- **`mo()` point estimates are not always the maximum.**",
    "CLOSED at 0.69.0 (lane optima): 199 of 200 seeds at the exact maximum to 1e-6."),
  c("- **Unseeded data in the ported brms suite**",
    "CLOSED at 0.69.0 (lane ciharden): each block is seeded."),
  c("- **`(cs(1) | g)` with brms attached after frmtmb dies",
    "CLOSED at 0.69.0 (lane surface), with the categorical message."),
  c("- **A sum-to-zero component of an `order = \"none\"` mixture lists no",
    "CLOSED at 0.69.0 as a documented divergence (lane surface, item 7)."),
  c("- **`test-perf.R` \"fit time grows with n",
    "CLOSED at 0.69.0 (lane ciharden): node and byte counts replace the clock."),
  c("- **vigport defect 3: `plot()` of a fit suggests `x`",
    "CLOSED at 0.69.0 (lane surface)."),
  c("- **vigport defect 5: `fixef()` has no `frm_multiple()` method**",
    "CLOSED at 0.69.0 (lane surface): Rubin's rules, brms's names."),
  c("- **The probit's log-odds form underflows past",
    "CLOSED at 0.69.0 (lane optima), with cloglog and softit."),
  c("- **The SE check's tier 1 accepts an inverse built on a Hessian row of",
    "CLOSED at 0.69.0 (lane setier): tier 3 is asked first."),
  c("- **A nonlinear ridge on a fit with random effects reports a huge SE",
    "CLOSED at 0.69.0 (lane setier)."),
  c("- **Separation is not named at the default budget**",
    "CLOSED at 0.69.0 (lane setier): named whatever the optimizer's code."),
  c("- **`ranef(condVar = TRUE)` still reads sdreport()'s own",
    "CLOSED at 0.69.0 (lane setier)."),
  c("- **`test-cumulative-cs.R:132` and `test-ordinal-mixture.R:751` are",
    "CLOSED at 0.69.0 (lane ciharden)."),
  c("- **frmtmb.sample lets two warnings escape on the Ubuntu runner**",
    "CLOSED at 0.69.0 (lane ciharden)."),
  c("- **Core's gp() position key is 15 significant digits, not exact**",
    "CLOSED at 0.69.0 by decision: the key stays brms's 15-digit key, which brms 2.23.0's new-data rule uses (lane ciharden, review B2)."),
  c("- **`gp_krig_cov()` forms `outer(w, w)` whole",
    "CLOSED at 0.69.0 (lane ciharden).")
)
for (c2 in cl) bl <- close_item(bl, c2[1], c2[2])
sec <- rd("dev/rel069-backlog-section.md")
ref <- which(bl == "## Reference")
stopifnot(length(ref) == 1L)
bl <- append(bl, sec, after = ref - 1L)
bl <- bl[!(c(FALSE, bl[-1] == "" & bl[-length(bl)] == ""))]
wr(bl, "dev/test-backlog.md")
cat("backlog: closed", length(cl), "items in place\n")
cat("RESOLVE DONE\n")
