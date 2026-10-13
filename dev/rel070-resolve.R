# Consolidation edits of the 0.70.0 release tree, applied after the
# textual merge (dev/rel070-merge.sh) and before the NEWS merge
# (dev/rel070-news-merge.R). Each edit checks that its old text is
# there, so a rerun on a fresh merge either applies all of them or
# stops on the first that no longer fits. A file no lane changed is
# restored from its HEAD blob (2884ccbd, 0.69.1) before it is edited,
# so every edit here is applied once however often this runs.
#
#   Rscript dev/rel070-resolve.R
rel <- "C:/Users/adf44/source/r/frmtmb-wt-release"
root <- "C:/Users/adf44/source/r"
setwd(rel)
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("show", paste0("HEAD:", f)), stdout = TRUE)
  wr(sub("\r$", "", x), f)
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

# ---- R/fit.R: the mo_search() tail (moparity against optrap) ---------
# moparity's step 3 (cox_face_leave()) first, then optrap's settle of
# the objective whether or not anything was refitted and its early
# return. moparity's own restore after the probes goes: obj_settle()
# restores the saved state when `best` is `opt` and the state is not at
# it, which is the probes' case.
edit_file("R/fit.R",
  paste0(
"<<<<<<< release
  # 3. a cox() baseline weight left at a face the likelihood rises from
  probed <- FALSE
  for (tm in terms) {
    if (!isTRUE(tm$cox)) next
    probed <- TRUE
    p <- cox_face_leave(obj, best, tm)
    if (is.null(p)) next
    o <- refit(pmin(pmax(p, lower), upper))
    if (better(o)) best <- o
  }
  if (!runs) {
    # the probes moved the objective's state; put it back
    if (probed) obj_state_restore(obj, saved)
    return(opt)
  }
=======
  # settled whether or not anything was refitted: the run that reached
  # `opt` need not have ended at it (obj_settle())
>>>>>>> optrap
  obj_settle(obj, opt, best, saved)
  if (!runs) return(opt)
"),
  paste0(
"  # 3. a cox() baseline weight left at a face the likelihood rises from
  for (tm in terms) {
    if (!isTRUE(tm$cox)) next
    p <- cox_face_leave(obj, best, tm)
    if (is.null(p)) next
    o <- refit(pmin(pmax(p, lower), upper))
    if (better(o)) best <- o
  }
  # settled whether or not anything was refitted: the run that reached
  # `opt` need not have ended at it, and the cox() probes above move
  # the objective's state (obj_settle())
  obj_settle(obj, opt, best, saved)
  if (!runs) return(opt)
"))
# the documentation: moparity's step 3 made the lead's two things
# three, the third for a cox() baseline only
edit_file("R/fit.R",
  paste0("#' `mo_simplex()` has no plateau at the faces, but two things ",
         "remain,\n#' and both are measured on brms_monotonic's ",
         "`ls ~ mo(income) * age`\n#' against its exact maximum ",
         "(dev/optima-findings.md):\n"),
  paste0("#' `mo_simplex()` has no plateau at the faces, but two things ",
         "remain for\n#' a `mo()` term, both measured on brms_monotonic's ",
         "`ls ~ mo(income) *\n#' age` against its exact maximum ",
         "(dev/optima-findings.md), and a third\n#' for a `cox()` ",
         "baseline:\n"))
edit_file("R/fit.R",
  paste0("#' The better optimum is kept. Silent, like ",
         "`escape_stationary()`;\n"),
  paste0("#' The better optimum is kept, and the objective is left at it\n",
         "#' (`obj_settle()`). Silent, like `escape_stationary()`;\n"))

# ---- versions ---------------------------------------------------------
for (f in c("DESCRIPTION", "codemeta.json",
            "extensions/frmtmb.sample/DESCRIPTION",
            "extensions/frmtmb.eam/DESCRIPTION",
            "extensions/frmtmb.eam/NEWS.md")) from_head(f)
edit_file("DESCRIPTION", "Version: 0.69.1\n", "Version: 0.70.0\n")
edit_file("codemeta.json", '"version": "0.69.1",', '"version": "0.70.0",')
edit_file("extensions/frmtmb.sample/DESCRIPTION", "Version: 0.17.0",
          "Version: 0.18.0")
edit_file("extensions/frmtmb.sample/DESCRIPTION",
          "frmtmb (>= 0.69.0)", "frmtmb (>= 0.70.0)")
# frmtmb.eam: test-defects.R requires the prediction warning that lane
# seflags' unit-diagonal metric gives (the 0.68.0 rule: a patch bump
# and a floor); rstan is called by test-sampling.R (since 0.69.0)
edit_file("extensions/frmtmb.eam/DESCRIPTION", "Version: 0.11.2",
          "Version: 0.11.3")
edit_file("extensions/frmtmb.eam/DESCRIPTION", "frmtmb (>= 0.66.0)",
          "frmtmb (>= 0.70.0)")
edit_file("extensions/frmtmb.eam/DESCRIPTION",
          "    Rmpfr,\n    rtdists,\n",
          "    Rmpfr,\n    rstan,\n    rtdists,\n")
x <- rd("extensions/frmtmb.eam/NEWS.md")
stopifnot(x[1] == "# frmtmb.eam 0.11.2")
wr(c(rd("dev/rel070-news-eam.md"), "", x), "extensions/frmtmb.eam/NEWS.md")
cat("edited extensions/frmtmb.eam/NEWS.md\n")

# ---- the release harness reads rellib-r8 ------------------------------
for (f in c("dev/release/run-tests.R", "dev/release/run-check.ps1",
            "dev/release/run-docs.ps1", "dev/warnleak-scan.R",
            "dev/brmsport-run.R", "dev/brmsport-guards.R",
            "dev/ciharden-run1.R", "dev/ciharden-scan.sh")) {
  from_head(f)
  x <- rd(f)
  y <- gsub("rellib-r7", "rellib-r8", x, fixed = TRUE)
  if (!identical(x, y)) { wr(y, f); cat("harness:", f, "\n") }
  stopifnot(any(grepl("rellib-r8", rd(f), fixed = TRUE)))
}

# ---- _pkgdown.yml -----------------------------------------------------
from_head("_pkgdown.yml")
pk_extra <- "dev/rel070-pkgdown-edits.R"
if (file.exists(pk_extra)) source(pk_extra, local = TRUE)

# ---- dev/test-backlog.md ----------------------------------------------
# The HEAD blob with each item this round closed marked in place, then
# one "Filed at the 0.70.0 release" section before "## Reference": the
# consolidation's own text (dev/rel070-backlog-section.md), then lane
# moparity's and lane seflags' sections as the lanes wrote them, one
# heading level down, read from the lanes' trees.
bl <- sub("\r$", "",
          system2("git", c("show", "HEAD:dev/test-backlog.md"),
                  stdout = TRUE))
close_item <- function(x, start, note) {
  i <- which(startsWith(x, start))
  if (length(i) != 1L) stop("backlog item not unique: ", start)
  j <- i + 1L
  # an item's indented code block follows a blank line; the note goes
  # after it, so the sentence that introduces the block stays with it
  # (the consolidation review, m5)
  while (j <= length(x) && (startsWith(x[j], "  ") ||
           (x[j] == "" && j < length(x) && startsWith(x[j + 1L], "  ")))) {
    j <- j + 1L
  }
  wrapped <- strwrap(note, width = 78, initial = "  ", prefix = "  ")
  append(x, wrapped, after = j - 1L)
}
cl <- list(
  c("- **A gaussian random-intercept fit stops short of the maximum near",
    "CLOSED at 0.70.0 (lane optrap): the identity-link mu of ten location-scale families has a unit, the residual MAD; 0 of 40 short at 1e-3, 1 and 1e3."),
  c("- **A prediction along a lost direction is flagged when it is",
    "CLOSED at 0.70.0 (lane seflags): the design alias is taken from the design (se_design_null()); the prediction test reads the unit-diagonal metric. data-helpers:7 holds under OpenBLAS LAPACK."),
  c("- **A restart that raises discards a finite first run** (the",
    "CLOSED at 0.70.0 (lane optrap): a restart that raises keeps the run before it; fuzz seed 20379118 fits under OpenBLAS LAPACK."),
  c("- **A boundary message on an sd the likelihood does not read** (the",
    "CLOSED at 0.70.0 (lane seflags): a boundary verdict needs the likelihood to change when the sd moves up, as far as +10 on the log scale."),
  c("- **One fit in `test-id-kron.R` ends at one of two optima from run to",
    "CLOSED for the test at 0.70.0 (lane optrap): the test no longer depends on which point it reaches; the two points by core type remain."),
  c("- **cox()'s baseline simplex (`sbhaz_raw`) is a softmax and meets the",
    "CLOSED at 0.70.0 (lane moparity): the baseline is in the mo() chart, searched at its faces (cox_face_leave()); frmtmb.sample samples it in softmax coordinates under dirichlet(1)."),
  c("- **`escape_stationary()`'s sibling paths**: the objective's state is",
    "CLOSED at 0.70.0 (lane optrap): every run the optimizer discards leaves the objective at the kept optimum (obj_settle())."),
  c("- **A dense block with several scales still has a `+Inf` or `NaN` log",
    "CLOSED at 0.70.0 (lane moparity): sd_floored() on us, us_t, cs, toep, ar1, hetar1 and gr(cov =, prec =)."),
  c("- **A degree-7 raw polynomial with `(1 | g)` keeps standard errors far",
    "CLOSED at 0.70.0 (lane seflags): all eight SEs within 0.5 percent of the QR reference (second differences of the objective)."),
  c("- **Curved ridges with an eigenvalue above 1e-9 keep silent finite",
    "CLOSED at 0.70.0 (lane seflags): the curvature is measured again on the ridge's floor."),
  c("- **Two results of the merged tree depend on OpenBLAS's LAPACK**,",
    "CLOSED at 0.70.0: the fuzz spec (lane optrap, the restart) and data-helpers:7 (lane seflags, the alias); see dev/round-20261013.md for the release build under OpenBLAS 0.3.26 and 0.3.32 with LAPACK."),
  c("- **`update(prior =)` replaces the stored prior; brms merges** the new",
    "CLOSED at 0.70.0 (lane moparity): merged by slot as brms; an empty table row resets its slot."),
  c("- **`mo()` seed 194** of `dev/optima-mo-study.R` stays 0.0489 below",
    "CLOSED at 0.70.0 (lane optrap): the search tries the next vertex while the coefficient ends on its bound; 400 of 400 at the maximum."),
  c("- **A `simo` prior row is refused by `frm()` and `frm_sample()`**:",
    "CLOSED at 0.70.0 (lane moparity): prior_dirichlet(), classes simo and sbhaz, in frm() (MAP) and frm_sample()."),
  c("- **`as_tmbstan(fit)` on every `mo()` fit**, ML or MAP, samples the",
    "CLOSED at 0.70.0 (lane moparity): contract kept; the message names only the simplexes with no density on the tape."),
  c("- **`mo(x) * f` with a factor and `mo()` in group-level terms**, both",
    "CLOSED at 0.70.0 (lane moparity), in brms 2.22.0's parameterization (brms 2.23.0 is broken here, brms#1828)."),
  c("- **The SE check misses an exactly flat simplex direction that is not",
    "NOT CLOSED at 0.70.0: much reduced, a residue remains (the consolidation review, m1; dev/relrev070-segap.R, 40 seeds by 4 designs): fits silent with every zeta SE finite, 64 of 160 on 0.69.1, 8 of 160 on the release under either BLAS (lane seflags). See \"Filed at the 0.70.0 release\", Open - low."),
  c("- **`parnames.frmtmb_fit()` and `parnames.frmtmb_draws()` decide",
    "CLOSED at 0.70.0 (lane moparity): decided from the dispatching frame; 26 of 26 load orders and routes warn once."),
  c("- **frmtmb.sample's exact `gp()` sampling mixes poorly from either",
    "CLOSED at 0.70.0 (lane moparity): the field is sampled non-centered when its sd and length scales have priors.")
)
for (c2 in cl) bl <- close_item(bl, c2[1], c2[2])
lane_sec <- function(l, head) {
  x <- rd(file.path(root, paste0("frmtmb-wt-", l), "dev/test-backlog.md"))
  i <- which(x == head)
  r <- which(x == "## Reference")
  stopifnot(length(i) == 1L, length(r) == 1L, i < r)
  y <- x[i:(r - 1L)]
  y <- sub("^### ", "#### ", y)
  sub("^## ", "### ", y)
}
sec <- c(rd("dev/rel070-backlog-section.md"), "",
         lane_sec("moparity",
                  "## Filed by lane moparity after punch round 1, 2026-10-07"),
         lane_sec("seflags",
                  "## Filed by lane seflags (round opened 2026-10-07)"))
ref <- which(bl == "## Reference")
stopifnot(length(ref) == 1L)
bl <- append(bl, sec, after = ref - 1L)
bl <- bl[!(c(FALSE, bl[-1] == "" & bl[-length(bl)] == ""))]
wr(bl, "dev/test-backlog.md")
cat("backlog: closed", length(cl), "items in place\n")
cat("RESOLVE DONE\n")
