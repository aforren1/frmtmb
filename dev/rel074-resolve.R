# Consolidation edits of the 0.74.0 release tree, applied after the
# textual merge (dev/rel074-merge.sh) and before the NEWS merge
# (dev/rel074-news-merge.R). Each edit checks that its old text is
# there exactly once, so a rerun on a fresh merge either applies all of
# them or stops on the first that no longer fits. A file no lane
# changed is restored from its HEAD blob (19b162b9, 0.73.1) before it
# is edited, so every edit here is applied once however often this
# runs. Adapted from dev/rel073-resolve.R.
#
#   Rscript dev/rel074-resolve.R
rel <- Sys.getenv("REL074_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-release")
git_dir <- Sys.getenv("REL074_GIT", rel)
setwd(rel)
lane_dir <- function(l) file.path("C:/Users/adf44/source/r",
                                  paste0("frmtmb-wt-", l))
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("-C", git_dir, "show", paste0("HEAD:", f)),
               stdout = TRUE)
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
merged <- readLines("dev/rel074-log/merge-work/files")
restore <- function(f) {
  if (f %in% merged) stop(f, " is a merged file; edit it, do not restore it")
  from_head(f)
}
# Replace the one conflict block of `f` whose two sides contain `ours`
# and `theirs` (fixed strings) with `new` (a character vector of lines).
# Stops unless exactly one block matches.
resolve_block <- function(f, ours, theirs, new) {
  x <- rd(f)
  a <- which(startsWith(x, "<<<<<<< "))
  m <- which(x == "=======")
  z <- which(startsWith(x, ">>>>>>> "))
  stopifnot(length(a) == length(z))
  hit <- integer(0)
  for (i in seq_along(a)) {
    mid <- m[m > a[i] & m < z[i]]
    stopifnot(length(mid) == 1L)
    o <- paste(x[(a[i] + 1L):(mid - 1L)], collapse = "\n")
    t <- if (z[i] > mid + 1L) paste(x[(mid + 1L):(z[i] - 1L)],
                                    collapse = "\n") else ""
    if (grepl(ours, o, fixed = TRUE) && grepl(theirs, t, fixed = TRUE)) {
      hit <- c(hit, i)
    }
  }
  if (length(hit) != 1L) stop(f, ": expected 1 conflict block, found ",
                              length(hit), ": ", ours)
  i <- hit
  wr(c(x[seq_len(a[i] - 1L)], new, x[-seq_len(z[i])]), f)
  cat("resolved a conflict in", f, "\n")
}

# ---- R/extension-home.R: one copy ----------------------------------------
# Both lanes created it. Lane spatial's copy is lane splitsurv's file
# (require_extension(), require_model_packages()) plus
# require_structure_packages(), as the spatial review measured (A3).
# Checked here line for line, then spatial's copy is taken
eh <- "R/extension-home.R"
ss <- rd(file.path(lane_dir("splitsurv"), eh))
sp <- rd(file.path(lane_dir("spatial"), eh))
stopifnot(length(sp) > length(ss), identical(sp[seq_along(ss)], ss),
          any(grepl("^require_structure_packages <- function", sp)),
          !any(grepl("^require_structure_packages", ss)))
wr(sp, eh)
cat("R/extension-home.R: spatial's copy (splitsurv's", length(ss),
    "lines plus", length(sp) - length(ss), "), one file\n")

# ---- splitsurv x spatial: the absent-package checks ----------------------
# splitsurv checks the family's package and spatial each block's
# structure package at the same two entry points: both run
resolve_block("R/fit.R", "require_model_packages(fit)",
              "require_structure_packages(fit)",
              c("  require_model_packages(fit)",
                "  require_structure_packages(fit)"))
resolve_block("R/objective.R", "require_model_packages(frame)",
              "require_structure_packages(frame)",
              c("  require_model_packages(frame)",
                "  require_structure_packages(frame)"))

# ---- R/compat.R: two adjacent rule deletions -----------------------------
# splitsurv cut the cens() x cox rule (frmtmb.survival registers it) and
# mvnlatent the cens() and trunc() x mixture_mvn rules (frmtmb.latent
# registers them): all three go
resolve_block("R/compat.R", 'r("cens()", "mixture_mvn", "refused",',
              'r("cens()", "cox", "works",', character(0))

# ---- R/structure.R: two new sections of ?frmtmb-extension-api ------------
# mvnlatent's "The mixture component interface" and spatial's "How to
# add a covariance structure", inserted at the same place: both, in lane
# order
local({
  x <- rd("R/structure.R")
  a <- which(startsWith(x, "<<<<<<< "))
  m <- which(x == "=======")
  z <- which(startsWith(x, ">>>>>>> "))
  stopifnot(length(a) == 1L, length(m) == 1L, length(z) == 1L,
            grepl("The mixture component interface", x[a + 1L]),
            grepl("How to add a covariance structure", x[m + 1L]))
  wr(c(x[seq_len(a - 1L)], x[(a + 1L):(m - 1L)], "#'",
       x[(m + 1L):(z - 1L)], x[-seq_len(z)]), "R/structure.R")
  cat("resolved R/structure.R: both sections\n")
})
# mvnlatent's F2 left the two seams of an unstructured block to link
# each other once both were in one tree
edit_file("R/covstruct.R",
"#' @seealso [frmtmb_family()] for `extra_pars`, and [frmtmb-robust-dpars]
#'   for the other helpers a density can call while it is taped.",
"#' @seealso [frmtmb_family()] for `extra_pars`, [frmtmb-robust-dpars]
#'   for the other helpers a density can call while it is taped, and
#'   [frmtmb_register_covstruct()] for a covariance structure whose
#'   `lpdf` or `from_natural` can use the same map.")
edit_file("R/covstruct-seam.R",
"#'   [frm_compat()] (the structure itself joins the vocabulary as a
#'   feature of kind `\"covstruct\"` when it is registered), and
#'   [frmtmb_register_aterm()] for the seam this one is modeled on.",
"#'   [frm_compat()] (the structure itself joins the vocabulary as a
#'   feature of kind `\"covstruct\"` when it is registered),
#'   [frmtmb_register_aterm()] for the seam this one is modeled on, and
#'   [us_sigma()] and [us_theta()] for the parameterization of an
#'   unstructured block, which an `lpdf` can reuse.")

# ---- package lists: frmtmb.survival and frmtmb.spatial, no frmtmb.spline -
# splitsurv replaced frmtmb.spline by frmtmb.survival in each list, and
# spatial appended frmtmb.spatial to the base list
resolve_block(".github/workflows/pkgdown.yaml", '"frmtmb.survival"',
              '"frmtmb.spatial"', c(
'                   "frmtmb.survival", "frmtmb.coupling", "frmtmb.sample",',
'                   "frmtmb.learn",',
'                   "frmtmb.spatial")'))
resolve_block("dev/release/build-docs.R", '"frmtmb.survival"',
              '"frmtmb.spatial"', c(
'                "frmtmb.survival", "frmtmb.coupling", "frmtmb.sample",',
'                "frmtmb.learn",',
'                "frmtmb.spatial")'))
resolve_block("dev/release/run-session.R", '"frmtmb.survival"',
              '"frmtmb.spatial"', c(
'          "frmtmb.learn", "frmtmb.ode", "frmtmb.sample", "frmtmb.survival",',
'          "frmtmb.spatial")'))
resolve_block("dev/release/run-sessions.sh", "frmtmb.survival",
              "frmtmb.spatial", c(
'                           frmtmb.survival',
'                           frmtmb.spatial)'))
edit_file("dev/release/run-sessions.sh",
          "# With no package named, all eight run, so 24 sessions.",
          "# With no package named, all nine run, so 27 sessions.")
resolve_block("tests/testthat/test-compat-extensions.R",
              '"frmtmb.survival"', '"frmtmb.spatial"', c(
'                       "frmtmb.survival",',
'                       "frmtmb.spatial")'))
resolve_block("_pkgdown.yml", "frmtmb.survival", "frmtmb.spatial", c(
"      - text: frmtmb.survival (survival families)",
"        href: https://aforren1.github.io/frmtmb/frmtmb.survival/",
"      - text: frmtmb.spatial (spatial random fields)",
"        href: https://aforren1.github.io/frmtmb/frmtmb.spatial/"))
# splitsurv's lists that spatial did not change (its check and suite
# drivers): frmtmb.spatial joins them
edit_file("dev/release/run-check.ps1",
'foreach ($n in @("frmtmb.eam", "frmtmb.latent", "frmtmb.ode",
                 "frmtmb.survival", "frmtmb.coupling", "frmtmb.sample",
                 "frmtmb.learn")) {',
'foreach ($n in @("frmtmb.eam", "frmtmb.latent", "frmtmb.ode",
                 "frmtmb.survival", "frmtmb.coupling", "frmtmb.sample",
                 "frmtmb.learn", "frmtmb.spatial")) {')
edit_file("dev/release/run-suite.ps1",
'                 "frmtmb.survival", "frmtmb.learn", "frmtmb.latent",
                 "frmtmb.ode")) {',
'                 "frmtmb.survival", "frmtmb.learn", "frmtmb.latent",
                 "frmtmb.ode", "frmtmb.spatial")) {')

# ---- the bhaz() refusal after cox() moved (spatial review, B2) ----------
# lane spatial wrote it against a core that had cox(); after the merge
# cox() is frmtmb.survival's, so the refusal says where it is and how to
# load it, in the review's words. The constant stays the one place the
# pointer lives (test-guardrail1.R reads it)
edit_file("R/brms-vocabulary.R",
'bhaz_family_pointer <- paste0(
  "family = cox(df = 5, degree = 3, intercept = TRUE), which is brms\'s ",
  "bhaz() default")',
'bhaz_family_pointer <- paste0(
  "cox() is in the extension package frmtmb.survival: install it, run ",
  "library(frmtmb.survival), then use family = cox(df = 5, degree = 3, ",
  "intercept = TRUE), which is brms\'s bhaz() default")')
edit_file("R/brms-vocabulary.R",
'    "model, and frmtmb sets the same basis in the family: ",
    bhaz_family_pointer, ". A group-specific baseline, bhaz(gr = ), has ",',
'    "model, and frmtmb sets the same basis in the family. ",
    bhaz_family_pointer, ". A group-specific baseline, bhaz(gr = ), has ",')
edit_file("R/brms-vocabulary.R",
"# Where the bhaz() refusal sends a user: the family that takes brms's
# baseline-hazard basis. The one place to change when cox() moves to an
# extension package (the refusal reads it).",
"# Where the bhaz() refusal sends a user: the family that takes brms's
# baseline-hazard basis, and the extension package that holds it. The
# one place to change when cox() moves again (the refusal reads it).")

# ---- versions and floors ---------------------------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
ver("DESCRIPTION", "0.73.1", "0.74.0")
# codemeta.json: mvnlatent took mclust out; splitsurv took splines2 out
# of Suggests and put gratia in; spatial took fmesher out. codemeta
# lists every Suggests
edit_file("codemeta.json", '"version": "0.73.1",', '"version": "0.74.0",')
cm_pkg <- function(p) paste0(
'    {
      "@type": "SoftwareApplication",
      "identifier": "', p, '",
      "name": "', p, '",
      "provider": {
        "@id": "https://cran.r-project.org",
        "@type": "Organization",
        "name": "Comprehensive R Archive Network (CRAN)",
        "url": "https://cran.r-project.org"
      },
      "sameAs": "https://CRAN.R-project.org/package=', p, '"
    },')
edit_file("codemeta.json", paste0(cm_pkg("fmesher"), "\n"), "")
edit_file("codemeta.json", paste0(cm_pkg("splines2"), "\n"), "")
edit_file("codemeta.json", cm_pkg("glmmTMB"),
          paste0(cm_pkg("glmmTMB"), "\n", cm_pkg("gratia")))
# each new package needs this release's seams: survival
# frmtmb_family(simplex_pars =, aterm_data =), spatial
# frmtmb_register_covstruct(), latent us_sigma() and the fieldterms
# spec slot, and a core without mixture_mvn (measured by
# dev/rel074-floor.sh)
for (p in c("survival", "spatial")) {
  edit_file(file.path("extensions", paste0("frmtmb.", p), "DESCRIPTION"),
            "frmtmb (>= 0.73.1)", "frmtmb (>= 0.74.0)")
}
dl <- "extensions/frmtmb.latent/DESCRIPTION"
ver(dl, "0.6.2", "0.7.0")
edit_file(dl, "frmtmb (>= 0.72.0)", "frmtmb (>= 0.74.0)")
ds <- "extensions/frmtmb.sample/DESCRIPTION"
ver(ds, "0.21.0", "0.21.1")
# its code loads on 0.73.0, but four of its test files need 0.74.0
# (frmtmb:::cs_entry(), frmtmb_family(simplex_pars =)) or
# frmtmb.latent 0.7.0 (mixture_mvn()), measured by dev/rel074-floor.sh;
# a file that fails only on the old core names the floor
edit_file(ds, "frmtmb (>= 0.73.0)", "frmtmb (>= 0.74.0)")
dlearn <- "extensions/frmtmb.learn/DESCRIPTION"
restore(dlearn)
ver(dlearn, "0.7.2", "0.7.3")

# ---- the release harness reads rellib-r12, and the nine suites -----------
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R", "dev/brmsport-bin2-all.sh",
            "dev/brmsport-bin2-record.sh", "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r11", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("rellib-r11", "rellib-r12", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}
cs <- "dev/ciharden-scan.sh"
edit_file(cs, "# Fragility scan: run every test file of all eight suites, one R process",
          "# Fragility scan: run every test file of all nine suites, one R process")
edit_file(cs, "# lines; default every test-*.R of core and the seven extensions.",
          paste0("# lines; default every test-*.R of core and the eight extensions.\n",
                 "# CIHARDEN_UNGATED=true sets NOT_CRAN only, for the plain tier."))
edit_file(cs, "  for e in coupling eam latent learn ode sample spline; do",
          "  for e in coupling eam latent learn ode sample spatial survival; do")
edit_file(cs,
'export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
export FRMTMB_SCALE_TESTS=true',
'export NOT_CRAN=true
# the plain tier is NOT_CRAN alone; the default sets every gate
if [ "${CIHARDEN_UNGATED:-}" = true ]; then
  unset FRMTMB_BRMS_FIT_TESTS FRMTMB_DRMTMB_FIT_TESTS FRMTMB_FUZZ \\
    FRMTMB_SCALE_TESTS
else
  export FRMTMB_BRMS_FIT_TESTS=true
  export FRMTMB_DRMTMB_FIT_TESTS=true
  export FRMTMB_FUZZ=true
  export FRMTMB_SCALE_TESTS=true
fi')

# ---- the ported brms suite with every extension loaded (guardrail 4) -----
# dev/core-boundary.md measures parity on core plus every extension. The
# two runners loaded frmtmb and the suite's own package only; with
# FRMTMB_PORT_EXTENSIONS=all they load every frmtmb.* package of the
# library under test first, and say so on the lib: line
port_ext <- '# guardrail 4 of dev/core-boundary.md: parity is measured on core plus
# every extension, so with FRMTMB_PORT_EXTENSIONS=all every frmtmb.*
# package of the library under test is loaded before the file runs
port_loaded <- character()
if (identical(Sys.getenv("FRMTMB_PORT_EXTENSIONS"), "all")) {
  port_loaded <- grep("^frmtmb[.]", rownames(utils::installed.packages(lib)),
                      value = TRUE)
  for (q in port_loaded) suppressMessages(loadNamespace(q))
}
cat("extensions loaded:", if (length(port_loaded)) port_loaded else "none",
    "\\n")'
for (f in c("dev/brmsport-run.R", "dev/brmsport-bin2-run.R")) {
  edit_file(f,
'suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))',
paste0('suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
', port_ext))
}

# ---- release fix: the saved-fit arms under R CMD check -------------------
# frmtmb.survival's test-saved-fits.R could not build its library stack
# under R CMD check (1 ERROR); the fix and its reason are
# dev/rel074-release-fix.R, approved by the coordinator as a release fix
source("dev/rel074-release-fix.R", local = TRUE)

# ---- dev/test-backlog.md and dev/round-handoff.md ------------------------
bl <- "dev/test-backlog.md"
restore(bl)
bsec <- "dev/rel074-backlog-section.md"
if (file.exists(bsec)) {
  edit_file(bl, "\n## Reference\n",
            paste0("\n", paste(rd(bsec), collapse = "\n"),
                   "\n\n## Reference\n"))
} else cat("NOTE: no", bsec, "yet\n")
hf <- "dev/round-handoff.md"
restore(hf)
hsec <- "dev/rel074-handoff-section.md"
if (file.exists(hsec)) {
  edit_file(hf,
"0.72.0 release and 2026-10-21 at the 0.73.0 release. Read this, then",
"0.72.0 release, 2026-10-21 at the 0.73.0 release and 2026-10-23 at the
0.74.0 release. Read this, then")
  edit_file(hf, "\n## The 0.73.0 round, 2026-10-20 to 2026-10-21\n",
            paste0("\n", paste(rd(hsec), collapse = "\n"),
                   "\n\n## The 0.73.0 round, 2026-10-20 to 2026-10-21\n"))
} else cat("NOTE: no", hsec, "yet\n")
cat("RESOLVE DONE\n")
