# Consolidation edits of the 0.75.0 release tree, applied after the
# textual merge (dev/rel075-merge.sh) and before the NEWS merge
# (dev/rel075-news-merge.R). Each edit checks that its old text is
# there exactly once (or the count it names), so a rerun on a fresh
# merge either applies all of them or stops on the first that no longer
# fits. A file no lane changed is restored from its HEAD blob (4fd134db)
# before it is edited, so every edit here is applied once however often
# this runs. Every hand resolution of the round is here, each with its
# reason. Adapted from dev/rel074-resolve.R.
#
#   Rscript dev/rel075-resolve.R
rel <- Sys.getenv("REL075_ROOT", "C:/Users/adf44/source/r/frmtmb-wt-rel075")
git_dir <- Sys.getenv("REL075_GIT", rel)
setwd(rel)
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
merged <- readLines("dev/rel075-log/merge-work/files")
copied <- unlist(lapply(list.files("dev/rel075-log", "^copied-.*[.]lst$",
                                   full.names = TRUE), readLines))
restore <- function(f) {
  if (f %in% c(merged, copied)) {
    stop(f, " is a merged or copied file; edit it, do not restore it")
  }
  from_head(f)
}
# Replace the one conflict block of `f` whose two sides contain `ours`
# and `theirs` (fixed strings) with `new` (a character vector of lines,
# or a function of the two sides' lines). Stops unless exactly one
# block matches.
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
  mid <- m[m > a[i] & m < z[i]]
  if (is.function(new)) {
    new <- new(x[seq_len(mid - a[i] - 1L) + a[i]],
               if (z[i] > mid + 1L) x[(mid + 1L):(z[i] - 1L)] else
                 character(0))
  }
  wr(c(x[seq_len(a[i] - 1L)], new, x[-seq_len(z[i])]), f)
  cat("resolved a conflict in", f, "\n")
}

# ---- 1. R/extension-home.R: smallfix x sarfcor ---------------------------
# smallfix rewrote the refusal for a saved fit's missing packages into
# one collector (family_missing_parts(), structure_missing_parts(),
# stop_extension_parts()); sarfcor appended a check of the residual
# correlation terms (autocor_missing_parts(), stop_autocor_missing()) to
# the old require_structure_packages(). The conflict is that append: the
# release side (smallfix's end of stop_extension_parts()) is kept, and
# the autocor parts join smallfix's combined parts in both checks, so
# one refusal names every missing package across family, structures and
# residual terms.
resolve_block("R/extension-home.R",
              'frm_stop(msg, class = "frmtmb_extension_missing"',
              "stop_autocor_missing(autocor_missing_parts(frame))",
              function(ours, theirs) ours)
edit_file("R/extension-home.R",
'#\' Check every family of a fit, a frame, a spec or a single family.
#\'
#\' @noRd
require_model_packages <- function(x) {
  parts <- family_missing_parts(x)
  if (!length(parts)) return(invisible(TRUE))
  stop_extension_parts(c(parts, structure_missing_parts(x)))
}',
'#\' The parts of a fit or a frame that this session cannot use because
#\' of a covariance structure or a residual correlation term
#\' (`autocor_missing_parts()`, R/autocor-seam.R): both are seams an
#\' extension registers, and a refusal names them with the families.
#\'
#\' @noRd
seam_missing_parts <- function(x) {
  c(structure_missing_parts(x), autocor_missing_parts(x))
}

#\' Check every family of a fit, a frame, a spec or a single family.
#\'
#\' @noRd
require_model_packages <- function(x) {
  parts <- family_missing_parts(x)
  if (!length(parts)) return(invisible(TRUE))
  stop_extension_parts(c(parts, seam_missing_parts(x)))
}')
edit_file("R/extension-home.R",
'#\' Check every random-effect block of a fit or a frame. Called where the
#\' families are checked (the objective, the standard-error machinery)
#\' and at the head of the methods that read a block\'s entry directly.
#\'
#\' @noRd
require_structure_packages <- function(x) {
  parts <- structure_missing_parts(x)',
'#\' Check every random-effect block and every residual correlation term
#\' of a fit or a frame. Called where the families are checked (the
#\' objective, the standard-error machinery) and at the head of the
#\' methods that read a block\'s entry directly.
#\'
#\' @noRd
require_structure_packages <- function(x) {
  parts <- seam_missing_parts(x)')
# The "session" sentence names the seam that registers the part: a
# covariance structure's is frmtmb_register_covstruct(), a residual
# correlation term's frmtmb_register_autocor() (sarfcor's own refusal
# said so; the part carries it in `register`)
edit_file("R/extension-home.R",
'  for (p in of("session")) {
    s <- c(s, paste0(p$what, " was registered in an R session, not by a ",
                     "package, and is not registered in this one. Register ",
                     "it again with frmtmb_register_covstruct() before using ",
                     "the fit."))
  }',
'  for (p in of("session")) {
    s <- c(s, paste0(p$what, " was registered in an R session, not by a ",
                     "package, and is not registered in this one. Register ",
                     "it again with ",
                     p[["register"]] %||% "frmtmb_register_covstruct()",
                     " before using the fit."))
  }')
edit_file("R/extension-home.R",
"# A check fires only for the part it guards, so a method that needs
# nothing of an absent package still answers. Once one fires, the
# refusal names every part of the model this session cannot use: a fit
# with a family from one extension and a structure from another named
# one package at a time before, so the user found the second package
# only after installing the first.",
"# A check fires only for the part it guards, so a method that needs
# nothing of an absent package still answers. Once one fires, the
# refusal names every part of the model this session cannot use
# (families, covariance structures, residual correlation terms): a fit
# with a family from one extension and a structure from another named
# one package at a time before, so the user found the second package
# only after installing the first.")
# R/autocor-seam.R (sarfcor): the parts say which seam registers them,
# and the block's own refusal goes through the same builder, so a
# session term, a stale package and a missing one read as a covariance
# structure's do and test-message-uniqueness.R sees one template each
edit_file("R/autocor-seam.R",
'    parts[[length(parts) + 1L]] <- list(
      kind = if (is.null(pkg)) "session" else if (!usable) "missing" else
        "stale",
      pkg = pkg, what = paste0("The ", nm, "() term of this model"))',
'    parts[[length(parts) + 1L]] <- list(
      kind = if (is.null(pkg)) "session" else if (!usable) "missing" else
        "stale",
      pkg = pkg, what = paste0("The ", nm, "() term of this model"),
      register = "frmtmb_register_autocor()")')
edit_file("R/autocor-seam.R",
'Kept apart from the
#\' covariance structures\' check so that a combined refusal can list
#\' both.',
'Kept apart from the
#\' covariance structures\' parts so that one refusal lists both
#\' (`seam_missing_parts()`, R/extension-home.R).')
edit_file("R/autocor-seam.R",
'#\' Stop on the first part `autocor_missing_parts()` reports, in the
#\' words the covariance structures\' refusal uses.
#\'
#\' @noRd
stop_autocor_missing <- function(parts) {
  if (!length(parts)) return(invisible(TRUE))
  p <- parts[[1L]]
  if (identical(p$kind, "missing")) require_extension(p$pkg, p$what)
  if (identical(p$kind, "session")) {
    frm_stop(p$what, " was registered in an R session, not by a package, ",
             "and is not registered in this one. Register it again with ",
             "frmtmb_register_autocor() before using the fit",
             class = "frmtmb_extension_missing", call. = FALSE)
  }
  frm_stop(p$what, " is from the extension package ", p$pkg, ", which ",
           "loads but does not register it. Install the version of ",
           p$pkg, " the fit was made with",
           class = "frmtmb_extension_missing", call. = FALSE)
}

#\' A block\'s entry, refused by its package\'s name when that package is
#\' not loaded: a saved fit with a registered term, read without the
#\' package.
#\'
#\' @noRd
ac_block_entry <- function(ac) {
  e <- ac_entry(ac[["struct"]])
  if (!is.null(e)) return(e)
  stop_autocor_missing(autocor_missing_parts(list(autocor = list(ac))))
  ac_entry(ac[["struct"]])
}',
'#\' A block\'s entry, refused by its package\'s name when that package is
#\' not loaded: a saved fit with a registered term, read without the
#\' package. The refusal is the one every missing part of a model gets
#\' (`stop_extension_parts()`, R/extension-home.R).
#\'
#\' @noRd
ac_block_entry <- function(ac) {
  e <- ac_entry(ac[["struct"]])
  if (!is.null(e)) return(e)
  parts <- autocor_missing_parts(list(autocor = list(ac)))
  if (length(parts)) stop_extension_parts(parts)
  ac_entry(ac[["struct"]])
}')

# ---- 1b. handlers that can swallow the refusal (smallfix x sarfcor) -----
# dev/smallfix-handlers.R on the merged build (dev/rel075-log/
# handlers-1.txt) lists 2 handlers whose expression can reach the
# refusal. sarfcor's frmtmb_covstruct_info() wraps cs_entry() in a
# handler that read the refusal as "no entry" and replaced its class
# (frmtmb_extension_missing) with a plain frmtmb_error: it now uses
# keep_missing(), as every such handler does, so a query for "car"
# without frmtmb.spatial gets the refusal that names the package and how
# to load it. The other, recover_data.frmtmb_fit() (R/interop.R), is
# pre-existing and kept: it returns every error's message as emmeans's
# own form of refusal (emmeans replaces a raised error by "Perhaps a
# 'data' or 'params' argument is needed"), so the refusal's text
# reaches the user, and require_fit_packages() at its head refuses first
edit_file("R/covstruct-seam.R",
"    e <- tryCatch(cs_entry(name), error = function(err) NULL)",
"    e <- tryCatch(cs_entry(name), error = keep_missing(NULL))")

# Two methods reach a residual term's entry with no fit-level check
# before it, so a saved car() and sar() fit read without frmtmb.spatial
# named the sar() term alone (the consolidation's test, item 6, found
# autocor_matrix(); frmtmb_autocor_pars() has the same shape). Both
# check the whole fit first, as smallfix's methods do
edit_file("R/autocor.R",
"autocor_matrix <- function(fit, resp = NULL) {
  acs <- fit$frame[[\"autocor\"]]
  if (!length(acs)) return(NULL)",
"autocor_matrix <- function(fit, resp = NULL) {
  acs <- fit$frame[[\"autocor\"]]
  if (!length(acs)) return(NULL)
  require_fit_packages(fit)")
edit_file("R/autocor-seam.R",
"                    stringsAsFactors = FALSE)
  for (r in names(fit$frame[[\"autocor\"]] %||% list())) {",
"                    stringsAsFactors = FALSE)
  if (length(fit$frame[[\"autocor\"]])) require_fit_packages(fit)
  for (r in names(fit$frame[[\"autocor\"]] %||% list())) {")

# ---- 2. dev/core-boundary.md: bhaz x sarfcor ------------------------------
# Both lanes rewrote the sentence on the tables of refusals: bhaz took
# bhaz() out of them, sarfcor took sar() and fcor() out. Both facts, in
# the present tense of a reference document
resolve_block("dev/core-boundary.md", "(`custom`, `sar()`, `fcor()`)",
              "(`custom` and `bhaz()` at 0.74.0)", c(
"(`custom` since 0.75.0), each with its reason. `bhaz()` is in",
"`extension_aterm_homes`, and frmtmb.survival refuses `bhaz(gr = )` by",
"name. `fcor()` is in core, and `sar()` is in `frmtmb.spatial`",
"(`extension_special_homes`)."))
# sarfcor wrote "(after 0.74.0)" where a lane cannot write the number
edit_file("dev/core-boundary.md", "(after 0.74.0)", "(0.75.0)", n = 4L)
edit_file("dev/core-boundary.md",
"| `frmtmb_register_aterm(option =)` | an addition term",
"| `frmtmb_register_aterm(option =)` (0.75.0) | an addition term")
edit_file("dev/core-boundary.md",
"A saved fit is a second way to meet a moved feature. Since 0.74.0 each
fit records the package of its family and of each covariance
structure. Read where that package is not installed, the fit stops
with the package's name (class `frmtmb_extension_missing`) at the
first method that needs the package (`R/extension-home.R`).",
"A saved fit is a second way to meet a moved feature. Since 0.74.0 each
fit records the package of its family and of each covariance
structure, and since 0.75.0 of each residual correlation term. Read
where such a package is not installed, the fit stops with the
package's name (class `frmtmb_extension_missing`) at the first method
that needs it (`R/extension-home.R`). Since 0.75.0 that one refusal
names every package the fit lacks.")
edit_file("dev/core-boundary.md",
"`frmtmb.spline` is retired. The section \"Moves made in 0.74.0\" has the
numbers.",
"`frmtmb.spline` is retired. The section \"Moves made in 0.74.0\" has the
numbers. 0.75.0 (2026-10-24) moved no feature. It added the seams
`frmtmb_register_aterm(option =)` and `frmtmb_register_autocor()`,
through which `frmtmb.survival` supplies brms's `bhaz()` and
`frmtmb.spatial` brms's `sar()`; `fcor()` is in core.")

# ---- 3. dev/upstream-bugs.md: two brms-36 ---------------------------------
# bhaz and sarfcor each appended entries after brms-35. bhaz's
# interval-censored cox() stays brms-36; sarfcor's three become brms-37
# (complex eigenvalues), brms-38 (stanc rejects sar with weights, se,
# cens, trunc) and brms-39 (error SAR posterior_epred()). The common
# trailing "- Status: not reported." the merge took out of both sides
# closes bhaz's entry
renum <- function(x) {
  for (k in c(38, 37, 36)) {
    x <- gsub(paste0("brms-", k, "(?![0-9])"), paste0("brms-", k + 1L), x,
              perl = TRUE)
  }
  x
}
resolve_block("dev/upstream-bugs.md",
              "### brms-36. An interval-censored `cox()` row",
              "### brms-36. `sar()` passes complex eigenvalues",
              function(ours, theirs) {
                stopifnot(sum(grepl("^### brms-", theirs)) == 3L)
                c(ours, "- Status: not reported.", "", renum(theirs))
              })
# every reference to sarfcor's three numbers in the lane's records and
# scripts (copied into this tree); no test, NEWS bullet or package file
# names one (grep of the lanes' changed and new files)
renum_file <- function(f, n) {
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("brms-3[678](?![0-9])", s,
                                      perl = TRUE)))
  if (k != n) stop(f, ": expected ", n, " brms-36..38, found ", k)
  wr(strsplit(renum(s), "\n", fixed = TRUE)[[1]], f)
  cat("renumbered", k, "in", f, "\n")
}
renum_file("dev/sarfcor-findings.md", 8L)
renum_file("dev/sarfcor-f-findings.R", 3L)
renum_file("dev/sarfcor-f7-upstream.R", 5L)
renum_file("dev/reviews/2026-10-24-sarfcor.md", 2L)
# (the renumbering above also moved this sentence's brms-36, which is
# bhaz's and stays)
edit_file("dev/sarfcor-findings.md",
"NULL). The bhaz lane added a brms-37 in its tree, so consolidation
  may renumber. Nothing posted.",
"NULL). The bhaz lane added a brms-36 in its tree, so the
  consolidation renumbered these three to brms-37, brms-38 and brms-39
  (`dev/rel075-resolve.R`). Nothing posted.")
s <- paste(rd("dev/upstream-bugs.md"), collapse = "\n")
h <- regmatches(s, gregexpr("(?m)^### brms-[0-9]+", s, perl = TRUE))[[1]]
stopifnot(!anyDuplicated(h), identical(tail(h, 4L), paste0("### brms-",
                                                          36:39)))
cat("dev/upstream-bugs.md: brms entries", length(h), ", last",
    tail(h, 1L), ", no number twice\n")

# ---- 4. versions and floors -----------------------------------------------
ver <- function(f, from, to) edit_file(f, paste0("Version: ", from),
                                       paste0("Version: ", to))
dcore <- "DESCRIPTION"
restore(dcore)
ver(dcore, "0.74.0", "0.75.0")
cm <- "codemeta.json"
restore(cm)
# no lane changed a dependency of core, so only the version moves
edit_file(cm, '"version": "0.74.0",', '"version": "0.75.0",')
# frmtmb.survival: .onLoad() calls frmtmb_register_aterm(option =)
dsv <- "extensions/frmtmb.survival/DESCRIPTION"
restore(dsv)
ver(dsv, "0.1.0", "0.2.0")
edit_file(dsv, "frmtmb (>= 0.74.0)", "frmtmb (>= 0.75.0)")
# frmtmb.spatial: imports frmtmb_register_autocor()
dsp <- "extensions/frmtmb.spatial/DESCRIPTION"
ver(dsp, "0.1.0", "0.2.0")
edit_file(dsp, "frmtmb (>= 0.74.0)", "frmtmb (>= 0.75.0)")
# frmtmb.sample: calls frmtmb_covstruct_info() and frmtmb_autocor_pars();
# its new tests read the names 0.75.0 gives. Minor: its NEWS has
# breaking changes (draw names, the order of fixef()) and new features
# (the lagsar column and default prior)
dsa <- "extensions/frmtmb.sample/DESCRIPTION"
restore(dsa)
ver(dsa, "0.21.1", "0.22.0")
edit_file(dsa, "frmtmb (>= 0.74.0)", "frmtmb (>= 0.75.0)")
# frmtmb.latent: a patch (tests only); its test of the predict(newdata =)
# text needs 0.75.0
dla <- "extensions/frmtmb.latent/DESCRIPTION"
restore(dla)
ver(dla, "0.7.0", "0.7.1")
edit_file(dla, "frmtmb (>= 0.74.0)", "frmtmb (>= 0.75.0)")

# ---- 5. the release harness reads rellib-r13 -------------------------------
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-guards.R",
            "dev/brmsport-run.R", "dev/ciharden-run1.R",
            "dev/release/run-check.ps1", "dev/release/run-docs.ps1",
            "dev/release/run-tests.R", "dev/warnleak-scan.R",
            "dev/release/run-session.R", "dev/brmsport-bin2-all.sh",
            "dev/brmsport-bin2-record.sh", "dev/ciharden-scan.sh")) {
  if (!f %in% merged) restore(f)
  s <- paste(rd(f), collapse = "\n")
  k <- lengths(regmatches(s, gregexpr("rellib-r12", s, fixed = TRUE)))
  stopifnot(k >= 1L)
  wr(strsplit(gsub("rellib-r12", "rellib-r13", s, fixed = TRUE), "\n",
              fixed = TRUE)[[1]], f)
  cat("edited", f, "(", k, ")\n")
}

# ---- 6. integration tests the lanes could not write ------------------------
# smallfix's combined refusal and sarfcor's residual term meet only in
# the merged tree: a fit with a family, a structure and a residual term
# from three absent packages, and a saved car() and sar() fit, each
# named in one refusal (dev/rel075-mixed-sar-test.R holds the blocks;
# cox() and sar() cannot be in one model, dev/rel075-mixed-probe.R)
blk <- "dev/rel075-mixed-sar-test.R"
if (file.exists(blk)) {
  f <- "tests/testthat/test-extension-home-mixed.R"
  x <- rd(f)
  stopifnot(!any(grepl("cox(), car() and sar()", x, fixed = TRUE)))
  wr(c(x, "", rd(blk)), f)
  cat("appended", blk, "to", f, "\n")
} else cat("NOTE: no", blk, "yet\n")

# ---- 6b. release fixes ---------------------------------------------------
# the release review (dev/reviews/2026-10-24-release.md, FIX) and the
# coordinator's instructions; the fixes and their reasons are
# dev/rel075-release-fix.R, sourced here so a merge rerun applies them
source("dev/rel075-release-fix.R", local = TRUE)

# ---- 7. dev/test-backlog.md and dev/round-handoff.md ----------------------
bl <- "dev/test-backlog.md"
restore(bl)
bsec <- "dev/rel075-backlog-section.md"
if (file.exists(bsec)) {
  edit_file(bl, "\n## Reference\n",
            paste0("\n", paste(rd(bsec), collapse = "\n"),
                   "\n\n## Reference\n"))
} else cat("NOTE: no", bsec, "yet\n")
hf <- "dev/round-handoff.md"
restore(hf)
hsec <- "dev/rel075-handoff-section.md"
if (file.exists(hsec)) {
  edit_file(hf,
"0.72.0 release, 2026-10-21 at the 0.73.0 release and 2026-10-23 at the
0.74.0 release. Read this, then",
"0.72.0 release, 2026-10-21 at the 0.73.0 release, 2026-10-23 at the
0.74.0 release and 2026-10-24 at the 0.75.0 release. Read this, then")
  edit_file(hf, "\n## The 0.74.0 round, ",
            paste0("\n", paste(rd(hsec), collapse = "\n"),
                   "\n\n## The 0.74.0 round, "))
} else cat("NOTE: no", hsec, "yet\n")
cat("RESOLVE DONE\n")
