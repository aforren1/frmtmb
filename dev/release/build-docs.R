#!/usr/bin/env Rscript
#
# Build the eight pkgdown sites of this monorepo into one tree: the core
# site at the root of the destination and each extension in a
# subdirectory named by the `destination` in its own `_pkgdown.yml`.
#
# This replaces a PowerShell-only driver so that the published site does
# not depend on one machine's state. The same file runs in
# .github/workflows/pkgdown.yaml and on a developer box.
#
# Usage:
#   Rscript dev/release/build-docs.R [options]
#
#   --dest=PATH        destination root. Default: the FRMTMB_DOCS_DEST
#                      environment variable, else `docs` under the
#                      repository root, which is the layout the
#                      committed site has today.
#   --pkgs=a,b,c       build only these packages, by package name.
#                      A partial build still has to satisfy
#                      --min-sites, so it fails unless you lower that
#                      deliberately. Deploying a partial tree would
#                      publish a site with dead navbar links.
#   --check-links      run the whole-tree checks, the search index and
#                      the navbar hrefs, even on a partial build. They
#                      are on by default for a full build and off for a
#                      partial one, because a partial tree has dead
#                      links by construction.
#   --min-sites=N      how many sites the run must produce. Default 8.
#   --no-clean         keep whatever is already in the destination.
#                      The default removes it, which is how a renamed
#                      package's stale directory disappears.
#   --log=PATH         where the per-site build output goes. Default
#                      dev/release/docs.log under the repository root.
#   --require-articles fail before building if a package that an article
#                      gates itself on is not installed. A build machine
#                      without it renders the article with every gated
#                      chunk dropped, which looks like a success and
#                      publishes a shorter page. CI passes this; a local
#                      build leaves it off and gets a warning instead.
#   --site=PATH        internal. Build exactly one site, in this
#                      process, and exit. The driver spawns one child
#                      per site so that a crash in one site is
#                      attributable and cannot poison the next one.
#   --reindex=PATH     internal. Rebuild the ROOT search index and
#                      sitemap over the finished tree, in this process,
#                      and exit. See the comment on that branch.
#
# The packages must already be installed. pkgdown runs the examples and
# the articles against the INSTALLED package, so building against a
# stale install silently publishes stale output. The driver checks the
# installed version against each DESCRIPTION and refuses on a mismatch.

# ---------------------------------------------------------------- args

args <- commandArgs(trailingOnly = TRUE)

opt_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[length(hit)])
}
opt_flag <- function(name) paste0("--", name) %in% args

# The script lives in dev/release/, so the repository root is two up.
# Resolved from the script's own path rather than the working directory
# so that a run from any directory finds the packages.
script_path <- local({
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(hit)) sub("^--file=", "", hit[1]) else "dev/release/build-docs.R"
})
script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                      winslash = "/", mustWork = TRUE)

# ------------------------------------------------- one site, one child

# The child branch runs first and exits, so nothing below it can run in
# a worker by accident.
site_arg <- opt_value("site")
if (!is.null(site_arg)) {
  dest <- opt_value("dest")
  if (is.null(dest)) stop("--site needs --dest")
  pkgdown::build_site(
    pkg = site_arg,
    override = list(destination = dest),
    preview = FALSE,
    install = FALSE,
    new_process = FALSE,
    lazy = FALSE
  )
  quit(save = "no", status = 0)
}

# ------------------------------------ the root index, after everything

# The root search index and sitemap are built by WALKING THE
# DESTINATION, not from the package alone, so they cover the extension
# sites only once those sites are on disk. The live site has that
# coverage by accident: the old driver built core into a docs/ that still
# held the previous run's subsites. Measured on the committed tree, the
# root index carries 2,088 entries of which 125 paths are subsite pages,
# and the root sitemap 497 URLs of which 268 are.
#
# Core cannot simply be built last: pkgdown's check_dest_is_pkgdown()
# refuses a non-empty destination whose root holds no pkgdown.yml, and
# seven subdirectories are not a pkgdown site. So core is built first as
# before, and the root index is regenerated here at the end.
reindex_arg <- opt_value("reindex")
if (!is.null(reindex_arg)) {
  dest <- opt_value("dest")
  if (is.null(dest)) stop("--reindex needs --dest")
  p <- pkgdown::as_pkgdown(reindex_arg,
                           override = list(destination = dest))
  pkgdown::build_search(p)
  # build_sitemap() is not exported. Looked up by name, so that a rename
  # upstream is a loud failure here rather than a root sitemap that
  # quietly covers the core site only.
  bs <- get0("build_sitemap", envir = asNamespace("pkgdown"),
             mode = "function", ifnotfound = NULL)
  if (is.null(bs)) {
    stop("this pkgdown has no build_sitemap(); the root sitemap would ",
         "cover only the core site")
  }
  bs(p)
  quit(save = "no", status = 0)
}

# ------------------------------------------------------- the packages

# Discovered from the tree rather than listed, so that a new extension
# joins the site by existing. The count guard below is what keeps a
# MISSING one from passing silently.
pkg_dirs <- c(root, sort(Sys.glob(file.path(root, "extensions", "*"))))
pkg_dirs <- pkg_dirs[file.exists(file.path(pkg_dirs, "_pkgdown.yml"))]
pkg_name <- vapply(pkg_dirs, function(d) {
  as.character(read.dcf(file.path(d, "DESCRIPTION"))[1, "Package"])
}, character(1), USE.NAMES = FALSE)

# Core first, because it is the one that turns the destination into a
# pkgdown site and pkgdown refuses a non-empty destination that is not
# one. The root search index and sitemap are regenerated after the
# extensions by the --reindex step above; building core LAST was tried
# instead and fails in check_dest_is_pkgdown().
#
# The extensions then run in dependency order, which matters because
# pkgdown loads each package and frmtmb.learn imports frmtmb.eam.
#
# Core's build does not remove the subdirectories it did not write. That
# is not an assumption: docs/frmtmb.ddm/ survived every core rebuild
# since the package was renamed.
order_hint <- c("frmtmb", "frmtmb.eam", "frmtmb.latent", "frmtmb.ode",
                "frmtmb.spline", "frmtmb.coupling", "frmtmb.sample",
                "frmtmb.learn")
ord <- order(match(pkg_name, order_hint, nomatch = length(order_hint) + 1L),
             pkg_name)
pkg_dirs <- pkg_dirs[ord]
pkg_name <- pkg_name[ord]

only <- opt_value("pkgs")
if (!is.null(only)) {
  keep <- trimws(strsplit(only, ",", fixed = TRUE)[[1]])
  unknown <- setdiff(keep, pkg_name)
  if (length(unknown)) {
    stop("--pkgs names no such package: ", paste(unknown, collapse = ", "))
  }
  pkg_dirs <- pkg_dirs[pkg_name %in% keep]
  pkg_name <- pkg_name[pkg_name %in% keep]
}

# The subdirectory of each site comes from its own `_pkgdown.yml`, so
# the tree this produces is the one the committed site has and the one
# the core navbar links to. Core has no `destination` and takes the
# root.
sub_of <- function(dir) {
  cfg <- yaml::read_yaml(file.path(dir, "_pkgdown.yml"))
  if (is.null(cfg$destination)) return("")
  basename(cfg$destination)
}
pkg_sub <- vapply(pkg_dirs, sub_of, character(1), USE.NAMES = FALSE)

# A `destination` whose last segment is not the package name would put
# the site somewhere the navbar does not link. Caught here rather than
# by a reader noticing a 404.
bad <- which(nzchar(pkg_sub) & pkg_sub != pkg_name)
if (length(bad)) {
  stop("destination does not end in the package name: ",
       paste(pkg_name[bad], "->", pkg_sub[bad], collapse = "; "))
}

# ---------------------------------------------------- the destination

dest_root <- opt_value("dest", Sys.getenv("FRMTMB_DOCS_DEST", ""))
if (!nzchar(dest_root)) dest_root <- file.path(root, "docs")
dir.create(dest_root, recursive = TRUE, showWarnings = FALSE)
dest_root <- normalizePath(dest_root, winslash = "/", mustWork = TRUE)

min_sites <- as.integer(opt_value("min-sites", "8"))
log_path <- opt_value("log", file.path(root, "dev", "release", "docs.log"))
whole_tree <- is.null(only) || opt_flag("check-links")

# ------------------------------------------------------- preflight
#
# Everything that can refuse the run happens BEFORE the clean below, so
# a refusal leaves the site that is already in the destination intact.
# The directory itself has been created by now; nothing in it has been
# removed.

for (p in c("pkgdown", "yaml")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("package missing: ", p)
}

# pkgdown renders against the INSTALLED package. A version older than
# the checkout means the published reference and articles describe code
# that is not in the tree, which is the failure this whole file exists
# to remove from the release.
for (i in seq_along(pkg_dirs)) {
  want <- as.character(read.dcf(
    file.path(pkg_dirs[i], "DESCRIPTION"))[1, "Version"])
  have <- tryCatch(as.character(utils::packageVersion(pkg_name[i])),
                   error = function(e) NA_character_)
  if (is.na(have)) {
    stop("not installed: ", pkg_name[i],
         ". Install every package before building the site.")
  }
  if (!identical(have, want)) {
    stop("installed ", pkg_name[i], " is ", have, " but the checkout is ",
         want, ". pkgdown renders the installed package.")
  }
}

# The articles gate their fitted output on `requireNamespace()`, so a
# missing Suggests turns a page into its own table of contents with no
# numbers in it and the build still exits 0. The list is READ OUT of the
# article sources rather than written down here: a typed list would go
# stale the first time an article adds a gate, and this guard exists
# because the failure it catches is silent.
gate_pkgs <- local({
  rmd <- c(Sys.glob(file.path(root, "vignettes", "*.Rmd")),
           Sys.glob(file.path(root, "extensions", "*", "vignettes", "*.Rmd")))
  txt <- unlist(lapply(rmd, readLines, warn = FALSE))
  hit <- unlist(regmatches(
    txt, gregexpr('requireNamespace[(]"[^"]+"', txt)))
  sort(unique(sub('"$', "", sub('^requireNamespace[(]"', "", hit))))
})
gate_missing <- gate_pkgs[!vapply(gate_pkgs, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(gate_missing)) {
  msg <- paste0("articles gate themselves on ", length(gate_missing),
                " package(s) that are not installed: ",
                paste(gate_missing, collapse = ", "))
  if (opt_flag("require-articles")) stop(msg)
  warning(msg, call. = FALSE, immediate. = TRUE)
}

# ------------------------------------------------- clean, after that

if (!opt_flag("no-clean")) {
  # Refuse to empty a directory that is not already a pkgdown site.
  # A typo in --dest would otherwise delete an unrelated tree.
  present <- list.files(dest_root, all.files = TRUE, no.. = TRUE)
  # A tree of subsites and no core site is a legitimate destination: a
  # partial build of the extensions alone produces exactly that, and the
  # first spelling of this check refused it. Found by the guard case for
  # the count arm, which failed here instead and so did not test what it
  # was for.
  subsite <- length(present) > 0 && all(vapply(present, function(e) {
    file.exists(file.path(dest_root, e, "pkgdown.yml"))
  }, logical(1)))
  looks_like_site <- any(c("pkgdown.yml", "index.html") %in% present) ||
    subsite
  if (length(present) && !looks_like_site) {
    stop("refusing to clean ", dest_root,
         ": not empty and holds no pkgdown.yml or index.html. ",
         "Pass --no-clean, or point --dest somewhere else.")
  }
  unlink(file.path(dest_root, present), recursive = TRUE, force = TRUE)
  left <- list.files(dest_root, all.files = TRUE, no.. = TRUE)
  if (length(left)) {
    stop("clean left ", length(left), " entries in ", dest_root)
  }
}

# Children get the parent's library paths. R_LIBS is the only channel
# that survives a fresh Rscript, because a profile-set .libPaths() does
# not.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
rscript <- file.path(R.home("bin"), "Rscript")

unlink(log_path)
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = log_path, append = TRUE)
}

say("DOCS root ", root)
say("DOCS dest ", dest_root)
say("DOCS sites ", length(pkg_dirs), ": ", paste(pkg_name, collapse = " "))
say("DOCS article gates ", length(gate_pkgs), ", missing ",
    length(gate_missing),
    if (length(gate_missing)) paste0(": ",
      paste(gate_missing, collapse = " ")) else "")
if (opt_flag("no-clean")) {
  # With the destination kept, a file left by an EARLIER run satisfies
  # the existence checks below. Said out loud rather than left for a
  # reader to work out from a green line.
  say("DOCS WARNING --no-clean: a file from an earlier run can satisfy ",
      "a check this run did not produce.")
}

# ----------------------------------------------------------- build

exit_of <- integer(length(pkg_dirs))
secs_of <- numeric(length(pkg_dirs))
# Wall clock, not proc.time(): the work happens in a child process, and
# the number is only ever used to size a CI timeout, where minutes are
# the unit. proc.time() also ticks at 10 ms on the development box.
t0 <- Sys.time()
for (i in seq_along(pkg_dirs)) {
  target <- if (nzchar(pkg_sub[i])) {
    file.path(dest_root, pkg_sub[i])
  } else {
    dest_root
  }
  say("===== SITE ", pkg_name[i], " -> ", target, " =====")
  ti <- Sys.time()
  out <- suppressWarnings(system2(
    rscript,
    args = c(shQuote(script_path),
             paste0("--site=", shQuote(pkg_dirs[i])),
             paste0("--dest=", shQuote(target))),
    stdout = TRUE, stderr = TRUE
  ))
  secs_of[i] <- as.numeric(difftime(Sys.time(), ti, units = "secs"))
  status <- attr(out, "status")
  exit_of[i] <- if (is.null(status)) 0L else as.integer(status)
  cat(out, sep = "\n", file = log_path, append = TRUE)
  cat("\n", file = log_path, append = TRUE)
  say("===== ", if (exit_of[i] == 0L) "ok" else
      paste0("FAILED exit ", exit_of[i]),
      " ", sprintf("%.0f", secs_of[i]), "s =====")
}

# The root index has to see the finished tree, so it is rebuilt here and
# only on a full build: a root index written over a partial tree would
# advertise a site that is not all there. Its own child process, like
# every site.
reindex_exit <- NA_integer_
core_i <- match("frmtmb", pkg_name)
if (whole_tree && !is.na(core_i) && exit_of[core_i] == 0L) {
  say("===== REINDEX root search and sitemap =====")
  out <- suppressWarnings(system2(
    rscript,
    args = c(shQuote(script_path),
             paste0("--reindex=", shQuote(pkg_dirs[core_i])),
             paste0("--dest=", shQuote(dest_root))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  reindex_exit <- if (is.null(status)) 0L else as.integer(status)
  cat(out, sep = "\n", file = log_path, append = TRUE)
  cat("\n", file = log_path, append = TRUE)
  say("===== ", if (reindex_exit == 0L) "ok" else
      paste0("FAILED exit ", reindex_exit), " =====")
}

total_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ----------------------------------------------------------- guards
#
# Every assertion below is on a FILE the run produced, not on the loop
# that produced it. A driver reporting success over zero work has
# happened three times in this project, and in each case the count came
# from launches rather than from results.

fail <- character(0)
note <- function(...) fail <<- c(fail, paste0(...))

for (i in seq_along(pkg_dirs)) {
  if (exit_of[i] != 0L) {
    note("site ", pkg_name[i], " exited ", exit_of[i])
  }
}

built <- 0L
for (i in seq_along(pkg_dirs)) {
  target <- if (nzchar(pkg_sub[i])) {
    file.path(dest_root, pkg_sub[i])
  } else {
    dest_root
  }
  index <- file.path(target, "index.html")
  refs <- file.path(target, "reference", "index.html")
  ok <- file.exists(index) && file.size(index) > 0
  if (!ok) {
    note("no usable index.html for ", pkg_name[i], " at ", index)
  } else if (!file.exists(refs)) {
    # A home page alone is what a site that died partway leaves behind.
    note("no reference/index.html for ", pkg_name[i])
  } else {
    built <- built + 1L
  }
}

if (built < min_sites) {
  note("built ", built, " sites, expected at least ", min_sites)
}

# The two checks below are about the WHOLE tree, so a deliberately
# partial build would fail them for the wrong reason and make --pkgs
# useless for a quick local rebuild. They run on a full build, and on
# demand with --check-links against a tree that is already there.
if (whole_tree) {
  if (!is.na(reindex_exit) && reindex_exit != 0L) {
    note("the root reindex exited ", reindex_exit)
  }
  search_index <- file.path(dest_root, "search.json")
  if (!file.exists(search_index)) {
    note("no search.json at the root of ", dest_root)
  } else {
    # THE ROOT SEARCH BOX MUST REACH THE EXTENSIONS. It does so only
    # because the index is rebuilt after the subsites exist, and the
    # first version of this script lost that silently: the site still
    # built, every link still worked, and the root search went from
    # 2,088 entries to 1,110. Existence of the file was not enough, so
    # the check is on its CONTENT, one extension at a time.
    idx <- readLines(search_index, warn = FALSE)
    for (s in pkg_sub[nzchar(pkg_sub)]) {
      if (!any(grepl(paste0("/", s, "/"), idx, fixed = TRUE))) {
        note("the root search index has no path under ", s)
      }
    }
    sm <- file.path(dest_root, "sitemap.xml")
    if (!file.exists(sm)) {
      note("no sitemap.xml at the root of ", dest_root)
    } else {
      smx <- readLines(sm, warn = FALSE)
      for (s in pkg_sub[nzchar(pkg_sub)]) {
        if (!any(grepl(paste0("/", s, "/"), smx, fixed = TRUE))) {
          note("the root sitemap has no URL under ", s)
        }
      }
    }
  }

  # The core navbar is what a reader uses to reach an extension site.
  # Its hrefs are absolute URLs into the deployed site, so a renamed or
  # unbuilt extension shows up only as a 404 unless it is checked here.
  core_cfg <- yaml::read_yaml(file.path(root, "_pkgdown.yml"))
  menu <- core_cfg$navbar$components$extensions$menu
  hrefs <- vapply(menu, function(x) {
    if (is.null(x$href)) NA_character_ else as.character(x$href)
  }, character(1))
  hrefs <- hrefs[!is.na(hrefs)]
  if (!length(hrefs)) note("the core navbar has no extensions menu")
  for (h in hrefs) {
    seg <- basename(sub("/+$", "", h))
    if (!file.exists(file.path(dest_root, seg, "index.html"))) {
      note("navbar href ", h, " has no built site at ", seg,
           "/index.html")
    }
  }
} else {
  say("DOCS whole-tree checks SKIPPED: this is a partial build. ",
      "Do not deploy it. Pass --check-links to run them anyway.")
}

# Anything in the destination that is not the core site or one of the
# built subdirectories. This is what left docs/frmtmb.ddm/ behind after
# the package was renamed to frmtmb.eam.
subdirs <- list.dirs(dest_root, recursive = FALSE, full.names = FALSE)
stray <- setdiff(subdirs, c(pkg_sub[nzchar(pkg_sub)],
                            "articles", "deps", "news", "reference"))
if (length(stray)) {
  say("DOCS stray directories: ", paste(stray, collapse = " "))
}

say("DOCS seconds per site: ",
    paste(pkg_name, sprintf("%.0f", secs_of), collapse = ", "))
say("DOCS total ", sprintf("%.0f", total_secs), "s (",
    sprintf("%.1f", total_secs / 60), " min)")
say("DOCS built ", built, " of ", length(pkg_dirs),
    " (minimum ", min_sites, ")")

if (length(fail)) {
  say("DOCS FAILED")
  for (f in fail) say("  - ", f)
  stop("the documentation build failed ", length(fail),
       " check(s); see ", log_path)
}
say("DOCS ok")
