# Compare local R libraries with CRAN for the round's dependency closure,
# and print what lags. It reads only: it installs and removes nothing.
#
# WHY. At the 0.79.0 release the local library held RTMBdist 1.1.0, while
# CI installed 1.2.0, published inside the round, and went red on it. The
# local suites had passed on a package version that CI no longer reads.
# Run this before the final suites of a round, and update the library
# (or say why not) before they run.
#
#   Rscript dev/release/cran-lag.R [libraries] [CRAN mirror]
#
# libraries    a ";"-separated stack, read in order as the suites read
#              it, for example "<release library>;<user library>".
#              Default: .libPaths() of the R that runs this, which for a
#              plain Rscript on this machine is the user library,
#              C:/Users/adf44/AppData/Local/R/win-library/4.6, and R's
#              own library.
# CRAN mirror  default: https://cloud.r-project.org
#
# The closure: Depends, Imports, LinkingTo and Suggests of core and of
# every extension under extensions/, plus the release tools covr,
# rcmdcheck, roxygen2, testthat, srr and pkgstats; then the recursive
# Depends, Imports and LinkingTo of all of those. Base packages are left
# out. A package that is not on CRAN is looked up in the
# Additional_repositories its DESCRIPTION names (RTMBode, from
# kaskr.r-universe.dev).
#
# A package resolves to its first copy in the stack, then in R's own
# library (.Library), which holds the recommended packages. The
# reference is the SOURCE version on the repository, which is what the
# Ubuntu runners install; on Windows the binary version is printed too,
# because a binary can trail its source by a few days.
#
# It fails closed. available.packages() only warns when it cannot read an
# index, and returns an empty table; that once made every package read
# "not in any repository" and the run end "nothing lags" with CRAN
# unreachable. So an index that cannot be read stops the run, an empty
# CRAN index stops it, and a package in no repository counts as a
# failure (a Suggests from GitHub only, say, is one that the runners
# cannot install either).
#
# Exit status: 0 when every package is current; 1 when a package lags,
# is missing or is in no repository, or when an index cannot be read.
args <- commandArgs(trailingOnly = TRUE)
libs <- if (length(args) >= 1L && nzchar(args[[1L]])) {
  strsplit(args[[1L]], ";", fixed = TRUE)[[1L]]
} else {
  .libPaths()
}
cran <- if (length(args) >= 2L) args[[2L]] else "https://cloud.r-project.org"
if (!all(dir.exists(libs))) {
  stop("no such library: ", paste(libs[!dir.exists(libs)], collapse = ", "))
}
libs <- unique(normalizePath(c(libs, .Library), winslash = "/"))

# an index that cannot be read is an error, not an empty table
read_index <- function(repos, type) {
  withCallingHandlers(
    available.packages(repos = repos, type = type),
    warning = function(w) stop(conditionMessage(w), call. = FALSE))
}

# the repository root: this file is dev/release/cran-lag.R. When it is
# sourced from another script, --file names that script, so fall back
# to the working directory
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                 value = TRUE))
root <- if (length(self) == 1L) {
  normalizePath(file.path(dirname(self), "..", ".."), winslash = "/",
                mustWork = FALSE)
} else {
  ""
}
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  root <- normalizePath(".", winslash = "/")
}
descs <- c(file.path(root, "DESCRIPTION"),
           Sys.glob(file.path(root, "extensions", "*", "DESCRIPTION")))
stopifnot(file.exists(descs[[1L]]), length(descs) >= 2L)

dep_names <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  p <- trimws(sub("[(].*", "", strsplit(x, ",", fixed = TRUE)[[1L]]))
  p[nzchar(p)]
}
fields <- c("Depends", "Imports", "LinkingTo", "Suggests")
seed <- character()
extra_repos <- character()
local_pkgs <- character()
for (d in descs) {
  dcf <- read.dcf(d)
  local_pkgs <- c(local_pkgs, dcf[1L, "Package"])
  for (f in intersect(fields, colnames(dcf))) {
    seed <- c(seed, dep_names(dcf[1L, f]))
  }
  if ("Additional_repositories" %in% colnames(dcf)) {
    extra_repos <- c(extra_repos,
                     dep_names(dcf[1L, "Additional_repositories"]))
  }
}
tools_pkgs <- c("covr", "rcmdcheck", "roxygen2", "testthat", "srr",
                "pkgstats")
base_pkgs <- c("R", rownames(installed.packages(priority = "base")))
seed <- setdiff(unique(c(seed, tools_pkgs)), c(base_pkgs, local_pkgs))

ap_cran <- read_index(cran, "source")
if (nrow(ap_cran) == 0L) stop("the CRAN index at ", cran, " is empty")
ap_extra <- if (length(extra_repos)) {
  read_index(unique(extra_repos), "source")
} else {
  ap_cran[0L, , drop = FALSE]
}
# CRAN first: an extra repository's copy counts only for a package CRAN
# does not carry
ap <- rbind(ap_cran, ap_extra[!rownames(ap_extra) %in% rownames(ap_cran), ,
                              drop = FALSE])
closure <- tools::package_dependencies(
  seed, db = ap, which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE)
pkgs <- sort(setdiff(unique(c(seed, unlist(closure, use.names = FALSE))),
                     c(base_pkgs, local_pkgs)))

# the binary column is information only, so its index may fail
bin <- if (.Platform$OS.type == "windows") {
  tryCatch(read_index(cran, "win.binary"), error = function(e) {
    cat("binary index not read:", conditionMessage(e), "\n")
    NULL
  })
}
ip <- installed.packages(lib.loc = libs)
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
rownames(ip) <- ip[, "Package"]

ver_of <- function(tab, p, col = "Version") {
  if (!is.null(tab) && p %in% rownames(tab)) tab[p, col] else NA_character_
}
rows <- lapply(pkgs, function(p) {
  have <- ver_of(ip, p)
  want <- ver_of(ap, p)
  data.frame(
    package = p,
    installed = have,
    where = if (is.na(have)) NA_character_ else
      normalizePath(ip[p, "LibPath"], winslash = "/"),
    repository = want,
    binary = ver_of(bin, p),
    from = if (is.na(want)) NA_character_ else
      if (p %in% rownames(ap_cran)) "CRAN" else
        sub("/src/contrib.*$", "", ap[p, "Repository"]),
    stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab$status <- ifelse(
  is.na(tab$repository), "not in any repository",
  ifelse(is.na(tab$installed), "MISSING",
         ifelse(package_version(tab$installed, strict = FALSE) <
                  package_version(tab$repository, strict = FALSE),
                "LAGS", "current")))
# a version string package_version() cannot read must not pass as current
tab$status[is.na(tab$status)] <- "unreadable version"

cat("libraries:", paste0("\n  ", libs), "\n")
cat("repositories:", cran, if (length(extra_repos)) extra_repos, "\n")
cat("closure:", nrow(tab), "packages from", length(seed), "direct ones",
    "of", length(descs), "DESCRIPTION files\n")
for (s in c("current", "LAGS", "MISSING", "not in any repository",
            "unreadable version")) {
  cat(sprintf("  %-22s %d\n", s, sum(tab$status == s)))
}
show <- tab[tab$status != "current", , drop = FALSE]
if (nrow(show)) {
  cat("\n")
  print(show[, c("package", "status", "installed", "where", "repository",
                 "binary", "from")], row.names = FALSE)
}
bad <- sum(tab$status != "current")
cat(if (bad) {
  sprintf("\n%d package(s) lag, are missing or are in no repository\n",
          bad)
} else {
  "\nnothing lags\n"
})
quit(status = if (bad) 1L else 0L)
