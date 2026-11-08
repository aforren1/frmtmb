# Reviewer check of dev/release/cran-lag.R, item 4.
#  (1) Run the lane's script with quit() trapped, keep its table.
#  (2) An independent closure: tools:::.split_dependencies() on each
#      DESCRIPTION, package_dependencies() level by level by hand.
#  (3) Version comparison on the forms that matter.
#  (4) The library's mtime and package count before and after, so
#      "installs nothing" is measured.
#
#   cd <worktree>; Rscript dev/review-ci-cranlag.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[[1L]] else
  "C:/Users/adf44/AppData/Local/R/win-library/4.6"
before <- list(n = length(list.dirs(lib, recursive = FALSE)),
               mt = max(file.mtime(list.dirs(lib, recursive = FALSE))))

# (1) the lane's script, its quit() trapped
lane <- new.env()
lane$quit <- function(status = 0L, ...) {
  assign("lane_status", status, envir = globalenv())
  invokeRestart("lane_done")
}
withRestarts(sys.source("dev/release/cran-lag.R", envir = lane),
             lane_done = function() NULL)
cat("lane exit status:", lane_status, "\n")
lane_set <- lane$tab$package

# (2) independent closure
cran <- "https://cloud.r-project.org"
descs <- c("DESCRIPTION", Sys.glob("extensions/*/DESCRIPTION"))
own <- vapply(descs, function(d) read.dcf(d, "Package")[1, 1], "")
first <- character()
for (d in descs) {
  x <- read.dcf(d)
  for (fl in intersect(c("Depends", "Imports", "LinkingTo", "Suggests"),
                       colnames(x))) {
    sp <- tools:::.split_dependencies(x[1, fl])
    first <- c(first, names(sp))
  }
}
first <- unique(c(first, "covr", "rcmdcheck", "roxygen2", "testthat", "srr",
                  "pkgstats"))
basep <- c("R", rownames(installed.packages(.Library, priority = "base")))
first <- setdiff(first, c(basep, own))
db <- available.packages(repos = c(CRAN = cran,
                                   kaskr = "https://kaskr.r-universe.dev"),
                         type = "source")
seen <- character(); frontier <- first
while (length(frontier)) {
  seen <- union(seen, frontier)
  nx <- character()
  for (p in intersect(frontier, rownames(db))) {
    for (fl in c("Depends", "Imports", "LinkingTo")) {
      v <- db[p, fl]
      if (!is.na(v)) nx <- c(nx, names(tools:::.split_dependencies(v)))
    }
  }
  frontier <- setdiff(setdiff(unique(nx), c(basep, own)), seen)
}
mine <- sort(seen)
cat("independent closure:", length(mine), "packages from", length(first),
    "direct; lane:", length(lane_set), "\n")
cat("  in mine, not the lane's:", paste(setdiff(mine, lane_set),
                                         collapse = " "), "\n")
cat("  in the lane's, not mine:", paste(setdiff(lane_set, mine),
                                         collapse = " "), "\n")
nr <- setdiff(mine, rownames(db))
cat("  in neither repository:", paste(nr, collapse = " "), "\n")
# recommended packages in the closure, and where they resolve
rec <- rownames(installed.packages(.Library, priority = "recommended"))
cat("  recommended packages in the closure:",
    paste(intersect(mine, rec), collapse = " "), "\n")

# (3) version comparisons
cmp <- function(a, b) package_version(a, strict = FALSE) <
  package_version(b, strict = FALSE)
cases <- rbind(c("3.1-171", "3.1-172"), c("3.1-172", "3.1-171"),
               c("1.7-6", "1.7-10"), c("2.39.1.1", "2.39.1"),
               c("0.9", "0.10"), c("1.0", "1.0-0"), c("1.0.0", "1.0"),
               c("2.3.3", "2.3.4"), c("abc", "1.0"))
for (i in seq_len(nrow(cases))) {
  cat(sprintf("  %s < %s : %s\n", cases[i, 1], cases[i, 2],
              format(cmp(cases[i, 1], cases[i, 2]))))
}

after <- list(n = length(list.dirs(lib, recursive = FALSE)),
              mt = max(file.mtime(list.dirs(lib, recursive = FALSE))))
cat("library dirs before/after:", before$n, after$n,
    " newest mtime unchanged:", identical(before$mt, after$mt), "\n")
