# Reviewer check: does the installed lane library hold the worktree
# source? Every function defined in R/ of core and of frmtmb.sample is
# sourced into a scratch environment and compared, by deparsed formals
# and body, with the binding in the installed namespace. The control is
# the same comparison against the base build rellib-r3, which must
# report exactly the functions the diff changed.
#   Rscript dev/postfit2-rev-libcheck.R > dev/postfit2-rev-log/libcheck.txt
args <- commandArgs(TRUE)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm", arm, ": frmtmb from", find.package("frmtmb"),
    "; frmtmb.sample from", find.package("frmtmb.sample"), "\n")
same_fn <- function(a, b) {
  identical(deparse(formals(a)), deparse(formals(b))) &&
    identical(deparse(body(a)), deparse(body(b)))
}
check <- function(pkg, rdir) {
  ns <- asNamespace(pkg)
  src <- new.env(parent = ns)
  files <- list.files(rdir, "[.]R$", full.names = TRUE)
  bad_files <- character(0)
  for (f in files) {
    ex <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
    if (is.null(ex)) {
      bad_files <- c(bad_files, basename(f))
      next
    }
    for (e in ex) {
      # only assignments of a function literal: a top-level call could
      # have side effects, and data objects are not the question here
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
          is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")) &&
          is.name(e[[2]])) {
        assign(as.character(e[[2]]), eval(e[[3]], src), envir = src)
      }
    }
  }
  nms <- ls(src, all.names = TRUE)
  missing <- nms[!vapply(nms, exists, NA, envir = ns, inherits = FALSE)]
  differ <- character(0)
  for (n in setdiff(nms, missing)) {
    inst <- get(n, envir = ns, inherits = FALSE)
    if (!is.function(inst)) next
    if (!same_fn(get(n, envir = src), inst)) differ <- c(differ, n)
  }
  cat(sprintf("%s: %d source functions; %d not in the installed namespace; %d differ\n",
              pkg, length(nms), length(missing), length(differ)))
  if (length(bad_files)) cat("  unparsed files:", bad_files, "\n")
  if (length(missing)) cat("  missing:", paste(missing, collapse = " "), "\n")
  if (length(differ)) cat("  differ:", paste(differ, collapse = " "), "\n")
}
check("frmtmb", file.path(wt, "R"))
check("frmtmb.sample", file.path(wt, "extensions/frmtmb.sample/R"))
