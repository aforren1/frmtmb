# Reviewer check (lane ceplot): does the installed lane library hold the
# worktree source? Every function literal assigned at top level in R/ of
# core and of frmtmb.sample is compared, by deparsed formals and body,
# with the installed namespace binding. The control is the same check
# against rellib-r4, which must list exactly the functions the diff
# changed or added.
#   Rscript dev/ceplot-rev-libcheck.R lane|base
args <- commandArgs(TRUE)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot"
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
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
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
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
  cat(sprintf("%s: %d source functions; %d not installed; %d differ\n",
              pkg, length(nms), length(missing), length(differ)))
  if (length(missing)) cat("  missing:", paste(missing, collapse = " "), "\n")
  if (length(differ)) cat("  differ:", paste(differ, collapse = " "), "\n")
}
check("frmtmb", file.path(wt, "R"))
check("frmtmb.sample", file.path(wt, "extensions/frmtmb.sample/R"))
