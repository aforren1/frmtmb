# The functions each lane removed or re-signed, against the base
# 9e902909, and every call site of them left in the merged release
# tree. A clean textual merge can still call a function another lane
# removed (dev/round-handoff.md, 0.66.0), so each one is read.
#
#   Rscript dev/rel068-removed.R
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
base <- "9e902909"
defs <- function(lines) {
  ex <- tryCatch(parse(text = lines, keep.source = FALSE),
                 error = function(e) NULL)
  out <- list()
  for (e in ex) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
        is.name(e[[2]]) && is.call(e[[3]]) &&
        identical(as.character(e[[3]][[1]]), "function")) {
      out[[as.character(e[[2]])]] <- paste(names(e[[3]][[2]]),
                                           collapse = ",")
    }
  }
  out
}
all_r <- function(dir) {
  c(list.files(file.path(dir, "R"), "[.]R$", full.names = TRUE),
    list.files(file.path(dir, "extensions"), "[.]R$", recursive = TRUE,
               full.names = TRUE))
}
call_sites <- function(fn) {
  hits <- character(0)
  for (f in all_r(rel)) {
    if (grepl("/dev/|/man/", f)) next
    l <- readLines(f, warn = FALSE)
    w <- grep(paste0("(^|[^A-Za-z0-9._])", gsub(".", "[.]", fn, fixed = TRUE),
                     "[(]"), l)
    w <- w[!grepl(paste0("^\\s*", fn, " <- function"), l[w])]
    if (length(w)) hits <- c(hits, paste0(sub(paste0(rel, "/"), "", f), ":",
                                          w))
  }
  hits
}
for (lane in c("fixes", "gpby", "ordmix")) {
  src <- file.path(root, paste0("frmtmb-wt-", lane))
  ch <- system2("git", c("-C", src, "diff", "--name-only", "HEAD"),
                stdout = TRUE)
  ch <- ch[grepl("[.]R$", ch) & grepl("(^|/)R/", ch)]
  cat("== lane", lane, ":", length(ch), "changed R files\n")
  for (f in ch) {
    b <- defs(system2("git", c("-C", src, "show", paste0(base, ":", f)),
                      stdout = TRUE))
    n <- defs(readLines(file.path(src, f), warn = FALSE))
    gone <- setdiff(names(b), names(n))
    resig <- Filter(function(x) !identical(b[[x]], n[[x]]),
                    intersect(names(b), names(n)))
    for (x in gone) {
      cs <- call_sites(x)
      cat(sprintf("  REMOVED %s (%s): %d call sites in the release tree%s\n",
                  x, f, length(cs),
                  if (length(cs)) paste0(": ", paste(cs, collapse = " ")) else
                    ""))
    }
    for (x in resig) {
      cat(sprintf("  RESIGNED %s (%s): (%s) -> (%s); %d call sites\n", x, f,
                  b[[x]], n[[x]], length(call_sites(x))))
    }
  }
}
cat("REMOVED DONE\n")
