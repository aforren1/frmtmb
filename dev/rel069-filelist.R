# The release tree's diff against 4f5ea39f, file by file, with the
# lanes that changed each file and whether the consolidation edited it
# after the merge, so the coordinator can commit per lane.
#
#   Rscript dev/rel069-filelist.R > dev/rel069-log/filelist.txt
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
git <- function(dir, ...) {
  suppressWarnings(system2("git", c("-C", dir, ...), stdout = TRUE,
                           stderr = FALSE))
}
lane_files <- function(l) {
  d <- file.path(root, paste0("frmtmb-wt-", l))
  c(git(d, "diff", "--name-only", "HEAD"),
    git(d, "ls-files", "--others", "--exclude-standard"))
}
lanes <- c("ciharden", "surface", "setier", "optima")
lf <- setNames(lapply(lanes, lane_files), lanes)
tracked <- git(rel, "diff", "--name-only", "4f5ea39f")
deleted <- git(rel, "diff", "--name-only", "--diff-filter=D", "4f5ea39f")
untracked <- git(rel, "ls-files", "--others", "--exclude-standard")
all <- sort(unique(c(tracked, untracked)))
# a file the consolidation changed after the merge: its content differs
# from every lane's copy of it (for a file one lane alone touched), or
# it is in no lane at all
same_as_lane <- function(f, l) {
  a <- file.path(rel, f)
  b <- file.path(root, paste0("frmtmb-wt-", l), f)
  # bytes without carriage returns, so a binary file compares too
  rb <- function(p) {
    r <- readBin(p, "raw", file.size(p))
    r[r != as.raw(13)]
  }
  file.exists(a) && file.exists(b) && identical(rb(a), rb(b))
}
out <- data.frame(file = all, status = ifelse(all %in% deleted, "D",
                                       ifelse(all %in% tracked, "M", "A")),
                  lanes = "", edited = "", stringsAsFactors = FALSE)
for (i in seq_along(all)) {
  f <- all[i]
  ls <- lanes[vapply(lanes, function(l) f %in% lf[[l]], NA)]
  out$lanes[i] <- if (length(ls)) paste(ls, collapse = "+") else "release"
  out$edited[i] <- if (length(ls) == 1L) {
    if (same_as_lane(f, ls)) "" else "edited"
  } else if (length(ls) > 1L) "merged" else "new"
}
for (g in unique(out$lanes)) {
  s <- out[out$lanes == g, ]
  cat("\n## ", g, " (", nrow(s), " files)\n", sep = "")
  for (j in seq_len(nrow(s))) {
    cat(sprintf("%s %s%s\n", s$status[j], s$file[j],
                if (nzchar(s$edited[j])) paste0("  [", s$edited[j], "]") else ""))
  }
}
cat("\nTOTAL", nrow(out), "files:", sum(out$status == "M"), "modified,",
    sum(out$status == "A"), "added,",
    sum(out$status == "D"), "deleted\n")
