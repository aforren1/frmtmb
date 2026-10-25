# Review of the 0.75.0 merge: line accounting per file.
# For every file any lane changed (lane tree vs 18b536fc), and every file
# the release tree changed (vs 4fd134db):
#  (a) lane-added lines missing from the release file (lost hunks);
#  (b) lane-removed lines still present in the release file, where no
#      other lane keeps them and HEAD had them (a lost deletion);
#  (c) release-added lines that no lane added (resolutions/integration);
#  (d) release-removed lines that no lane removed.
# Counts are multisets of exact lines (CR stripped).
#   Rscript dev/review-rel075-merge.R > dev/review-rel075-log/merge.txt
rel <- "C:/Users/adf44/source/r/frmtmb-wt-rel075"
lanes <- c("smallfix", "bhaz", "brmsnames", "sarfcor")
lt <- function(l) sprintf("C:/Users/adf44/source/r/frmtmb-wt-%s", l)
git <- function(dir, ...) suppressWarnings(system2("git", c("-C", dir, ...),
                                                  stdout = TRUE, stderr = FALSE))
blob <- function(rev, f) {
  x <- git(rel, "show", paste0(rev, ":", f))
  if (!is.null(attr(x, "status"))) return(NULL)
  sub("\r$", "", x)
}
rdf <- function(p) if (file.exists(p)) sub("\r$", "", readLines(p, warn = FALSE)) else NULL
gen <- "^(man/.*\\.Rd|NAMESPACE|extensions/[^/]*/man/.*\\.Rd|extensions/[^/]*/NAMESPACE)$"
lane_files <- list()
for (l in lanes) {
  st <- git(lt(l), "diff", "--name-status", "18b536fc")
  st <- do.call(rbind, strsplit(st, "\t"))
  lane_files[[l]] <- st[, 2]
  un <- git(lt(l), "ls-files", "--others", "--exclude-standard")
  un <- character(0)
  lane_files[[l]] <- c(lane_files[[l]], un)
}
# multiset difference a - b
msdiff <- function(a, b) {
  if (is.null(a)) return(character(0))
  if (is.null(b)) return(a)
  ua <- unique(a)
  ca <- tabulate(match(a, ua), length(ua))
  cb <- tabulate(match(b, ua), length(ua))
  d <- ca - cb
  rep(ua[d > 0], d[d > 0])
}
allf <- unique(unlist(lane_files))
relst <- git(rel, "status", "--porcelain", "--untracked-files=all")
relf <- substring(relst, 4)
relf <- character(0)
allf <- sort(unique(c(allf, relf)))
allf <- allf[grepl("^dev/", allf)]
cat("files considered (outside dev/):", length(allf), "\n")
for (f in allf) {
  base <- blob("18b536fc", f); head <- blob("4fd134db", f)
  relx <- rdf(file.path(rel, f))
  lx <- lapply(lanes, function(l) {
    if (f %in% lane_files[[l]]) rdf(file.path(lt(l), f)) else base
  })
  names(lx) <- lanes
  changed_by <- lanes[vapply(lanes, function(l) f %in% lane_files[[l]], NA)]
  msgs <- character(0)
  for (l in changed_by) {
    add <- msdiff(lx[[l]], base)
    miss <- msdiff(add, relx)
    if (length(miss)) msgs <- c(msgs, sprintf("  (a) %s added, missing in release: %d", l, length(miss)),
                                paste0("      | ", head(miss, 15)))
    rem <- msdiff(base, lx[[l]])
    keep <- msdiff(rem, msdiff(base, relx))
    # still present lines that the lane removed
    if (length(keep)) {
      others <- setdiff(changed_by, l)
      kept_by_other <- vapply(keep, function(z) any(vapply(others, function(o) z %in% msdiff(lx[[o]], base), NA)), NA)
      keep <- keep[!kept_by_other & nzchar(trimws(keep))]
      if (length(keep)) msgs <- c(msgs, sprintf("  (b) %s removed, still in release: %d", l, length(keep)),
                                  paste0("      | ", head(keep, 15)))
    }
  }
  radd <- msdiff(relx, head)
  laneadd <- unlist(lapply(lanes, function(l) msdiff(lx[[l]], base)))
  headadd <- msdiff(head, base)
  un <- msdiff(radd, c(laneadd, headadd))
  un <- un[nzchar(trimws(un))]
  if (length(un) && !grepl(gen, f)) msgs <- c(msgs, sprintf("  (c) release added, no lane added: %d", length(un)),
                            paste0("      | ", head(un, 40)))
  rrem <- msdiff(head, relx)
  lanerem <- unlist(lapply(lanes, function(l) msdiff(base, lx[[l]])))
  un2 <- msdiff(rrem, lanerem)
  un2 <- un2[nzchar(trimws(un2))]
  if (length(un2) && !grepl(gen, f)) msgs <- c(msgs, sprintf("  (d) release removed, no lane removed: %d", length(un2)),
                             paste0("      | ", head(un2, 40)))
  if (length(msgs)) {
    cat("==", f, "[", paste(changed_by, collapse = ","), "]\n")
    cat(msgs, sep = "\n")
  }
}
cat("DONE\n")
