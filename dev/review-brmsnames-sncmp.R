## Compare the reviewer's re-recorded suite names with the lane's.
for (a in c("base", "lane")) for (f in c("test-v15", "test-emmeans", "test-edgecases")) {
  mine <- readRDS(sprintf("dev/review-brmsnames-suite/%s/frmtmb/%s.rds", a, f))
  lane <- readRDS(sprintf("dev/brmsnames-1024-suite/%s/frmtmb/%s.rds", a, f))
  strip <- function(x) lapply(x$fits, function(r) r[setdiff(names(r), "call")])
  cat(a, f, "fits", length(mine$fits), length(lane$fits), "identical names:",
      identical(strip(mine), strip(lane)), "\n")
}
