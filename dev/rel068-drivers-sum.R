# The merged-tree reruns of dev/rel068-ordmix-drivers.sh against each
# lane's own final run of the same driver (read from the lane worktree,
# not edited): per log, the REP lines that differ, and the warning
# counts the reviews quoted.
#
#   Rscript dev/rel068-drivers-sum.R > dev/rel068-log/drivers-sum.txt
rel <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel068-log/drivers"
om <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev"
fx <- "C:/Users/adf44/source/r/frmtmb-wt-fixes/dev/fixes-log"
reps <- function(f) {
  if (!file.exists(f)) return(NA_character_)
  x <- readLines(f, warn = FALSE)
  x[grepl("^REP ", x)]
}
cmp <- function(tag, a, b) {
  ra <- reps(a)
  rb <- reps(b)
  if (identical(ra, NA_character_)) {
    cat(sprintf("%-28s lane log missing\n", tag))
    return(invisible())
  }
  same <- sum(ra %in% rb)
  cat(sprintf("%-28s lane %3d  release %3d  identical lines %3d\n", tag,
              length(ra), length(rb), same))
  for (d in setdiff(rb, ra)) cat("    release:", d, "\n")
  for (d in setdiff(ra, rb)) cat("    lane:   ", d, "\n")
}
cat("== round-1 false-alarm table (dev/ordmix-p1-fa.sh)\n")
for (f in list.files(file.path(rel, "fa"))) {
  cmp(paste0("fa/", f), file.path(om, "ordmix-p1-fa-log", f),
      file.path(rel, "fa", f))
}
cat("\n== margin designs (dev/ordmix-rev2-margin.R)\n")
sound <- function(x) {
  se <- grepl("se_finite=TRUE", x)
  re <- as.numeric(sub(".*slope_relerr=([0-9.]+).*", "\\1", x))
  se & re < 0.5
}
for (f in list.files(file.path(rel, "margin"))) {
  rb <- reps(file.path(rel, "margin", f))
  s <- sound(rb)
  w <- grepl("degwarn=TRUE", rb)
  cat(sprintf("%-20s fits %2d sound %2d sound-and-warn %d warn %d\n", f,
              length(rb), sum(s), sum(s & w), sum(w)))
  cmp(paste0("margin/", f), file.path(om, "ordmix-p2-margin-log", f),
      file.path(rel, "margin", f))
}
cat("\n== B1 edge cases (dev/ordmix-rev2-b1.R), WARN lines (one per warning fit)\n")
for (f in list.files(file.path(rel, "b1"))) {
  cnt <- function(p) {
    if (!file.exists(p)) return(NA)
    x <- readLines(p, warn = FALSE)
    sum(grepl("^   WARN", x))
  }
  cat(sprintf("%-20s lane %s release %s\n", f,
              cnt(file.path(om, "ordmix-p2-b1-log", f)),
              cnt(file.path(rel, "b1", f))))
  cmp(paste0("b1/", f), file.path(om, "ordmix-p2-b1-log", f),
      file.path(rel, "b1", f))
}
cat("\n== B3/B5 weight designs (dev/ordmix-rev3-b3.R)\n")
for (g in c("w_sum1", "w_tenth")) {
  cmp(paste0("b3/", g), file.path(om, "ordmix-p2b-log", paste0("b3-", g, ".txt")),
      file.path(rel, "b3", paste0(g, ".txt")))
}
cat("\n== lane fixes' flat-warning table (dev/fixes-p2-table.R)\n")
tl <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- x[grepl("conv [0-9]", x)]
  sub("[[:space:]]+$", "", x)
}
a <- tl(file.path(fx, "p2-table-wt-fixes-lib.txt"))
b <- tl(file.path(rel, "fixes-p2-table.txt"))
cat(sprintf("rows lane %d release %d identical %d\n", length(a), length(b),
            sum(a %in% b)))
for (d in setdiff(b, a)) cat("  release:", d, "\n")
for (d in setdiff(a, b)) cat("  lane:   ", d, "\n")
cat("flat warnings on the release build:",
    sum(grepl("flat warning YES", b)), "of", length(b), "rows\n")
