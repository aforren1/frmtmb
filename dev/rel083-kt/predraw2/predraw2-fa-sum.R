# Lane predraw2, fix round: the false-alarm table of the two-condition
# check (median excess above draw_check_tol AND its ratio to the median
# drop the law implies above draw_check_rel). Reads the K designs
# (dev/predraw2-kd-*.tsv), the corpus (dev/predraw2-meas3-*-lane.rds) and
# the mixture coverage runs (dev/predraw2-mc42-*-lane.tsv). Per design:
# fits with a verdict, fired, the largest excess and ratio of the unfired
# fits, the smallest of the fired ones, and the least distance of any
# verdict from the rule's boundary, max over the two conditions of
# |log(x / threshold)| for the condition that decides it.
# Usage: Rscript dev/predraw2-fa-sum.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
tol <- ns$draw_check_tol
rtol <- ns$draw_check_rel
cat("draw_check_tol", tol, "draw_check_rel", rtol, "\n")
rows <- list()
for (f in Sys.glob("dev/rel083-log/kt/predraw2/predraw2-kd-*.tsv")) {
  x <- utils::read.delim(f)
  rows[[f]] <- data.frame(design = x$design, excess = x$excess, rel = x$rel,
                          fire = x$fire)
}
for (f in Sys.glob("dev/rel083-log/kt/predraw2/predraw2-meas3-*-lane.rds")) {
  a <- readRDS(f)
  id <- vapply(a, `[[`, "", "id")
  g <- sub(":.*", "", id)
  g[startsWith(id, "sweep:")] <- paste0("sweep_", sub(".*:", "", id[
    startsWith(id, "sweep:")]))
  num <- function(k) vapply(a, function(z) z[[k]] %||% NA_real_, 0)
  rows[[f]] <- data.frame(design = paste0("corpus_", g), excess = num("excess"),
                          rel = num("rel"),
                          fire = vapply(a, function(z) isTRUE(z$fire), TRUE))
}
for (f in Sys.glob("dev/rel083-log/kt/predraw2/predraw2-mc42-*-lane.tsv")) {
  x <- utils::read.delim(f)
  rows[[f]] <- data.frame(design = paste0("mix40_", sub(".*mc42-(.)-.*", "\\1",
                                                       f)),
                          excess = x$excess, rel = x$rel, fire = x$fire)
}
d <- do.call(rbind, rows)
d <- d[!is.na(d$excess), ]
# the boundary distance: a verdict decided by one condition is as far
# from flipping as that condition is from its threshold
dist <- function(ex, rel) {
  a <- log(pmax(ex, 1e-12) / tol)
  b <- log(pmax(rel, 1e-12) / rtol)
  ifelse(a > 0 & b > 0, pmin(a, b), pmax(-a * (a <= 0), -b * (b <= 0)))
}
d$dist <- dist(d$excess, d$rel)
fmt <- function(v) if (length(v)) sprintf("%.4g", v) else "-"
cat(sprintf("%-22s %4s %5s | %-22s | %-22s | %s\n", "design", "fits", "fired",
            "unfired max excess/rel", "fired min excess/rel",
            "least boundary distance (log)"))
for (g in sort(unique(d$design))) {
  x <- d[d$design == g, ]
  u <- x[!x$fire, ]
  fi <- x[x$fire, ]
  cat(sprintf("%-22s %4d %5d | %-22s | %-22s | %s\n", g, nrow(x), nrow(fi),
              paste(fmt(max(c(-Inf, u$excess))), fmt(max(c(-Inf, u$rel))),
                    sep = " / "),
              if (nrow(fi)) paste(fmt(min(fi$excess)), fmt(min(fi$rel)),
                                  sep = " / ") else "-",
              fmt(min(x$dist))))
}
cat(sprintf("all: %d verdicts, %d fired; least boundary distance %.3f (%s)\n",
            nrow(d), sum(d$fire), min(d$dist), d$design[which.min(d$dist)]))
