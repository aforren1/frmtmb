# The merged compatibility table against the base, with every extension
# loaded (rule 4 of a move, dev/core-boundary.md). The base is rellib-r11
# (frmtmb 0.73.0 and its seven extensions), the merged build rellib-r12
# (frmtmb 0.74.0 and its eight). Every unordered pair must keep its
# status; its note may differ only on the pairs a lane accepted, and
# there it must be that lane's note: splitsurv's 6
# (dev/splitsurv-compare.R), spatial's 16 (dev/spatial-compat-compare.R).
# A pair's orientation (which feature is feature_a) follows the order of
# the vocabulary, which a contributed feature joins last, so mvnlatent's
# 19 flips (dev/mvnlatent-compat-compare.R) are counted apart; they are
# compared with the flips of every lane's own dumps.
#   Rscript dev/rel074-compat.R      (after the two dumps below)
#   Rscript dev/rel074-compat-dump.R dev/rel074-log/compat-base.rds \
#     C:/Users/adf44/source/r/rellib-r11 <the seven>
#   Rscript dev/rel074-compat-dump.R dev/rel074-log/compat-merged.rds \
#     C:/Users/adf44/source/r/rellib-r12 <the eight>
setwd("C:/Users/adf44/source/r/frmtmb-wt-release")
L <- "C:/Users/adf44/source/r/frmtmb-wt-"
tab_of <- function(x) {
  t <- if (!is.null(x$tab)) x$tab else x$tb
  t$pair <- ifelse(t$feature_a < t$feature_b,
                   paste(t$feature_a, t$feature_b, sep = " x "),
                   paste(t$feature_b, t$feature_a, sep = " x "))
  t$orient <- paste(t$feature_a, t$feature_b, sep = " x ")
  stopifnot(!anyDuplicated(t$pair))
  t[, c("pair", "orient", "status", "note")]
}
cmp <- function(b, a) {
  m <- merge(b, a, by = "pair", all = TRUE, suffixes = c(".b", ".a"))
  both <- m[!is.na(m$status.a) & !is.na(m$status.b), ]
  list(only_b = m$pair[is.na(m$status.a)], only_a = m$pair[is.na(m$status.b)],
       status = both$pair[both$status.b != both$status.a],
       note = both$pair[both$note.b != both$note.a],
       flip = both$pair[both$orient.b != both$orient.a], m = both)
}
base <- readRDS("dev/rel074-log/compat-base.rds")
mrg <- readRDS("dev/rel074-log/compat-merged.rds")
cat("base:  ", paste(base$libs, collapse = "\n        "), "\n")
cat("merged:", paste(mrg$libs, collapse = "\n        "), "\n")
B <- tab_of(base); M <- tab_of(mrg)
r <- cmp(B, M)
cat("pairs", nrow(B), "->", nrow(M), "| only base", length(r$only_b),
    "| only merged", length(r$only_a), "\n")
fb <- base$features[order(base$features$name), c("name", "kind")]
fm <- mrg$features[order(mrg$features$name), c("name", "kind")]
rownames(fb) <- rownames(fm) <- NULL
cat("features", nrow(fb), "->", nrow(fm), "| identical (name, kind):",
    identical(fb, fm), "\n")
if (!identical(fb, fm)) {
  mm <- merge(fb, fm, by = "name", all = TRUE)
  print(mm[is.na(mm$kind.x) | is.na(mm$kind.y) | mm$kind.x != mm$kind.y, ])
}
cat("status differs:", length(r$status), "\n")
if (length(r$status)) print(r$m[r$m$pair %in% r$status,
                                c("pair", "status.b", "status.a")])

# what each lane accepted, from its own dumps
lane <- list(
  splitsurv = cmp(tab_of(readRDS(paste0(L, "splitsurv/dev/splitsurv-out/compat-before.rds"))),
                  tab_of(readRDS(paste0(L, "splitsurv/dev/splitsurv-out/compat-after.rds")))),
  mvnlatent = cmp(tab_of(readRDS(paste0(L, "mvnlatent/dev/mvnlatent-out/compat-base-all.rds"))),
                  tab_of(readRDS(paste0(L, "mvnlatent/dev/mvnlatent-out/compat-lane-all.rds")))),
  spatial = cmp(tab_of(readRDS(paste0(L, "spatial/dev/spatial-log/compat-before-all.rds"))),
                tab_of(readRDS(paste0(L, "spatial/dev/spatial-log/compat-f1b-last.rds")))))
for (l in names(lane)) {
  x <- lane[[l]]
  cat(sprintf("lane %-9s status differs %d | note differs %d | flips %d\n",
              l, length(x$status), length(x$note), length(x$flip)))
}
# splitsurv also reworded one rule of frmtmb.learn ("core::cox() and
# frmtmb.spline::royston_parmar()" now "frmtmb.survival::cox() and
# royston_parmar()", dev/splitsurv-findings.md, "Overlaps"), which its
# dumps (core and frmtmb.spline or frmtmb.survival) do not hold: the
# pairs whose note differs by exactly that rewording are accepted too
lo <- "core::cox() and frmtmb.spline::royston_parmar() decline"
ln <- "frmtmb.survival::cox() and royston_parmar() decline"
lr <- r$m[r$m$pair %in% r$note & grepl(lo, r$m$note.b, fixed = TRUE), ]
learn_ok <- lr$pair[gsub(lo, ln, lr$note.b, fixed = TRUE) == lr$note.a]
cat("frmtmb.learn pairs whose note differs by exactly splitsurv's",
    "rewording:", length(learn_ok), "\n")
acc <- unique(c(lane$splitsurv$note, lane$spatial$note, lane$mvnlatent$note,
                learn_ok))
cat("notes accepted by the lanes (union):", length(acc), "\n")
cat("merged note differs:", length(r$note), "| of them accepted:",
    sum(r$note %in% acc), "| not accepted:", sum(!r$note %in% acc),
    "| accepted but unchanged in the merge:", sum(!acc %in% r$note), "\n")
if (any(!r$note %in% acc)) {
  for (p in setdiff(r$note, acc)) {
    z <- r$m[r$m$pair == p, ]
    cat("  NOT ACCEPTED:", p, "\n    base:  ", substr(z$note.b, 1, 150),
        "\n    merged:", substr(z$note.a, 1, 150), "\n")
  }
}
if (any(!acc %in% r$note)) cat("  unchanged:", setdiff(acc, r$note), sep = "\n   ")
# the merged note is the lane's own note on each accepted pair
same_as_lane <- function(l) {
  x <- lane[[l]]$m
  p <- intersect(lane[[l]]$note, r$note)
  z <- merge(x[x$pair %in% p, c("pair", "note.a")],
             r$m[r$m$pair %in% p, c("pair", "note.a")], by = "pair")
  c(n = nrow(z), equal = sum(z$note.a.x == z$note.a.y))
}
for (l in c("splitsurv", "spatial")) {
  s <- same_as_lane(l)
  cat(sprintf("  %s's accepted pairs: %d, merged note equal to the lane's: %d\n",
              l, s[["n"]], s[["equal"]]))
}
# A flip changes neither status nor note of the unordered pair; it says
# that the vocabulary order moved (a contributed feature joins it last,
# and curve inference is core's now). Each lane's dumps hold a
# different set of extensions, so their flips are only a guide; the
# merged flips are listed by feature and do not enter the verdict
fl <- unique(unlist(lapply(lane, `[[`, "flip")))
cat("orientation flips: merged", length(r$flip), "| lanes' union",
    length(fl), "| merged flips in the union:", sum(r$flip %in% fl),
    "| outside:", sum(!r$flip %in% fl), "\n")
if (any(!r$flip %in% fl)) {
  ff <- r$flip[!r$flip %in% fl]
  cat("  features of the flips outside the lanes' union:\n")
  print(sort(table(unlist(strsplit(ff, " x ", fixed = TRUE))),
             decreasing = TRUE)[1:min(10, length(ff))])
}
ok <- length(r$only_b) == 0 && length(r$only_a) == 0 && identical(fb, fm) &&
  length(r$status) == 0 && all(r$note %in% acc)
cat("COMPAT VERDICT:", if (ok) "PASS" else "FAIL", "\n")
