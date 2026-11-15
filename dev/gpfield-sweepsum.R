# Lane gpfield: summary of dev/gpfield-sweep.R logs, lane against base by row.
# For the 150 rows, also against lane gpstop's references (best and bestc
# of 15 theta starts, frmtmb-wt-gpstop/dev/gpstop-out/ng-lane-[ab].txt).
# Args: label log [log ...]
a <- commandArgs(TRUE)
lab <- a[1]
rows <- list()
for (f in a[-1]) {
  for (l in grep("^ROW ", readLines(f), value = TRUE)) {
    p <- strsplit(l, " [|] ")[[1]]
    tag <- sub("^ROW ", "", p[1])
    v <- strsplit(p[2], " +")[[1]]
    arm <- v[1]
    if ("ERROR" %in% v) {
      rows[[length(rows) + 1L]] <- data.frame(tag = tag, arm = arm, ll = NA,
        conv = NA, nwarn = NA, wall = NA, runs = NA, warn = "ERROR")
      next
    }
    g <- function(k) v[match(k, v) + 1L]
    rows[[length(rows) + 1L]] <- data.frame(
      tag = tag, arm = arm, ll = as.numeric(g("ll")),
      conv = as.integer(g("conv")), nwarn = as.integer(g("nwarn")),
      wall = as.numeric(g("wall")), runs = g("runs"),
      warn = if (length(p) > 2) trimws(p[3]) else "")
  }
}
d <- do.call(rbind, rows)
L <- d[d$arm == "lane", ]
B <- d[d$arm == "base", ]
m <- merge(L, B, by = "tag", suffixes = c(".l", ".b"))
dl <- m$ll.l - m$ll.b
cat(sprintf(paste0("SUM %s rows %d (lane %d base %d) | up >1e-6 %d >1e-3 %d | ",
                   "down < -1e-6 %d < -1e-3 %d | max up %.4f min %.4f\n"),
            lab, nrow(m), nrow(L), nrow(B), sum(dl > 1e-6, na.rm = TRUE),
            sum(dl > 1e-3, na.rm = TRUE), sum(dl < -1e-6, na.rm = TRUE),
            sum(dl < -1e-3, na.rm = TRUE), max(c(0, dl), na.rm = TRUE),
            min(c(0, dl), na.rm = TRUE)))
cc <- m$conv.l != m$conv.b
wl <- m$warn.l != m$warn.b
cat(sprintf(paste0("SUM %s conv changes %d | warning-set changes %d | ",
                   "wall lane %.1f base %.1f ratio %.3f\n"),
            lab, sum(cc, na.rm = TRUE), sum(wl, na.rm = TRUE),
            sum(m$wall.l, na.rm = TRUE), sum(m$wall.b, na.rm = TRUE),
            sum(m$wall.l, na.rm = TRUE) / sum(m$wall.b, na.rm = TRUE)))
for (i in which(abs(dl) > 1e-6 | cc | wl)) {
  cat(sprintf("MOVE %s %s | lane %.6f c%d w%d [%s] | base %.6f c%d w%d [%s] | d %+.6f\n",
              lab, m$tag[i], m$ll.l[i], m$conv.l[i], m$nwarn.l[i], m$warn.l[i],
              m$ll.b[i], m$conv.b[i], m$nwarn.b[i], m$warn.b[i], dl[i]))
}
# the 150 rows against lane gpstop's references
gs <- "C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-out"
ref <- NULL
for (f in file.path(gs, c("ng-lane-a.txt", "ng-lane-b.txt"))) {
  if (!file.exists(f)) next
  for (l in grep("^ANG ", readLines(f), value = TRUE)) {
    v <- strsplit(l, " +")[[1]]
    if ("ERROR" %in% v) next
    g <- function(k) as.numeric(v[match(k, v) + 1L])
    ref <- rbind(ref, data.frame(tag = paste(v[3], v[4], v[5], v[6]),
                                 best = g("best"), bestc = g("bestc")))
  }
}
if (!is.null(ref) && any(m$tag %in% ref$tag)) {
  r <- merge(m, ref, by = "tag")
  for (arm in c("l", "b")) {
    ll <- r[[paste0("ll.", arm)]]
    sh <- ll < r$best - 1e-3
    shc <- ll < r$bestc - 1e-3
    sil <- r[[paste0("nwarn.", arm)]] == 0
    cat(sprintf(paste0("REF %s arm %s rows %d | short of best %d (silent %d) ",
                       "worst %.3f | short of bestc %d (silent %d) worst %.3f\n"),
                lab, arm, nrow(r), sum(sh), sum(sh & sil),
                min(c(0, ll - r$best)), sum(shc), sum(shc & sil),
                min(c(0, ll - r$bestc))))
  }
  for (i in which(r$ll.l < r$bestc - 1e-3 | r$ll.b < r$bestc - 1e-3)) {
    cat(sprintf("REFSHORT %s %s | lane %.4f w%d [%s] | base %.4f | bestc %.4f best %.4f\n",
                lab, r$tag[i], r$ll.l[i], r$nwarn.l[i], r$warn.l[i], r$ll.b[i],
                r$bestc[i], r$best[i]))
  }
}
