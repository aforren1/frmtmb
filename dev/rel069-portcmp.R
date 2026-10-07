# dev/surface-portcmp.R for the 0.69.0 consolidation: the mechanical
# port's spell pass and the hand translation of two builds, row by row.
# "r6" is rellib-r6 and "lane" lane surface's build, both as lane surface
# recorded them in its worktree; "rel" is the release library rellib-r7
# (bash dev/surface-port.sh rel C:/Users/adf44/source/r/rellib-r7).
#
#   Rscript dev/rel069-portcmp.R r6 rel  > dev/rel069-log/portcmp-r6.txt
#   Rscript dev/rel069-portcmp.R lane rel > dev/rel069-log/portcmp-lane.txt
dirs <- list(
  r6 = c("C:/Users/adf44/source/r/frmtmb-wt-surface/dev/surface-port-out/r6", "C:/Users/adf44/source/r/frmtmb-wt-surface/dev/surface-bv-out/r6"),
  lane = c("C:/Users/adf44/source/r/frmtmb-wt-surface/dev/surface-port-out/lane", "C:/Users/adf44/source/r/frmtmb-wt-surface/dev/surface-bv-out/lane"),
  rel = c("dev/surface-port-out/rel", "dev/surface-bv-out/rel"))
arms <- commandArgs(trailingOnly = TRUE)
one_line <- function(s) substr(gsub("[\r\n]+", " ", s), 1, 110)
m <- lapply(stats::setNames(arms, c("r6", "lane")), function(t) {
  utils::read.csv(file.path(dirs[[t]][1], "results-merged.csv"),
                  stringsAsFactors = FALSE)
})
stopifnot(identical(m$r6$id, m$lane$id))
cat("## mechanical port, spell pass\n\n")
cat("rows:", nrow(m$r6), "\n")
for (k in c("model", "post", "other")) {
  i <- m$r6$kind == k
  cat(sprintf("%-5s spell OK: r6 %d -> lane %d of %d; class CLEAN: r6 %d -> lane %d\n",
              k, sum(m$r6$status_spell[i] == "OK"),
              sum(m$lane$status_spell[i] == "OK"), sum(i),
              sum(m$r6$class[i] == "CLEAN"), sum(m$lane$class[i] == "CLEAN")))
}
cat("\nheadline class counts (model and post rows), r6 -> lane:\n")
for (k in c("model", "post")) {
  i <- m$r6$kind == k
  for (cl in c("CLEAN", "SPELLING", "FAIL", "CASCADE")) {
    cat(sprintf("  %-5s %-8s %3d -> %3d\n", k, cl, sum(m$r6$class[i] == cl),
                sum(m$lane$class[i] == cl)))
  }
}
ch <- which(m$r6$status_spell != m$lane$status_spell |
              m$r6$msg_spell != m$lane$msg_spell)
cat("\nrows whose spell status or message changed:", length(ch), "\n")
for (i in ch) {
  cat(sprintf("%-26s %-5s %s -> %s  %s\n", m$r6$id[i], m$r6$kind[i],
              m$r6$status_spell[i], m$lane$status_spell[i],
              one_line(m$r6$src[i])))
  if (nzchar(m$r6$msg_spell[i])) {
    cat("    r6:   ", one_line(m$r6$msg_spell[i]), "\n")
  }
  if (nzchar(m$lane$msg_spell[i])) {
    cat("    lane: ", one_line(m$lane$msg_spell[i]), "\n")
  }
}
reg <- which(m$r6$status_spell == "OK" & m$lane$status_spell != "OK")
cat("\nregressions (spell OK on r6, not on the lane):", length(reg), "\n")

cat("\n## hand translation\n\n")
b <- lapply(stats::setNames(arms, c("r6", "lane")), function(t) {
  fs <- list.files(dirs[[t]][2], "^brms.*[.]csv$",
                   full.names = TRUE)
  d <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
  d[d$kind != "data", ]
})
key <- function(d) paste(d$vignette, d$label, sep = " | ")
# a label can repeat within a vignette, so the key carries its order
kk <- function(d) paste(key(d), stats::ave(seq_along(key(d)), key(d),
                                           FUN = seq_along))
stopifnot(setequal(kk(b$r6), kk(b$lane)))
l <- b$lane[match(kk(b$r6), kk(b$lane)), ]
cat("rows:", nrow(b$r6), "| ok: r6", sum(b$r6$ok), "-> lane", sum(l$ok), "\n")
chb <- which(b$r6$ok != l$ok | b$r6$msg != l$msg)
cat("rows whose outcome or message changed:", length(chb), "\n")
for (i in chb) {
  cat(sprintf("%s -> %s  [%s] %s | %s\n", b$r6$ok[i], l$ok[i],
              b$r6$edge[i], b$r6$vignette[i], one_line(b$r6$label[i])))
  if (nzchar(b$r6$msg[i])) cat("    r6:   ", one_line(b$r6$msg[i]), "\n")
  if (nzchar(l$msg[i])) cat("    lane: ", one_line(l$msg[i]), "\n")
}
cat("\nregressions (ok on r6, not on the lane):",
    sum(b$r6$ok & !l$ok), "\n")
cat("\n(arm \"r6\" above is", arms[1], "and \"lane\" is", arms[2], ")\n")
