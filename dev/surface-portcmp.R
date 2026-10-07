# Lane surface: the mechanical port's spell pass and the hand
# translation on rellib-r6 (0.68.1) against the lane build, row by row.
#
#   Rscript dev/surface-portcmp.R > dev/surface-out/portcmp.txt
#
# Inputs: dev/surface-port-out/{r6,lane}/results-merged.csv (written by
# dev/brms-port/summarize.R after run-all.R raw and spell) and
# dev/surface-bv-out/{r6,lane}/brms*.csv (dev/brms-vignettes/_run-all.sh).
one_line <- function(s) substr(gsub("[\r\n]+", " ", s), 1, 110)
m <- lapply(c(r6 = "r6", lane = "lane"), function(t) {
  utils::read.csv(file.path("dev/surface-port-out", t, "results-merged.csv"),
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
b <- lapply(c(r6 = "r6", lane = "lane"), function(t) {
  fs <- list.files(file.path("dev/surface-bv-out", t), "^brms.*[.]csv$",
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
