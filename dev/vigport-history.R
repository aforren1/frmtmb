# The outcome of every vignette call across three release builds, so a
# call that ran on an older build and fails on 0.67.0 is found by
# construction and not by reading.
#
#   Rscript dev/vigport-history.R
#
# Reads the mechanical port's merged records and the hand translation's
# CSVs written for rellib-r3 (0.65.0), rellib-r4 (0.66.0) and rellib-r5
# (0.67.0) under dev/vigport-port-out/<r> and dev/vigport-bv-out/<r>.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
builds <- c(r3 = "0.65.0", r4 = "0.66.0", r5 = "0.67.0")

cat("## mechanical port, raw pass status by build (model and post rows)\n\n")
M <- lapply(names(builds), function(r)
  readRDS(file.path(root, "vigport-port-out", r, "results-merged.rds")))
names(M) <- names(builds)
B <- readRDS(file.path(root, "brms-port", "results-merged.rds"))
ids <- M$r5$id[M$r5$kind %in% c("model", "post")]
st <- sapply(names(builds), function(r)
  M[[r]]$status_raw[match(ids, M[[r]]$id)])
st034 <- B$status_raw[match(ids, B$id)]
ok <- st == "OK"
cat("raw OK, model rows:",
    paste0(c("0.34.0", builds), "=",
           c(sum(st034[M$r5$kind[match(ids, M$r5$id)] == "model"] == "OK"),
             colSums(ok[M$r5$kind[match(ids, M$r5$id)] == "model", ])),
           collapse = "  "), "\n")
cat("raw OK, post rows: ",
    paste0(c("0.34.0", builds), "=",
           c(sum(st034[M$r5$kind[match(ids, M$r5$id)] == "post"] == "OK"),
             colSums(ok[M$r5$kind[match(ids, M$r5$id)] == "post", ])),
           collapse = "  "), "\n\n")
chg <- which(apply(cbind(st034 == "OK", ok), 1,
                   function(z) length(unique(z)) > 1))
cat("rows whose raw outcome differs between builds:", length(chg), "\n")
for (i in chg) {
  cat(sprintf("%-28s 0.34.0=%-6s 0.65.0=%-6s 0.66.0=%-6s 0.67.0=%-6s %s\n",
              ids[i], st034[i], st[i, "r3"], st[i, "r4"], st[i, "r5"],
              substr(M$r5$src[M$r5$id == ids[i]], 1, 60)))
}
lost <- ids[st034 == "OK" & !ok[, "r5"]]
cat("\nran raw at 0.34.0 and fail raw at 0.67.0:", length(lost), "\n")
for (id in lost) {
  first_bad <- names(builds)[which(!ok[ids == id, ])[1]]
  cat(sprintf("%-28s failing from %s or earlier: %s\n", id,
              builds[first_bad],
              substr(M$r5$msg_raw[M$r5$id == id], 1, 120)))
}
# The wider count: a row that ran at 0.34.0 in either pass (CLEAN, or
# SPELLING behind a patch) and fails in both passes now (FAIL or
# CASCADE). A row that ran only behind an upstream patch then is lost
# just the same.
b34 <- B$bucket[match(ids, B$id)]
b67 <- M$r5$bucket[match(ids, M$r5$id)]
lost2 <- ids[b34 %in% c("CLEAN", "SPELLING") & b67 %in% c("FAIL", "CASCADE")]
cat("\nran at 0.34.0 (raw or patched) and fail raw and patched at 0.67.0:",
    length(lost2), "\n")
for (id in lost2) {
  first_bad <- names(builds)[which(!ok[ids == id, ])[1]]
  cat(sprintf("%-28s %-18s failing raw from %s or earlier: %s\n", id,
              B$class[B$id == id], builds[first_bad],
              substr(M$r5$msg_spell[M$r5$id == id], 1, 100)))
}

cat("\n## hand translation, ok by build\n\n")
H <- lapply(names(builds), function(r) {
  d <- do.call(rbind, lapply(list.files(file.path(root, "vigport-bv-out", r),
                                        "^brms.*[.]csv$", full.names = TRUE),
                             utils::read.csv, stringsAsFactors = FALSE))
  d$edge[is.na(d$edge) | d$edge == ""] <- "CLEAN"
  d <- d[d$kind != "data", ]
  d$key <- paste(d$vignette, ave(seq_along(d$label), d$vignette,
                                 FUN = seq_along), d$label, sep = " | ")
  d
})
names(H) <- names(builds)
keys <- H$r5$key
okh <- sapply(names(builds), function(r) H[[r]]$ok[match(keys, H[[r]]$key)])
cat("rows:", length(keys), "  ok:",
    paste0(builds, "=", colSums(okh, na.rm = TRUE), collapse = "  "),
    "  missing keys:", sum(is.na(okh)), "\n")
flip <- which(apply(okh, 1, function(z) length(unique(z)) > 1))
cat("rows whose outcome differs between builds:", length(flip), "\n")
for (i in flip) {
  cat(sprintf("%-6s %-6s %-6s [%s] %s\n", okh[i, 1], okh[i, 2], okh[i, 3],
              H$r5$edge[i], substr(keys[i], 1, 110)))
  if (!okh[i, 3]) cat("       now: ", substr(H$r5$msg[i], 1, 160), "\n")
}
