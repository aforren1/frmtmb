# One-time: seed the tracked outcome tables that dev/remeasure-history.R
# extends, from the records that survive.
#
#   Rscript dev/remeasure-seed-history.R
#
# Why: the per-row logs of the 0.65.0, 0.66.0 and 0.67.0 measurements
# were deleted with the untracked dev logs, and the release libraries
# rellib-r3 to rellib-r5 no longer exist. What survives is (a) the
# tracked 0.34.0 records in dev/brms-port, (b) the generated history
# block in dev/vigport-findings.md, which lists EVERY row whose outcome
# differs between 0.34.0, 0.65.0, 0.66.0 and 0.67.0, and (c) a fresh
# 0.67.0 run (frmtmb 0.67.0 and frmtmb.sample 0.15.0 installed from
# 9e902909 into the lane library, tag v067). The 0.65.0 and 0.66.0
# columns are rebuilt from (b) and checked against the counts the block
# printed; the rebuilt 0.67.0 column must equal the fresh run row for
# row, which checks the rebuild and the rerun at once.
#
# Writes dev/brms-port/outcomes-mech.tsv, dev/brms-port/outcomes-hand.tsv
# and dev/brms-port/outcome-builds.tsv (tracked: dev/brms-port/** is).
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
bp <- file.path(root, "brms-port")
out067 <- file.path(root, "remeasure-out", "v067")
rec <- readLines(file.path(root, "vigport-findings.md"), warn = FALSE)
blk <- rec[(grep("BEGIN generated: history", rec) + 1):
             (grep("END generated", rec)[grep("END generated", rec) >
                                           grep("BEGIN generated: history",
                                                rec)][1] - 1)]

## mechanical port --------------------------------------------------------
B <- readRDS(file.path(bp, "results-merged.rds"))
V <- readRDS(file.path(out067, "port", "results-merged.rds"))
stopifnot(identical(sort(B$id), sort(V$id)))
ids <- V$id
same_src <- B$src[match(ids, B$id)] == V$src
cat("mechanical rows:", length(ids), " same transformed source as 0.34.0:",
    sum(same_src), "\n")
re <- paste0("^(brms_[A-Za-z_]+[.][0-9]+[.][0-9]+) +0[.]34[.]0=([A-Z-]+) +",
             "0[.]65[.]0=([A-Z-]+) +0[.]66[.]0=([A-Z-]+) +0[.]67[.]0=([A-Z-]+)")
hl <- grep(re, blk, value = TRUE)
hm <- regmatches(hl, regexec(re, hl))
H <- do.call(rbind, lapply(hm, function(m) m[2:6]))
colnames(H) <- c("id", "v034", "v065", "v066", "v067")
n_listed <- as.integer(sub(".*: ", "", grep("^rows whose raw outcome differs",
                                            blk, value = TRUE)))
cat("rows listed in the record:", nrow(H), " the record says:", n_listed,
    "\n")
stopifnot(nrow(H) == n_listed)
ok034 <- B$status_raw[match(ids, B$id)] == "OK"
rebuild <- function(col) {
  v <- ok034
  i <- match(H[, "id"], ids)
  v[i] <- H[, col] == "OK"
  v
}
ok065 <- rebuild("v065"); ok066 <- rebuild("v066"); ok067 <- rebuild("v067")
mp <- V$kind %in% c("model", "post")
# the counts the record printed, for the model and post rows
want <- grep("^raw OK, (model|post) rows", blk, value = TRUE)
cat(want, sep = "\n")
for (k in c("model", "post")) {
  cat(sprintf("rebuilt raw OK, %-5s 0.34.0=%d 0.65.0=%d 0.66.0=%d 0.67.0=%d\n",
              k, sum(ok034[V$kind == k]), sum(ok065[V$kind == k]),
              sum(ok066[V$kind == k]), sum(ok067[V$kind == k])))
}
fresh <- V$status_raw == "OK"
cat("rebuilt 0.67.0 against the fresh 0.67.0 run, model and post rows:",
    sum(ok067[mp] != fresh[mp]), "rows differ\n")
diff_other <- ids[!mp & ok067 != fresh]
cat("other rows (not in the record's block) where 0.34.0 and the fresh",
    "0.67.0 run differ:", length(diff_other), "\n")
stopifnot(sum(ok067[mp] != fresh[mp]) == 0)
okstr <- function(v) ifelse(v, "OK", "NOT-OK")
M <- data.frame(
  id = ids, vignette = V$vignette, kind = V$kind,
  src = substr(gsub("[\t\r\n]", " ", V$src), 1, 100),
  v034.raw = B$status_raw[match(ids, B$id)],
  v034.spell = B$status_spell[match(ids, B$id)],
  v034.class = B$class[match(ids, B$id)],
  v065.raw = ifelse(mp, okstr(ok065), NA),
  v066.raw = ifelse(mp, okstr(ok066), NA),
  v067.raw = V$status_raw, v067.spell = V$status_spell, v067.class = V$class,
  check.names = FALSE, stringsAsFactors = FALSE)

## hand translation -------------------------------------------------------
read_bv <- function(dir) {
  d <- do.call(rbind, lapply(list.files(dir, "^brms.*[.]csv$",
                                        full.names = TRUE),
                             utils::read.csv, stringsAsFactors = FALSE))
  d$edge[is.na(d$edge) | d$edge == ""] <- "CLEAN"
  d <- d[d$kind != "data", ]
  d$key <- paste(d$vignette, stats::ave(seq_along(d$label), d$vignette,
                                        FUN = seq_along), d$label,
                 sep = " | ")
  d
}
Hd <- read_bv(file.path(out067, "bv"))
cat("\nhand rows:", nrow(Hd), " ok in the fresh 0.67.0 run:", sum(Hd$ok),
    "\n")
reh <- "^(TRUE|FALSE) +(TRUE|FALSE) +(TRUE|FALSE) +\\[[A-Z]+\\] (.*)$"
fl <- grep(reh, blk, value = TRUE)
fm <- regmatches(fl, regexec(reh, fl))
cat("hand rows listed in the record as differing between builds:",
    length(fm), "\n")
h065 <- h066 <- h067 <- Hd$ok
for (m in fm) {
  i <- which(startsWith(Hd$key, trimws(m[5])))
  stopifnot(length(i) == 1L)
  h065[i] <- as.logical(m[2]); h066[i] <- as.logical(m[3])
  h067[i] <- as.logical(m[4])
}
cat(grep("^rows: 539", blk, value = TRUE), "\n")
cat(sprintf("rebuilt ok: 0.65.0=%d 0.66.0=%d 0.67.0=%d\n", sum(h065),
            sum(h066), sum(h067)))
cat("rebuilt 0.67.0 against the fresh run:", sum(h067 != Hd$ok),
    "rows differ\n")
stopifnot(sum(h067 != Hd$ok) == 0)
Ht <- data.frame(key = Hd$key, vignette = Hd$vignette, kind = Hd$kind,
                 edge = Hd$edge, v065.ok = h065, v066.ok = h066,
                 v067.ok = Hd$ok, stringsAsFactors = FALSE,
                 check.names = FALSE)

builds <- data.frame(
  tag = c("v034", "v065", "v066", "v067"),
  version = c("0.34.0", "0.65.0", "0.66.0", "0.67.0"),
  source = c("dev/brms-port/results-merged.rds (tracked)",
             "rebuilt from the dev/vigport-findings.md history block",
             "rebuilt from the dev/vigport-findings.md history block",
             paste("fresh run, frmtmb 0.67.0 and frmtmb.sample 0.15.0 from",
                   "9e902909, dev/remeasure-run.sh")),
  stringsAsFactors = FALSE)
utils::write.table(M, file.path(bp, "outcomes-mech.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
utils::write.table(Ht, file.path(bp, "outcomes-hand.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
utils::write.table(builds, file.path(bp, "outcome-builds.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("\nwrote outcomes-mech.tsv (", nrow(M), " rows), outcomes-hand.tsv (",
    nrow(Ht), " rows), outcome-builds.tsv\n", sep = "")
