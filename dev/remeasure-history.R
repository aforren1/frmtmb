# The outcome of every vignette call across builds, so a call that ran on
# an older build and fails on the newest is found by construction. The
# successor of dev/vigport-history.R, whose input logs no longer exist.
#
#   Rscript dev/remeasure-history.R <tag> [<tag> ...]
#
# Each <tag> is a run of dev/remeasure-run.sh (dev/remeasure-out/<tag>,
# after summarize.R). Its outcomes are added to the tracked tables
# dev/brms-port/outcomes-mech.tsv and outcomes-hand.tsv, replacing a
# column of the same tag, and the build is listed in outcome-builds.tsv.
# The tables were seeded by dev/remeasure-seed-history.R. The LAST tag
# named is the build judged: a row that ran on any build from 0.67.0 on
# and fails on it is a regression.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
bp <- file.path(root, "brms-port")
tags <- commandArgs(trailingOnly = TRUE)
stopifnot(length(tags) >= 1)
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE,
                                    check.names = FALSE, quote = "",
                                    na.strings = "")
M <- rd(file.path(bp, "outcomes-mech.tsv"))
H <- rd(file.path(bp, "outcomes-hand.tsv"))
Bt <- rd(file.path(bp, "outcome-builds.tsv"))
one <- function(s, n = 110) substr(gsub("[\r\n\t]+", " ", s), 1, n)
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
cur <- list()
for (t in tags) {
  d <- file.path(root, "remeasure-out", t)
  R <- readRDS(file.path(d, "port", "results-merged.rds"))
  stopifnot(setequal(R$id, M$id))
  R <- R[match(M$id, R$id), ]
  stopifnot(identical(substr(gsub("[\t\r\n]", " ", R$src), 1, 100),
                      M$src))
  hdr <- readLines(file.path(d, "port", "results", "brms_overview.log"),
                   n = 1)
  ver <- sub("^# frmtmb ([^ ]+) from ([^,]+),.*$", "\\1", hdr)
  lib <- sub("^# frmtmb ([^ ]+) from ([^,]+),.*$", "\\2", hdr)
  M[[paste0(t, ".raw")]] <- R$status_raw
  M[[paste0(t, ".spell")]] <- R$status_spell
  M[[paste0(t, ".class")]] <- R$class
  Hd <- read_bv(file.path(d, "bv"))
  stopifnot(setequal(Hd$key, H$key))
  Hd <- Hd[match(H$key, Hd$key), ]
  H[[paste0(t, ".ok")]] <- Hd$ok
  Bt <- Bt[Bt$tag != t, ]
  Bt <- rbind(Bt, data.frame(tag = t, version = ver,
                             source = paste("fresh run from", lib,
                                            "dev/remeasure-run.sh")))
  cur[[t]] <- list(R = R, H = Hd)
}
utils::write.table(M, file.path(bp, "outcomes-mech.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
utils::write.table(H, file.path(bp, "outcomes-hand.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE, na = "")
utils::write.table(Bt, file.path(bp, "outcome-builds.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)

# Columns in build order: the seeded ones, then the tags as named.
all_tags <- Bt$tag
new <- tags[length(tags)]
lab <- stats::setNames(paste0(Bt$version, "(", Bt$tag, ")"), Bt$tag)
cat("## builds\n\n")
print(Bt, row.names = FALSE, right = FALSE)
cat("\njudged build:", lab[[new]], "\n")

cat("\n## totals by build\n\n")
has <- function(t, f) paste0(t, f) %in% names(M)
cnt <- function(t, kind, f, val) {
  if (!has(t, f)) return(NA_integer_)
  x <- M[[paste0(t, f)]][M$kind == kind]
  if (identical(f, ".class")) sum(sub(":.*$", "", x) == val) else
    sum(x == val)
}
tab <- data.frame(build = unname(lab[all_tags]),
  model_raw_ok = vapply(all_tags, cnt, 0L, "model", ".raw", "OK"),
  model_clean = vapply(all_tags, cnt, 0L, "model", ".class", "CLEAN"),
  model_spell_ok = vapply(all_tags, cnt, 0L, "model", ".spell", "OK"),
  post_raw_ok = vapply(all_tags, cnt, 0L, "post", ".raw", "OK"),
  post_clean = vapply(all_tags, cnt, 0L, "post", ".class", "CLEAN"),
  post_spell_ok = vapply(all_tags, cnt, 0L, "post", ".spell", "OK"),
  hand_ok = vapply(all_tags, function(t) {
    if (paste0(t, ".ok") %in% names(H)) sum(H[[paste0(t, ".ok")]]) else
      NA_integer_
  }, 0L), row.names = NULL)
cat("model calls 42, post calls 102, hand rows", nrow(H), "\n\n")
print(tab, row.names = FALSE)

# From 0.67.0 on every build was run fresh, so these columns are complete.
since <- all_tags[seq(match("v067", all_tags), length(all_tags))]
cat("\n## mechanical port: rows whose raw or spell status changed since",
    "0.67.0\n\n")
mk <- function(f) vapply(since, function(t) M[[paste0(t, f)]], M$id)
sr <- mk(".raw"); ss <- mk(".spell")
chg <- which(apply(sr, 1, function(z) length(unique(z)) > 1) |
               apply(ss, 1, function(z) length(unique(z)) > 1))
cat("rows:", length(chg), "of", nrow(M), "\n")
cat("columns: raw/spell per build,", paste(lab[since], collapse = " "),
    "\n\n")
Rn <- cur[[new]]$R
for (i in chg) {
  cat(sprintf("%-24s %-5s %s  %s\n", M$id[i], M$kind[i],
              paste(sprintf("%s/%s", substr(sr[i, ], 1, 5),
                            substr(ss[i, ], 1, 5)), collapse = " "),
              one(M$src[i], 70)))
  if (ss[i, ncol(ss)] != "OK") {
    cat("      now: ", one(Rn$msg_spell[i], 150), "\n")
  }
}
cat("\n## hand translation: rows whose outcome changed since 0.67.0\n\n")
ho <- vapply(since, function(t) H[[paste0(t, ".ok")]], logical(nrow(H)))
hchg <- which(apply(ho, 1, function(z) length(unique(z)) > 1))
cat("rows:", length(hchg), "of", nrow(H), "\n")
cat("columns: ok per build,", paste(lab[since], collapse = " "), "\n\n")
Hn <- cur[[new]]$H
for (i in hchg) {
  cat(sprintf("%s  [%s] %s\n",
              paste(ifelse(ho[i, ], "ok ", "ERR"), collapse = " "),
              H$edge[i], one(H$key[i], 100)))
  if (!ho[i, ncol(ho)]) cat("      now: ", one(Hn$msg[i], 150), "\n")
}

cat("\n## regressions: ran on a build from 0.67.0 on, fails on",
    lab[[new]], "\n\n")
ran_m <- apply(sr == "OK" | ss == "OK", 1, any)
now_m <- sr[, ncol(sr)] == "OK" | ss[, ncol(ss)] == "OK"
reg_m <- which(ran_m & !now_m)
cat("mechanical (raw or spell pass):", length(reg_m), "\n")
for (i in reg_m) {
  cat(sprintf("  %-24s %s\n      cause: %s\n", M$id[i], one(M$src[i], 70),
              one(Rn$msg_spell[i], 200)))
}
# a raw pass that ran and now needs a patch is a regression of CLEAN
reg_c <- which(apply(sr == "OK", 1, any) & sr[, ncol(sr)] != "OK")
cat("mechanical raw pass only (CLEAN lost):", length(reg_c), "\n")
for (i in reg_c) {
  cat(sprintf("  %-24s %s\n      cause: %s\n", M$id[i], one(M$src[i], 70),
              one(Rn$msg_raw[i], 200)))
}
reg_h <- which(apply(ho, 1, any) & !ho[, ncol(ho)])
cat("hand translation:", length(reg_h), "\n")
for (i in reg_h) {
  cat(sprintf("  [%s] %s\n      cause: %s\n", H$edge[i], one(H$key[i], 100),
              one(Hn$msg[i], 200)))
}

cat("\n## against 0.34.0: ran then (raw or behind a patch), fail now\n\n")
b34 <- sub(":.*$", "", M$v034.class)
bnow <- sub(":.*$", "", M[[paste0(new, ".class")]])
lost <- which(M$kind %in% c("model", "post") &
                b34 %in% c("CLEAN", "SPELLING") &
                bnow %in% c("FAIL", "CASCADE"))
cat("rows:", length(lost), "\n")
for (i in lost) {
  cat(sprintf("  %-24s %-18s %s\n      now: %s\n", M$id[i], M$v034.class[i],
              one(M$src[i], 60), one(Rn$msg_spell[i], 150)))
}
cat("\n## end\n")
