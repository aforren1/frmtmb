# Reduce the runs to the scorecard tables.
#
#   PORT_OUT=<dir> Rscript summarize.R [baseline.rds]
#
# Reads <PORT_OUT>/results and <PORT_OUT>/results-spell and writes
# results-merged.rds and .csv next to them. With a baseline (the
# results-merged.rds of an earlier measurement, by default the 0.34.0
# one tracked in this directory), it also prints how every expression's
# class moved between the two.
HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
source(file.path(HERE, "env.R"))
source(file.path(HERE, "port-lib.R"))
source(file.path(HERE, "patches.R"))
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args)) args[1] else file.path(HERE, "results-merged.rds")
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")
`%||%` <- function(a, b) if (is.null(a)) b else a

load_dir <- function(d) {
  out <- list()
  for (v in VIGS) {
    f <- file.path(PORT_OUT, d, paste0(v, ".rds"))
    if (file.exists(f)) out <- c(out, readRDS(f))
  }
  out
}
raw <- load_dir("results")
spl <- load_dir("results-spell")
build <- vapply(c("results", "results-spell"), function(d) {
  f <- file.path(PORT_OUT, d, paste0(VIGS[1], ".log"))
  if (file.exists(f)) sub("^# ", "", readLines(f, n = 1)) else "?"
}, "")
cat("## provenance\n")
cat("raw:  ", build[1], "\nspell:", build[2], "\n")
cat("vignettes with a raw record:", sum(file.exists(file.path(
  PORT_OUT, "results", paste0(VIGS, ".rds")))), "of", length(VIGS), "\n\n")

as_df <- function(res) {
  do.call(rbind, lapply(res, function(r) data.frame(
    id = r$id, vignette = r$vignette, kind = r$kind, status = r$status,
    patch = paste(r$patch, collapse = "+"),
    msg = substr(gsub("\n", " ", r$msg %||% ""), 1, 220),
    src = gsub("\n *", " ", r$src),
    src_run = gsub("\n *", " ", r$src_run %||% r$src),
    secs = r$secs, stringsAsFactors = FALSE
  )))
}
R <- as_df(raw); S <- as_df(spl)
M <- merge(R[, c("id", "vignette", "kind", "status", "msg", "src", "secs")],
           S[, c("id", "status", "patch", "msg", "src_run", "secs")],
           by = "id", suffixes = c("_raw", "_spell"), all = TRUE)
M <- M[order(match(M$vignette, VIGS),
             as.integer(sub(".*\\.(\\d+)\\.\\d+$", "\\1", M$id)),
             as.integer(sub(".*\\.(\\d+)$", "\\1", M$id))), ]

M$class <- vapply(seq_len(nrow(M)), function(i) classify(M[i, ]), character(1))
M$bucket <- sub(":.*$", "", M$class)

tally <- function(d, k) {
  x <- d[d$kind == k & d$status_raw != "SETUP-SKIP", ]
  b <- factor(x$bucket, levels = c("CLEAN", "SPELLING", "FAIL", "CASCADE"))
  c(total = nrow(x), table(b))
}
cat("## headline\n")
for (k in c("model", "post", "other")) {
  z <- tally(M, k)
  cat(sprintf(paste0("%-6s total=%2d  CLEAN=%2d  SPELLING=%2d  FAIL=%2d",
                     "  CASCADE=%2d\n"),
              k, z[1], z[2], z[3], z[4], z[5]))
}
cat("\n## spelling changes by patch (model calls)\n")
print(table(M$class[M$kind == "model" & M$bucket == "SPELLING"]))
cat("\n## spelling changes by patch (post calls)\n")
print(table(M$class[M$kind == "post" & M$bucket == "SPELLING"]))

# A patch is measured, not assumed: an explicit PATCH entry is applied
# in the spell pass whether or not the raw line fails, so the raw pass
# alone says whether the porter still needs it. An AUTO_RETRY rule
# fires only on the error it was written for, so a rule that never
# fires is a gap that closed.
cat("\n## patch necessity (explicit per-expression patches)\n")
for (id in names(PATCH)) {
  r <- M[M$id == id, ]
  if (!nrow(r)) { cat(sprintf("%-26s no such expression\n", id)); next }
  v <- if (r$status_raw == "OK") "UNNEEDED: the raw line runs" else
    if (r$status_spell == "OK") "NEEDED" else "INSUFFICIENT: fails patched"
  cat(sprintf("%-26s raw=%-7s spell=%-7s %s\n", id, r$status_raw,
              r$status_spell, v))
  if (r$status_raw != "OK") cat("      raw! ", substr(r$msg_raw, 1, 160), "\n")
}
need <- load_dir("results-need")
if (length(need)) {
  N <- as_df(need)
  cat("\n## patch necessity, need pass (only the patches marked * applied)\n")
  for (id in intersect(names(PATCH), N$id)) {
    r <- N[N$id == id, ]
    cat(sprintf("%-26s %s need=%-7s %s\n", id,
                if (nzchar(r$patch)) "*" else " ", r$status,
                if (nzchar(r$patch)) "" else if (r$status == "OK")
                  "UNNEEDED once the patches upstream of it are applied"
                else "NEEDED"))
  }
}
cat("\n## patch necessity (generic retry rules)\n")
for (pn in names(AUTO_RETRY)) {
  hit <- grepl(pn, M$patch, fixed = TRUE)
  cat(sprintf("%-20s fired on %d expressions, %d of them then ran\n", pn,
              sum(hit), sum(hit & M$status_spell == "OK")))
}

cat("\n## per vignette",
    "(model calls: total | clean | spelling | fail | cascade)\n")
for (v in VIGS) {
  z <- tally(M[M$vignette == v, ], "model")
  zp <- tally(M[M$vignette == v, ], "post")
  cat(sprintf(paste0("%-20s model %2d | %2d %2d %2d %2d",
                     "    post %3d | %2d %2d %2d %2d\n"),
              v, z[1], z[2], z[3], z[4], z[5],
              zp[1], zp[2], zp[3], zp[4], zp[5]))
}

if (file.exists(BASE)) {
  B <- readRDS(BASE)
  cat("\n## moves against the baseline", basename(dirname(BASE)), "\n")
  key <- intersect(M$id, B$id)
  cat("expressions in both:", length(key), " only now:",
      length(setdiff(M$id, B$id)), " only then:",
      length(setdiff(B$id, M$id)), "\n")
  src_same <- B$src[match(key, B$id)] == M$src[match(key, M$id)]
  cat("of those, same transformed source:", sum(src_same), "\n")
  for (k in c("model", "post")) {
    kk <- key[M$kind[match(key, M$id)] == k]
    then <- factor(B$bucket[match(kk, B$id)],
                   levels = c("CLEAN", "SPELLING", "FAIL", "CASCADE"))
    now <- factor(M$bucket[match(kk, M$id)],
                  levels = c("CLEAN", "SPELLING", "FAIL", "CASCADE"))
    cat("\n", k, ": rows = 0.34.0 bucket, columns = now\n", sep = "")
    print(table(then, now))
  }
  # Rows never executed (BRMS-ONLY and the like) have no rank, so they
  # can be neither better nor worse.
  lv <- c("CLEAN", "SPELLING", "CASCADE", "FAIL")
  r_now <- match(M$bucket[match(key, M$id)], lv)
  r_then <- match(B$bucket[match(key, B$id)], lv)
  worse <- key[!is.na(r_now) & !is.na(r_then) & r_now > r_then]
  cat("\n## worse than the baseline (regression candidates)\n")
  if (!length(worse)) cat("none\n")
  for (id in worse) {
    cat(sprintf("%-28s %s -> %s  %s\n", id, B$class[B$id == id],
                M$class[M$id == id], substr(M$src[M$id == id], 1, 70)))
    cat("      now!  ", substr(M$msg_spell[M$id == id], 1, 200), "\n")
  }
}

cat("\n## every model and post expression\n")
sub <- M[M$kind %in% c("model", "post"), ]
for (i in seq_len(nrow(sub))) {
  r <- sub[i, ]
  cat(sprintf("%-28s %-5s %-10s %-9s %s\n", r$id, r$kind, r$status_raw,
              r$status_spell, substr(r$src, 1, 80)))
  if (nzchar(r$msg_raw)) cat("      raw!  ", substr(r$msg_raw, 1, 200), "\n")
  if (nzchar(r$msg_spell) && !identical(r$msg_spell, r$msg_raw))
    cat("      spell!", substr(r$msg_spell, 1, 200), "\n")
  if (nzchar(r$patch)) cat("      patch:", r$patch, "\n")
}
cat("\n## distinct failure messages, by kind (spell pass, FAIL only)\n")
f <- M[M$bucket == "FAIL", ]
if (nrow(f)) print(table(substr(f$msg_spell, 1, 52), f$kind)) else cat("none\n")
saveRDS(M, file.path(PORT_OUT, "results-merged.rds"))
utils::write.csv(M, file.path(PORT_OUT, "results-merged.csv"),
                 row.names = FALSE)
