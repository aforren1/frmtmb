# Reduce the runs to the scorecard tables. Rscript summarize.R
HERE <- "C:/Users/adf44/source/r/frmtmb-wt-audit/dev/brms-port"
source(file.path(HERE, "port-lib.R"))
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")

load_dir <- function(d) {
  out <- list()
  for (v in VIGS) {
    f <- file.path(HERE, d, paste0(v, ".rds"))
    if (file.exists(f)) out <- c(out, readRDS(f))
  }
  out
}
raw <- load_dir("results")
spl <- load_dir("results-spell")

as_df <- function(res, tag) {
  do.call(rbind, lapply(res, function(r) data.frame(
    id = r$id, vignette = r$vignette, kind = r$kind, status = r$status,
    patch = paste(r$patch, collapse = "+"),
    msg = substr(gsub("\n", " ", r$msg %||% ""), 1, 220),
    src = gsub("\n *", " ", r$src),
    src_run = gsub("\n *", " ", r$src_run %||% r$src),
    secs = r$secs, stringsAsFactors = FALSE
  )))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
R <- as_df(raw); S <- as_df(spl)
M <- merge(R[, c("id", "vignette", "kind", "status", "msg", "src", "secs")],
           S[, c("id", "status", "patch", "msg", "src_run", "secs")],
           by = "id", suffixes = c("_raw", "_spell"), all = TRUE)
M <- M[order(match(M$vignette, VIGS),
             as.integer(sub(".*\\.(\\d+)\\.\\d+$", "\\1", M$id)),
             as.integer(sub(".*\\.(\\d+)$", "\\1", M$id))), ]

M$class <- vapply(seq_len(nrow(M)), function(i) classify(M[i, ]), character(1))
M$bucket <- sub(":.*$", "", M$class)

# The projection pass is optional; count it only if it has been run.
v035 <- load_dir("results-v035")
if (length(v035)) {
  V <- as_df(v035)
  M$status_v035 <- V$status[match(M$id, V$id)]
}

tally <- function(d, k) {
  x <- d[d$kind == k & d$status_raw != "SETUP-SKIP", ]
  b <- factor(x$bucket, levels = c("CLEAN", "SPELLING", "FAIL", "CASCADE"))
  c(total = nrow(x), table(b))
}
cat("## headline\n")
for (k in c("model", "post", "other")) {
  z <- tally(M, k)
  cat(sprintf("%-6s total=%2d  CLEAN=%2d  SPELLING=%2d  FAIL=%2d  CASCADE=%2d\n",
              k, z[1], z[2], z[3], z[4], z[5]))
}
cat("\n## spelling changes by patch (model calls)\n")
print(table(M$class[M$kind == "model" & M$bucket == "SPELLING"]))
if (!is.null(M$status_v035)) {
  cat("\n## projection: only FN-1 (default family) and FN-10 (lf) fixed\n")
  for (k in c("model", "post")) {
    x <- M[M$kind == k & M$status_raw != "SETUP-SKIP", ]
    cat(sprintf("%-6s clean=%d of %d   (still failing: %d)\n", k,
                sum(x$status_v035 == "OK", na.rm = TRUE), nrow(x),
                sum(x$status_v035 != "OK", na.rm = TRUE)))
  }
}
cat("\n## per vignette (model calls: total | clean | spelling | fail | cascade)\n")
for (v in VIGS) {
  z <- tally(M[M$vignette == v, ], "model")
  zp <- tally(M[M$vignette == v, ], "post")
  cat(sprintf("%-20s model %2d | %2d %2d %2d %2d    post %3d | %2d %2d %2d %2d\n",
              v, z[1], z[2], z[3], z[4], z[5], zp[1], zp[2], zp[3], zp[4], zp[5]))
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
print(table(substr(f$msg_spell, 1, 52), f$kind))
saveRDS(M, file.path(HERE, "results-merged.rds"))
write.csv(M, file.path(HERE, "results-merged.csv"), row.names = FALSE)
