# Emit the generated tables of dev/remeasure-findings.md as markdown
# files under dev/remeasure-out/blocks/, from the records of two runs of
# dev/remeasure-run.sh (the "before" and the "after" build) and the
# tracked 0.34.0 baseline. dev/remeasure-fill.R pastes them.
#
#   Rscript dev/remeasure-blocks.R [before after]   (default v067 r8)
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
a <- commandArgs(trailingOnly = TRUE)
tb <- if (length(a) >= 2) a[1] else "v067"
ta <- if (length(a) >= 2) a[2] else "r8"
od <- file.path(root, "remeasure-out", "blocks")
dir.create(od, showWarnings = FALSE, recursive = TRUE)
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")
pm <- function(t) readRDS(file.path(root, "remeasure-out", t, "port",
                                    "results-merged.rds"))
ver <- function(t) sub("^# frmtmb ([^ ]+) from.*$", "\\1", readLines(
  file.path(root, "remeasure-out", t, "port", "results", "brms_overview.log"),
  n = 1))
B <- readRDS(file.path(root, "brms-port", "results-merged.rds"))
M0 <- pm(tb); M1 <- pm(ta)
v0 <- ver(tb); v1 <- ver(ta)
lv <- c("CLEAN", "SPELLING", "FAIL", "CASCADE")
tal <- function(d, k) {
  x <- d[d$kind == k & d$status_raw != "SETUP-SKIP", ]
  c(nrow(x), table(factor(x$bucket, levels = lv)))
}
pct <- function(a, n) sprintf("%d%%", round(100 * a / n))
w <- function(name, lines) writeLines(lines, file.path(od, paste0(name,
                                                                  ".md")))

x <- c(sprintf("| calls | outcome | 0.34.0 | %s | %s |", v0, v1),
       "|---|---|---|---|---|")
for (k in c("model", "post")) {
  a0 <- tal(B, k); b0 <- tal(M0, k); b1 <- tal(M1, k)
  for (j in 2:5) {
    x <- c(x, sprintf("| %s (%d) | %s | %d (%s) | %d (%s) | %d (%s) |", k,
                      b1[1], lv[j - 1], a0[j], pct(a0[j], a0[1]), b0[j],
                      pct(b0[j], b0[1]), b1[j], pct(b1[j], b1[1])))
  }
}
w("headline", x)

cell <- function(z) sprintf("%d \\| %d %d %d %d", z[1], z[2], z[3], z[4],
                            z[5])
x <- c(sprintf("| vignette | model %s | model %s | post %s | post %s |", v0,
               v1, v0, v1), "|---|---|---|---|---|")
for (v in VIGS) {
  x <- c(x, sprintf("| %s | %s | %s | %s | %s |", v,
                    cell(tal(M0[M0$vignette == v, ], "model")),
                    cell(tal(M1[M1$vignette == v, ], "model")),
                    cell(tal(M0[M0$vignette == v, ], "post")),
                    cell(tal(M1[M1$vignette == v, ], "post"))))
}
w("per-vignette", x)

# spell-pass OK counts, the second headline brms port users read
sp <- function(d, k) sum(d$status_spell[d$kind == k] == "OK")
x <- c(sprintf("| pass | calls | %s | %s |", v0, v1), "|---|---|---|---|")
for (k in c("model", "post", "other")) {
  x <- c(x, sprintf("| raw OK | %s (%d) | %d | %d |", k,
                    sum(M1$kind == k), sum(M0$status_raw[M0$kind == k] == "OK"),
                    sum(M1$status_raw[M1$kind == k] == "OK")),
         sprintf("| spell OK | %s (%d) | %d | %d |", k, sum(M1$kind == k),
                 sp(M0, k), sp(M1, k)))
}
# the keep-prior transform (run-vignette.R mode keepprior)
kp <- function(t) {
  fs <- list.files(file.path(root, "remeasure-out", t, "port",
                             "results-keepprior"), "[.]rds$",
                   full.names = TRUE)
  r <- do.call(c, lapply(fs, readRDS))
  data.frame(id = vapply(r, `[[`, "", "id"),
             kind = vapply(r, `[[`, "", "kind"),
             status = vapply(r, `[[`, "", "status"))
}
K0 <- kp(tb); K1 <- kp(ta)
for (k in c("model", "post")) {
  x <- c(x, sprintf("| keep-prior OK | %s (%d) | %d | %d |", k,
                    sum(K1$kind == k), sum(K0$status[K0$kind == k] == "OK"),
                    sum(K1$status[K1$kind == k] == "OK")))
}
w("passes", x)
kd <- merge(K1, M1[c("id", "status_raw")], by = "id")
kd <- kd[kd$kind %in% c("model", "post") &
           (kd$status == "OK") != (kd$status_raw == "OK"), ]
x <- sprintf("- %s (%s): raw %s, keep-prior %s", kd$id, kd$kind,
             kd$status_raw, kd$status)
w("keepprior-diff", c(sprintf("rows of %s whose outcome differs between the",
                              v1), "raw and the keep-prior transform:", "",
                      x))

# the hand translation's drift at the after build
D <- utils::read.csv(file.path(root, "remeasure-out", ta, "bv", "drift.csv"),
                     stringsAsFactors = FALSE)
D$path <- ifelse(grepl("^SAMPLE:", D$label), "SAMPLE", "ML")
lab <- c("CLEAN", "SPELLING", "BEHAVIOR", "MISSING", "REFUSAL")
x <- c(sprintf("| label at 0.42.0 | rows | ran at %s | erred at %s |", v1, v1),
       "|---|---|---|---|")
for (l in lab) {
  d <- D[D$edge == l, ]
  x <- c(x, sprintf("| %s | %d | %d | %d |", l, nrow(d), sum(d$ok),
                    sum(!d$ok)))
}
x <- c(x, sprintf("| all | %d | %d | %d |", nrow(D), sum(D$ok), sum(!D$ok)),
       "", "| drift | ML | SAMPLE | all |", "|---|---|---|---|")
for (k in c("AS-LABELED", "FAILS-NOW", "RUNS-NOW", "STILL-OPEN")) {
  x <- c(x, sprintf("| %s | %d | %d | %d |", k,
                    sum(D$drift == k & D$path == "ML"),
                    sum(D$drift == k & D$path == "SAMPLE"),
                    sum(D$drift == k)))
}
w("hand-drift", x)

# the bin-2P survey (dev/brmsport-bin2-survey.R)
gf <- file.path(root, "brmsport-bin2-groups.tsv")
if (file.exists(gf)) {
  G <- utils::read.delim(gf, quote = "", stringsAsFactors = FALSE)
  x <- c(paste("| group | blocks | assertions | ported | reads |",
               "frmtmb counterpart | cost | class if unported |"),
         "|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(G))) {
    x <- c(x, sprintf(paste("| %s | %d | %d | %d | %s | %s | %s |",
                            "cannot transfer: %s |"),
                      G$group[i], G$blocks[i], G$assertions[i], G$ported[i],
                      G$reads[i], G$frmtmb_counterpart[i], G$cost[i],
                      G$class[i]))
  }
  x <- c(x, sprintf("| **total** | **%d** | **%d** | **%d** | | | | |",
                    sum(G$blocks), sum(G$assertions), sum(G$ported)))
  w("bin2-groups", x)
}
cat("wrote", length(list.files(od)), "blocks to", od, "\n")
