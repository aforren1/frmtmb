# Emit the generated blocks of the 2026-10-05 re-measurement, as
# markdown, for dev/brms-vignette-port.md, dev/brms-vignette-audit.md
# and dev/vigport-findings.md. Every count in those sections is pasted
# from this script's output, never typed.
#
#   Rscript dev/vigport-blocks.R > dev/vigport-blocks.txt
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")
M <- readRDS(file.path(root, "vigport-port-out/r5/results-merged.rds"))
B <- readRDS(file.path(root, "brms-port/results-merged.rds"))
lv <- c("CLEAN", "SPELLING", "FAIL", "CASCADE")
tal <- function(d, k) {
  x <- d[d$kind == k & d$status_raw != "SETUP-SKIP", ]
  c(nrow(x), table(factor(x$bucket, levels = lv)))
}
pct <- function(a, n) sprintf("%d%%", round(100 * a / n))

cat("<!-- BEGIN generated: mechanical headline -->\n")
cat("| calls | outcome | 0.34.0 | 0.67.0 |\n|---|---|---|---|\n")
for (k in c("model", "post")) {
  a <- tal(B, k); b <- tal(M, k)
  for (j in 2:5) {
    cat(sprintf("| %s (%d) | %s | %d (%s) | %d (%s) |\n", k, b[1], lv[j - 1],
                a[j], pct(a[j], a[1]), b[j], pct(b[j], b[1])))
  }
}
cat("<!-- END generated -->\n\n")

cat("<!-- BEGIN generated: mechanical per vignette -->\n")
cat("| vignette | model 0.34.0 | model 0.67.0 | post 0.34.0 | post 0.67.0 |\n")
cat("|---|---|---|---|---|\n")
cell <- function(z) sprintf("%d \\| %d %d %d %d", z[1], z[2], z[3], z[4], z[5])
for (v in VIGS) {
  cat(sprintf("| %s | %s | %s | %s | %s |\n", v,
              cell(tal(B[B$vignette == v, ], "model")),
              cell(tal(M[M$vignette == v, ], "model")),
              cell(tal(B[B$vignette == v, ], "post")),
              cell(tal(M[M$vignette == v, ], "post"))))
}
cat("<!-- END generated -->\n\n")

cat("<!-- BEGIN generated: mechanical moves -->\n")
for (k in c("model", "post")) {
  key <- M$id[M$kind == k & M$status_raw != "SETUP-SKIP"]
  then <- factor(B$bucket[match(key, B$id)], levels = lv)
  now <- factor(M$bucket[match(key, M$id)], levels = lv)
  t <- table(then, now)
  cat("\n", k,
      " calls: rows are the 0.34.0 class, columns the 0.67.0 class\n\n",
      sep = "")
  cat("| 0.34.0 \\\\ 0.67.0 |", paste(lv, collapse = " | "), "|\n")
  cat("|---|---|---|---|---|\n")
  for (r in lv) cat("|", r, "|", paste(t[r, ], collapse = " | "), "|\n")
}
cat("<!-- END generated -->\n\n")

# The spell-pass edits that are still needed, by the need pass.
cat("<!-- BEGIN generated: spelling rows at 0.67.0 -->\n")
s <- M[M$bucket == "SPELLING", ]
cat("| id | kind | class | raw error |\n|---|---|---|---|\n")
for (i in seq_len(nrow(s))) {
  cat(sprintf("| %s | %s | %s | %s |\n", s$id[i], s$kind[i], s$class[i],
              gsub("\\|", "/", substr(s$msg_raw[i], 1, 70))))
}
cat("<!-- END generated -->\n\n")

cat("<!-- BEGIN generated: hand translation drift -->\n")
D <- utils::read.csv(file.path(root, "vigport-bv-out/r5/drift.csv"),
                     stringsAsFactors = FALSE)
D$path <- ifelse(grepl("^SAMPLE:", D$label), "SAMPLE", "ML")
cat("| label at 0.42.0 | rows | ran at 0.67.0 | erred at 0.67.0 |\n")
cat("|---|---|---|---|\n")
for (e in c("CLEAN", "SPELLING", "BEHAVIOR", "MISSING", "REFUSAL")) {
  x <- D[D$edge == e, ]
  cat(sprintf("| %s | %d | %d | %d |\n", e, nrow(x), sum(x$ok), sum(!x$ok)))
}
cat(sprintf("| all | %d | %d | %d |\n", nrow(D), sum(D$ok), sum(!D$ok)))
cat("\n| drift | ML | SAMPLE | all |\n|---|---|---|---|\n")
for (k in c("AS-LABELED", "FAILS-NOW", "RUNS-NOW", "STILL-OPEN")) {
  cat(sprintf("| %s | %d | %d | %d |\n", k, sum(D$drift == k & D$path == "ML"),
              sum(D$drift == k & D$path == "SAMPLE"), sum(D$drift == k)))
}
cat("\nPer vignette, model then post calls that RAN at 0.67.0, of all:\n\n")
cat("| vignette | ML model | ML post | SAMPLE model | SAMPLE post |\n")
cat("|---|---|---|---|---|\n")
for (v in sort(unique(D$vignette))) {
  f <- function(p, k) {
    x <- D[D$vignette == v & D$path == p & D$kind == k, ]
    if (!nrow(x)) "0" else sprintf("%d of %d", sum(x$ok), nrow(x))
  }
  cat(sprintf("| %s | %s | %s | %s | %s |\n", v, f("ML", "model"),
              f("ML", "post"), f("SAMPLE", "model"), f("SAMPLE", "post")))
}
cat("<!-- END generated -->\n\n")

cat("<!-- BEGIN generated: workaround necessity -->\n")
W <- utils::read.csv(file.path(root, "vigport-bv-out/r5-wa/workarounds.csv"),
                     stringsAsFactors = FALSE)
W <- W[W$kind != "data", ]
cat("| row (the brms spelling, run on 0.67.0) | runs |",
    "the 0.42.0 workaround | error now |\n")
cat("|---|---|---|---|\n")
for (i in seq_len(nrow(W))) {
  cat(sprintf("| %s | %s | %s | %s |\n", gsub("\\|", "/", W$label[i]),
              if (W$ok[i]) "yes" else "NO", gsub("\\|", "/", W$why[i]),
              gsub("\\|", "/", substr(gsub("\n", " ", W$msg[i]), 1, 90))))
}
w <- W[grepl("^WA ", W$label), ]
cat(sprintf(paste0("\nworkaround rows: %d; the brms spelling now runs on %d,",
                   " still fails on %d\n"),
            nrow(w), sum(w$ok), sum(!w$ok)))
cat("<!-- END generated -->\n")

cat("\n<!-- BEGIN generated: keep-prior transform -->\n")
kp <- list()
for (v in VIGS) {
  f <- file.path(root, "vigport-port-out/r5/results-keepprior",
                 paste0(v, ".rds"))
  kp <- c(kp, readRDS(f))
}
K <- data.frame(id = vapply(kp, `[[`, "", "id"),
                kind = vapply(kp, `[[`, "", "kind"),
                status = vapply(kp, `[[`, "", "status"),
                msg = vapply(kp, function(r) r$msg %||% "", ""),
                stringsAsFactors = FALSE)
for (k in c("model", "post")) {
  x <- K[K$kind == k & K$status != "SETUP-SKIP", ]
  r <- M[M$kind == k & M$status_raw != "SETUP-SKIP", ]
  cat(sprintf(paste0("%s calls: raw transform %d of %d run;",
                     " keep-prior transform %d of %d run\n"),
              k, sum(r$status_raw == "OK"), nrow(r), sum(x$status == "OK"),
              nrow(x)))
}
d <- merge(K, M[, c("id", "status_raw")], by = "id")
d <- d[d$kind %in% c("model", "post") & d$status != d$status_raw, ]
cat("\nrows whose outcome differs between the two transforms:\n\n")
for (i in seq_len(nrow(d))) {
  cat(sprintf("- %s (%s): raw %s, keep-prior %s%s\n", d$id[i], d$kind[i],
              d$status_raw[i], d$status[i],
              if (d$status[i] != "OK") paste0(": ", substr(d$msg[i], 1, 110))
              else ""))
}
cat("<!-- END generated -->\n")

# Blocks cut from the other reducers' outputs, so the records quote them
# verbatim. Each source file is regenerated by its own script first.
# `to` is the first line after the block; NULL runs to the end of file.
# fence = TRUE wraps plain-text output in a text fence; otherwise a
# "## " heading of the source becomes bold, so it does not break the
# heading structure of the record it lands in.
cut_block <- function(name, file, from, to = NULL, fence = FALSE) {
  x <- readLines(file.path(root, file), warn = FALSE)
  a <- grep(from, x)[1]
  if (is.na(a)) stop("block ", name, " not found in ", file)
  b <- if (is.null(to)) NA else a + grep(to, x[(a + 1):length(x)])[1]
  if (is.na(b)) b <- length(x) + 1
  body <- x[a:(b - 1)]
  while (length(body) && !nzchar(trimws(body[length(body)]))) {
    body <- body[-length(body)]
  }
  body <- if (fence) c("```text", body, "```") else
    sub("^## (.*)$", "**\\1**", body)
  cat("\n<!-- BEGIN generated: ", name, " -->\n", sep = "")
  cat(body, sep = "\n")
  cat("<!-- END generated -->\n")
}
cut_block("plausibility", "vigport-port-out/r5/plausibility.txt",
          "^## provenance", "^## detail")
cut_block("plausibility detail", "vigport-port-out/r5/plausibility.txt",
          "^## detail", "^## models not compared", fence = TRUE)
cut_block("plausibility not compared", "vigport-port-out/r5/plausibility.txt",
          "^## models not compared")
cut_block("gap ranking", "vigport-gaps.txt", "^## mechanical port")
cut_block("history", "vigport-history.txt", "^## mechanical port",
          fence = TRUE)
cut_block("families", "vigport-port-out/r5/families.txt", "^# frmtmb",
          "^frmtmb registry", fence = TRUE)
cut_block("repros", "vigport-port-out/repros-rellib-r5.txt", "^frmtmb 0",
          fence = TRUE)
# Punch round 1 evidence. sx.txt and sx0.txt are this lane's reruns of
# the reviewer's dev/vigport-rev-sx.R and dev/vigport-rev-sx0.R.
cut_block("pr1 spell change", "vigport-pr1-out/spellcmp.txt", "^rows:",
          fence = TRUE)
cut_block("pr1 sx factors", "vigport-pr1-out/sx.txt", "^## 1[.]",
          "^## 2[.]", fence = TRUE)
cut_block("pr1 sx fit_smooth1", "vigport-pr1-out/sx.txt", "^## 6[.]",
          fence = TRUE)
cut_block("pr1 nhanes", "vigport-pr1-out/nhanes.txt", "^observed",
          fence = TRUE)
cut_block("pr1 convergence", "vigport-pr1-out/conv-allfit.txt", "^frmtmb",
          fence = TRUE)
cut_block("pr1 mo sweep", "vigport-pr1-out/mo-sweep.txt", "^frmtmb",
          fence = TRUE)
cut_block("pr1 noseed", "vigport-pr1-out/noseed-cmp.txt", "^rows",
          fence = TRUE)
cut_block("pr1 sx0", "vigport-pr1-out/sx0.txt", "^linpreds", fence = TRUE)
