# Reduce the per-vignette CSV records into the audit scoreboard.
#
# Run after every dev/brms-vignettes/<name>.R has been run with BV_OUT
# pointing at one directory:
#   Rscript dev/brms-vignettes/_summarize.R
out <- Sys.getenv("BV_OUT", unset = tempdir())
files <- list.files(out, pattern = "[.]csv$", full.names = TRUE)
d <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
d$edge[is.na(d$edge) | d$edge == ""] <- "CLEAN"
d <- d[d$kind != "data", ]

# every label carries its inference path as a prefix, "ML: " or
# "SAMPLE: ", so the two paths are scored apart without a second column
# in every bv() call
d$path <- ifelse(grepl("^SAMPLE:", d$label), "SAMPLE",
                 ifelse(grepl("^ML:", d$label), "ML", "ML"))

cls <- c("CLEAN", "SPELLING", "BEHAVIOR", "MISSING", "REFUSAL")

score <- function(s) c(nrow(s), vapply(cls, function(cc) sum(s$edge == cc), 1L))

for (pth in c("ML", "SAMPLE")) {
  cat("## Per vignette, ", pth, " path\n\n", sep = "")
  cat("| vignette | model | clean | spell | behav | miss | refuse |",
      "post | clean | spell | behav | miss | refuse |\n")
  cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
  dp <- d[d$path == pth, ]
  for (v in sort(unique(d$vignette))) {
    s <- dp[dp$vignette == v, ]
    row <- c(v, score(s[s$kind == "model", ]), score(s[s$kind == "post", ]))
    cat("| ", paste(row, collapse = " | "), " |\n", sep = "")
  }
  cat("\n")
}

cat("\n## Totals\n\n")
print(table(d$path, d$kind, d$edge))
for (pth in c("ML", "SAMPLE")) {
  dp <- d[d$path == pth, ]
  cat("\n", pth, " model calls:", sum(dp$kind == "model"),
      " clean:", sum(dp$kind == "model" & dp$edge == "CLEAN"),
      " ran ok:", sum(dp$kind == "model" & dp$ok), "\n")
  cat(pth, " post calls:", sum(dp$kind == "post"),
      " clean:", sum(dp$kind == "post" & dp$edge == "CLEAN"),
      " ran ok:", sum(dp$kind == "post" & dp$ok), "\n")
}
cat("\nall calls:", nrow(d), " clean:", sum(d$edge == "CLEAN"),
    " ran ok:", sum(d$ok), "\n")
cat("total seconds:", round(sum(d$secs), 1), "\n")

cat("\n## Every non-clean edge\n\n")
e <- d[d$edge != "CLEAN",
       c("vignette", "path", "kind", "edge", "label", "why", "ok", "msg")]
e <- e[order(e$edge, e$path, e$vignette), ]
for (i in seq_len(nrow(e))) {
  cat("- [", e$edge[i], "|", e$path[i], "] ", e$vignette[i], " / ",
      e$label[i], "\n", "    ", e$why[i], "\n", sep = "")
  if (!e$ok[i]) cat("    error: ", e$msg[i], "\n", sep = "")
}

cat("\n## Failures with no edge label (translation bugs, should be zero)\n\n")
b <- d[!d$ok & d$edge == "CLEAN", c("vignette", "label", "msg")]
if (!nrow(b)) cat("none\n") else print(b)
