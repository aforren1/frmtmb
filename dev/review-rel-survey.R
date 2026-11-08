# Release review 0.80.0: the survey. Which of the 172 sampled pairs are
# among the core pairs whose compat status moved between 0.79.0's frame
# and this build's, and what each such pair's recorded verdict and status
# are; then 5 pairs that touch bigx, osa, readers or mfxse features,
# side by side with the 0.79.0 record, every column.
#
#   Rscript dev/review-rel-survey.R
rd <- function(f) utils::read.delim(f, colClasses = "character",
                                    na.strings = "NA")
Sys.setenv(GIT_OPTIONAL_LOCKS = "0")
old <- tempfile(fileext = ".tsv")
writeLines(system2("git", c("-c", "core.autocrlf=false", "show",
                            "3cd2a0c8:dev/equiv-data/results.tsv"),
                   stdout = TRUE), old)
a <- rd(old)
b <- rd("dev/equiv-data/results.tsv")
cat("columns:", paste(names(b), collapse = " "), "\n")
fr <- rd("dev/equiv-frame-out/compat.tsv")
r19 <- rd("dev/rel080-log/frame-r19/compat.tsv")
m2 <- merge(fr[fr$in_core_only == "TRUE", c("feature_a", "feature_b", "status")],
            r19[r19$in_core_only == "TRUE", c("feature_a", "feature_b", "status")],
            by = c("feature_a", "feature_b"))
mv <- m2[m2$status.x != m2$status.y, ]
cat("moved core pairs:", nrow(mv), "\n")
print(mv, row.names = FALSE)
key <- function(d) paste(d$feature_a, d$feature_b, sep = " x ")
key2 <- function(d) paste(d$feature_b, d$feature_a, sep = " x ")
hit <- b[key(b) %in% key(mv) | key2(b) %in% key(mv), ]
cat("\nsampled pairs among them:", nrow(hit), "\n")
if (nrow(hit)) print(hit[, intersect(c("id", "feature_a", "feature_b",
  "status", "verdict", "level", "agree", "reason"), names(b))],
  row.names = FALSE)
# 5 pairs whose features the round's fit and reader changes reach
pat <- "osa|resid|new|level|fitted|predict|marginal|me\\(|offset|poly|ce|conditional"
cand <- b[grepl(pat, paste(b$feature_a, b$feature_b), ignore.case = TRUE), ]
cat("\npairs naming those features:", nrow(cand), "\n")
print(cand[, c("id", "feature_a", "feature_b", "verdict")], row.names = FALSE)
pick <- utils::head(cand$id, 5L)
for (i in pick) {
  cat("\n==", i, "\n")
  x <- rbind(old = a[a$id == i, ], new = b[b$id == i, ])
  same <- vapply(names(x), function(n) identical(x$old[[n]], x$new[[n]]) ||
                   identical(x[1, n], x[2, n]), NA)
  for (n in names(x)) {
    cat(sprintf("  %-14s %s %s | %s\n", n, if (identical(x[1, n], x[2, n]))
      "=" else "*", substr(x[1, n], 1, 70), substr(x[2, n], 1, 70)))
  }
}
cat("SURVEY CHECK DONE\n")
