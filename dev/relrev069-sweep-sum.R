# Reviewer of the 0.69.0 consolidation: summarize the fuzz-family sweep
# (dev/relrev069-sweep.sh): per variant and BLAS, errors on 0.69.0, on
# 0.68.1 and on 0.69.0 with the candidate fix, and the seeds where
# 0.69.0 errors and 0.68.1 returns a fit.
d <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev069-log/sweep"
rd <- function(f) {
  x <- utils::read.delim(f, colClasses = c(seed = "integer",
                                          status = "character"),
                         quote = "")
  x
}
for (v in c("spec", "mo", "ml")) {
  for (b in c("ref", "ob26", "ob32")) {
    fs <- file.path(d, paste0(v, "-", b, "-", c("r7", "r6", "r7fix"),
                              ".tsv"))
    if (!all(file.exists(fs))) next
    r7 <- rd(fs[1]); r6 <- rd(fs[2]); fx <- rd(fs[3])
    s <- Reduce(intersect, list(r7$seed, r6$seed, fx$seed))
    r7 <- r7[match(s, r7$seed), ]; r6 <- r6[match(s, r6$seed), ]
    fx <- fx[match(s, fx$seed), ]
    reg <- s[r7$status == "error" & r6$status == "ok"]
    imp <- s[r7$status == "ok" & r6$status == "error"]
    cat(sprintf(paste0("%-4s %-4s seeds %3d nodata %d | errors: 0.69.0 %2d,",
                       " 0.68.1 %2d, fix %2d | 0.69.0 errors where 0.68.1",
                       " fits: %d (%s) | 0.68.1 errors where 0.69.0 fits: %d",
                       " | fix errors where 0.68.1 fits: %d\n"),
                v, b, length(s), sum(r7$status == "nodata"),
                sum(r7$status == "error"), sum(r6$status == "error"),
                sum(fx$status == "error"), length(reg),
                paste(head(reg, 6), collapse = " "), length(imp),
                sum(fx$status == "error" & r6$status == "ok")))
    em <- unique(substr(r7$detail[r7$status == "error"], 1, 90))
    if (length(em)) cat("     0.69.0 messages:", paste(em, collapse = " || "), "\n")
    cat(sprintf(paste0("     code 0 fits: 0.69.0 %d, 0.68.1 %d, fix %d;",
                       " median objective 0.69.0 - 0.68.1 on both-ok: %.3g\n"),
                sum(r7$code == 0, na.rm = TRUE),
                sum(r6$code == 0, na.rm = TRUE),
                sum(fx$code == 0, na.rm = TRUE),
                stats::median(as.numeric(r7$objective) -
                                as.numeric(r6$objective), na.rm = TRUE)))
  }
}
