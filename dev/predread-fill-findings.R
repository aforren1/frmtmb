# Paste the generated blocks into dev/predread-findings.md from the logs, so
# that no count is typed. Rerun after any new run.
#   Rscript dev/predread-fill-findings.R
f <- "dev/predread-findings.md"
s <- readLines(f, warn = FALSE)
put <- function(s, tag, body) {
  a <- grep(paste0("<!-- predread-", tag, " begin -->"), s, fixed = TRUE)
  b <- grep(paste0("<!-- predread-", tag, " end -->"), s, fixed = TRUE)
  stopifnot(length(a) == 1L, length(b) == 1L, b > a)
  c(s[seq_len(a)], body, s[b:length(s)])
}
seen <- function(tier) {
  d <- file.path("dev/release", paste0(tier, "-files"))
  out <- character(0)
  for (g in list.files(d, full.names = TRUE)) {
    x <- readLines(g, warn = FALSE)
    out <- c(out, grep("^RESULT", x, value = TRUE))
  }
  out
}
seen_body <- c("```", sub("^RESULT ", "", c(seen("predread-seen"),
                                            seen("predread-seen2"))), "```")
s <- put(s, "seen", seen_body)
s <- put(s, "audit", readLines("dev/predread-log/audit-table.md", warn = FALSE))
s <- put(s, "suite", readLines("dev/predread-log/suite-all.md", warn = FALSE))
writeLines(s, f)
