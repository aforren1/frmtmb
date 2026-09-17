# Lane wt-conditions, punch round 1: append the "Punch round 1" section
# to dev/conditions-findings.md from its template, pasting the logs into
# the generated blocks rather than typing any number. Run once, after
# the checks.
#   Rscript dev/conditions-punch-findings.R
lg <- function(f) readLines(file.path("dev/conditions-punch-log", f),
                            warn = FALSE)
probe <- c(paste("before:", tail(lg("probe-before.txt"), 1L)),
           paste("after: ", tail(lg("probe-after.txt"), 1L)))
census <- lg("census.txt")
census <- census[nzchar(census) & !grepl("^Warning|dir.create|already",
                                         census)]
before <- lg("summary-before.txt")
after <- lg("summary-after.txt")
block <- function(x, from, to) {
  i <- grep(from, x)[1L]
  j <- i + grep(to, x[i:length(x)])[1L] - 1L
  x[i:j]
}
sweep <- c("before (summary-before.txt):",
           block(before, "^== extension subclass", "^ALL "),
           "", "after (summary-after.txt):",
           block(after, "^== suites", "^ALL "))
checks <- unlist(lapply(c("check-core.txt", "check-sample.txt"), function(f) {
  x <- lg(f)
  c(paste0(f, ":"),
    grep("^[*] checking .*(NOTE|WARNING|ERROR)$|^Status:|^  Skipping|^Ran",
         x, value = TRUE), "")
}))
t <- readLines("dev/conditions-punch-findings-template.md")
fill <- function(t, key, val) {
  i <- which(t == key)
  stopifnot(length(i) == 1L)
  c(t[seq_len(i - 1L)], val, t[(i + 1L):length(t)])
}
t <- fill(t, "@PROBE@", probe)
t <- fill(t, "@CENSUS@", census)
t <- fill(t, "@SWEEP@", sweep)
t <- fill(t, "@CHECKS@", checks)
f <- "dev/conditions-findings.md"
old <- readLines(f)
stopifnot(!any(old == "## Punch round 1"))
writeLines(c(old, t), f)
