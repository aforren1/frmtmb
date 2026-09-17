# Lane wt-conditions, punch round 2: append the "Punch round 2" section
# to dev/conditions-findings.md, pasting the logs into its blocks.
#   Rscript dev/conditions-punch2-findings.R
lg <- function(f) readLines(file.path("dev/conditions-punch-log", f),
                            warn = FALSE)
probe <- c("before:", lg("p2-probe-before.txt"), "after:",
           lg("p2-probe-after.txt"))
s <- lg("p2-summary.txt")
suites <- s[seq_len(grep("^frmtmb[.]sample logs", s)[1L])]
t <- readLines("dev/conditions-punch2-findings-template.md")
fill <- function(t, key, val) {
  i <- which(t == key)
  stopifnot(length(i) == 1L)
  c(t[seq_len(i - 1L)], val, t[(i + 1L):length(t)])
}
t <- fill(t, "@PROBE@", probe)
t <- fill(t, "@OVERREACH@", lg("p2-census-overreach.txt"))
t <- fill(t, "@CENSUS@", lg("p2-census.txt"))
t <- fill(t, "@SUITES@", suites)
f <- "dev/conditions-findings.md"
old <- readLines(f)
stopifnot(!any(old == "## Punch round 2"))
writeLines(c(old, t), f)
