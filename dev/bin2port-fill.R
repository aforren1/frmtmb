# Fill the generated blocks of dev/bin2port-findings.md from the logs,
# so that no count in the record is typed by hand. A block is the text
# between "<!-- begin NAME -->" and "<!-- end NAME -->"; its content is
# replaced on every run.
#
#   Rscript dev/bin2port-fill.R
f <- "dev/bin2port-findings.md"
x <- readLines(f, warn = FALSE)
rd <- function(p) utils::read.delim(p, quote = "", colClasses = "character",
                                    na.strings = NULL)
strip <- function(p) {
  y <- readLines(p, warn = FALSE)
  y[!grepl("^<!-- generated", y)]
}
logl <- function(p) c("```", readLines(p, warn = FALSE), "```")
tier <- function(dir) {
  fs <- list.files(dir, "^tier-.*[.]txt$", full.names = TRUE)
  c("```", unlist(lapply(fs, function(p) {
    grep("^(lib:|RESULT|FAILING)", readLines(p, warn = FALSE), value = TRUE)
  })), "```")
}
rate709 <- function(p) {
  s <- readLines(p, warn = FALSE)
  m <- regmatches(s, regexpr("pass: [0-9]+ of 823", s))
  np <- as.integer(sub("pass: ([0-9]+) of 823", "\\1", m[1]))
  sprintf(paste("Over the 709 assertions of 0.71.0's denominator (every",
                "assertion but its 114 decisions): %d of 709 pass (%.1f%%)."),
          np, 100 * np / 709)
}
vtab <- function() {
  v <- rd("dev/brmsport-bin2-valuechecks.tsv")
  m <- rd("dev/brmsport-bin2-valuechecks-manual.tsv")
  L <- rd("dev/brmsport-bin2-ledger.tsv")
  row <- L$outcome[match(v$id, L$id)]
  vv <- ifelse(v$outcome == "holds", "pass",
               m$verdict[match(v$id, m$id)])
  vc <- m$class[match(v$id, m$id)]
  c("| row | brms row | value check | verdict | message |",
    "|---|---|---|---|---|",
    sprintf("| `%s` | %s | %s | %s | %s |", v$id, row, v$outcome,
            ifelse(is.na(vv), "(row does not pass)",
                   ifelse(is.na(vc), vv, paste0(vv, " (", vc, ")"))),
            gsub("|", "/", v$msg, fixed = TRUE)),
    "", sprintf("%d value checks: %s.", nrow(v),
                paste(sprintf("%s %d", names(table(v$outcome)),
                              table(v$outcome)), collapse = ", ")))
}
blocks <- list(
  before = strip("dev/bin2port-before-log/summary.md"),
  after = strip("dev/brmsport-bin2-summary.md"),
  rate709 = rate709("dev/brmsport-bin2-summary.md"),
  rate709before = rate709("dev/bin2port-before-log/summary.md"),
  diff = readLines("dev/bin2port-out/diff.md", warn = FALSE),
  values = vtab(),
  tierbefore = tier("dev/bin2port-before-log"),
  tierafter = tier("dev/brmsport-bin2-log"),
  probe0 = logl("dev/bin2port-out/probe0.txt"),
  probe1 = logl("dev/bin2port-out/probe1.txt"),
  probe2 = logl("dev/bin2port-out/probe2.txt"),
  probe3 = logl("dev/bin2port-out/probe3.txt"),
  probe4 = logl("dev/bin2port-out/probe4.txt"),
  probe4b = logl("dev/bin2port-out/probe4b.txt"),
  probe5 = logl("dev/bin2port-out/probe5.txt"),
  defect = logl("dev/bin2port-out/defect-cumulative.txt"),
  guard = c("```", grep("^(lib:|RESULT|FAILING|  )",
                        readLines("dev/bin2port-out/guard-tier.txt",
                                  warn = FALSE), value = TRUE), "```"),
  guardmsg = c("```", grep("^(lib:|RESULT|FAILING|  )",
                           readLines("dev/bin2port-out/guard-msg.txt",
                                     warn = FALSE), value = TRUE), "```"),
  guardledger = logl("dev/bin2port-out/guard-ledger.txt"),
  tol = logl("dev/bin2port-out/tol.txt"),
  probe6 = logl("dev/bin2port-out/probe6.txt"),
  probe6b = logl("dev/bin2port-out/probe6b.txt"),
  md5 = c("```", readLines("dev/bin2port-out/md5-1.txt", warn = FALSE),
          "```", "", "Second run:", "", "```",
          readLines("dev/bin2port-out/md5-2.txt", warn = FALSE), "```"))
for (nm in names(blocks)) {
  b <- grep(sprintf("^<!-- begin %s -->$", nm), x)
  e <- grep(sprintf("^<!-- end %s -->$", nm), x)
  if (length(b) != 1L || length(e) != 1L || e < b) {
    stop("block ", nm, ": markers not found once")
  }
  x <- c(x[seq_len(b)], blocks[[nm]], x[e:length(x)])
}
writeLines(x, f)
cat("filled", length(blocks), "blocks\n")
