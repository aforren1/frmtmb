# Lane osa: fill the generated blocks of dev/osa-findings.md from the logs.
#
#   Rscript dev/osa-findings-gen.R
#
# A block is the lines between "<!-- gen:<name> -->" and "<!-- /gen -->".
# <name> is a log under dev/osa-out/ without its ".log"; the block gets
# the log's lines from its "BEGIN" line to its "END" line, or the whole
# log where it has none, inside a code fence. Rerun after every log
# changes; never type a number that a block carries.
root <- "C:/Users/adf44/source/r/frmtmb-wt-osa/dev"
p <- file.path(root, "osa-findings.md")
x <- readLines(p)
out <- character(0)
i <- 1L
while (i <= length(x)) {
  m <- regmatches(x[i], regexec("^<!-- gen:([A-Za-z0-9_./-]+) -->$", x[i]))[[1]]
  if (!length(m)) {
    out <- c(out, x[i])
    i <- i + 1L
    next
  }
  j <- i + 1L
  while (j <= length(x) && x[j] != "<!-- /gen -->") j <- j + 1L
  if (j > length(x)) stop("unclosed block ", m[2])
  lg <- file.path(root, "osa-out", paste0(m[2], ".log"))
  body <- if (file.exists(lg)) {
    l <- readLines(lg, warn = FALSE)
    b <- grep("^BEGIN", l)
    e <- grep("^END", l)
    if (length(b) && length(e)) l[b[1]:e[length(e)]] else l
  } else paste("MISSING LOG", lg)
  out <- c(out, x[i], "```", body, "```", "<!-- /gen -->")
  i <- j + 1L
}
writeLines(out, p)
cat("filled", sum(grepl("^<!-- gen:", out)), "blocks\n")
