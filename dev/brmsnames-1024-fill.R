## Fill the generated blocks of dev/brmsnames-1024-findings.md from the
## lane's logs, so that every count in the record comes from the run
## that produced it.
##
##   Rscript dev/brmsnames-1024-fill.R
##
## A block is the text between `<!-- BEGIN <log> -->` and
## `<!-- END <log> -->`, where <log> is a file of dev/brmsnames-1024-log/
## or of dev/release/. A log that does not exist leaves its block alone
## and is named on stdout.
p <- "dev/brmsnames-1024-findings.md"
x <- readLines(p)
b <- grep("^<!-- BEGIN .* -->$", x)
for (i in rev(b)) {
  nm <- sub("^<!-- BEGIN (.*) -->$", "\\1", x[i])
  e <- which(x == paste0("<!-- END ", nm, " -->"))
  e <- e[e > i][1L]
  if (is.na(e)) stop("no END for ", nm)
  src <- c(file.path("dev/brmsnames-1024-log", nm),
           file.path("dev/release", nm))
  src <- src[file.exists(src)][1L]
  if (is.na(src)) {
    cat("missing:", nm, "\n")
    next
  }
  body <- readLines(src, warn = FALSE)
  body <- body[!grepl("^<!-- generated", body)]
  x <- c(x[seq_len(i)], "", body, "", x[e:length(x)])
}
writeLines(x, p)
cat("filled", length(b), "blocks\n")
