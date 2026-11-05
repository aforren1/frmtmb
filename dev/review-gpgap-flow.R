# Short check of fix round 1, item 4: all 172 pairs through the new flow, measured
# here on the lane's library (dev/review-gpgap-equiv-out), against the coordinator's
# baseline on 0.78.0. Which verdicts moved, before and after adjudication? Do the
# headline figures come from the flow's own verdicts?
#
#   Rscript dev/review-gpgap-flow.R
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
mine <- "dev/review-gpgap-equiv-out"
ids <- read.delim("dev/equiv-data/sample.tsv", stringsAsFactors = FALSE)$id
rd <- function(d, i) {
  p <- file.path(d, paste0(i, ".tsv"))
  if (!file.exists(p)) return(NULL)
  read.delim(p, stringsAsFactors = FALSE, na.strings = "NA")
}
ch <- character()
n <- 0L
for (i in ids) {
  a <- rd(base, i); b <- rd(mine, i)
  if (is.null(a) || is.null(b)) {
    cat("MISSING", i, is.null(a), is.null(b), "\n")
    next
  }
  n <- n + 1L
  same <- identical(as.character(a$verdict), as.character(b$verdict)) &&
    identical(as.character(a$level), as.character(b$level)) &&
    identical(as.logical(a$agree), as.logical(b$agree))
  if (!same) {
    ch <- c(ch, i)
    cat(sprintf("MOVED %s: %s (%s, agree %s) -> %s (%s, agree %s)\n    reason: %s\n", i,
                a$verdict, a$level, a$agree, b$verdict, b$level, b$agree,
                substr(b$reason, 1, 330)))
  }
}
cat(sprintf("FLOW pairs compared %d, moved %d: %s\n", n, length(ch), paste(ch, collapse = " ")))
b136 <- rd(mine, "p136")
cat(sprintf("FLOW p136: %s, level %s, l4 z %s, rhat %s, l1_const %s, reason: %s\n", b136$verdict, b136$level,
            b136$l4_z, signif(b136$l4_rhat, 4), signif(b136$l1_const, 4), substr(b136$reason, 1, 90)))
# the summary, on a scratch tree, adjudication less the four rows
root <- "dev/review-gpgap-rates"
unlink(root, recursive = TRUE)
dir.create(file.path(root, "dev"), recursive = TRUE)
file.copy("dev/equiv-data", file.path(root, "dev"), recursive = TRUE)
file.copy("dev/equiv-summary.R", file.path(root, "dev"))
adj <- read.delim("dev/equiv-data/adjudication.tsv", stringsAsFactors = FALSE,
                  colClasses = "character", na.strings = character())
utils::write.table(adj[!adj$id %in% c("p022", "p023", "p094", "p109"), ],
                   file.path(root, "dev/equiv-data/adjudication.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
o <- file.path(root, "out")
dir.create(o)
file.copy(Sys.glob(file.path(mine, "p*.tsv")), o)
wd <- getwd()
setwd(root)
system2("C:/Program Files/R/R-4.6.1/bin/Rscript.exe", c("dev/equiv-summary.R", "out"),
        stdout = FALSE, stderr = FALSE)
setwd(wd)
s <- readLines(file.path(o, "summary.md"), warn = FALSE)
for (p in c("the headline", "^- no defect", "sensitivity, weak", "among measured")) {
  cat(sprintf("RATE %s\n", substr(grep(p, s, value = TRUE)[1], 1, 90)))
}
cat(sprintf("RATE the rewrite named in the summary: %d lines\n",
            length(grep("rewritten", s))))
x <- s[grep("exact gp\\(\\) pairs read L1 twice", s) + 2:6]
cat(substr(x, 1, 230), sep = "\n")
unlink(root, recursive = TRUE)
