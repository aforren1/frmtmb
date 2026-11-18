# Q2 (coverage half) of the codesize lane: lines of core R/ that no core
# test runs, per file and subsystem, and the functions no test enters.
#   Rscript dev/codesize-covsum.R [cov.rds]
# Reads dev/codesize-cov/cov.rds (dev/codesize-covr-env.R) or the given
# snapshot, dev/codesize-cov/results.tsv, and dev/codesize-out/*.tsv.
a <- commandArgs(trailingOnly = TRUE)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-codesize"
cf <- if (length(a)) a[1] else file.path(wt, "dev/codesize-cov/cov.rds")
out <- file.path(wt, "dev/codesize-out")
.libPaths(c("C:/Users/adf44/source/r/wt-codesize-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
invisible(loadNamespace("covr"))
cov <- readRDS(cf)
df <- as.data.frame(cov)
df$file <- basename(df$filename)
# a line is covered when any expression on it ran
ln <- aggregate(value ~ file + first_line, df, max)
ln$hit <- ln$value > 0
pf <- aggregate(cbind(relevant = 1, covered = ln$hit) ~ file, ln, sum)
pf$uncovered <- pf$relevant - pf$covered
pf$pct <- round(100 * pf$covered / pf$relevant, 1)
P <- readRDS(file.path(out, "parsed.rds"))
pf$subsystem <- P$sub_of[sub("[.]R$", "", pf$file)]
pf <- pf[order(-pf$uncovered), ]
write.table(pf, file.path(out, "coverage-files.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
rt <- read.delim(file.path(wt, "dev/codesize-cov/results.tsv"))
cat("BEGIN codesize-covsum\n")
cat(sprintf("source: %s; test files run %d of %d; pass %d, fail %d, error files %d, skip %d; %.0f s summed under covr\n",
            basename(cf), nrow(rt),
            length(list.files(file.path(wt, "tests/testthat"), "^test-.*[.]R$")),
            sum(rt$pass, na.rm = TRUE), sum(rt$fail, na.rm = TRUE),
            sum(rt$err > 0, na.rm = TRUE), sum(rt$skip, na.rm = TRUE),
            sum(rt$sec, na.rm = TRUE)))
cat("files with an error:", paste(rt$file[rt$err > 0 & !is.na(rt$err)],
                                  collapse = " "), "\n")
cat(sprintf("relevant lines %d, covered %d (%.1f%%), uncovered %d\n",
            sum(pf$relevant), sum(pf$covered),
            100 * sum(pf$covered) / sum(pf$relevant), sum(pf$uncovered)))
ps <- aggregate(cbind(relevant, covered, uncovered) ~ subsystem, pf, sum)
ps$pct <- round(100 * ps$covered / ps$relevant, 1)
ps <- ps[order(-ps$uncovered), ]
cat("\n| subsystem | relevant | uncovered | covered % |\n|---|---|---|---|\n")
cat(sprintf("| %s | %d | %d | %.1f |\n", ps$subsystem, ps$relevant,
            ps$uncovered, ps$pct), sep = "")
cat("\nPer file (all files, most uncovered first):\n")
cat("| file | relevant | uncovered | covered % |\n|---|---|---|---|\n")
cat(sprintf("| %s | %d | %d | %.1f |\n", pf$file, pf$relevant, pf$uncovered,
            pf$pct), sep = "")

# functions: relevant lines inside each top-level function
fn <- read.delim(file.path(out, "callgraph.tsv"), stringsAsFactors = FALSE)
fn$relevant <- 0L; fn$covered <- 0L
for (f in unique(fn$file)) {
  l <- ln[ln$file == f, ]
  i <- which(fn$file == f)
  for (j in i) {
    k <- l$first_line >= fn$first[j] & l$first_line <= fn$last[j]
    fn$relevant[j] <- sum(k); fn$covered[j] <- sum(l$hit[k])
  }
}
z <- fn[fn$relevant > 0 & fn$covered == 0, ]
z <- z[order(-z$relevant), ]
write.table(z, file.path(out, "coverage-zero-functions.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("\nFunctions no core test enters: %d of %d with relevant lines, %d code lines (%d relevant lines); exported among them: %d\n",
            nrow(z), sum(fn$relevant > 0), sum(z$code), sum(z$relevant),
            sum(z$exported)))
zs <- aggregate(cbind(n = 1, code = z$code) ~ subsystem, z, sum)
zs <- zs[order(-zs$code), ]
cat(sprintf("  %s: %d functions, %d code lines\n", zs$subsystem, zs$n,
            zs$code), sep = "")
cat("Largest 30:\n")
x <- head(z, 30)
cat(sprintf("  %s (%s:%d) %d code lines%s\n", x$name, x$file, x$first,
            x$code, ifelse(x$exported, ", exported", "")), sep = "")
# A function reached only as a value in a registry list (family_registry,
# covstruct slots) is stored in the lazy-load database as its own copy, so
# covr, which traces the namespace binding, never sees it run. Those are
# not evidence of untested code.
nv <- read.delim(file.path(out, "callgraph-novalue.tsv"),
                 stringsAsFactors = FALSE)
zi <- z[!z$exported, ]
art <- zi$name %in% nv$name
cat(sprintf("\nInternal functions no test enters: %d (%d code lines); reached only as a value, so possibly a covr artifact: %d (%d code lines: %s); the rest: %d (%d code lines)\n",
            nrow(zi), sum(zi$code), sum(art), sum(zi$code[art]),
            paste(zi$name[art], collapse = " "), sum(!art),
            sum(zi$code[!art])))
y <- zi[!art, ]
cat(sprintf("  %s (%s:%d) %d code lines\n", y$name, y$file, y$first, y$code),
    sep = "")
cat("END codesize-covsum\n")
