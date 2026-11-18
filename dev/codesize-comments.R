# Part of Q1/Q7 of the codesize lane: what the 34,000 lines of internal
# commentary (roxygen on internal objects plus other comments) hold.
# Counts comment lines that carry project history or provenance (a
# version, "used to", a backlog, lane, round or dev/ reference), lines
# that cite a measured number or case, and the comment-to-code ratio per
# file.
#   Rscript dev/codesize-comments.R
# Reads dev/codesize-out/parsed.rds and files.tsv.
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
fl <- read.delim(file.path(out, "files.tsv"), stringsAsFactors = FALSE)
exported <- unique(c(P$exports, P$s3$method))
hist_pat <- paste0("\\b(v?0\\.[0-9]{2}(\\.[0-9])?)\\b|used to|no longer|",
                   "previously|\\bwas\\b.*\\b(now|before)\\b|backlog|",
                   "\\blane\\b|\\bround\\b of|dev/|item [0-9]|review")
case_pat <- "[0-9]\\.[0-9]{2,}|[0-9]e-[0-9]|seed"
res <- list()
for (f in names(P$core)) {
  x <- P$core[[f]]
  cm <- x$pd[x$pd$token == "COMMENT", ]
  rox <- grepl("^#'", cm$text)
  # internal roxygen = roxygen lines in blocks above a non-exported
  # function or with @noRd; approximate with the per-file split from
  # files.tsv and classify every comment line by pattern
  res[[f]] <- data.frame(
    file = paste0(f, ".R"),
    comment_lines = nrow(cm),
    history = sum(grepl(hist_pat, cm$text, perl = TRUE, ignore.case = TRUE)),
    numbers = sum(grepl(case_pat, cm$text, perl = TRUE)),
    rox_noRd_tags = sum(grepl("@noRd", cm$text)))
}
r <- do.call(rbind, res)
m <- merge(r, fl[c("file", "code", "rox_internal", "comment", "rox_exported")],
           by = "file")
m$internal_comment <- m$rox_internal + m$comment
m$ratio <- round(m$internal_comment / pmax(m$code, 1), 2)
cat("BEGIN codesize-comments\n")
cat(sprintf("history pattern: %s\n", hist_pat))
cat(sprintf("All comment lines %d (roxygen of exports %d, internal roxygen %d, other %d); with a history/provenance marker %d (%.1f%%); with a measured number %d\n",
            sum(m$comment_lines), sum(m$rox_exported), sum(m$rox_internal),
            sum(m$comment), sum(m$history),
            100 * sum(m$history) / sum(m$comment_lines), sum(m$numbers)))
cat(sprintf("Internal commentary (internal roxygen + other comments) per code line: %.2f\n",
            sum(m$internal_comment) / sum(m$code)))
m <- m[order(-m$internal_comment), ]
cat("Top 20 files by internal commentary:\n")
cat("| file | code | internal roxygen | other comments | per code line | history lines |\n|---|---|---|---|---|---|\n")
x <- head(m, 20)
cat(sprintf("| %s | %d | %d | %d | %.2f | %d |\n", x$file, x$code,
            x$rox_internal, x$comment, x$ratio, x$history), sep = "")
cat("END codesize-comments\n")
