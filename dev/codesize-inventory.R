# Q1 of the codesize lane: where the size of core R/ is, per file and per
# subsystem, split into code, roxygen of exported objects, roxygen of
# internal objects, other comments and blank lines; and how many Rd pages
# document internal functions.
#   Rscript dev/codesize-inventory.R
# Reads dev/codesize-out/parsed.rds (dev/codesize-parse.R). Writes
# dev/codesize-out/files.tsv, subsystems.tsv, functions.tsv, rd.tsv,
# and prints the summary blocks pasted into dev/codesize-findings.md.
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
exported <- unique(c(P$exports, P$s3$method))

line_kinds <- function(x) {
  n <- length(x$lines)
  kind <- rep("blank", n)
  cm <- x$pd[x$pd$token == "COMMENT", ]
  rox <- grepl("^#'", cm$text)
  kind[cm$line1[!rox]] <- "comment"
  kind[cm$line1[rox]] <- "roxygen"
  nc <- x$pd[x$pd$token != "COMMENT", ]
  # a token may span lines (a long string); every spanned line is code
  cl <- unlist(Map(seq, nc$line1, nc$line2))
  kind[unique(cl)] <- "code"
  kind
}

# One roxygen block = a run of roxygen lines; it documents the next
# top-level expression.
rox_blocks <- function(x, kind) {
  r <- rle(kind == "roxygen")
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  keep <- r$values
  if (!any(keep)) return(NULL)
  b <- data.frame(start = starts[keep], end = ends[keep])
  b$target <- vapply(b$end, function(e) {
    i <- which(x$top$first > e)[1]
    if (is.na(i)) NA_character_ else x$top$name[i] %||% NA_character_
  }, "")
  b$target_fun <- vapply(b$end, function(e) {
    i <- which(x$top$first > e)[1]
    !is.na(i) && x$top$is_fun[i]
  }, TRUE)
  txt <- vapply(seq_len(nrow(b)), function(i)
    paste(x$lines[b$start[i]:b$end[i]], collapse = "\n"), "")
  b$noRd <- grepl("@noRd", txt)
  b$export_tag <- grepl("@export", txt)
  b$kw_internal <- grepl("@keywords\\s+internal", txt)
  b$topic <- grepl("@name|@docType|_PACKAGE", txt) & !b$target_fun
  b$cls <- ifelse(b$noRd, "rox_internal",
           ifelse(b$kw_internal, "rox_internal",
           ifelse(b$export_tag | b$topic | b$target %in% exported,
                  "rox_exported", "rox_internal")))
  b
}
`%||%` <- function(a, b) if (is.null(a)) b else a

files <- list(); funs <- list()
for (f in names(P$core)) {
  x <- P$core[[f]]
  k <- line_kinds(x)
  b <- rox_blocks(x, k)
  rk <- k
  if (!is.null(b)) for (i in seq_len(nrow(b)))
    rk[b$start[i]:b$end[i]] <- b$cls[i]
  tb <- table(factor(rk, c("code", "rox_exported", "rox_internal",
                           "comment", "blank")))
  files[[f]] <- data.frame(file = paste0(f, ".R"), subsystem = P$sub_of[[f]],
                           lines = length(k), as.list(tb),
                           bytes = file.size(file.path(P$wt, "R",
                                                       paste0(f, ".R"))),
                           top_funs = sum(x$top$is_fun))
  t <- x$top[x$top$is_fun, ]
  if (nrow(t)) {
    t$file <- paste0(f, ".R")
    t$subsystem <- P$sub_of[[f]]
    t$code <- vapply(seq_len(nrow(t)), function(i)
      sum(k[t$first[i]:t$last[i]] == "code"), 1L)
    t$comment <- vapply(seq_len(nrow(t)), function(i)
      sum(k[t$first[i]:t$last[i]] == "comment"), 1L)
    # the roxygen block right above the function
    t$rox <- vapply(t$first, function(s) {
      if (is.null(b)) return(0)
      j <- which(b$end < s & b$end >= s - 3L)
      if (length(j)) as.numeric(b$end[j[1]] - b$start[j[1]] + 1) else 0
    }, 1)
    t$rox_noRd <- vapply(t$first, function(s) {
      if (is.null(b)) return(FALSE)
      j <- which(b$end < s & b$end >= s - 3L)
      length(j) > 0 && b$noRd[j[1]]
    }, TRUE)
    t$exported <- t$name %in% exported
    funs[[f]] <- t
  }
}
files <- do.call(rbind, files)
funs <- do.call(rbind, funs)
rownames(files) <- rownames(funs) <- NULL
write.table(files, file.path(out, "files.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
write.table(funs, file.path(out, "functions.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

num <- c("lines", "code", "rox_exported", "rox_internal", "comment", "blank",
         "bytes", "top_funs")
sb <- aggregate(files[num], files["subsystem"], sum)
sb$files <- as.vector(table(files$subsystem)[sb$subsystem])
sb <- sb[order(-sb$lines), ]
sb$share <- round(100 * sb$lines / sum(sb$lines), 1)
write.table(sb, file.path(out, "subsystems.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# Rd pages: internal when keyworded internal, or when no alias is exported
# and the page is not a package, data or topic page
man <- file.path(P$wt, "man")
rd <- do.call(rbind, lapply(list.files(man, "[.]Rd$", full.names = TRUE),
                            function(p) {
  r <- tools::parse_Rd(p)
  tags <- vapply(r, attr, "", "Rd_tag")
  al <- unlist(lapply(r[tags == "\\alias"], as.character))
  kw <- unlist(lapply(r[tags == "\\keyword"], as.character))
  dt <- unlist(lapply(r[tags == "\\docType"], as.character))
  us <- any(tags == "\\usage")
  data.frame(rd = basename(p), aliases = length(al),
             exported_alias = any(al %in% exported),
             kw_internal = "internal" %in% kw,
             doctype = if (length(dt)) dt[1] else "",
             has_usage = us, lines = length(readLines(p, warn = FALSE)))
}))
rd$internal_page <- rd$kw_internal |
  (!rd$exported_alias & rd$has_usage & rd$doctype == "")
write.table(rd, file.path(out, "rd.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

cat("BEGIN codesize-inventory\n")
tot <- colSums(files[num])
cat(sprintf("core R/: %d files, %d lines, %.0f KB; code %d, roxygen exported %d, roxygen internal %d, other comments %d, blank %d; top-level functions %d\n",
            nrow(files), tot[["lines"]], tot[["bytes"]] / 1024, tot[["code"]],
            tot[["rox_exported"]], tot[["rox_internal"]], tot[["comment"]],
            tot[["blank"]], tot[["top_funs"]]))
cat("\nPer subsystem (lines; share of all lines):\n")
cat("| subsystem | files | lines | share % | code | rox exported | rox internal | other comments | blank | functions |\n|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(sb))) with(sb[i, ], cat(sprintf(
  "| %s | %d | %d | %.1f | %d | %d | %d | %d | %d | %d |\n", subsystem,
  files, lines, share, code, rox_exported, rox_internal, comment, blank,
  top_funs)))
cat("\nLargest 30 files:\n")
cat("| file | subsystem | lines | code | rox exported | rox internal | other comments | functions |\n|---|---|---|---|---|---|---|---|\n")
fo <- files[order(-files$lines), ][1:30, ]
for (i in seq_len(nrow(fo))) with(fo[i, ], cat(sprintf(
  "| %s | %s | %d | %d | %d | %d | %d | %d |\n", file, subsystem, lines,
  code, rox_exported, rox_internal, comment, top_funs)))
cat(sprintf("\nFunctions: %d top-level, %d exported or registered S3, %d internal; internal functions with a roxygen block: %d (%d of them @noRd), roxygen lines on internal functions: %d\n",
    nrow(funs), sum(funs$exported), sum(!funs$exported),
    sum(!funs$exported & funs$rox > 0), sum(!funs$exported & funs$rox_noRd),
    sum(funs$rox[!funs$exported])))
q <- quantile(funs$code, c(.25, .5, .75, .9, .99))
cat("Code lines per function, quartiles and 90/99%:", q, " max", max(funs$code), "\n")
cat(sprintf("\nRd pages: %d; internal (keyword internal, or no exported alias): %d (%d keyword internal, %d with no exported alias), %d Rd lines\n",
    nrow(rd), sum(rd$internal_page), sum(rd$kw_internal),
    sum(!rd$exported_alias & rd$has_usage & rd$doctype == ""),
    sum(rd$lines[rd$internal_page])))
cat("Internal Rd pages:", paste(sub("[.]Rd$", "", rd$rd[rd$internal_page]),
                                collapse = " "), "\n")
cat("END codesize-inventory\n")
