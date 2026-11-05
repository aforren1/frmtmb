# Verdict changes of the rerun pairs against the baseline on 0.78.0
# (frmtmb-wt-r16base/dev/r16base-equiv-out, read-only). Writes
# dev/nlfix-equiv-compare.txt as a markdown table.
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
ids <- readLines("dev/nlfix-equiv-ids.txt")
rd <- function(d, id) {
  f <- file.path(d, paste0(id, ".tsv"))
  if (!file.exists(f)) return(NULL)
  read.delim(f, stringsAsFactors = FALSE)
}
rows <- lapply(ids, function(id) {
  a <- rd(base, id)
  b <- rd("dev/nlfix-equiv-out", id)
  data.frame(id = id, pair = paste(a$feature_a, a$feature_b, sep = " x "),
             before = a$verdict, before_level = a$level,
             after = if (is.null(b)) "MISSING" else b$verdict,
             after_level = if (is.null(b)) "" else b$level,
             l1_const = if (is.null(b)) NA else format(b$l1_const, digits = 3),
             l1_grad = if (is.null(b)) NA else format(b$l1_grad, digits = 3),
             frm_error = if (is.null(b)) "" else substr(ifelse(is.na(b$frm_error), "", b$frm_error), 1, 60),
             stringsAsFactors = FALSE)
})
tb <- do.call(rbind, rows)
md <- function(df) {
  df[] <- lapply(df, function(v) gsub("|", "\\|", ifelse(is.na(v), "", v), fixed = TRUE))
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
ch <- tb[tb$before != tb$after, ]
out <- c(sprintf("%d pairs rerun, %d changed verdict.", nrow(tb), nrow(ch)), "",
         "Changed:", "", md(ch), "", "All rerun pairs:", "", md(tb))
writeLines(out, "dev/nlfix-equiv-compare.txt")
print(ch, row.names = FALSE)
if (file.exists("dev/nlfix-equiv-out/summary.md")) {
  s <- readLines("dev/nlfix-equiv-out/summary.md")
  cat(grep("headline|verdicts by|^\\| (untested|conditional|works|refused) ", s,
           value = TRUE, ignore.case = TRUE), sep = "\n")
}
