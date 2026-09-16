# Counts test blocks and assertions in brms's own suite by PARSING each
# file, because a grep over expect_*() misses calls split across lines
# and hits them inside comments and strings.
dir <- "dev/brms-suite/brms/tests/testthat"
files <- sort(list.files(dir, pattern = "[.]R$", full.names = TRUE))

rows <- list()
kinds <- character(0)
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  # all.names walks the parsed AST, so it sees calls split over several
  # lines and never sees a comment or a string literal.
  nms <- all.names(ex, functions = TRUE, unique = FALSE)
  # expected_ar_mat and friends are data, not assertions, so the verb
  # pattern must not be a bare prefix match on "expect".
  which_exp <- grep("^expect_[a-z0-9_]+$", nms, value = TRUE)
  which_exp <- which_exp[!grepl("^expected", which_exp)]
  acc <- list(blocks = sum(nms == "test_that"),
              expects = length(which_exp))
  kinds <- c(kinds, which_exp)
  rows[[length(rows) + 1L]] <- data.frame(
    file = basename(f),
    lines = length(readLines(f, warn = FALSE)),
    blocks = acc$blocks,
    expects = acc$expects,
    stringsAsFactors = FALSE
  )
}
# A per-block dump, because three files mix bins inside one file and a
# file-level label would hide that.
blocks <- list()
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  for (i in seq_along(ex)) {
    e <- ex[[i]]
    if (!is.call(e) || !identical(e[[1]], as.name("test_that"))) next
    nms <- all.names(e, functions = TRUE, unique = FALSE)
    v <- grep("^expect_[a-z0-9_]+$", nms, value = TRUE)
    v <- v[!grepl("^expected", v)]
    blocks[[length(blocks) + 1L]] <- data.frame(
      file = basename(f),
      desc = gsub("[\r\n\t]+", " ", paste(as.character(e[[2]]),
                                          collapse = " ")),
      expects = length(v),
      match2 = sum(v == "expect_match2"),
      stringsAsFactors = FALSE
    )
  }
}
write.table(do.call(rbind, blocks), "dev/brmssuite-blocks.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

tab <- do.call(rbind, rows)
tab <- tab[order(-tab$expects), ]
write.table(tab, "dev/brmssuite-counts.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat(sprintf("%-32s %6s %7s %8s\n", "file", "lines", "blocks", "expects"))
for (i in seq_len(nrow(tab))) {
  cat(sprintf("%-32s %6d %7d %8d\n", tab$file[i], tab$lines[i],
              tab$blocks[i], tab$expects[i]))
}
cat(sprintf("%-32s %6d %7d %8d\n", "TOTAL", sum(tab$lines),
            sum(tab$blocks), sum(tab$expects)))
cat("\nassertion verbs used:\n")
print(sort(table(kinds), decreasing = TRUE))
