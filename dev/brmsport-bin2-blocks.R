# Locate brms 2.23.0's bin-2 blocks and every assertion inside them, by
# parsing, and check the counts against the audit's: 59 assertions in 14
# blocks of bin 2T, 270 in 68 of bin 2P (dev/brms-suite-audit.md).
# Sourced by dev/brmsport-bin2-gen.R and dev/brmsport-bin2-ledger.R;
# run alone it prints the census.
#
#   Rscript dev/brmsport-bin2-blocks.R
#
# dev/brmsport-blocks.R holds the tarball check and the assertion
# counter; this file only widens its block finder from bin 1 to any
# set of bins, so bin 1's finder stays as it is.
source("dev/brmsport-blocks.R")

brms_bin_blocks <- function(bins) {
  cls <- utils::read.delim("dev/brmssuite-classified.tsv",
                           stringsAsFactors = FALSE)
  cls$ord <- stats::ave(seq_len(nrow(cls)), cls$file, FUN = seq_along)
  out <- list()
  for (f in unique(cls$file[cls$bin %in% bins])) {
    ex <- parse(file.path(suite, f), keep.source = TRUE)
    srcs <- attr(ex, "srcref")
    k <- 0L
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      if (!is.call(e) || !identical(e[[1]], as.name("test_that"))) next
      k <- k + 1L
      row <- cls[cls$file == f & cls$ord == k, ]
      stopifnot(nrow(row) == 1L)
      if (!row$bin %in% bins) next
      label <- eval(e[[2]])
      stopifnot(identical(gsub("[\r\n\t]+", " ", label), row$desc) ||
                  grepl("^paste", row$desc))
      body <- e[[3]]
      bsrc <- attr(body, "srcref")
      stmts <- as.list(body)[-1]
      lines <- vapply(bsrc[-1], function(s) as.integer(s[1]), 1L)
      cols <- vapply(bsrc[-1], function(s) as.integer(s[5]), 1L)
      texts <- vapply(bsrc[-1], function(s) paste(as.character(s),
                                                  collapse = "\n"), "")
      out[[length(out) + 1L]] <- list(
        file = f, ord = k, label = label, tier = row$tier, bin = row$bin,
        line = as.integer(srcs[[i]][1]), expects = row$expects,
        stmts = stmts, lines = lines, cols = cols, texts = texts)
    }
  }
  out
}

if (sys.nframe() == 0L) {
  for (bin in c("2T", "2P")) {
    blocks <- brms_bin_blocks(bin)
    tot <- 0L
    for (b in blocks) {
      n <- sum(vapply(b$stmts, n_expect, 1L))
      stopifnot(n == b$expects)
      tot <- tot + n
    }
    cat(sprintf("bin %s: %d blocks, %d assertions\n", bin, length(blocks),
                tot))
  }
}
