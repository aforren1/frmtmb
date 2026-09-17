# Locate brms 2.23.0's bin-1 blocks and every assertion inside them, by
# parsing, and check the count against the audit's 494 in 95 blocks.
# Sourced by dev/brmsport-gen.R; run alone it prints the census.
#
# A block is matched to dev/brmssuite-classified.tsv by FILE and by its
# ORDINAL among that file's test_that() calls, because two files repeat
# a label (tests.stop2.R) and a label alone would be ambiguous.
suite <- "dev/brms-suite/brms/tests/testthat"
stopifnot(dir.exists(suite))
local({
  # the sha256 recorded in dev/brms-suite-audit.md section 1; the
  # extracted tree is trusted only because it came from this tarball
  tb <- "dev/brms-suite/brms_2.23.0.tar.gz"
  h <- digest::digest(file = tb, algo = "sha256")
  if (!identical(h, paste0("b5f5bb5604ec3f87b3ad99f0da4e6a37",
                           "d8c5488221e8abe459171340992c37a5"))) {
    stop("brms tarball sha256 mismatch: ", h)
  }
})

is_expect_name <- function(nm) {
  grepl("^expect_[a-z0-9_]+$", nm) & !grepl("^expected", nm)
}
n_expect <- function(e) {
  nms <- all.names(e, functions = TRUE, unique = FALSE)
  sum(is_expect_name(nms))
}

brms_bin1_blocks <- function() {
  cls <- utils::read.delim("dev/brmssuite-classified.tsv",
                           stringsAsFactors = FALSE)
  cls$ord <- stats::ave(seq_len(nrow(cls)), cls$file, FUN = seq_along)
  out <- list()
  for (f in unique(cls$file[cls$bin == "1"])) {
    ex <- parse(file.path(suite, f), keep.source = TRUE)
    srcs <- attr(ex, "srcref")
    k <- 0L
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      if (!is.call(e) || !identical(e[[1]], as.name("test_that"))) next
      k <- k + 1L
      row <- cls[cls$file == f & cls$ord == k, ]
      stopifnot(nrow(row) == 1L)
      if (row$bin != "1") next
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
        file = f, ord = k, label = label, tier = row$tier,
        line = as.integer(srcs[[i]][1]), expects = row$expects,
        stmts = stmts, lines = lines, cols = cols, texts = texts)
    }
  }
  out
}

if (sys.nframe() == 0L) {
  blocks <- brms_bin1_blocks()
  tot <- 0L
  for (b in blocks) {
    n <- sum(vapply(b$stmts, n_expect, 1L))
    stopifnot(n == b$expects)
    tot <- tot + n
    for (j in seq_along(b$stmts)) {
      s <- b$stmts[[j]]
      k <- n_expect(s)
      top <- is.call(s) && is_expect_name(deparse(s[[1]])[1])
      if (k > 0 && (!top || k > 1)) {
        cat(sprintf("NONTRIVIAL %s:%d k=%d top=%s: %s\n", b$file,
                    b$lines[j], k, top, substr(b$texts[j], 1, 70)))
      }
    }
  }
  cat(length(blocks), "blocks,", tot, "assertions\n")
}
