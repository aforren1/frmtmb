# Assertion verbs per file. The audit needs this because a file's bin
# often follows its verb mix: expect_match2 is Stan code, expect_error
# is validation, expect_equal on a number is the bin that needs an ML
# restatement.
dir <- "dev/brms-suite/brms/tests/testthat"
files <- sort(list.files(dir, pattern = "[.]R$", full.names = TRUE))
rows <- list()
for (f in files) {
  ex <- parse(f, keep.source = FALSE)
  nms <- all.names(ex, functions = TRUE, unique = FALSE)
  v <- grep("^expect_[a-z0-9_]+$", nms, value = TRUE)
  v <- v[!grepl("^expected", v)]
  tb <- table(v)
  rows[[basename(f)]] <- tb
}
verbs <- sort(unique(unlist(lapply(rows, names))))
m <- matrix(0L, nrow = length(rows), ncol = length(verbs),
            dimnames = list(names(rows), verbs))
for (i in names(rows)) m[i, names(rows[[i]])] <- as.integer(rows[[i]])
keep <- c("expect_equal", "expect_equivalent", "expect_identical",
          "expect_match2", "expect_error", "expect_warning",
          "expect_message", "expect_true", "expect_false", "expect_output")
m2 <- cbind(m[, intersect(keep, colnames(m)), drop = FALSE],
            other = rowSums(m) - rowSums(m[, intersect(keep, colnames(m)),
                                           drop = FALSE]))
print(m2)
cat("\ncolumn totals\n")
print(colSums(m2))
