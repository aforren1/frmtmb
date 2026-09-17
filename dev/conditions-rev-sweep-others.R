# Reviewer, lane wt-conditions (round 2): every error an expect_error()
# caught on the lane build that is not a frmtmb_error, with the test
# expression that raised it.
#   Rscript dev/conditions-rev-sweep-others.R <lane sweep dir>
fs <- list.files(commandArgs(TRUE)[1L], "rds$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, readRDS))
o <- x[grepl("^expect_error", x$key) & !grepl("frmtmb_error", x$class), ]
cat("rows", nrow(o), "\n\n")
for (i in seq_len(nrow(o))) {
  cat(sprintf("[%d] %s\n    %s\n    %s | %s\n", i,
              sub(".*/tests/testthat/", "", o$file[i]),
              substr(sub("^expect_error ", "", o$key[i]), 1, 150),
              o$class[i], substr(o$message[i], 1, 110)))
}
