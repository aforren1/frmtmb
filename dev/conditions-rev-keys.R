# Reviewer, lane wt-conditions (round 2): the test expressions that
# raised a classed warning or message in a sweep, to pick triggers for
# the printed-output harness.
#   Rscript dev/conditions-rev-keys.R <lane sweep dir>
fs <- list.files(commandArgs(TRUE)[1L], "rds$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, readRDS))
w <- x[grepl("^expect_(warning|message)", x$key) &
         grepl("frmtmb_", x$class), ]
w <- w[!duplicated(w$message), ]
for (i in seq_len(nrow(w))) {
  cat(sub(".*testthat/", "", w$file[i]), "|", substr(w$key[i], 1, 170),
      "|", w$class[i], "|", substr(w$message[i], 1, 50), "\n")
}
