# The nl x family sweep on this build against the record.
a <- utils::read.delim("dev/equiv-data/nl-sweep.tsv")
b <- utils::read.delim("dev/nlstart-out/nlsweep/nl-sweep.tsv")
b <- b[match(a$family, b$family), ]
num <- c("frm_linear", "frm_nl", "frm_intercept", "default_linear",
         "default_nl_default_start", "default_nl_positive_start")
for (k in num) {
  d <- abs(a[[k]] - b[[k]])
  cat(sprintf("%-26s rows %d, differ > 1e-6: %d, NA changed: %d\n", k,
              nrow(a), sum(d > 1e-6, na.rm = TRUE),
              sum(is.na(a[[k]]) != is.na(b[[k]]))))
}
for (k in c("outcome", "default_outcome", "brms_nl")) {
  cat(sprintf("%-26s changed: %d\n", k,
              sum(as.character(a[[k]]) != as.character(b[[k]]), na.rm = TRUE)))
}
