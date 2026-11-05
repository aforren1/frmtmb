# Release fix RF1 of the refuse review, on the merged build: the rule as the
# package has it, on the data of dev/review-refuse-f2-rule3.R (same seeds, same
# order). Each fit's own warning state, read from fit$mi_laplace, is joined by tag
# to the exact shift the reviewer measured (dev/rel079-rf1-score.R).
#   REV_LIB=<rellib-r17> Rscript dev/rel079-rf1.R table1 <family> | table2 <slope> | vm
args0 <- commandArgs(trailingOnly = TRUE)
txt <- readLines("dev/review-refuse-f2-rule.R")
i <- which(startsWith(txt, "mode <- args[1]"))
eval(parse(text = sub("args <- commandArgs(trailingOnly = TRUE)", "args <- \"none\"",
                      txt[seq_len(i - 1L)], fixed = TRUE)), envir = globalenv())
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
out <- function(tag, a) {
  ml <- a$mi_laplace
  mv <- if (is.null(ml$moves)) NA_real_ else max(ml$moves)
  w <- isTRUE(!is.null(ml$window)) || isTRUE(abs(ml$tk) > 1 && mv >= 0.5)
  pr("PKG %s tk %s mv %s warn %s window %s", tag, format(ml$tk, digits = 4),
     format(mv, digits = 4), w, format(ml$window %||% NA, digits = 3))
}
r3 <- readLines("dev/review-refuse-f2-rule3.R")
j <- which(startsWith(r3, "m <- args0[1]"))
eval(parse(text = r3[j:length(r3)]))
