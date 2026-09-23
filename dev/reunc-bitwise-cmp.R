# Lane wt-reunc: compare two dev/reunc-bitwise.R outputs with
# identical(), entry by entry.
#
#   Rscript dev/reunc-bitwise-cmp.R <base.rds> <lane.rds>
a <- commandArgs(trailingOnly = TRUE)
b <- readRDS(a[1])
l <- readRDS(a[2])
# these the task changes on purpose: the conditional predict()
expect_move <- c("predict_null")
same <- 0L
moved <- 0L
bad <- 0L
for (k in union(names(b), names(l))) {
  for (e in union(names(b[[k]]), names(l[[k]]))) {
    x <- b[[k]][[e]]
    y <- l[[k]][[e]]
    ok <- identical(x, y)
    tag <- if (ok) "IDENTICAL" else if (e %in% expect_move) "moved" else {
      "DIFFERS"
    }
    if (ok) same <- same + 1L else if (e %in% expect_move) {
      moved <- moved + 1L
    } else bad <- bad + 1L
    extra <- ""
    if (!ok && is.numeric(unlist(x)) && is.numeric(unlist(y)) &&
          length(unlist(x)) == length(unlist(y))) {
      extra <- sprintf(" max abs diff %.3e",
                       max(abs(unlist(x) - unlist(y)), na.rm = TRUE))
    }
    if (inherits(x, "grab_error") || inherits(y, "grab_error")) {
      extra <- paste0(extra, " [error: ",
                      substr(if (inherits(x, "grab_error")) x else y, 1, 60),
                      "]")
    }
    cat(sprintf("%-16s %-22s %s%s\n", k, e, tag, extra))
  }
}
cat(sprintf("\nidentical %d, moved as intended %d, differs %d\n", same,
            moved, bad))
