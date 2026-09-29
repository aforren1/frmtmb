# Does `[[` on a list return NULL or error for a missing name / an
# out-of-range index? cs_newdata_columns() relies on the answer.
p <- function(lab, ex) {
  cat(sprintf("%-34s %s\n", lab,
              tryCatch(paste(utils::capture.output(str(ex)),
                             collapse = " "),
                       error = function(e) paste("ERROR:",
                                                 conditionMessage(e)))))
}
l <- list(a = 1)
p('list(a=1)[["b"]]', l[["b"]])
p('list()[["b"]]', list()[["b"]])
p('list(a=1)[[3]]', l[[3]])
p('list(a=1)[[2]]', l[[2]])
e <- list()
e[["7"]] <- 5
p('after e[["7"]] <- 5, e[["7"]]', e[["7"]])
p('e[["8"]]', e[["8"]])
