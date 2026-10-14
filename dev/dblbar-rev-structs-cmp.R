# Reviewer, lane dblbar: base against lane for dev/dblbar-rev-structs.R.
# Output: dev/dblbar-rev-log/structs-cmp.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
b <- readRDS(file.path(wt, "dev/dblbar-rev-log/structs-base.rds"))
l <- readRDS(file.path(wt, "dev/dblbar-rev-log/structs-lane.rds"))
sink(file.path(wt, "dev/dblbar-rev-log/structs-cmp.txt"), split = TRUE)
for (nm in names(b)) {
  x <- b[[nm]]
  y <- l[[nm]]
  if (!is.null(x$error) || !is.null(y$error)) {
    cat(sprintf("%-8s error base %s | lane %s\n", nm,
                !is.null(x$error), !is.null(y$error)))
    next
  }
  same <- vapply(c("vars", "blocks", "rand", "vc", "cv", "hyp"),
                 function(k) identical(x[[k]], y[[k]]), NA)
  cat(sprintf("%-8s ll rel %.2e | identical: %s\n", nm,
              abs(x$ll - y$ll) / abs(x$ll),
              paste(names(same), ifelse(same, "yes", "NO"), sep = "=",
                    collapse = " ")))
  if (!all(same)) {
    for (k in names(same)[!same]) {
      cat("   ", k, "base only:", paste(setdiff(unlist(x[[k]]),
                                                unlist(y[[k]])),
                                        collapse = " "), "\n")
      cat("   ", k, "lane only:", paste(setdiff(unlist(y[[k]]),
                                                unlist(x[[k]])),
                                        collapse = " "), "\n")
    }
  }
}
sink()
